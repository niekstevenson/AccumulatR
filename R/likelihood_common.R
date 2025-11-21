`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.integrate_rel_tol <- function() getOption("uuber.integrate.rel.tol", 1e-5)
.integrate_abs_tol <- function() getOption("uuber.integrate.abs.tol", 1e-6)

.native_integrate <- function(fn, lower, upper,
                              rel.tol = .integrate_rel_tol(),
                              abs.tol = .integrate_abs_tol(),
                              max.depth = 10L) {
  if (!is.finite(lower) || !is.finite(upper)) {
    res <- tryCatch(
      stats::integrate(fn,
                       lower = lower,
                       upper = upper,
                       rel.tol = rel.tol,
                       abs.tol = abs.tol,
                       stop.on.error = FALSE),
      error = function(e) list(value = NA_real_)
    )
    val <- as.numeric(res[["value"]])
    if (!is.finite(val)) val <- 0.0
    return(val)
  }
  val <- boost_integrate_cpp(
    fn,
    as.numeric(lower),
    as.numeric(upper),
    rel.tol,
    abs.tol,
    as.integer(max.depth)
  )
  if (!is.finite(val)) val <- 0.0
  val
}

.expr_is_event <- function(expr) {
  !is.null(expr) && !is.null(expr[["kind"]]) && identical(expr[["kind"]], "event")
}

.collect_guards <- function(expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  if (identical(kind, "guard")) return(list(expr))
  if (kind %in% c("and", "or")) {
    return(unlist(lapply(expr[["args"]], .collect_guards), recursive = FALSE))
  }
  if (identical(kind, "not")) {
    return(.collect_guards(expr[["arg"]]))
  }
  list()
}

.guards_conflict <- function(primary_guard, other_expr) {
  primary_blocker <- primary_guard[["blocker"]]
  primary_reference <- primary_guard[["reference"]]
  if (!.expr_is_event(primary_blocker) || !.expr_is_event(primary_reference)) return(FALSE)
  if (length(primary_guard[["unless"]] %||% list()) > 0) return(FALSE)
  guards <- .collect_guards(other_expr)
  if (length(guards) == 0) return(FALSE)
  blk_src <- primary_blocker[["source"]]
  ref_src <- primary_reference[["source"]]
  for (g in guards) {
    if (length(g[["unless"]] %||% list()) > 0) next
    g_blocker <- g[["blocker"]]
    g_reference <- g[["reference"]]
    if (!.expr_is_event(g_blocker) || !.expr_is_event(g_reference)) next
    if (identical(blk_src, g_reference[["source"]]) &&
        identical(ref_src, g_blocker[["source"]])) {
      return(TRUE)
    }
  }
  FALSE
}

.expr_signature <- function(expr) {
  cached <- attr(expr, ".lik_signature", exact = TRUE)
  if (!is.null(cached)) return(cached)
  if (is.null(expr) || is.null(expr[["kind"]])) return("null")
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]] %||% ""
    kval <- expr[["k"]] %||% ""
    return(sprintf("event:%s:%s", as.character(src), as.character(kval)))
  }
  if (kind %in% c("and", "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0L) {
      return(sprintf("%s:empty", kind))
    }
    parts <- vapply(args, .expr_signature, character(1))
    return(sprintf("%s[%s]", kind, paste(parts, collapse = ",")))
  }
  if (identical(kind, "guard")) {
    blk_sig <- .expr_signature(expr[["blocker"]])
    ref_sig <- .expr_signature(expr[["reference"]])
    unl <- expr[["unless"]] %||% list()
    unl_sig <- if (length(unl) == 0L) "" else paste(vapply(unl, .expr_signature, character(1)), collapse = ";")
    return(sprintf("guard[%s|%s|%s]", blk_sig, ref_sig, unl_sig))
  }
  if (identical(kind, "not")) {
    return(sprintf("not(%s)", .expr_signature(expr[["arg"]])))
  }
  sprintf("other:%s", kind)
}

.compile_once <- function(fn) {
  if (!is.function(fn)) return(NULL)
  if (!isTRUE(getOption("uuber.compile.likelihood", TRUE))) return(fn)
  if (!requireNamespace("compiler", quietly = TRUE)) return(fn)
  compiler::cmpfun(fn)
}

.make_scenario_record <- function(prep, weight, forced_complete, forced_survive) {
  if (!is.finite(weight) || weight <= 0) return(NULL)
  list(
    weight = as.numeric(weight),
    forced_complete = .coerce_forced_ids(prep, forced_complete),
    forced_survive = .coerce_forced_ids(prep, forced_survive)
  )
}

.prep_native_payload <- function(prep) {
  keep <- c(
    "accumulators",
    "pools",
    "components",
    "outcomes",
    "shared_triggers",
    "default_deadline"
  )
  payload <- prep[intersect(keep, names(prep))]
  payload$default_deadline <- payload$default_deadline %||% prep$default_deadline %||% Inf

  expr_compiled <- prep[[".expr_compiled"]] %||% .prep_expr_compiled(prep)
  expr_compiled <- .trim_expr_compiled(expr_compiled)
  id_index <- prep[[".id_index"]] %||% .prep_id_index(prep)
  competitor_map <- prep[[".competitors"]] %||% .prep_competitors(prep)
  runtime <- prep[[".runtime"]]
  runtime_info <- list()
  if (!is.null(runtime) && !is.null(runtime$cache_bundle)) {
    cache_stub <- .cache_bundle_native_stub(runtime$cache_bundle)
    if (!is.null(cache_stub)) {
      runtime_info$cache_bundle <- cache_stub
    }
  }
  if (!is.null(expr_compiled)) runtime_info$expr_compiled <- expr_compiled
  if (!is.null(id_index)) runtime_info$id_index <- id_index
  if (!is.null(competitor_map)) runtime_info$competitor_map <- competitor_map
  if (length(runtime_info) == 0L) runtime_info <- NULL

  if (!is.null(expr_compiled)) payload[[".expr_compiled"]] <- expr_compiled
  if (!is.null(id_index)) payload[[".id_index"]] <- id_index
  if (!is.null(competitor_map)) payload[[".competitors"]] <- competitor_map
  if (!is.null(runtime_info)) {
    payload_runtime <- runtime_info
  } else {
    payload_runtime <- NULL
  }
  payload_no_runtime <- payload
  payload_no_runtime[[".runtime"]] <- NULL
  payload_no_runtime <- .strip_functions(payload_no_runtime)
  if (!is.null(payload_runtime)) {
    payload_no_runtime[[".runtime"]] <- payload_runtime
  }
  payload <- payload_no_runtime
  payload
}

.trim_expr_compiled <- function(expr_compiled) {
  if (is.null(expr_compiled)) return(NULL)
  nodes <- expr_compiled[["nodes"]] %||% list()
  if (length(nodes) == 0) return(expr_compiled)
  drop_fields <- c(
    "density_fn", "surv_fn", "cdf_fn", "scenario_fn",
    "density_fast_fn", "surv_fast_fn", "cdf_fast_fn",
    "ops"
  )
  expr_copy <- unserialize(serialize(expr_compiled, NULL))
  expr_copy$nodes <- lapply(expr_copy$nodes, function(node) {
    if (is.null(node)) return(node)
    node[setdiff(names(node), drop_fields)]
  })
  expr_copy
}

.strip_functions <- function(x) {
  if (is.function(x) || is.environment(x)) return(NULL)
  if (is.pairlist(x)) x <- as.list(x)
  if (is.list(x)) {
    atts <- attributes(x)
    if (length(x) == 0L) return(x)
    x <- lapply(x, .strip_functions)
    keep <- !vapply(x, is.null, logical(1))
    x <- x[keep]
    if (!is.null(atts)) {
      atts$names <- names(x)
      attributes(x) <- atts
    }
  }
  x
}

.cache_bundle_native_stub <- function(bundle) {
  if (is.null(bundle)) return(NULL)
  stub <- list()
  native_ctx <- bundle$native_ctx
  if (!is.environment(native_ctx)) {
    native_ctx <- new.env(parent = emptyenv())
    native_ctx$ptr <- NULL
    native_ctx$proto <- raw(0)
  }
  stub$native_ctx <- native_ctx
  stub$guard_quadrature_limit <- bundle$guard_quadrature_limit %||% 0L
  if (!is.null(bundle$version)) {
    stub$version <- bundle$version
  }
  class(stub) <- attr(bundle, "class")
  stub
}

.prep_runtime_get <- function(prep, key, default = NULL) {
  runtime <- prep[[".runtime"]]
  if (is.null(runtime)) return(default)
  val <- runtime[[key]]
  if (is.null(val)) default else val
}

.prep_expr_compiled <- function(prep) .prep_runtime_get(prep, "expr_compiled")
.prep_label_cache <- function(prep) .prep_runtime_get(prep, "label_cache")
.prep_competitors <- function(prep) .prep_runtime_get(prep, "competitor_map")
.prep_id_index <- function(prep) .prep_runtime_get(prep, "id_index")
.prep_cache_bundle <- function(prep) .prep_runtime_get(prep, "cache_bundle")

.prep_native_context <- function(prep) {
  bundle <- .prep_cache_bundle(prep)
  native_env <- NULL
  if (!is.null(bundle) && !is.null(bundle$native_ctx) && is.environment(bundle$native_ctx)) {
    native_env <- bundle$native_ctx
  }
  if (!is.null(bundle) && is.null(native_env)) {
    native_env <- new.env(parent = emptyenv())
    native_env$ptr <- NULL
    native_env$proto <- raw(0)
    bundle$native_ctx <- native_env
  }
  proto_blob <- if (!is.null(native_env)) native_env$proto %||% raw(0) else raw(0)
  native_ptr <- if (!is.null(native_env)) native_env$ptr %||% NULL else NULL
  touch_ctx <- function(ptr) {
    if (!inherits(ptr, "externalptr")) return(ptr)
    try(native_context_touch_cpp(ptr), silent = TRUE)
    ptr
  }
  if (inherits(native_ptr, "externalptr")) {
    return(touch_ctx(native_ptr))
  }
  if (length(proto_blob) > 0L) {
    native_ptr <- tryCatch(native_context_from_proto_cpp(proto_blob), error = function(e) NULL)
    if (inherits(native_ptr, "externalptr")) {
      if (!is.null(native_env)) {
        native_env$ptr <- native_ptr
      } else {
        attr(native_ptr, "native_proto") <- proto_blob
      }
      return(touch_ctx(native_ptr))
    }
  }
  payload <- .prep_native_payload(prep)
  native_ptr <- native_context_build(payload)
  if (inherits(native_ptr, "externalptr")) {
    proto_blob <- tryCatch(native_prep_serialize_cpp(payload), error = function(e) raw(0))
    if (!is.null(native_env)) {
      native_env$ptr <- native_ptr
      native_env$proto <- proto_blob
    } else {
      attr(native_ptr, "native_proto") <- proto_blob
    }
  }
  touch_ctx(native_ptr)
}

.cache_component_key <- function(component) {
  if (is.null(component) || length(component) == 0L) {
    "__default__"
  } else {
    comp_chr <- as.character(component)
    if (length(comp_chr) == 0L) "__default__" else comp_chr[[1]] %||% "__default__"
  }
}

.canonical_named_list <- function(x) {
  if (!is.list(x)) return(x)
  nm <- names(x)
  if (is.null(nm)) return(lapply(x, .canonical_named_list))
  ord <- order(nm)
  res <- lapply(x[ord], .canonical_named_list)
  names(res) <- nm[ord]
  res
}

.structure_hash_payload <- function(prep) {
  accs <- prep$accumulators %||% list()
  acc_ids <- sort(names(accs))
  acc_payload <- lapply(acc_ids, function(id) {
    acc <- accs[[id]] %||% list()
    list(
      id = acc$id %||% id,
      dist = acc$dist %||% "",
      components = sort(acc$components %||% character(0)),
      tags = .canonical_named_list(acc$tags %||% list()),
      attrs = .canonical_named_list(acc$attrs %||% list()),
      shared_trigger_id = acc$shared_trigger_id %||% NA_character_,
      shared_trigger_q = acc$shared_trigger_q %||% NA_real_
    )
  })
  pools <- prep$pools %||% list()
  pool_ids <- sort(names(pools))
  pool_payload <- lapply(pool_ids, function(id) {
    pool <- pools[[id]] %||% list()
    list(
      id = pool$id %||% id,
      members = sort(pool$members %||% character(0)),
      k = as.integer(pool$k %||% 1L),
      attrs = .canonical_named_list(pool$attrs %||% list())
    )
  })
  outcomes <- prep$outcomes %||% list()
  outcome_ids <- sort(names(outcomes))
  outcome_payload <- lapply(outcome_ids, function(id) {
    def <- outcomes[[id]] %||% list()
    sig <- .expr_signature(def$expr %||% list())
    list(
      label = id,
      signature = sig,
      options = .canonical_named_list(def$options %||% list())
    )
  })
  shared <- prep$shared_triggers %||% list()
  shared_ids <- sort(names(shared))
  shared_payload <- lapply(shared_ids, function(id) {
    trig <- shared[[id]] %||% list()
    list(
      id = trig$id %||% id,
      q = trig$q %||% NA_real_,
      members = sort(trig$members %||% character(0))
    )
  })
  list(
    accumulators = acc_payload,
    pools = pool_payload,
    outcomes = outcome_payload,
    components = .canonical_named_list(prep$components %||% list()),
    default_deadline = prep$default_deadline %||% NA_real_,
    special_outcomes = .canonical_named_list(prep$special_outcomes %||% list()),
    shared_triggers = shared_payload
  )
}

.structure_hash_value <- function(prep) {
  payload <- .structure_hash_payload(prep)
  raw_bytes <- serialize(payload, connection = NULL, ascii = TRUE)
  paste(as.character(raw_bytes), collapse = "")
}

.na_cache_limit_value <- function() {
  opt <- getOption("uuber.cache.na.max_per_trial", 128L)
  if (!is.numeric(opt) || length(opt) == 0L || is.na(opt[[1]])) return(128L)
  val <- as.integer(opt[[1]])
  if (!is.finite(val) || val < 0L) return(0L)
  val
}

.build_likelihood_cache_bundle <- function(prep) {
  native_ctx_env <- new.env(parent = emptyenv())
  native_ctx_env$ptr <- NULL
  native_ctx_env$proto <- raw(0)
  native_ctx_env$ptr <- tryCatch(native_context_build(prep), error = function(e) NULL)
  native_ctx_env$proto <- tryCatch(native_prep_serialize_cpp(prep), error = function(e) raw(0))
  structure(
    list(
      node_plan = prep[[".expr_compiled"]],
      precomputed_values = new.env(parent = emptyenv(), hash = TRUE),
      pool_templates = new.env(parent = emptyenv(), hash = TRUE),
      guard_quadrature = new.env(parent = emptyenv(), hash = TRUE),
      guard_quadrature_meta = new.env(parent = emptyenv(), hash = TRUE),
      guard_quadrature_limit = .na_cache_limit_value(),
      native_ctx = native_ctx_env,
      version = Sys.time()
    ),
    class = "likelihood_cache_bundle"
  )
}

.prep_set_cache_bundle <- function(prep, bundle) {
  runtime <- prep[[".runtime"]]
  if (is.null(runtime)) return(prep)
  runtime$cache_bundle <- bundle
  prep[[".runtime"]] <- runtime
  prep
}

likelihood_build_native_bundle <- function(model_spec = NULL, prep = NULL) {
  if (is.null(prep)) {
    if (is.null(model_spec)) {
      stop("Provide either a prepared likelihood object or a model_spec", call. = FALSE)
    }
    prep <- .prepare_model_for_likelihood(model_spec)
  }
  if (is.null(prep[[".runtime"]]) || is.null(prep$.runtime$cache_bundle)) {
    prep <- .prep_set_cache_bundle(prep, .build_likelihood_cache_bundle(prep))
  }
  proto_blob <- native_prep_serialize_cpp(prep)
  ctx_ptr <- native_context_from_proto_cpp(proto_blob)
  list(
    prep = prep,
    proto = proto_blob,
    context = ctx_ptr
  )
}

.native_node_batch_eval <- function(prep, node_ids, times, component = NULL,
                                    forced_complete = integer(0),
                                    forced_survive = integer(0)) {
  node_ids <- as.integer(node_ids)
  times <- as.numeric(times)
  if (length(node_ids) == 0L || length(times) == 0L) return(list())
  native_ctx <- .prep_native_context(prep)
  tasks <- lapply(node_ids, function(nid) {
    list(
      node_id = as.integer(nid),
      times = times,
      component = component,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  })
  native_likelihood_eval_cpp(native_ctx, tasks)
}

.likelihood_cache_bundle_clone <- function(bundle) {
  if (is.null(bundle)) return(NULL)
  structure(
    list(
      node_plan = bundle$node_plan,
      precomputed_values = {
        env <- new.env(parent = emptyenv(), hash = TRUE)
        if (!is.null(bundle$precomputed_values)) {
          keys <- ls(bundle$precomputed_values, all.names = TRUE)
          for (k in keys) {
            env[[k]] <- bundle$precomputed_values[[k]]
          }
        }
        env
      },
      pool_templates = {
        env <- new.env(parent = emptyenv(), hash = TRUE)
        if (!is.null(bundle$pool_templates)) {
          keys <- ls(bundle$pool_templates, all.names = TRUE)
          for (k in keys) {
            env[[k]] <- bundle$pool_templates[[k]]
          }
        }
        env
      },
      guard_quadrature = {
        env <- new.env(parent = emptyenv(), hash = TRUE)
        if (!is.null(bundle$guard_quadrature)) {
          keys <- ls(bundle$guard_quadrature, all.names = TRUE)
          for (k in keys) {
            env[[k]] <- bundle$guard_quadrature[[k]]
          }
        }
        env
      },
      guard_quadrature_meta = {
        env <- new.env(parent = emptyenv(), hash = TRUE)
        if (!is.null(bundle$guard_quadrature_meta)) {
          keys <- ls(bundle$guard_quadrature_meta, all.names = TRUE)
          for (k in keys) {
            env[[k]] <- bundle$guard_quadrature_meta[[k]]
          }
        }
        env
      },
      guard_quadrature_limit = bundle$guard_quadrature_limit %||% .na_cache_limit_value(),
      native_ctx = {
        env <- new.env(parent = emptyenv())
        if (!is.null(bundle$native_ctx) && is.environment(bundle$native_ctx)) {
          env$ptr <- bundle$native_ctx$ptr %||% NULL
          env$proto <- bundle$native_ctx$proto %||% raw(0)
        } else {
          env$ptr <- NULL
          env$proto <- raw(0)
        }
        env
      },
      version = bundle$version
    ),
    class = "likelihood_cache_bundle"
  )
}

.likelihood_outcome_cache_key <- function(bundle_key, params_hash, outcome_label, rt_val) {
  outcome_chr <- if (is.null(outcome_label) || length(outcome_label) == 0L || is.na(outcome_label[[1]])) {
    "NA"
  } else {
    as.character(outcome_label)[[1]]
  }
  hash_key <- params_hash %||% "__base__"
  if (length(hash_key) == 0L || is.na(hash_key[[1]]) || !nzchar(hash_key[[1]])) {
    hash_key <- "__base__"
  } else {
    hash_key <- as.character(hash_key[[1]])
  }
  rt_key <- .eval_state_time_key(rt_val)
  paste(bundle_key, hash_key, outcome_chr, rt_key, sep = "|")
}

.bundle_precomputed_get <- function(bundle, key) {
  if (is.null(bundle) || is.null(bundle$precomputed_values) || !nzchar(key)) return(NULL)
  bundle$precomputed_values[[key]]
}

.bundle_precomputed_set <- function(bundle, key, value) {
  if (is.null(bundle) || is.null(bundle$precomputed_values) || !nzchar(key)) return(value)
  bundle$precomputed_values[[key]] <- value
  value
}

.pool_template_cache_key <- function(pool_id, component, k) {
  comp_key <- .cache_component_key(component)
  paste(pool_id %||% "", comp_key, as.integer(k) %||% 0L, sep = "|")
}

.bundle_pool_templates_get <- function(bundle, key) {
  if (is.null(bundle) || is.null(bundle$pool_templates) || !nzchar(key)) return(NULL)
  bundle$pool_templates[[key]]
}

.bundle_pool_templates_set <- function(bundle, key, value) {
  if (is.null(bundle) || is.null(bundle$pool_templates) || !nzchar(key)) return(value)
  bundle$pool_templates[[key]] <- value
  value
}

.na_cache_component_key <- function(component_key) {
  if (is.null(component_key) || length(component_key) == 0L) return("__default__")
  comp_chr <- as.character(component_key)[[1]]
  if (!nzchar(comp_chr) || is.na(comp_chr)) "__default__" else comp_chr
}

.na_cache_limit <- function(bundle = NULL) {
  limit <- NULL
  if (!is.null(bundle)) limit <- bundle$guard_quadrature_limit
  if (is.null(limit) || !is.numeric(limit) || length(limit) == 0L || is.na(limit[[1]])) {
    limit <- .na_cache_limit_value()
    if (!is.null(bundle)) bundle$guard_quadrature_limit <- limit
  }
  limit <- as.integer(limit[[1]])
  if (!is.finite(limit) || limit < 0L) return(0L)
  limit
}

.na_cache_touch <- function(bundle, component_key, cache_key) {
  if (is.null(bundle) || is.null(bundle$guard_quadrature) || !nzchar(cache_key)) return(invisible(NULL))
  comp_key <- .na_cache_component_key(component_key)
  order_env <- bundle$guard_quadrature_meta
  if (is.null(order_env) || !is.environment(order_env)) {
    order_env <- new.env(parent = emptyenv(), hash = TRUE)
    bundle$guard_quadrature_meta <- order_env
  }
  order_vec <- order_env[[comp_key]]
  if (is.null(order_vec)) order_vec <- character(0)
  if (length(order_vec) > 0L) {
    order_vec <- order_vec[order_vec != cache_key]
  }
  order_vec <- c(order_vec, cache_key)
  limit <- .na_cache_limit(bundle)
  if (limit <= 0L) {
    if (length(order_vec) > 0L) {
      for (drop_key in order_vec) {
        if (exists(drop_key, envir = bundle$guard_quadrature, inherits = FALSE)) {
          rm(list = drop_key, envir = bundle$guard_quadrature)
        }
      }
    }
    order_env[[comp_key]] <- character(0)
    return(invisible(NULL))
  }
  while (length(order_vec) > limit) {
    drop_key <- order_vec[[1]]
    order_vec <- order_vec[-1]
    if (exists(drop_key, envir = bundle$guard_quadrature, inherits = FALSE)) {
      rm(list = drop_key, envir = bundle$guard_quadrature)
    }
  }
  order_env[[comp_key]] <- order_vec
  invisible(NULL)
}

.cache_bundle_ensure <- function(prep) {
  bundle <- .prep_cache_bundle(prep)
  if (is.null(bundle)) {
    .cache_metrics_inc("struct_misses")
    bundle <- .build_likelihood_cache_bundle(prep)
    prep <- .prep_set_cache_bundle(prep, bundle)
  } else {
    .cache_metrics_inc("struct_hits")
  }
  list(prep = prep, bundle = bundle)
}

likelihood_reset_cache <- function(prep) {
  if (is.null(prep)) return(prep)
  runtime <- prep[[".runtime"]]
  if (is.null(runtime)) return(prep)
  .prep_set_cache_bundle(prep, .build_likelihood_cache_bundle(prep))
}

.likelihood_outcome_cache_key <- function(component_key, params_hash, outcome_label, rt_val) {
  outcome_chr <- if (is.null(outcome_label) || length(outcome_label) == 0L || is.na(outcome_label[[1]])) {
    "NA"
  } else {
    as.character(outcome_label)[[1]]
  }
  hash_key <- params_hash %||% "__base__"
  if (length(hash_key) == 0L || is.na(hash_key[[1]]) || !nzchar(hash_key[[1]])) {
    hash_key <- "__base__"
  } else {
    hash_key <- as.character(hash_key[[1]])
  }
  rt_key <- .eval_state_time_key(rt_val)
  paste(component_key, hash_key, outcome_chr, rt_key, sep = "|")
}

.likelihood_outcome_cached <- function(prep, cache_key, compute_fn) {
  info <- .cache_bundle_ensure(prep)
  prep <- info$prep
  bundle <- info$bundle
  cached <- .bundle_precomputed_get(bundle, cache_key)
  if (is.null(cached)) {
    cached <- compute_fn()
    .bundle_precomputed_set(bundle, cache_key, cached)
  }
  list(value = cached, prep = prep)
}

.params_hash_from_rows <- function(rows) {
  if (is.null(rows) || !nrow(rows)) return("__base__")
  if ("params_hash" %in% names(rows)) {
    vals <- rows$params_hash
    vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
    if (length(vals) > 0L) {
      hash_val <- as.character(vals[[1]])
      if (nzchar(hash_val)) return(hash_val)
    }
  }
  tmp <- rows
  if ("params_hash" %in% names(tmp)) tmp$params_hash <- NULL
  hash_val <- .param_rows_hash_value(tmp)
  if (!nzchar(hash_val)) "__base__" else hash_val
}

.na_integral_cache_key <- function(source_signature, component_key, upper_limit, competitor_signature,
                                   params_hash = "__base__") {
  upper_tag <- .eval_state_time_key(upper_limit)
  hash_tag <- params_hash %||% "__base__"
  if (!nzchar(hash_tag) || is.na(hash_tag)) hash_tag <- "__base__"
  paste("na_map", source_signature, component_key, hash_tag, upper_tag, competitor_signature, sep = "|")
}

.na_integral_fetch <- function(prep, component_key, cache_key) {
  bundle <- .prep_cache_bundle(prep)
  if (is.null(bundle)) {
    .cache_metrics_inc("na_misses")
    return(NULL)
  }
  val <- bundle$guard_quadrature[[cache_key]]
  if (!is.null(val)) {
    .cache_metrics_inc("na_hits")
    .na_cache_touch(bundle, component_key, cache_key)
  } else {
    .cache_metrics_inc("na_misses")
  }
  val
}

.na_integral_store <- function(prep, component_key, cache_key, value) {
  info <- .cache_bundle_ensure(prep)
  bundle <- info$bundle
  if (.na_cache_limit(bundle) <= 0L) {
    if (exists(cache_key, envir = bundle$guard_quadrature, inherits = FALSE)) {
      rm(list = cache_key, envir = bundle$guard_quadrature)
    }
    .na_cache_touch(bundle, component_key, cache_key)
    return(invisible(NULL))
  }
  bundle$guard_quadrature[[cache_key]] <- value
  .na_cache_touch(bundle, component_key, cache_key)
  invisible(NULL)
}
.inspect_likelihood_plan <- function(prep, include_cache = FALSE) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) {
    if (!isTRUE(include_cache)) return(data.frame())
    bundle <- .prep_cache_bundle(prep)
    guard_meta <- if (is.null(bundle) || is.null(bundle$guard_quadrature_meta)) list() else {
      keys <- ls(bundle$guard_quadrature_meta, all.names = TRUE)
      stats::setNames(lapply(keys, function(k) bundle$guard_quadrature_meta[[k]]), keys)
    }
    return(list(nodes = data.frame(), cache = list(
      precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
      pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
      guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE),
      guard_quadrature_orders = guard_meta
    )))
  }
  nodes <- comp$nodes %||% list()
  node_df <- do.call(rbind, lapply(nodes, function(node) {
    data.frame(
      id = node$id,
      kind = node$kind,
      needs_forced = isTRUE(node$needs_forced),
      scenario_sensitive = isTRUE(node$scenario_sensitive),
      sources = paste(node$sources %||% integer(0), collapse = ","),
      args = paste(node$args %||% integer(0), collapse = ","),
      has_fast_density = is.function(node$density_fast_fn),
      has_fast_survival = is.function(node$surv_fast_fn),
      has_fast_cdf = is.function(node$cdf_fast_fn),
      stringsAsFactors = FALSE
    )
  }))
  if (!isTRUE(include_cache)) return(node_df)
  bundle <- .prep_cache_bundle(prep)
  guard_meta <- if (is.null(bundle) || is.null(bundle$guard_quadrature_meta)) list() else {
    keys <- ls(bundle$guard_quadrature_meta, all.names = TRUE)
    stats::setNames(lapply(keys, function(k) bundle$guard_quadrature_meta[[k]]), keys)
  }
  cache_summary <- list(
    precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
    pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
    guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE),
    guard_quadrature_orders = guard_meta
  )
  list(nodes = node_df, cache = cache_summary)
}
