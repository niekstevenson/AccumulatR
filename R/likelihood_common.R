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
  native <- .lik_native_fn("boost_integrate_cpp")
  val <- native(
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
  if (inherits(native_ptr, "externalptr")) {
    return(native_ptr)
  }
  if (length(proto_blob) > 0L) {
    from_proto <- .lik_native_fn("native_context_from_proto_cpp")
    native_ptr <- tryCatch(from_proto(proto_blob), error = function(e) NULL)
    if (inherits(native_ptr, "externalptr")) {
      if (!is.null(native_env)) {
        native_env$ptr <- native_ptr
      } else {
        attr(native_ptr, "native_proto") <- proto_blob
      }
      return(native_ptr)
    }
  }
  builder <- .lik_native_fn("native_context_build")
  native_ptr <- builder(prep)
  if (inherits(native_ptr, "externalptr")) {
    serialize_fn <- .lik_native_fn("native_prep_serialize_cpp")
    proto_blob <- tryCatch(serialize_fn(prep), error = function(e) raw(0))
    if (!is.null(native_env)) {
      native_env$ptr <- native_ptr
      native_env$proto <- proto_blob
    } else {
      attr(native_ptr, "native_proto") <- proto_blob
    }
  }
  native_ptr
}

.cache_component_key <- function(component) {
  if (is.null(component) || length(component) == 0L) {
    "__default__"
  } else {
    comp_chr <- as.character(component)
    if (length(comp_chr) == 0L) "__default__" else comp_chr[[1]] %||% "__default__"
  }
}

.build_likelihood_cache_bundle <- function(prep) {
  native_ctx_env <- new.env(parent = emptyenv())
  native_ctx_env$ptr <- NULL
  native_ctx_env$proto <- raw(0)
  native_builder <- NULL
  native_proto_builder <- NULL
  native_ctx_env$ptr <- tryCatch({
    native_builder <- .lik_native_fn("native_context_build")
    native_builder(prep)
  }, error = function(e) NULL)
  native_ctx_env$proto <- tryCatch({
    native_proto_builder <- .lik_native_fn("native_prep_serialize_cpp")
    native_proto_builder(prep)
  }, error = function(e) raw(0))
  structure(
    list(
      node_plan = prep[[".expr_compiled"]],
      precomputed_values = new.env(parent = emptyenv(), hash = TRUE),
      pool_templates = new.env(parent = emptyenv(), hash = TRUE),
      guard_quadrature = new.env(parent = emptyenv(), hash = TRUE),
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
  serialize_fn <- .lik_native_fn("native_prep_serialize_cpp")
  proto_blob <- serialize_fn(prep)
  ctx_fn <- .lik_native_fn("native_context_from_proto_cpp")
  ctx_ptr <- ctx_fn(proto_blob)
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
  native_fn <- .lik_native_fn("native_likelihood_eval_cpp")
  tasks <- lapply(node_ids, function(nid) {
    list(
      node_id = as.integer(nid),
      times = times,
      component = component,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  })
  native_fn(native_ctx, tasks)
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

.likelihood_outcome_cache_key <- function(bundle_key, outcome_label, rt_val) {
  outcome_chr <- if (is.null(outcome_label) || length(outcome_label) == 0L || is.na(outcome_label[[1]])) {
    "NA"
  } else {
    as.character(outcome_label)[[1]]
  }
  rt_key <- .eval_state_time_key(rt_val)
  paste(bundle_key, outcome_chr, rt_key, sep = "|")
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

.cache_bundle_ensure <- function(prep) {
  bundle <- .prep_cache_bundle(prep)
  if (is.null(bundle)) {
    bundle <- .build_likelihood_cache_bundle(prep)
    prep <- .prep_set_cache_bundle(prep, bundle)
  }
  list(prep = prep, bundle = bundle)
}

.likelihood_outcome_cache_key <- function(component_key, outcome_label, rt_val) {
  outcome_chr <- if (is.null(outcome_label) || length(outcome_label) == 0L || is.na(outcome_label[[1]])) {
    "NA"
  } else {
    as.character(outcome_label)[[1]]
  }
  rt_key <- .eval_state_time_key(rt_val)
  paste(component_key, outcome_chr, rt_key, sep = "|")
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

.guard_integral_cache_key <- function(guard_signature, component_key, upper_limit, competitor_signature) {
  upper_tag <- .eval_state_time_key(upper_limit)
  paste("guard_int", guard_signature, component_key, upper_tag, competitor_signature, sep = "|")
}

.guard_integral_fetch <- function(prep, cache_key) {
  bundle <- .prep_cache_bundle(prep)
  if (is.null(bundle)) return(NULL)
  bundle$guard_quadrature[[cache_key]]
}

.guard_integral_store <- function(prep, cache_key, value) {
  info <- .cache_bundle_ensure(prep)
  bundle <- info$bundle
  bundle$guard_quadrature[[cache_key]] <- value
  invisible(NULL)
}
.inspect_likelihood_plan <- function(prep, include_cache = FALSE) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) {
    if (!isTRUE(include_cache)) return(data.frame())
    bundle <- .prep_cache_bundle(prep)
    return(list(nodes = data.frame(), cache = list(
      precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
      pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
      guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE)
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
  cache_summary <- list(
    precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
    pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
    guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE)
  )
  list(nodes = node_df, cache = cache_summary)
}
