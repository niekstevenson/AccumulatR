`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.integrate_rel_tol <- function() getOption("uuber.integrate.rel.tol", 1e-5)
.integrate_abs_tol <- function() getOption("uuber.integrate.abs.tol", 1e-6)
.eval_state_create <- function(parent = emptyenv()) {
  new.env(parent = parent, hash = FALSE)
}

.expr_sources <- function(expr, prep) {
  gather_labels <- function(ex) {
    if (is.null(ex) || is.null(ex[["kind"]])) {
      return(character(0))
    }
    kind <- ex[["kind"]]
    if (identical(kind, "event")) {
      source_id <- ex[["source"]]
      if (!is.null(prep[["pools"]][[source_id]])) {
        members <- prep[["pools"]][[source_id]][["members"]] %||% character(0)
        unique(unlist(lapply(members, function(mid) {
          if (!is.null(prep[["accumulators"]][[mid]])) {
            return(mid)
          }
          if (!is.null(prep[["pools"]][[mid]])) {
            return(gather_labels(list(kind = "event", source = mid)))
          }
          character(0)
        })))
      } else {
        source_id
      }
    } else if (kind %in% c("and", "or")) {
      unique(unlist(lapply(ex[["args"]], gather_labels)))
    } else if (identical(kind, "guard")) {
      ref_src <- gather_labels(ex[["reference"]])
      blk_src <- gather_labels(ex[["blocker"]])
      unless_src <- unique(unlist(lapply(ex[["unless"]] %||% list(), gather_labels)))
      unique(c(ref_src, blk_src, unless_src))
    } else if (identical(kind, "not")) {
      gather_labels(ex[["arg"]])
    } else {
      character(0)
    }
  }
  .labels_to_ids(prep, gather_labels(expr))
}

.expr_is_event <- function(expr) {
  !is.null(expr) && !is.null(expr[["kind"]]) && identical(expr[["kind"]], "event")
}

.collect_guards <- function(expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) {
    return(list())
  }
  kind <- expr[["kind"]]
  if (identical(kind, "guard")) {
    return(list(expr))
  }
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
  if (!.expr_is_event(primary_blocker) || !.expr_is_event(primary_reference)) {
    return(FALSE)
  }
  if (length(primary_guard[["unless"]] %||% list()) > 0) {
    return(FALSE)
  }
  guards <- .collect_guards(other_expr)
  if (length(guards) == 0) {
    return(FALSE)
  }
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
  if (!is.null(cached)) {
    return(cached)
  }
  if (is.null(expr) || is.null(expr[["kind"]])) {
    return("null")
  }
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

.prep_native_payload <- function(prep) {
  keep <- c(
    "accumulators",
    "pools",
    "components",
    "outcomes",
    "shared_triggers"
  )
  payload <- prep[intersect(keep, names(prep))]

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
  if (is.null(expr_compiled)) {
    return(NULL)
  }
  nodes <- expr_compiled[["nodes"]] %||% list()
  if (length(nodes) == 0) {
    return(expr_compiled)
  }
  drop_fields <- c(
    "density_fn", "surv_fn", "cdf_fn", "scenario_fn",
    "density_fast_fn", "surv_fast_fn", "cdf_fast_fn",
    "ops"
  )
  expr_copy <- unserialize(serialize(expr_compiled, NULL))
  expr_copy$nodes <- lapply(expr_copy$nodes, function(node) {
    if (is.null(node)) {
      return(node)
    }
    node[setdiff(names(node), drop_fields)]
  })
  expr_copy
}

.strip_functions <- function(x) {
  if (is.function(x) || is.environment(x)) {
    return(NULL)
  }
  if (is.pairlist(x)) x <- as.list(x)
  if (is.list(x)) {
    atts <- attributes(x)
    if (length(x) == 0L) {
      return(x)
    }
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

.precompile_likelihood_expressions <- function(prep) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  if (length(outcome_defs) == 0L) {
    prep[[".expr_compiled"]] <- NULL
    return(prep)
  }
  sig_env <- new.env(parent = emptyenv(), hash = TRUE)
  nodes <- list()
  next_id <- 1L

  compile_expr <- function(expr) {
    if (is.null(expr) || is.null(expr[["kind"]])) {
      return(expr)
    }
    sig <- .expr_signature(expr)
    existing_id <- sig_env[[sig]]
    if (!is.null(existing_id)) {
      attr(expr, ".lik_id") <- as.integer(existing_id)
      return(expr)
    }
    kind <- expr[["kind"]]
    child_args <- integer(0)
    reference_id <- NA_integer_
    blocker_id <- NA_integer_
    unless_ids <- integer(0)
    arg_id <- NA_integer_

    if (kind %in% c("and", "or")) {
      args <- expr[["args"]] %||% list()
      if (length(args) > 0L) {
        expr[["args"]] <- lapply(args, compile_expr)
        child_args <- vapply(expr[["args"]], function(a) {
          attr(a, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      }
    } else if (identical(kind, "guard")) {
      expr[["reference"]] <- compile_expr(expr[["reference"]])
      expr[["blocker"]] <- compile_expr(expr[["blocker"]])
      unless_list <- expr[["unless"]] %||% list()
      if (length(unless_list) > 0L) {
        expr[["unless"]] <- lapply(unless_list, compile_expr)
        unless_ids <- vapply(expr[["unless"]], function(unl) {
          attr(unl, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      }
      reference_id <- attr(expr[["reference"]], ".lik_id", exact = TRUE) %||% NA_integer_
      blocker_id <- attr(expr[["blocker"]], ".lik_id", exact = TRUE) %||% NA_integer_
    } else if (identical(kind, "not")) {
      expr[["arg"]] <- compile_expr(expr[["arg"]])
      arg_id <- attr(expr[["arg"]], ".lik_id", exact = TRUE) %||% NA_integer_
    }

    node <- list(
      id = next_id,
      kind = kind,
      expr = expr,
      sources = .expr_sources(expr, prep),
      args = if (length(child_args) > 0L) child_args else NULL,
      reference_id = reference_id,
      blocker_id = blocker_id,
      unless_ids = if (length(unless_ids) > 0L) unless_ids else NULL,
      arg = arg_id
    )

    if (identical(kind, "event")) {
      node$source <- expr[["source"]] %||% NULL
      node$needs_forced <- TRUE
      node$scenario_sensitive <- !is.null(prep[["pools"]][[node$source %||% ""]])
    } else if (identical(kind, "guard")) {
      node$needs_forced <- TRUE
      node$scenario_sensitive <- TRUE
    } else {
      child_nodes <- if (length(child_args) > 0L) nodes[as.integer(child_args)] else list()
      child_needs <- if (length(child_nodes) > 0L) {
        vapply(child_nodes, function(n) isTRUE(n$needs_forced), logical(1))
      } else {
        logical(0)
      }
      child_scen <- if (length(child_nodes) > 0L) {
        vapply(child_nodes, function(n) isTRUE(n$scenario_sensitive), logical(1))
      } else {
        logical(0)
      }
      node$needs_forced <- any(child_needs)
      node$scenario_sensitive <- any(child_scen)
    }

    sig_env[[sig]] <- next_id
    attr(expr, ".lik_id") <- next_id
    nodes[[next_id]] <<- node
    next_id <<- next_id + 1L
    expr
  }

  for (i in seq_along(outcome_defs)) {
    out_expr <- outcome_defs[[i]][["expr"]]
    outcome_defs[[i]][["expr"]] <- compile_expr(out_expr)
  }

  prep[["outcomes"]] <- outcome_defs
  prep[[".expr_compiled"]] <- list(nodes = nodes, signatures = sig_env)
  prep
}

.expr_lookup_compiled <- function(expr, prep) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) {
    return(NULL)
  }
  node_id <- attr(expr, ".lik_id", exact = TRUE)
  if (is.null(node_id) || is.na(node_id)) {
    sig <- .expr_signature(expr)
    node_id <- comp[["signatures"]][[sig]]
  }
  if (is.null(node_id) || is.na(node_id)) {
    return(NULL)
  }
  node_idx <- as.integer(node_id)
  nodes <- comp[["nodes"]] %||% list()
  if (length(nodes) < node_idx || node_idx < 1L) {
    return(NULL)
  }
  nodes[[node_idx]]
}

.cache_bundle_native_stub <- function(bundle) {
  if (is.null(bundle)) {
    return(NULL)
  }
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
  if (is.null(runtime)) {
    return(default)
  }
  val <- runtime[[key]]
  if (is.null(val)) default else val
}

.prep_expr_compiled <- function(prep) .prep_runtime_get(prep, "expr_compiled", prep[[".expr_compiled"]])
.prep_label_cache <- function(prep) .prep_runtime_get(prep, "label_cache", prep[[".label_cache"]])
.prep_competitors <- function(prep) .prep_runtime_get(prep, "competitor_map", prep[[".competitors"]])
.prep_id_index <- function(prep) .prep_runtime_get(prep, "id_index", prep[[".id_index"]])
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
    native_ptr <- tryCatch(native_context_from_proto_cpp(proto_blob), error = function(e) NULL)
    if (inherits(native_ptr, "externalptr")) {
      if (!is.null(native_env)) {
        native_env$ptr <- native_ptr
      } else {
        attr(native_ptr, "native_proto") <- proto_blob
      }
      return(native_ptr)
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
  native_ptr
}

.canonical_named_list <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  nm <- names(x)
  if (is.null(nm)) {
    return(lapply(x, .canonical_named_list))
  }
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
  if (!is.numeric(opt) || length(opt) == 0L || is.na(opt[[1]])) {
    return(128L)
  }
  val <- as.integer(opt[[1]])
  if (!is.finite(val) || val < 0L) {
    return(0L)
  }
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
  if (is.null(runtime)) {
    return(prep)
  }
  runtime$cache_bundle <- bundle
  prep[[".runtime"]] <- runtime
  prep
}

.prepare_competitor_map <- function(prep) {
  outcomes <- prep[["outcomes"]] %||% list()
  if (length(outcomes) == 0L) {
    return(list())
  }
  comps <- lapply(seq_along(outcomes), function(idx) {
    .build_competitor_exprs_at_index(prep, idx)
  })
  names(comps) <- names(outcomes)
  attr(comps, "by_index") <- TRUE
  comps
}

.build_competitor_exprs_at_index <- function(prep, target_idx) {
  outcome_defs <- prep[["outcomes"]]
  if (is.null(outcome_defs) || length(outcome_defs) == 0L) {
    return(list())
  }
  if (is.na(target_idx) || target_idx < 1 || target_idx > length(outcome_defs)) {
    return(list())
  }
  target_def <- outcome_defs[[target_idx]]
  target_expr <- target_def[["expr"]]
  if (is.null(target_expr)) {
    return(list())
  }
  target_label <- names(outcome_defs)[target_idx]
  target_opts <- target_def[["options"]] %||% list()
  target_map <- target_opts[["map_outcome_to"]] %||% target_label
  comp_indices <- setdiff(seq_along(outcome_defs), target_idx)
  if (length(comp_indices) == 0) {
    return(list())
  }
  comp_indices <- Filter(function(idx) {
    comp_def <- outcome_defs[[idx]]
    comp_opts <- comp_def[["options"]] %||% list()
    if (!is.null(comp_opts[["alias_of"]])) {
      return(FALSE)
    }
    expr_lbl <- comp_def[["expr"]]
    if (!is.null(expr_lbl[["kind"]]) && identical(expr_lbl[["kind"]], "event")) {
      src <- expr_lbl[["source"]]
      if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) {
        return(FALSE)
      }
    }
    comp_label <- names(outcome_defs)[idx]
    comp_target <- comp_opts[["map_outcome_to"]] %||% comp_label
    if (identical(as.character(comp_target), as.character(target_map))) {
      return(FALSE)
    }
    TRUE
  }, comp_indices)
  if (length(comp_indices) == 0) {
    return(list())
  }
  comps <- lapply(comp_indices, function(idx) outcome_defs[[idx]][["expr"]])
  target_sources <- .expr_sources(target_expr, prep)
  if (length(comps) > 0) {
    comps <- lapply(comps, function(comp_expr) {
      if (!is.null(comp_expr[["kind"]]) && identical(comp_expr[["kind"]], "guard")) {
        blocker_sources <- .expr_sources(comp_expr[["blocker"]], prep)
        if (length(target_sources) > 0 && all(target_sources %in% blocker_sources)) {
          return(comp_expr[["reference"]])
        }
      }
      comp_expr
    })
  }
  if (!is.null(target_expr[["kind"]]) && identical(target_expr[["kind"]], "guard") && length(comps) > 0) {
    blocker_sources <- .expr_sources(target_expr[["blocker"]], prep)
    keep_idx <- vapply(comps, function(comp_expr) {
      comp_sources <- .expr_sources(comp_expr, prep)
      length(comp_sources) == 0 || any(!comp_sources %in% blocker_sources)
    }, logical(1))
    comps <- comps[keep_idx]
  }
  if (length(comps) > 0 && isTRUE(getOption("uuber.filter_guard_competitors", TRUE))) {
    primary_guards <- .collect_guards(target_expr)
    if (length(primary_guards) > 0) {
      guard_keep <- vapply(seq_along(comps), function(i) {
        comp_expr <- comps[[i]]
        !any(vapply(primary_guards, function(pg) .guards_conflict(pg, comp_expr), logical(1)))
      }, logical(1))
      comps <- comps[guard_keep]
    }
  }
  comps
}

.build_competitor_exprs <- function(prep, target_label, target_expr) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  target_idx <- match(target_label, names(outcome_defs))
  if (is.na(target_idx)) {
    return(list())
  }
  .build_competitor_exprs_at_index(prep, target_idx)
}
