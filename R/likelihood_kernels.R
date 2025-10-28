.precompile_likelihood_expressions <- function(prep) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  if (length(outcome_defs) == 0L) {
    prep[[".expr_compiled"]] <- NULL
    return(prep)
  }
  sig_env <- new.env(parent = emptyenv(), hash = TRUE)
  nodes_env <- new.env(parent = emptyenv())
  nodes_env$nodes <- list()
  
  compile_expr <- function(expr) {
    if (is.null(expr) || is.null(expr[["kind"]])) return(expr)
    sig <- .expr_signature(expr)
    existing_id <- sig_env[[sig]]
    if (!is.null(existing_id)) {
      attr(expr, ".lik_id") <- as.integer(existing_id)
      return(expr)
    }
    kind <- expr[["kind"]]
    child_info <- list()
    if (kind %in% c("and", "or")) {
      args <- expr[["args"]] %||% list()
      if (length(args) > 0L) {
        expr[["args"]] <- lapply(args, compile_expr)
        child_info$args <- vapply(expr[["args"]], function(arg) {
          attr(arg, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      } else {
        child_info$args <- integer(0)
      }
    } else if (identical(kind, "guard")) {
      expr[["reference"]] <- compile_expr(expr[["reference"]])
      expr[["blocker"]] <- compile_expr(expr[["blocker"]])
      unless_list <- expr[["unless"]] %||% list()
      if (length(unless_list) > 0L) {
        expr[["unless"]] <- lapply(unless_list, compile_expr)
        child_info$unless <- vapply(expr[["unless"]], function(unl) {
          attr(unl, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      } else {
        expr[["unless"]] <- list()
        child_info$unless <- integer(0)
      }
      child_info$reference <- attr(expr[["reference"]], ".lik_id", exact = TRUE) %||% NA_integer_
      child_info$blocker <- attr(expr[["blocker"]], ".lik_id", exact = TRUE) %||% NA_integer_
    } else if (identical(kind, "not")) {
      expr[["arg"]] <- compile_expr(expr[["arg"]])
      child_info$arg <- attr(expr[["arg"]], ".lik_id", exact = TRUE) %||% NA_integer_
    }
    node <- .build_compiled_node(expr, child_info, prep, nodes_env)
    node_id <- length(nodes_env$nodes) + 1L
    node$id <- node_id
    nodes_env$nodes[[node_id]] <- node
    sig_env[[sig]] <- node_id
    attr(expr, ".lik_id") <- node_id
    expr
  }
  
  for (lbl in names(outcome_defs)) {
    out_expr <- outcome_defs[[lbl]][["expr"]]
    outcome_defs[[lbl]][["expr"]] <- compile_expr(out_expr)
  }
  
  prep[["outcomes"]] <- outcome_defs
  prep[[".expr_compiled"]] <- list(
    nodes = nodes_env$nodes,
    signatures = sig_env
  )
  prep
}

.expr_lookup_compiled <- function(expr, prep) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) return(NULL)
  node_id <- attr(expr, ".lik_id", exact = TRUE)
  if (is.null(node_id) || is.na(node_id)) {
    sig <- .expr_signature(expr)
    node_id <- comp[["signatures"]][[sig]]
  }
  if (is.null(node_id) || is.na(node_id)) return(NULL)
  node <- comp[["nodes"]][[as.integer(node_id)]]
  if (is.null(node)) return(NULL)
  node
}

.build_compiled_node <- function(expr, child_info, prep, nodes_env) {
  kind <- expr[["kind"]]
  node <- list(
    kind = kind,
    sources = .expr_sources(expr, prep),
    expr = expr
  )
  if (identical(kind, "event")) {
    node$cdf_fn <- .compile_once(.make_event_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_event_survival_fn(expr, prep))
    node$scenario_fn <- .compile_once(.make_event_scenario_fn(expr, prep))
    return(node)
  }
  if (kind %in% c("and", "or")) {
    child_ids <- as.integer(child_info$args %||% integer(0))
    node$args <- child_ids
    if (length(child_ids) == 0L || any(is.na(child_ids))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$scenario_fn <- NULL
      return(node)
    }
    child_nodes <- lapply(child_ids, function(id) nodes_env$nodes[[id]])
    if (any(vapply(child_nodes, is.null, logical(1)))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$scenario_fn <- NULL
      return(node)
    }
    if (identical(kind, "and")) {
      have_child_cdf <- all(vapply(child_nodes, function(n) is.function(n$cdf_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      node$cdf_fn <- if (have_child_cdf) .compile_once(.make_and_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (is.function(node$cdf_fn)) .compile_once(.make_and_survival_fn(node$cdf_fn)) else NULL
      node$scenario_fn <- if (have_child_cdf && have_child_scen) .compile_once(.make_and_scenario_fn(child_nodes, prep)) else NULL
    } else {
      have_child_surv <- all(vapply(child_nodes, function(n) is.function(n$surv_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      node$cdf_fn <- if (have_child_surv) .compile_once(.make_or_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (have_child_surv) .compile_once(.make_or_survival_fn(child_nodes, prep)) else NULL
      node$scenario_fn <- if (have_child_scen) .compile_once(.make_or_scenario_fn(child_nodes, prep)) else NULL
    }
    return(node)
  }
  if (identical(kind, "guard")) {
    node$reference_id <- child_info$reference %||% NA_integer_
    node$blocker_id <- child_info$blocker %||% NA_integer_
    node$unless_ids <- as.integer(child_info$unless %||% integer(0))
    node$cdf_fn <- .compile_once(.make_guard_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_guard_survival_fn(expr, prep))
    guard_scen <- .make_guard_scenario_fn(expr, prep)
    node$scenario_fn <- if (is.function(guard_scen)) .compile_once(guard_scen) else NULL
    return(node)
  }
  node$cdf_fn <- NULL
  node$surv_fn <- NULL
  node$scenario_fn <- NULL
  node
}

.make_event_cdf_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    .event_cdf_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  }
}

.make_event_survival_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    .event_survival_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  }
}

.make_event_scenario_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  source_idx <- .label_to_id(prep, source_id)
  source_ids <- if (is.na(source_idx)) integer(0) else source_idx
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    weight <- .event_density_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    if (!is.finite(weight) || weight <= 0) return(list())
    pool_scenarios <- attr(weight, "scenarios")
    if (!is.null(pool_scenarios) && length(pool_scenarios) > 0L) {
      out <- list()
      for (psc in pool_scenarios) {
        if (is.null(psc)) next
        w <- psc$weight
        if (!is.finite(w) || w <= 0) next
        fcomp <- .forced_union(prep, forced_complete, c(source_ids, psc$forced_complete))
        fsurv <- .forced_union(prep, forced_survive, psc$forced_survive)
        sc <- .make_scenario_record(prep, w, fcomp, fsurv)
        if (!is.null(sc)) out[[length(out) + 1L]] <- sc
      }
      return(out)
    }
    sc <- .make_scenario_record(
      prep,
      weight,
      .forced_union(prep, forced_complete, source_ids),
      forced_survive
    )
    if (is.null(sc)) list() else list(sc)
  }
}

.make_and_cdf_fn <- function(child_nodes, prep) {
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fn")
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    if (length(child_cdf) == 0L) return(0.0)
    prod_val <- 1.0
    for (fn in child_cdf) {
      val <- fn(t, component, forced_complete, forced_survive)
      prod_val <- prod_val * val
      if (prod_val == 0) break
    }
    prod_val
  }
}

.make_and_survival_fn <- function(cdf_fn) {
  function(t, component, forced_complete, forced_survive) {
    val <- cdf_fn(t, component, forced_complete, forced_survive)
    if (!is.finite(val)) val <- 0.0
    out <- 1.0 - val
    if (out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    out
  }
}

.make_or_survival_fn <- function(child_nodes, prep) {
  child_surv <- lapply(child_nodes, `[[`, "surv_fn")
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    if (length(child_surv) == 0L) return(1.0)
    prod_val <- 1.0
    for (fn in child_surv) {
      val <- fn(t, component, forced_complete, forced_survive)
      prod_val <- prod_val * val
      if (prod_val == 0) break
    }
    prod_val
  }
}

.make_or_cdf_fn <- function(child_nodes, prep) {
  surv_fn <- .make_or_survival_fn(child_nodes, prep)
  function(t, component, forced_complete, forced_survive) {
    val <- surv_fn(t, component, forced_complete, forced_survive)
    out <- 1.0 - val
    if (!is.finite(out)) out <- 0.0
    if (out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    out
  }
}

.make_and_scenario_fn <- function(child_nodes, prep) {
  child_scen <- lapply(child_nodes, `[[`, "scenario_fn")
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fn")
  child_sources <- lapply(child_nodes, function(node) node$sources %||% integer(0))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    out <- list()
    n <- length(child_scen)
    if (n == 0L) return(out)
    for (i in seq_len(n)) {
      si <- child_scen[[i]](t, component, forced_complete, forced_survive)
      if (length(si) == 0L) next
      others <- seq_len(n)
      if (length(others) > 0L) others <- others[others != i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        ok <- TRUE
        if (length(others) > 0L) {
          for (j in others) {
            Fj <- child_cdf[[j]](t, component, fcomp, fsurv)
            if (!is.finite(Fj) || Fj <= 0) {
              ok <- FALSE
              break
            }
            weight <- weight * Fj
            if (!is.finite(weight) || weight <= 0) {
              ok <- FALSE
              break
            }
            fcomp <- .forced_union(prep, fcomp, child_sources[[j]])
          }
        }
        if (!ok) next
        sc_out <- .make_scenario_record(prep, weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    out
  }
}

.make_or_scenario_fn <- function(child_nodes, prep) {
  child_scen <- lapply(child_nodes, `[[`, "scenario_fn")
  child_sources <- lapply(child_nodes, function(node) node$sources %||% integer(0))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    out <- list()
    n <- length(child_scen)
    if (n == 0L) return(out)
    for (i in seq_len(n)) {
      si <- child_scen[[i]](t, component, forced_complete, forced_survive)
      if (length(si) == 0L) next
      others <- seq_len(n)
      if (length(others) > 0L) others <- others[others != i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        valid <- TRUE
        if (length(others) > 0L) {
          for (j in others) {
            req_j <- child_sources[[j]]
            if (length(req_j) == 0L) next
            if (all(req_j %in% fcomp)) {
              valid <- FALSE
              break
            }
            witness <- req_j[!(req_j %in% fcomp)]
            if (length(witness) >= 1L) {
              fsurv <- .forced_union(prep, fsurv, witness[[1]])
            } else {
              valid <- FALSE
              break
            }
          }
        }
        if (!valid) next
        sc_out <- .make_scenario_record(prep, weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    out
  }
}

.make_guard_cdf_fn <- function(expr, prep) {
  function(t, component, forced_complete, forced_survive) {
    if (t <= 0) return(0.0)
    if (!is.finite(t)) return(1.0)
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    local_cache <- lik_cache_create()
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        val <- .guard_density_at(
          expr, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = local_cache
        )
        as.numeric(val)
      }, numeric(1))
    }
    val <- tryCatch(
      stats::integrate(
        dens_fun,
        lower = 0,
        upper = t,
        rel.tol = .integrate_rel_tol(),
        abs.tol = .integrate_abs_tol(),
        stop.on.error = FALSE
      )[["value"]],
      error = function(e) 0.0
    )
    out <- as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
}

.make_guard_survival_fn <- function(expr, prep) {
  cdf_fn <- .make_guard_cdf_fn(expr, prep)
  function(t, component, forced_complete, forced_survive) {
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    val <- cdf_fn(t, component, forced_complete, forced_survive)
    out <- 1.0 - val
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
}

.guard_effective_survival <- function(guard_expr, t, prep, component,
                                      forced_complete, forced_survive,
                                      cache, block_node = NULL,
                                      unless_nodes = list()) {
  blocker <- guard_expr[["blocker"]]
  unless_list <- guard_expr[["unless"]] %||% list()
  if (is.null(blocker)) return(1.0)
  if (length(unless_list) == 0) {
    if (!is.null(block_node) && is.function(block_node$surv_fn)) {
      return(block_node$surv_fn(t, component, forced_complete, forced_survive))
    }
    return(.eval_expr_survival_cond(
      blocker, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    ))
  }
  blocker_density <- function(u) {
    vapply(u, function(ui) {
      if (!is.finite(ui) || ui < 0) return(0.0)
      .eval_expr_likelihood(
        blocker, ui, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        cache = cache
      )
    }, numeric(1))
  }
  protector_survival_product <- function(u) {
    if (length(unless_list) == 0) return(rep(1.0, length(u)))
    vapply(u, function(ui) {
      vals <- vapply(seq_along(unless_list), function(idx) {
        unl <- unless_list[[idx]]
        unl_node <- unless_nodes[[idx]]
        if (!is.null(unl_node) && is.function(unl_node$surv_fn)) {
          unl_node$surv_fn(ui, component, forced_complete, forced_survive)
        } else {
          .eval_expr_survival_cond(
            unl, ui, prep, component,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            cache = cache
          )
        }
      }, numeric(1))
      prod(vals)
    }, numeric(1))
  }
  tryCatch({
    if (!is.finite(t) || t <= 0) return(1.0)
    val <- stats::integrate(
      function(u) {
        blocker_density(u) * protector_survival_product(u)
      },
      lower = 0,
      upper = t,
      rel.tol = .integrate_rel_tol(),
      abs.tol = .integrate_abs_tol(),
      stop.on.error = FALSE
    )[["value"]]
    s <- 1.0 - as.numeric(val)
    if (!is.finite(s)) 0.0 else max(0.0, min(1.0, s))
  }, error = function(e) 1.0)
}

.guard_scenarios_slow <- function(expr, t, prep, component,
                                  forced_complete, forced_survive,
                                  cache) {
  reference <- expr[["reference"]]
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  ref_scen <- .expr_scenarios_at(
    reference, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    cache = cache
  )
  if (length(ref_scen) == 0) return(list())
  blocker_sources <- integer(0)
  if (!is.null(blocker)) {
    blocker_sources <- .expr_sources(blocker, prep)
    blocker_label <- .label_to_id(prep, blocker[["source"]])
    if (!is.na(blocker_label)) {
      blocker_sources <- .forced_union(prep, blocker_sources, blocker_label)
    }
  }
  block_node <- if (!is.null(blocker)) .expr_lookup_compiled(blocker, prep) else NULL
  unless_nodes <- lapply(unless_list, function(unl) .expr_lookup_compiled(unl, prep))
  out <- list()
  for (sc in ref_scen) {
    if (is.null(sc)) next
    fc_ctx <- .coerce_forced_ids(prep, sc$forced_complete)
    fs_ctx <- .coerce_forced_ids(prep, sc$forced_survive)
    S_eff <- .guard_effective_survival(
      expr, t, prep, component,
      fc_ctx, fs_ctx, cache,
      block_node = block_node,
      unless_nodes = unless_nodes
    )
    if (!is.finite(S_eff) || S_eff <= 0) next
    weight <- sc$weight * S_eff
    if (!is.finite(weight) || weight <= 0) next
    fcomp <- fc_ctx
    fsurv <- .forced_union(prep, fs_ctx, blocker_sources)
    out[[length(out) + 1L]] <- list(
      weight = as.numeric(weight),
      forced_complete = fcomp,
      forced_survive = .coerce_forced_ids(prep, fsurv)
    )
  }
  out
}

.guard_density_at <- function(expr, t, prep, component,
                              forced_complete, forced_survive,
                              cache = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  cache <- lik_cache_resolve(cache, prep)
  compiled <- .expr_lookup_compiled(expr, prep)
  if (!is.null(compiled) && is.function(compiled$scenario_fn)) {
    scenarios <- compiled$scenario_fn(t, component, forced_complete, forced_survive)
    total <- if (length(scenarios) == 0) 0.0 else sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total, "scenarios") <- scenarios %||% list()
    return(total)
  }
  scenarios <- .guard_scenarios_slow(
    expr, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    cache = cache
  )
  total <- if (length(scenarios) == 0) 0.0 else sum(vapply(scenarios, `[[`, numeric(1), "weight"))
  attr(total, "scenarios") <- scenarios %||% list()
  total
}

.make_guard_scenario_fn <- function(expr, prep) {
  reference <- expr[["reference"]]
  ref_node <- .expr_lookup_compiled(reference, prep)
  if (is.null(ref_node) || !is.function(ref_node$scenario_fn)) return(NULL)
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  block_node <- if (!is.null(blocker)) .expr_lookup_compiled(blocker, prep) else NULL
  unless_nodes <- lapply(unless_list, function(unl) .expr_lookup_compiled(unl, prep))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    cache <- lik_cache_create()
    ref_scen <- ref_node$scenario_fn(t, component, forced_complete, forced_survive)
    if (length(ref_scen) == 0) return(list())
    blocker_sources <- integer(0)
    if (!is.null(blocker)) {
      blocker_sources <- .expr_sources(blocker, prep)
      blocker_label <- .label_to_id(prep, blocker[["source"]])
      if (!is.na(blocker_label)) {
        blocker_sources <- .forced_union(prep, blocker_sources, blocker_label)
      }
    }
    out <- list()
    for (sc in ref_scen) {
      if (is.null(sc)) next
      fc_ctx <- .coerce_forced_ids(prep, sc$forced_complete)
      fs_ctx <- .coerce_forced_ids(prep, sc$forced_survive)
      S_eff <- .guard_effective_survival(
        expr, t, prep, component,
        fc_ctx, fs_ctx, cache,
        block_node = block_node,
        unless_nodes = unless_nodes
      )
      if (!is.finite(S_eff) || S_eff <= 0) next
      weight <- sc$weight * S_eff
      if (!is.finite(weight) || weight <= 0) next
      fcomp <- fc_ctx
      fsurv <- .forced_union(prep, fs_ctx, blocker_sources)
      out[[length(out) + 1L]] <- list(
        weight = as.numeric(weight),
        forced_complete = fcomp,
        forced_survive = .coerce_forced_ids(prep, fsurv)
      )
    }
    out
  }
}

.event_survival_at <- function(prep, id, component, t,
                               forced_complete = integer(0),
                               forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(1.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    return(.pool_survival(prep, id, component, t,
                          forced_complete = forced_complete,
                          forced_survive = forced_survive))
  }
  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) return(.acc_survival(acc_def, t))
  0.0
}

.event_cdf_at <- function(prep, id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  1.0 - .event_survival_at(prep, id, component, t,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive)
}

.event_density_at <- function(prep, id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(0.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    return(.pool_density(prep, id, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive))
  }
  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) return(.acc_density(acc_def, t))
  0.0
}

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  if (!is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  k <- as.integer(pool_def[["k"]] %||% 1L)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L || k < 1L || k > n) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  member_ids <- .labels_to_ids(prep, members)
  pool_idx <- .label_to_id(prep, pool_id)
  
  idx_seq <- seq_len(n)
  dens_vec <- numeric(n)
  cdf_vec <- numeric(n)
  surv_vec <- numeric(n)
  for (i in idx_seq) {
    mid <- members[[i]]
    dens_vec[[i]] <- .event_density_at(
      prep, mid, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    cdf_vec[[i]] <- .event_cdf_at(
      prep, mid, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    surv_vec[[i]] <- 1.0 - cdf_vec[[i]]
  }
  
  scenarios <- list()
  add_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return()
    scenarios[[length(scenarios) + 1L]] <<- list(
      weight = as.numeric(weight),
      forced_complete = .coerce_forced_ids(prep, fcomp),
      forced_survive = .coerce_forced_ids(prep, fsurv)
    )
  }
  
  fast_mode <- (length(forced_complete) == 0L && length(forced_survive) == 0L)
  fast_done <- FALSE
  if (fast_mode) {
    if (k == 1L) {
      if (n == 1L) {
        weight <- dens_vec[[1]]
        if (is.finite(weight) && weight > 0) {
          finisher_ids <- member_ids[1]
          if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
          finisher_ids <- finisher_ids[!is.na(finisher_ids)]
          add_scenario(weight, finisher_ids, integer(0))
          fast_done <- TRUE
        }
      } else {
        for (idx in idx_seq) {
          dens_mid <- dens_vec[[idx]]
          if (!is.finite(dens_mid) || dens_mid <= 0) next
          survivors_idx <- idx_seq[idx_seq != idx]
          weight <- dens_mid
          if (length(survivors_idx) > 0L) {
            weight <- weight * prod(surv_vec[survivors_idx])
          }
          if (!is.finite(weight) || weight <= 0) next
          finisher_ids <- member_ids[idx]
          if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
          finisher_ids <- finisher_ids[!is.na(finisher_ids)]
          survivor_ids <- if (length(survivors_idx) > 0L) member_ids[survivors_idx] else integer(0)
          add_scenario(weight, finisher_ids, survivor_ids)
          fast_done <- TRUE
        }
      }
    } else if (k == n) {
      for (idx in idx_seq) {
        dens_mid <- dens_vec[[idx]]
        if (!is.finite(dens_mid) || dens_mid <= 0) next
        others_idx <- idx_seq[idx_seq != idx]
        weight <- dens_mid
        if (length(others_idx) > 0L) {
          weight <- weight * prod(cdf_vec[others_idx])
        }
        if (!is.finite(weight) || weight <= 0) next
        finisher_ids <- member_ids
        if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
        finisher_ids <- finisher_ids[!is.na(finisher_ids)]
        add_scenario(weight, finisher_ids, integer(0))
        fast_done <- TRUE
      }
    } else {
      template_info <- .build_pool_templates(pool_id, members, member_ids, pool_idx, k)
      if (!is.null(template_info) && length(template_info$templates) > 0L) {
        templates <- template_info$templates
        finisher_map <- template_info$finisher_map
        for (idx in idx_seq) {
          dens_mid <- dens_vec[[idx]]
          if (!is.finite(dens_mid) || dens_mid <= 0) next
          template_ids <- finisher_map[[idx]]
          if (length(template_ids) == 0L) next
          for (tid in template_ids) {
            tmpl <- templates[[tid]]
            weight <- dens_mid
            if (length(tmpl$complete_idx) > 0L) {
              weight <- weight * prod(cdf_vec[tmpl$complete_idx])
            }
            if (length(tmpl$survivor_idx) > 0L) {
              weight <- weight * prod(surv_vec[tmpl$survivor_idx])
            }
            if (!is.finite(weight) || weight <= 0) next
            add_scenario(weight, tmpl$forced_complete_ids, tmpl$forced_survive_ids)
            fast_done <- TRUE
          }
        }
      }
    }
  }
  
  if (!fast_mode || !fast_done) {
    others_idx_list <- lapply(idx_seq, function(i) idx_seq[idx_seq != i])
    forced_complete_flags <- if (length(forced_complete) > 0L) member_ids %in% forced_complete else rep(FALSE, n)
    forced_survive_flags <- if (length(forced_survive) > 0L) member_ids %in% forced_survive else rep(FALSE, n)
    free_mask_base <- !(forced_complete_flags | forced_survive_flags)
    
    for (idx in idx_seq) {
      if (forced_survive_flags[[idx]]) next
      dens_mid <- dens_vec[[idx]]
      if (!is.finite(dens_mid) || dens_mid <= 0) next
      
      need <- k - 1L
      others_idx <- others_idx_list[[idx]]
      if (need > length(others_idx)) next
      
      fc_idx <- if (length(others_idx) > 0L) others_idx[forced_complete_flags[others_idx]] else integer(0)
      fs_idx <- if (length(others_idx) > 0L) others_idx[forced_survive_flags[others_idx]] else integer(0)
      forced_count <- length(fc_idx)
      if (forced_count > need) next
      remaining_need <- need - forced_count
      
      free_idx <- if (length(others_idx) > 0L) others_idx[free_mask_base[others_idx]] else integer(0)
      if (remaining_need > length(free_idx)) next
      
      base_factor <- dens_mid
      if (length(fc_idx) > 0L) base_factor <- base_factor * prod(cdf_vec[fc_idx])
      if (length(fs_idx) > 0L) base_factor <- base_factor * prod(surv_vec[fs_idx])
      if (!is.finite(base_factor) || base_factor <= 0) next
      
      finisher_ids <- member_ids[idx]
      if (length(fc_idx) > 0L) finisher_ids <- c(finisher_ids, member_ids[fc_idx])
      if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
      forced_complete_ids <- .forced_union(prep, forced_complete, finisher_ids)
      base_survive_ids <- .forced_union(prep, forced_survive, member_ids[fs_idx])
      
      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
          )
        }
        next
      }
      
      combos <- if (remaining_need == length(free_idx)) {
        list(free_idx)
      } else if (remaining_need == 1L) {
        lapply(free_idx, function(i) i)
      } else if (remaining_need <= 0L) {
        list(integer(0))
      } else {
        utils::combn(free_idx, remaining_need, simplify = FALSE)
      }
      
      if (length(combos) == 0L) next
      
      for (combo_idx in combos) {
        if (length(combo_idx) != remaining_need) next
        
        combo_weight <- base_factor
        if (length(combo_idx) > 0L) combo_weight <- combo_weight * prod(cdf_vec[combo_idx])
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        if (length(free_idx) > 0L) {
          drop_pos <- match(combo_idx, free_idx, nomatch = 0L)
          drop_pos <- drop_pos[drop_pos > 0L]
          survivors_idx <- if (length(drop_pos) > 0L) free_idx[-drop_pos] else free_idx
        } else {
          survivors_idx <- integer(0)
        }
        
        if (length(survivors_idx) > 0L) {
          combo_weight <- combo_weight * prod(surv_vec[survivors_idx])
        }
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        add_scenario(
          combo_weight,
          .forced_union(prep, forced_complete_ids, member_ids[combo_idx]),
          .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
        )
      }
    }
  }
  
  if (!fast_done || !fast_mode) {
    others_idx_list <- lapply(idx_seq, function(i) idx_seq[idx_seq != i])
    forced_complete_flags <- if (length(forced_complete) > 0L) member_ids %in% forced_complete else rep(FALSE, n)
    forced_survive_flags <- if (length(forced_survive) > 0L) member_ids %in% forced_survive else rep(FALSE, n)
    free_mask_base <- !(forced_complete_flags | forced_survive_flags)
    
    for (idx in idx_seq) {
      if (forced_survive_flags[[idx]]) next
      dens_mid <- dens_vec[[idx]]
      if (!is.finite(dens_mid) || dens_mid <= 0) next
      
      need <- k - 1L
      others_idx <- others_idx_list[[idx]]
      if (need > length(others_idx)) next
      
      fc_idx <- if (length(others_idx) > 0L) others_idx[forced_complete_flags[others_idx]] else integer(0)
      fs_idx <- if (length(others_idx) > 0L) others_idx[forced_survive_flags[others_idx]] else integer(0)
      forced_count <- length(fc_idx)
      if (forced_count > need) next
      remaining_need <- need - forced_count
      
      free_idx <- if (length(others_idx) > 0L) others_idx[free_mask_base[others_idx]] else integer(0)
      if (remaining_need > length(free_idx)) next
      
      base_factor <- dens_mid
      if (length(fc_idx) > 0L) base_factor <- base_factor * prod(cdf_vec[fc_idx])
      if (length(fs_idx) > 0L) base_factor <- base_factor * prod(surv_vec[fs_idx])
      if (!is.finite(base_factor) || base_factor <= 0) next
      
      finisher_ids <- member_ids[idx]
      if (length(fc_idx) > 0L) finisher_ids <- c(finisher_ids, member_ids[fc_idx])
      if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
      forced_complete_ids <- .forced_union(prep, forced_complete, finisher_ids)
      base_survive_ids <- .forced_union(prep, forced_survive, member_ids[fs_idx])
      
      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
          )
        }
        next
      }
      
      combos <- if (remaining_need == length(free_idx)) {
        list(free_idx)
      } else if (remaining_need == 1L) {
        lapply(free_idx, function(i) i)
      } else if (remaining_need <= 0L) {
        list(integer(0))
      } else {
        utils::combn(free_idx, remaining_need, simplify = FALSE)
      }
      
      if (length(combos) == 0L) next
      
      for (combo_idx in combos) {
        if (length(combo_idx) != remaining_need) next
        
        combo_weight <- base_factor
        if (length(combo_idx) > 0L) combo_weight <- combo_weight * prod(cdf_vec[combo_idx])
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        if (length(free_idx) > 0L) {
          drop_pos <- match(combo_idx, free_idx, nomatch = 0L)
          drop_pos <- drop_pos[drop_pos > 0L]
          survivors_idx <- if (length(drop_pos) > 0L) free_idx[-drop_pos] else free_idx
        } else {
          survivors_idx <- integer(0)
        }
        
        if (length(survivors_idx) > 0L) {
          combo_weight <- combo_weight * prod(surv_vec[survivors_idx])
        }
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        add_scenario(
          combo_weight,
          .forced_union(prep, forced_complete_ids, member_ids[combo_idx]),
          .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
        )
      }
    }
  }
  
  total_density <- if (length(scenarios) == 0L) 0.0 else sum(vapply(scenarios, `[[`, numeric(1), "weight"))
  if (length(scenarios) > 0L) {
    attr(total_density, "scenarios") <- scenarios
  } else {
    attr(total_density, "scenarios") <- list()
  }
  total_density
}

.pool_survival <- function(prep, pool_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)) {
  if (!is.finite(t)) return(0.0)
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  k <- as.integer(pool_def[["k"]] %||% 1L)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(1.0)
  if (k > n) return(1.0)
  if (k < 1L) return(0.0)
  
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  
  Fvec <- numeric(n)
  for (i in seq_len(n)) {
    id <- members[[i]]
    Fvec[[i]] <- .event_cdf_at(prep, id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive)
  }
  Svec <- 1.0 - Fvec
  coeffs <- pool_coeffs(Svec, Fvec)
  upto <- min(length(coeffs), k)
  sum(coeffs[seq_len(upto)])
}

# ==============================================================================
# PART 4: Expression Evaluation (Likelihood)
# ==============================================================================

# Scenario helpers for generic expressions at time t.
# Returns a list of scenario records:
#   list(weight = numeric(1),
#        forced_complete = integer(),
#        forced_survive  = integer())
.expr_scenarios_at <- function(expr, t, prep, component,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               cache = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  cache <- lik_cache_resolve(cache, prep)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- lik_cache_key("scen", node_id, component, t,
                             forced_complete, forced_survive)
  cached <- lik_cache_get(cache, cache_key)
  if (!is.null(cached)) return(cached)
  if (!is.null(compiled) && is.function(compiled$scenario_fn)) {
    res <- compiled$scenario_fn(t, component, forced_complete, forced_survive)
  } else {
    res <- .expr_scenarios_at_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) lik_cache_set(cache, cache_key, res) else res
}

.expr_scenarios_at_slow <- function(expr, t, prep, component,
                                    forced_complete = integer(0),
                                    forced_survive = integer(0),
                                    cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  
  make_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return(NULL)
    list(
      weight = as.numeric(weight),
      forced_complete = .coerce_forced_ids(prep, fcomp),
      forced_survive = .coerce_forced_ids(prep, fsurv)
    )
  }
  
  if (identical(kind, "event")) {
    source_id <- expr[["source"]]
    source_idx <- .label_to_id(prep, source_id)
    source_ids <- if (is.na(source_idx)) integer(0) else source_idx
    weight <- .event_density_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    if (!is.finite(weight) || weight <= 0) return(list())
    pool_scenarios <- attr(weight, "scenarios")
    if (!is.null(pool_scenarios) && length(pool_scenarios) > 0) {
      out <- list()
      for (psc in pool_scenarios) {
        if (is.null(psc)) next
        w <- psc$weight
        if (!is.finite(w) || w <= 0) next
        fcomp <- .forced_union(prep, forced_complete, c(source_ids, psc$forced_complete))
        fsurv <- .forced_union(prep, forced_survive, psc$forced_survive)
        sc <- make_scenario(w, fcomp, fsurv)
        if (!is.null(sc)) out[[length(out) + 1L]] <- sc
      }
      return(out)
    }
    sc <- make_scenario(
      weight,
      fcomp = .forced_union(prep, forced_complete, source_ids),
      fsurv = forced_survive
    )
    return(if (is.null(sc)) list() else list(sc))
  }
  
  if (identical(kind, "and")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(list())
    out <- list()
    for (i in seq_along(args)) {
      si <- .expr_scenarios_at(
        args[[i]], t, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        cache = cache)
      if (length(si) == 0) next
      others <- args[-i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        ok <- TRUE
        for (aj in others) {
          Fj <- .eval_expr_cdf_cond(
            aj, t, prep, component,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            cache = cache)
          if (!is.finite(Fj) || Fj <= 0) {
            ok <- FALSE
            break
          }
          weight <- weight * Fj
          if (!is.finite(weight) || weight <= 0) {
            ok <- FALSE
            break
          }
          fcomp <- .forced_union(prep, fcomp, .expr_sources(aj, prep))
        }
        if (!ok) next
        sc_out <- make_scenario(weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    return(out)
  }
  
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(list())
    out <- list()
    for (i in seq_along(args)) {
      si <- .expr_scenarios_at(
        args[[i]], t, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        cache = cache)
      if (length(si) == 0) next
      others <- args[-i]
      others_req <- lapply(others, function(aj) .expr_sources(aj, prep))
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        valid <- TRUE
        for (idx in seq_along(others)) {
          req_j <- others_req[[idx]]
          if (length(req_j) == 0L) next
          if (all(req_j %in% fcomp)) {
            valid <- FALSE
            break
          }
          witness <- req_j[!(req_j %in% fcomp)]
          if (length(witness) >= 1L) {
            fsurv <- .forced_union(prep, fsurv, witness[[1]])
          } else {
            valid <- FALSE
            break
          }
        }
        if (!valid) next
        sc_out <- make_scenario(weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    return(out)
  }
  
  if (identical(kind, "guard")) {
    return(.guard_scenarios_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    ))
  }
  
  weight <- .eval_expr_likelihood(
    expr, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    cache = cache
  )
  if (!is.finite(weight) || weight <= 0) return(list())
  sc <- make_scenario(
    weight,
    fcomp = .forced_union(prep, forced_complete, .expr_sources(expr, prep)),
    fsurv = forced_survive
  )
  if (is.null(sc)) list() else list(sc)
}

.eval_expr_cdf_cond <- function(expr, t, prep, component,
                                forced_complete = integer(0),
                                forced_survive = integer(0),
                                cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  cache <- lik_cache_resolve(cache, prep)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- lik_cache_key("cdf", node_id, component, t,
                             forced_complete, forced_survive)
  cached <- lik_cache_get(cache, cache_key)
  if (!is.null(cached)) return(cached)
  res <- if (!is.null(compiled) && is.function(compiled$cdf_fn)) {
    compiled$cdf_fn(t, component, forced_complete, forced_survive)
  } else {
    .eval_expr_cdf_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) lik_cache_set(cache, cache_key, res) else res
}

.eval_expr_cdf_cond_slow <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    return(.event_cdf_at(prep, src, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive))
  }
  if (identical(kind, "and")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(0.0)
    prod_val <- 1.0
    for (a in args) {
      Fa <- .eval_expr_cdf_cond(a, t, prep, component,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                cache = cache)
      prod_val <- prod_val * Fa
      if (prod_val == 0) break
    }
    return(prod_val)
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(0.0)
    prod_surv <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache)
      prod_surv <- prod_surv * Sa
      if (prod_surv == 0) break
    }
    return(1.0 - prod_surv)
  }
  if (identical(kind, "guard")) {
    if (t <= 0) return(0.0)
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        .eval_expr_likelihood(
          expr, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = cache
        )
      }, numeric(1))
    }
    val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                     rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                     stop.on.error = FALSE)$value,
                    error = function(e) 0.0)
    out <- as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    return(max(0.0, min(1.0, out)))
  }
  return(.eval_expr_cdf(expr, t, prep, component))
}

.eval_expr_survival_cond <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  cache <- lik_cache_resolve(cache, prep)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- lik_cache_key("surv", node_id, component, t,
                             forced_complete, forced_survive)
  cached <- lik_cache_get(cache, cache_key)
  if (!is.null(cached)) return(cached)
  res <- if (!is.null(compiled) && is.function(compiled$surv_fn)) {
    compiled$surv_fn(t, component, forced_complete, forced_survive)
  } else {
    .eval_expr_survival_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) lik_cache_set(cache, cache_key, res) else res
}

.eval_expr_survival_cond_slow <- function(expr, t, prep, component,
                                          forced_complete = integer(0),
                                          forced_survive = integer(0),
                                          cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    return(.event_survival_at(prep, src, component, t,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive))
  }
  if (identical(kind, "and")) {
    return(1.0 - .eval_expr_cdf_cond(expr, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache))
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(1.0)
    prod_val <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache)
      prod_val <- prod_val * Sa
      if (prod_val == 0) break
    }
    return(prod_val)
  }
  if (identical(kind, "guard")) {
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .eval_expr_likelihood(
          expr, uu, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = cache
        )
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                stop.on.error = FALSE)[["value"]],
               error = function(e) 0.0)
    }
    val <- compute_integral()
    out <- 1.0 - as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    return(max(0.0, min(1.0, out)))
  }
  return(.eval_expr_survival(expr, t, prep, component))
}

.scenario_density_with_competitors <- function(expr, t, prep, component,
                                               competitor_exprs = list(),
                                               cache = NULL) {
  if (is.na(t) || !is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  cache <- lik_cache_resolve(cache, prep)
  scenarios <- .expr_scenarios_at(
    expr, t, prep, component,
    forced_complete = integer(0),
    forced_survive = integer(0),
    cache = cache)
  if (length(scenarios) == 0) {
    return(0.0)
  }
  total <- 0.0
  for (sc in scenarios) {
    if (is.null(sc) || sc$weight <= 0) next
    weight <- sc$weight
    if (length(competitor_exprs) > 0) {
      surv_prod <- .compute_survival_product(
        outcome_expr = expr,
        competitor_exprs = competitor_exprs,
        prep = prep,
        component = component,
        t = t,
        cache = cache
      )
      weight <- weight * surv_prod
    }
    total <- total + weight
  }
  if (length(scenarios) > 0) attr(total, "scenarios") <- scenarios
  total
}

# ============================================================================== 
# Fast-path helpers for competitor survival products (ported from legacy API)
# ==============================================================================

.cluster_competitors <- function(exprs, prep) {
  n <- length(exprs)
  if (n == 0) return(list())
  src_sets <- lapply(exprs, .expr_sources, prep)
  visited <- rep(FALSE, n)
  clusters <- list()
  for (i in seq_len(n)) {
    if (visited[[i]]) next
    stack <- i
    visited[[i]] <- TRUE
    members <- integer(0)
    while (length(stack) > 0) {
      idx <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      members <- c(members, idx)
      for (j in seq_len(n)) {
        if (visited[[j]]) next
        if (length(intersect(src_sets[[idx]], src_sets[[j]])) > 0) {
          visited[[j]] <- TRUE
          stack <- c(stack, j)
        }
      }
    }
    clusters[[length(clusters) + 1L]] <- list(indices = members, exprs = exprs[members], sources = src_sets[members])
  }
  clusters
}

.get_event_source <- function(node) {
  if (is.null(node)) return(NULL)
  if (!is.null(node$kind) && identical(node$kind, "event")) {
    return(node$source %||% NULL)
  }
  NULL
}

.extract_guard_plus <- function(expr) {
  if (is.null(expr) || is.null(expr$kind)) return(NULL)
  if (identical(expr$kind, "guard")) {
    ref_id <- .get_event_source(expr$reference)
    blk_id <- .get_event_source(expr$blocker)
    if (is.null(ref_id) || is.null(blk_id)) return(NULL)
    return(list(guard = expr, ref = ref_id, block = blk_id, extras = character(0)))
  }
  if (identical(expr$kind, "and")) {
    args <- expr$args %||% list()
    guard_idx <- which(vapply(args, function(arg) !is.null(arg$kind) && identical(arg$kind, "guard"), logical(1)))
    if (length(guard_idx) != 1) return(NULL)
    guard <- args[[guard_idx[[1]]]]
    ref_id <- .get_event_source(guard$reference)
    blk_id <- .get_event_source(guard$blocker)
    if (is.null(ref_id) || is.null(blk_id)) return(NULL)
    extras <- args[-guard_idx[[1]]]
    extra_sources <- character(0)
    for (extra in extras) {
      src <- .get_event_source(extra)
      if (is.null(src)) return(NULL)
      extra_sources <- c(extra_sources, src)
    }
    return(list(guard = guard, ref = ref_id, block = blk_id, extras = unique(extra_sources)))
  }
  NULL
}

.event_density_fn <- function(source_id, prep, component) {
  event_expr <- list(kind = "event", source = source_id, k = NULL)
  function(t) vapply(
    t,
    function(tt) .eval_expr_likelihood(event_expr, tt, prep, component),
    numeric(1)
  )
}

.event_survival_fn <- function(source_id, prep, component) {
  event_expr <- list(kind = "event", source = source_id, k = NULL)
  function(t) vapply(
    t,
    function(tt) .eval_expr_survival(event_expr, tt, prep, component),
    numeric(1)
  )
}

.compute_survival_product <- function(outcome_expr, competitor_exprs, prep, component, t,
                                      cache = NULL) {
  if (length(competitor_exprs) == 0) return(1.0)
  clusters <- .cluster_competitors(competitor_exprs, prep)
  prod_val <- 1.0
  for (cluster in clusters) {
    exprs <- cluster$exprs
    cluster_val <- prod(vapply(exprs, function(ce) {
      .eval_expr_survival_cond(
        ce, t, prep, component,
        forced_complete = integer(0),
        forced_survive = integer(0),
        cache = cache
      )
    }, numeric(1)))
    prod_val <- prod_val * cluster_val
  }
  prod_val
}

# Get CDF (cumulative probability) for an expression
.eval_expr_cdf <- function(expr, t, prep, component) {
  kind <- expr[['kind']]
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(if (t >= 0) 1.0 else 0.0)
    }
    
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_cdf_at(prep, source_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)))
    }
    
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = instant completion
      }
      return(.event_cdf_at(prep, source_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)))
    }
    
    return(1.0)
  }
  
  if (identical(kind, "and")) {
    # CDF of max: product of CDFs
    args <- expr[['args']]
    prod_cdf <- 1.0
    for (arg in args) {
      cdf_arg <- .eval_expr_cdf(arg, t, prep, component)
      prod_cdf <- prod_cdf * cdf_arg
    }
    return(prod_cdf)
  }
  
  if (identical(kind, "or")) {
    # CDF of min: 1 - product of survivals
    args <- expr[['args']]
    prod_surv <- 1.0
    for (arg in args) {
      surv_arg <- .eval_expr_survival(arg, t, prep, component)
      prod_surv <- prod_surv * surv_arg
    }
    return(1 - prod_surv)
  }
  
  if (identical(kind, "guard")) {
    # Guard: integrate density to obtain CDF
    if (t <= 0) return(0.0)
    if (!is.finite(t)) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        .eval_expr_likelihood(expr, ui, prep, component)
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(
        stats::integrate(dens_fun, lower = 0, upper = t,
                         rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                         stop.on.error = FALSE)$value,
        error = function(e) 0.0
      )
    }
    val <- compute_integral()
    cdf <- as.numeric(val)
    if (!is.finite(cdf)) cdf <- 0.0
    if (cdf < 0) cdf <- 0.0
    if (cdf > 1) cdf <- 1.0
    return(cdf)
  }
  
  return(NA_real_)
}

# Get survival function for an expression
.eval_expr_survival <- function(expr, t, prep, component) {
  kind <- expr[['kind']]
  
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(if (t < 0) 1.0 else 0.0)
    }
    
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_survival_at(prep, source_id, component, t,
                                forced_complete = integer(0),
                                forced_survive = integer(0)))
    }
    
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = always survive (never fires)
      }
      return(.event_survival_at(prep, source_id, component, t,
                                forced_complete = integer(0),
                                forced_survive = integer(0)))
    }
    
    return(1.0)
  }
  
  if (identical(kind, "and")) {
    # Survival of max: at least one must be > t
    # S_max(t) = 1 - F_max(t) = 1 - Π F_i(t)
    return(1 - .eval_expr_cdf(expr, t, prep, component))
  }
  
  if (identical(kind, "or")) {
    # Survival of min: all must be > t
    # S_min(t) = Π S_i(t)
    args <- expr[['args']]
    prod_surv <- 1.0
    for (arg in args) {
      surv_arg <- .eval_expr_survival(arg, t, prep, component)
      prod_surv <- prod_surv * surv_arg
    }
    return(prod_surv)
  }
  
  if (identical(kind, "guard")) {
    # Survival of guard = 1 - ∫_0^t f_guard(u) du
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .eval_expr_likelihood(expr, uu, prep, component)
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                stop.on.error = FALSE)[['value']],
               error = function(e) 0.0)
    }
    val <- compute_integral()
    s <- 1.0 - as.numeric(val)
    if (!is.finite(s)) s <- 0.0
    if (s < 0) s <- 0.0
    if (s > 1) s <- 1.0
    return(s)
  }
  
  return(0.0)
}

# Compute likelihood (density or probability) for an expression
.eval_expr_likelihood <- function(expr, t, prep, component,
                                  forced_complete = integer(0),
                                  forced_survive = integer(0),
                                  cache = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  kind <- expr[['kind']]
  
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    # Handle special sources
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      # These don't contribute density directly
      return(0.0)
    }
    
    # Check if it's a pool
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_density_at(prep, source_id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive))
    }
    
    # Check if it's an accumulator
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      # Check component membership
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(0.0)  # Not active in this component
      }
      return(.event_density_at(prep, source_id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive))
    }
    
    return(0.0)
  }
  
  if (identical(kind, "and")) {
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    scenarios <- .expr_scenarios_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache)
    if (length(scenarios) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    total_density <- sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total_density, "scenarios") <- scenarios
    return(total_density)
  }
  
  if (identical(kind, "or")) {
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    scenarios <- .expr_scenarios_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache)
    if (length(scenarios) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    total_density <- sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total_density, "scenarios") <- scenarios
    return(total_density)
  }
  
  if (identical(kind, "not")) {
    # NOT typically doesn't have density, only affects probability
    # If we reach here, it means NOT is used in an outcome expression
    # This is unusual but could mean "outcome never happens"
    return(0.0)
  }
  
  if (identical(kind, "guard")) {
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    val <- .guard_density_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
    return(val)
  }
  
  stop(sprintf("Unsupported expression kind '%s'", kind))
}

# ==============================================================================
