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
    if (is.null(attr(expr, ".lik_signature", exact = TRUE))) {
      attr(expr, ".lik_signature") <- sig
    }
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
    node <- .compiled_node_attach_ops(node)
    node_id <- length(nodes_env$nodes) + 1L
    node$id <- node_id
    node$signature <- sig
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

.node_eval_entry <- function(state, node, component, t,
                             forced_complete = integer(0),
                             forced_survive = integer(0),
                             extra = NULL,
                             init = TRUE) {
  if (is.null(state) || is.null(node)) return(NULL)
  .eval_state_entry(
    state,
    node$id,
    component,
    t,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    extra = extra,
    init = init
  )
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

.compiled_nodes_table <- function(prep) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) return(list())
  comp$nodes %||% list()
}

.compiled_node_fetch <- function(prep, node_id) {
  if (is.null(node_id) || is.na(node_id)) return(NULL)
  nodes <- .compiled_nodes_table(prep)
  nodes[[as.integer(node_id)]]
}

.compiled_node_attach_ops <- function(node) {
  if (!is.list(node$ops)) node$ops <- list()
  node$ops$density <- node$density_fn
  node$ops$density_fast <- node$density_fast_fn
  node$ops$survival <- node$surv_fn
  node$ops$survival_fast <- node$surv_fast_fn
  node$ops$cdf <- node$cdf_fn
  node$ops$cdf_fast <- node$cdf_fast_fn
  node$ops$scenario <- node$scenario_fn
  node
}

.node_ops_get <- function(node, key, fallback = NULL) {
  if (is.null(node)) return(fallback)
  ops <- node$ops
  fn <- if (is.list(ops)) ops[[key]] else NULL
  if (is.function(fn)) return(fn)
  fallback
}

.fast_survival_product <- function(exprs, prep, component, t) {
  n <- length(exprs)
  if (n == 0) return(1.0)
  prod_val <- 1.0
  for (ex in exprs) {
    node <- .expr_lookup_compiled(ex, prep)
    if (is.null(node)) return(NULL)
    surv_fast <- .node_ops_get(node, "survival_fast", node$surv_fast_fn)
    if (!is.function(surv_fast)) return(NULL)
    val <- surv_fast(t, component)
    surv <- if (length(val) == 0L) 1.0 else as.numeric(val[[1]])
    if (!is.finite(surv)) return(NULL)
    if (surv < 0) surv <- 0.0 else if (surv > 1) surv <- 1.0
    prod_val <- prod_val * surv
    if (!is.finite(prod_val)) return(NULL)
    if (prod_val == 0) break
  }
  prod_val
}

.build_compiled_node <- function(expr, child_info, prep, nodes_env) {
  kind <- expr[["kind"]]
  node <- list(
    kind = kind,
    sources = .expr_sources(expr, prep),
    expr = expr,
    needs_forced = FALSE,
    scenario_sensitive = FALSE,
    density_fast_fn = NULL,
    surv_fast_fn = NULL,
    cdf_fast_fn = NULL
  )
  if (identical(kind, "event")) {
    node$cdf_fn <- .compile_once(.make_event_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_event_survival_fn(expr, prep))
    node$density_fn <- .compile_once(.make_event_density_fn(expr, prep))
    node$scenario_fn <- .compile_once(.make_event_scenario_fn(expr, prep))
    node$needs_forced <- TRUE
    source_id <- expr[["source"]]
    node$scenario_sensitive <- !is.null(prep[["pools"]][[source_id]])
    fast_cdf <- .make_event_cdf_fast_fn(expr, prep)
    fast_surv <- .make_event_survival_fast_fn(expr, prep)
    fast_dens <- .make_event_density_fast_fn(expr, prep)
    node$cdf_fast_fn <- if (is.function(fast_cdf)) .compile_once(fast_cdf) else NULL
    node$surv_fast_fn <- if (is.function(fast_surv)) .compile_once(fast_surv) else NULL
    node$density_fast_fn <- if (is.function(fast_dens)) .compile_once(fast_dens) else NULL
    return(.compiled_node_attach_ops(node))
  }
  if (kind %in% c("and", "or")) {
    child_ids <- as.integer(child_info$args %||% integer(0))
    node$args <- child_ids
    if (length(child_ids) == 0L || any(is.na(child_ids))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$scenario_fn <- NULL
      return(.compiled_node_attach_ops(node))
    }
    child_nodes <- lapply(child_ids, function(id) nodes_env$nodes[[id]])
    if (any(vapply(child_nodes, is.null, logical(1)))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$density_fn <- NULL
      node$scenario_fn <- NULL
      return(.compiled_node_attach_ops(node))
    }
    node$needs_forced <- any(vapply(child_nodes, function(n) isTRUE(n$needs_forced), logical(1)))
    node$scenario_sensitive <- any(vapply(child_nodes, function(n) isTRUE(n$scenario_sensitive), logical(1)))
    if (identical(kind, "and")) {
      have_child_cdf <- all(vapply(child_nodes, function(n) is.function(n$cdf_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      have_child_density <- all(vapply(child_nodes, function(n) is.function(n$density_fn), logical(1)))
      node$cdf_fn <- if (have_child_cdf) .compile_once(.make_and_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (is.function(node$cdf_fn)) .compile_once(.make_and_survival_fn(node$cdf_fn)) else NULL
      node$density_fn <- if (have_child_density && have_child_cdf) .compile_once(.make_and_density_fn(child_nodes, prep)) else NULL
      node$scenario_fn <- if (have_child_cdf && have_child_scen) .compile_once(.make_and_scenario_fn(child_nodes, prep)) else NULL
      if (!node$scenario_sensitive) {
        raw_cdf_fast <- .make_and_cdf_fast_fn(child_nodes, prep)
        node$cdf_fast_fn <- .compile_once(raw_cdf_fast)
        node$surv_fast_fn <- if (is.function(node$cdf_fast_fn)) .compile_once(.make_and_survival_fast_fn(node$cdf_fast_fn)) else NULL
        node$density_fast_fn <- .compile_once(.make_and_density_fast_fn(child_nodes, prep))
      }
    } else {
      have_child_surv <- all(vapply(child_nodes, function(n) is.function(n$surv_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      have_child_density <- all(vapply(child_nodes, function(n) is.function(n$density_fn), logical(1)))
      node$cdf_fn <- if (have_child_surv) .compile_once(.make_or_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (have_child_surv) .compile_once(.make_or_survival_fn(child_nodes, prep)) else NULL
      node$density_fn <- if (have_child_density && have_child_surv) .compile_once(.make_or_density_fn(child_nodes, prep)) else NULL
      node$scenario_fn <- if (have_child_scen) .compile_once(.make_or_scenario_fn(child_nodes, prep)) else NULL
      if (!node$scenario_sensitive) {
        node$surv_fast_fn <- .compile_once(.make_or_survival_fast_fn(child_nodes, prep))
        node$cdf_fast_fn <- if (is.function(node$surv_fast_fn)) .compile_once(.make_or_cdf_fast_fn(node$surv_fast_fn)) else NULL
        node$density_fast_fn <- .compile_once(.make_or_density_fast_fn(child_nodes, prep))
      }
    }
    return(.compiled_node_attach_ops(node))
  }
  if (identical(kind, "guard")) {
    node$reference_id <- child_info$reference %||% NA_integer_
    node$blocker_id <- child_info$blocker %||% NA_integer_
    node$unless_ids <- as.integer(child_info$unless %||% integer(0))
    node$cdf_fn <- .compile_once(.make_guard_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_guard_survival_fn(expr, prep))
    guard_scen <- .make_guard_scenario_fn(expr, prep)
    node$scenario_fn <- if (is.function(guard_scen)) .compile_once(guard_scen) else NULL
    node$density_fn <- .compile_once(.make_guard_density_fn(expr, prep))
    node$needs_forced <- TRUE
    node$scenario_sensitive <- TRUE
    ref_node_compiled <- if (!is.na(node$reference_id)) nodes_env$nodes[[as.integer(node$reference_id)]] else NULL
    block_node_compiled <- if (!is.na(node$blocker_id)) nodes_env$nodes[[as.integer(node$blocker_id)]] else NULL
    unless_nodes_compiled <- if (length(node$unless_ids) > 0L) {
      lapply(node$unless_ids, function(uid) nodes_env$nodes[[as.integer(uid)]])
    } else {
      list()
    }
    fast_guard <- .make_guard_fast_fns(
      expr, prep,
      ref_node_override = ref_node_compiled,
      block_node_override = block_node_compiled,
      unless_nodes_override = unless_nodes_compiled
    )
    node$density_fast_fn <- fast_guard$density
    node$surv_fast_fn <- fast_guard$survival
    node$cdf_fast_fn <- fast_guard$cdf
    return(.compiled_node_attach_ops(node))
  }
  node$cdf_fn <- NULL
  node$surv_fn <- NULL
  node$density_fn <- NULL
  node$scenario_fn <- NULL
  .compiled_node_attach_ops(node)
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

.make_event_density_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    .event_density_at(
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

.make_event_component_active <- function(component, components) {
  if (length(components) == 0L) return(TRUE)
  if (is.null(component) || identical(component, "__default__")) return(TRUE)
  comp_chr <- as.character(component)
  if (length(comp_chr) == 0L) return(TRUE)
  comp_val <- comp_chr[[1]]
  if (!nzchar(comp_val) || is.na(comp_val)) return(TRUE)
  comp_val %in% components
}

.make_event_density_fast_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  pool_def <- prep[["pools"]][[source_id]]
  if (!is.null(pool_def)) {
    return(function(t, component) {
      t_vals <- as.numeric(t)
      if (length(t_vals) == 0L) return(numeric(0))
      vapply(t_vals, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        .pool_density_fast_value(prep, source_id, component, tt)
      }, numeric(1))
    })
  }
  acc_def <- prep[["accumulators"]][[source_id]]
  if (is.null(acc_def)) return(NULL)
  components <- acc_def[["components"]] %||% character(0)
  function(t, component) {
    active <- .make_event_component_active(component, components)
    if (!active) return(rep(0.0, length(t)))
    vapply(t, function(tt) .acc_density(acc_def, tt), numeric(1))
  }
}

.make_event_survival_fast_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  pool_def <- prep[["pools"]][[source_id]]
  if (!is.null(pool_def)) {
    return(function(t, component) {
      t_vals <- as.numeric(t)
      if (length(t_vals) == 0L) return(numeric(0))
      vapply(t_vals, function(tt) {
        .pool_survival_fast_value(prep, source_id, component, tt)
      }, numeric(1))
    })
  }
  acc_def <- prep[["accumulators"]][[source_id]]
  if (is.null(acc_def)) return(NULL)
  components <- acc_def[["components"]] %||% character(0)
  function(t, component) {
    active <- .make_event_component_active(component, components)
    if (!active) return(rep(1.0, length(t)))
    vals <- vapply(t, function(tt) .acc_survival(acc_def, tt), numeric(1))
    out <- vals
    out[!is.finite(out)] <- 0.0
    out[out < 0] <- 0.0
    out[out > 1] <- 1.0
    out
  }
}

.make_event_cdf_fast_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  pool_def <- prep[["pools"]][[source_id]]
  if (!is.null(pool_def)) {
    surv_fast <- function(t, component) {
      t_vals <- as.numeric(t)
      if (length(t_vals) == 0L) return(numeric(0))
      vapply(t_vals, function(tt) .pool_survival_fast_value(prep, source_id, component, tt), numeric(1))
    }
    return(function(t, component) {
      surv_vals <- surv_fast(t, component)
      if (length(surv_vals) == 0L) return(numeric(0))
      out <- 1.0 - surv_vals
      out[!is.finite(out)] <- 0.0
      out[out < 0] <- 0.0
      out[out > 1] <- 1.0
      out
    })
  }
  acc_def <- prep[["accumulators"]][[source_id]]
  if (is.null(acc_def)) return(NULL)
  components <- acc_def[["components"]] %||% character(0)
  function(t, component) {
    active <- .make_event_component_active(component, components)
    if (!active) return(rep(0.0, length(t)))
    surv_vals <- vapply(t, function(tt) .acc_survival(acc_def, tt), numeric(1))
    out <- 1.0 - surv_vals
    out[!is.finite(out)] <- 0.0
    out[out < 0] <- 0.0
    out[out > 1] <- 1.0
    out
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

.make_and_cdf_fast_fn <- function(child_nodes, prep) {
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fast_fn")
  if (!all(vapply(child_cdf, is.function, logical(1)))) return(NULL)
  function(t, component) {
    vapply(t, function(tt) {
      prod_val <- 1.0
      for (fn in child_cdf) {
        val_vec <- fn(tt, component)
        val <- if (length(val_vec) > 0L) val_vec[[1]] else 0.0
        prod_val <- prod_val * val
        if (!is.finite(prod_val) || prod_val == 0) break
      }
      prod_val
    }, numeric(1))
  }
}

.make_and_survival_fast_fn <- function(cdf_fast_fn) {
  if (!is.function(cdf_fast_fn)) return(NULL)
  function(t, component) {
    vals <- cdf_fast_fn(t, component)
    out <- 1.0 - vals
    out[!is.finite(out)] <- 0.0
    out[out < 0] <- 0.0
    out[out > 1] <- 1.0
    out
  }
}

.make_and_density_fast_fn <- function(child_nodes, prep) {
  child_density <- lapply(child_nodes, `[[`, "density_fast_fn")
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fast_fn")
  if (!all(vapply(child_density, is.function, logical(1)))) return(NULL)
  if (!all(vapply(child_cdf, is.function, logical(1)))) return(NULL)
  function(t, component) {
    vapply(t, function(tt) {
      total <- 0.0
      n <- length(child_density)
      for (i in seq_len(n)) {
        dens_vec <- child_density[[i]](tt, component)
        dens_i <- if (length(dens_vec) > 0L) dens_vec[[1]] else 0.0
        if (!is.finite(dens_i) || dens_i <= 0) next
        prod_cdf <- 1.0
        for (j in seq_len(n)) {
          if (i == j) next
          Fj_vec <- child_cdf[[j]](tt, component)
          Fj <- if (length(Fj_vec) > 0L) Fj_vec[[1]] else 0.0
          prod_cdf <- prod_cdf * Fj
          if (!is.finite(prod_cdf) || prod_cdf == 0) break
        }
        if (!is.finite(prod_cdf) || prod_cdf == 0) next
        total <- total + dens_i * prod_cdf
      }
      if (!is.finite(total) || total < 0) total <- 0.0
      total
    }, numeric(1))
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

.make_and_density_fn <- function(child_nodes, prep) {
  child_density <- lapply(child_nodes, `[[`, "density_fn")
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fn")
  if (!all(vapply(child_density, is.function, logical(1))) ||
      !all(vapply(child_cdf, is.function, logical(1)))) {
    return(NULL)
  }
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    total <- 0.0
    n <- length(child_density)
    for (i in seq_len(n)) {
      dens_i <- child_density[[i]](t, component, forced_complete, forced_survive)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_cdf <- 1.0
      for (j in seq_len(n)) {
        if (i == j) next
        Fj <- child_cdf[[j]](t, component, forced_complete, forced_survive)
        prod_cdf <- prod_cdf * Fj
        if (!is.finite(prod_cdf) || prod_cdf == 0) break
      }
      if (!is.finite(prod_cdf) || prod_cdf == 0) next
      total <- total + dens_i * prod_cdf
    }
    if (!is.finite(total) || total < 0) total <- 0.0
    total
  }
}

.make_or_density_fn <- function(child_nodes, prep) {
  child_density <- lapply(child_nodes, `[[`, "density_fn")
  child_surv <- lapply(child_nodes, `[[`, "surv_fn")
  if (!all(vapply(child_density, is.function, logical(1))) ||
      !all(vapply(child_surv, is.function, logical(1)))) {
    return(NULL)
  }
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    total <- 0.0
    n <- length(child_density)
    for (i in seq_len(n)) {
      dens_i <- child_density[[i]](t, component, forced_complete, forced_survive)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_surv <- 1.0
      for (j in seq_len(n)) {
        if (i == j) next
        Sj <- child_surv[[j]](t, component, forced_complete, forced_survive)
        prod_surv <- prod_surv * Sj
        if (!is.finite(prod_surv) || prod_surv == 0) break
      }
      if (!is.finite(prod_surv) || prod_surv == 0) next
      total <- total + dens_i * prod_surv
    }
    if (!is.finite(total) || total < 0) total <- 0.0
    total
  }
}

.make_or_survival_fast_fn <- function(child_nodes, prep) {
  child_surv <- lapply(child_nodes, `[[`, "surv_fast_fn")
  if (!all(vapply(child_surv, is.function, logical(1)))) return(NULL)
  function(t, component) {
    vapply(t, function(tt) {
      prod_val <- 1.0
      for (fn in child_surv) {
        val_vec <- fn(tt, component)
        val <- if (length(val_vec) > 0L) val_vec[[1]] else 1.0
        prod_val <- prod_val * val
        if (!is.finite(prod_val) || prod_val == 0) break
      }
      prod_val
    }, numeric(1))
  }
}

.make_or_cdf_fast_fn <- function(surv_fast_fn) {
  if (!is.function(surv_fast_fn)) return(NULL)
  function(t, component) {
    surv_vals <- surv_fast_fn(t, component)
    out <- 1.0 - surv_vals
    out[!is.finite(out)] <- 0.0
    out[out < 0] <- 0.0
    out[out > 1] <- 1.0
    out
  }
}

.make_or_density_fast_fn <- function(child_nodes, prep) {
  child_density <- lapply(child_nodes, `[[`, "density_fast_fn")
  child_surv <- lapply(child_nodes, `[[`, "surv_fast_fn")
  if (!all(vapply(child_density, is.function, logical(1)))) return(NULL)
  if (!all(vapply(child_surv, is.function, logical(1)))) return(NULL)
  function(t, component) {
    vapply(t, function(tt) {
      total <- 0.0
      n <- length(child_density)
      for (i in seq_len(n)) {
        dens_vec <- child_density[[i]](tt, component)
        dens_i <- if (length(dens_vec) > 0L) dens_vec[[1]] else 0.0
        if (!is.finite(dens_i) || dens_i <= 0) next
        prod_surv <- 1.0
        for (j in seq_len(n)) {
          if (i == j) next
          Sj_vec <- child_surv[[j]](tt, component)
          Sj <- if (length(Sj_vec) > 0L) Sj_vec[[1]] else 1.0
          prod_surv <- prod_surv * Sj
          if (!is.finite(prod_surv) || prod_surv == 0) break
        }
        if (!is.finite(prod_surv) || prod_surv == 0) next
        total <- total + dens_i * prod_surv
      }
      if (!is.finite(total) || total < 0) total <- 0.0
      total
    }, numeric(1))
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

guard_scope_ids <- function(expr, prep) {
  prep_idx <- .prep_id_index(prep)
  if (is.null(prep_idx)) {
    idx_map <- prep[[".id_index"]] %||% NULL
    if (!is.null(idx_map)) {
      prep_tmp <- prep
      prep_tmp[[".runtime"]] <- list(id_index = idx_map)
      prep <- prep_tmp
    }
  }
  base_ids <- .expr_sources(expr, prep)
  gather_label_ids <- function(node) {
    if (is.null(node) || is.null(node[["kind"]])) return(integer(0))
    kind <- node[["kind"]]
    if (identical(kind, "event")) {
      lab <- node[["source"]]
      id <- .label_to_id(prep, lab)
      if (is.na(id)) integer(0) else id
    } else if (kind %in% c("and", "or")) {
      args <- node[["args"]] %||% list()
      if (length(args) == 0L) return(integer(0))
      unique(unlist(lapply(args, gather_label_ids)))
    } else if (identical(kind, "guard")) {
      ref_ids <- gather_label_ids(node[["reference"]])
      blk_ids <- gather_label_ids(node[["blocker"]])
      unl_ids <- unique(unlist(lapply(node[["unless"]] %||% list(), gather_label_ids)))
      unique(c(ref_ids, blk_ids, unl_ids))
    } else if (identical(kind, "not")) {
      gather_label_ids(node[["arg"]])
    } else {
      integer(0)
    }
  }
  label_ids <- gather_label_ids(expr)
  ids <- c(base_ids, label_ids)
  if (length(ids) == 0L) return(integer(0))
  ids <- unique(ids)
  ids[!is.na(ids)]
}

.guard_prepare_for_ids <- function(prep) {
  if (!is.null(.prep_id_index(prep))) return(prep)
  idx_map <- prep[[".id_index"]] %||% NULL
  if (is.null(idx_map)) return(prep)
  prep_stub <- prep
  prep_stub[[".runtime"]] <- list(id_index = idx_map)
  prep_stub
}

.guard_filter_forced_ids <- function(ids, scope_ids, prep) {
  prep_ids <- .guard_prepare_for_ids(prep)
  if (length(scope_ids) == 0L) {
    return(.coerce_forced_ids(prep_ids, integer(0)))
  }
  if (length(ids) == 0L) {
    return(.coerce_forced_ids(prep_ids, integer(0)))
  }
  raw <- as.integer(ids)
  keep <- raw[!is.na(raw) & raw %in% scope_ids]
  .coerce_forced_ids(prep_ids, keep)
}

.make_guard_cdf_fn <- function(expr, prep) {
  function(t, component, forced_complete, forced_survive) {
    if (t <= 0) return(0.0)
    if (!is.finite(t)) return(1.0)
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    scope_ids <- guard_scope_ids(expr, prep)
    forced_complete <- .guard_filter_forced_ids(forced_complete, scope_ids, prep)
    forced_survive <- .guard_filter_forced_ids(forced_survive, scope_ids, prep)
    local_state <- .eval_state_create()
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        val <- .guard_density_at(
          expr, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          state = local_state
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

.make_guard_density_fn <- function(expr, prep) {
  local_state <- .eval_state_create()
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    scope_ids <- guard_scope_ids(expr, prep)
    forced_complete <- .guard_filter_forced_ids(forced_complete, scope_ids, prep)
    forced_survive <- .guard_filter_forced_ids(forced_survive, scope_ids, prep)
    val <- .guard_density_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = local_state
    )
    as.numeric(val)
  }
}

.guard_density_fast_builder <- function(ref_density_fast, block_surv_fast, block_surv) {
  if (!is.function(ref_density_fast)) return(NULL)
  function(t, component) {
    t_vals <- as.numeric(t)
    if (length(t_vals) == 0L) return(numeric(0))
    vapply(t_vals, function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      dens_ref_vec <- ref_density_fast(tt, component)
      dens_ref <- if (length(dens_ref_vec) == 0L) 0.0 else dens_ref_vec[[1]]
      if (!is.finite(dens_ref) || dens_ref <= 0) return(0.0)
      surv_block <- if (is.function(block_surv_fast)) {
        surv_vec <- block_surv_fast(tt, component)
        if (length(surv_vec) == 0L) 1.0 else surv_vec[[1]]
      } else if (is.function(block_surv)) {
        block_surv(tt, component, integer(0), integer(0))
      } else {
        1.0
      }
      if (!is.finite(surv_block) || surv_block <= 0) return(0.0)
      total <- dens_ref * surv_block
      if (!is.finite(total) || total <= 0) return(0.0)
      total
    }, numeric(1))
  }
}

.guard_fast_survival_from_density <- function(density_fast) {
  if (!is.function(density_fast)) return(list(cdf = NULL, survival = NULL))
  guard_cdf <- function(t, component) {
    t_vals <- as.numeric(t)
    if (length(t_vals) == 0L) return(numeric(0))
    vapply(t_vals, function(tt) {
      if (tt <= 0) return(0.0)
      if (!is.finite(tt)) return(1.0)
      val <- tryCatch(
        stats::integrate(
          function(u) density_fast(u, component),
          lower = 0,
          upper = tt,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]],
        error = function(e) 0.0
      )
      out <- as.numeric(val)
      if (!is.finite(out)) out <- 0.0
      max(0.0, min(1.0, out))
    }, numeric(1))
  }
  guard_survival <- function(t, component) {
    vals <- guard_cdf(t, component)
    if (length(vals) == 0L) return(vals)
    out <- 1.0 - vals
    out[!is.finite(out)] <- 0.0
    out[out < 0] <- 0.0
    out[out > 1] <- 1.0
    out
  }
  list(cdf = guard_cdf, survival = guard_survival)
}

.make_guard_fast_fns <- function(expr, prep,
                                 ref_node_override = NULL,
                                 block_node_override = NULL,
                                 unless_nodes_override = NULL) {
  reference <- expr[["reference"]]
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  # Fast path is only valid when there are no protectors; otherwise we need
  # the full integral for effective blocker survival.
  if (length(unless_list) > 0L) {
    return(list(density = NULL, survival = NULL, cdf = NULL))
  }
  ref_node <- ref_node_override %||% .expr_lookup_compiled(reference, prep)
  block_node <- if (!is.null(blocker)) {
    block_node_override %||% .expr_lookup_compiled(blocker, prep)
  } else {
    NULL
  }
  density_fast <- .guard_density_fast_builder(
    .node_ops_get(ref_node, "density_fast", NULL),
    .node_ops_get(block_node, "survival_fast", NULL),
    .node_ops_get(block_node, "survival", NULL)
  )
  if (!is.function(density_fast)) {
    return(list(density = NULL, survival = NULL, cdf = NULL))
  }
  fast_surv_cdf <- .guard_fast_survival_from_density(density_fast)
  list(density = density_fast,
       survival = fast_surv_cdf$survival,
       cdf = fast_surv_cdf$cdf)
}

.node_cdf_cond <- function(node, t, prep, component,
                           forced_complete = integer(0),
                           forced_survive = integer(0),
                           state = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  cdf_fast_fn <- .node_ops_get(node, "cdf_fast", node$cdf_fast_fn)
  cdf_fn <- .node_ops_get(node, "cdf", node$cdf_fn)
  state_env <- if (is.environment(state)) state else NULL
  entry_extra <- if (!is.null(state_env)) list(time_id = .eval_state_time_id(state_env, t)) else NULL
  entry <- .node_eval_entry(state, node, component, t,
                            forced_complete = forced_complete,
                            forced_survive = forced_survive,
                            extra = entry_extra)
  cached <- .state_entry_get(entry, "cdf")
  if (!is.null(cached)) return(cached)
  if (length(forced_complete) == 0L && length(forced_survive) == 0L &&
      is.function(cdf_fast_fn)) {
    val <- cdf_fast_fn(t, component)
    out <- if (length(val) == 0L) 0.0 else as.numeric(val[[1]])
    .state_entry_set(entry, "cdf", out)
    return(out)
  }
  if (is.function(cdf_fn)) {
    val <- cdf_fn(t, component, forced_complete, forced_survive)
    out <- if (length(val) == 0L) 0.0 else as.numeric(val)
    .state_entry_set(entry, "cdf", out)
    return(out)
  }
  expr <- node$expr
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    out <- .event_cdf_at(
      prep, expr[["source"]], component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    .state_entry_set(entry, "cdf", out)
    return(out)
  }
  if (identical(kind, "and")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) {
      .state_entry_set(entry, "cdf", 0.0)
      return(0.0)
    }
    prod_val <- 1.0
    for (cid in child_ids) {
      child_node <- .compiled_node_fetch(prep, cid)
      if (is.null(child_node)) next
      Fc <- .node_cdf_cond(child_node, t, prep, component,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive,
                           state = state)
      prod_val <- prod_val * Fc
      if (!is.finite(prod_val) || prod_val == 0) break
    }
    if (!is.finite(prod_val) || prod_val < 0) prod_val <- 0.0
    if (prod_val > 1) prod_val <- 1.0
    .state_entry_set(entry, "cdf", prod_val)
    return(prod_val)
  }
  if (identical(kind, "or")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) {
      .state_entry_set(entry, "cdf", 0.0)
      return(0.0)
    }
    prod_surv <- 1.0
    for (cid in child_ids) {
      child_node <- .compiled_node_fetch(prep, cid)
      if (is.null(child_node)) next
      Sc <- .node_survival_cond(child_node, t, prep, component,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                state = state)
      prod_surv <- prod_surv * Sc
      if (!is.finite(prod_surv) || prod_surv == 0) break
    }
    if (!is.finite(prod_surv) || prod_surv < 0) prod_surv <- 0.0
    if (prod_surv > 1) prod_surv <- 1.0
    out <- 1.0 - prod_surv
    if (!is.finite(out) || out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    .state_entry_set(entry, "cdf", out)
    return(out)
  }
  if (identical(kind, "guard")) {
    surv <- .node_survival_cond(node, t, prep, component,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                state = state)
    out <- 1.0 - surv
    if (!is.finite(out) || out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    .state_entry_set(entry, "cdf", out)
    return(out)
  }
  .state_entry_set(entry, "cdf", 0.0)
  0.0
}

.node_survival_cond <- function(node, t, prep, component,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               state = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  surv_fast_fn <- .node_ops_get(node, "survival_fast", node$surv_fast_fn)
  surv_fn <- .node_ops_get(node, "survival", node$surv_fn)
  state_env <- if (is.environment(state)) state else NULL
  entry <- NULL
  entry_extra <- if (!is.null(state_env)) list(time_id = .eval_state_time_id(state_env, t)) else NULL
  if (!is.null(state_env)) {
    entry <- .node_eval_entry(
      state_env, node, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      extra = entry_extra,
      init = FALSE
    )
  } else {
    entry <- NULL
  }
  cached <- .state_entry_get(entry, "survival")
  if (!is.null(cached)) return(cached)
  expr <- node$expr
  kind <- expr[["kind"]]
  if (identical(kind, "guard")) {
    if (!is.null(state_env) && is.null(entry)) {
      entry <- .node_eval_entry(
        state_env, node, component, t,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        extra = entry_extra,
        init = TRUE
      )
    }
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .guard_density_at(
          node$expr, uu, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          state = state
        )
      }, numeric(1))
    }
    val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                     rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                     stop.on.error = FALSE)[["value"]],
                    error = function(e) 0.0)
    out <- 1.0 - as.numeric(val)
    if (!is.finite(out) || out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    .state_entry_set(entry, "survival", out)
    return(out)
  }
  fast_unforced <- length(forced_complete) == 0L && length(forced_survive) == 0L &&
    is.function(surv_fast_fn)
  if (fast_unforced) {
    val <- surv_fast_fn(t, component)
    out <- if (length(val) == 0L) 1.0 else as.numeric(val[[1]])
    .state_entry_set(entry, "survival", out)
    return(out)
  }
  if (!is.null(state_env) && is.null(entry)) {
    entry <- .node_eval_entry(
      state_env, node, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      extra = entry_extra,
      init = TRUE
    )
  }
  if (is.function(surv_fn)) {
    val <- surv_fn(t, component, forced_complete, forced_survive)
    out <- if (length(val) == 0L) 1.0 else as.numeric(val)
    .state_entry_set(entry, "survival", out)
    return(out)
  }
  if (identical(kind, "event")) {
    out <- .event_survival_at(
      prep, expr[["source"]], component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    .state_entry_set(entry, "survival", out)
    return(out)
  }
  if (identical(kind, "and")) {
    out <- 1.0 - .node_cdf_cond(
      node, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    if (!is.finite(out) || out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    .state_entry_set(entry, "survival", out)
    return(out)
  }
  if (identical(kind, "or")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) {
      .state_entry_set(entry, "survival", 1.0)
      return(1.0)
    }
    prod_val <- 1.0
    for (cid in child_ids) {
      child_node <- .compiled_node_fetch(prep, cid)
      if (is.null(child_node)) next
      Sv <- .node_survival_cond(child_node, t, prep, component,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                state = state)
      prod_val <- prod_val * Sv
      if (!is.finite(prod_val) || prod_val == 0) break
    }
    if (!is.finite(prod_val) || prod_val < 0) prod_val <- 0.0
    if (prod_val > 1) prod_val <- 1.0
    .state_entry_set(entry, "survival", prod_val)
    return(prod_val)
  }
  .state_entry_set(entry, "survival", 1.0)
  1.0
}

.node_density <- function(node, t, prep, component,
                          forced_complete = integer(0),
                          forced_survive = integer(0),
                          state = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  density_fast_fn <- .node_ops_get(node, "density_fast", node$density_fast_fn)
  density_fn <- .node_ops_get(node, "density", node$density_fn)
  state_env <- if (is.environment(state)) state else NULL
  entry_extra <- if (!is.null(state_env)) list(time_id = .eval_state_time_id(state_env, t)) else NULL
  if (!is.null(state_env)) {
    entry <- .node_eval_entry(
      state_env, node, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      extra = entry_extra,
      init = FALSE
    )
  } else {
    entry <- NULL
  }
  cached <- .state_entry_get(entry, "density")
  if (!is.null(cached)) return(cached)
  fast_unforced <- length(forced_complete) == 0L && length(forced_survive) == 0L &&
    is.function(density_fast_fn)
  if (fast_unforced) {
    val <- density_fast_fn(t, component)
    out <- if (length(val) == 0L) 0.0 else as.numeric(val[[1]])
    .state_entry_set(entry, "density", out)
    return(out)
  }
  if (!is.null(state_env) && is.null(entry)) {
    entry <- .node_eval_entry(
      state_env, node, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      extra = entry_extra,
      init = TRUE
    )
  }
  if (is.function(density_fn)) {
    val <- density_fn(t, component, forced_complete, forced_survive)
    out <- if (length(val) == 0L) 0.0 else as.numeric(val)
    .state_entry_set(entry, "density", out)
    return(out)
  }
  expr <- node$expr
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    out <- .event_density_at(
      prep, expr[["source"]], component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    .state_entry_set(entry, "density", out)
    return(out)
  }
  if (identical(kind, "and")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) {
      .state_entry_set(entry, "density", 0.0)
      return(0.0)
    }
    total <- 0.0
    for (i in seq_along(child_ids)) {
      cid <- child_ids[[i]]
      child_node <- .compiled_node_fetch(prep, cid)
      if (is.null(child_node)) next
      dens_i <- .node_density(child_node, t, prep, component,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive,
                              state = state)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_cdf <- 1.0
      for (j in seq_along(child_ids)) {
        if (i == j) next
        other_node <- .compiled_node_fetch(prep, child_ids[[j]])
        if (is.null(other_node)) next
        Fj <- .node_cdf_cond(other_node, t, prep, component,
                             forced_complete = forced_complete,
                             forced_survive = forced_survive,
                             state = state)
        prod_cdf <- prod_cdf * Fj
        if (!is.finite(prod_cdf) || prod_cdf == 0) break
      }
      if (!is.finite(prod_cdf) || prod_cdf == 0) next
      total <- total + dens_i * prod_cdf
    }
    if (!is.finite(total) || total < 0) total <- 0.0
    .state_entry_set(entry, "density", total)
    return(total)
  }
  if (identical(kind, "or")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) {
      .state_entry_set(entry, "density", 0.0)
      return(0.0)
    }
    total <- 0.0
    for (i in seq_along(child_ids)) {
      cid <- child_ids[[i]]
      child_node <- .compiled_node_fetch(prep, cid)
      if (is.null(child_node)) next
      dens_i <- .node_density(child_node, t, prep, component,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive,
                              state = state)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_surv <- 1.0
      for (j in seq_along(child_ids)) {
        if (i == j) next
        other_node <- .compiled_node_fetch(prep, child_ids[[j]])
        if (is.null(other_node)) next
        Sj <- .node_survival_cond(other_node, t, prep, component,
                                  forced_complete = forced_complete,
                                  forced_survive = forced_survive,
                                  state = state)
        prod_surv <- prod_surv * Sj
        if (!is.finite(prod_surv) || prod_surv == 0) break
      }
      if (!is.finite(prod_surv) || prod_surv == 0) next
      total <- total + dens_i * prod_surv
    }
    if (!is.finite(total) || total < 0) total <- 0.0
    .state_entry_set(entry, "density", total)
    return(total)
  }
  if (identical(kind, "guard")) {
    out <- .guard_density_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    .state_entry_set(entry, "density", out)
    return(out)
  }
  .state_entry_set(entry, "density", 0.0)
  0.0
}

.node_scenarios_at <- function(node, t, prep, component,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               state = NULL) {
  scenario_fn <- .node_ops_get(node, "scenario", node$scenario_fn)
  if (is.function(scenario_fn)) {
    return(scenario_fn(t, component, forced_complete, forced_survive))
  }
  expr <- node$expr
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  state_env <- if (is.environment(state)) state else NULL
  entry_extra <- if (!is.null(state_env)) list(time_id = .eval_state_time_id(state_env, t)) else NULL
  entry <- .node_eval_entry(state, node, component, t,
                            forced_complete = forced_complete,
                            forced_survive = forced_survive,
                            extra = entry_extra)
  cached <- .state_entry_get(entry, "scenarios")
  if (!is.null(cached)) return(cached)
  store <- function(res) {
    .state_entry_set(entry, "scenarios", res)
    res
  }
  if (identical(kind, "event")) {
    weight <- .event_density_at(
      prep, expr[["source"]], component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    if (!is.finite(weight) || weight <= 0) return(store(list()))
    pool_scenarios <- attr(weight, "scenarios")
    source_idx <- .label_to_id(prep, expr[["source"]])
    source_ids <- if (is.na(source_idx)) integer(0) else source_idx
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
      return(store(out))
    }
    sc <- .make_scenario_record(
      prep,
      weight,
      .forced_union(prep, forced_complete, source_ids),
      forced_survive
    )
    return(store(if (is.null(sc)) list() else list(sc)))
  }
  if (identical(kind, "and")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) return(store(list()))
    child_nodes <- lapply(child_ids, function(id) .compiled_node_fetch(prep, id))
    child_sources <- lapply(child_nodes, function(n) n$sources %||% integer(0))
    out <- list()
    for (i in seq_along(child_nodes)) {
      child_node <- child_nodes[[i]]
      if (is.null(child_node)) next
      si <- .node_scenarios_at(child_node, t, prep, component,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive,
                               state = state)
      if (length(si) == 0L) next
      others <- seq_along(child_nodes)
      if (length(others) > 0L) others <- others[others != i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        ok <- TRUE
        if (length(others) > 0L) {
          for (j in others) {
            other_node <- child_nodes[[j]]
            if (is.null(other_node)) next
            Fj <- .node_cdf_cond(other_node, t, prep, component,
                                 forced_complete = fcomp,
                                 forced_survive = fsurv,
                                 state = state)
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
    return(store(out))
  }
  if (identical(kind, "or")) {
    child_ids <- node$args %||% integer(0)
    if (length(child_ids) == 0L) return(store(list()))
    child_nodes <- lapply(child_ids, function(id) .compiled_node_fetch(prep, id))
    child_sources <- lapply(child_nodes, function(n) n$sources %||% integer(0))
    out <- list()
    for (i in seq_along(child_nodes)) {
      child_node <- child_nodes[[i]]
      if (is.null(child_node)) next
      si <- .node_scenarios_at(child_node, t, prep, component,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive,
                               state = state)
      if (length(si) == 0L) next
      others <- seq_along(child_nodes)
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
    return(store(out))
  }
  if (identical(kind, "guard")) {
    reference_node <- .compiled_node_fetch(prep, node$reference_id)
    blocker_node <- .compiled_node_fetch(prep, node$blocker_id)
  unless_nodes <- if (length(node$unless_ids %||% integer(0)) > 0L) {
    lapply(node$unless_ids, function(uid) .compiled_node_fetch(prep, uid))
  } else {
    list()
  }
    ref_scen <- if (is.null(reference_node)) list() else .node_scenarios_at(
      reference_node, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
    if (length(ref_scen) == 0) return(store(list()))
    blocker_sources <- integer(0)
    blocker_expr <- node$expr[["blocker"]]
    if (!is.null(blocker_expr)) {
      blocker_sources <- .expr_sources(blocker_expr, prep)
      blocker_label <- .label_to_id(prep, blocker_expr[["source"]])
      if (!is.na(blocker_label)) {
        blocker_sources <- .forced_union(prep, blocker_sources, blocker_label)
      }
    }
    out <- list()
    for (sc in ref_scen) {
      if (is.null(sc)) next
      fc_ctx <- sc$forced_complete %||% integer(0)
      fs_ctx <- sc$forced_survive %||% integer(0)
      S_eff <- .guard_effective_survival(
        node$expr, t, prep, component,
        fc_ctx, fs_ctx, state,
        block_node = blocker_node,
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
        forced_survive = fsurv
      )
    }
    return(store(out))
  }
  store(list())
}

.guard_effective_survival <- function(guard_expr, t, prep, component,
                                      forced_complete, forced_survive,
                                      state, block_node = NULL,
                                      unless_nodes = list()) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  scope_ids <- guard_scope_ids(guard_expr, prep)
  forced_complete <- .guard_filter_forced_ids(forced_complete, scope_ids, prep)
  forced_survive <- .guard_filter_forced_ids(forced_survive, scope_ids, prep)
  # Cache S_eff per (guard signature, component, t, forced sets)
  comp_key <- .eval_state_component_key(component)
  time_key <- .eval_state_time_key(t)
  fc_key <- .eval_state_ids_key(forced_complete)
  fs_key <- .eval_state_ids_key(forced_survive)
  sig <- .expr_signature(guard_expr)
  cache_tag <- paste("guard_S_eff", sig, comp_key, time_key, fc_key, fs_key, sep = "|")
  cached <- .eval_state_get_extra(state, cache_tag)
  if (!is.null(cached)) return(as.numeric(cached))
  blocker <- guard_expr[["blocker"]]
  unless_list <- guard_expr[["unless"]] %||% list()
  if (is.null(blocker)) return(1.0)
  block_surv_fast <- .node_ops_get(block_node, "survival_fast", NULL)
  block_surv <- .node_ops_get(block_node, "survival", NULL)
  if (length(unless_list) == 0) {
    if (!is.null(block_node)) {
      if (length(forced_complete) == 0L && length(forced_survive) == 0L &&
          is.function(block_surv_fast)) {
        val <- block_surv_fast(t, component)
        return(if (length(val) == 0L) 1.0 else as.numeric(val[[1]]))
      }
      if (is.function(block_surv)) {
        val <- block_surv(t, component, forced_complete, forced_survive)
        return(if (length(val) == 0L) 1.0 else as.numeric(val))
      }
    }
    return(.eval_expr_survival_cond(
      blocker, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    ))
  }
  block_density_fast <- .node_ops_get(block_node, "density_fast", NULL)
  block_density_fn <- .node_ops_get(block_node, "density", NULL)
  block_fast_unforced <- length(forced_complete) == 0L && length(forced_survive) == 0L &&
    is.function(block_density_fast)
  block_has_density_fn <- is.function(block_density_fn)
  blocker_density <- function(u) {
    vapply(u, function(ui) {
      if (!is.finite(ui) || ui < 0) return(0.0)
      if (block_fast_unforced) {
        val <- block_density_fast(ui, component)
        return(if (length(val) == 0L) 0.0 else as.numeric(val[[1]]))
      }
      if (block_has_density_fn) {
        val <- block_density_fn(ui, component, forced_complete, forced_survive)
        return(if (length(val) == 0L) 0.0 else as.numeric(val))
      }
      .eval_expr_likelihood(
        blocker, ui, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        state = state
      )
    }, numeric(1))
  }
  protector_survival_product <- function(u) {
    if (length(unless_list) == 0) return(rep(1.0, length(u)))
    vapply(u, function(ui) {
      vals <- vapply(seq_along(unless_list), function(idx) {
        unl <- unless_list[[idx]]
        unl_node <- unless_nodes[[idx]]
        unl_surv_fast <- .node_ops_get(unl_node, "survival_fast", NULL)
        unl_surv <- .node_ops_get(unl_node, "survival", NULL)
        if (!is.null(unl_node)) {
          if (length(forced_complete) == 0L && length(forced_survive) == 0L &&
              is.function(unl_surv_fast)) {
            val <- unl_surv_fast(ui, component)
            return(if (length(val) == 0L) 1.0 else as.numeric(val[[1]]))
          }
          if (is.function(unl_surv)) {
            val <- unl_surv(ui, component, forced_complete, forced_survive)
            return(if (length(val) == 0L) 1.0 else as.numeric(val))
          }
        }
        .eval_expr_survival_cond(
          unl, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          state = state
        )
      }, numeric(1))
      prod(vals)
    }, numeric(1))
  }
  if (!is.finite(t) || t <= 0) {
    .eval_state_set_extra(state, cache_tag, 1.0)
    return(1.0)
  }
  integrand <- function(u) blocker_density(u) * protector_survival_product(u)
  native_guard <- .lik_native_fn("guard_effective_survival_cpp")
  val <- tryCatch(
    native_guard(
      integrand,
      as.numeric(t),
      .integrate_rel_tol(),
      .integrate_abs_tol(),
      12L
    ),
    error = function(e) 1.0
  )
  if (!is.finite(val)) val <- 1.0
  val <- max(0.0, min(1.0, as.numeric(val)))
  .eval_state_set_extra(state, cache_tag, val)
  val
}

.guard_scenarios_slow <- function(expr, t, prep, component,
                                  forced_complete, forced_survive,
                                  state = NULL) {
  compiled_guard <- .expr_lookup_compiled(expr, prep)
  if (!is.null(compiled_guard)) {
    return(.node_scenarios_at(
      compiled_guard, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    ))
  }
  reference <- expr[["reference"]]
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  ref_scen <- .expr_scenarios_at(
    reference, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = state
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
    fc_ctx <- sc$forced_complete %||% integer(0)
    fs_ctx <- sc$forced_survive %||% integer(0)
    S_eff <- .guard_effective_survival(
      expr, t, prep, component,
      fc_ctx, fs_ctx, state,
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
      forced_survive = fsurv
    )
  }
  out
}

.guard_density_at <- function(expr, t, prep, component,
                              forced_complete, forced_survive,
                              state = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  scope_ids <- guard_scope_ids(expr, prep)
  forced_complete <- .guard_filter_forced_ids(forced_complete, scope_ids, prep)
  forced_survive <- .guard_filter_forced_ids(forced_survive, scope_ids, prep)
  # Direct formulation: f_guard(t) = f_ref(t) * S_eff_blocker(t)
  reference <- expr[["reference"]]
  if (is.null(reference)) return(0.0)
  dens_ref <- .eval_expr_likelihood(
    reference, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = state
  )
  if (!is.finite(dens_ref) || dens_ref <= 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  # Reuse compiled nodes when available for efficiency
  block_node <- NULL
  unless_nodes <- list()
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  if (!is.null(blocker)) block_node <- .expr_lookup_compiled(blocker, prep)
  if (length(unless_list) > 0L) unless_nodes <- lapply(unless_list, function(unl) .expr_lookup_compiled(unl, prep))
  S_eff <- .guard_effective_survival(
    expr, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = state,
    block_node = block_node,
    unless_nodes = unless_nodes
  )
  out <- as.numeric(dens_ref) * as.numeric(S_eff)
  if (!is.finite(out) || out < 0) out <- 0.0
  attr(out, "scenarios") <- list()
  out
}

.make_guard_scenario_fn <- function(expr, prep) {
  reference <- expr[["reference"]]
  ref_node <- .expr_lookup_compiled(reference, prep)
  ref_scen_fn <- .node_ops_get(ref_node, "scenario", NULL)
  if (!is.function(ref_scen_fn)) return(NULL)
  blocker <- expr[["blocker"]]
  unless_list <- expr[["unless"]] %||% list()
  block_node <- if (!is.null(blocker)) .expr_lookup_compiled(blocker, prep) else NULL
  unless_nodes <- lapply(unless_list, function(unl) .expr_lookup_compiled(unl, prep))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    state_local <- .eval_state_create()
    ref_scen <- ref_scen_fn(t, component, forced_complete, forced_survive)
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
      fc_ctx <- sc$forced_complete
      fs_ctx <- sc$forced_survive
      S_eff <- .guard_effective_survival(
        expr, t, prep, component,
        fc_ctx, fs_ctx, state_local,
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
        forced_survive = fsurv
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
      forced_complete = fcomp,
      forced_survive = fsurv
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
                               state = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  entry <- .eval_state_entry(
    state,
    node_id,
    component,
    t,
    forced_complete = forced_complete,
    forced_survive = forced_survive
  )
  if (!is.null(entry)) {
    cached <- entry$scenarios
    if (!is.null(cached)) return(cached)
  }
  res <- if (!is.null(compiled)) {
    .node_scenarios_at(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  } else {
    .expr_scenarios_at_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  }
  if (!is.null(entry)) entry$scenarios <- res
  res
}

.expr_scenarios_at_slow <- function(expr, t, prep, component,
                                    forced_complete = integer(0),
                                    forced_survive = integer(0),
                                    state = NULL) {
  compiled <- .expr_lookup_compiled(expr, prep)
  if (!is.null(compiled)) {
    return(.node_scenarios_at(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    ))
  }
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  
  make_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return(NULL)
    list(
      weight = as.numeric(weight),
      forced_complete = fcomp,
      forced_survive = fsurv
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
        state = state)
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
            state = state)
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
        state = state)
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
      state = state
    ))
  }
  
  weight <- .eval_expr_likelihood(
    expr, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = state
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
                                state = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  entry <- .eval_state_entry(
    state,
    node_id,
    component,
    t,
    forced_complete = forced_complete,
    forced_survive = forced_survive
  )
  if (!is.null(entry)) {
    cached <- entry$cdf
    if (!is.null(cached)) return(cached)
  }
  res <- if (!is.null(compiled)) {
    .node_cdf_cond(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  } else {
    .eval_expr_cdf_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  }
  if (!is.null(entry)) entry$cdf <- res
  res
}

.eval_expr_cdf_cond_slow <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     state = NULL) {
  compiled <- .expr_lookup_compiled(expr, prep)
  if (!is.null(compiled)) {
    return(.node_cdf_cond(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    ))
  }
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
                                state = state)
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
                                     state = state)
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
          state = state
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
                                     state = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  entry <- .eval_state_entry(
    state,
    node_id,
    component,
    t,
    forced_complete = forced_complete,
    forced_survive = forced_survive
  )
  if (!is.null(entry)) {
    cached <- entry$survival
    if (!is.null(cached)) return(cached)
  }
  res <- if (!is.null(compiled)) {
    .node_survival_cond(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  } else {
    .eval_expr_survival_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    )
  }
  if (!is.null(entry)) entry$survival <- res
  res
}

.eval_expr_survival_cond_slow <- function(expr, t, prep, component,
                                          forced_complete = integer(0),
                                          forced_survive = integer(0),
                                          state = NULL) {
  compiled <- .expr_lookup_compiled(expr, prep)
  if (!is.null(compiled)) {
    return(.node_survival_cond(
      compiled, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      state = state
    ))
  }
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
                                     state = state))
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(1.0)
    prod_val <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     state = state)
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
        .guard_density_at(
          expr, uu, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          state = state
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
                                               state = NULL) {
  if (is.na(t) || !is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  state <- state %||% .eval_state_create()
  # Guard: compute density directly to avoid stale compiled closures and ensure
  # correct protector handling.
  if (!is.null(expr) && !is.null(expr$kind) && identical(expr$kind, "guard")) {
    dens_scalar <- .guard_density_at(
      expr, t, prep, component,
      forced_complete = integer(0),
      forced_survive = integer(0),
      state = state
    )
    if (!is.finite(dens_scalar) || dens_scalar <= 0) {
      dens_scalar <- 0.0
    } else if (length(competitor_exprs) > 0) {
      surv_prod <- .compute_survival_product(
        outcome_expr = expr,
        competitor_exprs = competitor_exprs,
        prep = prep,
        component = component,
        t = t,
        state = state
      )
      if (!is.finite(surv_prod) || surv_prod <= 0) {
        dens_scalar <- 0.0
      } else {
        dens_scalar <- dens_scalar * surv_prod
      }
    }
    attr(dens_scalar, "scenarios") <- list()
    return(dens_scalar)
  }
  compiled <- .expr_lookup_compiled(expr, prep)
  density_fn <- if (!is.null(compiled)) compiled$density_fn else NULL
  if (!is.null(compiled) && is.function(compiled$density_fast_fn)) {
    dens <- compiled$density_fast_fn(t, component)
  } else if (is.function(density_fn)) {
    dens <- density_fn(t, component, integer(0), integer(0))
  } else {
    dens <- NULL
  }
  if (!is.null(dens)) {
    dens_val <- as.numeric(dens)
    dens_scalar <- if (length(dens_val) > 0L) dens_val[[1]] else 0.0
    if (!is.finite(dens_scalar) || dens_scalar <= 0) {
      dens_scalar <- 0.0
    } else if (length(competitor_exprs) > 0) {
      surv_prod <- .compute_survival_product(
        outcome_expr = expr,
        competitor_exprs = competitor_exprs,
        prep = prep,
        component = component,
        t = t,
        state = state
      )
      if (!is.finite(surv_prod) || surv_prod <= 0) {
        dens_scalar <- 0.0
      } else {
        dens_scalar <- dens_scalar * surv_prod
      }
    }
    attr(dens_scalar, "scenarios") <- list()
    return(dens_scalar)
  }
  if (!is.null(compiled)) {
    dens_scalar <- .node_density(
      compiled, t, prep, component,
      forced_complete = integer(0),
      forced_survive = integer(0),
      state = state
    )
    if (!is.finite(dens_scalar) || dens_scalar <= 0) {
      dens_scalar <- 0.0
    } else if (length(competitor_exprs) > 0) {
      surv_prod <- .compute_survival_product(
        outcome_expr = expr,
        competitor_exprs = competitor_exprs,
        prep = prep,
        component = component,
        t = t,
        state = state
      )
      if (!is.finite(surv_prod) || surv_prod <= 0) {
        dens_scalar <- 0.0
      } else {
        dens_scalar <- dens_scalar * surv_prod
      }
    }
    attr(dens_scalar, "scenarios") <- list()
    return(dens_scalar)
  }
  scenarios <- .expr_scenarios_at(
    expr, t, prep, component,
    forced_complete = integer(0),
    forced_survive = integer(0),
    state = state)
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
        state = state
      )
      weight <- weight * surv_prod
    }
    total <- total + weight
  }
  if (length(scenarios) > 0) attr(total, "scenarios") <- scenarios
  total
}

# ============================================================================== 
# Fast-path helpers for competitor survival products
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

.expr_guard_free <- function(expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(TRUE)
  kind <- expr[["kind"]]
  if (identical(kind, "guard")) return(FALSE)
  if (kind %in% c("and", "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0L) return(TRUE)
    return(all(vapply(args, .expr_guard_free, logical(1))))
  }
  if (identical(kind, "not")) {
    return(.expr_guard_free(expr[["arg"]]))
  }
  TRUE
}

.cluster_is_guard_free <- function(exprs) {
  if (length(exprs) == 0L) return(TRUE)
  all(vapply(exprs, .expr_guard_free, logical(1)))
}

.compute_survival_product <- function(outcome_expr, competitor_exprs, prep, component, t,
                                     state = NULL) {
  if (length(competitor_exprs) == 0) return(1.0)
  fast_val <- .fast_survival_product(competitor_exprs, prep, component, t)
  if (!is.null(fast_val)) {
    fast_vec <- as.numeric(fast_val)
    fast_scalar <- if (length(fast_vec) > 0L) fast_vec[[1]] else NA_real_
    if (!is.finite(fast_scalar) || is.na(fast_scalar) || fast_scalar <= 0) return(0.0)
    if (fast_scalar > 1) fast_scalar <- 1.0
    return(fast_scalar)
  }
  state <- state %||% .eval_state_create()
  clusters <- .cluster_competitors(competitor_exprs, prep)
  if (length(clusters) == 0) return(1.0)
  prod_val <- 1.0
  for (cluster in clusters) {
    cluster_exprs <- cluster$exprs %||% list()
    cluster_indices <- cluster$indices %||% seq_along(cluster_exprs)
    if (.cluster_is_guard_free(cluster_exprs)) {
      cluster_val <- .compute_survival_product_cluster(cluster_exprs, prep, component, t, state)
      prod_val <- prod_val * cluster_val
      if (!is.finite(prod_val) || prod_val == 0) return(0.0)
      next
    }
    component_key <- .eval_state_component_key(component)
    time_key <- .eval_state_time_key(t)
    cache_key <- paste(
      "joint_cluster",
      paste(cluster_indices, collapse = ","),
      component_key,
      time_key,
      sep = "|"
    )
    cached <- .eval_state_get_extra(state, cache_key)
    if (!is.null(cached)) {
      cluster_val <- cached
    } else {
      cluster_val <- .joint_survival_for_cluster(cluster_exprs, prep, component, t, state)
      if (!is.finite(cluster_val) || cluster_val < 0) cluster_val <- 0.0
      if (cluster_val > 1) cluster_val <- 1.0
      .eval_state_set_extra(state, cache_key, cluster_val)
    }
    prod_val <- prod_val * cluster_val
    if (!is.finite(prod_val) || prod_val == 0) return(0.0)
  }
  if (!is.finite(prod_val) || prod_val < 0) prod_val <- 0.0
  if (prod_val > 1) prod_val <- 1.0
  prod_val
}

# Compute joint survival for a cluster of dependent competitor expressions at time t
# using scenario-based conditioning. This is fully general (supports guards on guards).
.joint_survival_for_cluster <- function(cluster_exprs, prep, component, t, state) {
  n <- length(cluster_exprs)
  if (n == 0) return(1.0)
  if (n == 1) {
    return(.eval_expr_survival_cond(
      cluster_exprs[[1]], t, prep, component,
      forced_complete = integer(0),
      forced_survive = integer(0),
      state = state
    ))
  }
  # Chain-rule conditional product with heuristic ordering
  # Precompute metadata for ordering
  meta <- lapply(cluster_exprs, function(ex) {
    node <- .expr_lookup_compiled(ex, prep)
    sources_chr <- .expr_sources(ex, prep)
    sources_ids <- .labels_to_ids(prep, sources_chr)
    is_fast <- if (!is.null(node)) is.function(node$surv_fast_fn) else FALSE
    scenario_sensitive <- if (!is.null(node)) isTRUE(node$scenario_sensitive) else TRUE
    list(expr = ex, node = node, src = sources_ids, is_fast = is_fast, sens = scenario_sensitive)
  })
  # Overlap heuristic: fewer overlaps and smaller src sets first; fast and non-sensitive first
  all_src <- unique(unlist(lapply(meta, function(m) m$src)))
  overlap_counts <- vapply(meta, function(m) length(intersect(m$src, all_src)), integer(1))
  ord <- order(
    vapply(meta, function(m) if (m$is_fast) 0L else 1L, integer(1)),
    vapply(meta, function(m) if (m$sens) 1L else 0L, integer(1)),
    vapply(meta, function(m) length(m$src), integer(1)),
    overlap_counts,
    decreasing = FALSE
  )
  meta <- meta[ord]
  forced_surv <- integer(0)
  prod_val <- 1.0
  for (m in meta) {
    s_i <- .eval_expr_survival_cond(
      m$expr, t, prep, component,
      forced_complete = integer(0),
      forced_survive = forced_surv,
      state = state
    )
    if (!is.finite(s_i) || s_i <= 0) return(0.0)
    prod_val <- prod_val * s_i
    if (!is.finite(prod_val) || prod_val <= 0) return(0.0)
    if (length(m$src) > 0) {
      # Append-only when disjoint; fallback to union when overlap exists
      if (length(forced_surv) == 0L) {
        forced_surv <- as.integer(m$src)
      } else if (all(!(m$src %in% forced_surv))) {
        forced_surv <- c(forced_surv, as.integer(m$src))
      } else {
        forced_surv <- .forced_union(prep, forced_surv, m$src)
      }
    }
  }
  if (!is.finite(prod_val) || prod_val < 0) prod_val <- 0.0
  if (prod_val > 1) prod_val <- 1.0
  prod_val
}

.compute_survival_product_cluster <- function(exprs, prep, component, t, state) {
  if (length(exprs) == 0) return(1.0)
  state <- state %||% .eval_state_create()
  clusters <- .cluster_competitors(exprs, prep)
  prod_val <- 1.0
  for (cluster in clusters) {
    cluster_exprs <- cluster$exprs
    cluster_indices <- cluster$indices %||% seq_along(cluster_exprs)
    component_key <- if (is.null(component)) "__default__" else as.character(component)[[1]]
    time_key <- .eval_state_time_key(t)
    cache_key <- paste(
      "cluster",
      paste(cluster_indices, collapse = ","),
      component_key,
      time_key,
      sep = "|"
    )
    cached <- .eval_state_get_extra(state, cache_key)
    if (!is.null(cached)) {
      cluster_val <- cached
    } else {
      cluster_val <- prod(vapply(cluster_exprs, function(ce) {
        .eval_expr_survival_cond(
          ce, t, prep, component,
          forced_complete = integer(0),
          forced_survive = integer(0),
          state = state
        )
      }, numeric(1)))
      if (!is.finite(cluster_val) || cluster_val < 0) cluster_val <- 0.0
      .eval_state_set_extra(state, cache_key, cluster_val)
    }
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
    # Survival of guard = 1 - ∫_0^t f_guard(u) du (use direct guard density)
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .guard_density_at(expr, uu, prep, component, integer(0), integer(0), state = NULL)
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
                                  state = NULL) {
  compiled <- .expr_lookup_compiled(expr, prep)
  if (is.null(compiled)) {
    return(0.0)
  }
  .node_density(
    compiled, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = state
  )
}

# ==============================================================================
