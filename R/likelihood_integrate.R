# PART 5: Outcome Likelihood
# ==============================================================================

# Integrand for scenario-conditioned probability (density already encodes competitors)
.integrand_outcome_density <- function(t, expr, prep, component, competitor_exprs,
                                       state = NULL, trial_rows = NULL,
                                       trial_overrides = NULL) {
  state <- state %||% .eval_state_create()
  times <- as.numeric(t)
  if (!length(times)) return(numeric(0))
  res <- .scenario_density_with_competitors(expr, times, prep, component,
                                            competitor_exprs = competitor_exprs,
                                            state = state,
                                            trial_rows = trial_rows,
                                            trial_overrides = trial_overrides)
  as.numeric(res)
}

# Scenario-conditioned probability of an expression over [0, upper_limit].
# competitor_exprs retained for backward compatibility but ignored.
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf,
                                           competitor_exprs = NULL,
                                           state = NULL,
                                           trial_rows = NULL,
                                           trial_overrides = NULL) {
  if (!is.finite(upper_limit)) upper_limit <- Inf
  state <- state %||% .eval_state_create()
  guard_cache_key <- NULL
  component_key <- .cache_component_key(component)
  competitor_exprs <- competitor_exprs %||% list()
  expr_kind <- expr[['kind']]
  if (identical(expr_kind, "guard")) {
    guard_sig <- .expr_signature(expr)
    competitor_sig <- if (length(competitor_exprs) == 0L) {
      "none"
    } else {
      paste(sort(vapply(competitor_exprs, .expr_signature, character(1))), collapse = ",")
    }
    guard_cache_key <- .guard_integral_cache_key(guard_sig, component_key, upper_limit, competitor_sig)
    cached_guard <- .guard_integral_fetch(prep, guard_cache_key)
    if (!is.null(cached_guard)) {
      return(cached_guard)
    }
  }
  native_ctx <- NULL
  compiled <- .expr_lookup_compiled(expr, prep)
  competitor_exprs <- competitor_exprs %||% list()
  comp_ids <- integer(0)
  if (length(competitor_exprs) > 0L) {
    comp_nodes <- lapply(competitor_exprs, .expr_lookup_compiled, prep = prep)
    if (!any(vapply(comp_nodes, is.null, logical(1)))) {
      comp_ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
      if (any(is.na(comp_ids))) comp_ids <- NULL
    } else {
      comp_ids <- NULL
    }
  }
  trial_rows_df <- NULL
  if (!is.null(trial_rows) && inherits(trial_rows, "data.frame") && nrow(trial_rows) > 0L) {
    trial_rows_df <- as.data.frame(trial_rows)
  }
  use_native <- is.finite(upper_limit) && !is.null(compiled) && !is.null(comp_ids)
  if (use_native) {
    native_ctx <- .prep_native_context(prep)
    if (!is.null(trial_overrides) && inherits(trial_overrides, "externalptr")) {
      native_fn <- .lik_native_fn("native_outcome_probability_overrides_cpp")
      val <- native_fn(
        native_ctx,
        as.integer(compiled$id),
        as.numeric(upper_limit),
        component,
        integer(0),
        integer(0),
        as.integer(comp_ids),
        .integrate_rel_tol(),
        .integrate_abs_tol(),
        12L,
        trial_overrides
      )
    } else if (!is.null(trial_rows_df)) {
      native_fn <- .lik_native_fn("native_outcome_probability_params_cpp")
      val <- native_fn(
        native_ctx,
        as.integer(compiled$id),
        as.numeric(upper_limit),
        component,
        integer(0),
        integer(0),
        as.integer(comp_ids),
        .integrate_rel_tol(),
        .integrate_abs_tol(),
        12L,
        trial_rows_df
      )
    } else {
      native_fn <- .lik_native_fn("native_outcome_probability_cpp")
      val <- native_fn(
        native_ctx,
        as.integer(compiled$id),
        as.numeric(upper_limit),
        component,
        integer(0),
        integer(0),
        as.integer(comp_ids),
        .integrate_rel_tol(),
        .integrate_abs_tol(),
        12L
      )
    }
  } else {
    integrand <- function(t) {
      .integrand_outcome_density(
        t,
        expr = expr,
        prep = prep,
        component = component,
        competitor_exprs = competitor_exprs,
        state = state,
        trial_rows = trial_rows_df,
        trial_overrides = trial_overrides
      )
    }
    val <- .native_integrate(
      integrand,
      lower = 0,
      upper = upper_limit,
      rel.tol = .integrate_rel_tol(),
      abs.tol = .integrate_abs_tol(),
      max.depth = 12L
    )
  }
  if (!is.finite(val)) val <- 0.0
  if (val < 0) val <- 0.0
  if (val > 1) val <- 1.0
  if (!is.null(guard_cache_key)) {
    .guard_integral_store(prep, guard_cache_key, val)
  }
  val
}

.expr_sources <- function(expr, prep) {
  gather_labels <- function(ex) {
    if (is.null(ex) || is.null(ex[['kind']])) return(character(0))
    kind <- ex[['kind']]
    if (identical(kind, "event")) {
      source_id <- ex[['source']]
      if (!is.null(prep[["pools"]][[source_id]])) {
        members <- prep[["pools"]][[source_id]][['members']] %||% character(0)
        unique(unlist(lapply(members, function(mid) {
          if (!is.null(prep[["accumulators"]][[mid]])) return(mid)
          if (!is.null(prep[["pools"]][[mid]])) {
            return(gather_labels(list(kind = "event", source = mid)))
          }
          character(0)
        })))
      } else {
        source_id
      }
    } else if (identical(kind, "and") || identical(kind, "or")) {
      unique(unlist(lapply(ex[['args']], gather_labels)))
    } else if (identical(kind, "guard")) {
      ref_src <- gather_labels(ex[['reference']])
      blk_src <- gather_labels(ex[['blocker']])
      unless_src <- unique(unlist(lapply(ex[['unless']] %||% list(), gather_labels)))
      unique(c(ref_src, blk_src, unless_src))
    } else if (identical(kind, "not")) {
      gather_labels(ex[['arg']])
    } else {
      character(0)
    }
  }
  .labels_to_ids(prep, gather_labels(expr))
}

.build_competitor_exprs <- function(prep, target_label, target_expr) {
  outcome_defs <- prep[["outcomes"]]
  comp_labels <- setdiff(names(outcome_defs), target_label)
  if (length(comp_labels) == 0) return(list())
  target_opts <- outcome_defs[[target_label]][['options']] %||% list()
  target_map <- target_opts[['map_outcome_to']] %||% target_label
  comp_labels <- Filter(function(lbl) {
    comp_opts <- outcome_defs[[lbl]][['options']] %||% list()
    if (!is.null(comp_opts[['alias_of']])) return(FALSE)
    expr_lbl <- outcome_defs[[lbl]][['expr']]
    if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
      src <- expr_lbl[['source']]
      if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
    }
    comp_target <- comp_opts[['map_outcome_to']] %||% lbl
    if (identical(as.character(comp_target), as.character(target_map))) return(FALSE)
    TRUE
  }, comp_labels)
  if (length(comp_labels) == 0) return(list())
  comps <- lapply(comp_labels, function(lbl) outcome_defs[[lbl]][['expr']])
  target_sources <- .expr_sources(target_expr, prep)
  if (length(comps) > 0) {
    comps <- lapply(comps, function(comp_expr) {
      if (!is.null(comp_expr[['kind']]) && identical(comp_expr[['kind']], "guard")) {
        blocker_sources <- .expr_sources(comp_expr[['blocker']], prep)
        if (length(target_sources) > 0 && all(target_sources %in% blocker_sources)) {
          return(comp_expr[['reference']])
        }
      }
      comp_expr
    })
  }
  if (!is.null(target_expr[['kind']]) && identical(target_expr[['kind']], "guard") && length(comps) > 0) {
    blocker_sources <- .expr_sources(target_expr[['blocker']], prep)
    keep_idx <- vapply(comps, function(comp_expr) {
      comp_sources <- .expr_sources(comp_expr, prep)
      length(comp_sources) == 0 || any(!comp_sources %in% blocker_sources)
    }, logical(1))
    comps <- comps[keep_idx]
    comp_labels <- comp_labels[keep_idx]
  }
  if (length(comps) > 0 && isTRUE(getOption("uuber.filter_guard_competitors", TRUE))) {
    primary_guards <- .collect_guards(target_expr)
    if (length(primary_guards) > 0) {
      guard_keep <- vapply(seq_along(comps), function(i) {
        comp_expr <- comps[[i]]
        !any(vapply(primary_guards, function(pg) .guards_conflict(pg, comp_expr), logical(1)))
      }, logical(1))
      comps <- comps[guard_keep]
      comp_labels <- comp_labels[guard_keep]
    }
  }
  comps
}

.prepare_competitor_map <- function(prep) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  if (length(outcome_defs) == 0) return(list())
  comps <- lapply(names(outcome_defs), function(lbl) {
    exprs <- .build_competitor_exprs(prep, lbl, outcome_defs[[lbl]][['expr']])
    attr(exprs, ".lik_outcome_label") <- lbl
    exprs
  })
  names(comps) <- names(outcome_defs)
  .attach_competitor_fast_clusters(prep, comps)
}

.attach_competitor_fast_clusters <- function(prep, competitor_map) {
  if (length(competitor_map) == 0) return(competitor_map)
  cluster_counter <- 0L
  for (lbl in names(competitor_map)) {
    expr_list <- competitor_map[[lbl]] %||% list()
    if (!is.list(expr_list)) expr_list <- list(expr_list)
    if (length(expr_list) == 0) {
      competitor_map[[lbl]] <- expr_list
      next
    }
    clusters <- .cluster_competitors(expr_list, prep)
    if (length(clusters) == 0) {
      competitor_map[[lbl]] <- expr_list
      next
    }
    for (cl in clusters) {
      nodes <- lapply(cl$exprs, function(expr) .expr_lookup_compiled(expr, prep))
      if (length(nodes) == 0) next
      has_fast <- all(vapply(nodes, function(node) is.function(node$surv_fast_fn) && !isTRUE(node$scenario_sensitive), logical(1)))
      if (!has_fast) next
      cluster_counter <- cluster_counter + 1L
      key <- paste0("fast_comp_cluster_", cluster_counter)
      fast_fn <- .make_competitor_cluster_fast_closure(nodes)
      indices <- cl$indices %||% seq_along(cl$exprs)
      for (idx in indices) {
        if (is.na(idx) || idx < 1L || idx > length(expr_list)) next
        expr_obj <- expr_list[[idx]]
        attr(expr_obj, ".lik_cluster_fast_fn") <- fast_fn
        attr(expr_obj, ".lik_cluster_fast_key") <- key
        attr(expr_obj, ".lik_outcome_label") <- lbl
        expr_list[[idx]] <- expr_obj
      }
    }
    competitor_map[[lbl]] <- expr_list
  }
  competitor_map
}

.make_competitor_cluster_fast_closure <- function(nodes) {
  function(t, component) {
    t_vals <- as.numeric(t)
    if (length(t_vals) == 0L) return(numeric(0))
    vapply(t_vals, function(tt) {
      prod(vapply(nodes, function(node) {
        vals <- node$surv_fast_fn(tt, component)
        if (length(vals) == 0L) 1.0 else vals[[1]]
      }, numeric(1)))
    }, numeric(1))
  }
}

# Detect shared-gate pattern (A & C) vs (B & C) across outcome definitions
.detect_shared_gate <- function(outcome_defs, target_label) {
  out_def <- outcome_defs[[target_label]]
  if (is.null(out_def) || is.null(out_def[['expr']])) return(NULL)
  expr <- out_def[['expr']]
  if (is.null(expr[['kind']]) || !identical(expr[['kind']], "and") || length(expr[['args']]) != 2) return(NULL)
  get_evt_id <- function(node) {
    if (!is.null(node) && !is.null(node[['kind']]) && identical(node[['kind']], "event")) return(node[['source']])
    NULL
  }
  x_id <- get_evt_id(expr[['args']][[1]])
  c_id <- get_evt_id(expr[['args']][[2]])
  if (is.null(x_id) || is.null(c_id)) {
    x_id <- get_evt_id(expr[['args']][[2]])
    c_id <- get_evt_id(expr[['args']][[1]])
  }
  if (is.null(x_id) || is.null(c_id)) return(NULL)
  # Find another outcome with same gate c_id
  for (lbl in names(outcome_defs)) {
    if (identical(lbl, target_label)) next
    other_def <- outcome_defs[[lbl]]
    if (!is.null(other_def[['options']][['alias_of']])) next
    other <- other_def[['expr']]
    if (is.null(other) || is.null(other[['kind']]) || !identical(other[['kind']], "and") || length(other[['args']]) != 2) next
    a1 <- get_evt_id(other[['args']][[1]])
    a2 <- get_evt_id(other[['args']][[2]])
    ids_other <- c(a1, a2)
    if (c_id %in% ids_other) {
      y_id <- setdiff(ids_other, c_id)
      if (length(y_id) == 1 && !identical(y_id[[1]], x_id)) {
        return(list(x_id = x_id, y_id = y_id[[1]], c_id = c_id))
      }
    }
  }
  NULL
}

# Find a single shared-gate pair among non-alias outcomes, if present
.find_shared_gate_pair <- function(outcome_defs) {
  labels <- names(outcome_defs)
  core_labels <- Filter(function(lbl) is.null(outcome_defs[[lbl]][['options']][['alias_of']]), labels)
  # collect all ANDs with two events
  events <- lapply(core_labels, function(lbl) {
    ex <- outcome_defs[[lbl]][['expr']]
    if (!is.null(ex[['kind']]) && identical(ex[['kind']], "and") && length(ex[['args']]) == 2) {
      a1 <- if (!is.null(ex[['args']][[1]][['kind']]) && ex[['args']][[1]][['kind']] == "event") ex[['args']][[1]][['source']] else NULL
      a2 <- if (!is.null(ex[['args']][[2]][['kind']]) && ex[['args']][[2]][['kind']] == "event") ex[['args']][[2]][['source']] else NULL
      if (!is.null(a1) && !is.null(a2)) return(list(label = lbl, ids = c(a1, a2)))
    }
    NULL
  })
  events <- events[!vapply(events, is.null, logical(1))]
  if (length(events) < 2) return(NULL)
  # try to find two with exactly one shared id (the gate)
  for (i in seq_along(events)) {
    for (j in seq_along(events)) {
      if (j <= i) next
      ids_i <- events[[i]][['ids']]
      ids_j <- events[[j]][['ids']]
      shared <- intersect(ids_i, ids_j)
      if (length(shared) == 1) {
        c_id <- shared[[1]]
        x_id <- setdiff(ids_i, c_id)[[1]]
        y_id <- setdiff(ids_j, c_id)[[1]]
        return(list(
          label_x = events[[i]][['label']],
          label_y = events[[j]][['label']],
          x_id = x_id,
          y_id = y_id,
          c_id = c_id
        ))
      }
    }
  }
  NULL
}

# Compute probabilities for a shared-gate pair using the general likelihood machinery
.shared_gate_pair_probs <- function(prep, component, pair_info) {
  labels <- c(pair_info$label_x, pair_info$label_y)
  vapply(labels, function(lbl) {
    .outcome_likelihood(lbl, NA_real_, prep, component)
  }, numeric(1))
}

.shared_gate_fast_context <- function(prep, component, info, state) {
  if (is.null(info)) return(NULL)
  make_event_expr <- function(event_id) list(kind = "event", source = event_id)
  fetch_node <- function(event_id) {
    expr <- make_event_expr(event_id)
    .expr_lookup_compiled(expr, prep)
  }
  nodes <- list(
    x = fetch_node(info$x_id),
    y = fetch_node(info$y_id),
    c = fetch_node(info$c_id)
  )
  if (any(vapply(nodes, is.null, logical(1)))) return(NULL)
  if (any(vapply(nodes, function(node) {
    !is.function(node$density_fast_fn) || !is.function(node$surv_fast_fn)
  }, logical(1)))) {
    return(NULL)
  }
  wrap_density <- function(node) {
    fast <- node$density_fast_fn
    function(t_vec) {
      tv <- as.numeric(t_vec)
      if (length(tv) == 0L) return(numeric(0))
      vals <- fast(tv, component)
      out <- as.numeric(vals)
      out[!is.finite(out)] <- 0.0
      out
    }
  }
  wrap_survival <- function(node) {
    fast <- node$surv_fast_fn
    function(t_vec) {
      tv <- as.numeric(t_vec)
      if (length(tv) == 0L) return(numeric(0))
      vals <- fast(tv, component)
      out <- as.numeric(vals)
      out[!is.finite(out)] <- 0.0
      out[out < 0] <- 0.0
      out[out > 1] <- 1.0
      out
    }
  }
  fX <- wrap_density(nodes$x)
  fY <- wrap_density(nodes$y)
  fC <- wrap_density(nodes$c)
  SX <- wrap_survival(nodes$x)
  SY <- wrap_survival(nodes$y)
  SC <- wrap_survival(nodes$c)
  FX <- function(t_vec) {
    1.0 - SX(t_vec)
  }
  FY <- function(t_vec) {
    1.0 - SY(t_vec)
  }
  FC <- function(t_vec) {
    1.0 - SC(t_vec)
  }
  comp_key <- .eval_state_component_key(component)
  base_tag <- paste("shared_palloc", info$x_id, info$y_id, info$c_id, comp_key, sep = "|")
  shared_palloc <- function(limit) {
    limits <- as.numeric(limit)
    if (length(limits) == 0L) return(numeric(0))
    vapply(limits, function(lim) {
      if (!is.finite(lim) || lim <= 0) return(0.0)
      limit_key <- paste(base_tag, .eval_state_time_key(lim), sep = "|")
      cached <- .eval_state_get_extra(state, limit_key)
      if (!is.null(cached)) return(as.numeric(cached))
      denom_val <- as.numeric(FX(lim) * FY(lim))
      denom <- if (length(denom_val) == 0L) 0.0 else denom_val[[1]]
      if (!is.finite(denom) || denom <= 0) {
        out <- 0.0
      } else {
        integral_val <- tryCatch(
          stats::integrate(
            function(u) fX(u) * FY(u),
            lower = 0,
            upper = lim,
            rel.tol = .integrate_rel_tol(),
            abs.tol = .integrate_abs_tol(),
            stop.on.error = FALSE
          )[["value"]],
          error = function(e) 0.0
        )
        out <- 1.0 - (as.numeric(integral_val) / denom)
        if (!is.finite(out)) out <- 0.0
        if (out < 0) out <- 0.0
        if (out > 1) out <- 1.0
      }
      .eval_state_set_extra(state, limit_key, out)
      out
    }, numeric(1))
  }
  list(
    fX = fX,
    fY = fY,
    fC = fC,
    FX = FX,
    FY = FY,
    FC = FC,
    SY = SY,
    palloc = shared_palloc
  )
}

# Compute likelihood for a single trial/outcome
.outcome_likelihood <- function(outcome_label, rt, prep, component,
                                trial_rows = NULL,
                                trial_overrides = NULL) {
  outcome_defs <- prep[["outcomes"]]
  competitor_map <- .prep_competitors(prep) %||% list()
  state <- .eval_state_create()
  trial_rows_df <- NULL
  if (!is.null(trial_rows) && inherits(trial_rows, "data.frame") && nrow(trial_rows) > 0L) {
    trial_rows_df <- as.data.frame(trial_rows)
  }
  
  # FIRST: Check for GUESS outcomes (can be both defined and from guess policy)
  # GUESS outcomes arise from component-level guess policies
  if (identical(outcome_label, "GUESS")) {
    # Get the guess policy for this component
    guess_policy <- .get_component_attr(prep, component, "guess")
    
    if (!is.null(guess_policy) && !is.null(guess_policy[['weights']])) {
      # The guess policy has:
      # - outcome: the label for guess outcomes
      # - weights: named vector of keep probabilities
      
      guess_weights <- guess_policy[['weights']]
      
      # Probability of GUESS outcome is:
      # Sum over outcomes: P(outcome fires) * P(guess | outcome)
      # P(guess | outcome) = 1 - weight[outcome]
      
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep[["default_deadline"]]
      
      if (is.finite(deadline)) {
        if (is.na(rt)) {
          total_prob <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            expr <- outcome_defs[[label]][['expr']]
            comp_exprs_guess <- competitor_map[[label]] %||% list()
            prob_outcome <- .integrate_outcome_probability(expr, prep, component, deadline,
                                                           competitor_exprs = comp_exprs_guess,
                                                           state = state,
                                                           trial_rows = trial_rows_df,
                                                           trial_overrides = trial_overrides)
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_prob <- total_prob + prob_outcome * guess_prob
          }
          return(total_prob)
        } else {
          total_dens <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            expr <- outcome_defs[[label]][['expr']]
            comp_exprs_guess <- competitor_map[[label]] %||% list()
            dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                         competitor_exprs = comp_exprs_guess,
                                                         state = state,
                                                         trial_rows = trial_rows_df,
                                                         trial_overrides = trial_overrides)
            if (dens_r == 0) next
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_dens <- total_dens + guess_prob * dens_r
          }
          return(total_dens)
        }
      }
    }
    # If no guess policy or deadline, fall through to regular handling
  }
  
  # Find the outcome
  outcome_def <- outcome_defs[[outcome_label]]
  if (is.null(outcome_def)) {
    # Check if this is a special outcome
    if (is.na(outcome_label) || identical(outcome_label, "NA")) {
      # NA outcome - this could be from map_outcome_to = NA
      # Need to find which outcome maps to NA
      for (label in names(outcome_defs)) {
        out_def <- outcome_defs[[label]]
        if (!is.null(out_def[['options']][['map_outcome_to']])) {
          mapped <- out_def[['options']][['map_outcome_to']]
          if (is.na(mapped) || identical(mapped, "NA")) {
            # Found the outcome that maps to NA
            expr <- out_def[['expr']]
            if (is.na(rt)) {
              # Integrate race probability for the mapped source outcome
              deadline <- .get_component_attr(prep, component, "deadline")
              deadline <- deadline %||% prep[["default_deadline"]]
              comp_labels_map <- setdiff(names(outcome_defs), label)
              if (length(comp_labels_map) > 0) {
                comp_labels_map <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels_map)
              }
              comp_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              return(.integrate_outcome_probability(expr, prep, component, deadline,
                                                    competitor_exprs = comp_exprs_map,
                                                    state = state,
                                                    trial_rows = trial_rows_df,
                                                    trial_overrides = trial_overrides))
            } else {
              # Race density at rt for the mapped source outcome
              comp_labels_map <- setdiff(names(outcome_defs), label)
              if (length(comp_labels_map) > 0) {
                comp_labels_map <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels_map)
              }
              comp_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                           competitor_exprs = comp_exprs_map,
                                                           state = state,
                                                           trial_rows = trial_rows_df,
                                                           trial_overrides = trial_overrides)
              return(as.numeric(dens_r))
            }
          }
        }
      }
      # Couldn't find the NA mapping
      warning(sprintf("Unknown outcome '%s' and no NA mapping found", outcome_label))
      return(0.0)
    }
    
    # Check for special deadline outcomes
    if (identical(outcome_label, "NR_DEADLINE")) {
      # Deadline reached - compute probability that no outcome occurred by deadline
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep[["default_deadline"]]
      
      if (!is.finite(deadline)) return(0.0)
      
      # P(no outcome by deadline) = Π_i S_i(deadline)
      # where i ranges over all possible outcomes
      prob_none <- 1.0
      for (label in names(outcome_defs)) {
        out_def <- outcome_defs[[label]]
        expr <- out_def[['expr']]
        # Survival of this outcome at deadline
        surv <- .eval_expr_survival(expr, deadline, prep, component)
        prob_none <- prob_none * surv
      }
      return(prob_none)
    }
    
    warning(sprintf("Unknown outcome '%s'", outcome_label))
    return(0.0)
  }
  
  expr <- outcome_def[['expr']]
  options <- outcome_def[['options']] %||% list()
  
  # Check if outcome is mapped to something else
  if (!is.null(options[['map_outcome_to']])) {
    # This outcome is mapped, but we shouldn't reach here if data has the mapped label
    # Just proceed with the expression
  }
  
  # If this outcome has an outcome-level guess policy, it never appears as-is.
  # Its mass is redistributed to target labels via donors handled below.
  if (!is.null(options[['guess']])) {
    return(0.0)
  }
  
  competitor_exprs <- competitor_map[[outcome_label]] %||% list()
  
  # Build donor (guess) contributions: outcomes that guess into this label
  donors <- list()
  for (dlbl in setdiff(names(outcome_defs), outcome_label)) {
    dopt <- outcome_defs[[dlbl]][['options']][['guess']] %||% NULL
    if (is.null(dopt)) next
    glabels <- dopt[['labels']] %||% NULL
    gweights <- dopt[['weights']] %||% NULL
    if (is.null(glabels) || is.null(gweights)) next
    idx <- which(as.character(glabels) == outcome_label)
    if (length(idx) == 0) next
    w <- as.numeric(gweights[[idx[[1]]]])
    rt_policy <- dopt[['rt_policy']] %||% "keep"
    donors[[length(donors) + 1L]] <- list(
      expr = outcome_defs[[dlbl]][['expr']],
      weight = w,
      rt_policy = rt_policy,
      competitors = competitor_map[[dlbl]] %||% list()
    )
  }
  
  # Detect shared-gate AND pattern: (X & C) vs (Y & C)
  shared_gate_info <- .detect_shared_gate(outcome_defs, outcome_label)
  
  # Helper closures to get density/CDF for pools or accumulators
  get_f <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti < 0) return(0.0)
        .event_density_at(prep, id, component, ti,
                          forced_complete = integer(0),
                          forced_survive = integer(0))
      }, numeric(1))
    }
  }
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        .event_survival_at(prep, id, component, ti,
                           forced_complete = integer(0),
                           forced_survive = integer(0))
      }, numeric(1))
    }
  }
  get_F <- function(id) {
    Sfun <- get_S(id)
    function(tt) { 1.0 - Sfun(tt) }
  }
  
  shared_palloc <- NULL
  shared_gate_ctx <- NULL
  if (!is.null(shared_gate_info)) {
    shared_gate_ctx <- .shared_gate_fast_context(prep, component, shared_gate_info, state)
  }
  if (!is.null(shared_gate_ctx)) {
    fX <- shared_gate_ctx$fX
    fY <- shared_gate_ctx$fY
    fC <- shared_gate_ctx$fC
    FX <- shared_gate_ctx$FX
    FY <- shared_gate_ctx$FY
    FC <- shared_gate_ctx$FC
    SY <- shared_gate_ctx$SY
    shared_palloc <- shared_gate_ctx$palloc
  } else if (!is.null(shared_gate_info)) {
    fX <- get_f(shared_gate_info[['x_id']])
    fY <- get_f(shared_gate_info[['y_id']])
    fC <- get_f(shared_gate_info[['c_id']])
    FX <- get_F(shared_gate_info[['x_id']])
    FY <- get_F(shared_gate_info[['y_id']])
    FC <- get_F(shared_gate_info[['c_id']])
    SY <- get_S(shared_gate_info[['y_id']])
    shared_palloc <- function(limit) {
      if (!is.finite(limit) || limit <= 0) return(0.0)
      denom <- FX(limit) * FY(limit)
      if (!is.finite(denom) || denom <= 0) return(0.0)
      integral_val <- tryCatch(
        stats::integrate(
          function(u) fX(u) * FY(u),
          lower = 0,
          upper = limit,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]],
        error = function(e) 0.0
      )
      out <- 1.0 - (as.numeric(integral_val) / as.numeric(denom))
      if (!is.finite(out)) out <- 0.0
      max(0.0, min(1.0, out))
    }
  }
  
  # Handle NA/infinite RT by integration of race density
  if (is.na(rt) || !is.finite(rt) || rt < 0) {
    deadline <- .get_component_attr(prep, component, "deadline")
    deadline <- deadline %||% prep[["default_deadline"]]
    if (!is.null(shared_gate_info)) {
      integrand <- function(t) {
        vapply(t, function(tt) {
          term1 <- fX(tt) * FC(tt) * SY(tt)
          denom <- FX(tt) * FY(tt)
          termCother <- fC(tt) * FX(tt) * SY(tt)
          term2 <- if (!is.null(shared_palloc) && is.finite(denom) && denom > 0) {
            fC(tt) * FX(tt) * FY(tt) * shared_palloc(tt)
          } else {
            0.0
          }
          term1 + termCother + term2
        }, numeric(1))
      }
      keep_prob <- 1.0
      gp_comp <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp) && !is.null(gp_comp[['weights']])) {
        kp <- gp_comp[['weights']][[outcome_label]] %||% gp_comp[['weights']][[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob <- as.numeric(kp)
      }
      upper_lim <- if (is.finite(deadline)) deadline else Inf
      res <- tryCatch(
        stats::integrate(
          integrand,
          lower = 0,
          upper = upper_lim,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]],
        error = function(e) 0.0
      )
      res <- as.numeric(res)
      if (!is.finite(res) || res < 0) res <- 0.0
      return(keep_prob * res)
    }
    integrand <- function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        base <- .scenario_density_with_competitors(expr, tt, prep, component,
                                                   competitor_exprs = competitor_exprs,
                                                   state = state,
                                                   trial_rows = trial_rows_df,
                                                   trial_overrides = trial_overrides)
        if (base == 0) return(0.0)
        add <- 0.0
        if (length(donors) > 0) {
          for (d in donors) {
            dens_d <- .scenario_density_with_competitors(d[['expr']], tt, prep, component,
                                                         competitor_exprs = d[['competitors']] %||% list(),
                                                         state = state,
                                                         trial_rows = trial_rows_df,
                                                         trial_overrides = trial_overrides)
            if (dens_d == 0) next
            add <- add + d[['weight']] * dens_d
          }
        }
        base + add
      }, numeric(1))
    }
    upper_lim <- if (is.finite(deadline)) deadline else Inf
    res <- tryCatch(
      stats::integrate(
        integrand,
        lower = 0,
        upper = upper_lim,
        rel.tol = .integrate_rel_tol(),
        abs.tol = .integrate_abs_tol(),
        stop.on.error = FALSE
      )[["value"]],
      error = function(e) 0.0
    )
    res <- as.numeric(res)
    if (!is.finite(res) || res < 0) res <- 0.0
    return(res)
  }
  
  # Race density at observed RT
  if (!is.null(shared_gate_info)) {
    term1 <- fX(rt) * FC(rt) * SY(rt)
    denom <- FX(rt) * FY(rt)
    term2 <- 0.0
    termCother <- fC(rt) * FX(rt) * SY(rt)
    if (!is.null(shared_palloc) && is.finite(denom) && denom > 0) {
      palloc <- shared_palloc(rt)
      term2 <- fC(rt) * FX(rt) * FY(rt) * palloc
    }
    base_val <- term1 + termCother + term2
  } else {
      dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                   competitor_exprs = competitor_exprs,
                                                   state = state,
                                                   trial_rows = trial_rows_df,
                                                   trial_overrides = trial_overrides)
    if (dens_r == 0) return(0.0)
    base_val <- dens_r
    # Include donor mass (e.g., timeout guesses) that keeps the observed RT
    if (length(donors) > 0) {
      donor_add <- 0.0
      for (d in donors) {
        if (identical(d[['rt_policy']], "na")) next
        dens_d <- .scenario_density_with_competitors(d[['expr']], rt, prep, component,
                                                     competitor_exprs = d[['competitors']] %||% list(),
                                                     state = state,
                                                     trial_rows = trial_rows_df,
                                                     trial_overrides = trial_overrides)
        if (dens_d == 0) next
        donor_add <- donor_add + d[['weight']] * dens_d
      }
      base_val <- base_val + donor_add
    }
  }
  # Apply component-level keep probability to base density for observed labels
  keep_prob_rt <- 1.0
  gp_comp_rt <- .get_component_attr(prep, component, "guess")
  if (!is.null(gp_comp_rt) && !is.null(gp_comp_rt[['weights']])) {
    kp <- gp_comp_rt[['weights']][[outcome_label]] %||% gp_comp_rt[['weights']][[normalize_label(outcome_label)]] %||% NULL
    if (!is.null(kp)) keep_prob_rt <- as.numeric(kp)
  }
  base_val <- base_val * keep_prob_rt
  # Component-level guess mass is handled via the explicit GUESS outcome
  base_val
}

# ==============================================================================
# PART 6: Top-level Likelihood Function
# ==============================================================================

#' Compute log-likelihood for data given a model
#'
#' @param model Model specification from new_API.R (race_spec or race_model_spec)
#' @param data Data frame with columns: outcome, rt, (optionally: component)
#' @return Log-likelihood value with per-trial attributes
compute_loglik <- function(model, data) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  
  if (!"outcome" %in% names(data) || !"rt" %in% names(data)) {
    stop("data must contain 'outcome' and 'rt' columns")
  }
  
  prep <- .prepare_model_for_likelihood(model)
  
  comp_info <- prep[["components"]]
  comp_ids <- comp_info[['ids']] %||% list()
  weights <- comp_info[['weights']] %||% numeric(0)
  has_weight_param <- isTRUE(any(comp_info[['has_weight_param']] %||% FALSE))
  is_mixture <- length(comp_ids) > 1
  default_component <- if (length(comp_ids) > 0) comp_ids[[1]] else "__default__"
  
  n_rows <- nrow(data)
  if (n_rows == 0) {
    total_ll <- 0.0
    attr(total_ll, "contributions") <- numeric(0)
    attr(total_ll, "log_contributions") <- numeric(0)
    attr(total_ll, "keys") <- character(0)
    return(total_ll)
  }
  
  has_component_col <- "component" %in% names(data)
  raw_keys <- character(n_rows)
  lik_values <- numeric(n_rows)

  for (i in seq_len(n_rows)) {
    outcome <- as.character(data[['outcome']][[i]])
    rt_val <- as.numeric(data[['rt']][[i]])
    rt_key <- .eval_state_time_key(rt_val)
    component_key <- default_component
    if (is_mixture) {
      if (has_weight_param || !has_component_col) {
        component_key <- "__mixture__"
      } else {
        raw_component <- data[['component']][[i]]
        if (is.na(raw_component)) {
          component_key <- default_component
        } else {
          component_chr <- as.character(raw_component)
          component_key <- if (nzchar(component_chr)) component_chr else default_component
        }
      }
    }
    cache_key <- .likelihood_outcome_cache_key(component_key, outcome, rt_val)
    res <- .likelihood_outcome_cached(prep, cache_key, function() {
      if (is_mixture && identical(component_key, "__mixture__")) {
        base_weights <- if (length(comp_ids) > 0 && length(weights) == length(comp_ids)) {
          weights
        } else if (length(comp_ids) > 0) {
          rep(1 / length(comp_ids), length(comp_ids))
        } else {
          numeric(0)
        }
        total_mix <- 0.0
        if (length(comp_ids) > 0) {
          for (j in seq_along(comp_ids)) {
            sub_key <- .likelihood_outcome_cache_key(comp_ids[[j]], outcome, rt_val)
            sub_res <- .likelihood_outcome_cached(prep, sub_key, function() {
              .outcome_likelihood(outcome, rt_val, prep, comp_ids[[j]])
            })
            prep <<- sub_res$prep
            total_mix <- total_mix + base_weights[[j]] * as.numeric(sub_res$value)
          }
        }
        total_mix
      } else {
        .outcome_likelihood(outcome, rt_val, prep, component_key)
      }
    })
    prep <- res$prep
    lik <- res$value
    
    raw_keys[[i]] <- cache_key
    lik_values[[i]] <- lik
  }
  log_lik_values <- ifelse(is.finite(lik_values) & lik_values > 0, log(lik_values), -Inf)
  
  total_ll <- sum(log_lik_values)
  attr(total_ll, "contributions") <- lik_values
  attr(total_ll, "log_contributions") <- log_lik_values
  attr(total_ll, "keys") <- raw_keys
  
  total_ll
}

# ==============================================================================
# PART 7b: Response probability helpers (analytic, no simulation)
# ==============================================================================

#' Compute probability of an outcome label under the model
#' If the model has multiple components and no specific component is supplied,
#' returns the mixture-weighted probability.
#'
#' @param model Race model spec (race_spec or race_model_spec)
#' @param outcome_label Outcome label to evaluate (character)
#' @param component Optional component id (character) for mixture models
#' @return Probability (numeric in [0,1])
response_probability <- function(model, outcome_label, component = NULL) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  comp_ids <- prep[["components"]][['ids']]
  weights <- prep[["components"]][['weights']]
  
  # Helper for a single component
  .prob_for_component <- function(comp_id) {
    # Shared-gate pair fast-path (can be disabled globally)
    use_fastpath <- getOption("uuber.shared_gate_fastpath", default = TRUE)
    if (isTRUE(use_fastpath)) {
      pair <- .find_shared_gate_pair(prep[["outcomes"]])
      if (!is.null(pair) && outcome_label %in% c(pair[['label_x']], pair[['label_y']])) {
        vals <- .shared_gate_pair_probs(prep, comp_id, pair)
        if (identical(outcome_label, pair[['label_x']])) return(vals[[1]]) else return(vals[[2]])
      }
    }
    # Alias-of outcomes: sum of referenced labels
    out_def <- prep[["outcomes"]][[outcome_label]]
    if (!is.null(out_def[['options']][['alias_of']])) {
      refs <- out_def[['options']][['alias_of']]
      refs <- as.character(refs)
      vals <- vapply(refs, function(lbl) as.numeric(.outcome_likelihood(lbl, NA_real_, prep, comp_id)), numeric(1))
      return(sum(vals))
    }
    base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, comp_id))
    # Apply component-level keep for non-special labels (complement goes to GUESS)
    if (!identical(outcome_label, "GUESS")) {
      gp <- .get_component_attr(prep, comp_id, "guess")
      if (!is.null(gp) && !is.null(gp[['weights']])) {
        keep <- gp[['weights']][[outcome_label]] %||% gp[['weights']][[normalize_label(outcome_label)]] %||% 1.0
        base <- base * as.numeric(keep)
      }
    }
    base
  }
  
  if (!is.null(component)) {
    return(.prob_for_component(component))
  }
  
  if (length(comp_ids) <= 1L) {
    comp <- if (length(comp_ids) == 1L) comp_ids[[1]] else "__default__"
    return(.prob_for_component(comp))
  }
  
  # Mixture: average with weights
  vals <- vapply(seq_along(comp_ids), function(i) .prob_for_component(comp_ids[[i]]), numeric(1))
  sum(weights * vals)
}

#' Compute probabilities for all defined outcomes in the model
#' If component is NULL, returns mixture-weighted probabilities.
#'
#' @param model Race model spec
#' @param component Optional component id
#' @return Named numeric vector of probabilities
response_probabilities <- function(model, component = NULL) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  labels <- names(prep[["outcomes"]])
  vals <- vapply(labels, function(lbl) response_probability(model, lbl, component = component), numeric(1))
  stats::setNames(vals, labels)
}

#' observed_response_probabilities: Aggregate to observable outcomes
#'
#' For models where some outcomes are mapped (e.g., map_outcome_to = NA),
#' aggregate those mapped probabilities into an NA bucket (if include_na=TRUE)
#' and return only observable labels.
#'
#' @param model Race model spec
#' @param component Optional component id
#' @param include_na If TRUE, include an 'NA' entry with summed probability of
#'                   all outcomes mapped to NA. If FALSE, exclude mapped
#'                   outcomes entirely.
#' @return Named numeric vector of probabilities for observable outcomes (and NA)
observed_response_probabilities <- function(model, component = NULL, include_na = TRUE) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  labels <- names(prep[["outcomes"]])
  # Exclude alias-of outcomes from observed set
  labels <- Filter(function(lbl) {
    is.null(prep[["outcomes"]][[lbl]][['options']][['alias_of']])
  }, labels)
  # Compute base probabilities for all labels
  base <- vapply(labels, function(lbl) response_probability(model, lbl, component = component), numeric(1))
  names(base) <- labels
  # Partition into observable vs mapped-to-NA
  obs <- numeric(0)
  na_sum <- 0.0
  for (lbl in labels) {
    out_def <- prep[["outcomes"]][[lbl]]
    map_to <- out_def[['options']][['map_outcome_to']] %||% NULL
    if (is.null(map_to)) {
      # observable label
      obs[lbl] <- base[[lbl]]
    } else if (is.na(map_to)) {
      # mapped to NA
      na_sum <- na_sum + base[[lbl]]
    } else {
      # mapped to another label: aggregate into that label
      obs_lbl <- as.character(map_to)
      obs[obs_lbl] <- (obs[obs_lbl] %||% 0.0) + base[[lbl]]
    }
  }
  if (include_na && na_sum > 0) {
    obs <- c(obs, "NA" = na_sum)
  }
  if (include_na) {
    total <- sum(obs)
    resid <- 1.0 - total
    if (!is.finite(resid)) resid <- 0.0
    if (resid > .Machine$double.eps) {
      if ("NA" %in% names(obs)) {
        obs[["NA"]] <- obs[["NA"]] + resid
      } else {
        obs <- c(obs, "NA" = resid)
      }
    }
  }
  obs
}
