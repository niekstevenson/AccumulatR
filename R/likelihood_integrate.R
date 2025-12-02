# PART 5: Outcome Likelihood
# ==============================================================================

# Integrand for scenario-conditioned probability (density already encodes competitors)
.integrand_outcome_density <- function(t, expr, prep, component, competitor_exprs,
                                       state = NULL, trial_rows = NULL) {
  state <- state %||% .eval_state_create()
  times <- as.numeric(t)
  if (!length(times)) return(numeric(0))
  res <- .scenario_density_with_competitors(expr, times, prep, component,
                                            competitor_exprs = competitor_exprs,
                                            state = state,
                                            trial_rows = trial_rows)
  as.numeric(res)
}

# Scenario-conditioned probability of an expression over [0, upper_limit].
# competitor_exprs retained for backward compatibility but ignored.
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf,
                                           competitor_exprs = NULL,
                                           state = NULL,
                                           trial_rows = NULL) {
  if (!is.finite(upper_limit)) upper_limit <- Inf
  state <- state %||% .eval_state_create()
  component_key <- .cache_component_key(component)
  trial_rows_df <- NULL
  if (!is.null(trial_rows) && inherits(trial_rows, "data.frame") && nrow(trial_rows) > 0L) {
    trial_rows_df <- as.data.frame(trial_rows)
  }
  params_hash <- .params_hash_from_rows(trial_rows_df)
  competitor_exprs <- competitor_exprs %||% list()
  expr_kind <- expr[['kind']]
  if (identical(expr_kind, "guard")) {
    # guard caches removed; always compute directly
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
  upper_val <- as.numeric(upper_limit)[1]
  if (is.na(upper_val)) {
    use_native <- FALSE
  } else {
    use_native <- !is.null(compiled) && !is.null(comp_ids)
  }
  if (use_native) {
    native_ctx <- .prep_native_context(prep)
    if (!is.null(trial_rows_df) && nrow(trial_rows_df) > 0L) {
      val <- native_outcome_probability_params_cpp(
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
      val <- native_outcome_probability_cpp(
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
        trial_rows = trial_rows_df
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
.shared_gate_pair_probs <- function(prep, component, pair_info,
                                    trial_rows = NULL,
                                    trial_state = NULL) {
  labels <- c(pair_info$label_x, pair_info$label_y)
  vapply(labels, function(lbl) {
    .outcome_likelihood(
      outcome_label = lbl,
      rt = NA_real_,
      prep = prep,
      component = component,
      trial_rows = trial_rows,
      state = trial_state
    )
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
                                state = NULL) {
  outcome_defs <- prep[["outcomes"]]
  competitor_map <- .prep_competitors(prep) %||% list()
  state <- state %||% .eval_state_create()
  trial_rows_df <- NULL
  if (!is.null(trial_rows) && inherits(trial_rows, "data.frame") && nrow(trial_rows) > 0L) {
    trial_rows_df <- as.data.frame(trial_rows)
  }
  component_key <- .cache_component_key(component)
  params_hash <- .params_hash_from_rows(trial_rows_df)
  
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
                                                           trial_rows = trial_rows_df)
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
                                                         trial_rows = trial_rows_df)
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
        competitor_sig <- if (length(comp_exprs_map) == 0L) {
          "none"
        } else {
          paste(sort(vapply(comp_exprs_map, .expr_signature, character(1))), collapse = ",")
        }
        source_sig <- .expr_signature(expr)
        na_cache_key <- .na_integral_cache_key(
          source_signature = source_sig,
          component_key = component_key,
          upper_limit = deadline,
          competitor_signature = competitor_sig,
          params_hash = params_hash
        )
        cached_na <- .na_integral_fetch(prep, component_key, na_cache_key)
        if (!is.null(cached_na)) {
          return(as.numeric(cached_na))
        }
        val <- .integrate_outcome_probability(expr, prep, component, deadline,
                                              competitor_exprs = comp_exprs_map,
                                              state = state,
                                              trial_rows = trial_rows_df)
        .na_integral_store(prep, component_key, na_cache_key, val)
        return(val)
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
                                                           trial_rows = trial_rows_df)
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
                                                   trial_rows = trial_rows_df)
        if (base == 0) return(0.0)
        add <- 0.0
        if (length(donors) > 0) {
          for (d in donors) {
            dens_d <- .scenario_density_with_competitors(d[['expr']], tt, prep, component,
                                                         competitor_exprs = d[['competitors']] %||% list(),
                                                         state = state,
                                                         trial_rows = trial_rows_df)
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
                                                   trial_rows = trial_rows_df)
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
                                                     trial_rows = trial_rows_df)
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
