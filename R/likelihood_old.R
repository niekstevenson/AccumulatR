likelihood_old <- local({
  # Likelihood system for the new expression-based API (new_API.R)
  # Matches the generator logic in generator_new.R


  `%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

  prefix_sum <- function(coeff, K) {
    if (length(coeff) == 0) return(0.0)
    idx <- seq_len(min(K, length(coeff) - 1) + 1L)
    sum(coeff[idx])
  }

  .integrate_rel_tol <- function() getOption("uuber.integrate.rel.tol", 1e-5)
  .integrate_abs_tol <- function() getOption("uuber.integrate.abs.tol", 1e-6)

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

  .expr_requires_source <- function(expr, target_source) {
    if (is.null(expr) || is.null(expr[["kind"]])) return(FALSE)
    kind <- expr[["kind"]]
    if (identical(kind, "event")) {
      return(identical(expr[["source"]], target_source))
    }
    if (identical(kind, "and")) {
      args <- expr[["args"]]
      if (length(args) == 0) return(FALSE)
      return(any(vapply(args, .expr_requires_source, logical(1), target_source = target_source)))
    }
    if (identical(kind, "or")) {
      args <- expr[["args"]]
      if (length(args) == 0) return(FALSE)
      return(all(vapply(args, .expr_requires_source, logical(1), target_source = target_source)))
    }
    if (identical(kind, "guard")) {
      ref <- expr[["reference"]]
      return(.expr_requires_source(ref, target_source))
    }
    if (identical(kind, "not")) {
      return(FALSE)
    }
    FALSE
  }

  # ==============================================================================
  # PART 1: Model Preparation
  # ==============================================================================
  # Reuse the preparation logic from generator_new.R to ensure consistency

  .prepare_model_for_likelihood <- function(model) {
    if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
      model <- tables_to_model(model)
    }
    # Use the same preparation as the generator
    if (!exists("prepare_model", mode = "function")) {
      stop("prepare_model function not found - source generator_new.R first")
    }
    prepare_model(model)
  }

  # ==============================================================================
  # PART 2: Distribution Primitives
  # ==============================================================================

  # Get density function for an accumulator at time t
  .acc_density <- function(acc_def, t) {
    if (!is.finite(t) || t < 0) return(0.0)
    if (t < acc_def[['onset']]) return(0.0)

    # Account for q (probability of non-completion)
    success_prob <- 1 - acc_def[['q']]
    if (success_prob <= 0) return(0.0)

    reg <- dist_registry(acc_def[['dist']])
    if (is.null(reg) || is.null(reg[['d']])) {
      stop(sprintf("No density function for distribution '%s'", acc_def[['dist']]))
    }

    # Density of the underlying distribution at (t - onset)
    t_adj <- t - acc_def[['onset']]
    dens <- reg[['d']](t_adj, acc_def[['params']])

    # Scale by success probability
    success_prob * dens
  }

  # Get survival function for an accumulator at time t
  .acc_survival <- function(acc_def, t) {
    if (!is.finite(t)) return(0.0)
    if (t < 0) return(1.0)
    if (t < acc_def[['onset']]) return(1.0)

    # Survival = q + (1-q) * S(t - onset)
    reg <- dist_registry(acc_def[['dist']])
    if (is.null(reg) || is.null(reg[['p']])) {
      stop(sprintf("No CDF function for distribution '%s'", acc_def[['dist']]))
    }

    t_adj <- t - acc_def[['onset']]
    surv_underlying <- 1 - reg[['p']](t_adj, acc_def[['params']])

    acc_def[['q']] + (1 - acc_def[['q']]) * surv_underlying
  }

  # ==============================================================================
  # PART 3: Pool Likelihood
  # ==============================================================================

  # Density for k-th order statistic of a pool
  .pool_density <- function(prep, pool_id, component, t) {
    pool_defs <- prep[["pools"]]
    acc_defs <- prep[["accumulators"]]

    pool_def <- pool_defs[[pool_id]]
    if (is.null(pool_def)) {
      stop(sprintf("Unknown pool '%s'", pool_id))
    }

    # Get members active in this component
    members <- pool_def[['members']]
    active_members <- character(0)
    for (m in members) {
      if (!is.null(acc_defs[[m]])) {
        # Check if accumulator is active in this component
        acc_def <- acc_defs[[m]]
        comps <- acc_def[['components']]
        if (length(comps) == 0 || is.null(component) ||
            identical(component, "__default__") || component %in% comps) {
          active_members <- c(active_members, m)
        }
      } else if (!is.null(pool_defs[[m]])) {
        # Nested pool - include it
        active_members <- c(active_members, m)
      }
    }

    if (length(active_members) == 0) return(0.0)

    k <- pool_def[['k']]
    n <- length(active_members)

    if (k > n) return(0.0)
    if (k < 1) return(0.0)

    # Special case: k=1 (minimum/first-to-finish)
    if (k == 1) {
      return(.pool_k1_density(prep, active_members, component, t))
    }

    # General k-of-n order statistic density
    # Density of k-th order statistic:
    # f_{(k)}(t) = Σ_i f_i(t) * P(exactly k-1 of others have finished by t)
    # where "others" = all members except i

    total_density <- 0.0

    for (i in seq_along(active_members)) {
      mid <- active_members[[i]]

      # Density of member i at time t
      if (!is.null(acc_defs[[mid]])) {
        dens_i <- .acc_density(acc_defs[[mid]], t)
      } else if (!is.null(pool_defs[[mid]])) {
        dens_i <- .pool_density(prep, mid, component, t)
      } else {
        next
      }

      if (dens_i == 0) next

      # Get survival and CDF vectors for all others
      others <- active_members[-i]
      if (length(others) == 0) {
        # Only one member, k must be 1
        if (k == 1) {
          total_density <- total_density + dens_i
        }
        next
      }

      Svec <- numeric(length(others))
      Fvec <- numeric(length(others))

      for (j in seq_along(others)) {
        mjd <- others[[j]]
        if (!is.null(acc_defs[[mjd]])) {
          Svec[[j]] <- .acc_survival(acc_defs[[mjd]], t)
          Fvec[[j]] <- 1 - Svec[[j]]
        } else if (!is.null(pool_defs[[mjd]])) {
          Svec[[j]] <- .pool_survival(prep, mjd, component, t)
          Fvec[[j]] <- 1 - Svec[[j]]
        }
      }

      # Compute pool coefficients (probability of exactly m finished by t)
      coeffs <- pool_coeffs(Svec, Fvec)

      # We need P(exactly k-1 of others finished by t)
      # coeffs[m+1] = P(exactly m finished), so we want coeffs[k]
      if (k <= length(coeffs)) {
        prob_k_minus_1 <- coeffs[[k]]
        total_density <- total_density + dens_i * prob_k_minus_1
      }
    }

    total_density
  }

  # Density for first-to-finish (k=1) pool
  .pool_k1_density <- function(prep, member_ids, component, t) {
    acc_defs <- prep[["accumulators"]]
    pool_defs <- prep[["pools"]]

    total_density <- 0.0

    # For each member: f_i(t) * prod_{j != i} S_j(t)
    for (i in seq_along(member_ids)) {
      mid <- member_ids[[i]]

      # Get density of this member
      if (!is.null(acc_defs[[mid]])) {
        dens_i <- .acc_density(acc_defs[[mid]], t)
      } else if (!is.null(pool_defs[[mid]])) {
        # Nested pool
        dens_i <- .pool_density(prep, mid, component, t)
      } else {
        next
      }

      if (dens_i == 0) next

      # Get survival of all other members
      surv_prod <- 1.0
      for (j in seq_along(member_ids)) {
        if (i == j) next
        mjd <- member_ids[[j]]

        if (!is.null(acc_defs[[mjd]])) {
          surv_j <- .acc_survival(acc_defs[[mjd]], t)
        } else if (!is.null(pool_defs[[mjd]])) {
          # Nested pool survival
          surv_j <- .pool_survival(prep, mjd, component, t)
        } else {
          next
        }

        surv_prod <- surv_prod * surv_j
      }

      total_density <- total_density + dens_i * surv_prod
    }

    total_density
  }

  # Survival function for a pool (prob that less than k members have finished by t)
  .pool_survival <- function(prep, pool_id, component, t) {
    pool_defs <- prep[["pools"]]
    acc_defs <- prep[["accumulators"]]

    pool_def <- pool_defs[[pool_id]]
    if (is.null(pool_def)) {
      stop(sprintf("Unknown pool '%s'", pool_id))
    }

    members <- pool_def[['members']]
    k <- pool_def[['k']]

    # Get active members for this component
    active_members <- character(0)
    for (mid in members) {
      if (!is.null(acc_defs[[mid]])) {
        acc_def <- acc_defs[[mid]]
        comps <- acc_def[['components']]
        if (length(comps) == 0 || is.null(component) ||
            identical(component, "__default__") || component %in% comps) {
          active_members <- c(active_members, mid)
        }
      } else if (!is.null(pool_defs[[mid]])) {
        active_members <- c(active_members, mid)
      }
    }

    if (length(active_members) == 0) return(1.0)
    if (k > length(active_members)) return(1.0)
    if (k < 1) return(0.0)

    # Survival = P(< k finished by t) = P(≤ k-1 finished by t)
    # = sum of coeffs[0] through coeffs[k-1]
    # coeffs[m+1] = P(exactly m finished)

    # Get Svec and Fvec for all active members
    Svec <- numeric(length(active_members))
    Fvec <- numeric(length(active_members))

    for (i in seq_along(active_members)) {
      mid <- active_members[[i]]
      if (!is.null(acc_defs[[mid]])) {
        Svec[[i]] <- .acc_survival(acc_defs[[mid]], t)
        Fvec[[i]] <- 1 - Svec[[i]]
      } else if (!is.null(pool_defs[[mid]])) {
        Svec[[i]] <- .pool_survival(prep, mid, component, t)
        Fvec[[i]] <- 1 - Svec[[i]]
      }
    }

    # Compute pool coefficients
    coeffs <- pool_coeffs(Svec, Fvec)

    # Survival = P(< k finished) = sum(coeffs[1:k])
    # coeffs[m+1] = P(exactly m finished), so we want coeffs[1] + ... + coeffs[k]
    prefix_sum(coeffs, k - 1L)
  }

  # ==============================================================================
  # PART 4: Expression Evaluation (Likelihood)
  # ==============================================================================

  # Get CDF (cumulative probability) for an expression
  .eval_expr_cdf <- function(expr, t, prep, component) {
    kind <- expr[['kind']]

    if (identical(kind, "event")) {
      source_id <- expr[['source']]

      if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
        return(if (t >= 0) 1.0 else 0.0)
      }

      if (!is.null(prep[["pools"]][[source_id]])) {
        # Pool CDF: 1 - survival
        return(1 - .pool_survival(prep, source_id, component, t))
      }

      acc_def <- prep[["accumulators"]][[source_id]]
      if (!is.null(acc_def)) {
        comps <- acc_def[['components']]
        if (length(comps) > 0 && !is.null(component) &&
            !identical(component, "__default__") && !(component %in% comps)) {
          return(1.0)  # Not active = instant completion
        }
        return(1 - .acc_survival(acc_def, t))
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
      val <- tryCatch(
        stats::integrate(dens_fun, lower = 0, upper = t,
                         rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                         stop.on.error = FALSE)$value,
        error = function(e) 0.0
      )
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
        return(.pool_survival(prep, source_id, component, t))
      }

      acc_def <- prep[["accumulators"]][[source_id]]
      if (!is.null(acc_def)) {
        comps <- acc_def[['components']]
        if (length(comps) > 0 && !is.null(component) &&
            !identical(component, "__default__") && !(component %in% comps)) {
          return(1.0)  # Not active = always survive (never fires)
        }
        return(.acc_survival(acc_def, t))
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
      val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                       rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                       stop.on.error = FALSE)[['value']],
                      error = function(e) 0.0)
      s <- 1.0 - as.numeric(val)
      if (!is.finite(s)) s <- 0.0
      if (s < 0) s <- 0.0
      if (s > 1) s <- 1.0
      return(s)
    }

    return(0.0)
  }

  # Compute likelihood (density or probability) for an expression
  .eval_expr_likelihood <- function(expr, t, prep, component) {
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
    return(.pool_density(prep, source_id, component, t))
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
    return(.acc_density(acc_def, t))
  }

  return(0.0)
  }

  if (identical(kind, "and")) {
    # AND: all must finish, density at max
    # f_max(t) = Σ_i f_i(t) * Π_{j≠i} F_j(t)
    # Each term: i finishes last at t, all others finished by t
    args <- expr[['args']]
    n <- length(args)

    total_density <- 0.0
    for (i in seq_along(args)) {
      # Density of arg i at t
      dens_i <- .eval_expr_likelihood(args[[i]], t, prep, component)
      if (dens_i == 0) next

      # CDF of all others at t
      prod_cdf <- 1.0
      for (j in seq_along(args)) {
        if (i == j) next
        cdf_j <- .eval_expr_cdf(args[[j]], t, prep, component)
        prod_cdf <- prod_cdf * cdf_j
      }

      total_density <- total_density + dens_i * prod_cdf
    }

    return(total_density)
  }

  if (identical(kind, "or")) {
    # OR: first to finish
    # f_min(t) = Σ_i f_i(t) * Π_{j≠i} S_j(t)
    # Each term: i finishes first at t, all others survive past t
    args <- expr[['args']]

    total_density <- 0.0
    for (i in seq_along(args)) {
      # Density of arg i at t
      dens_i <- .eval_expr_likelihood(args[[i]], t, prep, component)
      if (dens_i == 0) next

      # Survival of all others at t
      prod_surv <- 1.0
      for (j in seq_along(args)) {
        if (i == j) next
        surv_j <- .eval_expr_survival(args[[j]], t, prep, component)
        prod_surv <- prod_surv * surv_j
      }

      total_density <- total_density + dens_i * prod_surv
    }

    return(total_density)
  }

  if (identical(kind, "not")) {
    # NOT typically doesn't have density, only affects probability
    # If we reach here, it means NOT is used in an outcome expression
    # This is unusual but could mean "outcome never happens"
    return(0.0)
  }

  if (identical(kind, "guard")) {
    # GUARD: reference fires only if blocker doesn't fire first
    # inhibit(reference, by = blocker, unless = protectors)
    #
    # Density: f_ref(t) * S_effective_blocker(t)
    # where S_effective_blocker(t) = 1 - ∫_0^t f_blocker(u) * Π_p S_protector_p(u) du

    reference <- expr[['reference']]
    blocker <- expr[['blocker']]
    unless_list <- expr[['unless']] %||% list()

    # Density of reference at t
    dens_ref <- .eval_expr_likelihood(reference, t, prep, component)
    if (dens_ref == 0) return(0.0)

    # Fast path: no protectors -> blocker survival is sufficient
    S_eff <- NULL
    if (length(unless_list) == 0) {
      S_eff <- .eval_expr_survival(blocker, t, prep, component)
    } else {
      # Build closures for f_blocker and S_protectors
      f_blocker <- function(tt) {
        .eval_expr_likelihood(blocker, tt, prep, component)
      }
      S_protector_prod <- function(tt) {
        if (length(unless_list) == 0) return(rep(1.0, length(tt)))
        vapply(tt, function(ui) {
          vals <- vapply(unless_list, function(unl) .eval_expr_survival(unl, ui, prep, component), numeric(1))
          prod(vals)
        }, numeric(1))
      }
      # Effective blocker survival up to time t
      S_eff <- tryCatch({
        if (!is.finite(t) || t <= 0) 1.0 else {
          val <- stats::integrate(function(u) {
            vapply(u, function(ui) f_blocker(ui) * S_protector_prod(ui), numeric(1))
          }, lower = 0, upper = t,
          rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(), stop.on.error = FALSE)[['value']]
          s <- 1.0 - as.numeric(val)
          if (!is.finite(s)) 0.0 else max(0.0, min(1.0, s))
        }
      }, error = function(e) 1.0)
    }
    if (!is.finite(S_eff)) S_eff <- 0.0
    if (S_eff < 0) S_eff <- 0.0
    if (S_eff > 1) S_eff <- 1.0
    return(dens_ref * S_eff)
  }

  stop(sprintf("Unsupported expression kind '%s'", kind))
  }

  # ==============================================================================
  # PART 5: Outcome Likelihood
  # ==============================================================================

  # Integrand function for probability integration
  .integrand_outcome_density <- function(t, expr, prep, component) {
    vapply(t, function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      .eval_expr_likelihood(expr, tt, prep, component)
    }, numeric(1))
  }

  # Compute probability by integrating race density
  # If competitor_exprs is provided, integrates:
  #   f_r(t) * prod_j S_j(t)  over t in [0, upper]
  .integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf, competitor_exprs = NULL) {
    # Set integration upper bound
    if (!is.finite(upper_limit)) {
      upper_limit <- Inf
    }

    result <- tryCatch({
      if (is.null(competitor_exprs) || length(competitor_exprs) == 0) {
        stats::integrate(
          .integrand_outcome_density,
          lower = 0,
          upper = upper_limit,
          expr = expr,
          prep = prep,
          component = component,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )
      } else {
        integrand <- function(t) {
          vapply(t, function(tt) {
            if (!is.finite(tt) || tt < 0) return(0.0)
            dens_r <- .eval_expr_likelihood(expr, tt, prep, component)
            if (dens_r == 0) return(0.0)
            surv_comp <- 1.0
            for (cexpr in competitor_exprs) {
              surv_comp <- surv_comp * .eval_expr_survival(cexpr, tt, prep, component)
              if (surv_comp == 0) break
            }
            dens_r * surv_comp
          }, numeric(1))
        }
        stats::integrate(
          integrand,
          lower = 0,
          upper = upper_limit,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )
      }
    }, error = function(e) {
      list(value = 0.0, message = as.character(e))
    })

    val <- as.numeric(result[['value']])
    if (!is.finite(val)) val <- 0.0
    if (val < 0) val <- 0.0
    if (val > 1) val <- 1.0
    val
  }

  .expr_sources <- function(expr, prep) {
    if (is.null(expr) || is.null(expr[['kind']])) return(character(0))
    kind <- expr[['kind']]
    if (identical(kind, "event")) {
      source_id <- expr[['source']]
      if (!is.null(prep[["pools"]][[source_id]])) {
        members <- prep[["pools"]][[source_id]][['members']] %||% character(0)
        unique(unlist(lapply(members, function(mid) {
          if (!is.null(prep[["accumulators"]][[mid]])) return(mid)
          if (!is.null(prep[["pools"]][[mid]])) {
            return(.expr_sources(list(kind = "event", source = mid), prep))
          }
          character(0)
        })))
      } else {
        source_id
      }
    } else if (identical(kind, "and") || identical(kind, "or")) {
      unique(unlist(lapply(expr[['args']], .expr_sources, prep = prep)))
    } else if (identical(kind, "guard")) {
      ref_src <- .expr_sources(expr[['reference']], prep)
      blk_src <- .expr_sources(expr[['blocker']], prep)
      unless_src <- unique(unlist(lapply(expr[['unless']] %||% list(), .expr_sources, prep = prep)))
      unique(c(ref_src, blk_src, unless_src))
    } else if (identical(kind, "not")) {
      .expr_sources(expr[['arg']], prep)
    } else {
      character(0)
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

  # Vectorized integrands for shared-gate pair; returns list(prob_x, prob_y)
  .shared_gate_pair_probs <- function(prep, component, pair_info,
                                      rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol()) {
    get_f <- function(id) {
      function(tt) {
        vapply(tt, function(ti) {
          if (!is.finite(ti) || ti < 0) return(0.0)
          if (!is.null(prep[["pools"]][[id]])) return(.pool_density(prep, id, component, ti))
          if (!is.null(prep[["accumulators"]][[id]])) return(.acc_density(prep[["accumulators"]][[id]], ti))
          0.0
        }, numeric(1))
      }
    }
    get_S <- function(id) {
      function(tt) {
        vapply(tt, function(ti) {
          if (!is.finite(ti)) return(0.0)
          if (ti < 0) return(1.0)
          if (!is.null(prep[["pools"]][[id]])) return(.pool_survival(prep, id, component, ti))
          if (!is.null(prep[["accumulators"]][[id]])) return(.acc_survival(prep[["accumulators"]][[id]], ti))
          0.0
        }, numeric(1))
      }
    }
    get_F <- function(id) {
      Sfun <- get_S(id)
      function(tt) { 1.0 - Sfun(tt) }
    }

    fX <- get_f(pair_info[['x_id']]); fY <- get_f(pair_info[['y_id']]); fC <- get_f(pair_info[['c_id']])
    FX <- get_F(pair_info[['x_id']]); FY <- get_F(pair_info[['y_id']]); FC <- get_F(pair_info[['c_id']])
    SX <- get_S(pair_info[['x_id']]); SY <- get_S(pair_info[['y_id']])

    palloc_x <- function(t) {
      denom <- FX(t) * FY(t)
      num <- tryCatch(stats::integrate(function(u) { fX(u) * FY(u) }, lower = 0, upper = t,
                                       rel.tol = rel.tol, abs.tol = abs.tol,
                                       stop.on.error = FALSE)[['value']], error = function(e) 0.0)
      out <- 1.0 - (as.numeric(num) / as.numeric(denom))
      if (!is.finite(out)) 0.0 else max(0.0, min(1.0, out))
    }

    integrand_x <- function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        term1 <- fX(tt) * FC(tt) * SY(tt)
        termCother <- fC(tt) * FX(tt) * SY(tt)
        den <- FX(tt) * FY(tt)
        tie <- if (is.finite(den) && den > 0) fC(tt) * FX(tt) * FY(tt) * palloc_x(tt) else 0.0
        term1 + termCother + tie
      }, numeric(1))
    }

    integrand_y <- function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        term1 <- fY(tt) * FC(tt) * SX(tt)
        termCother <- fC(tt) * FY(tt) * SX(tt)
        den <- FX(tt) * FY(tt)
        # symmetry: palloc_y = 1 - palloc_x
        tie <- if (is.finite(den) && den > 0) fC(tt) * FX(tt) * FY(tt) * (1.0 - palloc_x(tt)) else 0.0
        term1 + termCother + tie
      }, numeric(1))
    }

    upper <- Inf
    px <- tryCatch(stats::integrate(integrand_x, lower = 0, upper = upper,
                                    rel.tol = rel.tol, abs.tol = abs.tol,
                                    stop.on.error = FALSE)[['value']], error = function(e) 0.0)
    py <- tryCatch(stats::integrate(integrand_y, lower = 0, upper = upper,
                                    rel.tol = rel.tol, abs.tol = abs.tol,
                                    stop.on.error = FALSE)[['value']], error = function(e) 0.0)
    c(as.numeric(px), as.numeric(py))
  }

  # Compute likelihood for a single trial/outcome
  .outcome_likelihood <- function(outcome_label, rt, prep, component) {
    outcome_defs <- prep[["outcomes"]]

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
          # GUESS probability/density comes from diverting mass from base labels
          if (is.na(rt)) {
            total_prob <- 0.0
            for (label in names(guess_weights)) {
              if (!label %in% names(outcome_defs)) next
              out_def <- outcome_defs[[label]]
              expr <- out_def[['expr']]
              comp_labels <- setdiff(names(outcome_defs), label)
              if (length(comp_labels) > 0) {
                comp_labels <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels)
              }
              competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              prob_outcome <- .integrate_outcome_probability(expr, prep, component, deadline, competitor_exprs = competitor_exprs)
              keep_prob <- as.numeric(guess_weights[[label]])
              guess_prob <- 1.0 - keep_prob
              total_prob <- total_prob + prob_outcome * guess_prob
            }
            return(total_prob)
          } else {
            total_dens <- 0.0
            for (label in names(guess_weights)) {
              if (!label %in% names(outcome_defs)) next
              out_def <- outcome_defs[[label]]
              expr <- out_def[['expr']]
              # density of label at rt times survival of all competitors excluding this label
              dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
              if (dens_r == 0) next
              comp_labels <- setdiff(names(outcome_defs), label)
              if (length(comp_labels) > 0) {
                comp_labels <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels)
              }
              competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              surv_prod <- .compute_survival_product(expr, competitor_exprs, prep, component, rt)
              keep_prob <- as.numeric(guess_weights[[label]])
              guess_prob <- 1.0 - keep_prob
              total_dens <- total_dens + guess_prob * dens_r * surv_prod
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
                # competitors exclude the mapped source label
                comp_labels_map <- setdiff(names(outcome_defs), label)
                competitor_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
                return(.integrate_outcome_probability(expr, prep, component, deadline, competitor_exprs = competitor_exprs_map))
              } else {
                # Race density at rt for the mapped source outcome
                dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
                if (dens_r == 0) return(0.0)
                comp_labels_map <- setdiff(names(outcome_defs), label)
                competitor_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
                surv_comp <- .compute_survival_product(expr, competitor_exprs_map, prep, component, rt)
                return(dens_r * surv_comp)
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

    # Build competitor exprs for this outcome (exclude alias-of, special outcomes, and siblings mapping to the same top-level label)
    comp_labels <- setdiff(names(outcome_defs), outcome_label)
    if (length(comp_labels) > 0) {
      # Determine this outcome's top-level target (itself or map_outcome_to)
      this_target <- options[['map_outcome_to']] %||% outcome_label
      comp_labels <- Filter(function(lbl) {
        # skip alias-of
        if (!is.null(outcome_defs[[lbl]][['options']][['alias_of']])) return(FALSE)
        # skip special sources
        expr_lbl <- outcome_defs[[lbl]][['expr']]
        if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
          src <- expr_lbl[['source']]
          if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
        }
        # skip siblings that map to the same top-level target label
        comp_target <- outcome_defs[[lbl]][['options']][['map_outcome_to']] %||% lbl
        if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
        TRUE
      }, comp_labels)
    }
    competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
    if (length(competitor_exprs) > 0) {
      expr_sources <- .expr_sources(expr, prep)
      competitor_exprs <- lapply(competitor_exprs, function(comp_expr) {
        if (!is.null(comp_expr[['kind']]) && identical(comp_expr[['kind']], "guard")) {
          blocker_sources <- .expr_sources(comp_expr[['blocker']], prep)
          if (length(expr_sources) > 0 && all(expr_sources %in% blocker_sources)) {
            return(comp_expr[['reference']])
          }
        }
        comp_expr
      })
    }

    if (!is.null(expr[['kind']]) && identical(expr[['kind']], "guard") && length(competitor_exprs) > 0) {
      blocker_sources <- .expr_sources(expr[['blocker']], prep)
      keep_idx <- vapply(seq_along(competitor_exprs), function(i) {
        comp_sources <- .expr_sources(competitor_exprs[[i]], prep)
        # keep if competitor involves sources outside the blocker set
        length(comp_sources) == 0 || any(!comp_sources %in% blocker_sources)
      }, logical(1))
      competitor_exprs <- competitor_exprs[keep_idx]
      comp_labels <- comp_labels[keep_idx]
    }

    if (length(competitor_exprs) > 0) {
      primary_guards <- .collect_guards(expr)
      if (length(primary_guards) > 0) {
        guard_keep <- vapply(seq_along(competitor_exprs), function(i) {
          comp_expr <- competitor_exprs[[i]]
          any_conflict <- any(vapply(primary_guards, function(pg) .guards_conflict(pg, comp_expr), logical(1)))
          !any_conflict
        }, logical(1))
        competitor_exprs <- competitor_exprs[guard_keep]
        comp_labels <- comp_labels[guard_keep]
      }
    }

    if (!is.null(expr[['kind']]) && identical(expr[['kind']], "guard") && length(competitor_exprs) > 0) {
      blocker <- expr[['blocker']]
      blocker_source <- if (.expr_is_event(blocker)) blocker[['source']] else NULL
      if (!is.null(blocker_source)) {
        req_keep <- vapply(seq_along(competitor_exprs), function(i) {
          !.expr_requires_source(competitor_exprs[[i]], blocker_source)
        }, logical(1))
        competitor_exprs <- competitor_exprs[req_keep]
        comp_labels <- comp_labels[req_keep]
      }
    }

    # Build donor (guess) contributions: outcomes that guess into this label
    donors <- list()
    if (length(comp_labels) > 0) {
      for (dlbl in comp_labels) {
        dopt <- outcome_defs[[dlbl]][['options']][['guess']] %||% NULL
        if (is.null(dopt)) next
        glabels <- dopt[['labels']] %||% NULL
        gweights <- dopt[['weights']] %||% NULL
        if (is.null(glabels) || is.null(gweights)) next
        idx <- which(as.character(glabels) == outcome_label)
        if (length(idx) == 0) next
        w <- as.numeric(gweights[[idx[[1]]]])
        rt_policy <- dopt[['rt_policy']] %||% "keep"
        donors[[length(donors) + 1L]] <- list(expr = outcome_defs[[dlbl]][['expr']], weight = w, rt_policy = rt_policy)
      }
    }

    # Detect shared-gate AND pattern: (X & C) vs (Y & C)
    shared_gate_info <- .detect_shared_gate(outcome_defs, outcome_label)

    # Helper closures to get density/CDF for pools or accumulators
    get_f <- function(id) {
      function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        if (!is.null(prep[["pools"]][[id]])) return(.pool_density(prep, id, component, tt))
        if (!is.null(prep[["accumulators"]][[id]])) return(.acc_density(prep[["accumulators"]][[id]], tt))
        0.0
      }
    }
    get_S <- function(id) {
      function(tt) {
        if (!is.finite(tt)) return(0.0)
        if (tt < 0) return(1.0)
        if (!is.null(prep[["pools"]][[id]])) return(.pool_survival(prep, id, component, tt))
        if (!is.null(prep[["accumulators"]][[id]])) return(.acc_survival(prep[["accumulators"]][[id]], tt))
        0.0
      }
    }
    get_F <- function(id) {
      Sfun <- get_S(id)
      function(tt) { 1.0 - Sfun(tt) }
    }

    # Handle NA/infinite RT by integration of race density
    if (is.na(rt) || !is.finite(rt) || rt < 0) {
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep[["default_deadline"]]
      if (!is.null(shared_gate_info)) {
        fX <- get_f(shared_gate_info[['x_id']]); fC <- get_f(shared_gate_info[['c_id']])
        FX <- get_F(shared_gate_info[['x_id']]); FY <- get_F(shared_gate_info[['y_id']])
        SX <- get_S(shared_gate_info[['x_id']]); SY <- get_S(shared_gate_info[['y_id']]); FC <- get_F(shared_gate_info[['c_id']])
        integrand <- function(t) {
          vapply(t, function(tt) {
            # X last term
            term1 <- fX(tt) * FC(tt) * SY(tt)
            # C last tie term with proportional allocation
            denom <- FX(tt) * FY(tt)
            # C last with other > t: fC * F_X * S_Y
            termCother <- fC(tt) * FX(tt) * SY(tt)
            if (!is.finite(denom) || denom <= 0) return(term1 + termCother)
            # p(X<Y | both <= t) = 1 - [ ∫ fX(u) F_Y(u) du / (F_X(t) F_Y(t)) ]
            palloc <- tryCatch({
              val <- stats::integrate(function(u) {
                vapply(u, function(ui) fX(ui) * FY(ui), numeric(1))
              }, lower = 0, upper = tt,
              rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(), stop.on.error = FALSE)[['value']]
              1.0 - (as.numeric(val) / denom)
            }, error = function(e) 0.0)
            term2 <- fC(tt) * FX(tt) * FY(tt) * palloc
            term1 + termCother + term2
          }, numeric(1))
        }
        # Apply component-level keep (if any) to observed probability
        keep_prob <- 1.0
        gp_comp <- .get_component_attr(prep, component, "guess")
        if (!is.null(gp_comp) && !is.null(gp_comp[['weights']])) {
          kp <- gp_comp[['weights']][[outcome_label]] %||% gp_comp[['weights']][[normalize_label(outcome_label)]] %||% NULL
          if (!is.null(kp)) keep_prob <- as.numeric(kp)
        }
        res <- tryCatch(stats::integrate(integrand, lower = 0, upper = if (is.finite(deadline)) deadline else Inf,
                                         rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(), stop.on.error = FALSE)[['value']],
                        error = function(e) 0.0)
        return(as.numeric(keep_prob * res))
      } else {
        # Generic race integral with donor contributions
        integrand <- function(t) {
          vapply(t, function(tt) {
            if (!is.finite(tt) || tt < 0) return(0.0)
            base <- .eval_expr_likelihood(expr, tt, prep, component)
            if (base > 0 && length(competitor_exprs) > 0) {
              base <- base * .compute_survival_product(expr, competitor_exprs, prep, component, tt)
            }
            add <- 0.0
            if (length(donors) > 0) {
              for (d in donors) {
                dens_d <- .eval_expr_likelihood(d[['expr']], tt, prep, component)
                if (dens_d == 0) next
                donor_src <- NULL
                for (lblx in names(outcome_defs)) { if (identical(outcome_defs[[lblx]][['expr']], d[['expr']])) { donor_src <- lblx; break } }
                d_comp_labels <- setdiff(names(outcome_defs), donor_src)
                d_comp_exprs <- if (length(d_comp_labels) > 0) lapply(d_comp_labels, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
                surv_prod_d <- if (length(d_comp_exprs) > 0) prod(vapply(d_comp_exprs, function(ce) .eval_expr_survival(ce, tt, prep, component), numeric(1))) else 1.0
                add <- add + d[['weight']] * dens_d * surv_prod_d
              }
            }
            base + add
          }, numeric(1))
        }
        res <- tryCatch(stats::integrate(integrand, lower = 0, upper = if (is.finite(deadline)) deadline else Inf,
                                         rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(), stop.on.error = FALSE)[['value']],
                        error = function(e) 0.0)
        return(as.numeric(res))
      }
    }

    # Race density at observed RT
    if (!is.null(shared_gate_info)) {
      fX <- get_f(shared_gate_info[['x_id']]); fC <- get_f(shared_gate_info[['c_id']])
      FX <- get_F(shared_gate_info[['x_id']]); FY <- get_F(shared_gate_info[['y_id']])
      SX <- get_S(shared_gate_info[['x_id']]); SY <- get_S(shared_gate_info[['y_id']]); FC <- get_F(shared_gate_info[['c_id']])
      # X last term: A finishes at t, C finished before t => fA * F_C * (1 - F_B)
      term1 <- fX(rt) * FC(rt) * SY(rt)
      # C last tie term with proportional allocation
      denom <- FX(rt) * FY(rt)
      term2 <- 0.0
      termCother <- fC(rt) * FX(rt) * SY(rt)
      if (is.finite(denom) && denom > 0) {
        palloc <- tryCatch({
          val <- stats::integrate(function(u) {
            vapply(u, function(ui) fX(ui) * FY(ui), numeric(1))
          }, lower = 0, upper = rt,
          rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(), stop.on.error = FALSE)[['value']]
          1.0 - (as.numeric(val) / denom)
        }, error = function(e) 0.0)
        term2 <- fC(rt) * FX(rt) * FY(rt) * palloc
      }
      base_val <- term1 + termCother + term2
    } else {
      dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
      if (dens_r == 0) return(0.0)
      surv_comp <- .compute_survival_product(expr, competitor_exprs, prep, component, rt)
      base_val <- dens_r * surv_comp
      # Include donor mass (e.g., timeout guesses) that keeps the observed RT
      if (length(donors) > 0) {
        donor_add <- 0.0
        for (d in donors) {
          if (identical(d[['rt_policy']], "na")) next
          dens_d <- .eval_expr_likelihood(d[['expr']], rt, prep, component)
          if (dens_d == 0) next
          donor_src <- NULL
          for (lblx in names(outcome_defs)) {
            if (identical(outcome_defs[[lblx]][['expr']], d[['expr']])) {
              donor_src <- lblx
              break
            }
          }
          d_comp_labels <- setdiff(names(outcome_defs), donor_src)
          d_comp_exprs <- if (length(d_comp_labels) > 0) lapply(d_comp_labels, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
          surv_prod_d <- if (length(d_comp_exprs) > 0) prod(vapply(d_comp_exprs, function(ce) .eval_expr_survival(ce, rt, prep, component), numeric(1))) else 1.0
          donor_add <- donor_add + d[['weight']] * dens_d * surv_prod_d
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
  #' @param cache_env Optional environment for memoising per-trial likelihoods
  #' @return Log-likelihood value with per-trial attributes
  compute_loglik <- function(model, data, cache_env = NULL) {
    if (!is.data.frame(data)) {
      stop("data must be a data frame")
    }

    if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
      model <- tables_to_model(model)
    }

    if (!"outcome" %in% names(data) || !"rt" %in% names(data)) {
      stop("data must contain 'outcome' and 'rt' columns")
    }

    # Shared cache for memoisation (allows reuse across repeated evaluations)
    if (is.null(cache_env)) {
      cache_env <- new.env(parent = emptyenv())
    } else if (!is.environment(cache_env)) {
      stop("cache_env must be an environment or NULL")
    }

    # Prepare model
    prep <- .prepare_model_for_likelihood(model)

    # Check if mixture model and whether weights are parameters (auto-estimate)
    comp_info <- prep[["components"]]
    is_mixture <- length(comp_info[['ids']]) > 1
    has_weight_param <- FALSE
    if (!is.null(comp_info[['has_weight_param']])) {
      has_weight_param <- any(comp_info[['has_weight_param']])
    }

    # Helper to compute likelihood for a single data row (with memoisation key)
    eval_row_likelihood <- function(idx) {
      outcome <- as.character(data[['outcome']][[idx]])
      rt_val <- as.numeric(data[['rt']][[idx]])
      rt_key <- if (is.na(rt_val)) "NA" else format(rt_val, digits = 22, scientific = FALSE, trim = TRUE)
      component_key <- ""

      if (is_mixture) {
        if (has_weight_param || !("component" %in% names(data))) {
          # Mixture-averaged likelihood (weights may be parameters)
          comp_ids <- comp_info[['ids']]
          weights <- comp_info[['weights']]
          total_lik <- 0.0
          for (j in seq_along(comp_ids)) {
            comp_lik <- .outcome_likelihood(outcome, rt_val, prep, comp_ids[[j]])
            total_lik <- total_lik + weights[[j]] * comp_lik
          }
          component_key <- "__mixture__"
          lik <- total_lik
        } else {
          component <- as.character(data[['component']][[idx]])
          component_key <- component
          lik <- .outcome_likelihood(outcome, rt_val, prep, component)
        }
      } else {
        component <- if (length(comp_info[['ids']]) > 0) comp_info[['ids']][[1]] else "__default__"
        component_key <- component
        lik <- .outcome_likelihood(outcome, rt_val, prep, component)
      }

      key <- paste(component_key, outcome, rt_key, sep = "|")

      if (!is.null(cache_env[[key]])) {
        return(list(key = key, value = cache_env[[key]]))
      }

      cache_env[[key]] <- lik
      list(key = key, value = lik)
    }

    # Pre-evaluate unique rows to make reuse/caching explicit
    n_rows <- nrow(data)
    if (n_rows == 0) {
      total_ll <- 0.0
      attr(total_ll, "contributions") <- numeric(0)
      attr(total_ll, "log_contributions") <- numeric(0)
      attr(total_ll, "cache") <- cache_env
      attr(total_ll, "keys") <- character(0)
      return(total_ll)
    }

    raw_keys <- character(n_rows)
    for (i in seq_len(n_rows)) {
      outcome <- as.character(data[['outcome']][[i]])
      rt_val <- as.numeric(data[['rt']][[i]])
      rt_key <- if (is.na(rt_val)) "NA" else format(rt_val, digits = 22, scientific = FALSE, trim = TRUE)
      component_key <- ""
      if (is_mixture) {
        if (has_weight_param || !("component" %in% names(data))) {
          component_key <- "__mixture__"
        } else {
          component_key <- as.character(data[['component']][[i]])
        }
      } else {
        component_key <- if (length(comp_info[['ids']]) > 0) comp_info[['ids']][[1]] else "__default__"
      }
      raw_keys[[i]] <- paste(component_key, outcome, rt_key, sep = "|")
    }

    unique_keys <- unique(raw_keys)
    key_to_indices <- split(seq_len(n_rows), raw_keys)

    for (key in unique_keys) {
      idx <- key_to_indices[[key]][[1]]
      res <- eval_row_likelihood(idx)
      # Ensure cache is populated even if likelihood computed downstream encountered Inf/NaN
      if (is.null(cache_env[[res[['key']]]])) {
        cache_env[[res[['key']]]] <- res[['value']]
      }
    }

    lik_values <- vapply(raw_keys, function(key) cache_env[[key]], numeric(1))
    log_likes <- ifelse(is.finite(lik_values) & lik_values > 0, log(lik_values), -Inf)

    total_ll <- sum(log_likes)
    attr(total_ll, "contributions") <- lik_values
    attr(total_ll, "log_contributions") <- log_likes
    attr(total_ll, "keys") <- raw_keys
    attr(total_ll, "cache") <- cache_env

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
    function(t) vapply(t, function(tt) .eval_expr_likelihood(event_expr, tt, prep, component), numeric(1))
  }

  .event_survival_fn <- function(source_id, prep, component) {
    event_expr <- list(kind = "event", source = source_id, k = NULL)
    function(t) vapply(t, function(tt) .eval_expr_survival(event_expr, tt, prep, component), numeric(1))
  }

  .guard_pair_joint_survival <- function(t, prep, component, ref_a, block_a, extras) {
    if (t <= 0) return(1.0)
    f_a <- .event_density_fn(ref_a, prep, component)
    S_a <- .event_survival_fn(ref_a, prep, component)
    f_b <- .event_density_fn(block_a, prep, component)
    S_b <- .event_survival_fn(block_a, prep, component)
    S_extra_t <- 1.0
    if (length(extras) > 0) {
      S_extra_t <- prod(vapply(extras, function(src) .event_survival_fn(src, prep, component)(t), numeric(1)))
    }
    rel_tol <- .integrate_rel_tol()
    abs_tol <- .integrate_abs_tol()
    term_case2 <- if (is.finite(t)) {
      stats::integrate(function(v) f_a(v) * S_b(v), lower = t, upper = Inf,
                       rel.tol = rel_tol, abs.tol = abs_tol,
                       stop.on.error = FALSE)$value
    } else 0.0
    term_s1_gt_t <- stats::integrate(function(u) f_b(u) * S_a(u), lower = t, upper = Inf,
                                     rel.tol = rel_tol, abs.tol = abs_tol,
                                     stop.on.error = FALSE)$value
    if (!is.finite(term_s1_gt_t)) term_s1_gt_t <- 0.0
    term_s1_le_t <- stats::integrate(function(u) f_b(u) * S_a(u), lower = 0, upper = t,
                                     rel.tol = rel_tol, abs.tol = abs_tol,
                                     stop.on.error = FALSE)$value
    if (!is.finite(term_s1_le_t)) term_s1_le_t <- 0.0
    val <- term_case2 + term_s1_gt_t + S_extra_t * term_s1_le_t
    if (!is.finite(val)) val <- 0.0
    if (val < 0) val <- 0.0
    if (val > 1) val <- 1.0
    val
  }

  .compute_survival_product <- function(outcome_expr, competitor_exprs, prep, component, t) {
    if (length(competitor_exprs) == 0) return(1.0)
    clusters <- .cluster_competitors(competitor_exprs, prep)
    prod_val <- 1.0
    for (cluster in clusters) {
      exprs <- cluster$exprs
      cluster_val <- NULL
      if (length(exprs) == 2) {
        info1 <- .extract_guard_plus(exprs[[1]])
        info2 <- .extract_guard_plus(exprs[[2]])
        if (!is.null(info1) && !is.null(info2)) {
          if (identical(info1$ref, info2$block) && identical(info1$block, info2$ref)) {
            guard_go1 <- info1
            guard_stop <- info2
          } else if (identical(info2$ref, info1$block) && identical(info2$block, info1$ref)) {
            guard_go1 <- info2
            guard_stop <- info1
          } else {
            guard_go1 <- NULL
            guard_stop <- NULL
          }
          if (!is.null(guard_go1) && length(guard_stop$extras) == length(unique(guard_stop$extras))) {
            cluster_val <- .guard_pair_joint_survival(t, prep, component,
                                                      ref_a = guard_go1$ref,
                                                      block_a = guard_go1$block,
                                                      extras = guard_stop$extras)
          }
        }
      }
      if (is.null(cluster_val)) {
        cluster_val <- prod(vapply(exprs, function(ce) .eval_expr_survival(ce, t, prep, component), numeric(1)))
      }
      prod_val <- prod_val * cluster_val
    }
    prod_val
  }
  list(
    compute_loglik = compute_loglik,
    response_probability = response_probability,
    response_probabilities = response_probabilities
  )
})

compute_loglik_old <- likelihood_old$compute_loglik
response_probability_old <- likelihood_old$response_probability
response_probabilities_old <- likelihood_old$response_probabilities
