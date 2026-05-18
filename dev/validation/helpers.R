`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

validation_source_rebuild <- function(repo_root) {
  source(file.path(repo_root, "R", "helpers.R"))
  source(file.path(repo_root, "R", "model_definition.R"))
  source(file.path(repo_root, "R", "semantic_bridge.R"))
  source(file.path(repo_root, "R", "likelihood_param_interface.R"))
}

acc_parts <- function(prefix, params) {
  t0_name <- paste0(prefix, ".t0")
  list(
    m = unname(params[[paste0(prefix, ".m")]]),
    s = unname(params[[paste0(prefix, ".s")]]),
    q = 0.0,
    t0 = if (t0_name %in% names(params)) unname(params[[t0_name]]) else 0.0
  )
}

validation_separate_all_parameters <- function(spec) {
  lookup <- .default_parameter_lookup(spec)
  groups <- unique(unname(lookup))
  groups <- groups[nzchar(groups)]
  if (length(groups) == 0L) {
    return(spec)
  }
  spec$parameters <- .normalize_parameter_spec(
    separate = stats::setNames(rep(list(TRUE), length(groups)), groups)
  )
  spec
}

validation_spec_for_params <- function(spec, params) {
  if (all(names(params) %in% par_names(spec))) {
    return(spec)
  }
  validation_separate_all_parameters(spec)
}

acc_pdf_scalar <- function(t, p) {
  if (!is.finite(t) || t < p$t0) {
    return(0.0)
  }
  out <- (1.0 - p$q) * dlnorm(t - p$t0, meanlog = p$m, sdlog = p$s)
  if (!is.finite(out) || out < 0.0) {
    return(0.0)
  }
  out
}

acc_cdf_scalar <- function(t, p) {
  if (is.infinite(t) && t > 0.0) {
    return(1.0 - p$q)
  }
  if (!is.finite(t) || t < p$t0) {
    return(0.0)
  }
  out <- (1.0 - p$q) * plnorm(t - p$t0, meanlog = p$m, sdlog = p$s)
  min(max(out, 0.0), 1.0)
}

acc_survival_scalar <- function(t, p) {
  1.0 - acc_cdf_scalar(t, p)
}

integrate_scalar <- function(fn,
                             lower,
                             upper,
                             rel_tol = 1e-10,
                             abs_tol = 0,
                             subdivisions = 400L) {
  if (!is.finite(upper) || upper <= lower) {
    return(0.0)
  }
  integrate(
    function(x) vapply(x, fn, numeric(1)),
    lower = lower,
    upper = upper,
    rel.tol = rel_tol,
    abs.tol = abs_tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )$value
}

inhibit_density_scalar <- function(ref_pdf, blocker_cdf, t) {
  ref_pdf(t) * max(0.0, 1.0 - blocker_cdf(t))
}

inhibit_cdf_scalar <- function(ref_pdf,
                               blocker_cdf,
                               t,
                               rel_tol = 1e-10,
                               abs_tol = 0,
                               subdivisions = 400L) {
  if (!is.finite(t) || t <= 0.0) {
    return(0.0)
  }
  integrate_scalar(
    function(u) inhibit_density_scalar(ref_pdf, blocker_cdf, u),
    0.0,
    t,
    rel_tol = rel_tol,
    abs_tol = abs_tol,
    subdivisions = subdivisions
  )
}

engine_loglik <- function(structure, params, data_df, min_ll = -1e12) {
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  params_df <- build_param_matrix(
    validation_spec_for_params(structure$model_spec, params),
    params,
    trial_df = prepared
  )
  as.numeric(log_likelihood(ctx, prepared, params_df, min_ll = min_ll))
}

engine_density_or_mass <- function(structure, params, data_df, min_ll = -1e12) {
  exp(engine_loglik(structure, params, data_df, min_ll = min_ll))
}

check_row <- function(model_id,
                      check_id,
                      engine_value,
                      manual_value,
                      tolerance,
                      description) {
  data.frame(
    model_id = model_id,
    check_id = check_id,
    description = description,
    engine = engine_value,
    manual = manual_value,
    abs_diff = abs(engine_value - manual_value),
    tolerance = tolerance,
    passed = is.finite(engine_value) &&
      is.finite(manual_value) &&
      abs(engine_value - manual_value) <= tolerance,
    stringsAsFactors = FALSE
  )
}

shared_gate_pair_density_r <- function(eval_x,
                                       eval_y,
                                       eval_c,
                                       t,
                                       rel_tol = 1e-8,
                                       abs_tol = 1e-10) {
  if (!is.finite(t) || t < 0.0) {
    return(0.0)
  }

  ex <- eval_x(t)
  ey <- eval_y(t)
  ec <- eval_c(t)

  f_x <- ex$density
  f_c <- ec$density
  s_y <- min(max(ey$survival, 0.0), 1.0)
  f_y <- min(max(1.0 - s_y, 0.0), 1.0)
  s_c <- min(max(ec$survival, 0.0), 1.0)
  f_cdf <- min(max(1.0 - s_c, 0.0), 1.0)
  s_x <- min(max(ex$survival, 0.0), 1.0)
  f_xcdf <- min(max(1.0 - s_x, 0.0), 1.0)

  term1 <- if (is.finite(f_x) && f_x > 0.0) f_x * f_cdf * s_y else 0.0
  term_c_other <- if (is.finite(f_c) && f_c > 0.0) f_c * f_xcdf * s_y else 0.0
  term2 <- 0.0

  denom <- f_xcdf * f_y
  if (is.finite(f_c) && f_c > 0.0 && is.finite(denom) && denom > 0.0 && t > 0.0) {
    integral <- integrate(
      function(u) {
        vapply(
          u,
          function(uu) {
            ex_u <- eval_x(uu)
            ey_u <- eval_y(uu)
            fx_u <- ex_u$density
            if (!is.finite(fx_u) || fx_u <= 0.0) {
              return(0.0)
            }
            fy_u <- min(max(1.0 - min(max(ey_u$survival, 0.0), 1.0), 0.0), 1.0)
            out <- fx_u * fy_u
            if (!is.finite(out) || out <= 0.0) {
              return(0.0)
            }
            out
          },
          numeric(1)
        )
      },
      lower = 0.0,
      upper = t,
      rel.tol = rel_tol,
      abs.tol = abs_tol,
      subdivisions = 200L
    )$value
    term <- denom - integral
    if (!is.finite(term) || term < 0.0) {
      term <- 0.0
    }
    term2 <- f_c * term
  }

  out <- term1 + term_c_other + term2
  if (!is.finite(out) || out < 0.0) {
    return(0.0)
  }
  out
}

shared_gate_pair_probability_r <- function(eval_x,
                                           eval_y,
                                           eval_c,
                                           upper,
                                           rel_tol = 1e-8,
                                           abs_tol = 1e-10) {
  if (!is.finite(upper)) {
    upper <- Inf
  }
  if (upper <= 0.0) {
    return(0.0)
  }

  coeff <- 1.0
  if (is.finite(upper)) {
    ec_upper <- eval_c(upper)
    coeff <- min(max(1.0 - min(max(ec_upper$survival, 0.0), 1.0), 0.0), 1.0)
  }

  total <- integrate(
    function(tt) {
      vapply(
        tt,
        function(t) {
          ex <- eval_x(t)
          ec <- eval_c(t)

          s_x <- min(max(ex$survival, 0.0), 1.0)
          f_xcdf <- min(max(1.0 - s_x, 0.0), 1.0)
          f_c_dens <- ec$density

          out <- 0.0
          if (is.finite(f_c_dens) && f_c_dens > 0.0 && f_xcdf > 0.0) {
            out <- out + f_c_dens * f_xcdf
          }

          f_x_dens <- ex$density
          if (is.finite(f_x_dens) && f_x_dens > 0.0) {
            ey <- eval_y(t)
            f_ycdf <- min(max(1.0 - min(max(ey$survival, 0.0), 1.0), 0.0), 1.0)
            f_ccdf <- min(max(1.0 - min(max(ec$survival, 0.0), 1.0), 0.0), 1.0)
            out <- out + f_x_dens * (f_ccdf - coeff * f_ycdf)
          }
          if (!is.finite(out)) {
            return(0.0)
          }
          out
        },
        numeric(1)
      )
    },
    lower = 0.0,
    upper = upper,
    rel.tol = rel_tol,
    abs.tol = abs_tol,
    subdivisions = 200L,
    stop.on.error = FALSE
  )$value

  if (!is.finite(total) || total < 0.0) {
    return(0.0)
  }
  min(max(total, 0.0), 1.0)
}

shared_gate_many_density_r <- function(eval_x,
                                       eval_others,
                                       eval_gate,
                                       t,
                                       rel_tol = 1e-8,
                                       abs_tol = 1e-10) {
  if (!is.finite(t) || t < 0.0) {
    return(0.0)
  }
  ex <- eval_x(t)
  eg <- eval_gate(t)
  other_survival_at_t <- prod(vapply(
    eval_others,
    function(eval_other) min(max(eval_other(t)$survival, 0.0), 1.0),
    numeric(1)
  ))
  gate_cdf_at_t <- min(max(1.0 - min(max(eg$survival, 0.0), 1.0), 0.0), 1.0)
  active_target <- ex$density * gate_cdf_at_t * other_survival_at_t

  active_gate <- 0.0
  if (eg$density > 0.0 && t > 0.0) {
    active_gate <- eg$density * integrate_scalar(
      function(u) {
        eu <- eval_x(u)
        if (!is.finite(eu$density) || eu$density <= 0.0) {
          return(0.0)
        }
        other_survival <- prod(vapply(
          eval_others,
          function(eval_other) min(max(eval_other(u)$survival, 0.0), 1.0),
          numeric(1)
        ))
        eu$density * other_survival
      },
      0.0,
      t,
      rel_tol = rel_tol,
      abs_tol = abs_tol,
      subdivisions = 300L
    )
  }

  out <- active_target + active_gate
  if (!is.finite(out) || out < 0.0) {
    return(0.0)
  }
  out
}

oracle_source <- function(name) {
  list(kind = "source", name = name)
}

oracle_pool <- function(members, k = 1L) {
  list(kind = "pool", members = members, k = as.integer(k))
}

oracle_all <- function(...) {
  list(kind = "all", children = list(...))
}

oracle_first <- function(...) {
  list(kind = "first", children = list(...))
}

oracle_inhibit <- function(ref, by) {
  list(kind = "inhibit", ref = ref, blocker = by)
}

oracle_none <- function(expr) {
  list(kind = "none", child = expr)
}

oracle_no_event <- function() {
  list(ok = FALSE, time = Inf, active = NA_character_, readiness = Inf)
}

oracle_event <- function(time, active, readiness = 0.0) {
  list(ok = TRUE, time = time, active = active, readiness = readiness)
}

oracle_active_sources <- function(expr) {
  if (expr$kind == "source") {
    return(expr$name)
  }
  if (expr$kind == "pool") {
    return(expr$members)
  }
  if (expr$kind == "none") {
    return(character())
  }
  if (expr$kind == "inhibit") {
    return(oracle_active_sources(expr$ref))
  }
  if (expr$kind == "all" || expr$kind == "first") {
    out <- unlist(lapply(expr$children, oracle_active_sources), use.names = FALSE)
    return(unique(out))
  }
  stop("unsupported oracle expression kind: ", expr$kind, call. = FALSE)
}

oracle_eval_transition <- function(expr, times, eps = 1e-12) {
  if (expr$kind == "source") {
    return(oracle_event(unname(times[[expr$name]]), expr$name, 0.0))
  }

  if (expr$kind == "pool") {
    member_times <- unname(times[expr$members])
    if (any(!is.finite(member_times)) || length(member_times) < expr$k) {
      return(oracle_no_event())
    }
    ord <- order(member_times, expr$members)
    active_pos <- ord[[expr$k]]
    readiness <- if (expr$k <= 1L) 0.0 else member_times[[ord[[expr$k - 1L]]]]
    return(oracle_event(
      member_times[[active_pos]],
      expr$members[[active_pos]],
      readiness
    ))
  }

  if (expr$kind == "none") {
    return(oracle_no_event())
  }

  if (expr$kind == "inhibit") {
    ref <- oracle_eval_transition(expr$ref, times, eps)
    if (!ref$ok) {
      return(oracle_no_event())
    }
    blocker <- oracle_eval_transition(expr$blocker, times, eps)
    if (blocker$ok && blocker$time <= ref$time + eps) {
      return(oracle_no_event())
    }
    return(ref)
  }

  if (expr$kind == "first") {
    child_events <- lapply(expr$children, oracle_eval_transition, times = times, eps = eps)
    child_events <- child_events[vapply(child_events, function(x) x$ok, logical(1))]
    if (length(child_events) == 0L) {
      return(oracle_no_event())
    }
    child_times <- vapply(child_events, function(x) x$time, numeric(1))
    transition <- min(child_times)
    candidates <- child_events[abs(child_times - transition) <= eps]
    readiness <- vapply(candidates, function(x) x$readiness, numeric(1))
    chosen <- candidates[[which.min(readiness)]]
    return(oracle_event(transition, chosen$active, chosen$readiness))
  }

  if (expr$kind == "all") {
    absent <- Filter(function(x) identical(x$kind, "none"), expr$children)
    positive <- Filter(function(x) !identical(x$kind, "none"), expr$children)
    if (length(positive) == 0L) {
      return(oracle_no_event())
    }
    child_events <- lapply(positive, oracle_eval_transition, times = times, eps = eps)
    if (!all(vapply(child_events, function(x) x$ok, logical(1)))) {
      return(oracle_no_event())
    }
    child_times <- vapply(child_events, function(x) x$time, numeric(1))
    transition <- max(child_times)
    for (none_expr in absent) {
      blocked <- oracle_eval_transition(none_expr$child, times, eps)
      if (blocked$ok && blocked$time <= transition + eps) {
        return(oracle_no_event())
      }
    }
    candidates <- which(abs(child_times - transition) <= eps)
    candidate_readiness <- vapply(candidates, function(idx) {
      other_times <- child_times[-idx]
      max(c(child_events[[idx]]$readiness, other_times, 0.0))
    }, numeric(1))
    chosen_idx <- candidates[[which.min(candidate_readiness)]]
    return(oracle_event(
      transition,
      child_events[[chosen_idx]]$active,
      candidate_readiness[[which.min(candidate_readiness)]]
    ))
  }

  stop("unsupported oracle expression kind: ", expr$kind, call. = FALSE)
}

oracle_gauss_legendre <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (exists(key, cache, inherits = FALSE)) {
      return(get(key, cache, inherits = FALSE))
    }
    i <- seq_len(n - 1L)
    beta <- i / sqrt(4.0 * i * i - 1.0)
    jacobi <- matrix(0.0, n, n)
    jacobi[cbind(i, i + 1L)] <- beta
    jacobi[cbind(i + 1L, i)] <- beta
    eig <- eigen(jacobi, symmetric = TRUE)
    ord <- order(eig$values)
    nodes <- (eig$values[ord] + 1.0) / 2.0
    weights <- eig$vectors[1L, ord]^2
    out <- list(nodes = nodes, weights = weights)
    assign(key, out, cache)
    out
  }
})

oracle_source_quantile <- function(u, p) {
  live_prob <- 1.0 - p$q
  if (u >= live_prob) {
    return(Inf)
  }
  p$t0 + qlnorm(u / live_prob, meanlog = p$m, sdlog = p$s)
}

oracle_integrate_source_times <- function(source_names,
                                          source_parts,
                                          fixed_times,
                                          fn,
                                          nodes = 7L) {
  varying <- setdiff(source_names, names(fixed_times))
  rule <- oracle_gauss_legendre(nodes)
  times <- rep(NA_real_, length(source_names))
  names(times) <- source_names
  for (name in names(fixed_times)) {
    times[[name]] <- fixed_times[[name]]
  }

  walk <- function(pos, weight) {
    if (pos > length(varying)) {
      return(weight * fn(times))
    }
    name <- varying[[pos]]
    p <- source_parts[[name]]
    total <- 0.0
    for (i in seq_along(rule$nodes)) {
      times[[name]] <<- oracle_source_quantile(rule$nodes[[i]], p)
      total <- total + walk(pos + 1L, weight * rule$weights[[i]])
    }
    total
  }

  walk(1L, 1.0)
}

oracle_target_wins <- function(target_event,
                               competitor_event,
                               eps = 1e-10) {
  if (!competitor_event$ok) {
    return(TRUE)
  }
  if (competitor_event$time < target_event$time - eps) {
    return(FALSE)
  }
  if (competitor_event$time > target_event$time + eps) {
    return(TRUE)
  }
  if (!identical(competitor_event$active, target_event$active)) {
    return(FALSE)
  }
  target_event$readiness < competitor_event$readiness - eps
}

oracle_outcome_density <- function(outcomes,
                                   target,
                                   source_names,
                                   params,
                                   rt,
                                   nodes = 7L) {
  source_parts <- stats::setNames(
    lapply(source_names, acc_parts, params = params),
    source_names
  )
  total <- 0.0
  target_sources <- intersect(oracle_active_sources(outcomes[[target]]), source_names)
  for (active_source in target_sources) {
    density <- acc_pdf_scalar(rt, source_parts[[active_source]])
    if (!is.finite(density) || density <= 0.0) {
      next
    }
    win_probability <- oracle_integrate_source_times(
      source_names,
      source_parts,
      stats::setNames(list(rt), active_source),
      function(times) {
        target_event <- oracle_eval_transition(outcomes[[target]], times)
        if (!target_event$ok ||
            !identical(target_event$active, active_source) ||
            abs(target_event$time - rt) > 1e-10) {
          return(0.0)
        }
        for (name in setdiff(names(outcomes), target)) {
          competitor_event <- oracle_eval_transition(outcomes[[name]], times)
          if (!oracle_target_wins(target_event, competitor_event)) {
            return(0.0)
          }
        }
        1.0
      },
      nodes = nodes
    )
    total <- total + density * win_probability
  }
  total
}

oracle_density_rows <- function(model_id,
                                structure,
                                params,
                                outcomes,
                                source_names,
                                response,
                                rts,
                                tolerance,
                                description,
                                nodes = 7L,
                                component = NULL) {
  context <- make_context(structure)
  engine_for_data <- function(data_df) {
    prepared <- prepare_data(structure, data_df)
    params_df <- build_param_matrix(
      validation_spec_for_params(structure$model_spec, params),
      params,
      trial_df = prepared
    )
    exp(as.numeric(log_likelihood(context, prepared, params_df)))
  }
  rows <- vector("list", length(rts))
  for (i in seq_along(rts)) {
    rt <- rts[[i]]
    data_df <- data.frame(
      trials = 1L,
      R = response,
      rt = rt,
      stringsAsFactors = FALSE
    )
    if (!is.null(component)) {
      data_df$component <- component
    }
    engine <- engine_for_data(data_df)
    manual <- oracle_outcome_density(
      outcomes, response, source_names, params, rt, nodes = nodes)
    rows[[i]] <- check_row(
      model_id,
      paste0(response, "_rt_", format(rt, nsmall = 2)),
      engine,
      manual,
      tolerance,
      description
    )
  }
  do.call(rbind, rows)
}
