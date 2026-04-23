`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

validation_source_rebuild <- function(repo_root) {
  source(file.path(repo_root, "R", "helpers.R"))
  source(file.path(repo_root, "R", "model_definition.R"))
  source(file.path(repo_root, "R", "likelihood_prep.R"))
  source(file.path(repo_root, "R", "likelihood_common.R"))
  source(file.path(repo_root, "R", "semantic_bridge.R"))
  source(file.path(repo_root, "R", "likelihood_param_interface.R"))
}

acc_parts <- function(prefix, params) {
  q_name <- paste0(prefix, ".q")
  t0_name <- paste0(prefix, ".t0")
  list(
    m = unname(params[[paste0(prefix, ".m")]]),
    s = unname(params[[paste0(prefix, ".s")]]),
    q = if (q_name %in% names(params)) unname(params[[q_name]]) else 0.0,
    t0 = if (t0_name %in% names(params)) unname(params[[t0_name]]) else 0.0
  )
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
  params_df <- build_param_matrix(structure$model_spec, params, trial_df = prepared)
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
