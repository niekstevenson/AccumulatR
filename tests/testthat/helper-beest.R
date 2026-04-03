# Standalone likelihood for the specific BEEST-style stop/change model:
#
#   model <- race_spec() |>
#     add_accumulator("S", "exgauss") |>
#     add_accumulator("stop", "exgauss") |>
#     add_accumulator("change", "exgauss") |>
#     add_outcome("S", inhibit("S", by = "stop")) |>
#     add_outcome("X", all_of("change", "stop")) |>
#     add_component("go_only", members = "S", weight = .75) |>
#     add_component("go_stop", members = c("S", "stop", "change"), weight = .25) |>
#     add_trigger(
#       "stop_trigger",
#       members = c("stop", "change"),
#       q = 0.05,
#       param = "stop_trigger"
#     ) |>
#     set_mixture_options(mode = "fixed") |>
#     finalize_model()
#
# Assumptions:
# - The likelihood is hard-coded to the model above.
# - Input data are in the long format shown in the prompt, with one row per
#   active accumulator per trial.
# - Parameter vector order:
#     c(
#       S.mu, S.sigma, S.tau,
#       stop.mu, stop.sigma, stop.tau,
#       change.mu, change.sigma, change.tau,
#       stop_trigger
#     )
# - `onset` is treated as a hard lower bound, with ex-Gaussian latent times
#   truncated to `x > 0` and renormalized (no onset atom).

.beest_clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

.beest_exgauss_pdf <- function(x, mu, sigma, tau) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  if (!is.finite(mu) || !is.finite(sigma) || !is.finite(tau) ||
      sigma <= 0 || tau <= 0) {
    return(out)
  }

  na_idx <- is.na(x)
  out[na_idx] <- NA_real_

  finite_idx <- which(!na_idx & is.finite(x))
  if (length(finite_idx) > 0L) {
    xx <- x[finite_idx]
    z <- (xx - mu) / sigma
    log_pdf <- -log(tau) +
      (sigma^2) / (2 * tau^2) -
      (xx - mu) / tau +
      stats::pnorm(z - sigma / tau, log.p = TRUE)
    val <- exp(log_pdf)
    val[!is.finite(val) | val < 0] <- 0
    out[finite_idx] <- val
  }

  out[!na_idx & !is.finite(x)] <- 0
  out
}

.beest_exgauss_cdf <- function(x, mu, sigma, tau) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  if (!is.finite(mu) || !is.finite(sigma) || !is.finite(tau) ||
      sigma <= 0 || tau <= 0) {
    return(out)
  }

  na_idx <- is.na(x)
  out[na_idx] <- NA_real_

  finite_idx <- which(!na_idx & is.finite(x))
  if (length(finite_idx) > 0L) {
    xx <- x[finite_idx]
    z <- (xx - mu) / sigma
    correction <- exp(
      (sigma^2) / (2 * tau^2) -
        (xx - mu) / tau +
        stats::pnorm(z - sigma / tau, log.p = TRUE)
    )
    val <- stats::pnorm(z) - correction
    val[!is.finite(val)] <- NA_real_
    out[finite_idx] <- .beest_clamp01(val)
  }

  nonfinite_idx <- !na_idx & !is.finite(x)
  if (any(nonfinite_idx)) {
    out[nonfinite_idx & x < 0] <- 0
    out[nonfinite_idx & x > 0] <- 1
  }

  out
}

.beest_exgauss_pos_mass <- function(mu, sigma, tau) {
  cdf0 <- .beest_exgauss_cdf(0, mu, sigma, tau)
  out <- 1 - cdf0
  if (!is.finite(out) || out <= 0) {
    return(NA_real_)
  }
  out
}

.beest_exgauss_trunc_pdf <- function(x, mu, sigma, tau) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  norm_const <- .beest_exgauss_pos_mass(mu, sigma, tau)
  if (!is.finite(norm_const)) {
    return(out)
  }

  na_idx <- is.na(x)
  out[na_idx] <- NA_real_

  before_zero <- !na_idx & (is.infinite(x) | x < 0)
  out[before_zero] <- 0

  ok_idx <- which(!na_idx & is.finite(x) & x >= 0)
  if (length(ok_idx) > 0L) {
    out[ok_idx] <- .beest_exgauss_pdf(x[ok_idx], mu, sigma, tau) / norm_const
    bad_idx <- ok_idx[!is.finite(out[ok_idx]) | out[ok_idx] < 0]
    out[bad_idx] <- 0
  }

  out
}

.beest_exgauss_trunc_cdf <- function(x, mu, sigma, tau) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  norm_const <- .beest_exgauss_pos_mass(mu, sigma, tau)
  cdf0 <- .beest_exgauss_cdf(0, mu, sigma, tau)
  if (!is.finite(norm_const) || !is.finite(cdf0)) {
    return(out)
  }

  na_idx <- is.na(x)
  out[na_idx] <- NA_real_

  before_zero <- !na_idx & (x < 0)
  out[before_zero] <- 0

  ok_idx <- which(!na_idx & is.finite(x) & x >= 0)
  if (length(ok_idx) > 0L) {
    numer <- .beest_exgauss_cdf(x[ok_idx], mu, sigma, tau) - cdf0
    out[ok_idx] <- .beest_clamp01(numer / norm_const)
  }

  pos_inf_idx <- !na_idx & is.infinite(x) & x > 0
  out[pos_inf_idx] <- 1

  .beest_clamp01(out)
}

.beest_shifted_pdf <- function(t, mu, sigma, tau, onset = 0) {
  t <- as.numeric(t)
  out <- rep(NA_real_, length(t))

  na_idx <- is.na(t)
  out[na_idx] <- NA_real_

  before_onset <- !na_idx & (is.infinite(t) | t < 0 | t < onset)
  out[before_onset] <- 0

  ok_idx <- which(!na_idx & is.finite(t) & t >= 0 & t >= onset)
  if (length(ok_idx) > 0L) {
    out[ok_idx] <- .beest_exgauss_trunc_pdf(t[ok_idx] - onset, mu, sigma, tau)
    bad_idx <- ok_idx[!is.finite(out[ok_idx]) | out[ok_idx] < 0]
    out[bad_idx] <- 0
  }

  out
}

.beest_shifted_cdf <- function(t, mu, sigma, tau, onset = 0) {
  t <- as.numeric(t)
  out <- rep(NA_real_, length(t))

  na_idx <- is.na(t)
  out[na_idx] <- NA_real_

  before_onset <- !na_idx & (t < 0 | t < onset)
  out[before_onset] <- 0

  ok_idx <- which(!na_idx & is.finite(t) & t >= 0 & t >= onset)
  if (length(ok_idx) > 0L) {
    out[ok_idx] <- .beest_exgauss_trunc_cdf(t[ok_idx] - onset, mu, sigma, tau)
  }

  pos_inf_idx <- !na_idx & is.infinite(t) & t > 0
  out[pos_inf_idx] <- 1

  .beest_clamp01(out)
}

.beest_shifted_survival <- function(t, mu, sigma, tau, onset = 0) {
  1 - .beest_shifted_cdf(t, mu, sigma, tau, onset = onset)
}

.beest_parse_par <- function(par) {
  par <- as.numeric(par)
  if (length(par) != 10L) {
    stop("`par` must be a numeric vector of length 10.", call. = FALSE)
  }
  if (any(!is.finite(par))) {
    stop("All 10 parameters must be finite.", call. = FALSE)
  }

  names(par) <- c(
    "S.mu", "S.sigma", "S.tau",
    "stop.mu", "stop.sigma", "stop.tau",
    "change.mu", "change.sigma", "change.tau",
    "stop_trigger"
  )

  if (par[["S.sigma"]] <= 0 || par[["stop.sigma"]] <= 0 || par[["change.sigma"]] <= 0) {
    stop("All sigma parameters must be > 0.", call. = FALSE)
  }
  if (par[["S.tau"]] <= 0 || par[["stop.tau"]] <= 0 || par[["change.tau"]] <= 0) {
    stop("All tau parameters must be > 0.", call. = FALSE)
  }
  if (par[["stop_trigger"]] < 0 || par[["stop_trigger"]] > 1) {
    stop("`stop_trigger` must be in [0, 1].", call. = FALSE)
  }

  list(
    S = list(mu = par[["S.mu"]], sigma = par[["S.sigma"]], tau = par[["S.tau"]]),
    stop = list(mu = par[["stop.mu"]], sigma = par[["stop.sigma"]], tau = par[["stop.tau"]]),
    change = list(mu = par[["change.mu"]], sigma = par[["change.sigma"]], tau = par[["change.tau"]]),
    stop_trigger = par[["stop_trigger"]]
  )
}

.beest_trial_key <- function(data) {
  if ("trial" %in% names(data) && "subject" %in% names(data)) {
    return(paste(data$subject, data$trial, sep = "::"))
  }
  if ("trials" %in% names(data) && "subjects" %in% names(data)) {
    return(paste(data$subjects, data$trials, sep = "::"))
  }
  if ("trial" %in% names(data)) {
    return(as.character(data$trial))
  }
  if ("trials" %in% names(data)) {
    return(as.character(data$trials))
  }
  stop("Data must contain `trial` or `trials`.", call. = FALSE)
}

.beest_row_onset <- function(data) {
  if ("onset" %in% names(data)) {
    return(as.numeric(data$onset))
  }
  if ("SSD" %in% names(data)) {
    acc <- as.character(data$accumulator)
    ssd <- as.numeric(data$SSD)
    out <- ifelse(acc == "S", 0, ssd)
    out[is.na(out)] <- 0
    return(out)
  }
  rep(0, nrow(data))
}

.beest_prepare_trials <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  required <- c("R", "rt", "accumulator")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf("Data are missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  data$accumulator <- as.character(data$accumulator)
  data$R <- as.character(data$R)
  data$rt <- as.numeric(data$rt)
  data$onset_internal <- .beest_row_onset(data)
  data$trial_key <- .beest_trial_key(data)

  split_rows <- split(data, data$trial_key, drop = TRUE)

  lapply(split_rows, function(df) {
    R_vals <- unique(df$R[!is.na(df$R)])
    if (length(R_vals) != 1L) {
      stop("Each trial must have a single observed outcome label.", call. = FALSE)
    }

    rt_vals <- unique(df$rt[!is.na(df$rt)])
    if (length(rt_vals) != 1L) {
      stop("Each trial must have a single observed RT.", call. = FALSE)
    }

    accs <- df$accumulator
    if (anyDuplicated(accs)) {
      stop("Each accumulator may appear at most once per trial.", call. = FALSE)
    }

    onset_map <- stats::setNames(df$onset_internal, accs)
    component <- if (setequal(accs, "S")) {
      "go_only"
    } else if (setequal(accs, c("S", "stop", "change"))) {
      "go_stop"
    } else {
      stop(
        sprintf(
          "Trial '%s' does not match either component: accumulators were [%s].",
          df$trial_key[[1]],
          paste(sort(accs), collapse = ", ")
        ),
        call. = FALSE
      )
    }

    list(
      trial_key = df$trial_key[[1]],
      R = R_vals[[1]],
      rt = rt_vals[[1]],
      component = component,
      onset = onset_map
    )
  })
}

.beest_event_pdf <- function(t, acc, onset = 0, q = 0) {
  (1 - q) * .beest_shifted_pdf(t, acc$mu, acc$sigma, acc$tau, onset = onset)
}

.beest_event_cdf <- function(t, acc, onset = 0, q = 0) {
  (1 - q) * .beest_shifted_cdf(t, acc$mu, acc$sigma, acc$tau, onset = onset)
}

.beest_event_survival <- function(t, acc, onset = 0, q = 0) {
  q + (1 - q) * .beest_shifted_survival(t, acc$mu, acc$sigma, acc$tau, onset = onset)
}

.beest_stop_prefix_against_S_upto <- function(t, pars, s_onset, stop_onset,
                                              cache, rel.tol, abs.tol) {
  if (!is.finite(t) || t <= 0) {
    return(0)
  }

  key <- paste(
    "stop_prefix_against_S",
    formatC(t, digits = 16, format = "fg"),
    formatC(s_onset, digits = 16, format = "fg"),
    formatC(stop_onset, digits = 16, format = "fg"),
    sep = "|"
  )

  cached <- cache[[key]]
  if (!is.null(cached)) {
    return(cached)
  }

  lower <- max(0, stop_onset)
  if (!(t > lower)) {
    cache[[key]] <- 0
    return(0)
  }

  integrand <- function(u) {
    .beest_event_pdf(u, pars$stop, onset = stop_onset, q = 0) *
      .beest_event_survival(u, pars$S, onset = s_onset, q = 0)
  }
  out <- stats::integrate(
    integrand,
    lower = lower,
    upper = t,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    stop.on.error = FALSE
  )$value
  if (!is.finite(out) || out < 0) {
    out <- 0
  }
  cache[[key]] <- out
  out
}

.beest_trial_density <- function(trial, pars,
                                 include_component_weight,
                                 component_weights,
                                 cache,
                                 rel.tol,
                                 abs.tol) {
  t <- trial$rt
  R_obs <- trial$R
  s_onset <- unname(trial$onset[["S"]])
  if (is.na(s_onset)) {
    s_onset <- 0
  }

  if (trial$component == "go_only") {
    dens <- switch(
      R_obs,
      "S" = .beest_event_pdf(t, pars$S, onset = s_onset, q = 0),
      "X" = 0,
      0
    )
  } else {
    stop_onset <- unname(trial$onset[["stop"]])
    change_onset <- unname(trial$onset[["change"]])
    if (is.na(stop_onset)) {
      stop_onset <- 0
    }
    if (is.na(change_onset)) {
      change_onset <- 0
    }
    q_fail <- pars$stop_trigger

    f_S <- .beest_event_pdf(t, pars$S, onset = s_onset, q = 0)
    S_S <- .beest_event_survival(t, pars$S, onset = s_onset, q = 0)

    f_stop <- .beest_event_pdf(t, pars$stop, onset = stop_onset, q = 0)
    F_stop <- .beest_event_cdf(t, pars$stop, onset = stop_onset, q = 0)
    S_stop <- .beest_event_survival(t, pars$stop, onset = stop_onset, q = 0)

    f_change <- .beest_event_pdf(t, pars$change, onset = change_onset, q = 0)
    F_change <- .beest_event_cdf(t, pars$change, onset = change_onset, q = 0)

    stop_prefix_against_S <- .beest_stop_prefix_against_S_upto(
      t = t,
      pars = pars,
      s_onset = s_onset,
      stop_onset = stop_onset,
      cache = cache,
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )

    exact_S_density <- f_S * S_stop
    exact_X_density <- f_stop * F_change * S_S +
      f_change * stop_prefix_against_S

    dens <- switch(
      R_obs,
      "S" = q_fail * f_S +
        (1 - q_fail) * exact_S_density,
      "X" = (1 - q_fail) * exact_X_density,
      0
    )
  }

  dens <- as.numeric(dens)
  dens[!is.finite(dens) | dens < 0] <- 0

  if (isTRUE(include_component_weight)) {
    dens <- dens * component_weights[[trial$component]]
  }

  dens
}

# Returns a function that takes only the 10-parameter vector.
make_beest_loglik <- function(data,
                              include_component_weight = TRUE,
                              component_weights = c(go_only = 0.75, go_stop = 0.25),
                              return_log = TRUE,
                              rel.tol = 1e-8,
                              abs.tol = 0) {
  trials <- .beest_prepare_trials(data)

  if (length(component_weights) != 2L) {
    stop("`component_weights` must contain `go_only` and `go_stop`.", call. = FALSE)
  }
  if (is.null(names(component_weights))) {
    names(component_weights) <- c("go_only", "go_stop")
  }
  component_weights <- component_weights[c("go_only", "go_stop")]
  if (any(is.na(component_weights))) {
    stop("`component_weights` must be named `go_only` and `go_stop`.", call. = FALSE)
  }
  component_weights <- as.numeric(component_weights)
  names(component_weights) <- c("go_only", "go_stop")

  function(par) {
    pars <- .beest_parse_par(par)
    cache <- new.env(parent = emptyenv(), hash = TRUE)

    trial_dens <- vapply(
      trials,
      .beest_trial_density,
      numeric(1),
      pars = pars,
      include_component_weight = include_component_weight,
      component_weights = component_weights,
      cache = cache,
      rel.tol = rel.tol,
      abs.tol = abs.tol
    )

    if (isTRUE(return_log)) {
      if (any(!is.finite(trial_dens) | trial_dens <= 0)) {
        return(-Inf)
      }
      return(sum(log(trial_dens)))
    }

    if (any(!is.finite(trial_dens) | trial_dens < 0)) {
      return(0)
    }
    prod(trial_dens)
  }
}

# Convenience wrapper if you do not want a closure.
beest_loglik <- function(par, data, ...) {
  make_beest_loglik(data, ...)(par)
}

beest_nll <- function(par, data, ...) {
  -beest_loglik(par, data, ...)
}

# Example:
# ll_fun <- make_beest_loglik(emcs[[1]][[1]]$data$`8953`)
# par <- c(
#   0.35, 0.05, 0.08,
#   0.15, 0.04, 0.06,
#   0.25, 0.05, 0.09,
#   0.05
# )
# ll_fun(par)
