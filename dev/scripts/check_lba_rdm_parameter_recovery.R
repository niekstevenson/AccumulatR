#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(AccumulatR)
  }
})

args <- commandArgs(trailingOnly = TRUE)
n_trials <- if (length(args) >= 1L) as.integer(args[[1]]) else 1500L
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260224L

if (!is.finite(n_trials) || n_trials < 200L) {
  stop("n_trials must be >= 200")
}
if (!is.finite(seed)) {
  stop("seed must be finite")
}

set.seed(seed)

sanity_check_distribution <- function(dist, pars) {
  x <- seq(1e-4, 4.0, length.out = 400L)
  if (dist == "lba") {
    pdf <- dist_lba_pdf(x, pars[["v"]], pars[["sv"]], pars[["B"]], pars[["A"]])
    cdf <- dist_lba_cdf(x, pars[["v"]], pars[["sv"]], pars[["B"]], pars[["A"]])
  } else if (dist == "rdm") {
    pdf <- dist_rdm_pdf(x, pars[["v"]], pars[["B"]], pars[["A"]], pars[["s"]])
    cdf <- dist_rdm_cdf(x, pars[["v"]], pars[["B"]], pars[["A"]], pars[["s"]])
  } else {
    stop("Unsupported dist: ", dist)
  }
  if (any(!is.finite(pdf)) || any(pdf < -1e-12)) {
    stop(sprintf("%s sanity check failed: pdf invalid", toupper(dist)))
  }
  if (any(!is.finite(cdf)) || any(cdf < -1e-8) || any(cdf > 1 + 1e-8)) {
    stop(sprintf("%s sanity check failed: cdf out of [0,1]", toupper(dist)))
  }
  if (any(diff(cdf) < -1e-7)) {
    stop(sprintf("%s sanity check failed: cdf not monotone", toupper(dist)))
  }
}

build_single_acc_model <- function(dist) {
  race_spec() |>
    add_accumulator("a", dist) |>
    add_outcome("A", "a") |>
    finalize_model()
}

run_recovery <- function(dist,
                         true_core,
                         t0 = 0.05,
                         n_trials = 1500L,
                         fit_core = names(true_core),
                         fixed_core = c(),
                         n_starts = 4L,
                         maxit = 1200L,
                         rel_tol = 0.35) {
  spec <- build_single_acc_model(dist)
  core_nm <- names(true_core)
  true_core <- as.numeric(true_core)
  names(true_core) <- core_nm
  fit_core <- as.character(fit_core)
  if (!all(fit_core %in% names(true_core))) {
    stop("fit_core contains unknown parameter names")
  }
  fixed_nm <- names(fixed_core)
  fixed_core <- as.numeric(fixed_core)
  names(fixed_core) <- fixed_nm
  if (length(fixed_core) > 0L && is.null(names(fixed_core))) {
    stop("fixed_core must be a named numeric vector")
  }
  if (length(fixed_core) > 0L && any(!names(fixed_core) %in% names(true_core))) {
    stop("fixed_core contains unknown parameter names")
  }
  if (length(intersect(fit_core, names(fixed_core))) > 0L) {
    stop("fit_core and fixed_core overlap")
  }
  core_names <- names(true_core)
  full_true <- c(stats::setNames(true_core, paste0("a.", core_names)), a.q = 0, a.t0 = t0)

  params_true <- build_param_matrix(spec, full_true, n_trials = n_trials)
  sim <- simulate(spec, params_true, seed = sample.int(.Machine$integer.max, 1L))
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )

  ctx <- build_likelihood_context(spec, data_df)

  objective <- function(theta_log) {
    fit_vals <- exp(theta_log)
    names(fit_vals) <- fit_core
    all_core_vals <- true_core
    all_core_vals[fit_core] <- fit_vals
    if (length(fixed_core) > 0L) {
      all_core_vals[names(fixed_core)] <- fixed_core
    }
    full_vals <- c(stats::setNames(all_core_vals, paste0("a.", core_names)), a.q = 0, a.t0 = t0)
    pm <- build_param_matrix(
      spec,
      full_vals,
      n_trials = n_trials,
      layout = ctx$param_layout
    )
    ll <- as.numeric(log_likelihood(ctx, pm))
    if (!is.finite(ll)) {
      return(1e12)
    }
    -ll
  }

  starts <- vector("list", n_starts)
  starts[[1L]] <- log(unname(true_core[fit_core])) + stats::rnorm(length(fit_core), sd = 0.20)
  if (n_starts > 1L) {
    for (i in 2:n_starts) {
      starts[[i]] <- log(unname(true_core[fit_core])) + stats::rnorm(length(fit_core), sd = 0.45)
    }
  }

  fits <- lapply(starts, function(st) {
    stats::optim(
      par = st,
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = 1e-9)
    )
  })

  vals <- vapply(fits, `[[`, numeric(1), "value")
  best_idx <- which.min(vals)
  fit <- fits[[best_idx]]

  est_fit <- exp(fit$par)
  names(est_fit) <- fit_core

  ll_true <- -objective(log(unname(true_core[fit_core])))
  ll_est <- -fit$value

  rel_err <- abs(est_fit - true_core[fit_core]) / pmax(abs(true_core[fit_core]), 1e-8)
  out_tbl <- data.frame(
    parameter = fit_core,
    true = as.numeric(true_core[fit_core]),
    recovered = as.numeric(est_fit),
    abs_error = as.numeric(abs(est_fit - true_core[fit_core])),
    rel_error = as.numeric(rel_err),
    stringsAsFactors = FALSE
  )

  pass <- is.finite(ll_est) &&
    is.finite(ll_true) &&
    fit$convergence == 0L &&
    max(rel_err) <= rel_tol

  list(
    dist = dist,
    pass = pass,
    fit = fit,
    ll_true = ll_true,
    ll_est = ll_est,
    fit_core = fit_core,
    fixed_core = fixed_core,
    table = out_tbl
  )
}

lba_true <- c(v = 2.0, sv = 0.7, B = 1.15, A = 0.35)
rdm_true <- c(v = 1.80, B = 1.20, A = 0.80, s = 1.10)

sanity_check_distribution("lba", lba_true)
sanity_check_distribution("rdm", rdm_true)

lba_res <- run_recovery(
  dist = "lba",
  true_core = lba_true,
  t0 = 0.06,
  n_trials = n_trials,
  fit_core = c("v", "B", "A"),
  fixed_core = c(sv = lba_true[["sv"]])
)
rdm_res <- run_recovery(
  dist = "rdm",
  true_core = rdm_true,
  t0 = 0.03,
  n_trials = n_trials,
  fit_core = c("v", "B", "A"),
  fixed_core = c(s = rdm_true[["s"]]),
  rel_tol = 0.60
)

print_result <- function(res) {
  cat("\n==============================\n")
  cat(sprintf("Distribution: %s\n", toupper(res$dist)))
  if (length(res$fixed_core) > 0L) {
    cat("Fixed during fit:\n")
    print(res$fixed_core)
  }
  cat(sprintf("Convergence code: %d\n", res$fit$convergence))
  cat(sprintf("LogLik(true): %.6f\n", res$ll_true))
  cat(sprintf("LogLik(best): %.6f\n", res$ll_est))
  print(res$table, row.names = FALSE)
  cat(sprintf("PASS: %s\n", if (isTRUE(res$pass)) "YES" else "NO"))
}

print_result(lba_res)
print_result(rdm_res)

if (!lba_res$pass || !rdm_res$pass) {
  stop("Parameter recovery check failed for one or more distributions")
}

cat("\nAll LBA/RDM recovery checks passed.\n")
