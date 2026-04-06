#!/usr/bin/env Rscript

pkgload::load_all(".", quiet = TRUE)

dir.create("vignettes/data", recursive = TRUE, showWarnings = FALSE)

save_simple_model <- function() {
  model_spec <- race_spec() |>
    add_accumulator("R1_A", "LBA") |>
    add_accumulator("R1_B", "RDM") |>
    add_accumulator("R2", "lognormal") |>
    add_pool("R1", c("R1_A", "R1_B")) |>
    add_outcome("R1", "R1") |>
    add_outcome("R2", "R2")

  structure <- finalize_model(model_spec)

  true_params <- c(
    R1_A.v = 2,
    R1_A.B = 1,
    R1_A.A = 0.3,
    R1_A.sv = 1,
    R1_A.q = 0,
    R1_A.t0 = 0,
    R1_B.v = 3,
    R1_B.B = 1,
    R1_B.A = 0.3,
    R1_B.s = 1,
    R1_B.q = 0,
    R1_B.t0 = 0,
    R2.m = log(0.4),
    R2.s = 0.18,
    R2.q = 0,
    R2.t0 = 0
  )

  set.seed(123456)
  params_df <- build_param_matrix(model_spec, true_params, n_trials = 2000)
  sim <- simulate(structure, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(structure, data_df)

  neg_loglik <- function(theta) {
    est <- true_params
    est["R1_A.v"] <- theta[["R1_A.v"]]
    est["R1_B.v"] <- theta[["R1_B.v"]]
    est["R1_A.B"] <- exp(theta[["log_B_shared"]])
    est["R1_B.B"] <- exp(theta[["log_B_shared"]])
    est["R2.m"] <- theta[["R2.m"]]
    est["R2.s"] <- exp(theta[["log_R2.s"]])
    params_df <- build_param_matrix(
      model_spec,
      est,
      n_trials = max(data_df$trial),
      layout = ctx$param_layout
    )
    -as.numeric(log_likelihood(ctx, params_df))
  }

  start <- c(
    R1_A.v = 1.5,
    R1_B.v = 1.5,
    log_B_shared = log(1.0),
    R2.m = log(0.32),
    log_R2.s = log(0.15)
  )

  set.seed(123456)
  fit <- optim(
    start,
    neg_loglik,
    method = "Nelder-Mead",
    control = list(maxit = 4000, reltol = 1e-9)
  )

  save(fit, file = "vignettes/data/simple_model.RData")
}

save_multi_outcome <- function() {
  model_spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B")

  structure <- finalize_model(model_spec)

  true_params <- c(
    A.m = log(0.30),
    A.s = 0.18,
    A.q = 0,
    A.t0 = 0,
    B.m = log(0.38),
    B.s = 0.22,
    B.q = 0,
    B.t0 = 0
  )

  set.seed(123456)
  params_df <- build_param_matrix(model_spec, true_params, n_trials = 500)
  sim <- simulate(structure, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    R2 = factor(sim$R2),
    rt2 = sim$rt2,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(structure, data_df)

  neg_loglik <- function(theta) {
    est <- true_params
    est[c("A.m", "A.s", "B.m", "B.s")] <- theta[c("A.m", "A.s", "B.m", "B.s")]
    est[c("A.s", "B.s")] <- exp(est[c("A.s", "B.s")])
    params_df <- build_param_matrix(
      model_spec,
      est,
      n_trials = max(data_df$trial),
      layout = ctx$param_layout
    )
    -as.numeric(log_likelihood(ctx, params_df))
  }

  start <- c(
    A.m = log(0.24),
    A.s = log(0.12),
    B.m = log(0.48),
    B.s = log(0.12)
  )

  set.seed(123456)
  fit <- optim(start, neg_loglik, method = "Nelder-Mead")

  save(fit, file = "vignettes/data/multi_outcome.RData")
}

save_chained_onset <- function() {
  model_spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("C", "lognormal", onset = after("B")) |>
    add_outcome("A", "A") |>
    add_outcome("C", "C")

  structure <- finalize_model(model_spec)

  true_params <- c(
    A.m = log(0.28),
    A.s = 0.14,
    A.q = 0,
    A.t0 = 0,
    B.m = log(0.1),
    B.s = 0.1,
    B.q = 0,
    B.t0 = 0,
    C.m = log(0.15),
    C.s = 0.1,
    C.q = 0,
    C.t0 = 0
  )

  set.seed(123456)
  params_df <- build_param_matrix(model_spec, true_params, n_trials = 2000)
  sim <- simulate(structure, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(structure, data_df)

  neg_loglik <- function(theta) {
    est <- true_params
    est[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")] <- theta[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")]
    est[c("A.s", "B.s", "C.s")] <- exp(est[c("A.s", "B.s", "C.s")])
    params_df <- build_param_matrix(
      model_spec,
      est,
      n_trials = max(data_df$trial),
      layout = ctx$param_layout
    )
    -as.numeric(log_likelihood(ctx, params_df))
  }

  start <- c(
    A.m = log(0.22),
    A.s = log(0.10),
    B.m = log(0.28),
    B.s = log(0.10),
    C.m = log(0.28),
    C.s = log(0.10)
  )

  set.seed(123456)
  fit <- optim(start, neg_loglik, method = "Nelder-Mead")

  save(fit, file = "vignettes/data/chained_onset.RData")
}

save_simple_model()
save_multi_outcome()
save_chained_onset()
