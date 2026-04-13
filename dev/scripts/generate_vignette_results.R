#!/usr/bin/env Rscript

pkgload::load_all(".", quiet = TRUE)

dir.create("vignettes/data", recursive = TRUE, showWarnings = FALSE)

save_simple_model <- function() {
  model <- race_spec() |>
    add_accumulator("R1_A", "LBA") |>
    add_accumulator("R1_B", "RDM") |>
    add_accumulator("R2", "lognormal") |>
    add_pool("R1", c("R1_A", "R1_B")) |>
    add_outcome("R1", "R1") |>
    add_outcome("R2", "R2") |>
    set_parameters(list(
      B_shared = c("R1_A.B", "R1_B.B"),
      A_shared = c("R1_A.A", "R1_B.A"),
      noise_shared = c("R1_A.sv", "R1_B.s"),
      t0_shared = c("R1_A.t0", "R1_B.t0", "R2.t0")
    )) |>
    finalize_model()

  true_params <- c(
    R1_A.v = 2,
    R1_B.v = 3,
    B_shared = 1,
    A_shared = 0.3,
    noise_shared = 1,
    R2.m = log(0.4),
    R2.s = 0.18,
    t0_shared = 0.05
  )

  set.seed(123456)
  params_df <- build_param_matrix(model, true_params, n_trials = 2000)
  sim <- simulate(model, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(model, data_df)
  ctx <- make_context(model)

  neg_loglik <- function(theta) {
    est <- true_params
    est["R1_A.v"] <- theta[["R1_A.v"]]
    est["R1_B.v"] <- theta[["R1_B.v"]]
    est["B_shared"] <- exp(theta[["log_B_shared"]])
    est["R2.m"] <- theta[["R2.m"]]
    est["R2.s"] <- exp(theta[["log_R2.s"]])
    est["t0_shared"] <- exp(theta[["log_t0_shared"]])
    params_df <- build_param_matrix(
      model,
      est,
      trial_df = prepared
    )
    -as.numeric(log_likelihood(ctx, prepared, params_df))
  }

  start <- c(
    R1_A.v = 1.5,
    R1_B.v = 1.5,
    log_B_shared = log(1.0),
    R2.m = log(0.32),
    log_R2.s = log(0.15),
    log_t0_shared = log(0.03)
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

save_pool_model <- function() {
  model_spec <- race_spec() |>
    add_accumulator("A1", "lognormal") |>
    add_accumulator("A2", "lognormal") |>
    add_accumulator("A3", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_pool("A_pool", c("A1", "A2", "A3"), k = 2L) |>
    add_outcome("A", "A_pool") |>
    add_outcome("B", "B")

  structure <- finalize_model(model_spec)

  true_params <- c(
    A1.m = log(0.28),
    A1.s = 0.16,
    A2.m = log(0.28),
    A2.s = 0.16,
    A3.m = log(0.28),
    A3.s = 0.16,
    B.m = log(0.28),
    B.s = 0.18
  )

  set.seed(123456)
  params_df <- build_param_matrix(model_spec, true_params, n_trials = 1500)
  sim <- simulate(structure, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  neg_loglik <- function(theta) {
    est <- true_params
    est[c("A1.m", "A2.m", "A3.m")] <- theta[["A_pool.m"]]
    est[c("A1.s", "A2.s", "A3.s")] <- exp(theta[["log_A_pool.s"]])
    est["B.m"] <- theta[["B.m"]]
    est["B.s"] <- exp(theta[["log_B.s"]])
    params_df <- build_param_matrix(
      model_spec,
      est,
      trial_df = prepared
    )
    -as.numeric(log_likelihood(ctx, prepared, params_df))
  }

  start <- c(
    A_pool.m = log(0.24),
    log_A_pool.s = log(0.12),
    B.m = log(0.34),
    log_B.s = log(0.12)
  )

  set.seed(123456)
  fit <- optim(
    start,
    neg_loglik,
    method = "Nelder-Mead",
    control = list(maxit = 4000, reltol = 1e-9)
  )

  save(fit, file = "vignettes/data/pool_model.RData")
}

save_trigger_model <- function() {
  model_spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_outcome("R1", "go1") |>
    add_outcome("R2", "go2") |>
    add_trigger("shared_trigger",
      members = c("go1", "go2"),
      q = 0.15,
      draw = "shared"
    )

  structure <- finalize_model(model_spec)

  true_params <- c(
    go1.m = log(0.30),
    go1.s = 0.18,
    go2.m = log(0.35),
    go2.s = 0.18
  )

  set.seed(123456)
  params_df <- build_param_matrix(model_spec, true_params, n_trials = 1500)
  sim <- simulate(structure, params_df)
  data_df <- data.frame(
    trial = sim$trial,
    R = factor(sim$R),
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  neg_loglik <- function(theta) {
    est <- true_params
    est["go1.m"] <- theta[["go1.m"]]
    est["go1.s"] <- exp(theta[["log_go1.s"]])
    est["go2.m"] <- theta[["go2.m"]]
    est["go2.s"] <- exp(theta[["log_go2.s"]])
    params_df <- build_param_matrix(
      model_spec,
      est,
      trial_df = prepared
    )
    -as.numeric(log_likelihood(ctx, prepared, params_df))
  }

  start <- c(
    go1.m = log(0.25),
    log_go1.s = log(0.12),
    go2.m = log(0.42),
    log_go2.s = log(0.12)
  )

  set.seed(123456)
  fit <- optim(
    start,
    neg_loglik,
    method = "Nelder-Mead",
    control = list(maxit = 4000, reltol = 1e-9)
  )

  save(fit, file = "vignettes/data/trigger_model.RData")
}

save_mixtures <- function() {
  sampled_spec <- race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("TARGET", c("target_fast", "target_slow")) |>
    add_outcome("R1", "TARGET") |>
    add_outcome("R2", "competitor") |>
    add_component("fast", members = c("target_fast", "competitor"), weight_param = "p_fast") |>
    add_component("slow", members = c("target_slow", "competitor")) |>
    set_mixture_options(mode = "sample", reference = "slow")

  sampled_structure <- finalize_model(sampled_spec)

  true_params_sampled <- c(
    target_fast.m = log(0.25),
    target_fast.s = 0.15,
    target_slow.m = log(0.45),
    target_slow.s = 0.20,
    competitor.m = log(0.35),
    competitor.s = 0.18,
    p_fast = 0.35
  )

  set.seed(123456)
  params_df_sampled <- build_param_matrix(
    sampled_spec,
    true_params_sampled,
    n_trials = 1500
  )
  sim_sampled <- simulate(sampled_structure, params_df_sampled)
  data_sampled <- sim_sampled[, c("trial", "R", "rt")]

  prepared_sampled <- prepare_data(sampled_structure, data_sampled)
  ctx_sampled <- make_context(sampled_structure)

  neg_loglik <- function(theta) {
    est <- true_params_sampled
    est["p_fast"] <- plogis(theta[["logit_p_fast"]])
    params_df <- build_param_matrix(
      sampled_spec,
      est,
      trial_df = prepared_sampled
    )
    -as.numeric(log_likelihood(ctx_sampled, prepared_sampled, params_df))
  }

  start <- c(logit_p_fast = qlogis(0.60))

  set.seed(123456)
  fit <- optim(start, neg_loglik, method = "BFGS")

  save(fit, file = "vignettes/data/mixtures.RData")
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
    B.m = log(0.38),
    B.s = 0.22
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
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  neg_loglik <- function(theta) {
    est <- true_params
    est[c("A.m", "A.s", "B.m", "B.s")] <- theta[c("A.m", "A.s", "B.m", "B.s")]
    est[c("A.s", "B.s")] <- exp(est[c("A.s", "B.s")])
    params_df <- build_param_matrix(
      model_spec,
      est,
      trial_df = prepared
    )
    -as.numeric(log_likelihood(ctx, prepared, params_df))
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
    B.m = log(0.1),
    B.s = 0.1,
    C.m = log(0.15),
    C.s = 0.1
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
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  neg_loglik <- function(theta) {
    est <- true_params
    est[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")] <- theta[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")]
    est[c("A.s", "B.s", "C.s")] <- exp(est[c("A.s", "B.s", "C.s")])
    params_df <- build_param_matrix(
      model_spec,
      est,
      trial_df = prepared
    )
    -as.numeric(log_likelihood(ctx, prepared, params_df))
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
save_pool_model()
save_trigger_model()
save_mixtures()
save_multi_outcome()
save_chained_onset()
