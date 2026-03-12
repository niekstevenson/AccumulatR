race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
add_trigger <- AccumulatR::add_trigger
add_component <- AccumulatR::add_component
all_of <- AccumulatR::all_of
inhibit <- AccumulatR::inhibit
set_mixture_options <- AccumulatR::set_mixture_options
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood

testthat::test_that("build_param_matrix applies named shared-trigger q to all members", {
  spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_outcome("R1", "go1") |>
    add_outcome("R2", "go2") |>
    add_trigger("tg", members = c("go1", "go2"), q = 0.10, param = "tg_q", draw = "shared")

  params <- c(
    go1.meanlog = log(0.30),
    go1.sdlog = 0.18,
    go2.meanlog = log(0.32),
    go2.sdlog = 0.19,
    tg_q = 0.27
  )

  pm <- build_param_matrix(spec, params, n_trials = 2L)

  testthat::expect_equal(unname(pm[, "q"]), rep(0.27, 4L))
})

testthat::test_that("shared-trigger likelihood is stable across multi-trial onset changes", {
  spec <- race_spec() |>
    add_accumulator("S", "exgauss") |>
    add_accumulator("stop", "exgauss") |>
    add_accumulator("change", "exgauss") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_component("go_only", members = "S", weight = 0.75) |>
    add_component("go_stop", members = c("S", "stop", "change"), weight = 0.25) |>
    add_trigger("tg", members = c("stop", "change"), q = 0.10, param = "tg_q", draw = "shared") |>
    set_mixture_options(mode = "fixed")

  params <- c(
    S.mu = 0.38,
    S.sigma = 0.05,
    S.tau = 0.06,
    stop.mu = 0.08,
    stop.sigma = 0.16,
    stop.tau = 0.09,
    change.mu = 0.33,
    change.sigma = 0.09,
    change.tau = 0.05,
    tg_q = 0.22
  )

  batch_data <- data.frame(
    trial = c(1L, 1L, 1L, 2L, 2L, 2L),
    R = factor(c("X", "X", "X", "X", "X", "X"), levels = c("S", "X")),
    rt = c(0.3833333, 0.3833333, 0.3833333, 0.4666667, 0.4666667, 0.4666667),
    component = factor(
      c("go_stop", "go_stop", "go_stop", "go_stop", "go_stop", "go_stop"),
      levels = c("go_only", "go_stop")
    ),
    accumulator = c("S", "stop", "change", "S", "stop", "change"),
    onset = c(0.0, 0.1, 0.1, 0.0, 0.3, 0.3)
  )

  ctx_batch <- build_likelihood_context(spec, batch_data)
  pm_batch <- build_param_matrix(spec, params, n_trials = 2L)
  ll_batch <- as.numeric(log_likelihood(ctx_batch, batch_data, pm_batch))

  ll_split <- 0
  for (trial_id in 1:2) {
    trial_data <- subset(batch_data, trial == trial_id)
    trial_data$trial <- 1L
    ctx_trial <- build_likelihood_context(spec, trial_data)
    pm_trial <- build_param_matrix(spec, params, n_trials = 1L)
    ll_split <- ll_split + as.numeric(log_likelihood(ctx_trial, trial_data, pm_trial))
  }

  testthat::expect_equal(ll_batch, ll_split, tolerance = 1e-10)
})
