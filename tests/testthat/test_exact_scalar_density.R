race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_component <- AccumulatR::add_component
add_outcome <- AccumulatR::add_outcome
add_trigger <- AccumulatR::add_trigger
all_of <- AccumulatR::all_of
build_likelihood_context <- AccumulatR::build_likelihood_context
build_param_matrix <- AccumulatR::build_param_matrix
finalize_model <- AccumulatR::finalize_model
inhibit <- AccumulatR::inhibit
log_likelihood <- AccumulatR::log_likelihood
set_mixture_options <- AccumulatR::set_mixture_options

testthat::test_that("BEEST single-step likelihood matches the exact oracle", {
  source(testthat::test_path("../../dev/test_beest.R"), local = TRUE)

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

  pkg_params <- c(
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
  oracle_params <- c(
    S.mu = 0.38,
    S.sigma = 0.05,
    S.tau = 0.06,
    stop.mu = 0.08,
    stop.sigma = 0.16,
    stop.tau = 0.09,
    change.mu = 0.33,
    change.sigma = 0.09,
    change.tau = 0.05,
    stop_trigger = 0.22
  )

  data_df <- data.frame(
    trial = c(1L, 2L, 2L, 2L, 3L, 3L, 3L),
    R = factor(c("S", "S", "S", "S", "X", "X", "X"), levels = c("S", "X")),
    rt = c(0.43, 0.47, 0.47, 0.47, 0.39, 0.39, 0.39),
    component = factor(
      c("go_only", "go_stop", "go_stop", "go_stop", "go_stop", "go_stop", "go_stop"),
      levels = c("go_only", "go_stop")
    ),
    accumulator = c("S", "S", "stop", "change", "S", "stop", "change"),
    onset = c(0.00, 0.00, 0.12, 0.12, 0.00, 0.10, 0.10)
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, pkg_params, n_trials = 3L)

  ll_cpp <- as.numeric(log_likelihood(ctx, data_df, params_df))

  ll_oracle <- make_beest_loglik(data_df, include_component_weight = FALSE)(oracle_params)

  testthat::expect_true(is.finite(ll_cpp))
  testthat::expect_equal(ll_cpp, ll_oracle, tolerance = 1e-3)
})

testthat::test_that("generic coupling integration returns a finite probability for shared-source composites", {
  spec <- race_spec() |>
    add_accumulator("S", "exgauss") |>
    add_accumulator("stop", "exgauss") |>
    add_accumulator("change", "exgauss") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop"))

  structure <- finalize_model(spec)
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  ctx <- AccumulatR:::native_context_build(prep)
  x_compiled <- AccumulatR:::.expr_lookup_compiled(prep$outcomes[["X"]]$expr, prep)
  s_compiled <- AccumulatR:::.expr_lookup_compiled(prep$outcomes[["S"]]$expr, prep)

  prob <- AccumulatR:::native_outcome_probability_params_cpp_idx(
    ctx,
    as.integer(x_compiled$id),
    0.60,
    -1L,
    integer(),
    integer(),
    as.integer(s_compiled$id),
    1e-4,
    1e-6,
    8L,
    NULL
  )
  testthat::expect_true(is.finite(prob))
})

testthat::test_that("plain independent races use the direct competing fast path", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B")

  data_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("A", "B")),
    rt = 0.31
  )
  params <- c(
    A.meanlog = log(0.28), A.sdlog = 0.15, A.q = 0.00,
    B.meanlog = log(0.34), B.sdlog = 0.17, B.q = 0.00
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll <- as.numeric(log_likelihood(ctx, data_df, params_df))
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  testthat::expect_true(is.finite(ll))
  testthat::expect_gt(
    stats$direct_node_density_batch_simple_competing_fastpath_calls,
    0
  )
})

testthat::test_that("timeout guess races use the simple competing fast path", {
  spec <- race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_accumulator("timeout", "lognormal", onset = 0.05) |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
    ))

  data_df <- data.frame(
    trial = 1L,
    R = factor("Left", levels = c("Left", "Right", "TIMEOUT")),
    rt = 0.31
  )
  params <- c(
    go_left.meanlog = log(0.30),
    go_left.sdlog = 0.18,
    go_right.meanlog = log(0.325),
    go_right.sdlog = 0.18,
    timeout.meanlog = log(0.25),
    timeout.sdlog = 0.10
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll <- as.numeric(log_likelihood(ctx, data_df, params_df))
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  testthat::expect_true(is.finite(ll))
  testthat::expect_gt(
    stats$direct_node_density_batch_simple_competing_fastpath_calls,
    0
  )
  testthat::expect_equal(stats$direct_node_density_batch_kernel_eval_calls, 0)
  testthat::expect_equal(stats$kernel_event_batch_calls, 0)
})

testthat::test_that("disjoint guard competitors evaluate finitely", {
  spec <- race_spec() |>
    add_accumulator("plain", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("PLAIN", "plain") |>
    add_outcome("GUARD", inhibit("a", by = "b"))

  data_df <- data.frame(
    trial = 1L,
    R = factor("PLAIN", levels = c("PLAIN", "GUARD")),
    rt = 0.41
  )
  params <- c(
    plain.meanlog = log(0.38), plain.sdlog = 0.16,
    a.meanlog = log(0.31), a.sdlog = 0.15,
    b.meanlog = log(0.27), b.sdlog = 0.15
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  ll <- as.numeric(log_likelihood(ctx, data_df, params_df))

  testthat::expect_true(is.finite(ll))
})

testthat::test_that("component filtering yields finite per-component likelihoods", {
  spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("stop", "exgauss", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_outcome("R1", inhibit("go1", by = "stop")) |>
    add_outcome("R2", all_of("go2", "stop")) |>
    add_component("go_only", members = "go1", weight = 0.5) |>
    add_component("go_stop", members = c("go1", "stop", "go2"), weight = 0.5) |>
    set_mixture_options(mode = "fixed")

  params <- c(
    go1.meanlog = log(0.35), go1.sdlog = 0.20,
    stop.mu = 0.1, stop.sigma = 0.04, stop.tau = 0.1,
    go2.meanlog = log(0.60), go2.sdlog = 0.18
  )

  data_go_only <- data.frame(
    trial = 1L,
    R = factor("R1", levels = c("R1", "R2")),
    rt = 0.42,
    component = factor("go_only", levels = c("go_only", "go_stop"))
  )
  ctx_go_only <- build_likelihood_context(spec, data_go_only)
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  ll_go_only <- as.numeric(log_likelihood(ctx_go_only, data_go_only, params_df))

  data_go_stop <- data.frame(
    trial = 1L,
    R = factor("R1", levels = c("R1", "R2")),
    rt = 0.42,
    component = factor("go_stop", levels = c("go_only", "go_stop"))
  )
  ctx_go_stop <- build_likelihood_context(spec, data_go_stop)

  ll_go_stop <- as.numeric(log_likelihood(ctx_go_stop, data_go_stop, params_df))

  testthat::expect_true(is.finite(ll_go_only))
  testthat::expect_true(is.finite(ll_go_stop))
})

testthat::test_that("exact batch log-likelihood matches the one-trial exact kernel", {
  spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("stop", "exgauss", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_outcome("R1", inhibit("go1", by = "stop"),
      options = list(component = c("go_only", "go_stop"))
    ) |>
    add_outcome("R2", all_of("go2", "stop"), options = list(component = "go_stop")) |>
    add_component("go_only", members = c("go1"), weight = 0.5) |>
    add_component("go_stop", members = c("go1", "stop", "go2"), weight = 0.5) |>
    set_mixture_options(mode = "fixed")

  data_df <- data.frame(
    trial = 1:3,
    R = factor(rep("R1", 3), levels = c("R1", "R2")),
    rt = c(0.36, 0.42, 0.51),
    component = factor(rep("go_stop", 3), levels = c("go_only", "go_stop"))
  )
  params <- c(
    go1.meanlog = log(0.35), go1.sdlog = 0.20,
    stop.mu = 0.1, stop.sigma = 0.04, stop.tau = 0.1,
    go2.meanlog = log(0.60), go2.sdlog = 0.18
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, params, n_trials = nrow(data_df))

  ll_batch <- as.numeric(log_likelihood(ctx, data_df, params_df))

  ll_scalar_sum <- 0.0
  for (trial_id in data_df$trial) {
    ok <- data_df$trial == trial_id
    ll_scalar_sum <- ll_scalar_sum +
      as.numeric(log_likelihood(
        ctx,
        data_df,
        params_df,
        ok = ok,
        expand = trial_id
      ))
  }

  testthat::expect_equal(ll_batch, ll_scalar_sum, tolerance = 1e-12)
})

testthat::test_that("guard-heavy exact batches match one-trial scalar evaluation", {
  spec <- race_spec() |>
    add_accumulator("go_fast", "lognormal") |>
    add_accumulator("go_slow", "lognormal") |>
    add_accumulator("gate_shared", "lognormal") |>
    add_accumulator("stop_control", "lognormal") |>
    add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
    add_outcome("Slow", all_of("go_slow", "gate_shared"))

  data_df <- data.frame(
    trial = 1:4,
    R = factor(c("Fast", "Fast", "Slow", "Fast"), levels = c("Fast", "Slow")),
    rt = c(0.31, 0.38, 0.41, 0.47)
  )
  params <- c(
    go_fast.meanlog = log(0.28), go_fast.sdlog = 0.18,
    go_slow.meanlog = log(0.34), go_slow.sdlog = 0.18,
    gate_shared.meanlog = log(0.30), gate_shared.sdlog = 0.16,
    stop_control.meanlog = log(0.27), stop_control.sdlog = 0.15
  )

  ctx <- build_likelihood_context(spec, data_df)
  params_df <- build_param_matrix(spec, params, n_trials = nrow(data_df))

  ll_batch <- as.numeric(log_likelihood(ctx, data_df, params_df))

  ll_scalar_sum <- 0.0
  for (trial_id in data_df$trial) {
    ok <- data_df$trial == trial_id
    ll_scalar_sum <- ll_scalar_sum +
      as.numeric(log_likelihood(
        ctx,
        data_df,
        params_df,
        ok = ok,
        expand = trial_id
      ))
  }

  testthat::expect_equal(ll_batch, ll_scalar_sum, tolerance = 1e-12)
})
