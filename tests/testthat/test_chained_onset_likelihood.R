testthat::test_that("single-outcome chained onset matches numeric convolution", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.m = log(0.25), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  rt_obs <- 0.85
  data_df <- data.frame(trial = 1L, R = "B", rt = rt_obs, stringsAsFactors = FALSE)
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  ll <- as.numeric(log_likelihood(ctx, prepared, params_df))
  expected <- stats::integrate(
    function(u) {
      stats::dlnorm(u, params[["a.m"]], params[["a.s"]]) *
        stats::dlnorm(rt_obs - u, params[["b.m"]], params[["b.s"]])
    },
    lower = 0, upper = rt_obs, rel.tol = 1e-8
  )$value
  testthat::expect_equal(ll, log(expected), tolerance = 2e-3)
})

testthat::test_that("ranked chained onset uses observed source time (A then B)", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.m = log(0.30), a.s = 0.12, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.25), b.s = 0.12, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  t1 <- 0.30
  t2 <- 0.62
  data_df <- data.frame(
    trial = 1L,
    R = "A", rt = t1,
    R2 = "B", rt2 = t2,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)

  ll <- as.numeric(log_likelihood(ctx, prepared, params_df))
  expected <- stats::dlnorm(t1, params[["a.m"]], params[["a.s"]]) *
    stats::dlnorm(t2 - t1, params[["b.m"]], params[["b.s"]])
  testthat::expect_equal(ll, log(expected), tolerance = 2e-3)
})

testthat::test_that("pool-source chained onset likelihood is finite", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_accumulator("c", "lognormal", onset = after("p", lag = 0.03)) |>
    add_outcome("C", "c")

  structure <- finalize_model(spec)
  params <- c(
    a.m = log(0.20), a.s = 0.10, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.24), b.s = 0.12, b.q = 0.00, b.t0 = 0.00,
    c.m = log(0.18), c.s = 0.11, c.q = 0.00, c.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(trial = 1L, R = "C", rt = 0.65, stringsAsFactors = FALSE)
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  ll <- as.numeric(log_likelihood(ctx, prepared, params_df))
  testthat::expect_true(is.finite(ll))
})

testthat::test_that("missing trigger implies non-start in chained likelihood", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.m = log(0.25), a.s = 0.10, a.q = 1.00, a.t0 = 0.00,
    b.m = log(0.25), b.s = 0.10, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  min_ll <- -1e6

  observed_df <- data.frame(trial = 1L, R = "B", rt = 0.8, stringsAsFactors = FALSE)
  observed_prepared <- prepare_data(structure, observed_df)
  ctx <- make_context(structure)
  ll_obs <- as.numeric(log_likelihood(
    ctx,
    observed_prepared,
    params_df,
    min_ll = min_ll
  ))
  testthat::expect_equal(ll_obs, min_ll)

  nonresponse_df <- data.frame(
    trial = 1L,
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  nonresponse_prepared <- prepare_data(structure, nonresponse_df)
  ll_nr <- as.numeric(log_likelihood(
    ctx,
    nonresponse_prepared,
    params_df,
    min_ll = min_ll
  ))
  testthat::expect_true(is.finite(ll_nr))
  testthat::expect_gt(ll_nr, min_ll)
})

testthat::test_that("component with inactive onset source yields non-start", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("B", "b") |>
    add_component("with_source", members = c("a", "b"), weight = 0.5) |>
    add_component("without_source", members = c("b"), weight = 0.5)

  structure <- finalize_model(spec)
  params <- c(
    a.m = log(0.22), a.s = 0.10, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.20), b.s = 0.10, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  min_ll <- -1e6
  observed_df <- data.frame(
    trial = 1L,
    component = "without_source",
    R = "B",
    rt = 0.70,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, observed_df)
  ctx <- make_context(structure)
  params_df <- build_param_matrix(spec, params, trial_df = prepared)
  ll <- as.numeric(log_likelihood(
    ctx,
    prepared,
    params_df,
    min_ll = min_ll
  ))
  testthat::expect_equal(ll, min_ll)
})
