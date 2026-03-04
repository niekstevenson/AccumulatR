race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
inhibit <- AccumulatR::inhibit
add_pool <- AccumulatR::add_pool
add_trigger <- AccumulatR::add_trigger
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood

testthat::test_that("ranked single-step truncated path is deterministic", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("alt", "lognormal") |>
    add_outcome("GO", "go") |>
    add_outcome("ALT", "alt")

  structure <- finalize_model(spec)
  params <- c(
    go.meanlog = log(0.34), go.sdlog = 0.20,
    alt.meanlog = log(0.47), alt.sdlog = 0.19
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  ranked_single_df <- data.frame(
    trial = 1L,
    R = factor("GO", levels = c("GO", "ALT")),
    rt = 0.46,
    R2 = factor(NA_character_, levels = c("GO", "ALT")),
    rt2 = NA_real_
  )

  ctx <- build_likelihood_context(structure, ranked_single_df)
  ll1 <- as.numeric(log_likelihood(ctx, ranked_single_df, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, ranked_single_df, params_df))

  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})

testthat::test_that("ranked stopper two-outcome path is finite and deterministic", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("alt", "lognormal") |>
    add_outcome("GO", inhibit("go", by = "stop")) |>
    add_outcome("ALT", "alt")

  structure <- finalize_model(spec)
  params <- c(
    go.meanlog = log(0.34), go.sdlog = 0.20,
    stop.meanlog = log(0.41), stop.sdlog = 0.18,
    alt.meanlog = log(0.47), alt.sdlog = 0.19
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  ranked_single_df <- data.frame(
    trial = 1L,
    R = factor("GO", levels = c("GO", "ALT")),
    rt = 0.46,
    R2 = factor(NA_character_, levels = c("GO", "ALT")),
    rt2 = NA_real_
  )

  ctx <- build_likelihood_context(structure, ranked_single_df)
  ll1 <- as.numeric(log_likelihood(ctx, ranked_single_df, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, ranked_single_df, params_df))

  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})

testthat::test_that("ranked shared-trigger multi-rank path is deterministic", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("A", c("a1", "a2"), k = 1L) |>
    add_outcome("A", "A") |>
    add_outcome("B", "b") |>
    add_trigger("tg", members = c("a1", "a2"), q = 0.10, draw = "shared")

  structure <- finalize_model(spec)
  params <- c(
    a1.meanlog = log(0.30), a1.sdlog = 0.20, a1.q = 0.10, a1.t0 = 0,
    a2.meanlog = log(0.33), a2.sdlog = 0.20, a2.q = 0.10, a2.t0 = 0,
    b.meanlog = log(0.52), b.sdlog = 0.18, b.q = 0.00, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  ranked_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("A", "B")),
    rt = 0.36,
    R2 = factor("B", levels = c("A", "B")),
    rt2 = 0.61
  )

  ctx <- build_likelihood_context(structure, ranked_df)
  ll1 <- as.numeric(log_likelihood(ctx, ranked_df, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, ranked_df, params_df))

  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})

testthat::test_that("invalid ranked ordering returns min_ll", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.34), a.sdlog = 0.20,
    b.meanlog = log(0.48), b.sdlog = 0.18
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  invalid_ranked_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("A", "B")),
    rt = 0.50,
    R2 = factor("B", levels = c("A", "B")),
    rt2 = 0.40
  )

  min_ll <- log(1e-9)
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, invalid_ranked_df),
    invalid_ranked_df,
    params_df,
    min_ll = min_ll
  ))

  testthat::expect_equal(ll, min_ll, tolerance = 1e-12)
})

testthat::test_that("nonresponse is stable with ranked columns present", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.34), a.sdlog = 0.20,
    b.meanlog = log(0.48), b.sdlog = 0.18
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  plain_nonresp_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("A", "B")),
    rt = NA_real_
  )
  ranked_nonresp_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("A", "B")),
    rt = NA_real_,
    R2 = factor(NA_character_, levels = c("A", "B")),
    rt2 = NA_real_
  )

  ll_plain <- as.numeric(log_likelihood(
    build_likelihood_context(structure, plain_nonresp_df),
    plain_nonresp_df,
    params_df
  ))
  ll_ranked <- as.numeric(log_likelihood(
    build_likelihood_context(structure, ranked_nonresp_df),
    ranked_nonresp_df,
    params_df
  ))

  testthat::expect_true(is.finite(ll_plain))
  testthat::expect_true(is.finite(ll_ranked))
  testthat::expect_equal(ll_ranked, ll_plain, tolerance = 1e-10)
})
