testthat::test_that("ranked likelihood matches simple independent two-outcome formula", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.45), b.sdlog = 0.25, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.30,
    R2 = "B", rt2 = 0.55,
    stringsAsFactors = FALSE
  )

  ctx <- build_likelihood_context(structure, data_df)
  ll <- as.numeric(log_likelihood(ctx, params_df))

  expected <- dlnorm(0.30, meanlog = log(0.30), sdlog = 0.20) *
    dlnorm(0.55, meanlog = log(0.45), sdlog = 0.25)
  testthat::expect_equal(ll, log(expected), tolerance = 1e-4)
})

testthat::test_that("ranked likelihood treats missing later ranks as truncation", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.45), b.sdlog = 0.25, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  ranked_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.30,
    R2 = NA_character_, rt2 = NA_real_,
    stringsAsFactors = FALSE
  )
  single_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.30,
    stringsAsFactors = FALSE
  )

  ll_ranked <- as.numeric(log_likelihood(build_likelihood_context(structure, ranked_df), params_df))
  ll_single <- as.numeric(log_likelihood(build_likelihood_context(structure, single_df), params_df))
  testthat::expect_equal(ll_ranked, ll_single, tolerance = 1e-8)
})

testthat::test_that("ranked likelihood rejects mismatched rank pairs with min_ll", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.45), b.sdlog = 0.25, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  min_ll <- -1e6

  bad_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.30,
    R2 = "B", rt2 = NA_real_,
    stringsAsFactors = FALSE
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, bad_df),
    params_df,
    min_ll = min_ll
  ))
  testthat::expect_equal(ll, min_ll)
})

testthat::test_that("likelihood context enforces contiguous ranked columns", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  bad_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.30,
    R3 = "B", rt3 = 0.55,
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    build_likelihood_context(structure, bad_df),
    "contiguous"
  )
})

testthat::test_that("ranked likelihood rejects non-increasing times with min_ll", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.45), b.sdlog = 0.25, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  min_ll <- -1e6

  bad_time_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.50,
    R2 = "B", rt2 = 0.40,
    stringsAsFactors = FALSE
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, bad_time_df),
    params_df,
    min_ll = min_ll
  ))
  testthat::expect_equal(ll, min_ll)
})

testthat::test_that("ranked likelihood handles pool + shared trigger models", {
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
    a2.meanlog = log(0.32), a2.sdlog = 0.20, a2.q = 0.10, a2.t0 = 0,
    b.meanlog = log(0.50), b.sdlog = 0.20, b.q = 0.00, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.35,
    R2 = "B", rt2 = 0.60,
    stringsAsFactors = FALSE
  )

  ll <- as.numeric(log_likelihood(build_likelihood_context(structure, data_df), params_df))
  testthat::expect_true(is.finite(ll))
})
