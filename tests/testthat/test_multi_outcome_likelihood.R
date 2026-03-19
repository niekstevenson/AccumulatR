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
    R = factor("A", levels = c("A", "B")), rt = 0.30,
    R2 = factor("B", levels = c("A", "B")), rt2 = 0.55
  )

  ctx <- build_likelihood_context(structure, data_df)
  ll <- as.numeric(log_likelihood(ctx, data_df, params_df))

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
    R = factor("A", levels = c("A", "B")), rt = 0.30,
    R2 = factor(NA_character_, levels = c("A", "B")), rt2 = NA_real_
  )
  single_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("A", "B")), rt = 0.30
  )

  ll_ranked <- as.numeric(log_likelihood(build_likelihood_context(structure, ranked_df), ranked_df, params_df))
  ll_single <- as.numeric(log_likelihood(build_likelihood_context(structure, single_df), single_df, params_df))
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
    R = factor("A", levels = c("A", "B")), rt = 0.30,
    R2 = factor("B", levels = c("A", "B")), rt2 = NA_real_
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, bad_df),
    bad_df,
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
    R = factor("A", levels = c("A", "B")), rt = 0.30,
    R3 = factor("B", levels = c("A", "B")), rt3 = 0.55
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
    R = factor("A", levels = c("A", "B")), rt = 0.50,
    R2 = factor("B", levels = c("A", "B")), rt2 = 0.40
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, bad_time_df),
    bad_time_df,
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
    R = factor("A", levels = c("A", "B")), rt = 0.35,
    R2 = factor("B", levels = c("A", "B")), rt2 = 0.60
  )

  ll <- as.numeric(log_likelihood(build_likelihood_context(structure, data_df), data_df, params_df))
  testthat::expect_true(is.finite(ll))
})

testthat::test_that("ranked chained-onset trials match per-trial evaluation", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.16, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.22), b.sdlog = 0.16, b.q = 0, b.t0 = 0
  )
  data_df <- data.frame(
    trial = 1:4,
    R = factor(rep("A", 4), levels = c("A", "B")),
    rt = c(0.29, 0.31, 0.34, 0.37),
    R2 = factor(rep("B", 4), levels = c("A", "B")),
    rt2 = c(0.51, 0.56, 0.60, 0.65)
  )
  params_df <- build_param_matrix(spec, params, n_trials = nrow(data_df))

  ctx <- build_likelihood_context(structure, data_df)
  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll_batched <- as.numeric(log_likelihood(ctx, data_df, params_df))
  stats_batched <- AccumulatR:::unified_outcome_stats_cpp()

  ll_split <- 0.0
  for (i in seq_len(nrow(data_df))) {
    row_df <- data_df[i, , drop = FALSE]
    row_params <- build_param_matrix(spec, params, n_trials = 1L)
    ll_split <- ll_split + as.numeric(log_likelihood(
      build_likelihood_context(structure, row_df),
      row_df,
      row_params
    ))
  }

  testthat::expect_gt(stats_batched$exact_node_batch_calls, 0)
  testthat::expect_equal(ll_batched, ll_split, tolerance = 1e-10)
})

testthat::test_that("ranked bounded-source chained-onset batches avoid scalar fallback", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_accumulator("c", "lognormal") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.16, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.22), b.sdlog = 0.16, b.q = 0, b.t0 = 0,
    c.meanlog = log(0.35), c.sdlog = 0.18, c.q = 0, c.t0 = 0
  )
  data_df <- data.frame(
    trial = 1:4,
    R = factor(rep("B", 4), levels = c("B", "C")),
    rt = c(0.42, 0.48, 0.55, 0.63),
    R2 = factor(rep("C", 4), levels = c("B", "C")),
    rt2 = c(0.70, 0.76, 0.83, 0.92)
  )
  params_df <- build_param_matrix(spec, params, n_trials = nrow(data_df))

  ctx <- build_likelihood_context(structure, data_df)
  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll_batched <- as.numeric(log_likelihood(ctx, data_df, params_df))
  stats_batched <- AccumulatR:::unified_outcome_stats_cpp()

  ll_split <- 0.0
  for (i in seq_len(nrow(data_df))) {
    row_df <- data_df[i, , drop = FALSE]
    row_params <- build_param_matrix(spec, params, n_trials = 1L)
    ll_split <- ll_split + as.numeric(log_likelihood(
      build_likelihood_context(structure, row_df),
      row_df,
      row_params
    ))
  }

  testthat::expect_gt(stats_batched$exact_node_batch_calls, 0)
  testthat::expect_equal(ll_batched, ll_split, tolerance = 1e-10)
})

testthat::test_that("ranked interval-conditioned self constraints batch without scalar fallback", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_accumulator("c", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.16, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.22), b.sdlog = 0.16, b.q = 0, b.t0 = 0,
    c.meanlog = log(0.35), c.sdlog = 0.18, c.q = 0, c.t0 = 0
  )

  check_case <- function(data_df) {
    params_df <- build_param_matrix(spec, params, n_trials = nrow(data_df))
    ctx <- build_likelihood_context(structure, data_df)
    AccumulatR:::unified_outcome_stats_reset_cpp()
    ll_batched <- as.numeric(log_likelihood(ctx, data_df, params_df))
    stats_batched <- AccumulatR:::unified_outcome_stats_cpp()

    ll_split <- 0.0
    for (i in seq_len(nrow(data_df))) {
      row_df <- data_df[i, , drop = FALSE]
      row_params <- build_param_matrix(spec, params, n_trials = 1L)
      ll_split <- ll_split + as.numeric(log_likelihood(
        build_likelihood_context(structure, row_df),
        row_df,
        row_params
      ))
    }

    testthat::expect_gt(stats_batched$exact_node_batch_calls, 0)
    testthat::expect_equal(ll_batched, ll_split, tolerance = 1e-10)
  }

  check_case(data.frame(
    trial = 1:4,
    R = factor(rep("B", 4), levels = c("A", "B", "C")),
    rt = c(0.42, 0.48, 0.55, 0.63),
    R2 = factor(rep("C", 4), levels = c("A", "B", "C")),
    rt2 = c(0.70, 0.76, 0.83, 0.92)
  ))

  check_case(data.frame(
    trial = 1:4,
    R = factor(rep("C", 4), levels = c("A", "B", "C")),
    rt = c(0.42, 0.48, 0.55, 0.63),
    R2 = factor(rep("A", 4), levels = c("A", "B", "C")),
    rt2 = c(0.70, 0.76, 0.83, 0.92)
  ))
})
