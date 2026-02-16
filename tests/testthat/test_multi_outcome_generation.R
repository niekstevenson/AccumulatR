testthat::test_that("simulate adds R2/rt2 when n_outcomes = 2", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.35), b.sdlog = 0.20, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 200L)
  sim <- simulate(structure, params_df, seed = 123)

  testthat::expect_true(all(c("trial", "R", "rt", "R2", "rt2") %in% names(sim)))
  testthat::expect_false(any(is.na(sim$rt2)))
  testthat::expect_true(all(sim$rt <= sim$rt2))
  testthat::expect_true(all(sim$R %in% c("A", "B")))
  testthat::expect_true(all(sim$R2 %in% c("A", "B")))
})

testthat::test_that("simulate pads missing follow-up outcomes with NA", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.15, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.35), b.sdlog = 0.15, b.q = 1, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 50L)
  sim <- simulate(structure, params_df, seed = 456)

  testthat::expect_true(all(sim$R == "A"))
  testthat::expect_true(all(is.finite(sim$rt)))
  testthat::expect_true(all(is.na(sim$R2)))
  testthat::expect_true(all(is.na(sim$rt2)))
})

testthat::test_that("simulate remains backward compatible when n_outcomes = 1", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.35), b.sdlog = 0.20, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 25L)
  sim <- simulate(structure, params_df, seed = 789)

  testthat::expect_true(all(c("trial", "R", "rt") %in% names(sim)))
  testthat::expect_false(any(c("R2", "rt2") %in% names(sim)))
})

testthat::test_that("simulate extends naming beyond second response", {
  spec <- race_spec(n_outcomes = 3L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.35), b.sdlog = 0.20, b.q = 0, b.t0 = 0,
    c.meanlog = log(0.40), c.sdlog = 0.20, c.q = 0, c.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 100L)
  sim <- simulate(structure, params_df, seed = 321)

  testthat::expect_true(all(c("R2", "rt2", "R3", "rt3") %in% names(sim)))
  testthat::expect_false(any(is.na(sim$rt3)))
  testthat::expect_true(all(sim$rt <= sim$rt2))
  testthat::expect_true(all(sim$rt2 <= sim$rt3))
})

testthat::test_that("component n_outcomes override controls structural rank NAs", {
  spec <- race_spec(n_outcomes = 1L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_component("single", members = c("a", "b"), n_outcomes = 1L, weight = 0.5) |>
    add_component("double", members = c("a", "b"), n_outcomes = 2L, weight = 0.5)

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0, a.t0 = 0,
    b.meanlog = log(0.35), b.sdlog = 0.20, b.q = 0, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 4L)
  trial_df <- data.frame(
    trial = 1:4,
    component = c("single", "double", "single", "double"),
    stringsAsFactors = FALSE
  )
  sim <- simulate(structure, params_df, trial_df = trial_df, seed = 42)

  testthat::expect_true(all(c("R2", "rt2") %in% names(sim)))
  single_rows <- sim$component == "single"
  double_rows <- sim$component == "double"
  testthat::expect_true(all(is.na(sim$R2[single_rows])))
  testthat::expect_true(all(is.na(sim$rt2[single_rows])))
  testthat::expect_true(all(!is.na(sim$R2[double_rows])))
  testthat::expect_true(all(is.finite(sim$rt2[double_rows])))
})
