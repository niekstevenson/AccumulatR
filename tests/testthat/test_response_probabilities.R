testthat::test_that("response_probabilities returns deterministic mass for a single outcome", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_outcome("A_win", "A")

  structure <- finalize_model(spec)
  params <- build_param_matrix(
    spec,
    c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
    n_trials = 1L
  )

  probs <- response_probabilities(structure, params)

  testthat::expect_equal(probs, c(A_win = 1), tolerance = 1e-4)
})

testthat::test_that("response_probabilities respects mixture weights and component filtering", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Left", "A") |>
    add_outcome("Right", "B") |>
    add_component("left_only", members = "A", weight = 0.25) |>
    add_component("right_only", members = "B", weight = 0.75) |>
    set_mixture_options(mode = "fixed")

  structure <- finalize_model(spec)
  params <- build_param_matrix(
    spec,
    c(
      A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0,
      B.m = 0, B.s = 0.1, B.q = 0, B.t0 = 0
    ),
    n_trials = 1L
  )

  probs <- response_probabilities(structure, params)
  testthat::expect_equal(
    probs[c("Left", "Right")],
    c(Left = 0.25, Right = 0.75),
    tolerance = 1e-4
  )

  rows <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows$component <- "left_only"
  probs_left <- response_probabilities(structure, rows)

  testthat::expect_equal(unname(probs_left["Left"]), 1, tolerance = 1e-4)
  testthat::expect_equal(unname(probs_left["Right"]), 0, tolerance = 1e-8)
})

testthat::test_that("response_probabilities returns residual NA mass for mapped outcomes", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Seen", "A") |>
    add_outcome("Miss", "B", options = list(map_outcome_to = NA_character_)) |>
    add_component("seen", members = "A", weight = 0.7) |>
    add_component("missing", members = "B", weight = 0.3) |>
    set_mixture_options(mode = "fixed")

  structure <- finalize_model(spec)
  params <- build_param_matrix(
    spec,
    c(
      A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0,
      B.m = 0, B.s = 0.1, B.q = 0, B.t0 = 0
    ),
    n_trials = 1L
  )

  probs <- response_probabilities(structure, params, include_na = TRUE)
  testthat::expect_equal(
    probs[c("Seen", "NA")],
    stats::setNames(c(0.7, 0.3), c("Seen", "NA")),
    tolerance = 1e-4
  )

  probs_no_na <- response_probabilities(structure, params, include_na = FALSE)
  testthat::expect_equal(probs_no_na, c(Seen = 0.7), tolerance = 1e-4)
})
