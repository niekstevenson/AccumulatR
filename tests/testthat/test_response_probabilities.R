testthat::test_that("response_probabilities respects mixture weights and component filtering", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Left", "A") |>
    add_outcome("Right", "B") |>
    add_component("left_only", members = "A") |>
    add_component("right_only", members = "B") |>
    set_mixture(mode = "fixed", weights = c(left_only = 0.25, right_only = 0.75)) |>
    test_separate_all_parameters()

  structure <- finalize_model(spec)
  params <- build_param_matrix(
    structure,
    c(
      A.m = 0, A.s = 0.1, A.t0 = 0,
      B.m = 0, B.s = 0.1, B.t0 = 0
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
    add_component("seen", members = "A") |>
    add_component("missing", members = "B") |>
    set_mixture(mode = "fixed", weights = c(seen = 0.7, missing = 0.3)) |>
    test_separate_all_parameters()

  structure <- finalize_model(spec)
  params <- build_param_matrix(
    structure,
    c(
      A.m = 0, A.s = 0.1, A.t0 = 0,
      B.m = 0, B.s = 0.1, B.t0 = 0
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

testthat::test_that("response_probabilities marginalizes latent sampled mixtures", {
  structure <- latent_sampled_mixture_spec()
  params <- build_param_matrix(
    structure,
    latent_sampled_mixture_params(0.2),
    n_trials = 1L
  )

  probs_latent <- response_probabilities(structure, params)

  rows_na <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows_na$component <- NA_character_
  probs_explicit_na <- response_probabilities(structure, rows_na)

  rows_fast <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows_fast$component <- "fast"
  probs_fast <- response_probabilities(structure, rows_fast)
  rows_slow <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows_slow$component <- "slow"
  probs_slow <- response_probabilities(structure, rows_slow)

  testthat::expect_equal(probs_explicit_na, probs_latent, tolerance = 1e-8)
  testthat::expect_equal(
    probs_latent,
    0.2 * probs_fast + 0.8 * probs_slow,
    tolerance = 1e-8
  )
})
