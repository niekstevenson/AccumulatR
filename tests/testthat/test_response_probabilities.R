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

latent_sampled_response_spec <- function() {
  race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("Target", c("target_fast", "target_slow")) |>
    add_outcome("Target", "Target") |>
    add_outcome("Competitor", "competitor") |>
    add_component(
      "fast",
      members = c("target_fast", "competitor"),
      weight_param = "p_fast"
    ) |>
    add_component("slow", members = c("target_slow", "competitor")) |>
    set_mixture_options(mode = "sample", reference = "slow") |>
    finalize_model()
}

latent_sampled_response_params <- function(p_fast) {
  c(
    target_fast.m = log(0.25),
    target_fast.s = 0.15,
    target_slow.m = log(0.45),
    target_slow.s = 0.20,
    competitor.m = log(0.35),
    competitor.s = 0.18,
    p_fast = p_fast
  )
}

testthat::test_that("response_probabilities marginalizes sampled mixtures when component is latent", {
  structure <- latent_sampled_response_spec()

  probs_lo <- response_probabilities(
    structure,
    build_param_matrix(
      structure$model_spec,
      latent_sampled_response_params(0.2),
      n_trials = 1L
    )
  )
  probs_hi <- response_probabilities(
    structure,
    build_param_matrix(
      structure$model_spec,
      latent_sampled_response_params(0.8),
      n_trials = 1L
    )
  )

  testthat::expect_gt(
    abs(unname(probs_hi["Target"]) - unname(probs_lo["Target"])),
    1e-6
  )
})

testthat::test_that("response_probabilities treats explicit NA component like latent and explicit labels like observed", {
  structure <- latent_sampled_response_spec()
  params <- build_param_matrix(
    structure$model_spec,
    latent_sampled_response_params(0.2),
    n_trials = 1L
  )

  probs_latent <- response_probabilities(structure, params)

  rows_na <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows_na$component <- NA_character_
  probs_explicit_na <- response_probabilities(structure, rows_na)

  rows_fast <- AccumulatR:::.param_matrix_to_rows(structure, params)
  rows_fast$component <- "fast"
  probs_fast <- response_probabilities(structure, rows_fast)

  testthat::expect_equal(probs_explicit_na, probs_latent, tolerance = 1e-8)
  testthat::expect_gt(
    abs(unname(probs_fast["Target"]) - unname(probs_latent["Target"])),
    1e-6
  )
})
