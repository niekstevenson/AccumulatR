testthat::test_that("response_probabilities respects fixed mixture weights", {
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

  probs <- response_probabilities(make_context(structure), params)
  testthat::expect_equal(
    probs[c("Left", "Right")],
    c(Left = 0.25, Right = 0.75),
    tolerance = 1e-4
  )
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

  context <- make_context(structure)
  probs <- response_probabilities(context, params, include_na = TRUE)
  testthat::expect_equal(
    probs[c("Seen", "NA")],
    stats::setNames(c(0.7, 0.3), c("Seen", "NA")),
    tolerance = 1e-4
  )

  probs_no_na <- response_probabilities(context, params, include_na = FALSE)
  testthat::expect_equal(probs_no_na, c(Seen = 0.7), tolerance = 1e-4)
})

testthat::test_that("response_probabilities marginalizes latent sampled mixtures", {
  structure <- latent_sampled_mixture_spec()
  context <- make_context(structure)
  values <- latent_sampled_mixture_params(0.2)
  params <- build_param_matrix(
    structure,
    values,
    n_trials = 1L
  )
  probabilities <- response_probabilities(context, params)
  fast_target <- pnorm(
    (values[["competitor.m"]] - values[["target_fast.m"]]) /
      sqrt(values[["competitor.s"]]^2 + values[["target_fast.s"]]^2)
  )
  slow_target <- pnorm(
    (values[["competitor.m"]] - values[["target_slow.m"]]) /
      sqrt(values[["competitor.s"]]^2 + values[["target_slow.s"]]^2)
  )
  target <- 0.2 * fast_target + 0.8 * slow_target

  testthat::expect_equal(
    probabilities[c("Target", "Competitor")],
    c(Target = target, Competitor = 1 - target),
    tolerance = 1e-4
  )
})

testthat::test_that("response_probabilities averages heterogeneous parameter blocks", {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    test_separate_all_parameters() |>
    finalize_model()
  params <- build_param_matrix(
    structure,
    c(A.m = 0, A.s = 0.2, A.t0 = 0, B.m = 0, B.s = 0.2, B.t0 = 0),
    n_trials = 2L
  )
  params[c(1L, 3L), "p1"] <- log(c(0.2, 0.6))
  params[c(2L, 4L), "p1"] <- log(c(0.5, 0.3))

  p_a <- pnorm(
    (params[c(2L, 4L), "p1"] - params[c(1L, 3L), "p1"]) /
      sqrt(params[c(1L, 3L), "p2"]^2 + params[c(2L, 4L), "p2"]^2)
  )
  expected <- c(A = mean(p_a), B = 1 - mean(p_a))

  testthat::expect_equal(
    response_probabilities(make_context(structure), params)[c("A", "B")],
    expected,
    tolerance = 1e-4
  )
})
