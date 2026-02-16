testthat::test_that("race_spec defaults n_outcomes to 1", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  finalized <- finalize_model(spec)

  testthat::expect_equal(finalized$model_spec$metadata$observation$mode, "top_k")
  testthat::expect_equal(finalized$model_spec$metadata$observation$n_outcomes, 1L)
  testthat::expect_equal(finalized$prep$observation$mode, "top_k")
  testthat::expect_equal(finalized$prep$observation$n_outcomes, 1L)
})

testthat::test_that("race_spec stores custom n_outcomes", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  finalized <- finalize_model(spec)
  testthat::expect_equal(finalized$model_spec$metadata$observation$n_outcomes, 2L)
  testthat::expect_equal(finalized$prep$observation$n_outcomes, 2L)
})

testthat::test_that("invalid n_outcomes values fail fast", {
  testthat::expect_error(race_spec(n_outcomes = 0), "n_outcomes must be >= 1")
  testthat::expect_error(race_spec(n_outcomes = 1.5), "n_outcomes must be an integer value")
  testthat::expect_error(race_spec(n_outcomes = c(1, 2)), "n_outcomes must be a single non-missing value")
  testthat::expect_error(race_spec(n_outcomes = NA_real_), "n_outcomes must be a single non-missing value")
})

testthat::test_that("n_outcomes cannot exceed number of outcomes", {
  spec <- race_spec(n_outcomes = 3L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  testthat::expect_error(
    finalize_model(spec),
    "cannot exceed number of declared outcomes"
  )
})

testthat::test_that("n_outcomes > 1 rejects non-direct and remapped outcomes", {
  spec_expr <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", all_of("a", "b")) |>
    add_outcome("B", "b")
  testthat::expect_error(finalize_model(spec_expr), "direct event outcomes")

  spec_guess <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a", options = list(
      guess = list(labels = c("B"), weights = c(1))
    )) |>
    add_outcome("B", "b")
  testthat::expect_error(finalize_model(spec_guess), "guess option not supported")

  spec_alias <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b", options = list(alias_of = "A"))
  testthat::expect_error(finalize_model(spec_alias), "alias_of option not supported")

  spec_map <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a", options = list(map_outcome_to = NA_character_)) |>
    add_outcome("B", "b")
  testthat::expect_error(finalize_model(spec_map), "map_outcome_to option not supported")
})

testthat::test_that("n_outcomes > 1 rejects deterministic overlapping outcomes", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("p_a", members = "a", k = 1L) |>
    add_outcome("A", "a") |>
    add_outcome("A_dup", "p_a") |>
    add_outcome("B", "b")

  testthat::expect_error(
    finalize_model(spec),
    "deterministic overlapping outcomes"
  )
})

testthat::test_that("legacy outcome flexibility remains when n_outcomes == 1", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_outcome("A", inhibit("a", by = "stop"), options = list(
      guess = list(labels = c("B"), weights = c(1))
    )) |>
    add_outcome("B", "b", options = list(map_outcome_to = "B_OBS"))

  finalized <- finalize_model(spec)
  testthat::expect_s3_class(finalized, "model_structure")
  testthat::expect_equal(finalized$model_spec$metadata$observation$n_outcomes, 1L)
})

testthat::test_that("component n_outcomes override extends observation width", {
  spec <- race_spec(n_outcomes = 1L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_component("single", members = c("a", "b"), n_outcomes = 1L, weight = 0.5) |>
    add_component("double", members = c("a", "b"), n_outcomes = 2L, weight = 0.5)

  finalized <- finalize_model(spec)
  testthat::expect_equal(finalized$prep$observation$global_n_outcomes, 1L)
  testthat::expect_equal(finalized$prep$observation$component_n_outcomes$double, 2L)
  testthat::expect_equal(finalized$prep$observation$n_outcomes, 2L)
})

testthat::test_that("component n_outcomes values are validated", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  testthat::expect_error(
    add_component(spec, "bad", members = c("a", "b"), n_outcomes = 0),
    "n_outcomes must be >= 1"
  )
  testthat::expect_error(
    add_component(spec, "bad2", members = c("a", "b"), attrs = list(n_outcomes = 1.5)),
    "n_outcomes must be an integer value"
  )
})
