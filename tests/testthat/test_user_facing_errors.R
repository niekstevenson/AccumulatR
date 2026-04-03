testthat::test_that("model-building functions report readable spec errors", {
  bad_spec <- list()

  testthat::expect_error(
    add_accumulator(bad_spec, "a", "lognormal"),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    add_pool(bad_spec, "p", members = "a"),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    add_outcome(bad_spec, "A", "a"),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    add_component(bad_spec, "fast", members = "a"),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    add_trigger(bad_spec, "tg", members = "a"),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    add_shared_params(bad_spec, members = "a", meanlog = 0),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    set_mixture_options(bad_spec),
    "expects a model specification created with race_spec"
  )
  testthat::expect_error(
    set_metadata(bad_spec, rel_tol = 1e-4),
    "expects a model specification created with race_spec"
  )
})
