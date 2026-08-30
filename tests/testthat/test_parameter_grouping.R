testthat::test_that("compatible accumulator parameters are shared by default", {
  spec <- race_spec() |>
    add_accumulator("left", "lognormal") |>
    add_accumulator("right", "lognormal") |>
    add_outcome("left", "left") |>
    add_outcome("right", "right")

  testthat::expect_equal(par_names(spec), c("m", "s", "t0"))

  parameters <- build_param_matrix(
    finalize_model(spec),
    c(m = log(0.30), s = 0.17, t0 = 0.02),
    n_trials = 2L
  )
  testthat::expect_equal(dim(parameters), c(4L, 4L))
  testthat::expect_equal(unique(parameters[, "p1"]), log(0.30))
  testthat::expect_equal(unique(parameters[, "p2"]), 0.17)
  testthat::expect_equal(unique(parameters[, "t0"]), 0.02)
})
