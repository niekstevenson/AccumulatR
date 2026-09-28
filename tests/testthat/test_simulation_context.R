testthat::test_that("simulation survives model copies and serialization before and after use", {
  model <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    finalize_model()
  params <- build_param_matrix(model, c(m = log(0.3), s = 0.2), n_trials = 20L)
  unused <- serialize(model, NULL)
  expected <- simulate(model, params, seed = 47, keep_detail = TRUE)
  used <- serialize(model, NULL)

  copy <- model
  attr(copy, "label") <- "copy"
  rm(model)
  invisible(gc())
  testthat::expect_identical(simulate(copy, params, seed = 47, keep_detail = TRUE), expected)
  rm(copy)
  invisible(gc())

  for (saved in list(unused, used)) {
    restored <- unserialize(saved)
    testthat::expect_identical(simulate(restored, params, seed = 47, keep_detail = TRUE), expected)
    invisible(gc())
    testthat::expect_identical(simulate(restored, params, seed = 47, keep_detail = TRUE), expected)
  }
})
