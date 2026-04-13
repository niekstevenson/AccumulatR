pkgload::load_all(".", quiet = TRUE)

testthat::test_that("test_model fills missing parameters and respects set_parameters", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_outcome("go", "go") |>
    add_outcome("stop", "stop") |>
    set_parameters(list(
      drift = c("go.m", "stop.m")
    ))

  testthat::expect_message(
    res <- test_model(
      spec,
      c(
        drift = log(0.30),
        go.s = 0.15,
        stop.s = 0.18
      ),
      n_trials = 25L,
      seed = 1,
      profile_points = 5L,
      plot = FALSE
    ),
    "Assuming 0 for missing parameter"
  )

  testthat::expect_true(all(c("go", "stop") %in% res$comparison$response))
  testthat::expect_true(all(c("go.q", "go.t0", "stop.q", "stop.t0") %in% names(res$parameters)))
  testthat::expect_equal(names(res$profiles), c("drift", "go.s", "stop.s"))
})
