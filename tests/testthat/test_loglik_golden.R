fixture_path <- testthat::test_path("fixtures", "loglik_golden_v1.rds")

testthat::test_that("log-likelihood goldens stay stable across representative models", {
  testthat::expect_true(
    file.exists(fixture_path),
    info = "Missing tests/testthat/fixtures/loglik_golden_v1.rds"
  )

  fixture <- readRDS(fixture_path)
  registry <- loglik_golden_case_registry()
  fixture_cases <- setNames(fixture$cases, vapply(fixture$cases, `[[`, character(1), "label"))

  testthat::expect_setequal(names(fixture_cases), names(registry))

  for (label in names(registry)) {
    expected <- fixture_cases[[label]]
    values <- evaluate_loglik_golden_case(
      registry[[label]]$build(),
      expected$data_df,
      expected$params
    )
    tol <- expected$tolerance %||% 1e-4

    testthat::expect_equal(
      values$total_loglik,
      expected$total_loglik,
      tolerance = tol,
      info = paste(label, "total_loglik")
    )
    testthat::expect_equal(
      values$trial_loglik,
      expected$trial_loglik,
      tolerance = tol,
      info = paste(label, "trial_loglik")
    )
  }
})
