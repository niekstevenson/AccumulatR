testthat::test_that("set_parameters renames and shares parameters", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_outcome("go", "go") |>
    add_outcome("stop", "stop") |>
    set_parameters(list(
      drift = c("go.m", "stop.m"),
      spread = c("go.s", "stop.s"),
      onset = "go.t0"
    ))

  testthat::expect_equal(
    sampled_pars(spec),
    c("drift", "spread", "go.q", "onset", "stop.q", "stop.t0")
  )

  params_df <- build_param_matrix(
    spec,
    c(
      drift = log(0.30),
      spread = 0.18,
      go.q = 0.05,
      onset = 0.02,
      stop.q = 0.10,
      stop.t0 = 0.04
    ),
    n_trials = 1L
  )

  testthat::expect_equal(params_df[, "p1"], c(log(0.30), log(0.30)))
  testthat::expect_equal(params_df[, "p2"], c(0.18, 0.18))
  testthat::expect_equal(params_df[, "q"], c(0.05, 0.10))
  testthat::expect_equal(params_df[, "t0"], c(0.02, 0.04))
})

testthat::test_that("set_parameters rejects unknown and duplicated targets", {
  base_spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_outcome("go", "go") |>
    add_outcome("stop", "stop")

  unknown_spec <- set_parameters(base_spec, list(shared = c("go.m", "missing.m")))
  testthat::expect_error(
    sampled_pars(unknown_spec),
    "unknown parameter"
  )

  duplicate_spec <- set_parameters(base_spec, list(a = "go.m", b = "go.m"))
  testthat::expect_error(
    sampled_pars(duplicate_spec),
    "cannot be assigned more than once"
  )
})
