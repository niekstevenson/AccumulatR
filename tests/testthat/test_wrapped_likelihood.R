testthat::test_that("guess outcomes are aggregated in observed finite-label likelihood", {
  spec <- race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_accumulator("timeout", "lognormal", onset = 0.05) |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
    ))

  structure <- finalize_model(spec)
  params <- c(
    go_left.m = log(0.30), go_left.s = 0.18, go_left.q = 0.00, go_left.t0 = 0.00,
    go_right.m = log(0.325), go_right.s = 0.18, go_right.q = 0.00, go_right.t0 = 0.00,
    timeout.m = log(0.25), timeout.s = 0.10, timeout.q = 0.00, timeout.t0 = 0.00
  )
  rt_obs <- 0.35
  data_df <- data.frame(
    trial = 1L,
    R = "Right",
    rt = rt_obs,
    stringsAsFactors = FALSE
  )

  prepared <- prepare_data(structure, data_df)
  params_df <- build_param_matrix(spec, params, trial_df = prepared)
  ctx <- make_context(structure)
  ll <- as.numeric(log_likelihood(ctx, prepared, params_df))

  right_term <- stats::dlnorm(rt_obs, params[["go_right.m"]], params[["go_right.s"]]) *
    (1 - stats::plnorm(rt_obs, params[["go_left.m"]], params[["go_left.s"]])) *
    (1 - stats::plnorm(rt_obs - 0.05, params[["timeout.m"]], params[["timeout.s"]]))
  timeout_term <- stats::dlnorm(rt_obs - 0.05, params[["timeout.m"]], params[["timeout.s"]]) *
    (1 - stats::plnorm(rt_obs, params[["go_left.m"]], params[["go_left.s"]])) *
    (1 - stats::plnorm(rt_obs, params[["go_right.m"]], params[["go_right.s"]]))
  expected <- right_term + 0.8 * timeout_term

  testthat::expect_equal(ll, log(expected), tolerance = 2e-3)
})
