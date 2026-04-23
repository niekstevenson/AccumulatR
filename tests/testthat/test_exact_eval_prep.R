run_public_loglik <- function(spec, trial_df, params, min_ll = log(1e-10)) {
  structure <- finalize_model(spec)
  context <- make_context(structure)
  prepared <- prepare_data(structure, trial_df)
  params_mat <- build_param_matrix(spec, params, trial_df = prepared)
  as.numeric(log_likelihood(context, prepared, params_mat, min_ll = min_ll))
}

testthat::test_that("exact kernel carries positive-mass tie terms for the guarded shared-gate example", {
  spec <- race_spec() |>
    add_accumulator("go_fast", "lognormal") |>
    add_accumulator("go_slow", "lognormal") |>
    add_accumulator("gate_shared", "lognormal") |>
    add_accumulator("stop_control", "lognormal") |>
    add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
    add_outcome("Slow", all_of("go_slow", "gate_shared"))

  t_obs <- 0.30
  trial_df <- data.frame(
    trial = 1L,
    R = "Fast",
    rt = t_obs,
    stringsAsFactors = FALSE
  )
  params <- c(
    go_fast.m = log(0.28),
    go_fast.s = 0.18,
    go_fast.q = 0.00,
    go_fast.t0 = 0.00,
    go_slow.m = log(0.34),
    go_slow.s = 0.18,
    go_slow.q = 0.00,
    go_slow.t0 = 0.00,
    gate_shared.m = log(0.30),
    gate_shared.s = 0.16,
    gate_shared.q = 0.00,
    gate_shared.t0 = 0.00,
    stop_control.m = log(0.27),
    stop_control.s = 0.15,
    stop_control.q = 0.00,
    stop_control.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  f_fast <- function(x) dlnorm(x, params[["go_fast.m"]], params[["go_fast.s"]])
  F_fast <- function(x) plnorm(x, params[["go_fast.m"]], params[["go_fast.s"]])
  f_slow <- function(x) dlnorm(x, params[["go_slow.m"]], params[["go_slow.s"]])
  F_slow <- function(x) plnorm(x, params[["go_slow.m"]], params[["go_slow.s"]])
  f_gate <- function(x) dlnorm(x, params[["gate_shared.m"]], params[["gate_shared.s"]])
  F_gate <- function(x) plnorm(x, params[["gate_shared.m"]], params[["gate_shared.s"]])
  S_stop <- function(x) 1 - plnorm(x, params[["stop_control.m"]], params[["stop_control.s"]])

  expected <- S_stop(t_obs) * (
    f_fast(t_obs) * F_gate(t_obs) * (1 - F_slow(t_obs)) +
      f_gate(t_obs) * F_fast(t_obs) * (1 - F_slow(t_obs)) +
      f_gate(t_obs) * integrate(
        function(u) f_fast(u) * (F_slow(t_obs) - F_slow(u)),
        lower = 0,
        upper = t_obs,
        rel.tol = 1e-10,
        abs.tol = 0
      )$value
  )

  testthat::expect_equal(out, log(expected), tolerance = 2e-3)
})

testthat::test_that("exact kernel handles nested logical-guard outcomes directly", {
  spec <- race_spec() |>
    add_accumulator("s1", "lognormal") |>
    add_accumulator("is", "lognormal") |>
    add_accumulator("s2", "lognormal") |>
    add_outcome("STOP", all_of("s1", inhibit("s2", by = "is")))

  t_obs <- 0.42
  trial_df <- data.frame(
    trial = 1L,
    R = "STOP",
    rt = t_obs,
    stringsAsFactors = FALSE
  )
  params <- c(
    s1.m = log(0.35), s1.s = 0.16, s1.q = 0.00, s1.t0 = 0.00,
    is.m = log(0.30), is.s = 0.14, is.q = 0.00, is.t0 = 0.00,
    s2.m = log(0.28), s2.s = 0.14, s2.q = 0.00, s2.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  f_s1 <- function(x) dlnorm(x, params[["s1.m"]], params[["s1.s"]])
  F_s1 <- function(x) plnorm(x, params[["s1.m"]], params[["s1.s"]])
  f_s2 <- function(x) dlnorm(x, params[["s2.m"]], params[["s2.s"]])
  S_is <- function(x) 1 - plnorm(x, params[["is.m"]], params[["is.s"]])
  guard_cdf <- integrate(
    function(u) f_s2(u) * S_is(u),
    lower = 0,
    upper = t_obs,
    rel.tol = 1e-10,
    abs.tol = 0
  )$value
  guard_density <- f_s2(t_obs) * S_is(t_obs)
  expected <- f_s1(t_obs) * guard_cdf + F_s1(t_obs) * guard_density

  testthat::expect_equal(out, log(expected), tolerance = 2e-3)
})

testthat::test_that("exact kernel handles none_of as an absence condition inside conjunctions", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("other", "lognormal") |>
    add_outcome("RESPOND", all_of("go", none_of("stop"))) |>
    add_outcome("OTHER", "other")

  t_obs <- 0.41
  trial_df <- data.frame(
    trial = 1L,
    R = "RESPOND",
    rt = t_obs,
    stringsAsFactors = FALSE
  )
  params <- c(
    go.m = log(0.31), go.s = 0.15, go.q = 0.00, go.t0 = 0.00,
    stop.m = log(0.27), stop.s = 0.13, stop.q = 0.00, stop.t0 = 0.00,
    other.m = log(0.45), other.s = 0.17, other.q = 0.00, other.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  expected <- dlnorm(t_obs, params[["go.m"]], params[["go.s"]]) *
    (1 - plnorm(t_obs, params[["stop.m"]], params[["stop.s"]])) *
    (1 - plnorm(t_obs, params[["other.m"]], params[["other.s"]]))

  testthat::expect_equal(out, log(expected), tolerance = 2e-3)
})

testthat::test_that("exact kernel rejects logical-not branches inside first_of/or outcomes", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_outcome("RESPOND", first_of("go", none_of("stop")))

  trial_df <- data.frame(
    trial = 1L,
    R = "RESPOND",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  params <- c(
    go.m = log(0.31), go.s = 0.15, go.q = 0.00, go.t0 = 0.00,
    stop.m = log(0.27), stop.s = 0.13, stop.q = 0.00, stop.t0 = 0.00
  )
  testthat::expect_error(
    run_public_loglik(spec, trial_df, params),
    "logical-not branches inside first_of\\(\\)/or outcomes"
  )
})

testthat::test_that("exact kernel handles overlapping guarded competitors without factorizing them", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("C1", inhibit("a", by = "stop")) |>
    add_outcome("C2", inhibit("b", by = "stop")) |>
    add_outcome("D", "d")

  t_obs <- 0.43
  trial_df <- data.frame(
    trial = 1L,
    R = "D",
    rt = t_obs,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.29), a.s = 0.17, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
    stop.m = log(0.24), stop.s = 0.14, stop.q = 0.00, stop.t0 = 0.00,
    d.m = log(0.46), d.s = 0.19, d.q = 0.00, d.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  expected <- dlnorm(t_obs, params[["d.m"]], params[["d.s"]]) * (
    integrate(
      function(u) {
        dlnorm(u, params[["stop.m"]], params[["stop.s"]]) *
          (1 - plnorm(u, params[["a.m"]], params[["a.s"]])) *
          (1 - plnorm(u, params[["b.m"]], params[["b.s"]]))
      },
      lower = 0,
      upper = t_obs,
      rel.tol = 1e-10,
      abs.tol = 0
    )$value +
      (1 - plnorm(t_obs, params[["stop.m"]], params[["stop.s"]])) *
      (1 - plnorm(t_obs, params[["a.m"]], params[["a.s"]])) *
      (1 - plnorm(t_obs, params[["b.m"]], params[["b.s"]]))
  )

  testthat::expect_equal(out, log(expected), tolerance = 2e-3)
})
