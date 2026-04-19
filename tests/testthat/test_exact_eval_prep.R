.repo_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)

source(file.path(.repo_root, "R", "helpers.R"))
source(file.path(.repo_root, "R", "model_definition.R"))
source(file.path(.repo_root, "R", "semantic_bridge.R"))

testthat::test_that("exact kernel conditions shared-trigger competitors on the realized transition", {
  spec <- race_spec() |>
    add_accumulator("s", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("change", "lognormal") |>
    add_outcome("S", inhibit("s", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_trigger("tg", members = c("stop", "change"), q = 0.05, draw = "shared")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "S",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  params <- c(
    s.m = log(0.28), s.s = 0.12, s.q = 0.00, s.t0 = 0.00,
    stop.m = log(0.35), stop.s = 0.15, stop.t0 = 0.00,
    change.m = log(0.40), change.s = 0.18, change.t0 = 0.00,
    tg = 0.05
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, rebuild = TRUE, root = .repo_root)

  t <- trial_df$rt[[1]]
  fs <- dlnorm(t, params[["s.m"]], params[["s.s"]])
  surv_stop <- 1 - plnorm(t, params[["stop.m"]], params[["stop.s"]])
  expected <- log(params[["tg"]] * fs + (1 - params[["tg"]]) * fs * surv_stop)

  testthat::expect_equal(as.numeric(out$loglik), expected, tolerance = 1e-6)
  testthat::expect_equal(out$total_loglik, expected, tolerance = 1e-6)
})

testthat::test_that("exact kernel matches numeric convolution for single-outcome chained onset", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "B",
    rt = 0.85,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.25), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

  expected <- stats::integrate(
    function(u) {
      stats::dlnorm(u, params[["a.m"]], params[["a.s"]]) *
        stats::dlnorm(trial_df$rt[[1]] - u, params[["b.m"]], params[["b.s"]])
    },
    lower = 0,
    upper = trial_df$rt[[1]],
    rel.tol = 1e-8
  )$value

  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 2e-3)
})

testthat::test_that("exact kernel carries positive-mass tie terms for the guarded shared-gate example", {
  spec <- race_spec() |>
    add_accumulator("go_fast", "lognormal") |>
    add_accumulator("go_slow", "lognormal") |>
    add_accumulator("gate_shared", "lognormal") |>
    add_accumulator("stop_control", "lognormal") |>
    add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
    add_outcome("Slow", all_of("go_slow", "gate_shared"))

  prep <- prepare_model(spec)
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
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

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

  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 2e-3)
})

testthat::test_that("exact kernel handles nested logical-guard outcomes directly", {
  spec <- race_spec() |>
    add_accumulator("s1", "lognormal") |>
    add_accumulator("is", "lognormal") |>
    add_accumulator("s2", "lognormal") |>
    add_outcome("STOP", all_of("s1", inhibit("s2", by = "is")))

  prep <- prepare_model(spec)
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
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

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

  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 2e-3)
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

  prep <- prepare_model(spec)
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
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

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

  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 2e-3)
})

testthat::test_that("exact ranked kernel matches simple independent two-outcome formula", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.30,
    R2 = "B",
    rt2 = 0.55,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.30), a.s = 0.20, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.45), b.s = 0.25, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

  expected <- dlnorm(0.30, params[["a.m"]], params[["a.s"]]) *
    dlnorm(0.55, params[["b.m"]], params[["b.s"]])
  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 1e-4)
})

testthat::test_that("exact ranked kernel treats missing later ranks as truncation", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  ranked_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.30,
    R2 = NA_character_,
    rt2 = NA_real_,
    stringsAsFactors = FALSE
  )
  single_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.30), a.s = 0.20, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.45), b.s = 0.25, b.q = 0.00, b.t0 = 0.00
  )
  ranked_params <- build_param_matrix(spec, params, trial_df = ranked_df)
  single_params <- build_param_matrix(spec, params, trial_df = single_df)

  ranked_out <- .exact_loglik_prep(prep, ranked_params, ranked_df, root = .repo_root)
  single_out <- .exact_loglik_prep(prep, single_params, single_df, root = .repo_root)

  testthat::expect_equal(
    as.numeric(ranked_out$loglik),
    as.numeric(single_out$loglik),
    tolerance = 1e-8
  )
})

testthat::test_that("exact ranked kernel uses observed source time for chained onset", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.30,
    R2 = "B",
    rt2 = 0.60,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.25), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)

  expected <- dlnorm(0.30, params[["a.m"]], params[["a.s"]]) *
    dlnorm(0.60 - 0.30, params[["b.m"]], params[["b.s"]])
  testthat::expect_equal(as.numeric(out$loglik), log(expected), tolerance = 2e-3)
})

testthat::test_that("exact ranked kernel returns min_ll for malformed rank pairs", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.30,
    R2 = "B",
    rt2 = NA_real_,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.30), a.s = 0.20, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.45), b.s = 0.25, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)
  min_ll <- -1e6

  out <- .exact_loglik_prep(
    prep,
    params_mat,
    trial_df,
    min_ll = min_ll,
    root = .repo_root
  )
  testthat::expect_equal(as.numeric(out$loglik), min_ll)
})

testthat::test_that("exact ranked kernel is finite for pool plus shared trigger models", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("A_pool", c("a1", "a2"), k = 1L) |>
    add_outcome("A", "A_pool") |>
    add_outcome("B", "b") |>
    add_trigger("tg", members = c("a1", "a2"), q = 0.10, draw = "shared")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.35,
    R2 = "B",
    rt2 = 0.60,
    stringsAsFactors = FALSE
  )
  params <- c(
    a1.m = log(0.30), a1.s = 0.20, a1.q = 0.10, a1.t0 = 0.00,
    a2.m = log(0.32), a2.s = 0.20, a2.q = 0.10, a2.t0 = 0.00,
    b.m = log(0.50), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)

  out <- .exact_loglik_prep(prep, params_mat, trial_df, root = .repo_root)
  testthat::expect_true(is.finite(as.numeric(out$loglik)))
})

testthat::test_that("exact ranked kernel returns min_ll for non-increasing times", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  trial_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.50,
    R2 = "B",
    rt2 = 0.40,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.25), a.s = 0.12, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.45), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
  )
  params_mat <- build_param_matrix(spec, params, trial_df = trial_df)
  min_ll <- -1e6

  out <- .exact_loglik_prep(
    prep,
    params_mat,
    trial_df,
    min_ll = min_ll,
    root = .repo_root
  )
  testthat::expect_equal(as.numeric(out$loglik), min_ll)
})
