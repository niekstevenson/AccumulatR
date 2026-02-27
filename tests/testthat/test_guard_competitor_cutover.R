race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
all_of <- AccumulatR::all_of
inhibit <- AccumulatR::inhibit
add_pool <- AccumulatR::add_pool
add_trigger <- AccumulatR::add_trigger
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood

testthat::test_that("guard-heavy competitor likelihood remains finite", {
  spec <- race_spec() |>
    add_accumulator("plain", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("PLAIN", "plain") |>
    add_outcome(
      "GUARD",
      inhibit("a", by = inhibit("b", by = inhibit("c", by = "d")))
    )

  structure <- finalize_model(spec)
  params <- c(
    plain.meanlog = log(0.42), plain.sdlog = 0.18,
    a.meanlog = log(0.34), a.sdlog = 0.20,
    b.meanlog = log(0.30), b.sdlog = 0.20,
    c.meanlog = log(0.28), c.sdlog = 0.18,
    d.meanlog = log(0.26), d.sdlog = 0.16
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "PLAIN",
    rt = 0.55,
    stringsAsFactors = FALSE
  )

  ll <- as.numeric(log_likelihood(build_likelihood_context(structure, data_df), params_df))
  testthat::expect_true(is.finite(ll))
})

testthat::test_that("competitor with non-guard wrapper around guard is supported", {
  spec <- race_spec() |>
    add_accumulator("x", "lognormal") |>
    add_accumulator("y", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("R1", "x") |>
    add_outcome("R2", all_of("y", inhibit("a", by = "b")))

  structure <- finalize_model(spec)
  params <- c(
    x.meanlog = log(0.30), x.sdlog = 0.16,
    y.meanlog = log(0.46), y.sdlog = 0.18,
    a.meanlog = log(0.33), a.sdlog = 0.18,
    b.meanlog = log(0.29), b.sdlog = 0.16
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "R1",
    rt = 0.52,
    stringsAsFactors = FALSE
  )

  ll <- as.numeric(log_likelihood(build_likelihood_context(structure, data_df), params_df))
  testthat::expect_true(is.finite(ll))
})

testthat::test_that("ranked shared-trigger path remains finite after cutover", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("A", c("a1", "a2"), k = 1L) |>
    add_outcome("A", "A") |>
    add_outcome("B", "b") |>
    add_trigger("tg", members = c("a1", "a2"), q = 0.10, draw = "shared")

  structure <- finalize_model(spec)
  params <- c(
    a1.meanlog = log(0.30), a1.sdlog = 0.20, a1.q = 0.10, a1.t0 = 0,
    a2.meanlog = log(0.32), a2.sdlog = 0.20, a2.q = 0.10, a2.t0 = 0,
    b.meanlog = log(0.50), b.sdlog = 0.20, b.q = 0.00, b.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "A", rt = 0.35,
    R2 = "B", rt2 = 0.60,
    stringsAsFactors = FALSE
  )

  ll <- as.numeric(log_likelihood(build_likelihood_context(structure, data_df), params_df))
  testthat::expect_true(is.finite(ll))
})
