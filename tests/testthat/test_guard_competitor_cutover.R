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

make_nested_guard_expr <- function(ids) {
  stopifnot(length(ids) >= 2L)
  if (length(ids) == 2L) {
    return(inhibit(ids[[1L]], by = ids[[2L]]))
  }
  inhibit(ids[[1L]], by = make_nested_guard_expr(ids[-1L]))
}

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

testthat::test_that("deep nested guard chain (>8) builds and evaluates", {
  chain_ids <- sprintf("g%02d", seq_len(11L))
  spec <- race_spec()
  for (id in chain_ids) {
    spec <- add_accumulator(spec, id, "lognormal")
  }
  spec <- spec |>
    add_accumulator("plain", "lognormal") |>
    add_outcome("CHAIN", make_nested_guard_expr(chain_ids)) |>
    add_outcome("PLAIN", "plain")

  structure <- testthat::expect_no_error(finalize_model(spec))
  params <- c(
    plain.meanlog = log(0.58), plain.sdlog = 0.18
  )
  for (i in seq_along(chain_ids)) {
    id <- chain_ids[[i]]
    params[paste0(id, ".meanlog")] <- log(0.24 + 0.02 * i)
    params[paste0(id, ".sdlog")] <- 0.14 + 0.002 * i
  }
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "PLAIN",
    rt = 0.72,
    stringsAsFactors = FALSE
  )

  ctx <- testthat::expect_no_error(build_likelihood_context(structure, data_df))
  ll <- as.numeric(log_likelihood(ctx, params_df))
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

testthat::test_that("guard competitor wrapper path is deterministic", {
  spec <- race_spec() |>
    add_accumulator("x", "lognormal") |>
    add_accumulator("y", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("R1", "x") |>
    add_outcome("R2", all_of("y", inhibit("a", by = "b")))

  structure <- finalize_model(spec)
  params <- c(
    x.meanlog = log(0.32), x.sdlog = 0.16,
    y.meanlog = log(0.44), y.sdlog = 0.18,
    a.meanlog = log(0.31), a.sdlog = 0.18,
    b.meanlog = log(0.27), b.sdlog = 0.16
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "R1",
    rt = 0.54,
    stringsAsFactors = FALSE
  )

  ctx <- build_likelihood_context(structure, data_df)
  ll1 <- as.numeric(log_likelihood(ctx, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, params_df))
  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
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

testthat::test_that("guard+competitor shared-trigger path remains finite", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("comp", "lognormal") |>
    add_outcome("R1", inhibit("go", by = "stop")) |>
    add_outcome("R2", "comp") |>
    add_trigger("tg", members = c("go", "comp"), q = 0.08, draw = "shared")

  structure <- finalize_model(spec)
  params <- c(
    go.meanlog = log(0.31), go.sdlog = 0.19, go.q = 0.08, go.t0 = 0,
    stop.meanlog = log(0.40), stop.sdlog = 0.17, stop.q = 0.00, stop.t0 = 0,
    comp.meanlog = log(0.49), comp.sdlog = 0.18, comp.q = 0.08, comp.t0 = 0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "R1",
    rt = 0.34,
    stringsAsFactors = FALSE
  )

  ctx <- build_likelihood_context(structure, data_df)
  ll1 <- as.numeric(log_likelihood(ctx, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, params_df))
  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})
