testthat::test_that("after() onset is accepted and metadata is recorded", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal", onset = 0.02) |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  finalized <- finalize_model(spec)

  testthat::expect_equal(finalized$prep$onset_specs$a$kind, "absolute")
  testthat::expect_equal(finalized$prep$onset_specs$a$value, 0.02)
  testthat::expect_equal(finalized$prep$onset_specs$b$kind, "after")
  testthat::expect_equal(finalized$prep$onset_specs$b$source, "a")
  testthat::expect_equal(finalized$prep$onset_specs$b$source_kind, "accumulator")
  testthat::expect_equal(finalized$prep$onset_specs$b$lag, 0)
  testthat::expect_equal(finalized$prep$onset_dependencies$b, "a")
  testthat::expect_equal(finalized$prep$onset_topology, c("a", "b"))
  testthat::expect_true(finalized$prep$onset_has_dependencies)
})

testthat::test_that("after(pool, lag) is accepted and pool dependencies expand", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_accumulator("c", "lognormal", onset = after("p", lag = 0.05)) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  finalized <- finalize_model(spec)

  testthat::expect_equal(finalized$prep$onset_specs$c$kind, "after")
  testthat::expect_equal(finalized$prep$onset_specs$c$source, "p")
  testthat::expect_equal(finalized$prep$onset_specs$c$source_kind, "pool")
  testthat::expect_equal(finalized$prep$onset_specs$c$lag, 0.05)
  testthat::expect_setequal(finalized$prep$onset_dependencies$c, c("a", "b"))
})

testthat::test_that("onset character shorthand is explicitly rejected", {
  testthat::expect_error(
    race_spec() |> add_accumulator("a", "lognormal", onset = "b"),
    "Character onset shorthand is not supported"
  )
})

testthat::test_that("invalid chained onset sources are rejected", {
  unknown_source <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("missing")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  testthat::expect_error(
    finalize_model(unknown_source),
    "Unknown onset source 'missing'"
  )

  outcome_source <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("RESP")) |>
    add_outcome("RESP", "a") |>
    add_outcome("B", "b")

  testthat::expect_error(
    finalize_model(outcome_source),
    "references an outcome label"
  )
})

testthat::test_that("after() validates lag and malformed literals", {
  testthat::expect_error(after("a", lag = -0.01), ">= 0")
  testthat::expect_error(after("a", lag = Inf), ">= 0")

  malformed <- structure(
    list(kind = "after", source = "a", lag = -1),
    class = c("race_onset_after", "list")
  )
  testthat::expect_error(
    race_spec() |> add_accumulator("a", "lognormal", onset = malformed),
    "single finite numeric value >= 0"
  )
})

testthat::test_that("chained onset cycles are rejected", {
  self_cycle <- race_spec() |>
    add_accumulator("a", "lognormal", onset = after("a")) |>
    add_outcome("A", "a")
  testthat::expect_error(
    finalize_model(self_cycle),
    "Chained onset dependency cycle detected"
  )

  two_node_cycle <- race_spec() |>
    add_accumulator("a", "lognormal", onset = after("b")) |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")
  testthat::expect_error(
    finalize_model(two_node_cycle),
    "Chained onset dependency cycle detected"
  )

  pool_cycle <- race_spec() |>
    add_accumulator("a", "lognormal", onset = after("p")) |>
    add_accumulator("b", "lognormal") |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")
  testthat::expect_error(
    finalize_model(pool_cycle),
    "Chained onset dependency cycle detected"
  )
})

testthat::test_that("numeric onset models remain compatible", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal", onset = 0.01) |>
    add_accumulator("b", "lognormal", onset = 0.03) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  finalized <- finalize_model(spec)

  testthat::expect_false(finalized$prep$onset_has_dependencies)
  testthat::expect_equal(finalized$prep$onset_specs$a$kind, "absolute")
  testthat::expect_equal(finalized$prep$onset_specs$b$kind, "absolute")
  testthat::expect_equal(finalized$prep$onset_specs$a$value, 0.01)
  testthat::expect_equal(finalized$prep$onset_specs$b$value, 0.03)
})

testthat::test_that("simulate honors accumulator chained onset delays", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal", onset = after("a", lag = 0.05)) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  finalized <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.20), a.sdlog = 0.05, a.q = 0.0, a.t0 = 0.0,
    b.meanlog = log(0.25), b.sdlog = 0.05, b.q = 0.0, b.t0 = 0.0,
    c.meanlog = log(0.10), c.sdlog = 0.05, c.q = 0.0, c.t0 = 0.0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 64L)
  sim <- simulate(finalized, params_df, seed = 42, keep_detail = TRUE)
  details <- attr(sim, "details")

  acc_time <- function(detail, id) {
    val <- detail$acc_times[[id]]
    if (is.null(val) || length(val) == 0L) return(NA_real_)
    as.numeric(val[[1]])
  }
  lag_ok <- vapply(details, function(detail) {
    a_time <- acc_time(detail, "a")
    c_time <- acc_time(detail, "c")
    is.finite(a_time) && is.finite(c_time) && (c_time >= (a_time + 0.05 - 1e-12))
  }, logical(1))

  testthat::expect_equal(nrow(sim), 64L)
  testthat::expect_true(all(lag_ok))
})

testthat::test_that("simulate honors pool chained onset and missing source completion", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_accumulator("c", "lognormal", onset = after("p", lag = 0.02)) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")

  finalized <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.18), a.sdlog = 0.05, a.q = 0.0, a.t0 = 0.0,
    b.meanlog = log(0.22), b.sdlog = 0.05, b.q = 0.0, b.t0 = 0.0,
    c.meanlog = log(0.10), c.sdlog = 0.05, c.q = 0.0, c.t0 = 0.0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 64L)
  sim <- simulate(finalized, params_df, seed = 99, keep_detail = TRUE)
  details <- attr(sim, "details")

  acc_time <- function(detail, id) {
    val <- detail$acc_times[[id]]
    if (is.null(val) || length(val) == 0L) return(NA_real_)
    as.numeric(val[[1]])
  }
  pool_lag_ok <- vapply(details, function(detail) {
    a_time <- acc_time(detail, "a")
    b_time <- acc_time(detail, "b")
    c_time <- acc_time(detail, "c")
    src <- min(c(a_time, b_time))
    is.finite(src) && is.finite(c_time) && (c_time >= (src + 0.02 - 1e-12))
  }, logical(1))
  testthat::expect_true(all(pool_lag_ok))

  spec_missing <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")
  finalized_missing <- finalize_model(spec_missing)
  params_missing <- c(
    a.meanlog = log(0.20), a.sdlog = 0.05, a.q = 1.0, a.t0 = 0.0,
    b.meanlog = log(0.25), b.sdlog = 0.05, b.q = 0.0, b.t0 = 0.0,
    c.meanlog = log(0.10), c.sdlog = 0.05, c.q = 0.0, c.t0 = 0.0
  )
  params_df_missing <- build_param_matrix(spec_missing, params_missing, n_trials = 32L)
  sim_missing <- simulate(finalized_missing, params_df_missing, seed = 7, keep_detail = TRUE)
  details_missing <- attr(sim_missing, "details")
  c_times <- vapply(details_missing, function(detail) acc_time(detail, "c"), numeric(1))
  testthat::expect_true(all(is.infinite(c_times)))
  testthat::expect_true(all(is.na(sim_missing$R) | sim_missing$R != "C"))
})

testthat::test_that("chained onset likelihood context and evaluation are available", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")
  finalized <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.25), a.sdlog = 0.10, a.q = 0.0, a.t0 = 0.0,
    b.meanlog = log(0.20), b.sdlog = 0.12, b.q = 0.0, b.t0 = 0.0
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  data_df <- data.frame(
    trial = 1L,
    R = "B",
    rt = 0.75,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(finalized, data_df)
  ll <- as.numeric(log_likelihood(ctx, params_df))
  testthat::expect_true(is.finite(ll))
})
