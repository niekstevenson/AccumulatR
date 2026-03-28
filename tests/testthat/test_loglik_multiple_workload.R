race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
add_pool <- AccumulatR::add_pool
add_trigger <- AccumulatR::add_trigger
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood

testthat::test_that("multi-parameter direct evaluation matches repeated scalar calls", {
  spec <- race_spec() |>
    add_accumulator("go", "lognormal") |>
    add_accumulator("alt", "lognormal") |>
    add_outcome("GO", "go") |>
    add_outcome("ALT", "alt")

  structure <- finalize_model(spec)
  data_df <- data.frame(
    trial = 1:3,
    R = factor(c("GO", "ALT", "GO"), levels = c("GO", "ALT")),
    rt = c(0.34, 0.47, 0.39)
  )
  ctx <- build_likelihood_context(structure, data_df)

  params1 <- c(
    go.meanlog = log(0.30), go.sdlog = 0.18,
    alt.meanlog = log(0.45), alt.sdlog = 0.19
  )
  params2 <- c(
    go.meanlog = log(0.33), go.sdlog = 0.16,
    alt.meanlog = log(0.49), alt.sdlog = 0.17
  )
  pm1 <- build_param_matrix(spec, params1, n_trials = 3L)
  pm2 <- build_param_matrix(spec, params2, n_trials = 3L)

  ok <- c(TRUE, FALSE, TRUE)
  expand <- c(1L, 3L)

  ll_multi <- as.numeric(log_likelihood(
    ctx,
    data_df,
    list(pm1, pm2),
    ok = ok,
    expand = expand
  ))
  ll_split <- c(
    as.numeric(log_likelihood(ctx, data_df, pm1, ok = ok, expand = expand)),
    as.numeric(log_likelihood(ctx, data_df, pm2, ok = ok, expand = expand))
  )

  testthat::expect_equal(ll_multi, ll_split, tolerance = 1e-12)
})

testthat::test_that("multi-parameter ranked evaluation matches repeated scalar calls", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("A", c("a1", "a2"), k = 1L) |>
    add_outcome("A", "A") |>
    add_outcome("B", "b") |>
    add_trigger("tg", members = c("a1", "a2"), q = 0.10, draw = "shared")

  structure <- finalize_model(spec)
  ranked_df <- data.frame(
    trial = 1:2,
    R = factor(c("A", "A"), levels = c("A", "B")),
    rt = c(0.36, 0.38),
    R2 = factor(c("B", "B"), levels = c("A", "B")),
    rt2 = c(0.61, 0.64)
  )
  ctx <- build_likelihood_context(structure, ranked_df)

  params1 <- c(
    a1.meanlog = log(0.30), a1.sdlog = 0.20, a1.q = 0.10, a1.t0 = 0.00,
    a2.meanlog = log(0.33), a2.sdlog = 0.20, a2.q = 0.10, a2.t0 = 0.00,
    b.meanlog = log(0.52), b.sdlog = 0.18, b.q = 0.00, b.t0 = 0.00
  )
  params2 <- c(
    a1.meanlog = log(0.31), a1.sdlog = 0.18, a1.q = 0.10, a1.t0 = 0.00,
    a2.meanlog = log(0.35), a2.sdlog = 0.18, a2.q = 0.10, a2.t0 = 0.00,
    b.meanlog = log(0.56), b.sdlog = 0.17, b.q = 0.00, b.t0 = 0.00
  )
  pm1 <- build_param_matrix(spec, params1, n_trials = 2L)
  pm2 <- build_param_matrix(spec, params2, n_trials = 2L)

  ll_multi <- as.numeric(log_likelihood(ctx, ranked_df, list(pm1, pm2)))
  ll_split <- c(
    as.numeric(log_likelihood(ctx, ranked_df, pm1)),
    as.numeric(log_likelihood(ctx, ranked_df, pm2))
  )

  testthat::expect_equal(ll_multi, ll_split, tolerance = 1e-10)
})
