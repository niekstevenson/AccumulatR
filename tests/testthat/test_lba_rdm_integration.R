testthat::test_that("sampled_pars exposes canonical LBA/RDM parameter names", {
  spec <- race_spec() |>
    add_accumulator("a", "LBA") |>
    add_accumulator("b", "RDM") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  pars <- sampled_pars(spec)
  expected <- c(
    "a.v", "a.sv", "a.B", "a.A", "a.q", "a.t0",
    "b.v", "b.B", "b.A", "b.s", "b.q", "b.t0"
  )
  testthat::expect_setequal(pars, expected)
})

testthat::test_that("build_param_matrix maps LBA/RDM fourth parameter into p4", {
  spec <- race_spec() |>
    add_accumulator("a", "lba") |>
    add_accumulator("b", "rdm") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  vals <- c(
    a.v = 2.0, a.sv = 0.6, a.B = 1.2, a.A = 0.4, a.q = 0.0, a.t0 = 0.05,
    b.v = 1.5, b.B = 1.0, b.A = 0.3, b.s = 1.1, b.q = 0.0, b.t0 = 0.02
  )
  pm <- build_param_matrix(spec, vals, n_trials = 1L)

  testthat::expect_true(all(c("p1", "p2", "p3", "p4") %in% colnames(pm)))
  testthat::expect_equal(as.numeric(pm[1, "p4"]), vals[["a.A"]])
  testthat::expect_equal(as.numeric(pm[2, "p4"]), vals[["b.s"]])
})

testthat::test_that("LBA model supports simulate + likelihood end-to-end", {
  spec <- race_spec() |>
    add_accumulator("a", "lba") |>
    add_outcome("A", "a")
  structure <- finalize_model(spec)

  vals <- c(a.v = 2.0, a.sv = 0.6, a.B = 1.2, a.A = 0.4, a.q = 0.0, a.t0 = 0.05)
  params_df <- build_param_matrix(spec, vals, n_trials = 12L)
  data_df <- simulate(structure, params_df, seed = 1)
  ctx <- build_likelihood_context(structure, data_df)
  ll <- as.numeric(log_likelihood(ctx, params_df))

  testthat::expect_true(length(ll) == 1L)
  testthat::expect_true(is.finite(ll))
  probs <- response_probabilities(structure, build_param_matrix(spec, vals, n_trials = 1L))
  testthat::expect_true(all(is.finite(probs)))
})

testthat::test_that("RDM model supports simulate + likelihood end-to-end", {
  spec <- race_spec() |>
    add_accumulator("a", "rdm") |>
    add_outcome("A", "a")
  structure <- finalize_model(spec)

  vals <- c(a.v = 1.5, a.B = 1.0, a.A = 0.3, a.s = 1.1, a.q = 0.0, a.t0 = 0.02)
  params_df <- build_param_matrix(spec, vals, n_trials = 12L)
  data_df <- simulate(structure, params_df, seed = 2)
  ctx <- build_likelihood_context(structure, data_df)
  ll <- as.numeric(log_likelihood(ctx, params_df))

  testthat::expect_true(length(ll) == 1L)
  testthat::expect_true(is.finite(ll))
  probs <- response_probabilities(structure, build_param_matrix(spec, vals, n_trials = 1L))
  testthat::expect_true(all(is.finite(probs)))
})

testthat::test_that("RDM requires parameter s", {
  spec <- race_spec() |>
    add_accumulator("a", "rdm") |>
    add_outcome("A", "a")
  vals <- c(a.v = 1.5, a.B = 1.0, a.A = 0.3, a.q = 0.0, a.t0 = 0.02)

  testthat::expect_error(
    build_param_matrix(spec, vals, n_trials = 1L),
    "Missing required parameter 's'"
  )
})
