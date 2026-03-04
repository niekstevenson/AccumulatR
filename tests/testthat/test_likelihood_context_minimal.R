testthat::test_that("likelihood context stores only compact runtime fields", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B")
  structure <- finalize_model(spec)

  params <- build_param_matrix(
    spec,
    c(
      A.meanlog = log(0.30), A.sdlog = 0.20, A.q = 0, A.t0 = 0,
      B.meanlog = log(0.45), B.sdlog = 0.25, B.q = 0, B.t0 = 0
    ),
    n_trials = 3L
  )
  data_df <- simulate(structure, params, seed = 1)

  ctx <- build_likelihood_context(structure, data_df)
  expected_names <- c(
    "native_ctx", "accumulator_ids", "accumulator_onset_map", "component_ids",
    "accumulators_by_component",
    "rel_tol", "abs_tol", "max_depth"
  )
  testthat::expect_setequal(names(ctx), expected_names)
  testthat::expect_false(any(c(
    "data_df", "data_df_cpp", "prep", "structure", "trial_ids",
    "rank_width", "data", "param_layout", "n_trials"
  ) %in% names(ctx)))

  ll <- as.numeric(log_likelihood(ctx, data_df, params))
  testthat::expect_true(is.finite(ll))
})

testthat::test_that("ensure_native_ctx requires model_spec for rebuilds", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_outcome("A", "A")
  structure <- finalize_model(spec)

  params <- build_param_matrix(
    spec,
    c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
    n_trials = 1L
  )
  data_df <- simulate(structure, params, seed = 1)
  ctx <- build_likelihood_context(structure, data_df)

  ctx$native_ctx <- NULL
  testthat::expect_error(
    ensure_native_ctx(ctx),
    "model_spec required"
  )
  rebuilt <- ensure_native_ctx(ctx, model_spec = structure$model_spec)
  testthat::expect_true(inherits(rebuilt$native_ctx, "externalptr"))
})

testthat::test_that("build_likelihood_context validates factor levels once at creation", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B")
  structure <- finalize_model(spec)

  good_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("A", "B")),
    rt = 0.35
  )
  testthat::expect_no_error(build_likelihood_context(structure, good_df))

  bad_order_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = c("B", "A")),
    rt = 0.35
  )
  testthat::expect_error(
    build_likelihood_context(structure, bad_order_df),
    "exact model level order"
  )

  bad_type_df <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.35
  )
  testthat::expect_error(
    build_likelihood_context(structure, bad_type_df),
    "must be factor"
  )
})
