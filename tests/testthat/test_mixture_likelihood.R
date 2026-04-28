latent_sampled_mixture_spec <- function() {
  race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("TARGET", c("target_fast", "target_slow")) |>
    add_outcome("Target", "TARGET") |>
    add_outcome("Competitor", "competitor") |>
    add_component(
      "fast",
      members = c("target_fast", "competitor"),
      weight_param = "p_fast"
    ) |>
    add_component("slow", members = c("target_slow", "competitor")) |>
    set_mixture_options(mode = "sample", reference = "slow") |>
    finalize_model()
}

latent_sampled_mixture_params <- function(p_fast) {
  c(
    target_fast.m = log(0.25),
    target_fast.s = 0.15,
    target_slow.m = log(0.45),
    target_slow.s = 0.20,
    competitor.m = log(0.35),
    competitor.s = 0.18,
    p_fast = p_fast
  )
}

testthat::test_that("latent sampled mixtures marginalize on the C++ likelihood path", {
  structure <- latent_sampled_mixture_spec()
  ctx <- make_context(structure)

  latent_df <- data.frame(
    trial = 1L,
    R = "Target",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  latent_prepared <- prepare_data(structure, latent_df)

  explicit_na_df <- data.frame(
    trial = 1L,
    component = NA_character_,
    R = "Target",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  explicit_na_prepared <- prepare_data(structure, explicit_na_df)

  fast_df <- data.frame(
    trial = 1L,
    component = "fast",
    R = "Target",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  slow_df <- data.frame(
    trial = 1L,
    component = "slow",
    R = "Target",
    rt = 0.30,
    stringsAsFactors = FALSE
  )
  fast_prepared <- prepare_data(structure, fast_df)
  slow_prepared <- prepare_data(structure, slow_df)

  params_lo <- latent_sampled_mixture_params(0.2)
  latent_mat_lo <- build_param_matrix(
    structure$model_spec,
    params_lo,
    trial_df = latent_prepared
  )
  fast_mat_lo <- build_param_matrix(
    structure$model_spec,
    params_lo,
    trial_df = fast_prepared
  )
  slow_mat_lo <- build_param_matrix(
    structure$model_spec,
    params_lo,
    trial_df = slow_prepared
  )

  ll_latent_lo <- as.numeric(log_likelihood(ctx, latent_prepared, latent_mat_lo))
  ll_explicit_na_lo <- as.numeric(log_likelihood(
    ctx,
    explicit_na_prepared,
    build_param_matrix(structure$model_spec, params_lo, trial_df = explicit_na_prepared)
  ))
  ll_fast_lo <- as.numeric(log_likelihood(ctx, fast_prepared, fast_mat_lo))
  ll_slow_lo <- as.numeric(log_likelihood(ctx, slow_prepared, slow_mat_lo))

  expected_lo <- log(
    0.2 * exp(ll_fast_lo) +
      0.8 * exp(ll_slow_lo)
  )

  testthat::expect_equal(ll_latent_lo, expected_lo, tolerance = 1e-8)
  testthat::expect_equal(ll_explicit_na_lo, ll_latent_lo, tolerance = 1e-8)

  params_hi <- latent_sampled_mixture_params(0.8)
  ll_latent_hi <- as.numeric(log_likelihood(
    ctx,
    latent_prepared,
    build_param_matrix(structure$model_spec, params_hi, trial_df = latent_prepared)
  ))
  ll_fast_hi <- as.numeric(log_likelihood(
    ctx,
    fast_prepared,
    build_param_matrix(structure$model_spec, params_hi, trial_df = fast_prepared)
  ))

  testthat::expect_gt(abs(ll_latent_hi - ll_latent_lo), 1e-6)
  testthat::expect_equal(ll_fast_hi, ll_fast_lo, tolerance = 1e-10)
})

testthat::test_that("latent fixed mixtures use fixed component weights in likelihood", {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Left", "A") |>
    add_outcome("Right", "B") |>
    add_component("left_only", members = "A", weight = 0.25) |>
    add_component("right_only", members = "B", weight = 0.75) |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()

  params <- c(
    A.m = log(0.30), A.s = 0.15, A.q = 0, A.t0 = 0,
    B.m = log(0.40), B.s = 0.15, B.q = 0, B.t0 = 0
  )
  ctx <- make_context(structure)

  latent_prepared <- prepare_data(
    structure,
    data.frame(
      trial = 1L,
      R = "Left",
      rt = 0.32,
      stringsAsFactors = FALSE
    )
  )
  left_prepared <- prepare_data(
    structure,
    data.frame(
      trial = 1L,
      component = "left_only",
      R = "Left",
      rt = 0.32,
      stringsAsFactors = FALSE
    )
  )

  ll_latent <- as.numeric(log_likelihood(
    ctx,
    latent_prepared,
    build_param_matrix(structure$model_spec, params, trial_df = latent_prepared)
  ))
  ll_left <- as.numeric(log_likelihood(
    ctx,
    left_prepared,
    build_param_matrix(structure$model_spec, params, trial_df = left_prepared)
  ))

  testthat::expect_equal(ll_latent, log(0.25) + ll_left, tolerance = 1e-8)
})

testthat::test_that("latent missing-all mixtures use compiled terminal no-response when available", {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Left", "A") |>
    add_outcome("Right", "B") |>
    add_component("left_only", members = "A", weight = 0.25) |>
    add_component("right_only", members = "B", weight = 0.75) |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()

  prepared <- prepare_data(
    structure,
    data.frame(
      trial = 1L,
      R = NA_character_,
      rt = NA_real_,
      stringsAsFactors = FALSE
    )
  )
  params <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.20,
    B.m = log(0.32), B.s = 0.18, B.q = 0.40
  )

  out <- as.numeric(log_likelihood(
    make_context(structure),
    prepared,
    build_param_matrix(structure$model_spec, params, trial_df = prepared)
  ))

  testthat::expect_equal(out, log(0.25 * 0.20 + 0.75 * 0.40), tolerance = 1e-10)
})

testthat::test_that("observed labels with missing RT marginalize over finite response", {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Seen", "A") |>
    add_outcome("HiddenSeen", "B", options = list(
      guess = list(labels = "Seen", weights = 1, rt_policy = "na")
    )) |>
    finalize_model()

  prepared <- prepare_data(
    structure,
    data.frame(
      trial = 1L,
      R = "Seen",
      rt = NA_real_,
      stringsAsFactors = FALSE
    )
  )
  params <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.20, A.t0 = 0,
    B.m = log(0.32), B.s = 0.18, B.q = 0.30, B.t0 = 0
  )

  out <- as.numeric(log_likelihood(
    make_context(structure),
    prepared,
    build_param_matrix(structure$model_spec, params, trial_df = prepared)
  ))

  testthat::expect_equal(out, log(1 - 0.20 * 0.30), tolerance = 1e-10)
})
