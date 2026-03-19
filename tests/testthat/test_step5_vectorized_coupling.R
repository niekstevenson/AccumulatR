race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
add_pool <- AccumulatR::add_pool
all_of <- AccumulatR::all_of
add_trigger <- AccumulatR::add_trigger
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood
response_probabilities <- AccumulatR::response_probabilities

build_shared_gate_pair_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP", all_of("x2", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_))
}

params_shared_gate_pair <- c(
  x1.meanlog = log(0.32), x1.sdlog = 0.18,
  x2.meanlog = log(0.36), x2.sdlog = 0.18,
  gate.meanlog = log(0.24), gate.sdlog = 0.14
)

build_shared_gate_nway_trigger_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
    add_trigger("tg_x12", members = c("x1", "x2"), q = 0.10, draw = "shared") |>
    add_trigger("tg_x3g", members = c("x3", "gate"), q = 0.16, draw = "shared")
}

params_shared_gate_nway_trigger <- c(
  x1.meanlog = log(0.31), x1.sdlog = 0.18, x1.q = 0.10,
  x2.meanlog = log(0.34), x2.sdlog = 0.18, x2.q = 0.10,
  x3.meanlog = log(0.37), x3.sdlog = 0.18, x3.q = 0.16,
  gate.meanlog = log(0.23), gate.sdlog = 0.14, gate.q = 0.16
)

build_shared_gate_guess_trigger_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_accumulator("timeout", "lognormal") |>
    add_outcome("RESP", all_of("x2", "gate")) |>
    add_outcome("ALT", "x1") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(labels = c("RESP"), weights = c(1), rt_policy = "keep")
    )) |>
    add_trigger("tg_all", members = c("x1", "x2", "gate", "timeout"),
      q = 0.10, param = "tg_q", draw = "shared"
    )
}

params_shared_gate_guess_trigger <- c(
  x1.meanlog = log(0.31), x1.sdlog = 0.18,
  x2.meanlog = log(0.34), x2.sdlog = 0.18,
  gate.meanlog = log(0.23), gate.sdlog = 0.14,
  timeout.meanlog = log(0.46), timeout.sdlog = 0.12,
  tg_q = 0.25
)

build_shared_gate_pool_trigger_spec <- function() {
  race_spec() |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_pool("A", c("a1", "a2"), k = 1L) |>
    add_outcome("RESP", all_of("A", "gate")) |>
    add_outcome("NR_RAW", all_of("b", "gate"),
      options = list(map_outcome_to = NA_character_)
    ) |>
    add_trigger("tg_pool", members = c("a1", "a2", "gate"),
      q = 0.10, param = "tg_q", draw = "shared"
    )
}

params_shared_gate_pool_trigger <- c(
  a1.meanlog = log(0.28), a1.sdlog = 0.16,
  a2.meanlog = log(0.33), a2.sdlog = 0.16,
  b.meanlog = log(0.31), b.sdlog = 0.17,
  gate.meanlog = log(0.24), gate.sdlog = 0.14,
  tg_q = 0.22
)

testthat::test_that("shared-gate pair NA mass is deterministic and normalized", {
  spec <- build_shared_gate_pair_spec()
  structure <- finalize_model(spec)
  params_df <- build_param_matrix(spec, params_shared_gate_pair, n_trials = 1L)

  probs1 <- response_probabilities(structure, params_df, include_na = TRUE)
  probs2 <- response_probabilities(structure, params_df, include_na = TRUE)

  testthat::expect_equal(probs1, probs2, tolerance = 1e-12)
  testthat::expect_true(all(is.finite(as.numeric(probs1))))
  testthat::expect_true("NA" %in% names(probs1))
  testthat::expect_equal(sum(as.numeric(probs1)), 1.0, tolerance = 1e-8)
})

testthat::test_that("shared-gate nway with shared triggers is deterministic", {
  spec <- build_shared_gate_nway_trigger_spec()
  structure <- finalize_model(spec)
  params_df <- build_param_matrix(spec, params_shared_gate_nway_trigger, n_trials = 1L)

  probs1 <- response_probabilities(structure, params_df, include_na = TRUE)
  probs2 <- response_probabilities(structure, params_df, include_na = TRUE)

  testthat::expect_equal(probs1, probs2, tolerance = 1e-12)
  testthat::expect_true(all(is.finite(as.numeric(probs1))))
  testthat::expect_true("NA" %in% names(probs1))
  testthat::expect_equal(sum(as.numeric(probs1)), 1.0, tolerance = 1e-8)
})

testthat::test_that("shared-gate pair nonresponse log-likelihood is stable", {
  spec <- build_shared_gate_pair_spec()
  structure <- finalize_model(spec)
  params_df <- build_param_matrix(spec, params_shared_gate_pair, n_trials = 1L)
  nonresp_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("RESP", "NR_RAW")),
    rt = NA_real_
  )
  ctx <- build_likelihood_context(structure, nonresp_df)

  ll1 <- as.numeric(log_likelihood(ctx, nonresp_df, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, nonresp_df, params_df))

  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})

testthat::test_that("repeated nonresponse trials use outcome-mass batch in cpp_loglik", {
  spec <- build_shared_gate_pair_spec()
  structure <- finalize_model(spec)

  one_trial_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("RESP", "NR_RAW")),
    rt = NA_real_
  )
  repeated_df <- data.frame(
    trial = 1:3,
    R = factor(rep(NA_character_, 3), levels = c("RESP", "NR_RAW")),
    rt = rep(NA_real_, 3)
  )

  ll_one <- as.numeric(log_likelihood(
    build_likelihood_context(structure, one_trial_df),
    one_trial_df,
    build_param_matrix(spec, params_shared_gate_pair, n_trials = 1L)
  ))

  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll_repeated <- as.numeric(log_likelihood(
    build_likelihood_context(structure, repeated_df),
    repeated_df,
    build_param_matrix(spec, params_shared_gate_pair, n_trials = 3L)
  ))
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  testthat::expect_true(is.finite(ll_one))
  testthat::expect_true(is.finite(ll_repeated))
  testthat::expect_equal(ll_repeated, 3 * ll_one, tolerance = 1e-12)
  testthat::expect_gt(stats$cpp_loglik_outcome_mass_group_batch_calls, 0)
  testthat::expect_gt(
    stats$cpp_loglik_outcome_mass_group_batch_trials_total,
    stats$cpp_loglik_outcome_mass_group_batch_calls
  )
  testthat::expect_equal(stats$cpp_loglik_fallback_group_calls, 0)
  testthat::expect_equal(stats$cpp_loglik_outcome_mass_group_batch_exec_failures, 0)
})

testthat::test_that("non-shared-trigger outcome mass batches trials with different params", {
  spec <- build_shared_gate_pair_spec()
  structure <- finalize_model(spec)

  nonresp_df <- data.frame(
    trial = 1:3,
    R = factor(rep(NA_character_, 3), levels = c("RESP", "NR_RAW")),
    rt = rep(NA_real_, 3)
  )
  params_df <- build_param_matrix(spec, params_shared_gate_pair, n_trials = 3L)
  rows_per_trial <- nrow(params_df) %/% 3L

  testthat::expect_equal(nrow(params_df), rows_per_trial * 3L)

  p1_offsets <- c(0.00, 0.08, -0.05)
  for (trial_idx in seq_along(p1_offsets)) {
    row_idx <- ((trial_idx - 1L) * rows_per_trial + 1L):(trial_idx * rows_per_trial)
    params_df[row_idx, "p1"] <- params_df[row_idx, "p1"] + p1_offsets[trial_idx]
  }

  ctx_batch <- build_likelihood_context(structure, nonresp_df)
  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll_batch <- as.numeric(log_likelihood(ctx_batch, nonresp_df, params_df))
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  split_ll <- 0
  one_trial_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("RESP", "NR_RAW")),
    rt = NA_real_
  )
  for (trial_idx in seq_along(p1_offsets)) {
    row_idx <- ((trial_idx - 1L) * rows_per_trial + 1L):(trial_idx * rows_per_trial)
    split_ll <- split_ll + as.numeric(log_likelihood(
      build_likelihood_context(structure, one_trial_df),
      one_trial_df,
      params_df[row_idx, , drop = FALSE]
    ))
  }

  testthat::expect_true(is.finite(ll_batch))
  testthat::expect_equal(ll_batch, split_ll, tolerance = 1e-10)
  testthat::expect_gt(stats$cpp_loglik_outcome_mass_group_batch_calls, 0)
  testthat::expect_gt(
    stats$cpp_loglik_outcome_mass_group_batch_trials_total,
    stats$cpp_loglik_outcome_mass_group_batch_calls
  )
  testthat::expect_equal(stats$cpp_loglik_outcome_mass_group_batch_exec_failures, 0)
  testthat::expect_equal(stats$cpp_loglik_fallback_group_calls, 0)
})

testthat::test_that("shared-trigger nonresponse batches outcome mass without scalar labelref fallback", {
  spec <- build_shared_gate_nway_trigger_spec()
  structure <- finalize_model(spec)
  params_df <- build_param_matrix(spec, params_shared_gate_nway_trigger, n_trials = 1L)
  nonresp_df <- data.frame(
    trial = 1L,
    R = factor(NA_character_, levels = c("RESP2", "RESP3", "NR_RAW")),
    rt = NA_real_
  )
  ctx <- build_likelihood_context(structure, nonresp_df)

  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll <- as.numeric(log_likelihood(ctx, nonresp_df, params_df))
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  testthat::expect_true(is.finite(ll))
  testthat::expect_gt(stats$evaluate_outcome_coupling_unified_calls, 0)
  testthat::expect_gt(stats$shared_trigger_mask_batch_calls, 0)
})

testthat::test_that("shared-trigger guess-donor response probabilities batch masks without scalar fallback", {
  spec <- build_shared_gate_guess_trigger_spec()
  structure <- finalize_model(spec)

  params_shared <- params_shared_gate_guess_trigger
  params_success <- params_shared
  params_fail <- params_shared
  params_success[["tg_q"]] <- 0
  params_fail[["tg_q"]] <- 1

  probs_success <- response_probabilities(
    structure,
    build_param_matrix(spec, params_success, n_trials = 1L),
    include_na = TRUE
  )
  probs_fail <- response_probabilities(
    structure,
    build_param_matrix(spec, params_fail, n_trials = 1L),
    include_na = TRUE
  )

  AccumulatR:::unified_outcome_stats_reset_cpp()
  probs_shared <- response_probabilities(
    structure,
    build_param_matrix(spec, params_shared, n_trials = 1L),
    include_na = TRUE
  )
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  observed_names <- c("RESP", "ALT", "TIMEOUT")
  expected <- (1 - params_shared[["tg_q"]]) * probs_success[observed_names] +
    params_shared[["tg_q"]] * probs_fail[observed_names]

  testthat::expect_equal(probs_shared[observed_names], expected, tolerance = 1e-8)
  testthat::expect_gt(stats$evaluate_outcome_coupling_unified_calls, 0)
  testthat::expect_gt(stats$shared_trigger_mask_batch_calls, 0)
})

testthat::test_that("shared-trigger pooled coupling uses pool labelref batch without kernel-event pointwise work", {
  spec <- build_shared_gate_pool_trigger_spec()
  structure <- finalize_model(spec)

  params_shared <- params_shared_gate_pool_trigger
  params_success <- params_shared
  params_fail <- params_shared
  params_success[["tg_q"]] <- 0
  params_fail[["tg_q"]] <- 1

  probs_success <- response_probabilities(
    structure,
    build_param_matrix(spec, params_success, n_trials = 1L),
    include_na = TRUE
  )
  probs_fail <- response_probabilities(
    structure,
    build_param_matrix(spec, params_fail, n_trials = 1L),
    include_na = TRUE
  )

  AccumulatR:::unified_outcome_stats_reset_cpp()
  probs_shared <- response_probabilities(
    structure,
    build_param_matrix(spec, params_shared, n_trials = 1L),
    include_na = TRUE
  )
  stats <- AccumulatR:::unified_outcome_stats_cpp()

  expected <- (1 - params_shared[["tg_q"]]) * probs_success + params_shared[["tg_q"]] * probs_fail

  testthat::expect_equal(probs_shared, expected, tolerance = 1e-7)
  testthat::expect_gt(stats$evaluate_outcome_coupling_unified_calls, 0)
  testthat::expect_gt(stats$shared_trigger_mask_batch_calls, 0)
  testthat::expect_gt(stats$generic_poolref_batch_fastpath_calls, 0)
  testthat::expect_equal(stats$kernel_event_batch_calls, 0)
})
