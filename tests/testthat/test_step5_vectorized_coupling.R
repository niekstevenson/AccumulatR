race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
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
