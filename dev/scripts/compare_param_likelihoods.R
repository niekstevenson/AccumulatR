#!/usr/bin/env Rscript

# Compare legacy vs buffer likelihood pipelines on a synthetic data set.

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  source("examples/new_API.R")
  library(AccumulatR)
})

build_param_table <- function(structure, n_trials) {
  acc_lookup <- structure$accumulators
  param_table <- data.frame(
    trial = rep(seq_len(n_trials), each = nrow(acc_lookup)),
    accumulator_id = rep(acc_lookup$accumulator_id, times = n_trials),
    accumulator = rep(acc_lookup$accumulator_index, times = n_trials),
    component = ifelse(length(structure$components$component_id) == 1L,
                       structure$components$component_id,
                       NA_character_),
    role = rep(acc_lookup$role, times = n_trials),
    onset = rep(acc_lookup$onset, times = n_trials),
    q = rep(acc_lookup$q, times = n_trials),
    stringsAsFactors = FALSE
  )
  param_names <- unique(unlist(lapply(acc_lookup$params, names)))
  if (length(param_names) > 0L) {
    for (param_name in param_names) {
      base_vals <- vapply(acc_lookup$params, function(p) {
        val <- p[[param_name]]
        if (is.null(val)) NA_real_ else as.numeric(val)
      }, numeric(1))
      param_table[[param_name]] <- rep(base_vals, times = n_trials)
    }
  }
  param_table
}

run_speed_test <- function(example_idx = 6,
                           n_trials = 2000L,
                           seed = 2025) {
  set.seed(seed)
  model_spec <- new_api_examples[[example_idx]]
  structure <- build_generator_structure(model_spec)
  model_tables <- model_to_tables(model_spec)
  sim_data <- simulate_model(model_tables, n_trials = n_trials)
  param_table <- build_param_table(structure, n_trials)

  native_bundle <- likelihood_build_native_bundle(model_spec = model_spec)
  prep <- native_bundle$prep

  plan_legacy <- .likelihood_build_trial_plan(structure, param_table, prep)
  plan_buffer <- .likelihood_build_trial_plan(
    structure = structure,
    params_df = param_table,
    prep = prep
  )

  cat(sprintf("Example %d, %d trials\n", example_idx, n_trials))

  legacy_time <- system.time({
    ll_legacy <- log_likelihood_from_params(
      structure = structure,
      params_df = param_table,
      data_df = sim_data,
      prep = prep,
      trial_plan = plan_legacy
    )
  })

  buffer_time <- system.time({
    ll_buffer <- log_likelihood_from_params_buffer(
      structure = structure,
      params_df = param_table,
      data_df = sim_data,
      prep = prep,
      trial_plan = plan_buffer
    )
  })

  if (!isTRUE(all.equal(ll_legacy$loglik, ll_buffer$loglik, tolerance = 5e-3))) {
    warning(sprintf("Log-likelihood mismatch (legacy %.6f vs buffer %.6f)",
                    ll_legacy$loglik, ll_buffer$loglik))
  }

  timings <- rbind(
    data.frame(method = "legacy", elapsed = legacy_time[["elapsed"]],
               user = legacy_time[["user.self"]],
               system = legacy_time[["sys.self"]]),
    data.frame(method = "buffer", elapsed = buffer_time[["elapsed"]],
               user = buffer_time[["user.self"]],
               system = buffer_time[["sys.self"]])
  )

  print(timings)
  invisible(timings)
}

if (identical(environment(), globalenv())) {
  run_speed_test()
}
