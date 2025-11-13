#!/usr/bin/env Rscript

# Minimal example showing how an MCMC sampler could drive the native buffer
# likelihood for Example 6 without ever calling the R fallback path.

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

# 1. Prep model/data once ------------------------------------------------------
example_idx <- 6
n_trials <- 1500
set.seed(2025)

model_spec <- new_api_examples[[example_idx]]
structure <- build_generator_structure(model_spec)
model_tables <- model_to_tables(model_spec)
sim_data <- simulate_model(model_tables, n_trials = n_trials)
param_table <- build_param_table(structure, n_trials)

native_bundle <- likelihood_build_native_bundle(model_spec = model_spec)
prep <- native_bundle$prep

# Build a buffer-style plan: only row indices, no per-trial data frames.
buffer_plan <- .likelihood_build_trial_plan(
  structure = structure,
  params_df = param_table,
  prep = prep,
  build_override_ptr = FALSE,
  keep_trial_rows = FALSE,
  keep_component_rows = FALSE
)

# Helper the sampler would call each iteration ---------------------------------
evaluate_loglik <- function(params_df) {
  log_likelihood_from_params_buffer(
    structure = structure,
    params_df = params_df,
    data_df = sim_data,
    prep = prep,
    trial_plan = buffer_plan
  )$loglik
}

# 2. Demonstrate a tiny "sampling" loop ----------------------------------------
# Pretend the sampler proposes new lognormal means for TaskA/TaskB, while all
# other columns in param_table come from design matrices.
current_params <- param_table
current_ll <- evaluate_loglik(current_params)
cat(sprintf("Initial log-likelihood: %.3f\n", current_ll))

proposal_scales <- list(acc_taskA = 0.02, acc_taskB = 0.02)
n_iter <- 5

for (iter in seq_len(n_iter)) {
  proposal <- current_params
  for (acc in names(proposal_scales)) {
    mask <- proposal$accumulator_id == acc
    proposal$meanlog[mask] <- proposal$meanlog[mask] + rnorm(1, sd = proposal_scales[[acc]])
  }
  proposal_ll <- evaluate_loglik(proposal)
  cat(sprintf("Iter %02d: proposal log-likelihood %.3f\n", iter, proposal_ll))
  # In a real sampler we would compute an acceptance ratio here; for this demo
  # accept if the proposal improved the likelihood.
  if (proposal_ll > current_ll) {
    current_params <- proposal
    current_ll <- proposal_ll
    cat("  accepted\n")
  } else {
    cat("  rejected\n")
  }
}

cat(sprintf("Final log-likelihood after %d iterations: %.3f\n", n_iter, current_ll))
