#!/usr/bin/env Rscript

# Minimal example showing how an MCMC sampler could drive the native buffer
# likelihood for Example 6 without ever calling the R fallback path.

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  source("examples/new_API.R")
  source("R/dist.R")
  source("R/utils.R")
  source("R/pool_math.R")
  source("R/model_tables.R")
  source("R/generator_new.R")
  source("R/likelihood_cache.R")
  source("R/likelihood_common.R")
  source("R/likelihood_prep.R")
  source("R/likelihood_primitives.R")
  source("R/likelihood_kernels.R")
  source("R/likelihood_integrate.R")
  source("R/likelihood_param_interface.R")
  source("R/super_large_likelihood.R")
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
example_idx <- 3
n_trials <- 1500
set.seed(2025)

model_spec <- new_api_examples[[example_idx]]
structure <- build_generator_structure(model_spec)
model_tables <- model_to_tables(model_spec)
sim_data <- simulate_model(model_tables, n_trials = n_trials)
param_table <- build_param_table(structure, n_trials)

native_bundle <- likelihood_build_native_bundle(model_spec = model_spec)
prep <- native_bundle$prep
default_deadline <- prep$default_deadline %||% Inf

# Build a buffer-style plan: only row indices, no per-trial data frames.
buffer_plan <- .likelihood_build_trial_plan(
  structure = structure,
  params_df = param_table,
  prep = prep,
  build_override_ptr = FALSE,
  keep_trial_rows = FALSE,
  keep_component_rows = FALSE
)

buffer_bundle <- build_buffer_trial_entries(
  structure = structure,
  params_df = param_table,
  data_df = sim_data,
  prep = prep,
  trial_plan = buffer_plan
)
trial_entries <- buffer_bundle$entries

Rcpp::sourceCpp("scripts/demo_buffer_batch_cpp.cpp")

# Build a few example parameter frames to evaluate in one C++ call ----------
param_proposals <- list(
  base = param_table,
  taskA_plus = {
    df <- param_table
    mask <- df$accumulator_id == "go_left"
    df$meanlog[mask] <- df$meanlog[mask] + 0.05
    df
  },
  taskB_minus = {
    df <- param_table
    mask <- df$accumulator_id == "go_right"
    df$meanlog[mask] <- df$meanlog[mask] - 0.04
    df
  }
)

lik_values <- demo_buffer_batch_loglik_cpp(
  ctx_ptr = native_bundle$context,
  structure = structure,
  trial_entries = trial_entries,
  params_list = param_proposals,
  component_weights = buffer_bundle$component_weights,
  default_deadline = buffer_bundle$default_deadline,
  rel_tol = buffer_bundle$rel_tol,
  abs_tol = buffer_bundle$abs_tol,
  max_depth = buffer_bundle$max_depth
)

cat("Log-likelihoods (base, go_left_plus, go_right_minus):\n")
print(lik_values)
