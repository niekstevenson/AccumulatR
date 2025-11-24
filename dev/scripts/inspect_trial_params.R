#!/usr/bin/env Rscript
# Scratch helper to inspect the native trial entry inputs and the param table
# schema to understand what must change per proposal.

rm(list = ls())

if (!"AccumulatR" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(AccumulatR)
  }
}

source("dev/examples/new_API.R")

example_id <- "example_1_simple"
n_trials <- 5L
seed <- 2025L

model_spec <- new_api_examples[[example_id]]
core_params <- new_api_example_params[[example_id]]
structure <- build_generator_structure(model_spec)
params_df <- build_params_df(model_spec, core_params, n_trials = n_trials)
data_df <- simulate_data_from_params_table(structure, params_df, seed = seed)

ctx <- build_likelihood_context(
  structure = structure,
  params_df = params_df,
  data_df = data_df
)

entries <- ctx$plan$.native_cache$entries
cat(sprintf("Total entries cached: %d\n", length(entries)))

# Summarize parameter table schema
num_cols <- names(Filter(is.numeric, params_df))
other_cols <- setdiff(names(params_df), num_cols)
cat("Numeric columns in params_df:\n")
print(num_cols)
cat("Non-numeric columns in params_df:\n")
print(other_cols)

# Per-entry trial_rows summary (names/types)
for (i in seq_along(entries)) {
  entry <- entries[[i]]
  trial_rows <- entry$trial_rows
  cat(sprintf("\nEntry %d trial_rows (n=%d):\n", i, nrow(trial_rows)))
  cat("  Columns:\n")
  for (nm in names(trial_rows)) {
    cls <- paste(class(trial_rows[[nm]]), collapse = ",")
    cat(sprintf("    %s : %s\n", nm, cls))
  }
  comp_plan <- entry$component_plan
  comps <- comp_plan$components
  cat(sprintf("  Components: %s\n", paste(comps, collapse = ", ")))
  if (!is.null(entry$trial_params_ptr)) {
    cat("  trial_params_ptr: externalptr present\n")
  } else {
    cat("  trial_params_ptr: <NULL>\n")
  }
}

# Show a single trial_rows data frame for reference
cat("\nFirst trial_rows snapshot:\n")
print(entries[[1]]$trial_rows)

# Show the first component plan
cat("\nFirst component plan:\n")
print(entries[[1]]$component_plan)
