rm(list = ls())

source("examples/new_API.R")
source("R/generator_new.R")

# Use the third API example (stop task with NA outcome mapping)
model_spec <- new_api_examples[[3]]
generator_spec <- build_generator_structure(model_spec)

# Static lookup: accumulator metadata that stays constant across trials
accumulator_lookup <- generator_spec$accumulators[, c(
  "accumulator_index", "accumulator_id", "role", "dist", "onset", "q"
)]
accumulator_lookup$params <- generator_spec$accumulators$params

print("Accumulator lookup")
print(accumulator_lookup)

# Trial-level metadata: this would usually come from the experimental design
trial_layout <- data.frame(
  trial = 1:6,
  condition = rep(c("BASE", "STOP"), each = 3),
  stringsAsFactors = FALSE
)

# Build the parameter table: one row per trial × accumulator
n_acc <- nrow(accumulator_lookup)
n_trials <- nrow(trial_layout)
param_table <- data.frame(
  trial = rep(trial_layout$trial, each = n_acc),
  condition = rep(trial_layout$condition, each = n_acc),
  accumulator_id = rep(accumulator_lookup$accumulator_id, times = n_trials),
  accumulator = rep(accumulator_lookup$accumulator_index, times = n_trials),
  component = "__default__",
  role = rep(accumulator_lookup$role, times = n_trials),
  onset = rep(accumulator_lookup$onset, times = n_trials),
  q = rep(accumulator_lookup$q, times = n_trials),
  stringsAsFactors = FALSE
)

# Copy the distribution parameters into explicit columns
param_names <- unique(unlist(lapply(accumulator_lookup$params, names)))
if (length(param_names) > 0L) {
  for (param_name in param_names) {
    base_vals <- vapply(accumulator_lookup$params, function(p) {
      val <- p[[param_name]]
      if (is.null(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
    trial_adjustment <- if (identical(param_name, "meanlog")) {
      rep(seq_len(n_trials) - 1, each = n_acc) * 0.05
    } else {
      rep(0, n_acc * n_trials)
    }
    param_table[[param_name]] <- rep(base_vals, times = n_trials) + trial_adjustment
  }
}

print("Trial parameter table (first 12 rows)")
print(head(param_table, 12))

# Simulate data directly from the parameter table
set.seed(42)
simulated_trials <- simulate_trials_from_params(generator_spec, param_table, seed = 42)

simulated_trials <- merge(
  trial_layout,
  simulated_trials,
  by = "trial",
  sort = TRUE
)

print("Simulated trials")
print(simulated_trials)

# --- Mixture stop example ------------------------------------------------------
