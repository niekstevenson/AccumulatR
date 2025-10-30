rm(list = ls())

source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("examples/new_API.R")

# Build the shared-trigger stop model
shared_stop_model <- race_spec() |>
  add_accumulator("go_left", "lognormal",
                  params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("go_right", "lognormal",
                  params = list(meanlog = log(0.32), sdlog = 0.20)) |>
  add_accumulator("stop_fast", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.16), sdlog = 0.15)) |>
  add_accumulator("stop_slow", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.18), sdlog = 0.15)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("STOP", c("stop_fast", "stop_slow")) |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("STOP", "STOP", options = list(map_outcome_to = NA_character_)) |>
  add_group("stop_fast_component", members = "stop_fast", attrs = list(component = "fast")) |>
  add_group("stop_slow_component", members = "stop_slow", attrs = list(component = "slow")) |>
  add_group("shared_stop_trigger",
            members = c("stop_fast", "stop_slow"),
            attrs = list(shared_trigger = list(id = "stop_gate", q = 0.25))) |>
  set_metadata(mixture = list(
    components = list(
      component("fast", weight = 0.25),
      component("slow", weight = 0.75)
    )
  )) |>
  build_model()

# Prepare generator structure and per-trial parameter table
generator_struct <- build_generator_structure(shared_stop_model)
acc_lookup <- generator_struct$accumulators
trial_count <- 300L

base_table <- data.frame(
  trial = 1L,
  accumulator_id = acc_lookup$accumulator_id,
  accumulator = acc_lookup$accumulator_index,
  component = ifelse(acc_lookup$accumulator_id == "stop_fast", "fast",
                     ifelse(acc_lookup$accumulator_id == "stop_slow", "slow", NA_character_)),
  role = acc_lookup$role,
  onset = acc_lookup$onset,
  q = acc_lookup$q,
  stringsAsFactors = FALSE
)

param_names <- unique(unlist(lapply(acc_lookup$params, names)))
if (length(param_names) > 0L) {
  for (param_name in param_names) {
    base_vals <- vapply(acc_lookup$params, function(p) {
      val <- p[[param_name]]
      if (is.null(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
    base_table[[param_name]] <- base_vals
  }
}

param_table <- base_table[rep(seq_len(nrow(base_table)), times = trial_count), , drop = FALSE]
param_table$trial <- rep(seq_len(trial_count), each = nrow(base_table))

theoretical_probs <- observed_response_probabilities_from_params(generator_struct, param_table, include_na = TRUE)

simulated <- simulate_trials_from_params(generator_struct, param_table, seed = 2025)

total_trials <- nrow(simulated)
sim_probs <- c(
  Left = sum(simulated$outcome == "Left", na.rm = TRUE) / total_trials,
  Right = sum(simulated$outcome == "Right", na.rm = TRUE) / total_trials,
  `NA` = sum(is.na(simulated$outcome)) / total_trials
)

cat("Empirical probabilities (generator):\n")
print(sim_probs)
cat("Analytic probabilities (likelihood):\n")
print(theoretical_probs)

diffs <- abs(sim_probs - theoretical_probs[names(sim_probs)])
cat("Max absolute difference:", max(diffs), "\n")

tol <- 0.02
if (any(diffs > tol)) {
  stop(sprintf("Shared-trigger likelihood mismatch: max diff %.4f exceeds tolerance %.4f",
               max(diffs), tol))
}

cat("Shared-trigger generator vs likelihood check passed.\n")
