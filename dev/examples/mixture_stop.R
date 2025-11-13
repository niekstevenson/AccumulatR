rm(list = ls())

source("examples/new_API.R")
source("R/generator_new.R")

mixture_spec <- race_spec() |>
  add_accumulator("go_left", "lognormal",
                  params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("go_right", "lognormal",
                  params = list(meanlog = log(0.32), sdlog = 0.20)) |>
  add_accumulator("stop_fast", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.18), sdlog = 0.15)) |>
  add_accumulator("stop_slow", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.22), sdlog = 0.15)) |>
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
            attrs = list(shared_trigger = list(id = "stop_gate", q = 0.05))) |>
  set_metadata(mixture = list(
    components = list(
      component("fast", weight = 0.25),
      component("slow", weight = 0.75)
    )
  )) |>
  build_model()

mixture_structure <- build_generator_structure(mixture_spec)

mixture_acc_lookup <- mixture_structure$accumulators[, c(
  "accumulator_index", "accumulator_id", "role", "dist", "onset", "q", "shared_trigger_id"
)]
mixture_acc_lookup$params <- mixture_structure$accumulators$params
mixture_acc_lookup$components <- mixture_structure$accumulators$components

print("Mixture accumulator lookup")
print(mixture_acc_lookup)

mixture_trial_layout <- data.frame(
  trial = 1:20,
  condition = rep(c("BASE", "STOP"), each = 10),
  stringsAsFactors = FALSE
)

n_mix_acc <- nrow(mixture_acc_lookup)
n_mix_trials <- nrow(mixture_trial_layout)

mixture_param_table <- data.frame(
  trial = rep(mixture_trial_layout$trial, each = n_mix_acc),
  condition = rep(mixture_trial_layout$condition, each = n_mix_acc),
  accumulator_id = rep(mixture_acc_lookup$accumulator_id, times = n_mix_trials),
  accumulator = rep(mixture_acc_lookup$accumulator_index, times = n_mix_trials),
  component = rep(NA_character_, n_mix_acc * n_mix_trials),
  role = rep(mixture_acc_lookup$role, times = n_mix_trials),
  onset = rep(mixture_acc_lookup$onset, times = n_mix_trials),
  q = NA_real_,
  stringsAsFactors = FALSE
)

mixture_param_table$component[mixture_param_table$accumulator_id == "stop_fast"] <- "fast"
mixture_param_table$component[mixture_param_table$accumulator_id == "stop_slow"] <- "slow"

mixture_param_names <- unique(unlist(lapply(mixture_acc_lookup$params, names)))
if (length(mixture_param_names) > 0L) {
  for (param_name in mixture_param_names) {
    base_vals <- vapply(mixture_acc_lookup$params, function(p) {
      val <- p[[param_name]]
      if (is.null(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
    condition_shift <- if (identical(param_name, "meanlog")) {
      rep(ifelse(mixture_trial_layout$condition == "STOP", 0.05, 0), each = n_mix_acc)
    } else {
      rep(0, n_mix_acc * n_mix_trials)
    }
    mixture_param_table[[param_name]] <- rep(base_vals, times = n_mix_trials) + condition_shift
  }
}

# Demonstrate condition-specific trigger failure overrides (shared across members)
stop_rows <- mixture_param_table$accumulator_id %in% c("stop_fast", "stop_slow")
mixture_param_table$q[stop_rows & mixture_param_table$condition == "BASE"] <- 0.05
mixture_param_table$q[stop_rows & mixture_param_table$condition == "STOP"] <- 0.10

print("Mixture parameter table (first 16 rows)")
print(head(mixture_param_table, 16))

set.seed(123)
mixture_sim <- simulate_trials_from_params(mixture_structure, mixture_param_table, seed = 123)
mixture_sim <- merge(
  mixture_trial_layout,
  mixture_sim,
  by = "trial",
  sort = TRUE
)

print("Mixture simulation summary")
print(mixture_sim)
print(table(component = mixture_sim$component))
