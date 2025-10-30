rm(list = ls())
options(uuber.param_cache_across_trials = TRUE)
source("examples/new_API.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")

# Model: example_3_stop_na
model_spec <- new_api_examples[[3]]

set.seed(2025)
n_trials <- 1000L

# Generate synthetic data via generator
structure <- build_generator_structure(model_spec)

acc_lookup <- structure$accumulators
base_table <- data.frame(
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
    base_table[[param_name]] <- rep(base_vals, times = n_trials)
  }
}

set.seed(2025)
data_sim <- simulate_trials_from_params(structure, base_table, seed = 2025)
data_df <- data.frame(
  trial = data_sim$trial,
  outcome = data_sim$outcome,
  rt = data_sim$rt,
  stringsAsFactors = FALSE
)

# Old interface profile
model_tables <- model_to_tables(model_spec)

library(profvis)
profvis({
  ll_old <- compute_loglik(model_tables, data_df)
}, prof_output = "ll_old_table.out")

prof_old <- summaryRprof("ll_old_table.out")
write.csv(prof_old$by.total, file = "ll_old_table_total.csv")
write.csv(prof_old$by.self, file = "ll_old_table_self.csv")

# New interface using parameter table
profvis({
  res_new <- log_likelihood_from_params(structure, base_table, data_df)
}, prof_output = "ll_new_table.out")

prof_new <- summaryRprof("ll_new_table.out")
write.csv(prof_new$by.total, file = "ll_new_table_total.csv")
write.csv(prof_new$by.self, file = "ll_new_table_self.csv")

# Compare key metrics
cat("Old interface total time (seconds):", prof_old$sampling.time, "\n")
cat("New interface total time (seconds):", prof_new$sampling.time, "\n")

library(dplyr)

load_profile <- function(path_total, path_self) {
  total_df <- read.csv(path_total, stringsAsFactors = FALSE)
  self_df <- read.csv(path_self, stringsAsFactors = FALSE)
  total_df <- rename(total_df, total_time = total.time, total_pct = total.pct, fn = X)
  self_df <- rename(self_df, self_time = self.time, self_pct = self.pct, fn = X)
  full_join(total_df, self_df, by = "fn")
}

old_profile <- load_profile("ll_old_table_total.csv", "ll_old_table_self.csv")
new_profile <- load_profile("ll_new_table_total.csv", "ll_new_table_self.csv")

compare_profile <- full_join(old_profile, new_profile, by = "fn", suffix = c("_old", "_new")) %>%
  mutate(
    total_time_diff = total_time_new - total_time_old,
    total_pct_diff = total_pct_new - total_pct_old,
    self_time_diff = self_time_new - self_time_old,
    self_pct_diff = self_pct_new - self_pct_old
  ) %>%
  arrange(desc(abs(total_pct_diff)))

write.csv(compare_profile, file = "ll_table_profile_comparison.csv", row.names = FALSE)

cat("Profile comparison written to ll_table_profile_comparison.csv\n")
