rm(list = ls())

base::source("R/super_large_likelihood.R")
base::source("R/generator_new.R")
base::source("examples/stim_selective_versions.R")

model_spec <- stim_selective_versions[[1]]
model_tables <- model_to_tables(model_spec)

analytic <- observed_response_probabilities(model_tables, include_na = TRUE)
if (!isTRUE(all.equal(sum(analytic), 1.0, tolerance = 1e-8))) {
  stop(sprintf("Analytic probabilities sum to %.12f", sum(analytic)))
}

go2_prob <- analytic[["GO2"]] %||% 0
if (!is.finite(go2_prob) || go2_prob <= 0) {
  stop("GO2 probability should be positive for guarded stimulus-selective model")
}

set.seed(2025)
sim_df <- simulate_model(model_tables, n_trials = 4000)
sim_tab <- prop.table(table(sim_df$outcome, useNA = "ifany"))
sim_go2 <- sim_tab[["GO2"]] %||% 0

if (!is.finite(sim_go2)) {
  stop("Simulated GO2 probability is not finite")
}

if (abs(go2_prob - sim_go2) > 0.05) {
  stop(sprintf(
    "Analytic GO2 probability %.4f differs from simulation %.4f",
    go2_prob, sim_go2
  ))
}
