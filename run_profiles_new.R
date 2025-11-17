rm(list = ls())
devtools::load_all()
source("examples/new_API.R")
set.seed(123456)
example_id <- 1
model_spec <- new_api_examples[[example_id]]
core_params <- new_api_example_params[[example_id]]
structure <- build_generator_structure(model_spec)
n_trials <- 2000L
params_df <- build_params_df(model_spec, core_params, n_trials = n_trials)

# Simulate data under the parameter table
data <- simulate_trials_from_params(structure, params_df, seed = 123456)

cat("\nSimulation counts:\n")
print(table(data$outcome, useNA = "ifany"))
cat("\nSimulation probabilities:\n")
print(round(prop.table(table(data$outcome, useNA = "ifany")), 6))

# Analytic probabilities from the parameter table
single_params <- build_params_df(model_spec, core_params, n_trials = 1L)
probs <- observed_response_probabilities_from_params(structure, single_params, include_na = TRUE)
cat("\nAnalytic outcome probabilities (single-trial params):\n")
print(round(probs, 6))
cat(sprintf("Sum of probabilities: %.6f\n", sum(probs)))

# Log-likelihood check on simulated data
ll_res <- log_likelihood_from_params(structure, params_df, data)
cat(sprintf("\nLog-likelihood on test data: %.4f\n", ll_res$loglik))

source("dev/R_extra/profile_plot_new.R")
# Profile likelihood using the new parameter-table pipeline
profile_res <- profile_likelihood(
  structure = structure,
  model_spec = model_spec,
  base_params = core_params,
  data = data,
  n_points = 10,
  percent = 0.10,
  n_cores = 10
)

print(profile_res)
plot_profile(profile_res)
