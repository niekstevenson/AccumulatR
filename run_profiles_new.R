rm(list = ls())
source("examples/new_API.R")
source("R/generator_new.R")
source("R/likelihood_new.R")
source("R/profile_plot_new.R")

set.seed(123)
model <- new_api_examples[[5]]

# More trials for curvature
data <- simulate_model(model, n_trials = 5000)


# ---- Simple sanity checks (similar to examples/1_simple_two_response_race.R) ----

# 2) Compute and print simulation counts
cat("\nSimulation counts:\n")
print(table(data$outcome, useNA = "ifany"))
cat("\nSimulation probabilities:\n")
print(table(data$outcome, useNA = "ifany")/nrow(data))


# 3) Analytic probability check via likelihood_new.R helpers
cat("\nAnalytic outcome probabilities:\n")
probs <- observed_response_probabilities(model, include_na = TRUE)
print(round(probs, 6))
cat(sprintf("Sum of probabilities: %.6f\n", sum(probs)))
#
# 4) Ensure likelihood is finite on provided data
ll <- compute_loglik(model, data)
cat(sprintf("\nLog-likelihood on test data: %.4f\n", ll))


profile_res <- profile_likelihood(
  model = model,
  data = data,
  n_points = 10,
  percent = 0.10,
  n_cores = 10
)

print(profile_res)

plot_profile(profile_res)
#
#
#
