rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/likelihood_new.R")
source("R/profile_plot_new.R")
source("R/processing_tree.R")

set.seed(123456)
model_spec <- new_api_examples[[6]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

# processing_math_tree(model_tables, "A")


# More trials for curvature
data <- simulate_model(model_tables, n_trials = 1000)

# ---- Simple sanity checks (similar to examples/1_simple_two_response_race.R) ----

# 2) Compute and print simulation counts
cat("\nSimulation counts:\n")
print(table(data$outcome, useNA = "ifany"))
# print(table(data$component, data$outcome, useNA = "ifany"))
cat("\nSimulation probabilities:\n")
print(table(data$outcome, useNA = "ifany")/nrow(data))



#
# 3) Analytic probability check via likelihood helpers
# cat("\nAnalytic outcome probabilities:\n")
# probs <- observed_response_probabilities(model_tables, include_na = TRUE)
# print(round(probs, 6))
# cat(sprintf("Sum of probabilities: %.6f\n", sum(probs)))
#
# # #
# 4) Ensure likelihood is finite on provided data
library(profvis)
profvis(
  ll <- compute_loglik(model_tables, data)
)


cat(sprintf("\nLog-likelihood on test data: %.4f\n", ll))
#
# profile_res <- profile_likelihood(
#   model = model_tables,
#   data = data,
#   n_points = 10,
#   percent = 0.10,
#   n_cores = 10
# )
#
# print(profile_res)
#
# plot_profile(profile_res)
