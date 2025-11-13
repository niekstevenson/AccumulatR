rm(list = ls())

source("examples/new_API.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/likelihood_new.R")

# Choose an example specification
model_spec <- new_api_examples[[5]]

# Translate specification to table representation
model_tables <- model_to_tables(model_spec)

cat("struct_nodes:\n")
print(head(model_tables$struct_nodes))
cat("\nparam_values:\n")
print(head(model_tables$param_values))
cat("\ncomponents:\n")
print(model_tables$components)

# Simulate a small data set via tables
set.seed(123)
sim_data <- simulate_model(model_tables, n_trials = 1000)

cat("\nSimulation outcome proportions:\n")
print(prop.table(table(sim_data$outcome)))

# Compute likelihood directly from tables
ll_value <- compute_loglik(model_tables, sim_data)
cat(sprintf("\nLog-likelihood: %.4f\n", ll_value))
