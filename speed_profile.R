rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")

set.seed(123456)
model_spec <- new_api_examples[[2]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 1000)
# ll <- compute_loglik(model_tables, data)
# sum(ll)
# New likelihood
library(profvis)
profvis({
  ll <- compute_loglik(model_tables, data)
}, prof_output = "ll_new.out"

)
df_new <- summaryRprof("ll_new.out")
write.csv(df_new$by.total, file = "ll_new_total.csv")
write.csv(df_new$by.self, file = "ll_new_self.csv")

# Old likelihood
source("R/likelihood_old.R")
library(profvis)
profvis({
  ll <- compute_loglik(model_tables, data)
}, prof_output = "ll_old.out"

)
df_old <- summaryRprof("ll_old.out")
write.csv(df_old$by.total, file = "ll_old_total.csv")
write.csv(df_old$by.self, file = "ll_old_self.csv")


