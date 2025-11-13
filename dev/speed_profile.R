rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
source("R/likelihood_old.R")

set.seed(123456)

# Example 1 ---------------------------------------------------------------


model_spec <- new_api_examples[[1]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 1000)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_total1.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_self1.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_total1.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_self1.csv")


# Example 2 ---------------------------------------------------------------
rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
set.seed(123456)

model_spec <- new_api_examples[[2]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 500)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_total2.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_self2.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_total2.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_self2.csv")

# Example 3 ---------------------------------------------------------------
rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
set.seed(123456)

model_spec <- new_api_examples[[3]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 500)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_total3.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_self3.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_total3.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_self3.csv")




# Example 6 ---------------------------------------------------------------
rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
set.seed(123456)

model_spec <- new_api_examples[[6]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 500)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_total6.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_self6.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_total6.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_self6.csv")



# Example 12 --------------------------------------------------------------
rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
set.seed(123456)

model_spec <- new_api_examples[[12]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 50)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_total12.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_self12.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_total12.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_self12.csv")

# Example SS1 ---------------------------------------------------------------
rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/super_large_likelihood.R")
source("R/likelihood_old.R")
source("R/likelihood_old.R")

set.seed(123456)

model_spec <- stim_selective_versions[[1]]

# Translate model specification to table representation
model_tables <- model_to_tables(model_spec)

data <- simulate_model(model_tables, n_trials = 100)
ll <- compute_loglik_old(model_tables, data)
sum(ll)
library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_new.out"

)
df_new <- summaryRprof("profiles/ll_new.out")
write.csv(df_new$by.total, file = "profiles/ll_new_totalSS1.csv")
write.csv(df_new$by.self, file = "profiles/ll_new_selfSS1.csv")

# Old likelihood
source("R/likelihood_old.R")
ll <- compute_loglik_old(model_tables, data)
sum(ll)

library(profvis)
profvis({
  ll <- compute_loglik_old(model_tables, data)
}, prof_output = "profiles/ll_old.out"

)
df_old <- summaryRprof("profiles/ll_old.out")
write.csv(df_old$by.total, file = "profiles/ll_old_totalSS1.csv")
write.csv(df_old$by.self, file = "profiles/ll_old_selfSS1.csv")


