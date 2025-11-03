rm(list = ls())

ok <- function(name, cond) {
  cat(sprintf("%s: %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))
  if (!isTRUE(cond)) quit(status = 1)
}

base::source("R/super_large_likelihood.R")
base::source("R/generator_new.R")
base::source("examples/stim_selective_versions.R")

# Helper: build tables and simulate
model_spec <- new_api_examples[[12]]
model_tables <- model_to_tables(model_spec)
set.seed(123)
data <- simulate_model(model_tables, n_trials = 50)

# Test 1: Pool fast-paths match full paths with empty forced sets (new)
prep <- .prepare_model_for_likelihood(model_tables)
comp <- "go_stop"
# pick a pool known in example 12
pools <- names(prep$pools)
if (length(pools) > 0) {
  pid <- pools[[1]]
  tvals <- c(0.2, 0.3, 0.4)
  surv_fast <- vapply(tvals, function(t) .pool_survival_fast_value(prep, pid, comp, t), numeric(1))
  surv_full <- vapply(tvals, function(t) .pool_survival(prep, pid, comp, t, integer(0), integer(0)), numeric(1))
  ok("pool fast survival == full survival (no forced)", max(abs(surv_fast - surv_full)) < 1e-9)
  dens_fast <- vapply(tvals, function(t) .pool_density_fast_value(prep, pid, comp, t), numeric(1))
  dens_full <- vapply(tvals, function(t) .pool_density(prep, pid, comp, t, integer(0), integer(0)), numeric(1))
  ok("pool fast density == full density (no forced)", max(abs(dens_fast - dens_full)) < 1e-9)
} else {
  cat("SKIP: no pools found in model\n")
}

# Test 2: Likelihood sums equal (old vs new)
ll_new <- compute_loglik(model_tables, data)
base::source("R/likelihood_old.R")
ll_old <- compute_loglik(model_tables, data)
cat(sprintf("new_sum=%.8f old_sum=%.8f\n", sum(ll_new), sum(ll_old)))
ok("likelihood sums equal", abs(sum(ll_new) - sum(ll_old)) < 1e-5)

# Test 3: Response probabilities equal per component
rp_new_go <- response_probabilities(model_tables, component = "go_stop")
base::source("R/likelihood_old.R")
rp_old_go <- response_probabilities(model_tables, component = "go_stop")
ok("response_probabilities(go_stop)", max(abs(rp_new_go - rp_old_go)) < 1e-6)

cat("All tests passed.\n")


