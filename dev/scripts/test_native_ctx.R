#!/usr/bin/env Rscript

# Scratch test for ensure_native_ctx() and pointer rebuilds
suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(AccumulatR)
  }
})

# Build a minimal model
spec <- race_spec()
spec <- add_accumulator(spec, "go1", "lognormal", params = list(meanlog = 0, sdlog = 0.1))
spec <- add_accumulator(spec, "go2", "lognormal", params = list(meanlog = 0.1, sdlog = 0.15))
spec <- add_outcome(spec, "go1_win", "go1")
spec <- add_outcome(spec, "go2_win", "go2")
model <- finalize_model(spec)

# Params and data
params <- c(
  go1.meanlog = 0, go1.sdlog = 0.1, go1.q = 0.6, go1.t0 = 0.1,
  go2.meanlog = 0.1, go2.sdlog = 0.15, go2.q = 0.4, go2.t0 = 0.1
)
param_mat <- build_param_matrix(model, params, n_trials = 10)
data <- simulate(model, param_mat, seed = 123, keep_detail = FALSE)

# Build context and compute likelihood
ctx <- build_likelihood_context(model, data)
ll1 <- log_likelihood(ctx, param_mat)
cat("Initial log-likelihood:", ll1, "\n")

# Create a transportable context list (as described)
context <- list(
  native_ctx = ctx$native_ctx,
  data_df = ctx$data_df,
  rel_tol = ctx$rel_tol,
  abs_tol = ctx$abs_tol,
  max_depth = ctx$max_depth,
  model_spec = spec
)
class(context) <- "likelihood_context"

# Simulate a mangled pointer and rebuild
context$native_ctx <- new("externalptr") # intentionally invalid
context <- ensure_native_ctx(context, model_spec = context$model_spec)
ll2 <- log_likelihood(context, param_mat)
cat("Rebuilt log-likelihood:", ll2, "\n")

if (!isTRUE(all.equal(ll1, ll2))) {
  warning("Rebuilt likelihood differs from initial value")
} else {
  cat("Pointer rebuild succeeded; likelihood matches.\n")
}
