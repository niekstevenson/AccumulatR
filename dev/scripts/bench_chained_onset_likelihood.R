#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AccumulatR)
})

run_bench <- function(label, expr, n_rep = 8L) {
  timings <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    gc(FALSE)
    timings[[i]] <- system.time(force(expr))[["elapsed"]]
  }
  c(
    label = label,
    median_sec = stats::median(timings),
    min_sec = min(timings),
    max_sec = max(timings)
  )
}

set.seed(123)
n_trials <- 2000L

# Baseline: independent single-outcome likelihood path.
spec_base <- race_spec() |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal") |>
  add_outcome("A", "a") |>
  add_outcome("B", "b")
structure_base <- finalize_model(spec_base)
params_base <- c(
  a.meanlog = log(0.30), a.sdlog = 0.20, a.q = 0.00, a.t0 = 0.00,
  b.meanlog = log(0.40), b.sdlog = 0.20, b.q = 0.00, b.t0 = 0.00
)
params_df_base <- build_param_matrix(spec_base, params_base, n_trials = n_trials)
sim_base <- simulate(structure_base, params_df_base, seed = 11)
ctx_base <- build_likelihood_context(structure_base, sim_base)

# Chained single-outcome path.
spec_chain <- race_spec() |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal", onset = after("a")) |>
  add_outcome("B", "b")
structure_chain <- finalize_model(spec_chain)
params_chain <- c(
  a.meanlog = log(0.26), a.sdlog = 0.18, a.q = 0.00, a.t0 = 0.00,
  b.meanlog = log(0.28), b.sdlog = 0.18, b.q = 0.00, b.t0 = 0.00
)
params_df_chain <- build_param_matrix(spec_chain, params_chain, n_trials = n_trials)
sim_chain <- simulate(structure_chain, params_df_chain, seed = 22)
ctx_chain <- build_likelihood_context(structure_chain, sim_chain)

# Chained ranked path (k = 2).
spec_chain_ranked <- race_spec(n_outcomes = 2L) |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal", onset = after("a")) |>
  add_outcome("A", "a") |>
  add_outcome("B", "b")
structure_chain_ranked <- finalize_model(spec_chain_ranked)
params_chain_ranked <- c(
  a.meanlog = log(0.30), a.sdlog = 0.16, a.q = 0.00, a.t0 = 0.00,
  b.meanlog = log(0.22), b.sdlog = 0.16, b.q = 0.00, b.t0 = 0.00
)
params_df_chain_ranked <- build_param_matrix(
  spec_chain_ranked,
  params_chain_ranked,
  n_trials = n_trials
)
sim_chain_ranked <- simulate(structure_chain_ranked, params_df_chain_ranked, seed = 33)
ctx_chain_ranked <- build_likelihood_context(structure_chain_ranked, sim_chain_ranked)

bench <- rbind(
  run_bench(
    "single_outcome_baseline",
    log_likelihood(ctx_base, params_df_base)
  ),
  run_bench(
    "single_outcome_chained",
    log_likelihood(ctx_chain, params_df_chain)
  ),
  run_bench(
    "ranked_chained_k2",
    log_likelihood(ctx_chain_ranked, params_df_chain_ranked)
  )
)

bench <- as.data.frame(bench, stringsAsFactors = FALSE)
bench$median_sec <- as.numeric(bench$median_sec)
bench$min_sec <- as.numeric(bench$min_sec)
bench$max_sec <- as.numeric(bench$max_sec)
print(bench, row.names = FALSE)

baseline <- bench$median_sec[bench$label == "single_outcome_baseline"]
chained <- bench$median_sec[bench$label == "single_outcome_chained"]
if (length(baseline) == 1L && length(chained) == 1L) {
  ratio <- chained / baseline
  cat(sprintf("\nSingle-outcome chained/baseline ratio: %.3f\n", ratio))
}
