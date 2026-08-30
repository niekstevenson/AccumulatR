library(pkgload)

load_all(".", quiet = TRUE, helpers = FALSE)
source("dev/examples/benchmark_models.R")

n_trials <- 100L
samples <- 3L
evaluations_per_sample <- 1000L

make_case <- function(model) {
  parameters <- build_param_matrix(
    model$structure, model$pars, n_trials = n_trials
  )
  data <- simulate(
    model$structure, parameters, seed = 123L, keep_component = TRUE
  )
  data <- prepare_data(model$structure, data)
  parameters <- build_param_matrix(
    model$structure, model$pars, trial_df = data
  )
  list(
    context = make_context(model$structure),
    data = data,
    parameters = parameters
  )
}

time_case <- function(case) {
  evaluate <- function() {
    log_likelihood(case$context, case$data, case$parameters)
  }
  invisible(evaluate())
  elapsed <- replicate(samples, system.time({
    for (i in seq_len(evaluations_per_sample)) evaluate()
  })[["elapsed"]])
  c(
    median_ms = 1000 * median(elapsed) / evaluations_per_sample,
    min_ms = 1000 * min(elapsed) / evaluations_per_sample,
    max_ms = 1000 * max(elapsed) / evaluations_per_sample
  )
}

cases <- lapply(benchmark_models, make_case)
timings <- t(vapply(cases, time_case, numeric(3)))
results <- data.frame(
  model = names(cases),
  n_trials = n_trials,
  timings,
  row.names = NULL,
  check.names = FALSE
)

print(results, row.names = FALSE, digits = 4)
dir.create("dev/scripts/scratch_outputs", showWarnings = FALSE)
write.csv(
  results,
  "dev/scripts/scratch_outputs/benchmark_speed.csv",
  row.names = FALSE
)
