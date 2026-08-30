library(AccumulatR)
source("dev/examples/benchmark_models.R")

case_labels <- c(
  "example_1_simple",
  "example_5_timeout_guess",
  "example_21_simple_q",
  "example_22_shared_q",
  "example_2_stop_mixture",
  "example_7_mixture",
  "example_10_exclusion",
  "example_23_ranked_chain",
  "stop_change_shared_trigger"
)
n_trials <- 50L
duration <- 22

make_case <- function(model) {
  parameters <- build_param_matrix(
    model$structure, model$pars, n_trials = n_trials
  )
  data <- simulate(
    model$structure, parameters, seed = 123L, keep_component = TRUE
  )
  data <- prepare_data(model$structure, data)
  list(
    context = make_context(model$structure),
    data = data,
    parameters = parameters
  )
}

cases <- lapply(benchmark_models[case_labels], make_case)
counts <- setNames(integer(length(cases)), names(cases))
marker <- Sys.getenv("ACCUMULATR_PROFILE_START_FILE")
if (nzchar(marker)) writeLines("ready", marker)
deadline <- proc.time()[["elapsed"]] + duration

repeat {
  for (name in names(cases)) {
    case <- cases[[name]]
    log_likelihood(case$context, case$data, case$parameters)
    counts[[name]] <- counts[[name]] + 1L
  }
  if (proc.time()[["elapsed"]] >= deadline) break
}

print(counts)
