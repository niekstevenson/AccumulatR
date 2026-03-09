#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_PROFILE_SEED", "20260308")))
library(AccumulatR)

run_sec <- as.numeric(Sys.getenv("ACCUMULATR_PROFILE_RUN_SEC", "45"))
if (!is.finite(run_sec) || run_sec <= 0) {
  run_sec <- 45
}
n_trials <- as.integer(Sys.getenv("ACCUMULATR_PROFILE_TRIALS", "250"))
if (!is.finite(n_trials) || n_trials < 50L) {
  n_trials <- 250L
}

example_2_stop_mixture <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("stop", "exgauss", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_outcome("R1", inhibit("go1", by = "stop")) |>
    add_outcome("R2", all_of("go2", "stop")) |>
    add_component("go_only", members = c("go1"),
                  attrs = list(component = "go_only"), weight = 0.5) |>
    add_component("go_stop", members = c("go1", "stop", "go2"),
                  attrs = list(component = "go_stop"), weight = 0.5) |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()
}

params_example_2_stop_mixture <- c(
  go1.meanlog = log(0.35),
  go1.sdlog = 0.2,
  stop.mu = 0.1,
  stop.sigma = 0.04,
  stop.tau = 0.1,
  go2.meanlog = log(0.60),
  go2.sdlog = 0.18
)

make_case <- function(spec_builder, params) {
  spec_obj <- spec_builder()
  structure <- finalize_model(spec_obj)
  params_df <- build_param_matrix(spec_obj, params, n_trials = n_trials)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  ctx <- build_likelihood_context(structure, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  params_df_slim <- build_param_matrix(
    spec_obj,
    params,
    n_trials = max(data_df$trial),
    layout = ctx$param_layout
  )
  invisible(log_likelihood(ctx, data_prepped, params_df_slim))
  list(ctx = ctx, data_prepped = data_prepped, params_df = params_df_slim)
}

cases <- list(
  example_2_stop_mixture = make_case(
    example_2_stop_mixture,
    params_example_2_stop_mixture
  )
)

case_filter <- trimws(Sys.getenv("ACCUMULATR_PROFILE_CASES", ""))
if (nzchar(case_filter)) {
  selected <- strsplit(case_filter, ",", fixed = TRUE)[[1]]
  selected <- trimws(selected)
  selected <- selected[nzchar(selected)]
  selected <- intersect(selected, names(cases))
  if (length(selected) == 0L) {
    stop("ACCUMULATR_PROFILE_CASES did not match any known example cases")
  }
  cases <- cases[selected]
}

invisible(lapply(cases, function(x) {
  log_likelihood(x$ctx, x$data_prepped, x$params_df)
}))

start <- Sys.time()
checksum <- 0.0
iters <- integer(length(cases))
names(iters) <- names(cases)

while (as.numeric(difftime(Sys.time(), start, units = "secs")) < run_sec) {
  for (nm in names(cases)) {
    val <- log_likelihood(cases[[nm]]$ctx, cases[[nm]]$data_prepped,
                          cases[[nm]]$params_df)
    checksum <- checksum + sum(val)
    iters[[nm]] <- iters[[nm]] + 1L
  }
}

cat("profile_workload_examples complete\n")
print(iters)
cat("checksum:", format(checksum, digits = 16), "\n")
