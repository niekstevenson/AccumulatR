rm(list = ls())

devtools::load_all(quiet = TRUE)


source("dev/examples/new_API.R")
source("dev/scripts/compare_likelihoods.R")

example_ids <- names(new_api_examples)[c(1:11, 13:16, 18)]
response_summary <- compare_response_suite(new_api_examples, example_ids)
print(response_summary[c("model", "max_abs_diff")])

param_example <- "example_3_stop_na"
param_results <- run_param_table_benchmark(
  model_spec = new_api_examples[[param_example]],
  n_trials = 1000L,
  seed = 2025,
  bench_reps = 3
)

cat("\nLog-likelihood difference (native - R):", param_results$loglik$loglik_diff, "\n")
cat("Native batch used:", param_results$loglik$native_batch_used, "\n")
print(param_results$speed$timings)
