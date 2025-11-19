rm(list = ls())
devtools::load_all()


if (!"AccumulatR" %in% loadedNamespaces()) {
  skip_devtools <- identical(Sys.getenv("UUBER_SKIP_DEVTOOLS"), "1")
  if (!skip_devtools && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(AccumulatR)
  }
}

source("dev/examples/new_API.R")
source("dev/scripts/compare_likelihoods.R")

# example_ids <- names(new_api_examples)[c(1:11, 13:16, 18)]
# example_ids <- names(new_api_examples)[c(5, 15)]
example_ids <- names(new_api_examples)[1]

summary_table <- summarize_example_models(
  model_list = new_api_examples,
  params_list = new_api_example_params,
  example_ids = example_ids,
  n_trials = 1000L,
  seed = 2025
)

print(summary_table)
