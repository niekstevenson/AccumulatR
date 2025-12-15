rm(list = ls())

# Load the package (devtools if available, otherwise the installed build)
if (!"AccumulatR" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(AccumulatR)
  }
}

# Bring in example definitions and profiling helpers
source("dev/examples/stim_selective_versions.R")
source("dev/R_extra/profile_plot_new.R")
source("dev/R_extra/processing_tree.R")



# Simple helper to force native likelihood code paths
with_native_only <- function(expr) {
  old_opts <- options(
    uuber.use.native.node.eval = TRUE,
    uuber.use.native.param.rows = TRUE
  )
  on.exit(options(old_opts), add = TRUE)
  force(expr)
}

# 1. Which examples to run
example_ids <- names(new_api_examples)[23]

# Optional overrides via env vars
if (nzchar(Sys.getenv("UUBER_EXAMPLES"))) {
  example_ids <- strsplit(Sys.getenv("UUBER_EXAMPLES"), ",", fixed = TRUE)[[1]]
  example_ids <- trimws(example_ids)
}
n_trials <- as.integer(Sys.getenv("UUBER_N_TRIALS", "3000"))
seed <- as.integer(Sys.getenv("UUBER_SEED", "2025"))
profile_percent <- as.numeric(Sys.getenv("UUBER_PROFILE_PERCENT", "0.10"))
profile_points <- as.integer(Sys.getenv("UUBER_PROFILE_POINTS", "10"))
plot_profiles <- identical(Sys.getenv("UUBER_PROFILE_PLOT", "0"), "1")

for (example_id in example_ids) {
  cat(sprintf("\n=== Example %s ===\n", example_id))
  model_spec <- new_api_examples[[example_id]]
  core_params <- new_api_example_params[[example_id]]
  # processing_math_tree(model_spec, "A")
  
  # 2. Parameter table for this example
  params_df <- build_params_df(model_spec, core_params, n_trials = n_trials)

  # 3. Simulate data and compute observed probabilities
  sim <- simulate(model_spec, params_df, seed = seed)
  data_df <- data.frame(
    trial = sim$trial,
    outcome = sim$outcome,
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
  obs_counts <- table(data_df$outcome, useNA = "ifany")
  obs_probs <- prop.table(obs_counts)
  cat("Observed outcome probabilities (simulated data):\n")
  print(round(obs_probs, 6))

  # 4. Analytic response probabilities (native)
  single_params <- build_params_df(model_spec, core_params, n_trials = 1L)
  analytic_probs <- with_native_only({
    response_probabilities(model_spec, single_params, include_na = TRUE)
  })
  cat("Analytic outcome probabilities (native):\n")
  print(round(analytic_probs, 6))
  cat(sprintf("Sum: %.6f\n", sum(analytic_probs)))

  # 5. Native profiling over parameter grids and repeat log-likelihoods
  ctx <- build_likelihood_context(
    structure = model_spec,
    params_df = params_df,
    data_df = data_df
  )
  profile_res <- profile_likelihood(
    structure = model_spec,
    model_spec = model_spec,
    base_params = core_params,
    data = data_df,
    percent = profile_percent,
    n_points = profile_points,
    n_cores = 10
  )
  print(utils::head(profile_res))
  print(plot_profile(profile_res))
}
