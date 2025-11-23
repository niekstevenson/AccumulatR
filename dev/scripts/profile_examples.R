rm(list = ls())

repo_root <- Sys.getenv("UUBER_REPO_ROOT", unset = NA_character_)
if (!is.na(repo_root) && nzchar(repo_root)) {
  setwd(repo_root)
}

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

build_native_plan <- function(structure, param_table) {
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  if (is.null(prep[[".runtime"]]) || is.null(prep$.runtime$cache_bundle)) {
    prep <- AccumulatR:::.prep_set_cache_bundle(
      prep,
      AccumulatR:::.build_likelihood_cache_bundle(prep)
    )
  }
  plan <- AccumulatR:::.likelihood_build_trial_plan(
    structure,
    param_table,
    prep
  )
  list(
    prep = prep,
    prep_native = AccumulatR:::.prep_native_payload(prep),
    plan = plan
  )
}

profile_example <- function(example_id,
                            n_trials = 200L,
                            seed = 2025L) {
  spec <- new_api_examples[[example_id]]
  core_params <- new_api_example_params[[example_id]]
  if (is.null(spec) || is.null(core_params)) {
    stop(sprintf("Example '%s' not found", example_id))
  }

  message(sprintf("Profiling example %s", example_id))
  patched_model <- apply_core_params_to_spec(spec, core_params)
  structure <- build_generator_structure(patched_model)
  param_table <- build_params_df(spec, core_params, n_trials = n_trials)
  data_df <- simulate_data_from_params_table(structure, param_table, seed = seed)
  plan_bundle <- build_native_plan(structure, param_table)

  AccumulatR:::likelihood_cache_reset_stats(plan_bundle$prep)
  elapsed <- system.time({
    with_native_flags(TRUE, TRUE, {
      trial_keys <- as.character(data_df$trial)
      trial_keys[is.na(trial_keys)] <- NA_character_
      data_row_indices <- split(seq_len(nrow(data_df)), trial_keys)
      trial_ids <- unique(data_df$trial)
      AccumulatR:::.native_loglikelihood_batch(
        structure = structure,
        prep = plan_bundle$prep_native,
        plan = plan_bundle$plan,
        trial_ids = trial_ids,
        data_df = data_df,
        data_row_indices = data_row_indices,
        component_weights = NULL
      )
    })
  })[["elapsed"]]

  message(sprintf("Example %s native log-likelihood took %.3f seconds", example_id, elapsed))
  invisible(elapsed)
}

example_ids <- names(new_api_examples)[1]
trial_count <- as.integer(Sys.getenv("UUBER_PROFILE_TRIALS", "200"))
seed <- as.integer(Sys.getenv("UUBER_PROFILE_SEED", "2025"))

timings <- vapply(example_ids, profile_example, numeric(1), n_trials = trial_count, seed = seed)

print(timings)
