rm(list = ls())

repeat_count <- as.integer(Sys.getenv("UUBER_PROFILE_REPEATS", "10"))
n_rep <- 10

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
  if ("params_hash" %in% names(param_table)) {
    param_table$params_hash <- NULL
  }
  data_df <- simulate_data_from_params_table(structure, param_table, seed = seed)
  ctx <- build_likelihood_context(
    structure = structure,
    params_df = param_table,
    data_df = data_df
  )

  AccumulatR:::likelihood_cache_reset_stats(ctx$prep)
  param_tables <- lapply(seq_len(repeat_count), function(i) {
    df <- param_table
    if ("params_hash" %in% names(df)) {
      df$params_hash <- NULL
    }
    num_cols <- vapply(df, is.numeric, logical(1))
    if ("trial" %in% names(df)) num_cols[match("trial", names(df))] <- FALSE
    if (any(num_cols)) {
      delta <- (i - 1) * 1e-4
      for (nm in names(df)[num_cols]) {
        df[[nm]] <- df[[nm]] + delta
      }
    }
    df
  })
  elapsed <- system.time({
    with_native_flags(TRUE, TRUE, {
      AccumulatR:::native_loglikelihood_param_repeat(
        ctx,
        param_tables
      )
    })
  })[["elapsed"]]

  message(sprintf(
    "Example %s native log-likelihood repeat (%d proposals) took %.3f seconds",
    example_id, repeat_count, elapsed
  ))
  invisible(elapsed)
}

example_ids <- names(new_api_examples)[1]
trial_count <- as.integer(Sys.getenv("UUBER_PROFILE_TRIALS", "200"))
seed <- as.integer(Sys.getenv("UUBER_PROFILE_SEED", "2025"))

for(i in 1:n_rep){
  timings <- vapply(example_ids, profile_example, numeric(1), n_trials = trial_count, seed = seed)
}

print(timings)
