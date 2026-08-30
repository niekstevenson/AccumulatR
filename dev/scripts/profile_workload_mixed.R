#!/usr/bin/env Rscript

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) ""
)
repo_root <- normalizePath(file.path(dirname(script_path %||% "."), "..", ".."),
                           mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "dev", "examples"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

load_accumulatr <- function() {
  if (requireNamespace("AccumulatR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(AccumulatR))
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(pkgload))
  load_all(repo_root, quiet = TRUE, helpers = FALSE)
}

load_accumulatr()
source(file.path("dev", "examples", "benchmark_models.R"))

parse_int_env <- function(name, default, min_value = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < min_value) {
    return(as.integer(default))
  }
  value
}

parse_num_env <- function(name, default, min_value = 0) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value <= min_value) {
    return(default)
  }
  value
}

models <- benchmark_models

balanced_labels <- c(
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

case_env <- Sys.getenv("ACCUMULATR_PROFILE_CASES", "")
case_labels <- if (nzchar(case_env)) {
  trimws(strsplit(case_env, ",", fixed = TRUE)[[1]])
} else {
  balanced_labels
}
case_labels <- case_labels[nzchar(case_labels)]
unknown <- setdiff(case_labels, names(models))
if (length(unknown) > 0L) {
  stop("Unknown profile case(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

n_trials <- parse_int_env("ACCUMULATR_PROFILE_TRIALS", 50L)
workload_seconds <- parse_num_env(
  "ACCUMULATR_PROFILE_WORKLOAD_SECONDS",
  parse_num_env("ACCUMULATR_PROFILE_DURATION", 20)
)

make_case <- function(mod, n_trials) {
  params_df <- build_param_matrix(mod$structure, mod$pars, n_trials = n_trials)
  data_df <- simulate(mod$structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(mod$structure, data_df)
  ctx <- make_context(mod$structure)
  params_slim <- build_param_matrix(mod$structure, mod$pars, trial_df = prepared)
  value <- as.numeric(log_likelihood(ctx, prepared, params_slim))
  if (length(value) != 1L || !is.finite(value)) {
    stop("Non-finite log-likelihood during profile setup", call. = FALSE)
  }
  list(ctx = ctx, prepared = prepared, params = params_slim)
}

cases <- lapply(models[case_labels], make_case, n_trials = n_trials)
counts <- setNames(integer(length(cases)), names(cases))

start_file <- Sys.getenv("ACCUMULATR_PROFILE_START_FILE", "")
end_file <- Sys.getenv("ACCUMULATR_PROFILE_END_FILE", "")
if (nzchar(start_file)) {
  writeLines("start", start_file)
}
on.exit({
  if (nzchar(end_file)) {
    writeLines("end", end_file)
  }
}, add = TRUE)

deadline <- proc.time()[["elapsed"]] + workload_seconds
repeat {
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    value <- as.numeric(log_likelihood(case$ctx, case$prepared, case$params))
    if (length(value) != 1L || !is.finite(value)) {
      stop("Non-finite log-likelihood during profile loop", call. = FALSE)
    }
    counts[[case_name]] <- counts[[case_name]] + 1L
  }
  if (proc.time()[["elapsed"]] >= deadline) {
    break
  }
}

cat("Mixed profile workload completed\n")
cat("Trials per case:", n_trials, "\n")
print(counts)
