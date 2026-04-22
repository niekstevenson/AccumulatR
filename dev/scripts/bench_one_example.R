#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: bench_one_example.R <label> <n_trials> <cap_ms> <target_batch_ms> <out_csv>")
}

label <- args[[1]]
n_trials <- as.integer(args[[2]])
cap_ms <- as.numeric(args[[3]])
target_batch_ms <- as.numeric(args[[4]])
out_csv <- args[[5]]

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

script_file <- {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- all_args[grepl("^--file=", all_args)]
  if (length(file_arg) == 0L) {
    stop("This benchmark script must be run with Rscript --file")
  }
  sub("^--file=", "", file_arg[[1]])
}

repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

suppressPackageStartupMessages({
  library(pkgload)
})

load_all(repo_root, quiet = TRUE, helpers = FALSE)
source(file.path("dev", "examples", "new_API.R"))

bench_one <- function(label, structure, pars) {
  params_df <- build_param_matrix(structure, pars, n_trials = n_trials)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  params_slim <- build_param_matrix(structure, pars, trial_df = prepared)

  invisible(log_likelihood(ctx, prepared, params_slim))

  pilot_times <- numeric(1)
  pilot_times[[1]] <- system.time(
    invisible(log_likelihood(ctx, prepared, params_slim))
  )[["elapsed"]] * 1000
  if (pilot_times[[1]] < cap_ms / 2) {
    pilot_times <- c(
      pilot_times,
      replicate(
        2,
        system.time(invisible(log_likelihood(ctx, prepared, params_slim)))[["elapsed"]] * 1000
      )
    )
  }
  pilot_ms <- median(pilot_times)

  capped <- pilot_ms >= cap_ms
  inner_reps <- if (capped) {
    1L
  } else {
    max(1L, min(250L, as.integer(ceiling(target_batch_ms / max(pilot_ms, 0.1)))))
  }
  outer_reps <- if (capped) 1L else 5L

  times <- numeric(outer_reps)
  for (i in seq_len(outer_reps)) {
    gc(FALSE)
    times[[i]] <- system.time({
      for (j in seq_len(inner_reps)) {
        invisible(log_likelihood(ctx, prepared, params_slim))
      }
    })[["elapsed"]] * 1000 / inner_reps
  }

  data.frame(
    label = label,
    per_eval_ms = median(times),
    batch_ms = median(times) * inner_reps,
    inner_reps = inner_reps,
    capped = capped,
    status = "ok",
    error = "",
    stringsAsFactors = FALSE
  )
}

if (!label %in% names(new_api_examples)) {
  stop("Unknown example label: ", label)
}

row <- tryCatch(
  bench_one(label, new_api_examples[[label]], new_api_example_params[[label]]),
  error = function(e) data.frame(
    label = label,
    per_eval_ms = NA_real_,
    batch_ms = NA_real_,
    inner_reps = NA_integer_,
    capped = NA,
    status = "error",
    error = conditionMessage(e),
    stringsAsFactors = FALSE
  )
)

utils::write.csv(row, out_csv, row.names = FALSE)
