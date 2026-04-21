#!/usr/bin/env Rscript

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

script_file <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
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

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

old_csv <- file.path(out_dir, "benchmark_fast_reference_old.csv")
new_csv <- file.path(out_dir, "benchmark_fast_reference_new.csv")

selected_labels <- c(
  "example_1_simple",
  "example_7_mixture",
  "example_22_shared_q",
  "example_5_timeout_guess",
  "example_10_exclusion",
  "example_3_stop_na",
  "example_23_ranked_chain",
  "example_2_stop_mixture"
)

old_rows <- utils::read.csv(
  file.path(out_dir, "benchmark_centralized.csv"),
  stringsAsFactors = FALSE
)
old_rows <- subset(old_rows, section == "speed" & suite == "examples_models")
old_metrics <- split(old_rows, old_rows$label)
old_metrics <- lapply(old_metrics, function(df) setNames(df$value, df$metric))

models <- list(
  example_1_simple = list(
    structure = new_api_examples[["example_1_simple"]],
    pars = new_api_example_params[["example_1_simple"]]
  ),
  example_2_stop_mixture = list(
    structure = new_api_examples[["example_2_stop_mixture"]],
    pars = new_api_example_params[["example_2_stop_mixture"]]
  ),
  example_3_stop_na = list(
    structure = new_api_examples[["example_3_stop_na"]],
    pars = new_api_example_params[["example_3_stop_na"]]
  ),
  example_5_timeout_guess = list(
    structure = new_api_examples[["example_5_timeout_guess"]],
    pars = new_api_example_params[["example_5_timeout_guess"]]
  ),
  example_7_mixture = list(
    structure = new_api_examples[["example_7_mixture"]],
    pars = new_api_example_params[["example_7_mixture"]]
  ),
  example_10_exclusion = list(
    structure = new_api_examples[["example_10_exclusion"]],
    pars = new_api_example_params[["example_10_exclusion"]]
  ),
  example_22_shared_q = list(
    structure = new_api_examples[["example_22_shared_q"]],
    pars = new_api_example_params[["example_22_shared_q"]]
  ),
  example_23_ranked_chain = local({
    structure <- race_spec(n_outcomes = 2L) |>
      add_accumulator("a", "lognormal") |>
      add_accumulator("b", "lognormal", onset = after("a")) |>
      add_outcome("A", "a") |>
      add_outcome("B", "b") |>
      finalize_model()
    pars <- c(a.m = log(0.30), a.s = 0.16, b.m = log(0.22), b.s = 0.16)
    list(structure = structure, pars = pars)
  })
)

bench_one <- function(label, structure, pars, n_trials = 250L) {
  params_df <- build_param_matrix(structure, pars, n_trials = n_trials)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  params_slim <- build_param_matrix(structure, pars, trial_df = prepared)
  invisible(log_likelihood(ctx, prepared, params_slim))

  pilot <- median(replicate(
    3,
    system.time(invisible(log_likelihood(ctx, prepared, params_slim)))[["elapsed"]]
  ))
  inner_reps <- max(1L, min(250L, as.integer(ceiling(0.4 / max(pilot, 1e-4)))))

  times <- numeric(5)
  for (i in seq_along(times)) {
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
    stringsAsFactors = FALSE
  )
}

old_out <- do.call(
  rbind,
  lapply(seq_along(selected_labels), function(i) {
    label <- selected_labels[[i]]
    metrics <- old_metrics[[label]]
    data.frame(
      rank = i,
      label = label,
      per_eval_ms = as.numeric(metrics[["median_per_eval_sec"]]) * 1000,
      batch_ms = as.numeric(metrics[["median_sec"]]) * 1000,
      inner_reps = as.integer(metrics[["inner_reps"]]),
      stringsAsFactors = FALSE
    )
  })
)

new_out <- do.call(
  rbind,
  lapply(seq_along(selected_labels), function(i) {
    label <- selected_labels[[i]]
    row <- bench_one(label, models[[label]]$structure, models[[label]]$pars)
    row$rank <- i
    row[, c("rank", "label", "per_eval_ms", "batch_ms", "inner_reps")]
  })
)

utils::write.csv(old_out, old_csv, row.names = FALSE)
utils::write.csv(new_out, new_csv, row.names = FALSE)

cat("Wrote old benchmark reference to", old_csv, "\n")
cat("Wrote new benchmark reference to", new_csv, "\n")
