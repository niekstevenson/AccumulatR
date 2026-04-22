#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

take_arg <- function(flag, default = NULL) {
  prefix <- paste0(flag, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

model_label <- take_arg("--model")
trial_count <- as.integer(take_arg("--trials", "250"))
params_text <- take_arg("--params")
out_path <- take_arg("--out")
frequency <- as.integer(take_arg("--frequency", "500"))
target_seconds <- as.numeric(take_arg("--seconds", "5"))
min_ll <- as.numeric(take_arg("--min-ll", as.character(log(1e-10))))

if (is.null(model_label) || !nzchar(model_label)) {
  stop("missing --model=<label>")
}
if (!is.finite(trial_count) || trial_count <= 0L) {
  stop("--trials must be a positive integer")
}
if (!is.finite(frequency) || frequency <= 0L) {
  stop("--frequency must be a positive integer")
}
if (!is.finite(target_seconds) || target_seconds <= 0) {
  stop("--seconds must be positive")
}

script_file <- {
  full <- commandArgs(trailingOnly = FALSE)
  file_arg <- full[grepl("^--file=", full)]
  if (length(file_arg) == 0L) {
    stop("this script must be run with Rscript --file")
  }
  sub("^--file=", "", file_arg[[1]])
}
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

suppressPackageStartupMessages({
  library(AccumulatR)
})
source(file.path("dev", "examples", "new_API.R"))

if (!exists("new_api_examples", inherits = TRUE) ||
    !exists("new_api_example_params", inherits = TRUE)) {
  stop("example registry not found after sourcing dev/examples/new_API.R")
}

structure <- new_api_examples[[model_label]]
default_params <- new_api_example_params[[model_label]]
if (is.null(structure) && identical(model_label, "example_23_ranked_chain")) {
  structure <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    finalize_model()
  default_params <- c(a.m = log(0.30), a.s = 0.16, b.m = log(0.22), b.s = 0.16)
}
if (is.null(structure)) {
  stop(sprintf("unknown model '%s'", model_label))
}

params <- if (!is.null(params_text) && nzchar(params_text)) {
  eval(parse(text = params_text), envir = baseenv())
} else {
  default_params
}
if (is.null(params)) {
  stop(sprintf("no default parameters found for model '%s'; pass --params", model_label))
}

profile_dir <- file.path("dev", "scripts", "scratch_outputs", "profiles")
if (!dir.exists(profile_dir)) {
  dir.create(profile_dir, recursive = TRUE)
}
if (is.null(out_path) || !nzchar(out_path)) {
  out_path <- file.path(
    profile_dir,
    sprintf("%s_%s.prof", model_label, format(Sys.time(), "%Y%m%d-%H%M%S"))
  )
}

params_df <- build_param_matrix(structure, params, n_trials = trial_count)
data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
prepared <- prepare_data(structure, data_df)
ctx <- make_context(structure)
params_slim <- build_param_matrix(structure, params, trial_df = prepared)

if (!AccumulatR:::`.profiler_available`()) {
  stop("native profiler bridge is not available")
}

invisible(log_likelihood(ctx, prepared, params_slim, min_ll = min_ll))
pilot <- median(replicate(
  3,
  system.time(invisible(log_likelihood(ctx, prepared, params_slim, min_ll = min_ll)))[["elapsed"]]
))
reps <- max(1L, as.integer(ceiling(target_seconds / max(pilot, 1e-4))))

if (file.exists(out_path)) {
  invisible(file.remove(out_path))
}

invisible(AccumulatR:::`.profiler_start`(out_path, frequency = frequency))
on.exit(invisible(AccumulatR:::`.profiler_stop`()), add = TRUE)

for (i in seq_len(reps)) {
  invisible(log_likelihood(ctx, prepared, params_slim, min_ll = min_ll))
}

invisible(AccumulatR:::`.profiler_flush`())
invisible(AccumulatR:::`.profiler_stop`())

cat(sprintf("profile=%s\n", normalizePath(out_path, mustWork = FALSE)))
cat(sprintf("model=%s\n", model_label))
cat(sprintf("trials=%d\n", trial_count))
cat(sprintf("frequency=%d\n", frequency))
cat(sprintf("iterations=%d\n", reps))
