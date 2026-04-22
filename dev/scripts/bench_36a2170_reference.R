#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BENCH_SEED", "20260217")))

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
old_csv <- file.path(out_dir, "benchmark_36a2170_reference_old.csv")
new_csv <- file.path(out_dir, "benchmark_36a2170_reference_new.csv")

n_trials <- as.integer(Sys.getenv("ACCUMULATR_BENCH_TRIALS", "250"))
if (!is.finite(n_trials) || n_trials < 50L) {
  n_trials <- 250L
}
n_rep <- as.integer(Sys.getenv("ACCUMULATR_BENCH_N_REP", "4"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 4L
}
target_sample_sec <- as.numeric(Sys.getenv("ACCUMULATR_BENCH_TARGET_SEC", "0.08"))
if (!is.finite(target_sample_sec) || target_sample_sec <= 0) {
  target_sample_sec <- 0.08
}
min_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MIN_INNER_REPS", "5"))
if (!is.finite(min_inner_reps) || min_inner_reps < 1L) {
  min_inner_reps <- 5L
}
max_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MAX_INNER_REPS", "200"))
if (!is.finite(max_inner_reps) || max_inner_reps < min_inner_reps) {
  max_inner_reps <- max(1000L, min_inner_reps)
}

choose_inner_reps <- function(fn, target_sec, min_reps, max_reps) {
  calib_reps <- 1L
  elapsed <- 0.0
  while ((elapsed <= 0 || !is.finite(elapsed)) && calib_reps <= 4096L) {
    elapsed <- system.time({
      for (i in seq_len(calib_reps)) {
        fn()
      }
    })[["elapsed"]]
    if (elapsed <= 0 || !is.finite(elapsed)) {
      calib_reps <- calib_reps * 2L
    }
  }
  if (!is.finite(elapsed) || elapsed <= 0) {
    return(as.integer(min_reps))
  }
  per_eval <- elapsed / calib_reps
  if (!is.finite(per_eval) || per_eval <= 0) {
    return(as.integer(min_reps))
  }
  est <- ceiling(target_sec / per_eval)
  est <- max(min_reps, est)
  est <- min(max_reps, est)
  as.integer(est)
}

run_bench <- function(label, fn, n_rep, inner_reps) {
  timings <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    gc(FALSE)
    timings[[i]] <- system.time({
      for (j in seq_len(inner_reps)) {
        fn()
      }
    })[["elapsed"]]
  }
  data.frame(
    label = label,
    median_sec = stats::median(timings),
    min_sec = min(timings),
    max_sec = max(timings),
    inner_reps = inner_reps,
    median_per_eval_sec = stats::median(timings) / inner_reps,
    mode = "ir",
    stringsAsFactors = FALSE
  )
}

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
  example_6_dual_path = list(
    structure = new_api_examples[["example_6_dual_path"]],
    pars = new_api_example_params[["example_6_dual_path"]]
  ),
  example_7_mixture = list(
    structure = new_api_examples[["example_7_mixture"]],
    pars = new_api_example_params[["example_7_mixture"]]
  ),
  example_10_exclusion = list(
    structure = new_api_examples[["example_10_exclusion"]],
    pars = new_api_example_params[["example_10_exclusion"]]
  ),
  example_16_guard_tie_simple = list(
    structure = new_api_examples[["example_15_guard_tie_simple"]],
    pars = new_api_example_params[["example_15_guard_tie_simple"]]
  ),
  example_21_simple_q = list(
    structure = new_api_examples[["example_21_simple_q"]],
    pars = new_api_example_params[["example_21_simple_q"]]
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

make_case <- function(mod, n_trials) {
  params_df <- build_param_matrix(mod$structure, mod$pars, n_trials = n_trials)
  data_df <- simulate(mod$structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(mod$structure, data_df)
  ctx <- make_context(mod$structure)
  params_slim <- build_param_matrix(mod$structure, mod$pars, trial_df = prepared)
  ll <- as.numeric(log_likelihood(ctx, prepared, params_slim))
  if (length(ll) != 1L || !is.finite(ll)) {
    stop("Non-finite log-likelihood during setup")
  }
  list(ctx = ctx, prepared = prepared, params = params_slim)
}

cases <- lapply(models, make_case, n_trials = n_trials)
case_inner_reps <- setNames(integer(length(cases)), names(cases))
for (case_name in names(cases)) {
  case <- cases[[case_name]]
  eval_fn <- function() {
    log_likelihood(case$ctx, case$prepared, case$params)
  }
  case_inner_reps[[case_name]] <- choose_inner_reps(
    fn = eval_fn,
    target_sec = target_sample_sec,
    min_reps = min_inner_reps,
    max_reps = max_inner_reps
  )
}

new_out <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    case <- cases[[case_name]]
    eval_fn <- function() {
      log_likelihood(case$ctx, case$prepared, case$params)
    }
    invisible(eval_fn())
    run_bench(
      label = case_name,
      fn = eval_fn,
      n_rep = n_rep,
      inner_reps = case_inner_reps[[case_name]]
    )
  })
)

old_text <- system2(
  "git",
  c("show", "36a2170:dev/scripts/scratch_outputs/bench_examples_models_likelihood.csv"),
  stdout = TRUE,
  stderr = TRUE
)
old_status <- attr(old_text, "status") %||% 0L
if (old_status != 0L) {
  stop("Failed to read historical benchmark CSV from git")
}
old_con <- textConnection(old_text)
on.exit(close(old_con), add = TRUE)
old_out <- utils::read.csv(old_con, stringsAsFactors = FALSE)

utils::write.csv(old_out, old_csv, row.names = FALSE)
utils::write.csv(new_out, new_csv, row.names = FALSE)

cat("Wrote historical reference to", old_csv, "\n")
cat("Wrote current reference to", new_csv, "\n")
