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

args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args) >= 1L && nzchar(args[[1]])) {
  normalizePath(args[[1]])
} else {
  repo_root
}
source_label <- if (length(args) >= 2L && nzchar(args[[2]])) {
  args[[2]]
} else {
  basename(source_dir)
}

bench_lib <- tempfile("accumulatr_milestone2_lib_")
dir.create(bench_lib, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(bench_lib, recursive = TRUE, force = TRUE), add = TRUE)

install_out <- system2(
  "R",
  args = c(
    "CMD", "INSTALL", "--preclean", "--no-multiarch",
    paste0("--library=", bench_lib),
    source_dir
  ),
  stdout = TRUE,
  stderr = TRUE
)
install_status <- attr(install_out, "status") %||% 0L
if (!identical(install_status, 0L)) {
  cat(paste(install_out, collapse = "\n"), "\n")
  stop("Failed to install benchmark source tree")
}

suppressPackageStartupMessages(library(AccumulatR, lib.loc = bench_lib))

n_rep <- as.integer(Sys.getenv("ACCUMULATR_M2_BENCH_N_REP", "7"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 7L
}

run_bench <- function(label, fn, value_fn = NULL, inner_reps) {
  timings <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    gc(FALSE)
    timings[[i]] <- system.time({
      for (j in seq_len(inner_reps)) {
        fn()
      }
    })[["elapsed"]]
  }
  value <- if (is.null(value_fn)) NA_character_ else {
    paste(capture.output(print(value_fn())), collapse = " ")
  }
  data.frame(
    source_label = source_label,
    case = label,
    median_sec = stats::median(timings),
    min_sec = min(timings),
    max_sec = max(timings),
    inner_reps = inner_reps,
    median_per_eval_sec = stats::median(timings) / inner_reps,
    throughput_eval_per_sec = inner_reps / stats::median(timings),
    value = value,
    stringsAsFactors = FALSE
  )
}

build_ranked_spec <- function() {
  race_spec(n_outcomes = 2L) |>
    add_accumulator("a1", "lognormal") |>
    add_accumulator("a2", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("A", c("a1", "a2"), k = 1L) |>
    add_outcome("A", "A") |>
    add_outcome("B", "b") |>
    add_trigger("tg", members = c("a1", "a2"), q = 0.10, draw = "shared") |>
    finalize_model()
}

ranked_params <- c(
  a1.meanlog = log(0.30), a1.sdlog = 0.20, a1.q = 0.10, a1.t0 = 0.00,
  a2.meanlog = log(0.33), a2.sdlog = 0.20, a2.q = 0.10, a2.t0 = 0.00,
  b.meanlog = log(0.52), b.sdlog = 0.18, b.q = 0.00, b.t0 = 0.00
)

make_ranked_case <- function(label, rt1, rt2, inner_reps = 30L) {
  structure <- build_ranked_spec()
  trial_count <- length(rt1)
  ranked_df <- data.frame(
    trial = seq_len(trial_count),
    R = factor(rep("A", trial_count), levels = c("A", "B")),
    rt = rt1,
    R2 = factor(rep("B", trial_count), levels = c("A", "B")),
    rt2 = rt2
  )
  params_df <- build_param_matrix(structure, ranked_params, n_trials = trial_count)
  ctx <- build_likelihood_context(structure, ranked_df)
  fn <- function() log_likelihood(ctx, ranked_df, params_df)
  list(label = label, fn = fn, value_fn = fn, inner_reps = inner_reps)
}

build_beest_spec <- function() {
  race_spec() |>
    add_accumulator("S", "exgauss") |>
    add_accumulator("stop", "exgauss") |>
    add_accumulator("change", "exgauss") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_component("go_only", members = "S", weight = 0.75) |>
    add_component("go_stop", members = c("S", "stop", "change"), weight = 0.25) |>
    add_trigger("tg", members = c("stop", "change"), q = 0.10, param = "tg_q", draw = "shared") |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()
}

exact_params <- c(
  S.mu = 0.38, S.sigma = 0.05, S.tau = 0.06,
  stop.mu = 0.08, stop.sigma = 0.16, stop.tau = 0.09,
  change.mu = 0.33, change.sigma = 0.09, change.tau = 0.05,
  tg_q = 0.22
)

make_exact_case <- function(label, rt_shift = 0, inner_reps = 20L) {
  structure <- build_beest_spec()
  base <- data.frame(
    trial = c(1L, 1L, 1L),
    R = factor(c("X", "X", "X"), levels = c("S", "X")),
    rt = c(0.39, 0.39, 0.39) + rt_shift,
    component = factor(c("go_stop", "go_stop", "go_stop"),
      levels = c("go_only", "go_stop")
    ),
    accumulator = c("S", "stop", "change"),
    onset = c(0.00, 0.10, 0.10)
  )
  data_df <- do.call(rbind, lapply(seq_len(90L), function(i) {
    out <- base
    out$trial <- i
    out
  }))
  params_df <- build_param_matrix(structure, exact_params, n_trials = 90L)
  ctx <- build_likelihood_context(structure, data_df)
  fn <- function() log_likelihood(ctx, data_df, params_df)
  list(label = label, fn = fn, value_fn = fn, inner_reps = inner_reps)
}

ranked_duplicate_case <- make_ranked_case(
  "ranked_duplicate_loglik",
  rt1 = rep(0.36, 180L),
  rt2 = rep(0.61, 180L)
)

ranked_unique_case <- make_ranked_case(
  "ranked_unique_loglik",
  rt1 = seq(0.36, 0.54, length.out = 180L),
  rt2 = seq(0.61, 0.89, length.out = 180L)
)

exact_duplicate_case <- make_exact_case("exact_guarded_duplicate_loglik", rt_shift = 0)
exact_shifted_case <- make_exact_case("exact_guarded_shifted_loglik", rt_shift = 0.015)

results <- do.call(
  rbind,
  list(
    run_bench(
      ranked_duplicate_case$label,
      ranked_duplicate_case$fn,
      ranked_duplicate_case$value_fn,
      ranked_duplicate_case$inner_reps
    ),
    run_bench(
      ranked_unique_case$label,
      ranked_unique_case$fn,
      ranked_unique_case$value_fn,
      ranked_unique_case$inner_reps
    ),
    run_bench(
      exact_duplicate_case$label,
      exact_duplicate_case$fn,
      exact_duplicate_case$value_fn,
      exact_duplicate_case$inner_reps
    ),
    run_bench(
      exact_shifted_case$label,
      exact_shifted_case$fn,
      exact_shifted_case$value_fn,
      exact_shifted_case$inner_reps
    )
  )
)

print(results)

output_dir <- file.path(repo_root, "dev", "scripts", "scratch_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_dir,
  sprintf("benchmark_milestone2_%s.csv", source_label)
)
utils::write.csv(results, output_file, row.names = FALSE)
cat(sprintf("Wrote %s\n", output_file))
