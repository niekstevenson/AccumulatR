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

bench_lib <- tempfile("accumulatr_interval_self_lib_")
dir.create(bench_lib, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(bench_lib, recursive = TRUE, force = TRUE), add = TRUE)

install_out <- system2(
  "R",
  args = c(
    "CMD", "INSTALL", "--preclean", "--no-multiarch",
    paste0("--library=", bench_lib),
    repo_root
  ),
  stdout = TRUE,
  stderr = TRUE
)
install_status <- attr(install_out, "status") %||% 0L
if (!identical(install_status, 0L)) {
  cat(paste(install_out, collapse = "\n"), "\n")
  stop("Failed to install local AccumulatR into temporary library")
}

library(AccumulatR, lib.loc = bench_lib)

n_trials <- as.integer(Sys.getenv("ACCUMULATR_INTERVAL_SELF_TRIALS", "250"))
if (!is.finite(n_trials) || n_trials < 4L) {
  n_trials <- 250L
}
n_rep <- as.integer(Sys.getenv("ACCUMULATR_INTERVAL_SELF_N_REP", "6"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 6L
}
inner_reps <- as.integer(Sys.getenv("ACCUMULATR_INTERVAL_SELF_INNER_REPS", "40"))
if (!is.finite(inner_reps) || inner_reps < 1L) {
  inner_reps <- 40L
}

build_interval_self_spec <- function() {
  race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_accumulator("c", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_outcome("C", "c")
}

params_interval_self <- c(
  a.meanlog = log(0.30), a.sdlog = 0.16, a.q = 0, a.t0 = 0,
  b.meanlog = log(0.22), b.sdlog = 0.16, b.q = 0, b.t0 = 0,
  c.meanlog = log(0.35), c.sdlog = 0.18, c.q = 0, c.t0 = 0
)

make_ranked_df <- function(first, second, n_trials) {
  idx <- seq_len(n_trials)
  rt1 <- rep(c(0.42, 0.48, 0.55, 0.63), length.out = n_trials)
  rt2 <- rep(c(0.70, 0.76, 0.83, 0.92), length.out = n_trials)
  data.frame(
    trial = idx,
    R = factor(rep(first, n_trials), levels = c("A", "B", "C")),
    rt = rt1,
    R2 = factor(rep(second, n_trials), levels = c("A", "B", "C")),
    rt2 = rt2,
    stringsAsFactors = FALSE
  )
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
    stringsAsFactors = FALSE
  )
}

make_case <- function(first, second) {
  spec <- build_interval_self_spec()
  structure <- finalize_model(spec)
  data_df <- make_ranked_df(first, second, n_trials)
  ctx <- build_likelihood_context(structure, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  params_slim <- build_param_matrix(
    spec,
    params_interval_self,
    n_trials = n_trials,
    layout = ctx$param_layout
  )
  AccumulatR:::unified_outcome_stats_reset_cpp()
  ll <- as.numeric(log_likelihood(ctx, data_prepped, params_slim))
  stats <- AccumulatR:::unified_outcome_stats_cpp()
  if (!is.finite(ll)) {
    stop("Non-finite log-likelihood for interval-self case")
  }
  list(
    ctx = ctx,
    data_prepped = data_prepped,
    params_df = params_slim,
    stats = stats
  )
}

cases <- list(
  interval_self_B_then_C = make_case("B", "C"),
  interval_self_C_then_A = make_case("C", "A")
)

bench_rows <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    case <- cases[[case_name]]
    row <- run_bench(
      case_name,
      function() log_likelihood(case$ctx, case$data_prepped, case$params_df),
      n_rep = n_rep,
      inner_reps = inner_reps
    )
    row$exact_node_batch_calls <- case$stats$exact_node_batch_calls
    row$pointwise_fallback_calls <- case$stats$exact_node_batch_pointwise_fallback_calls
    row
  })
)

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, "benchmark_interval_self_ranked.csv")
utils::write.csv(bench_rows, out_file, row.names = FALSE)

print(bench_rows, row.names = FALSE)
cat("\nWrote interval-self benchmark output to", out_file, "\n")
