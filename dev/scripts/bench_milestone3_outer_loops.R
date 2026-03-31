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

bench_lib <- tempfile("accumulatr_milestone3_lib_")
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

n_rep <- as.integer(Sys.getenv("ACCUMULATR_M3_BENCH_N_REP", "7"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 7L
}
inner_reps <- as.integer(Sys.getenv("ACCUMULATR_M3_BENCH_INNER_REPS", "80"))
if (!is.finite(inner_reps) || inner_reps < 1L) {
  inner_reps <- 80L
}
nonresponse_trials <- as.integer(Sys.getenv("ACCUMULATR_M3_BENCH_NONRESP_TRIALS", "180"))
if (!is.finite(nonresponse_trials) || nonresponse_trials < 4L) {
  nonresponse_trials <- 180L
}

run_bench <- function(label, fn, value_fn = NULL, n_rep, inner_reps) {
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

build_example_5_timeout_guess <- function() {
  race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_accumulator("timeout", "lognormal", onset = 0.05) |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
    )) |>
    finalize_model()
}

params_example_5_timeout_guess <- c(
  go_left.meanlog = log(0.30),
  go_left.sdlog = 0.18,
  go_right.meanlog = log(0.325),
  go_right.sdlog = 0.18,
  timeout.meanlog = log(0.25),
  timeout.sdlog = 0.10
)

build_example_7_mixture <- function() {
  race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("TARGET", c("target_fast", "target_slow")) |>
    add_outcome("R1", "TARGET") |>
    add_outcome("R2", "competitor") |>
    add_component("fast", members = c("target_fast", "competitor"),
      weight_param = "p_fast"
    ) |>
    add_component("slow", members = c("target_slow", "competitor")) |>
    set_mixture_options(mode = "sample", reference = "slow") |>
    finalize_model()
}

params_example_7_mixture <- c(
  target_fast.meanlog = log(0.25),
  target_fast.sdlog = 0.15,
  target_slow.meanlog = log(0.45),
  target_slow.sdlog = 0.20,
  competitor.meanlog = log(0.35),
  competitor.sdlog = 0.18,
  p_fast = 0.2
)

build_shared_gate_pair_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP", all_of("x2", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"),
      options = list(map_outcome_to = NA_character_)
    )
}

params_shared_gate_pair <- c(
  x1.meanlog = log(0.32), x1.sdlog = 0.18,
  x2.meanlog = log(0.36), x2.sdlog = 0.18,
  gate.meanlog = log(0.24), gate.sdlog = 0.14
)

build_shared_trigger_guess_spec <- function() {
  race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_accumulator("timeout", "lognormal", onset = 0.05) |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
    )) |>
    add_trigger("tg",
      members = c("go_left", "go_right", "timeout"),
      q = 0.10, param = "tg_q", draw = "shared"
    )
}

params_shared_trigger_guess <- c(
  go_left.meanlog = log(0.30),
  go_left.sdlog = 0.18,
  go_right.meanlog = log(0.325),
  go_right.sdlog = 0.18,
  timeout.meanlog = log(0.25),
  timeout.sdlog = 0.10,
  tg_q = 0.21
)

make_mixture_native_case <- function() {
  structure <- build_example_7_mixture()
  params_df <- build_param_matrix(structure, params_example_7_mixture, n_trials = 1L)
  rows_df <- getFromNamespace(".param_matrix_to_rows", "AccumulatR")(structure, params_df)
  prep <- getFromNamespace(".prepare_model_for_likelihood", "AccumulatR")(
    getFromNamespace(".model_spec_with_params", "AccumulatR")(structure$model_spec, rows_df)
  )
  out_def <- prep[["outcomes"]][[1]]
  compiled <- getFromNamespace(".expr_lookup_compiled", "AccumulatR")(out_def[["expr"]], prep)
  competitor_map <- getFromNamespace(".prep_competitors", "AccumulatR")(prep) %||% list()
  comp_exprs <- competitor_map[["R1"]] %||% list()
  comp_nodes <- lapply(comp_exprs, function(ex) {
    getFromNamespace(".expr_lookup_compiled", "AccumulatR")(ex, prep)
  })
  comp_ids <- if (!any(vapply(comp_nodes, is.null, logical(1)))) {
    vapply(comp_nodes, function(node) as.integer(node[["id"]] %||% NA_integer_), integer(1))
  } else {
    integer(0)
  }
  native_ctx <- getFromNamespace(".prep_native_context", "AccumulatR")(prep)
  native_trial_mixture_cpp <- getFromNamespace("native_trial_mixture_cpp", "AccumulatR")
  mix_weights <- c(params_example_7_mixture[["p_fast"]],
                   1.0 - params_example_7_mixture[["p_fast"]])
  fn <- function() {
    native_trial_mixture_cpp(
      native_ctx,
      as.integer(compiled[["id"]]),
      0.35,
      c("fast", "slow"),
      mix_weights,
      NULL,
      as.integer(comp_ids),
      rows_df,
      NULL
    )
  }
  list(label = "mixture_native_helper", fn = fn, value_fn = fn, inner_reps = 4000L)
}

make_mixture_response_case <- function() {
  structure <- build_example_7_mixture()
  params_df <- build_param_matrix(structure, params_example_7_mixture, n_trials = 1L)
  fn <- function() response_probabilities(structure, params_df, include_na = TRUE)
  list(label = "mixture_response_probabilities", fn = fn, value_fn = fn, inner_reps = 120L)
}

make_guess_response_case <- function() {
  structure <- build_example_5_timeout_guess()
  params_df <- build_param_matrix(structure, params_example_5_timeout_guess, n_trials = 1L)
  fn <- function() response_probabilities(structure, params_df, include_na = TRUE)
  list(label = "guess_response_probabilities", fn = fn, value_fn = fn, inner_reps = 120L)
}

make_nonresponse_loglik_case <- function() {
  structure <- finalize_model(build_shared_gate_pair_spec())
  nonresp_df <- data.frame(
    trial = seq_len(nonresponse_trials),
    R = factor(rep(NA_character_, nonresponse_trials), levels = c("RESP", "NR_RAW")),
    rt = rep(NA_real_, nonresponse_trials)
  )
  params_df <- build_param_matrix(
    structure,
    params_shared_gate_pair,
    n_trials = nonresponse_trials
  )
  ctx <- build_likelihood_context(structure, nonresp_df)
  fn <- function() log_likelihood(ctx, nonresp_df, params_df)
  list(label = "nonresponse_loglik", fn = fn, value_fn = fn, inner_reps = 40L)
}

make_shared_trigger_guess_loglik_case <- function() {
  spec <- build_shared_trigger_guess_spec()
  structure <- finalize_model(spec)
  trial_count <- max(40L, nonresponse_trials)
  trial_data <- data.frame(
    trial = seq_len(trial_count),
    R = factor(rep("Left", trial_count), levels = c("Left", "Right", "TIMEOUT")),
    rt = rep(c(0.34, 0.35, 0.36, 0.37), length.out = trial_count)
  )
  params_df <- build_param_matrix(
    spec,
    params_shared_trigger_guess,
    n_trials = trial_count
  )
  ctx <- build_likelihood_context(structure, trial_data)
  fn <- function() log_likelihood(ctx, trial_data, params_df)
  list(label = "shared_trigger_guess_loglik", fn = fn, value_fn = fn, inner_reps = 40L)
}

cases <- list(
  make_mixture_native_case(),
  make_mixture_response_case(),
  make_guess_response_case(),
  make_nonresponse_loglik_case(),
  make_shared_trigger_guess_loglik_case()
)

rows <- do.call(
  rbind,
  lapply(cases, function(cs) {
    run_bench(
      cs$label,
      cs$fn,
      cs$value_fn,
      n_rep = n_rep,
      inner_reps = cs$inner_reps %||% inner_reps
    )
  })
)

out_dir <- file.path(repo_root, "dev", "scripts", "scratch_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("benchmark_milestone3_%s.csv", source_label))
utils::write.csv(rows, out_file, row.names = FALSE)

print(rows, row.names = FALSE)
cat("\nWrote milestone-3 benchmark output to", out_file, "\n")
