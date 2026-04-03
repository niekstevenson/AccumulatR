#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BEEST_BENCH_SEED", "20260326")))

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

bench_lib <- tempfile("accumulatr_beest_bench_lib_")
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

suppressPackageStartupMessages(library(AccumulatR, lib.loc = bench_lib))
source(file.path("tests", "testthat", "helper-beest.R"), local = TRUE)

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, "benchmark_centralized_shared_blocker_beest.csv")

n_rep <- as.integer(Sys.getenv("ACCUMULATR_BEEST_BENCH_N_REP", "5"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 5L
}
target_sample_sec <- as.numeric(Sys.getenv("ACCUMULATR_BEEST_BENCH_TARGET_SEC", "0.08"))
if (!is.finite(target_sample_sec) || target_sample_sec <= 0) {
  target_sample_sec <- 0.08
}
min_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BEEST_BENCH_MIN_INNER_REPS", "3"))
if (!is.finite(min_inner_reps) || min_inner_reps < 1L) {
  min_inner_reps <- 3L
}
max_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BEEST_BENCH_MAX_INNER_REPS", "200"))
if (!is.finite(max_inner_reps) || max_inner_reps < min_inner_reps) {
  max_inner_reps <- max(200L, min_inner_reps)
}
inner_reps_override <- as.integer(Sys.getenv("ACCUMULATR_BEEST_BENCH_INNER_REPS", NA_character_))
if (is.na(inner_reps_override) || !is.finite(inner_reps_override) || inner_reps_override < 1L) {
  inner_reps_override <- NA_integer_
}
fixture_reps <- as.integer(Sys.getenv("ACCUMULATR_BEEST_PATTERN_REPS", "80"))
if (!is.finite(fixture_reps) || fixture_reps < 1L) {
  fixture_reps <- 80L
}
example_trials <- as.integer(Sys.getenv("ACCUMULATR_BEEST_MODEL2_TRIALS", "250"))
if (!is.finite(example_trials) || example_trials < 10L) {
  example_trials <- 250L
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

choose_inner_reps <- function(fn, target_sec, min_reps, max_reps, override = NA_integer_) {
  if (!is.na(override) && is.finite(override) && override >= 1L) {
    return(as.integer(override))
  }
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

flatten_speed <- function(df, suite_name) {
  metrics <- c("median_sec", "min_sec", "max_sec", "inner_reps", "median_per_eval_sec")
  do.call(
    rbind,
    lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, , drop = FALSE]
      data.frame(
        section = "speed",
        suite = suite_name,
        label = as.character(row$label[[1]]),
        metric = metrics,
        value = as.numeric(unlist(row[1, metrics], use.names = FALSE)),
        mode = as.character(row$mode[[1]] %||% ""),
        stringsAsFactors = FALSE
      )
    })
  )
}

build_example_2_case <- function(n_trials) {
  spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("stop", "exgauss", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_outcome("R1", inhibit("go1", by = "stop"),
      options = list(component = c("go_only", "go_stop"))
    ) |>
    add_outcome("R2", all_of("go2", "stop"),
      options = list(component = "go_stop")
    ) |>
    add_component("go_only", members = "go1", weight = 0.5) |>
    add_component("go_stop", members = c("go1", "stop", "go2"), weight = 0.5) |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()

  params <- c(
    go1.meanlog = log(0.35),
    go1.sdlog = 0.20,
    stop.mu = 0.10,
    stop.sigma = 0.04,
    stop.tau = 0.10,
    go2.meanlog = log(0.60),
    go2.sdlog = 0.18
  )

  params_df <- build_param_matrix(spec, params, n_trials = n_trials)
  data_df <- simulate(spec, params_df, seed = 123, keep_component = TRUE)
  ctx <- build_likelihood_context(spec, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  params_df_slim <- build_param_matrix(
    spec,
    params,
    n_trials = max(data_df$trial),
    layout = ctx$param_layout
  )

  eval_fn <- function() log_likelihood(ctx, data_prepped, params_df_slim)
  ll_native <- as.numeric(eval_fn())
  if (length(ll_native) != 1L || !is.finite(ll_native)) {
    stop("Non-finite model-2 benchmark log-likelihood during setup")
  }

  list(
    label = "example_2_stop_mixture",
    eval_fn = eval_fn,
    native_loglik = ll_native
  )
}

make_beest_fixture_data <- function(n_reps) {
  base_rows <- data.frame(
    pattern_trial = c(1L, 2L, 2L, 2L, 3L, 3L, 3L),
    R = c("S", "S", "S", "S", "X", "X", "X"),
    rt = c(0.43, 0.47, 0.47, 0.47, 0.39, 0.39, 0.39),
    component = c("go_only", "go_stop", "go_stop", "go_stop",
                  "go_stop", "go_stop", "go_stop"),
    accumulator = c("S", "S", "stop", "change", "S", "stop", "change"),
    onset = c(0.00, 0.00, 0.12, 0.12, 0.00, 0.10, 0.10),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", n_reps)
  for (rep_idx in seq_len(n_reps)) {
    df <- base_rows
    trial_offset <- (rep_idx - 1L) * 3L
    rt_offset <- ((rep_idx - 1L) %% 7L) * 0.004
    onset_offset <- ((rep_idx - 1L) %% 5L) * 0.01
    df$trial <- df$pattern_trial + trial_offset
    df$rt <- df$rt + rt_offset
    is_stop_branch <- df$accumulator %in% c("stop", "change")
    df$onset[is_stop_branch] <- df$onset[is_stop_branch] + onset_offset
    rows[[rep_idx]] <- df[, c("trial", "R", "rt", "component", "accumulator", "onset")]
  }

  out <- do.call(rbind, rows)
  out$R <- factor(out$R, levels = c("S", "X"))
  out$component <- factor(out$component, levels = c("go_only", "go_stop"))
  out
}

build_beest_case <- function(shared_trigger, n_reps) {
  spec <- race_spec() |>
    add_accumulator("S", "exgauss") |>
    add_accumulator("stop", "exgauss") |>
    add_accumulator("change", "exgauss") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_component("go_only", members = "S", weight = 0.75) |>
    add_component("go_stop", members = c("S", "stop", "change"), weight = 0.25)
  if (shared_trigger) {
    spec <- spec |>
      add_trigger("tg", members = c("stop", "change"), q = 0.10,
                  param = "tg_q", draw = "shared")
  }
  spec <- spec |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()

  pkg_params <- c(
    S.mu = 0.38,
    S.sigma = 0.05,
    S.tau = 0.06,
    stop.mu = 0.08,
    stop.sigma = 0.16,
    stop.tau = 0.09,
    change.mu = 0.33,
    change.sigma = 0.09,
    change.tau = 0.05
  )
  oracle_params <- c(
    S.mu = 0.38,
    S.sigma = 0.05,
    S.tau = 0.06,
    stop.mu = 0.08,
    stop.sigma = 0.16,
    stop.tau = 0.09,
    change.mu = 0.33,
    change.sigma = 0.09,
    change.tau = 0.05,
    stop_trigger = if (shared_trigger) 0.22 else 0.0
  )
  if (shared_trigger) {
    pkg_params <- c(pkg_params, tg_q = 0.22)
  }

  data_df <- make_beest_fixture_data(n_reps)
  ctx <- build_likelihood_context(spec, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  params_df <- build_param_matrix(spec, pkg_params, n_trials = max(data_df$trial))

  eval_fn <- function() log_likelihood(ctx, data_prepped, params_df)
  ll_native <- as.numeric(eval_fn())
  if (length(ll_native) != 1L || !is.finite(ll_native)) {
    stop("Non-finite BEEST benchmark log-likelihood during setup")
  }

  oracle_fn <- make_beest_loglik(data_df, include_component_weight = FALSE)
  ll_oracle <- as.numeric(oracle_fn(oracle_params))
  if (length(ll_oracle) != 1L || !is.finite(ll_oracle)) {
    stop("Non-finite BEEST oracle benchmark log-likelihood during setup")
  }

  list(
    label = if (shared_trigger) {
      "beest_shared_blocker_shared_trigger"
    } else {
      "beest_shared_blocker_no_trigger"
    },
    eval_fn = eval_fn,
    oracle_eval_fn = function() oracle_fn(oracle_params),
    native_loglik = ll_native,
    oracle_loglik = ll_oracle
  )
}

bench_case <- function(label, fn, mode = "native_cpp") {
  wrapped_fn <- fn
  inner_reps <- choose_inner_reps(
    fn = wrapped_fn,
    target_sec = target_sample_sec,
    min_reps = min_inner_reps,
    max_reps = max_inner_reps,
    override = inner_reps_override
  )
  invisible(wrapped_fn())
  row <- run_bench(label = label, fn = wrapped_fn, n_rep = n_rep, inner_reps = inner_reps)
  row$mode <- mode
  row
}

cat("Preparing focused shared-blocker benchmark cases...\n")
model2_case <- build_example_2_case(example_trials)
beest_no_trigger_case <- build_beest_case(FALSE, fixture_reps)
beest_shared_case <- build_beest_case(TRUE, fixture_reps)

cat("Benchmarking native cases...\n")
native_df <- do.call(
  rbind,
  list(
    bench_case(model2_case$label, model2_case$eval_fn, "native_cpp"),
    bench_case(beest_no_trigger_case$label, beest_no_trigger_case$eval_fn, "native_cpp"),
    bench_case(beest_shared_case$label, beest_shared_case$eval_fn, "native_cpp"),
    bench_case("beest_oracle_shared_trigger", beest_shared_case$oracle_eval_fn, "oracle_r")
  )
)

native_df$suite <- "shared_blocker_beest"

speed_rows <- flatten_speed(native_df, "shared_blocker_beest")

speed_value <- function(label, mode) {
  native_df$median_per_eval_sec[
    native_df$label == label & native_df$mode == mode
  ]
}

comparison_rows <- data.frame(
  section = "comparison",
  suite = "shared_blocker_beest",
  label = c(
    "beest_shared_vs_no_trigger_runtime_ratio",
    "beest_shared_vs_example_2_runtime_ratio",
    "beest_shared_native_vs_oracle_abs_diff"
  ),
  metric = c(
    "median_per_eval_ratio",
    "median_per_eval_ratio",
    "total_loglik_abs_diff"
  ),
  value = c(
    speed_value("beest_shared_blocker_shared_trigger", "native_cpp") /
      speed_value("beest_shared_blocker_no_trigger", "native_cpp"),
    speed_value("beest_shared_blocker_shared_trigger", "native_cpp") /
      speed_value("example_2_stop_mixture", "native_cpp"),
    abs(beest_shared_case$native_loglik - beest_shared_case$oracle_loglik)
  ),
  mode = "",
  stringsAsFactors = FALSE
)

git_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)[1]),
  error = function(e) ""
)
run_ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

combined <- rbind(speed_rows, comparison_rows)
combined$run_timestamp_utc <- run_ts
combined$git_commit <- git_commit
combined <- combined[, c(
  "run_timestamp_utc", "git_commit",
  "section", "suite", "label", "metric", "value", "mode"
)]

utils::write.csv(combined, out_file, row.names = FALSE)

cat("\nWrote focused benchmark output to", out_file, "\n")
cat("Rows:", nrow(combined), "\n")
cat(sprintf(
  "BEEST shared/native oracle abs diff: %.6g\n",
  abs(beest_shared_case$native_loglik - beest_shared_case$oracle_loglik)
))
