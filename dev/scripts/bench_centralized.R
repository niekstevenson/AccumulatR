#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BENCH_SEED", "20260220")))

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

bench_lib <- tempfile("accumulatr_bench_lib_")
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

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, "benchmark_centralized.csv")

golden_script <- file.path("dev", "scripts", "check_loglik_golden.R")

n_trials <- as.integer(Sys.getenv("ACCUMULATR_BENCH_TRIALS", "250"))
if (!is.finite(n_trials) || n_trials < 50L) {
  n_trials <- 250L
}
n_rep <- as.integer(Sys.getenv("ACCUMULATR_BENCH_N_REP", "4"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 4L
}
inner_reps_override <- as.integer(Sys.getenv("ACCUMULATR_BENCH_INNER_REPS", NA_character_))
if (is.na(inner_reps_override) || !is.finite(inner_reps_override) || inner_reps_override < 1L) {
  inner_reps_override <- NA_integer_
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

run_script <- function(path) {
  out <- system2("Rscript", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status") %||% 0L
  list(status = as.integer(status), output = out)
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

parse_golden_output <- function(lines) {
  case_rows <- list()
  max_abs_diff <- NA_real_
  tolerance <- NA_real_
  passed <- FALSE

  for (ln in lines) {
    case_m <- regexec(
      "^([^:]+): total_diff=([0-9.eE+-]+) trial_max_diff=([0-9.eE+-]+)$",
      ln
    )
    case_hit <- regmatches(ln, case_m)[[1]]
    if (length(case_hit) == 4L) {
      case_rows[[length(case_rows) + 1L]] <- data.frame(
        section = "correctness",
        suite = "golden",
        label = case_hit[[2]],
        metric = c("total_diff", "trial_max_diff"),
        value = c(as.numeric(case_hit[[3]]), as.numeric(case_hit[[4]])),
        mode = "",
        stringsAsFactors = FALSE
      )
      next
    }

    overall_m <- regexec(
      "^max_abs_diff=([0-9.eE+-]+) tolerance=([0-9.eE+-]+)$",
      ln
    )
    overall_hit <- regmatches(ln, overall_m)[[1]]
    if (length(overall_hit) == 3L) {
      max_abs_diff <- as.numeric(overall_hit[[2]])
      tolerance <- as.numeric(overall_hit[[3]])
      next
    }

    if (identical(trimws(ln), "Golden check passed")) {
      passed <- TRUE
    }
  }

  rows <- do.call(rbind, case_rows)
  overall_rows <- data.frame(
    section = "correctness",
    suite = "golden",
    label = "overall",
    metric = c("max_abs_diff", "tolerance", "passed"),
    value = c(max_abs_diff, tolerance, if (passed) 1.0 else 0.0),
    mode = "",
    stringsAsFactors = FALSE
  )

  if (is.null(rows)) {
    overall_rows
  } else {
    rbind(rows, overall_rows)
  }
}

build_depth3_guard_spec <- function() {
  race_spec() |>
    add_accumulator("plain", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("PLAIN", "plain") |>
    add_outcome("GUARD", inhibit("a", by = inhibit("b", by = inhibit("c", by = "d")))) |>
    finalize_model()
}
params_depth3_guard <- c(
  plain.meanlog = log(0.42), plain.sdlog = 0.18,
  a.meanlog = log(0.34), a.sdlog = 0.20,
  b.meanlog = log(0.30), b.sdlog = 0.20,
  c.meanlog = log(0.28), c.sdlog = 0.18,
  d.meanlog = log(0.26), d.sdlog = 0.16
)

build_shared_gate_pair_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP", all_of("x2", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
}
params_shared_gate_pair <- c(
  x1.meanlog = log(0.32), x1.sdlog = 0.18,
  x2.meanlog = log(0.36), x2.sdlog = 0.18,
  gate.meanlog = log(0.24), gate.sdlog = 0.14
)

build_shared_gate_nway_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
}
params_shared_gate_nway <- c(
  x1.meanlog = log(0.31), x1.sdlog = 0.18,
  x2.meanlog = log(0.34), x2.sdlog = 0.18,
  x3.meanlog = log(0.37), x3.sdlog = 0.18,
  gate.meanlog = log(0.23), gate.sdlog = 0.14
)

build_shared_gate_nway_trigger_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
    add_trigger("tg_x12", members = c("x1", "x2"), q = 0.10, param = "q_x12", draw = "shared") |>
    add_trigger("tg_x3g", members = c("x3", "gate"), q = 0.16, param = "q_x3g", draw = "shared") |>
    finalize_model()
}
params_shared_gate_nway_trigger <- c(
  x1.meanlog = log(0.31), x1.sdlog = 0.18,
  x2.meanlog = log(0.34), x2.sdlog = 0.18,
  x3.meanlog = log(0.37), x3.sdlog = 0.18,
  gate.meanlog = log(0.23), gate.sdlog = 0.14,
  q_x12 = 0.10,
  q_x3g = 0.16
)

make_nested_case <- function(spec_builder, params, observed_label, rt_mean, n_trials = 6000L, inner_reps = 8L) {
  spec <- spec_builder()
  outcome_levels <- names((spec$prep %||% list())$outcomes %||% list())
  outcome_levels <- as.character(outcome_levels)
  outcome_levels <- outcome_levels[!is.na(outcome_levels) & nzchar(outcome_levels)]
  outcome_levels <- unique(outcome_levels)
  if (length(outcome_levels) == 0L) {
    stop("Could not determine outcome level order from model structure")
  }
  params_df <- build_param_matrix(spec, params, n_trials = n_trials)
  rt_col <- if (is.na(rt_mean)) {
    rep(NA_real_, n_trials)
  } else {
    stats::rlnorm(n_trials, meanlog = log(rt_mean), sdlog = 0.12)
  }
  data_df <- data.frame(
    trial = seq_len(n_trials),
    R = factor(rep(observed_label, n_trials), levels = outcome_levels),
    rt = rt_col,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(spec, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  ll <- as.numeric(log_likelihood(ctx, data_prepped, params_df))
  if (!is.finite(ll)) {
    stop("Non-finite log-likelihood for nested case: ", observed_label)
  }
  list(
    ctx = ctx,
    data_prepped = data_prepped,
    params_df = params_df,
    inner_reps = inner_reps
  )
}

nested_cases <- list(
  depth3_guard_competitor = make_nested_case(
    build_depth3_guard_spec,
    params_depth3_guard,
    observed_label = "PLAIN",
    rt_mean = 0.48,
    inner_reps = 8L
  ),
  shared_gate_pair = make_nested_case(
    build_shared_gate_pair_spec,
    params_shared_gate_pair,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  ),
  shared_gate_nway = make_nested_case(
    build_shared_gate_nway_spec,
    params_shared_gate_nway,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  ),
  shared_gate_nway_shared_triggers = make_nested_case(
    build_shared_gate_nway_trigger_spec,
    params_shared_gate_nway_trigger,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  )
)

run_nested_suite <- function(cases) {
  do.call(
    rbind,
    lapply(names(cases), function(case_name) {
      case <- cases[[case_name]]
      invisible(log_likelihood(case$ctx, case$data_prepped, case$params_df))
      row <- run_bench(
        label = case_name,
        fn = function() log_likelihood(case$ctx, case$data_prepped, case$params_df),
        n_rep = 6L,
        inner_reps = case$inner_reps
      )
      row$mode <- "ir"
      row
    })
  )
}

example_1_simple <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_outcome("R1", "go1") |>
    add_outcome("R2", "go2") |>
    finalize_model()
}
params_example_1_simple <- c(
  go1.meanlog = log(0.30),
  go1.sdlog = 0.18,
  go2.meanlog = log(0.32),
  go2.sdlog = 0.18
)

example_2_stop_mixture <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("stop", "exgauss", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_outcome("R1", inhibit("go1", by = "stop")) |>
    add_outcome("R2", all_of("go2", "stop")) |>
    add_component("go_only", members = c("go1"), attrs = list(component = "go_only"), weight = 0.5) |>
    add_component("go_stop", members = c("go1", "stop", "go2"), attrs = list(component = "go_stop"), weight = 0.5) |>
    set_mixture_options(mode = "fixed") |>
    finalize_model()
}
params_example_2_stop_mixture <- c(
  go1.meanlog = log(0.35),
  go1.sdlog = 0.2,
  stop.mu = 0.1,
  stop.sigma = 0.04,
  stop.tau = 0.1,
  go2.meanlog = log(0.60),
  go2.sdlog = 0.18
)

example_3_stop_na <- function() {
  race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_accumulator("stop", "lognormal", onset = 0.15) |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_outcome("STOP", "stop", options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
}
params_example_3_stop_na <- c(
  go_left.meanlog = log(0.30),
  go_left.sdlog = 0.20,
  go_right.meanlog = log(0.32),
  go_right.sdlog = 0.20,
  stop.meanlog = log(0.15),
  stop.sdlog = 0.18
)

example_5_timeout_guess <- function() {
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

example_6_dual_path <- function() {
  race_spec() |>
    add_accumulator("acc_taskA", "lognormal") |>
    add_accumulator("acc_taskB", "lognormal") |>
    add_accumulator("acc_gateC", "lognormal") |>
    add_outcome("Outcome_via_A", all_of("acc_taskA", "acc_gateC")) |>
    add_outcome("Outcome_via_B", all_of("acc_taskB", "acc_gateC")) |>
    finalize_model()
}
params_example_6_dual_path <- c(
  acc_taskA.meanlog = log(0.28),
  acc_taskA.sdlog = 0.18,
  acc_taskB.meanlog = log(0.32),
  acc_taskB.sdlog = 0.18,
  acc_gateC.meanlog = log(0.30),
  acc_gateC.sdlog = 0.18
)

example_7_mixture <- function() {
  race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("TARGET", c("target_fast", "target_slow")) |>
    add_outcome("R1", "TARGET") |>
    add_outcome("R2", "competitor") |>
    add_component("fast", members = c("target_fast", "competitor"), weight_param = "p_fast") |>
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

example_10_exclusion <- function() {
  race_spec() |>
    add_accumulator("R1_acc", "lognormal") |>
    add_accumulator("R2_acc", "lognormal") |>
    add_accumulator("X_acc", "lognormal") |>
    add_outcome("R1", inhibit("R1_acc", by = "X_acc")) |>
    add_outcome("R2", "R2_acc") |>
    finalize_model()
}
params_example_10_exclusion <- c(
  R1_acc.meanlog = log(0.35),
  R1_acc.sdlog = 0.18,
  R2_acc.meanlog = log(0.45),
  R2_acc.sdlog = 0.18,
  X_acc.meanlog = log(0.35),
  X_acc.sdlog = 0.18
)

example_16_guard_tie_simple <- function() {
  race_spec() |>
    add_accumulator("go_fast", "lognormal") |>
    add_accumulator("go_slow", "lognormal") |>
    add_accumulator("gate_shared", "lognormal") |>
    add_accumulator("stop_control", "lognormal") |>
    add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
    add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
    finalize_model()
}
params_example_16_guard_tie_simple <- c(
  go_fast.meanlog = log(0.28),
  go_fast.sdlog = 0.18,
  go_slow.meanlog = log(0.34),
  go_slow.sdlog = 0.18,
  gate_shared.meanlog = log(0.30),
  gate_shared.sdlog = 0.16,
  stop_control.meanlog = log(0.27),
  stop_control.sdlog = 0.15
)

example_21_simple_q <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_outcome("R1", "go1") |>
    add_outcome("R2", "go2") |>
    add_trigger(
      "shared_trigger",
      members = c("go1", "go2"),
      q = 0.10,
      param = "go_trigger_q",
      draw = "shared"
    ) |>
    finalize_model()
}
params_example_21_simple_q <- c(
  go1.meanlog = log(0.30),
  go1.sdlog = 0.18,
  go2.meanlog = log(0.32),
  go2.sdlog = 0.18,
  go_trigger_q = 0.10
)

example_22_shared_q <- function() {
  race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_trigger(
      "q_shared",
      members = c("go_left", "go_right"),
      q = 0.10,
      param = "q_shared",
      draw = "independent"
    ) |>
    finalize_model()
}
params_example_22_shared_q <- c(
  go_left.meanlog = log(0.30),
  go_left.sdlog = 0.18,
  go_right.meanlog = log(0.32),
  go_right.sdlog = 0.18,
  q_shared = 0.10
)

example_23_ranked_chain <- function() {
  race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    finalize_model()
}
params_example_23_ranked_chain <- c(
  a.meanlog = log(0.30),
  a.sdlog = 0.16,
  b.meanlog = log(0.22),
  b.sdlog = 0.16
)

example_models <- list(
  example_1_simple = list(spec = example_1_simple, params = params_example_1_simple),
  example_2_stop_mixture = list(spec = example_2_stop_mixture, params = params_example_2_stop_mixture),
  example_3_stop_na = list(spec = example_3_stop_na, params = params_example_3_stop_na),
  example_5_timeout_guess = list(spec = example_5_timeout_guess, params = params_example_5_timeout_guess),
  example_6_dual_path = list(spec = example_6_dual_path, params = params_example_6_dual_path),
  example_7_mixture = list(spec = example_7_mixture, params = params_example_7_mixture),
  example_10_exclusion = list(spec = example_10_exclusion, params = params_example_10_exclusion),
  example_16_guard_tie_simple = list(spec = example_16_guard_tie_simple, params = params_example_16_guard_tie_simple),
  example_21_simple_q = list(spec = example_21_simple_q, params = params_example_21_simple_q),
  example_22_shared_q = list(spec = example_22_shared_q, params = params_example_22_shared_q),
  example_23_ranked_chain = list(spec = example_23_ranked_chain, params = params_example_23_ranked_chain)
)

model_filter <- trimws(Sys.getenv("ACCUMULATR_BENCH_MODELS", ""))
if (nzchar(model_filter)) {
  keep <- trimws(strsplit(model_filter, ",", fixed = TRUE)[[1]])
  keep <- keep[nzchar(keep)]
  example_models <- example_models[intersect(keep, names(example_models))]
}
if (length(example_models) == 0L) {
  stop("No models selected for benchmark")
}

make_example_case <- function(mod, n_trials) {
  spec_obj <- mod$spec()
  structure <- finalize_model(spec_obj)
  params_vec <- mod$params

  params_df <- build_param_matrix(spec_obj, params_vec, n_trials = n_trials)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  ctx <- build_likelihood_context(structure, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  params_df_slim <- build_param_matrix(
    spec_obj,
    params_vec,
    n_trials = max(data_df$trial),
    layout = ctx$param_layout
  )

  ll <- as.numeric(log_likelihood(ctx, data_prepped, params_df_slim))
  if (length(ll) != 1L || !is.finite(ll)) {
    stop("Non-finite log-likelihood during setup")
  }

  list(
    ctx = ctx,
    data_prepped = data_prepped,
    params_df = params_df_slim
  )
}

run_examples_suite <- function(models, n_trials, n_rep, target_sec, min_reps, max_reps, inner_override) {
  cases <- lapply(models, make_example_case, n_trials = n_trials)
  case_inner_reps <- setNames(integer(length(cases)), names(cases))
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    eval_fn <- function() {
      log_likelihood(case$ctx, case$data_prepped, case$params_df)
    }
    case_inner_reps[[case_name]] <- choose_inner_reps(
      fn = eval_fn,
      target_sec = target_sec,
      min_reps = min_reps,
      max_reps = max_reps,
      override = inner_override
    )
  }

  bench_rows <- vector("list", length(cases))
  case_names <- names(cases)
  for (i in seq_along(case_names)) {
    case_name <- case_names[[i]]
    case <- cases[[case_name]]
    eval_fn <- function() {
      log_likelihood(case$ctx, case$data_prepped, case$params_df)
    }
    invisible(eval_fn())
    row <- run_bench(
      label = case_name,
      fn = eval_fn,
      n_rep = n_rep,
      inner_reps = case_inner_reps[[case_name]]
    )
    row$mode <- "ir"
    bench_rows[[i]] <- row
  }
  do.call(rbind, bench_rows)
}

cat("Running nested suite...\n")
nested_df <- run_nested_suite(nested_cases)

cat("Running examples suite...\n")
examples_df <- run_examples_suite(
  models = example_models,
  n_trials = n_trials,
  n_rep = n_rep,
  target_sec = target_sample_sec,
  min_reps = min_inner_reps,
  max_reps = max_inner_reps,
  inner_override = inner_reps_override
)

cat("Running golden correctness check...\n")
golden_run <- run_script(golden_script)
if (golden_run$status != 0L) {
  cat(paste(golden_run$output, collapse = "\n"), "\n")
  stop("Golden correctness check failed")
}

speed_rows <- rbind(
  flatten_speed(nested_df, "nested"),
  flatten_speed(examples_df, "examples_models")
)
correctness_rows <- parse_golden_output(golden_run$output)

git_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)[1]),
  error = function(e) ""
)
run_ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

combined <- rbind(speed_rows, correctness_rows)
combined$run_timestamp_utc <- run_ts
combined$git_commit <- git_commit
combined <- combined[, c(
  "run_timestamp_utc", "git_commit",
  "section", "suite", "label", "metric", "value", "mode"
)]

utils::write.csv(combined, out_file, row.names = FALSE)

cat("\nWrote canonical benchmark output to", out_file, "\n")
cat("Rows:", nrow(combined), "\n")
