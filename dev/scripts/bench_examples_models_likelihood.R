#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BENCH_SEED", "20260217")))

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

library(AccumulatR)

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

# Model definitions from tests/testthat/test_examples_models.R
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
    add_trigger("shared_trigger",
      members = c("go1", "go2"),
      q = 0.10, param = "go_trigger_q", draw = "shared"
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
    add_trigger("q_shared",
      members = c("go_left", "go_right"),
      q = 0.10, param = "q_shared", draw = "independent"
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

models <- list(
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
  models <- models[intersect(keep, names(models))]
}
if (length(models) == 0L) {
  stop("No models selected for benchmark")
}

make_case <- function(mod, n_trials) {
  spec_obj <- mod$spec()
  structure <- finalize_model(spec_obj)
  params_vec <- mod$params

  params_df <- build_param_matrix(spec_obj, params_vec, n_trials = n_trials)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  ctx <- build_likelihood_context(structure, data_df)
  params_df_slim <- build_param_matrix(
    spec_obj,
    params_vec,
    n_trials = max(data_df$trial),
    layout = ctx$param_layout
  )

  ll <- as.numeric(log_likelihood(ctx, params_df_slim))
  if (length(ll) != 1L || !is.finite(ll)) {
    stop("Non-finite log-likelihood during setup")
  }

  list(ctx = ctx, params_df = params_df_slim)
}

cases <- lapply(models, make_case, n_trials = n_trials)
case_inner_reps <- setNames(integer(length(cases)), names(cases))
for (case_name in names(cases)) {
  case <- cases[[case_name]]
  eval_fn <- function() {
    log_likelihood(case$ctx, case$params_df)
  }
  case_inner_reps[[case_name]] <- choose_inner_reps(
    fn = eval_fn,
    target_sec = target_sample_sec,
    min_reps = min_inner_reps,
    max_reps = max_inner_reps,
    override = inner_reps_override
  )
}

bench <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    case <- cases[[case_name]]
    eval_fn <- function() {
      log_likelihood(case$ctx, case$params_df)
    }
    invisible(eval_fn())
    local_inner_reps <- case_inner_reps[[case_name]]
    row <- run_bench(
      label = case_name,
      fn = eval_fn,
      n_rep = n_rep,
      inner_reps = local_inner_reps
    )
    row$mode <- "ir"
    row
  })
)

print(bench, row.names = FALSE)

out_dir <- "dev/scripts/scratch_outputs"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
utils::write.csv(bench, file.path(out_dir, "bench_examples_models_likelihood.csv"), row.names = FALSE)
cat("\nWrote benchmark outputs to", out_dir, "\n")
