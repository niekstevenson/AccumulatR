#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BENCH_SEED", "20260424")))

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
  if (identical(Sys.getenv("ACCUMULATR_BENCH_INSTALLED"), "true")) {
    library(AccumulatR)
  } else {
    library(pkgload)
    load_all(repo_root, quiet = TRUE, helpers = FALSE)
  }
})
source(file.path("dev", "examples", "new_API.R"))

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- Sys.getenv("ACCUMULATR_BENCH_OUT", "")
if (!nzchar(out_file)) {
  out_file <- file.path(out_dir, "benchmark_speed.csv")
}

n_trials <- as.integer(Sys.getenv("ACCUMULATR_BENCH_TRIALS", "50"))
if (!is.finite(n_trials) || n_trials < 1L) {
  n_trials <- 50L
}
n_rep <- as.integer(Sys.getenv("ACCUMULATR_BENCH_N_REP", "3"))
if (!is.finite(n_rep) || n_rep < 1L) {
  n_rep <- 3L
}
target_sample_sec <- as.numeric(Sys.getenv("ACCUMULATR_BENCH_TARGET_SEC", "0.05"))
if (!is.finite(target_sample_sec) || target_sample_sec <= 0) {
  target_sample_sec <- 0.05
}
min_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MIN_INNER_REPS", "1"))
if (!is.finite(min_inner_reps) || min_inner_reps < 1L) {
  min_inner_reps <- 1L
}
max_inner_reps <- as.integer(Sys.getenv("ACCUMULATR_BENCH_MAX_INNER_REPS", "200"))
if (!is.finite(max_inner_reps) || max_inner_reps < min_inner_reps) {
  max_inner_reps <- max(1000L, min_inner_reps)
}

git_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]]),
  error = function(e) ""
)
git_dirty <- tryCatch(
  length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
  error = function(e) NA
)
run_timestamp_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

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
    n_trials = n_trials,
    n_rep = n_rep,
    median_sec = stats::median(timings),
    min_sec = min(timings),
    max_sec = max(timings),
    inner_reps = inner_reps,
    median_per_eval_sec = stats::median(timings) / inner_reps,
    median_per_eval_ms = 1000.0 * stats::median(timings) / inner_reps,
    mode = "ir",
    run_timestamp_utc = run_timestamp_utc,
    git_commit = git_commit,
    git_dirty = git_dirty,
    stringsAsFactors = FALSE
  )
}

stop_change_model <- function() {
  structure <- race_spec() |>
    add_accumulator("S", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("change", "lognormal") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_component("go_only", members = "S", weight = .75) |>
    add_component("go_stop", members = c("S", "stop", "change"), weight = .25) |>
    add_trigger("stop_trigger", members = c("stop", "change"), q = 0.05) |>
    set_mixture_options(mode = "fixed") |>
    set_parameters(list(
      m_go = "S.m",
      m_stop = "stop.m",
      m_change = "change.m",
      s_go = "S.s",
      s_stop = "stop.s",
      s_change = "change.s",
      t0_go = "S.t0",
      t0_change = "change.t0",
      q = c("stop.q", "change.q")
    )) |>
    finalize_model()
  pars <- c(
    m_go = log(0.30), s_go = 0.18, t0_go = 0.00,
    m_stop = log(0.22), s_stop = 0.18,
    m_change = log(0.40), s_change = 0.18, t0_change = 0.00,
    q = 0.05, S.q = 0.00
  )
  list(structure = structure, pars = pars)
}

stim_selective_model <- function() {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("S1", "lognormal") |>
    add_accumulator("IS", "lognormal") |>
    add_accumulator("S2", "lognormal") |>
    add_outcome(
      "A",
      first_of(
        inhibit("A", by = "S1"),
        all_of("A", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "B",
      first_of(
        inhibit("B", by = "S1"),
        all_of("B", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "STOP",
      all_of("S1", inhibit("S2", by = "IS")),
      options = list(map_outcome_to = NA_character_)
    ) |>
    add_component("go", members = c("A", "B")) |>
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    set_parameters(list(
      m_go = c("A.m", "B.m"),
      s_go = c("A.s", "B.s"),
      t0_go = c("A.t0", "B.t0")
    )) |>
    finalize_model()
  pars <- c(
    m_go = log(0.30), s_go = 0.18, t0_go = 0.05,
    S1.m = log(0.26), S1.s = 0.18, S1.t0 = 0.00,
    IS.m = log(0.35), IS.s = 0.18, IS.t0 = 0.00,
    S2.m = log(0.32), S2.s = 0.18, S2.t0 = 0.00
  )
  list(structure = structure, pars = pars)
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
  }),
  stop_change_shared_trigger = stop_change_model(),
  stim_selective_stop = stim_selective_model()
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

out <- do.call(
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

utils::write.csv(out, out_file, row.names = FALSE)

cat("Wrote benchmark speed output to", out_file, "\n")
cat("Rows:", nrow(out), "\n")
