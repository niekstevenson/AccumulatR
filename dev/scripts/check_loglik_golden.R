#!/usr/bin/env Rscript

set.seed(20260218L)
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(AccumulatR))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

trial_loglik_vector <- function(ctx, eval_inputs, params_df, min_ll = log(1e-10)) {
  n_trials <- eval_inputs$n_trials
  ok <- rep(TRUE, n_trials)
  rel_tol <- ctx$rel_tol %||% 1e-5
  abs_tol <- ctx$abs_tol %||% 1e-6
  max_depth <- ctx$max_depth %||% 12L
  pm <- AccumulatR:::.align_param_matrix_to_layout(
    as.matrix(params_df),
    eval_inputs$param_layout
  )
  out <- numeric(n_trials)
  for (i in seq_len(n_trials)) {
    out[i] <- AccumulatR:::cpp_loglik(
      ctx$native_ctx,
      pm,
      eval_inputs$data,
      ok,
      as.integer(i),
      min_ll,
      rel_tol,
      abs_tol,
      max_depth
    )
  }
  out
}

build_depth3_case <- function(n_trials = 40L) {
  structure <- race_spec() |>
    add_accumulator("plain", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("PLAIN", "plain") |>
    add_outcome("GUARD", inhibit("a", by = inhibit("b", by = inhibit("c", by = "d")))) |>
    finalize_model()
  outcome_levels <- names(structure$prep$outcomes)
  params_df <- build_param_matrix(
    structure,
    c(
      plain.meanlog = log(0.42), plain.sdlog = 0.18,
      a.meanlog = log(0.34), a.sdlog = 0.20,
      b.meanlog = log(0.30), b.sdlog = 0.20,
      c.meanlog = log(0.28), c.sdlog = 0.18,
      d.meanlog = log(0.26), d.sdlog = 0.16
    ),
    n_trials = n_trials
  )
  data_df <- data.frame(
    trial = seq_len(n_trials),
    R = factor(rep("PLAIN", n_trials), levels = outcome_levels),
    rt = stats::rlnorm(n_trials, meanlog = log(0.48), sdlog = 0.12)
  )
  ctx <- build_likelihood_context(structure, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  list(
    label = "depth3_guard_competitor",
    ctx = ctx,
    data_df = data_df,
    data_prepped = data_prepped,
    params_df = params_df
  )
}

build_shared_nway_trigger_case <- function(n_trials = 40L) {
  structure <- race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"),
                options = list(map_outcome_to = NA_character_)) |>
    add_trigger("tg_x12", members = c("x1", "x2"), q = 0.10, param = "q_x12",
                draw = "shared") |>
    add_trigger("tg_x3g", members = c("x3", "gate"), q = 0.16, param = "q_x3g",
                draw = "shared") |>
    finalize_model()
  outcome_levels <- names(structure$prep$outcomes)
  params_df <- build_param_matrix(
    structure,
    c(
      x1.meanlog = log(0.31), x1.sdlog = 0.18,
      x2.meanlog = log(0.34), x2.sdlog = 0.18,
      x3.meanlog = log(0.37), x3.sdlog = 0.18,
      gate.meanlog = log(0.23), gate.sdlog = 0.14,
      q_x12 = 0.10,
      q_x3g = 0.16
    ),
    n_trials = n_trials
  )
  data_df <- data.frame(
    trial = seq_len(n_trials),
    R = factor(rep(NA_character_, n_trials), levels = outcome_levels),
    rt = rep(NA_real_, n_trials)
  )
  ctx <- build_likelihood_context(structure, data_df)
  data_prepped <- prepare_likelihood_data(ctx, data_df)
  list(
    label = "shared_gate_nway_shared_triggers",
    ctx = ctx,
    data_df = data_df,
    data_prepped = data_prepped,
    params_df = params_df
  )
}

fixture_path <- file.path("dev", "scripts", "scratch_outputs", "loglik_golden_v1.rds")
if (!file.exists(fixture_path)) {
  stop("Golden fixture missing: ", fixture_path)
}
golden <- readRDS(fixture_path)
tolerance <- golden$tolerance %||% 1e-4

cases <- list(
  build_depth3_case(),
  build_shared_nway_trigger_case()
)

case_map <- setNames(cases, vapply(cases, `[[`, character(1), "label"))
golden_map <- setNames(golden$cases, vapply(golden$cases, `[[`, character(1), "label"))

max_diff <- 0.0
for (label in names(case_map)) {
  if (!label %in% names(golden_map)) {
    stop("Golden fixture is missing case: ", label)
  }
  cs <- case_map[[label]]
  gl <- golden_map[[label]]
  total <- as.numeric(log_likelihood(cs$ctx, cs$data_prepped, cs$params_df))
  trial <- trial_loglik_vector(cs$ctx, cs$data_prepped, cs$params_df)
  if (length(trial) != gl$n_trials) {
    stop("Case ", label, " trial count mismatch: ", length(trial), " vs ", gl$n_trials)
  }
  diff_total <- abs(total - gl$total_loglik)
  diff_trial <- max(abs(trial - gl$trial_loglik))
  case_max <- max(diff_total, diff_trial)
  max_diff <- max(max_diff, case_max)
  cat(sprintf(
    "%s: total_diff=%.6g trial_max_diff=%.6g\n",
    label, diff_total, diff_trial
  ))
}

cat(sprintf("max_abs_diff=%.6g tolerance=%.6g\n", max_diff, tolerance))
if (!is.finite(max_diff) || max_diff > tolerance) {
  stop("Golden check failed")
}
cat("Golden check passed\n")
