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

trial_loglik_vector <- function(ctx, params_df, min_ll = log(1e-10)) {
  n_trials <- ctx$n_trials %||% length(unique(ctx$data_df$trial))
  ok <- rep(TRUE, n_trials)
  rel_tol <- ctx$rel_tol %||% 1e-5
  abs_tol <- ctx$abs_tol %||% 1e-6
  max_depth <- ctx$max_depth %||% 12L
  pm <- as.matrix(params_df)
  out <- numeric(n_trials)
  for (i in seq_len(n_trials)) {
    out[i] <- AccumulatR:::cpp_loglik(
      ctx$native_ctx,
      pm,
      ctx$data_df_cpp,
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
  spec <- race_spec() |>
    add_accumulator("plain", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("PLAIN", "plain") |>
    add_outcome("GUARD", inhibit("a", by = inhibit("b", by = inhibit("c", by = "d")))) |>
    finalize_model()
  params_df <- build_param_matrix(
    spec,
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
    R = rep("PLAIN", n_trials),
    rt = stats::rlnorm(n_trials, meanlog = log(0.48), sdlog = 0.12),
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(spec, data_df)
  list(label = "depth3_guard_competitor", ctx = ctx, params_df = params_df)
}

build_shared_nway_trigger_case <- function(n_trials = 40L) {
  spec <- race_spec() |>
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
  params_df <- build_param_matrix(
    spec,
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
    R = rep(NA_character_, n_trials),
    rt = rep(NA_real_, n_trials),
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(spec, data_df)
  list(label = "shared_gate_nway_shared_triggers", ctx = ctx, params_df = params_df)
}

cases <- list(
  build_depth3_case(),
  build_shared_nway_trigger_case()
)

golden <- lapply(cases, function(cs) {
  total <- as.numeric(log_likelihood(cs$ctx, cs$params_df))
  per_trial <- trial_loglik_vector(cs$ctx, cs$params_df)
  list(
    label = cs$label,
    n_trials = length(per_trial),
    total_loglik = total,
    trial_loglik = per_trial
  )
})

out <- list(
  created_at = as.character(Sys.time()),
  seed = 20260218L,
  tolerance = 1e-4,
  cases = golden
)

out_file <- file.path("dev", "scripts", "scratch_outputs", "loglik_golden_v1.rds")
saveRDS(out, out_file)
cat("Wrote golden fixture to", out_file, "\n")
