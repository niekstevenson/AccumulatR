#!/usr/bin/env Rscript

set.seed(20260212)

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
bench_lib <- tempfile("accumulatr_bench_lib_")
dir.create(bench_lib, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(bench_lib, recursive = TRUE, force = TRUE), add = TRUE)
install_out <- system2(
  "R",
  args = c(
    "CMD", "INSTALL", "--preclean", "--no-multiarch",
    paste0("--library=", bench_lib), repo_root
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

run_bench <- function(label, fn, n_rep = 6L, inner_reps = 8L) {
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

make_case <- function(spec_builder, params, observed_label, rt_mean, n_trials = 6000L, inner_reps = 8L) {
  spec <- spec_builder()
  params_df <- build_param_matrix(spec, params, n_trials = n_trials)
  rt_col <- if (is.na(rt_mean)) {
    rep(NA_real_, n_trials)
  } else {
    stats::rlnorm(n_trials, meanlog = log(rt_mean), sdlog = 0.12)
  }
  data_df <- data.frame(
    trial = seq_len(n_trials),
    R = rep(observed_label, n_trials),
    rt = rt_col,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(spec, data_df)
  ll <- as.numeric(log_likelihood(ctx, params_df))
  if (!is.finite(ll)) {
    stop("Non-finite log-likelihood for case: ", observed_label)
  }
  list(ctx = ctx, params_df = params_df, inner_reps = inner_reps)
}

cases <- list(
  depth3_guard_competitor = make_case(
    build_depth3_guard_spec,
    params_depth3_guard,
    observed_label = "PLAIN",
    rt_mean = 0.48,
    inner_reps = 8L
  ),
  shared_gate_pair = make_case(
    build_shared_gate_pair_spec,
    params_shared_gate_pair,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  ),
  shared_gate_nway = make_case(
    build_shared_gate_nway_spec,
    params_shared_gate_nway,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  ),
  shared_gate_nway_shared_triggers = make_case(
    build_shared_gate_nway_trigger_spec,
    params_shared_gate_nway_trigger,
    observed_label = NA_character_,
    rt_mean = NA_real_,
    inner_reps = 80L
  )
)

bench <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    case <- cases[[case_name]]
    # Warmup run.
    invisible(log_likelihood(case$ctx, case$params_df))
    row <- run_bench(
      case_name,
      function() log_likelihood(case$ctx, case$params_df),
      n_rep = 6L,
      inner_reps = case$inner_reps
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
out_file <- file.path(out_dir, "bench_nested_integral_baseline.csv")
utils::write.csv(bench, out_file, row.names = FALSE)
cat("\nWrote baseline benchmark to", out_file, "\n")
