#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_PROFILE_SEED", "20260213")))
library(AccumulatR)

run_sec <- as.numeric(Sys.getenv("ACCUMULATR_PROFILE_RUN_SEC", "45"))
if (!is.finite(run_sec) || run_sec <= 0) {
  run_sec <- 45
}
n_trials <- as.integer(Sys.getenv("ACCUMULATR_PROFILE_TRIALS", "6000"))
if (!is.finite(n_trials) || n_trials < 100) {
  n_trials <- 6000L
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

make_case <- function(spec_builder, params, observed_label, rt_mean) {
  spec <- spec_builder()
  params_df <- build_param_matrix(spec, params, n_trials = n_trials)
  rt_col <- if (is.na(rt_mean)) rep(NA_real_, n_trials) else stats::rlnorm(n_trials, meanlog = log(rt_mean), sdlog = 0.12)
  data_df <- data.frame(
    trial = seq_len(n_trials),
    R = rep(observed_label, n_trials),
    rt = rt_col,
    stringsAsFactors = FALSE
  )
  ctx <- build_likelihood_context(spec, data_df)
  list(ctx = ctx, params_df = params_df)
}

cases <- list(
  depth3_guard_competitor = make_case(build_depth3_guard_spec, params_depth3_guard, "PLAIN", 0.48),
  shared_gate_pair = make_case(build_shared_gate_pair_spec, params_shared_gate_pair, NA_character_, NA_real_),
  shared_gate_nway = make_case(build_shared_gate_nway_spec, params_shared_gate_nway, NA_character_, NA_real_),
  shared_gate_nway_shared_triggers = make_case(build_shared_gate_nway_trigger_spec, params_shared_gate_nway_trigger, NA_character_, NA_real_)
)

case_filter <- trimws(Sys.getenv("ACCUMULATR_PROFILE_CASES", ""))
if (nzchar(case_filter)) {
  selected <- strsplit(case_filter, ",", fixed = TRUE)[[1]]
  selected <- trimws(selected)
  selected <- selected[nzchar(selected)]
  selected <- intersect(selected, names(cases))
  if (length(selected) == 0L) {
    stop("ACCUMULATR_PROFILE_CASES did not match any known cases")
  }
  cases <- cases[selected]
}

# Warmup JIT/lookup paths before sampling window.
invisible(lapply(cases, function(x) {
  log_likelihood(x$ctx, x$params_df)
}))

start <- Sys.time()
checksum <- 0.0
iters <- integer(length(cases))
names(iters) <- names(cases)

while (as.numeric(difftime(Sys.time(), start, units = "secs")) < run_sec) {
  for (nm in names(cases)) {
    val <- log_likelihood(cases[[nm]]$ctx, cases[[nm]]$params_df)
    checksum <- checksum + sum(val)
    iters[[nm]] <- iters[[nm]] + 1L
  }
}

cat("profile_workload_nested complete\n")
print(iters)
cat("checksum:", format(checksum, digits = 16), "\n")
