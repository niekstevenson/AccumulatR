#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(AccumulatR))

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

resolve_competitor_ids <- function(prep, outcome_label, outcome_idx) {
  competitor_map <- AccumulatR:::.prep_competitors(prep) %||% list()
  use_indexed <- isTRUE(attr(competitor_map, "by_index")) &&
    length(competitor_map) == length(prep$outcomes)
  comp_exprs <- if (use_indexed) {
    competitor_map[[outcome_idx]] %||% list()
  } else {
    competitor_map[[outcome_label]] %||% list()
  }
  if (length(comp_exprs) == 0L) {
    return(integer(0))
  }
  comp_nodes <- lapply(comp_exprs, function(expr) AccumulatR:::.expr_lookup_compiled(expr, prep))
  if (any(vapply(comp_nodes, is.null, logical(1)))) {
    return(integer(0))
  }
  as.integer(vapply(comp_nodes, function(node) node$id %||% NA_integer_, integer(1)))
}

native_prob <- function(prep, native_ctx, params_df, outcome_label, upper) {
  labels <- vapply(prep$outcomes, function(x) as.character(x$label %||% ""), character(1))
  idx <- match(outcome_label, labels)
  if (is.na(idx)) {
    stop("Outcome not found in prep: ", outcome_label)
  }
  compiled <- AccumulatR:::.expr_lookup_compiled(prep$outcomes[[idx]]$expr, prep)
  if (is.null(compiled) || is.null(compiled$id)) {
    stop("Missing compiled node id for outcome: ", outcome_label)
  }
  comp_ids <- resolve_competitor_ids(prep, outcome_label, idx)
  AccumulatR:::native_outcome_probability_params_cpp_idx(
    native_ctx,
    as.integer(compiled$id),
    upper,
    0L,
    integer(0),
    integer(0),
    comp_ids,
    AccumulatR:::.integrate_rel_tol(),
    AccumulatR:::.integrate_abs_tol(),
    params_df
  )
}

build_shared_gate_pair <- function() {
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

build_shared_gate_nway <- function() {
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

check_case <- function(structure, params, outcome_label) {
  params_mat <- build_param_matrix(structure$model_spec, params, n_trials = 1L)
  params_rows <- AccumulatR:::.param_matrix_to_rows(structure, params_mat)
  trial_rows <- params_rows[params_rows$trial == 1L, , drop = FALSE]
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  native_ctx <- AccumulatR:::native_context_build(prep)

  p4 <- as.numeric(native_prob(prep, native_ctx, trial_rows, outcome_label, 4.0))
  p8 <- as.numeric(native_prob(prep, native_ctx, trial_rows, outcome_label, 8.0))
  p_inf <- as.numeric(native_prob(prep, native_ctx, trial_rows, outcome_label, Inf))

  cat(sprintf(
    "%s: p(4)=%.9g p(8)=%.9g p(inf)=%.9g\n",
    outcome_label, p4, p8, p_inf
  ))

  if (!all(is.finite(c(p4, p8, p_inf)))) {
    stop("Non-finite probabilities for outcome ", outcome_label)
  }
  monotonic_tol <- 5e-4
  if (!((p4 <= p8 + monotonic_tol) && (p8 <= p_inf + monotonic_tol))) {
    stop("Monotonicity failed for outcome ", outcome_label)
  }
  if (abs(p_inf - p8) > 2e-3) {
    stop("Tail approximation too loose for outcome ", outcome_label, ": p(inf)-p(8)=", p_inf - p8)
  }
}

pair_structure <- build_shared_gate_pair()
check_case(pair_structure, params_shared_gate_pair, "RESP")

nway_structure <- build_shared_gate_nway()
check_case(nway_structure, params_shared_gate_nway, "RESP2")
check_case(nway_structure, params_shared_gate_nway, "RESP3")

cat("Step 5 infinite-tail checks passed\n")
