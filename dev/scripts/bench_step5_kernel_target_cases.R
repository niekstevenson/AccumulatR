#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(AccumulatR))

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
all_of <- AccumulatR::all_of
inhibit <- AccumulatR::inhibit
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix

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
  comp_nodes <- lapply(comp_exprs, function(expr) {
    AccumulatR:::.expr_lookup_compiled(expr, prep)
  })
  if (any(vapply(comp_nodes, is.null, logical(1)))) {
    stop("Failed to resolve competitor nodes for outcome ", outcome_label)
  }
  as.integer(vapply(
    comp_nodes,
    function(node) node$id %||% NA_integer_,
    integer(1)
  ))
}

resolve_outcome_setup <- function(structure, params, target_label, forced_label,
                                  upper, forced_mode = c("complete", "survive"),
                                  resolve_competitors = TRUE) {
  forced_mode <- match.arg(forced_mode)
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  native_ctx <- AccumulatR:::native_context_build(prep)
  params_mat <- build_param_matrix(structure$model_spec, params, n_trials = 1L)
  params_rows <- AccumulatR:::.param_matrix_to_rows(structure, params_mat)
  trial_rows <- params_rows[params_rows$trial == 1L, , drop = FALSE]

  outcome_labels <- vapply(
    prep$outcomes,
    function(x) as.character(x$label %||% ""),
    character(1)
  )
  target_idx <- match(target_label, outcome_labels)
  if (is.na(target_idx)) {
    stop("Target outcome not found: ", target_label)
  }
  compiled <- AccumulatR:::.expr_lookup_compiled(prep$outcomes[[target_idx]]$expr, prep)
  if (is.null(compiled) || is.null(compiled$id)) {
    stop("Missing compiled node id for target outcome ", target_label)
  }
  label_df <- AccumulatR:::native_outcome_labels_cpp(native_ctx)
  forced_rows <- label_df[label_df$label == forced_label, , drop = FALSE]
  if (nrow(forced_rows) != 1L) {
    stop("Forced label not found or duplicated: ", forced_label)
  }

  competitor_ids <- if (isTRUE(resolve_competitors)) {
    resolve_competitor_ids(prep, target_label, target_idx)
  } else {
    integer(0)
  }
  forced_complete <- integer(0)
  forced_survive <- integer(0)
  forced_label_id <- as.integer(forced_rows$label_id[[1]])
  if (forced_mode == "complete") {
    forced_complete <- forced_label_id
  } else {
    forced_survive <- forced_label_id
  }

  list(
    native_ctx = native_ctx,
    trial_rows = trial_rows,
    target_node_id = as.integer(compiled$id),
    competitor_ids = competitor_ids,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    upper = upper,
    forced_label = forced_label,
    forced_mode = forced_mode
  )
}

eval_case_once <- function(setup) {
  AccumulatR:::native_outcome_probability_params_cpp_idx(
    setup$native_ctx,
    setup$target_node_id,
    setup$upper,
    0L,
    setup$forced_complete,
    setup$forced_survive,
    setup$competitor_ids,
    AccumulatR:::.integrate_rel_tol(),
    AccumulatR:::.integrate_abs_tol(),
    12L,
    setup$trial_rows
  )
}

benchmark_case <- function(label, surface, setup, inner, n_rep) {
  AccumulatR:::unified_outcome_stats_reset_cpp()
  probability <- as.numeric(eval_case_once(setup))
  single_stats <- AccumulatR:::unified_outcome_stats_cpp()

  per_eval_sec <- numeric(n_rep)
  noderef_per_eval <- numeric(n_rep)
  scalar_fallback_per_eval <- numeric(n_rep)
  coupling_per_eval <- numeric(n_rep)
  kernel_density_per_eval <- numeric(n_rep)

  for (rep_idx in seq_len(n_rep)) {
    AccumulatR:::unified_outcome_stats_reset_cpp()
    elapsed <- system.time({
      for (inner_idx in seq_len(inner)) {
        eval_case_once(setup)
      }
    })[["elapsed"]]
    stats <- AccumulatR:::unified_outcome_stats_cpp()
    per_eval_sec[[rep_idx]] <- elapsed / inner
    noderef_per_eval[[rep_idx]] <-
      as.numeric(stats$generic_noderef_batch_calls %||% NA_real_) / inner
    scalar_fallback_per_eval[[rep_idx]] <-
      as.numeric(stats$generic_scalar_fallback_calls %||% NA_real_) / inner
    coupling_per_eval[[rep_idx]] <-
      as.numeric(stats$evaluate_outcome_coupling_unified_calls %||% NA_real_) / inner
    kernel_density_per_eval[[rep_idx]] <-
      as.numeric(stats$kernel_node_density_entry_idx_calls %||% NA_real_) / inner
  }

  data.frame(
    label = label,
    surface = surface,
    forced_mode = setup$forced_mode,
    forced_label = setup$forced_label,
    competitor_count = length(setup$competitor_ids),
    upper = setup$upper,
    probability = probability,
    single_generic_noderef_batch_calls =
      as.numeric(single_stats$generic_noderef_batch_calls %||% NA_real_),
    single_generic_scalar_fallback_calls =
      as.numeric(single_stats$generic_scalar_fallback_calls %||% NA_real_),
    single_evaluate_outcome_coupling_unified_calls =
      as.numeric(single_stats$evaluate_outcome_coupling_unified_calls %||% NA_real_),
    single_kernel_node_density_entry_idx_calls =
      as.numeric(single_stats$kernel_node_density_entry_idx_calls %||% NA_real_),
    median_per_eval_sec = stats::median(per_eval_sec),
    mean_per_eval_sec = mean(per_eval_sec),
    generic_noderef_batch_per_eval = stats::median(noderef_per_eval),
    generic_scalar_fallback_per_eval = stats::median(scalar_fallback_per_eval),
    coupling_calls_per_eval = stats::median(coupling_per_eval),
    kernel_density_calls_per_eval = stats::median(kernel_density_per_eval),
    inner_reps = inner,
    n_rep = n_rep,
    stringsAsFactors = FALSE
  )
}

build_pair_guard_case <- function() {
  structure <- race_spec() |>
    add_accumulator("x", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_accumulator("other", "lognormal") |>
    add_outcome("RESP", all_of("x", "gate")) |>
    add_outcome("NR_RAW", "other", options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
  params <- c(
    x.meanlog = log(0.34), x.sdlog = 0.17,
    gate.meanlog = log(0.25), gate.sdlog = 0.14,
    other.meanlog = log(0.58), other.sdlog = 0.18
  )
  setup <- resolve_outcome_setup(
    structure = structure,
    params = params,
    target_label = "RESP",
    forced_label = "NR_RAW",
    upper = 0.60,
    forced_mode = "complete",
    resolve_competitors = FALSE
  )
  list(
    label = "pair_guard_generic_noderef",
    surface = "direct_guarded_generic_noderef",
    setup = setup
  )
}

build_nested_guard_case <- function() {
  structure <- race_spec() |>
    add_accumulator("x", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("other", "lognormal") |>
    add_outcome("RESP", all_of("x", inhibit("a", by = "b"))) |>
    add_outcome("NR_RAW", "other", options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
  params <- c(
    x.meanlog = log(0.41), x.sdlog = 0.18,
    a.meanlog = log(0.31), a.sdlog = 0.18,
    b.meanlog = log(0.28), b.sdlog = 0.16,
    other.meanlog = log(0.57), other.sdlog = 0.18
  )
  setup <- resolve_outcome_setup(
    structure = structure,
    params = params,
    target_label = "RESP",
    forced_label = "NR_RAW",
    upper = 0.62,
    forced_mode = "complete",
    resolve_competitors = FALSE
  )
  list(
    label = "nested_guard_generic_noderef",
    surface = "direct_guarded_nested_noderef",
    setup = setup
  )
}

build_competitor_guard_case <- function() {
  structure <- race_spec() |>
    add_accumulator("x", "lognormal") |>
    add_accumulator("y", "lognormal") |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("other", "lognormal") |>
    add_outcome("R1", "x") |>
    add_outcome("R2", all_of("y", inhibit("a", by = "b"))) |>
    add_outcome("NR_RAW", "other", options = list(map_outcome_to = NA_character_)) |>
    finalize_model()
  params <- c(
    x.meanlog = log(0.30), x.sdlog = 0.16,
    y.meanlog = log(0.46), y.sdlog = 0.18,
    a.meanlog = log(0.33), a.sdlog = 0.18,
    b.meanlog = log(0.29), b.sdlog = 0.16,
    other.meanlog = log(0.63), other.sdlog = 0.19
  )
  setup <- resolve_outcome_setup(
    structure = structure,
    params = params,
    target_label = "R2",
    forced_label = "NR_RAW",
    upper = 0.56,
    forced_mode = "complete"
  )
  list(
    label = "competitor_guard_generic_noderef",
    surface = "competitor_bearing_generic_noderef",
    setup = setup
  )
}

cases <- list(
  build_pair_guard_case(),
  build_nested_guard_case(),
  build_competitor_guard_case()
)

inner <- as.integer(Sys.getenv("ACCUMULATR_FP5K_INNER", unset = "600"))
n_rep <- as.integer(Sys.getenv("ACCUMULATR_FP5K_N_REP", unset = "12"))
if (!is.finite(inner) || inner <= 0L) {
  stop("ACCUMULATR_FP5K_INNER must be a positive integer")
}
if (!is.finite(n_rep) || n_rep <= 0L) {
  stop("ACCUMULATR_FP5K_N_REP must be a positive integer")
}

results <- do.call(
  rbind,
  lapply(cases, function(cs) {
    row <- benchmark_case(cs$label, cs$surface, cs$setup, inner, n_rep)
    if (!is.finite(row$probability[[1]]) || row$probability[[1]] < 0.0 ||
        row$probability[[1]] > 1.0) {
      stop("Invalid probability for case ", cs$label)
    }
    if (!is.finite(row$generic_noderef_batch_per_eval[[1]]) ||
        row$generic_noderef_batch_per_eval[[1]] < 1.0) {
      stop("Case did not exercise generic NodeRef batch path: ", cs$label)
    }
    if (!is.finite(row$generic_scalar_fallback_per_eval[[1]]) ||
        row$generic_scalar_fallback_per_eval[[1]] != 0.0) {
      stop("Scalar fallback observed for case ", cs$label)
    }
    row
  })
)

utils::write.csv(results, stdout(), row.names = FALSE, quote = TRUE)
cat("bench_step5_kernel_target_cases complete\n", file = stderr())
