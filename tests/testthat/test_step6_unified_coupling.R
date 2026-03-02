race_spec <- AccumulatR::race_spec
add_accumulator <- AccumulatR::add_accumulator
add_outcome <- AccumulatR::add_outcome
all_of <- AccumulatR::all_of
add_trigger <- AccumulatR::add_trigger
after <- AccumulatR::after
finalize_model <- AccumulatR::finalize_model
build_param_matrix <- AccumulatR::build_param_matrix
build_likelihood_context <- AccumulatR::build_likelihood_context
log_likelihood <- AccumulatR::log_likelihood
response_probabilities <- AccumulatR::response_probabilities

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
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

build_shared_gate_nway_trigger_spec <- function() {
  race_spec() |>
    add_accumulator("x1", "lognormal") |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
    add_trigger("tg_x12", members = c("x1", "x2"), q = 0.10, draw = "shared") |>
    add_trigger("tg_x3g", members = c("x3", "gate"), q = 0.16, draw = "shared")
}

params_shared_gate_nway_trigger <- c(
  x1.meanlog = log(0.31), x1.sdlog = 0.18, x1.q = 0.10,
  x2.meanlog = log(0.34), x2.sdlog = 0.18, x2.q = 0.10,
  x3.meanlog = log(0.37), x3.sdlog = 0.18, x3.q = 0.16,
  gate.meanlog = log(0.23), gate.sdlog = 0.14, gate.q = 0.16
)

testthat::test_that("forced bits + coupling remains finite under native probability calls", {
  spec <- build_shared_gate_nway_trigger_spec()
  structure <- finalize_model(spec)
  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  native_ctx <- AccumulatR:::native_context_build(prep)

  params_mat <- build_param_matrix(spec, params_shared_gate_nway_trigger, n_trials = 1L)
  params_rows <- AccumulatR:::.param_matrix_to_rows(structure, params_mat)
  trial_rows <- params_rows[params_rows$trial == 1L, , drop = FALSE]

  outcome_labels <- vapply(prep$outcomes, function(x) as.character(x$label %||% ""), character(1))
  resp2_idx <- match("RESP2", outcome_labels)
  testthat::expect_true(is.finite(resp2_idx))
  compiled <- AccumulatR:::.expr_lookup_compiled(prep$outcomes[[resp2_idx]]$expr, prep)
  testthat::expect_true(!is.null(compiled) && !is.null(compiled$id))

  comp_ids <- resolve_competitor_ids(prep, "RESP2", resp2_idx)
  labels_df <- AccumulatR:::native_outcome_labels_cpp(native_ctx)
  forced_label <- labels_df$label_id[labels_df$label == "RESP3"]
  testthat::expect_true(length(forced_label) >= 1L)
  forced_label <- as.integer(forced_label[[1]])

  p_forced_complete <- AccumulatR:::native_outcome_probability_params_cpp_idx(
    native_ctx,
    as.integer(compiled$id),
    4.0,
    0L,
    as.integer(forced_label),
    integer(0),
    comp_ids,
    AccumulatR:::.integrate_rel_tol(),
    AccumulatR:::.integrate_abs_tol(),
    12L,
    trial_rows
  )
  p_forced_survive <- AccumulatR:::native_outcome_probability_params_cpp_idx(
    native_ctx,
    as.integer(compiled$id),
    4.0,
    0L,
    integer(0),
    as.integer(forced_label),
    comp_ids,
    AccumulatR:::.integrate_rel_tol(),
    AccumulatR:::.integrate_abs_tol(),
    12L,
    trial_rows
  )

  testthat::expect_true(is.finite(p_forced_complete))
  testthat::expect_true(is.finite(p_forced_survive))
  testthat::expect_gte(p_forced_complete, 0.0)
  testthat::expect_lte(p_forced_complete, 1.0)
  testthat::expect_gte(p_forced_survive, 0.0)
  testthat::expect_lte(p_forced_survive, 1.0)
})

testthat::test_that("chained onset coupling path stays deterministic and normalized", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP", all_of("b", "gate")) |>
    add_outcome("NR_RAW", all_of("a", "gate"), options = list(map_outcome_to = NA_character_))

  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.24), a.sdlog = 0.13, a.q = 0.00, a.t0 = 0.00,
    b.meanlog = log(0.31), b.sdlog = 0.15, b.q = 0.00, b.t0 = 0.00,
    gate.meanlog = log(0.27), gate.sdlog = 0.14, gate.q = 0.00, gate.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)

  probs1 <- response_probabilities(structure, params_df, include_na = TRUE)
  probs2 <- response_probabilities(structure, params_df, include_na = TRUE)

  testthat::expect_equal(probs1, probs2, tolerance = 1e-12)
  testthat::expect_true(all(is.finite(as.numeric(probs1))))
  testthat::expect_true("NA" %in% names(probs1))
  testthat::expect_equal(sum(as.numeric(probs1)), 1.0, tolerance = 1e-8)
})

testthat::test_that("ranked coupling path remains deterministic after unified cutover", {
  spec <- race_spec() |>
    add_accumulator("x2", "lognormal") |>
    add_accumulator("x3", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_outcome("RESP2", all_of("x2", "gate")) |>
    add_outcome("RESP3", all_of("x3", "gate")) |>
    add_trigger("tg_x23", members = c("x2", "x3"), q = 0.08, draw = "shared")

  structure <- finalize_model(spec)
  params <- c(
    x2.meanlog = log(0.33), x2.sdlog = 0.16, x2.q = 0.08, x2.t0 = 0.00,
    x3.meanlog = log(0.37), x3.sdlog = 0.16, x3.q = 0.08, x3.t0 = 0.00,
    gate.meanlog = log(0.24), gate.sdlog = 0.13, gate.q = 0.00, gate.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  ranked_df <- data.frame(
    trial = 1L,
    R = "RESP2", rt = 0.35,
    R2 = "RESP3", rt2 = 0.62,
    stringsAsFactors = FALSE
  )

  ctx <- build_likelihood_context(structure, ranked_df)
  ll1 <- as.numeric(log_likelihood(ctx, params_df))
  ll2 <- as.numeric(log_likelihood(ctx, params_df))

  testthat::expect_true(is.finite(ll1))
  testthat::expect_equal(ll1, ll2, tolerance = 1e-12)
})
