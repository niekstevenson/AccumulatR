generated_region_metrics <- function(structure) {
  complexity_metrics(make_context(structure, diagnostics = TRUE))$total
}

generated_loglik <- function(structure, params, rt) {
  data <- data.frame(
    trial = 1L,
    R = "R",
    rt = rt,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data)
  params_df <- build_param_matrix(
    structure$model_spec, params, trial_df = prepared)
  as.numeric(log_likelihood(make_context(structure), prepared, params_df))
}

generated_shape_model <- function(target_expr, competitor_expr) {
  race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_accumulator("s", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_outcome("R", target_expr) |>
    add_outcome("C", competitor_expr) |>
    finalize_model()
}

testthat::test_that("bounded generated grammar shapes compile within symbolic budgets", {
  targets <- list(
    source = function() "a",
    all_pair = function() all_of("a", "b"),
    first_pair = function() first_of("a", "b"),
    absence = function() all_of("a", none_of("s")),
    guarded_source = function() inhibit("a", by = "s"),
    guarded_all = function() inhibit(all_of("a", "b"), by = "s"),
    composite_blocker = function() inhibit("a", by = all_of("s", "g")),
    nested_absence_choice = function() {
      first_of(
        all_of("a", none_of("s")),
        inhibit("b", by = "c")
      )
    }
  )
  competitors <- list(
    source = function() "d",
    shared_gate = function() all_of("b", "g"),
    guarded = function() inhibit("c", by = "s")
  )
  params <- c(
    a.m = log(0.31), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.36), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
    c.m = log(0.29), c.s = 0.17, c.q = 0.00, c.t0 = 0.00,
    d.m = log(0.48), d.s = 0.18, d.q = 0.00, d.t0 = 0.00,
    s.m = log(0.27), s.s = 0.14, s.q = 0.00, s.t0 = 0.00,
    g.m = log(0.34), g.s = 0.15, g.q = 0.00, g.t0 = 0.00
  )
  for (target_name in names(targets)) {
    for (competitor_name in names(competitors)) {
      shape_id <- paste(target_name, "vs", competitor_name)
      structure <- generated_shape_model(
        targets[[target_name]](),
        competitors[[competitor_name]]()
      )
      metrics <- generated_region_metrics(structure)
      ll <- generated_loglik(structure, params, 0.43)

      testthat::expect_true(is.finite(ll), info = shape_id)
      testthat::expect_equal(
        metrics[["negative_symbolic_cells"]], 0L, info = shape_id)
      testthat::expect_equal(
        metrics[["overlapping_symbolic_cell_pairs"]], 0L, info = shape_id)
      testthat::expect_equal(
        metrics[["generic_integral_kernels"]], 0L, info = shape_id)
      testthat::expect_true(
        metrics[["max_integral_depth"]] <= 4L, info = shape_id)
      testthat::expect_true(
        metrics[["integral_nodes"]] <= 80L, info = shape_id)
      testthat::expect_true(
        metrics[["symbolic_cells"]] <= 140L, info = shape_id)
    }
  }
})
