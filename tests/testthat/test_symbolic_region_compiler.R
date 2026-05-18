compiler_total_metrics <- function(structure) {
  complexity_metrics(make_context(structure, diagnostics = TRUE))$total
}

expect_clean_region_cells <- function(metrics) {
  testthat::expect_equal(metrics[["negative_symbolic_cells"]], 0L)
  testthat::expect_equal(metrics[["overlapping_symbolic_cell_pairs"]], 0L)
}

expect_no_generic_integrals <- function(metrics) {
  testthat::expect_equal(metrics[["generic_integral_kernels"]], 0L)
}

expect_complexity_budget <- function(metrics,
                                     max_integral_depth,
                                     max_integral_nodes,
                                     max_symbolic_cells) {
  testthat::expect_true(metrics[["max_integral_depth"]] <= max_integral_depth)
  testthat::expect_true(metrics[["integral_nodes"]] <= max_integral_nodes)
  testthat::expect_true(metrics[["symbolic_cells"]] <= max_symbolic_cells)
}

build_simple_first_of_model <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    finalize_model()
}

build_all_of_competitor_model <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("C", "lognormal") |>
    add_outcome("AB", all_of("A", "B")) |>
    add_outcome("C", "C") |>
    finalize_model()
}

build_shared_gate_competitor_model <- function(reverse = FALSE) {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_accumulator("d", "lognormal")

  if (reverse) {
    spec |>
      add_outcome("D", "d") |>
      add_outcome("C2", all_of("gate", "b")) |>
      add_outcome("C1", all_of("gate", "a")) |>
      test_separate_all_parameters() |>
      finalize_model()
  } else {
    spec |>
      add_outcome("C1", all_of("a", "gate")) |>
      add_outcome("C2", all_of("b", "gate")) |>
      add_outcome("D", "d") |>
      test_separate_all_parameters() |>
      finalize_model()
  }
}

build_guarded_overlap_model <- function() {
  race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("C1", inhibit("a", by = "stop")) |>
    add_outcome("C2", inhibit("b", by = "stop")) |>
    add_outcome("D", "d") |>
    finalize_model()
}

build_expr_distribution_model <- function(expr) {
  race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", expr) |>
    add_outcome("D", "d") |>
    finalize_model()
}

testthat::test_that("compiler complexity metrics expose structural acceptance fields", {
  metrics <- compiler_total_metrics(build_simple_first_of_model())
  testthat::expect_named(
    metrics,
    c(
      "symbolic_regions",
      "symbolic_cells",
      "max_symbolic_cells_per_region",
      "negative_symbolic_cells",
      "overlapping_symbolic_cell_pairs",
      "expr_relation_atoms",
      "compiled_roots",
      "compiled_nodes",
      "integral_nodes",
      "integral_kernels",
      "source_product_integral_kernels",
      "generic_integral_kernels",
      "max_integral_depth"
    )
  )
})

testthat::test_that("simple independent and max-transition regions stay closed form", {
  simple_metrics <- compiler_total_metrics(build_simple_first_of_model())
  expect_clean_region_cells(simple_metrics)
  expect_no_generic_integrals(simple_metrics)
  expect_complexity_budget(
    simple_metrics,
    max_integral_depth = 0L,
    max_integral_nodes = 0L,
    max_symbolic_cells = 4L
  )

  all_of_metrics <- compiler_total_metrics(build_all_of_competitor_model())
  expect_clean_region_cells(all_of_metrics)
  expect_no_generic_integrals(all_of_metrics)
  expect_complexity_budget(
    all_of_metrics,
    max_integral_depth = 0L,
    max_integral_nodes = 0L,
    max_symbolic_cells = 7L
  )
})

testthat::test_that("shared-gate and guarded-overlap competitors materialize without generic fallback", {
  shared_metrics <- compiler_total_metrics(build_shared_gate_competitor_model())
  expect_clean_region_cells(shared_metrics)
  expect_no_generic_integrals(shared_metrics)
  expect_complexity_budget(
    shared_metrics,
    max_integral_depth = 1L,
    max_integral_nodes = 4L,
    max_symbolic_cells = 31L
  )

  guarded_metrics <- compiler_total_metrics(build_guarded_overlap_model())
  expect_clean_region_cells(guarded_metrics)
  expect_no_generic_integrals(guarded_metrics)
  expect_complexity_budget(
    guarded_metrics,
    max_integral_depth = 2L,
    max_integral_nodes = 14L,
    max_symbolic_cells = 21L
  )
})

testthat::test_that("first_of expression distributions preserve cheap union cases", {
  cases <- list(
    independent_child_union = list(
      expr = first_of("a", "b"),
      max_roots = 11L,
      max_nodes = 23L,
      max_cells = 6L
    ),
    overlapping_child_union = list(
      expr = first_of(all_of("a", "g"), all_of("b", "g")),
      max_roots = 15L,
      max_nodes = 34L,
      max_cells = 8L
    ),
    absorbed_union = list(
      expr = first_of(all_of("a", "g"), "g"),
      max_roots = 6L,
      max_nodes = 14L,
      max_cells = 4L
    ),
    multi_child_union = list(
      expr = first_of("a", "b", "c"),
      max_roots = 14L,
      max_nodes = 30L,
      max_cells = 8L
    ),
    nested_first_of = list(
      expr = first_of("a", first_of("b", "c")),
      max_roots = 14L,
      max_nodes = 30L,
      max_cells = 8L
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    metrics <- compiler_total_metrics(
      build_expr_distribution_model(case$expr))
    expect_clean_region_cells(metrics)
    expect_no_generic_integrals(metrics)
    testthat::expect_equal(metrics[["integral_nodes"]], 0L, info = case_name)
    testthat::expect_equal(metrics[["integral_kernels"]], 0L, info = case_name)
    testthat::expect_true(
      metrics[["compiled_roots"]] <= case$max_roots, info = case_name)
    testthat::expect_true(
      metrics[["compiled_nodes"]] <= case$max_nodes, info = case_name)
    testthat::expect_true(
      metrics[["symbolic_cells"]] <= case$max_cells, info = case_name)
  }
})

testthat::test_that("common conjunct first_of canonicalizes to factored all_of cost", {
  unfactored <- compiler_total_metrics(
    build_expr_distribution_model(
      first_of(all_of("a", "g"), all_of("b", "g"))))
  factored <- compiler_total_metrics(
    build_expr_distribution_model(
      all_of(first_of("a", "b"), "g")))
  fields <- c(
    "symbolic_regions",
    "symbolic_cells",
    "expr_relation_atoms",
    "compiled_roots",
    "compiled_nodes",
    "integral_nodes",
    "generic_integral_kernels"
  )
  testthat::expect_equal(unfactored[fields], factored[fields])
})

testthat::test_that("all_of and simple guard expression distributions stay compact", {
  cases <- list(
    all_of_pair = list(
      expr = all_of("a", "b"),
      max_roots = 11L,
      max_nodes = 24L,
      max_integrals = 0L,
      max_depth = 0L,
      max_cells = 6L
    ),
    all_of_three = list(
      expr = all_of("a", "b", "c"),
      max_roots = 14L,
      max_nodes = 31L,
      max_integrals = 0L,
      max_depth = 0L,
      max_cells = 8L
    ),
    simple_guard = list(
      expr = inhibit("a", by = "g"),
      max_roots = 8L,
      max_nodes = 19L,
      max_integrals = 1L,
      max_depth = 1L,
      max_cells = 5L
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    metrics <- compiler_total_metrics(
      build_expr_distribution_model(case$expr))
    expect_clean_region_cells(metrics)
    expect_no_generic_integrals(metrics)
    testthat::expect_true(
      metrics[["compiled_roots"]] <= case$max_roots, info = case_name)
    testthat::expect_true(
      metrics[["compiled_nodes"]] <= case$max_nodes, info = case_name)
    testthat::expect_true(
      metrics[["integral_nodes"]] <= case$max_integrals, info = case_name)
    testthat::expect_true(
      metrics[["max_integral_depth"]] <= case$max_depth, info = case_name)
    testthat::expect_true(
      metrics[["symbolic_cells"]] <= case$max_cells, info = case_name)
  }
})

testthat::test_that("shared-gate likelihood is invariant to child and outcome order", {
  params <- c(
    a.m = log(0.29), a.s = 0.17, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.16, b.t0 = 0.00,
    gate.m = log(0.24), gate.s = 0.14, gate.t0 = 0.00,
    d.m = log(0.46), d.s = 0.19, d.t0 = 0.00
  )
  data <- data.frame(
    trials = 1L,
    R = "D",
    rt = 0.43,
    stringsAsFactors = FALSE
  )

  base <- build_shared_gate_competitor_model(reverse = FALSE)
  reversed <- build_shared_gate_competitor_model(reverse = TRUE)

  base_data <- prepare_data(base, data)
  reversed_data <- prepare_data(reversed, data)
  base_ll <- log_likelihood(
    make_context(base),
    base_data,
    build_param_matrix(base, params, trial_df = base_data)
  )
  reversed_ll <- log_likelihood(
    make_context(reversed),
    reversed_data,
    build_param_matrix(reversed, params, trial_df = reversed_data)
  )

  testthat::expect_equal(
    as.numeric(base_ll),
    as.numeric(reversed_ll),
    tolerance = 1e-10
  )
})
