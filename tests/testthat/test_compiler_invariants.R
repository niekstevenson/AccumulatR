compiler_metrics <- function(expr) {
  model <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", expr) |>
    add_outcome("D", "d") |>
    finalize_model()
  complexity_metrics(make_context(model, diagnostics = TRUE))$total
}

expect_compiler_budget <- function(metrics, budget, case) {
  testthat::expect_equal(metrics[["negative_symbolic_cells"]], 0L, info = case)
  testthat::expect_equal(metrics[["overlapping_symbolic_cell_pairs"]], 0L, info = case)
  testthat::expect_equal(metrics[["generic_integral_kernels"]], 0L, info = case)
  for (metric in names(budget)) {
    testthat::expect_true(metrics[[metric]] <= budget[[metric]], info = case)
  }
}

testthat::test_that("first_of distributions retain compact closed-form programs", {
  cases <- list(
    independent = list(
      first_of("a", "b"),
      c(integral_nodes = 0L, integral_kernels = 0L,
        compiled_roots = 12L, compiled_nodes = 24L, symbolic_cells = 6L)
    ),
    overlapping = list(
      first_of(all_of("a", "g"), all_of("b", "g")),
      c(integral_nodes = 0L, integral_kernels = 0L,
        compiled_roots = 16L, compiled_nodes = 35L, symbolic_cells = 8L)
    ),
    absorbed = list(
      first_of(all_of("a", "g"), "g"),
      c(integral_nodes = 0L, integral_kernels = 0L,
        compiled_roots = 7L, compiled_nodes = 15L, symbolic_cells = 4L)
    ),
    multi_child = list(
      first_of("a", "b", "c"),
      c(integral_nodes = 0L, integral_kernels = 0L,
        compiled_roots = 15L, compiled_nodes = 31L, symbolic_cells = 8L)
    )
  )
  for (case in names(cases)) {
    expect_compiler_budget(compiler_metrics(cases[[case]][[1L]]), cases[[case]][[2L]], case)
  }
})

testthat::test_that("all_of and guard distributions retain compact programs", {
  cases <- list(
    all_of_three = list(
      all_of("a", "b", "c"),
      c(compiled_roots = 15L, compiled_nodes = 32L, integral_nodes = 0L,
        max_integral_depth = 0L, symbolic_cells = 8L)
    ),
    simple_guard = list(
      inhibit("a", by = "g"),
      c(compiled_roots = 9L, compiled_nodes = 20L, integral_nodes = 1L,
        max_integral_depth = 1L, symbolic_cells = 5L)
    )
  )
  for (case in names(cases)) {
    expect_compiler_budget(compiler_metrics(cases[[case]][[1L]]), cases[[case]][[2L]], case)
  }
})

testthat::test_that("likelihood is invariant to outcome declaration order", {
  build <- function(reverse) {
    spec <- race_spec() |>
      add_accumulator("a", "lognormal") |>
      add_accumulator("b", "lognormal") |>
      add_accumulator("gate", "lognormal") |>
      add_accumulator("d", "lognormal")
    outcomes <- list(
      C1 = all_of("a", "gate"),
      C2 = all_of("b", "gate"),
      D = "d"
    )
    if (reverse) outcomes <- rev(outcomes)
    for (label in names(outcomes)) {
      spec <- add_outcome(spec, label, outcomes[[label]])
    }
    spec |> test_separate_all_parameters() |> finalize_model()
  }
  params <- c(
    a.m = log(0.29), a.s = 0.17, a.t0 = 0,
    b.m = log(0.35), b.s = 0.16, b.t0 = 0,
    gate.m = log(0.24), gate.s = 0.14, gate.t0 = 0,
    d.m = log(0.46), d.s = 0.19, d.t0 = 0
  )
  evaluate <- function(model) {
    data <- prepare_data(model, data.frame(trials = 1L, R = "D", rt = 0.43))
    log_likelihood(
      make_context(model), data,
      build_param_matrix(model, params, n_trials = 1L)
    )
  }
  testthat::expect_equal(evaluate(build(FALSE)), evaluate(build(TRUE)), tolerance = 1e-10)
})
