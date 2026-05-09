metamorphic_loglik <- function(structure, params, response, rt) {
  data <- data.frame(
    trial = 1L,
    R = response,
    rt = rt,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data)
  params_df <- build_param_matrix(
    structure$model_spec, params, trial_df = prepared)
  as.numeric(log_likelihood(make_context(structure), prepared, params_df))
}

testthat::test_that("all_of lowering is order and association invariant", {
  base <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", all_of("a", all_of("b", "g"))) |>
    add_outcome("D", "d") |>
    finalize_model()
  reordered <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", all_of("g", "b", "a")) |>
    add_outcome("D", "d") |>
    finalize_model()
  params <- c(
    a.m = log(0.29), a.s = 0.16, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.35), b.s = 0.17, b.q = 0.00, b.t0 = 0.00,
    g.m = log(0.24), g.s = 0.14, g.q = 0.00, g.t0 = 0.00,
    d.m = log(0.47), d.s = 0.18, d.q = 0.00, d.t0 = 0.00
  )

  testthat::expect_equal(
    metamorphic_loglik(base, params, "R", 0.42),
    metamorphic_loglik(reordered, params, "R", 0.42),
    tolerance = 1e-10
  )
})

testthat::test_that("first_of lowering is order and association invariant", {
  base <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", first_of("a", first_of("b", "c"))) |>
    add_outcome("D", "d") |>
    finalize_model()
  reordered <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", first_of("c", "b", "a")) |>
    add_outcome("D", "d") |>
    finalize_model()
  params <- c(
    a.m = log(0.28), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.33), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
    c.m = log(0.38), c.s = 0.17, c.q = 0.00, c.t0 = 0.00,
    d.m = log(0.46), d.s = 0.18, d.q = 0.00, d.t0 = 0.00
  )

  testthat::expect_equal(
    metamorphic_loglik(base, params, "R", 0.31),
    metamorphic_loglik(reordered, params, "R", 0.31),
    tolerance = 1e-10
  )
})

testthat::test_that("repeated shared subexpressions equal their flattened form", {
  repeated <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", all_of(all_of("a", "g"), all_of("b", "g"))) |>
    add_outcome("D", "d") |>
    finalize_model()
  flattened <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("g", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("R", all_of("a", "b", "g")) |>
    add_outcome("D", "d") |>
    finalize_model()
  params <- c(
    a.m = log(0.30), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.36), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
    g.m = log(0.25), g.s = 0.13, g.q = 0.00, g.t0 = 0.00,
    d.m = log(0.48), d.s = 0.18, d.q = 0.00, d.t0 = 0.00
  )

  testthat::expect_equal(
    metamorphic_loglik(repeated, params, "R", 0.41),
    metamorphic_loglik(flattened, params, "R", 0.41),
    tolerance = 2e-3
  )
})
