testthat::test_that("exact variant lowers to a dense runtime program", {
  spec <- race_spec() |>
    add_accumulator("s", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("change", "lognormal") |>
    add_outcome("S", inhibit("s", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_trigger("tg", members = c("stop", "change"), q = 0.2, draw = "shared")

  prep <- prepare_model(spec)
  lowered <- .lower_exact_prep(prep)

  testthat::expect_length(lowered$variants, 1L)
  variant <- lowered$variants[[1]]
  program <- variant$program

  testthat::expect_equal(variant$component_id, "__default__")
  testthat::expect_equal(variant$leaf_ids, c("s", "stop", "change"))
  testthat::expect_equal(variant$outcome_labels, c("S", "X"))

  testthat::expect_equal(program$layout$n_leaves, 3L)
  testthat::expect_equal(program$layout$n_triggers, 1L)
  testthat::expect_equal(program$trigger_kind, 1L)
  testthat::expect_equal(program$trigger_member_offsets, c(0L, 2L))
  testthat::expect_equal(program$trigger_member_indices, c(1L, 2L))
  testthat::expect_equal(program$parameter_layout$leaf_q_slots, c(2L, 6L, 6L))
  testthat::expect_length(program$expr_kind, 6L)
  testthat::expect_equal(program$outcome_expr_root, c(2L, 5L))
})

testthat::test_that("direct variants are rejected instead of half-lowering to exact", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)

  testthat::expect_error(
    .lower_exact_prep(prep),
    "cannot lower direct variant"
  )
})
