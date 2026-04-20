testthat::test_that("simple direct variant lowers to a dense runtime program", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "gamma") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  lowered <- .lower_direct_prep(prep)

  testthat::expect_length(lowered$variants, 1L)
  variant <- lowered$variants[[1]]
  program <- variant$program

  testthat::expect_equal(variant$component_id, "__default__")
  testthat::expect_equal(variant$leaf_ids, c("a", "b"))
  testthat::expect_equal(variant$outcome_labels, c("A", "B"))
  testthat::expect_equal(
    variant$param_keys,
    c("a.m", "a.s", "a.q", "a.t0", "b.shape", "b.rate", "b.q", "b.t0")
  )

  testthat::expect_equal(program$layout$n_leaves, 2L)
  testthat::expect_equal(program$layout$n_pools, 0L)
  testthat::expect_equal(program$layout$n_outcomes, 2L)
  testthat::expect_equal(program$layout$n_params, 8L)
  testthat::expect_equal(program$leaf_dist_kind, c(0L, 1L))
  testthat::expect_equal(program$parameter_layout$leaf_param_offsets, c(0L, 2L, 4L))
  testthat::expect_equal(program$parameter_layout$leaf_param_slots, c(0L, 1L, 4L, 5L))
  testthat::expect_equal(program$parameter_layout$leaf_q_slots, c(2L, 6L))
  testthat::expect_equal(program$parameter_layout$leaf_t0_slots, c(3L, 7L))
  testthat::expect_equal(program$outcome_source_kind, c(0L, 0L))
  testthat::expect_equal(program$outcome_source_index, c(0L, 1L))
  testthat::expect_equal(program$observed_label_index, c(0L, 1L))
})

testthat::test_that("projected shared trigger lowers as one independent runtime trigger", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_component("one", members = c("a")) |>
    add_trigger("tg", members = c("a", "b"), q = 0.2, draw = "shared")

  prep <- prepare_model(spec)
  lowered <- .lower_direct_prep(prep)

  variant <- lowered$variants[[1]]
  program <- variant$program

  testthat::expect_equal(program$layout$n_triggers, 1L)
  testthat::expect_equal(program$leaf_trigger_index, 0L)
  testthat::expect_equal(program$trigger_kind, 0L)
  testthat::expect_equal(program$trigger_has_fixed_q, 1L)
  testthat::expect_equal(program$trigger_fixed_q, 0.2)
  testthat::expect_equal(program$trigger_member_offsets, c(0L, 1L))
  testthat::expect_equal(program$trigger_member_indices, 0L)
  testthat::expect_equal(program$parameter_layout$leaf_q_slots, 2L)
  testthat::expect_equal(variant$param_keys, c("a.m", "a.s", "tg", "a.t0"))
})

testthat::test_that("pooled direct outcomes lower to pool topology without semantic trees", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_pool("ab", members = c("a", "b"), k = 1L) |>
    add_outcome("AB", "ab") |>
    add_outcome("C", "c")

  prep <- prepare_model(spec)
  lowered <- .lower_direct_prep(prep)

  variant <- lowered$variants[[1]]
  program <- variant$program

  testthat::expect_equal(variant$pool_ids, "ab")
  testthat::expect_equal(program$layout$n_pools, 1L)
  testthat::expect_equal(program$pool_k, 1L)
  testthat::expect_equal(program$pool_member_offsets, c(0L, 2L))
  testthat::expect_equal(program$pool_member_indices, c(0L, 1L))
  testthat::expect_equal(program$pool_member_kind, c(0L, 0L))
  testthat::expect_equal(program$outcome_source_kind, c(1L, 0L))
  testthat::expect_equal(program$outcome_source_index, c(0L, 2L))
})

testthat::test_that("exact variants are rejected instead of half-lowering", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)

  testthat::expect_error(
    .lower_direct_prep(prep),
    "cannot lower exact variant.*ranked observation"
  )
})
