testthat::test_that("semantic compiler handles direct top-k prep", {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "gamma") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  compiled <- .compile_semantic_prep(prep)

  testthat::expect_equal(vapply(compiled$leaves, `[[`, character(1), "id"), c("a", "b"))
  testthat::expect_equal(vapply(compiled$leaves, `[[`, character(1), "dist"), c("lognormal", "gamma"))
  testthat::expect_equal(compiled$observation$n_outcomes, 2L)
  testthat::expect_length(compiled$validation_issues, 0L)
})

testthat::test_that("semantic compiler preserves onset, pools, guards, and triggers", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a", lag = 0.1)) |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_outcome("A", "a") |>
    add_outcome("B", inhibit(all_of("b", "p"), by = "a")) |>
    add_trigger("tg", members = c("a", "b"), q = 0.2, draw = "shared")

  prep <- prepare_model(spec)
  compiled <- .compile_semantic_prep(prep)

  testthat::expect_equal(compiled$leaves[[2]]$onset$kind, "after")
  testthat::expect_equal(compiled$leaves[[2]]$onset$source_kind, "leaf")
  testthat::expect_equal(compiled$leaves[[2]]$onset$source_id, "a")
  testthat::expect_equal(compiled$leaves[[1]]$q_key, "tg")
  testthat::expect_equal(compiled$triggers[[1]]$id, "tg")
  testthat::expect_equal(compiled$triggers[[1]]$members, c("a", "b"))
  testthat::expect_equal(compiled$pools[[1]]$members, c("a", "b"))
  testthat::expect_true(any(vapply(compiled$expr_nodes, `[[`, character(1), "kind") == "guard"))
  testthat::expect_length(compiled$validation_issues, 0L)
})
