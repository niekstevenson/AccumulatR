.repo_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)

source(file.path(.repo_root, "R", "helpers.R"))
source(file.path(.repo_root, "R", "model_definition.R"))
source(file.path(.repo_root, "R", "semantic_bridge.R"))

testthat::test_that("projected simple race is classified as direct", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "gamma") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  prep <- prepare_model(spec)
  compiled <- .project_semantic_prep(prep, rebuild = TRUE, root = .repo_root)

  testthat::expect_length(compiled$variants, 1L)
  testthat::expect_equal(compiled$variants[[1]]$component_id, "__default__")
  testthat::expect_equal(compiled$variants[[1]]$backend, "direct")
  testthat::expect_equal(
    vapply(compiled$variants[[1]]$direct_outcomes, `[[`, character(1), "label"),
    c("A", "B")
  )
  testthat::expect_equal(
    vapply(compiled$variants[[1]]$leaves, `[[`, character(1), "id"),
    c("a", "b")
  )
  testthat::expect_length(compiled$variants[[1]]$backend_reasons, 0L)
})

testthat::test_that("projection removes dead guarded branches inside a component", {
  spec <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("S1", "lognormal") |>
    add_accumulator("IS", "lognormal") |>
    add_accumulator("S2", "lognormal") |>
    add_outcome(
      "A",
      first_of(
        inhibit("A", by = "S1"),
        all_of("A", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "B",
      first_of(
        inhibit("B", by = "S1"),
        all_of("B", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "STOP",
      all_of("S1", inhibit("S2", by = "IS")),
      options = list(map_outcome_to = NA_character_)
    ) |>
    add_component("go", members = c("A", "B")) |>
    add_component("stop", members = c("A", "B", "S1", "IS", "S2"))

  prep <- prepare_model(spec)
  compiled <- .project_semantic_prep(prep, root = .repo_root)

  by_id <- setNames(compiled$variants, vapply(compiled$variants, `[[`, character(1), "component_id"))

  testthat::expect_equal(by_id$go$backend, "direct")
  testthat::expect_equal(
    vapply(by_id$go$leaves, `[[`, character(1), "id"),
    c("A", "B")
  )
  testthat::expect_equal(
    vapply(by_id$go$outcomes, `[[`, character(1), "label"),
    c("A", "B")
  )
  testthat::expect_true(all(vapply(by_id$go$expr_nodes, `[[`, character(1), "kind") == "event"))

  testthat::expect_equal(by_id$stop$backend, "exact")
  testthat::expect_true("non-direct outcome" %in% by_id$stop$backend_reasons)
  testthat::expect_true("outcome remapping" %in% by_id$stop$backend_reasons)
})

testthat::test_that("shared triggers are projected away when only one member survives", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_component("one", members = c("a")) |>
    add_component("both", members = c("a", "b")) |>
    add_trigger("tg", members = c("a", "b"), q = 0.2, draw = "shared")

  prep <- prepare_model(spec)
  compiled <- .project_semantic_prep(prep, root = .repo_root)

  by_id <- setNames(compiled$variants, vapply(compiled$variants, `[[`, character(1), "component_id"))

  testthat::expect_equal(by_id$one$backend, "direct")
  testthat::expect_false("shared trigger" %in% by_id$one$backend_reasons)
  testthat::expect_equal(vapply(by_id$one$leaves, `[[`, character(1), "id"), "a")

  testthat::expect_equal(by_id$both$backend, "exact")
  testthat::expect_true("shared trigger" %in% by_id$both$backend_reasons)
})

testthat::test_that("component observation overrides classify only the widened variant as exact", {
  spec <- race_spec(n_outcomes = 1L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_component("single", members = c("a", "b"), n_outcomes = 1L, weight = 0.5) |>
    add_component("double", members = c("a", "b"), n_outcomes = 2L, weight = 0.5)

  prep <- prepare_model(spec)
  compiled <- .project_semantic_prep(prep, root = .repo_root)

  by_id <- setNames(compiled$variants, vapply(compiled$variants, `[[`, character(1), "component_id"))

  testthat::expect_equal(by_id$single$observation$n_outcomes, 1L)
  testthat::expect_equal(by_id$single$backend, "direct")

  testthat::expect_equal(by_id$double$observation$n_outcomes, 2L)
  testthat::expect_equal(by_id$double$backend, "exact")
  testthat::expect_true("ranked observation" %in% by_id$double$backend_reasons)
})
