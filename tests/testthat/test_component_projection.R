testthat::test_that("go-component probabilities match the simple A/B race", {
  guarded <- race_spec() |>
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
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    finalize_model()

  simple <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    finalize_model()

  params_guarded <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00,
    B.m = log(0.33), B.s = 0.18, B.q = 0.00, B.t0 = 0.00,
    S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
    IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
    S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
  )
  params_simple <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00,
    B.m = log(0.33), B.s = 0.18, B.q = 0.00, B.t0 = 0.00
  )

  data_guarded <- data.frame(
    trial = 1:2,
    component = "go",
    R = c("A", "B"),
    rt = c(0.42, 0.47),
    stringsAsFactors = FALSE
  )
  data_simple <- data.frame(
    trial = 1:2,
    R = c("A", "B"),
    rt = c(0.42, 0.47),
    stringsAsFactors = FALSE
  )

  params_df_guarded <- data.frame(
    trial = 1L,
    component = "go",
    accumulator = 1:5,
    q = 0,
    t0 = 0,
    m = c(log(0.30), log(0.33), log(0.26), log(0.35), log(0.32)),
    s = c(0.18, 0.18, 0.18, 0.18, 0.18),
    stringsAsFactors = FALSE
  )
  params_df_simple <- data.frame(
    trial = 1L,
    accumulator = 1:2,
    q = 0,
    t0 = 0,
    m = c(log(0.30), log(0.33)),
    s = c(0.18, 0.18),
    stringsAsFactors = FALSE
  )

  probs_guarded <- response_probabilities(guarded, params_df_guarded, include_na = TRUE)
  probs_simple <- response_probabilities(simple, params_df_simple, include_na = TRUE)

  testthat::expect_equal(probs_guarded, probs_simple, tolerance = 1e-8)
})

testthat::test_that("go-component probabilities drop impossible chained-onset branches", {
  guarded <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("S", "lognormal") |>
    add_accumulator("C", "lognormal", onset = after("S")) |>
    add_outcome("A", "A") |>
    add_outcome("C", "C") |>
    add_component("go", members = c("A")) |>
    add_component("stop", members = c("A", "S", "C")) |>
    finalize_model()

  simple <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_outcome("A", "A") |>
    finalize_model()

  params_guarded <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00,
    S.m = log(0.26), S.s = 0.18, S.q = 0.00, S.t0 = 0.00,
    C.m = log(0.35), C.s = 0.18, C.q = 0.00, C.t0 = 0.00
  )
  params_simple <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00
  )

  data_guarded <- data.frame(
    trial = 1,
    component = "go",
    R = "A",
    rt = 0.42,
    stringsAsFactors = FALSE
  )
  data_simple <- data.frame(
    trial = 1,
    R = "A",
    rt = 0.42,
    stringsAsFactors = FALSE
  )

  params_df_guarded <- data.frame(
    trial = 1L,
    component = "go",
    accumulator = 1:3,
    q = 0,
    t0 = 0,
    m = c(log(0.30), log(0.26), log(0.35)),
    s = c(0.18, 0.18, 0.18),
    stringsAsFactors = FALSE
  )
  params_df_simple <- data.frame(
    trial = 1L,
    accumulator = 1L,
    q = 0,
    t0 = 0,
    m = log(0.30),
    s = 0.18,
    stringsAsFactors = FALSE
  )

  probs_guarded <- response_probabilities(guarded, params_df_guarded, include_na = TRUE)
  probs_simple <- response_probabilities(simple, params_df_simple, include_na = TRUE)

  testthat::expect_equal(unname(probs_guarded[["A"]]), unname(probs_simple[["A"]]), tolerance = 1e-8)
  testthat::expect_equal(unname(probs_guarded[["NA"]]), unname(probs_simple[["NA"]]), tolerance = 1e-8)
  testthat::expect_equal(unname(probs_guarded[["C"]]), 0, tolerance = 1e-12)
})

testthat::test_that("go-component likelihood matches the simple A/B race", {
  guarded <- race_spec() |>
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
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    finalize_model()

  simple <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    finalize_model()

  params_guarded <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00,
    B.m = log(0.33), B.s = 0.18, B.q = 0.00, B.t0 = 0.00,
    S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
    IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
    S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
  )
  params_simple <- c(
    A.m = log(0.30), A.s = 0.18, A.q = 0.00, A.t0 = 0.00,
    B.m = log(0.33), B.s = 0.18, B.q = 0.00, B.t0 = 0.00
  )

  data_guarded <- data.frame(
    trial = 1:3,
    component = "go",
    R = c("A", "B", "A"),
    rt = c(0.42, 0.47, 0.51),
    stringsAsFactors = FALSE
  )
  data_simple <- data.frame(
    trial = 1:3,
    R = c("A", "B", "A"),
    rt = c(0.42, 0.47, 0.51),
    stringsAsFactors = FALSE
  )

  prepared_guarded <- prepare_data(guarded, data_guarded)
  prepared_simple <- prepare_data(simple, data_simple)
  ctx_guarded <- make_context(guarded)
  ctx_simple <- make_context(simple)
  params_df_guarded <- build_param_matrix(
    guarded$model_spec,
    params_guarded,
    trial_df = prepared_guarded
  )
  params_df_simple <- build_param_matrix(
    simple$model_spec,
    params_simple,
    trial_df = prepared_simple
  )

  ll_guarded <- as.numeric(log_likelihood(ctx_guarded, prepared_guarded, params_df_guarded))
  ll_simple <- as.numeric(log_likelihood(ctx_simple, prepared_simple, params_df_simple))

  testthat::expect_equal(ll_guarded, ll_simple, tolerance = 1e-8)
})
