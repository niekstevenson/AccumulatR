stop_na_model <- function() {
  race_spec() |>
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
    set_parameters(list(
      m_go = c("A.m", "B.m"),
      s_go = c("A.s", "B.s"),
      t0_go = c("A.t0", "B.t0")
    )) |>
    finalize_model()
}

stop_na_params <- c(
  m_go = log(0.30), s_go = 0.18, t0_go = 0.00,
  S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
  IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
  S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
)

testthat::test_that("prepare_data compresses repeated successful-stop trials", {
  model <- stop_na_model()
  data_df <- data.frame(
    trial = 1:6,
    component = "stop",
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )

  prepared_full <- prepare_data(model, data_df, compress = FALSE)
  prepared_comp <- prepare_data(model, data_df, compress = TRUE)

  testthat::expect_null(attr(prepared_full, "expand", exact = TRUE))
  testthat::expect_equal(attr(prepared_comp, "expand", exact = TRUE), rep(1L, 6L))
  testthat::expect_equal(max(prepared_comp$trial), 1L)

  ctx <- make_context(model)
  pm_full <- build_param_matrix(model$model_spec, stop_na_params, trial_df = prepared_full)
  pm_comp <- build_param_matrix(model$model_spec, stop_na_params, trial_df = prepared_comp)

  ll_full <- as.numeric(log_likelihood(ctx, prepared_full, pm_full))
  ll_comp <- as.numeric(log_likelihood(ctx, prepared_comp, pm_comp))

  testthat::expect_equal(ll_comp, ll_full, tolerance = 1e-10)
})

testthat::test_that("compressed prepared data honors original-length ok masks", {
  model <- stop_na_model()
  data_df <- data.frame(
    trial = 1:6,
    component = "stop",
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  ok <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)

  prepared_full <- prepare_data(model, data_df, compress = FALSE)
  prepared_comp <- prepare_data(model, data_df, compress = TRUE)

  ctx <- make_context(model)
  pm_full <- build_param_matrix(model$model_spec, stop_na_params, trial_df = prepared_full)
  pm_comp <- build_param_matrix(model$model_spec, stop_na_params, trial_df = prepared_comp)

  ll_full <- as.numeric(log_likelihood(ctx, prepared_full, pm_full, ok = ok))
  ll_comp <- as.numeric(log_likelihood(ctx, prepared_comp, pm_comp, ok = ok))

  testthat::expect_equal(ll_comp, ll_full, tolerance = 1e-10)
})
