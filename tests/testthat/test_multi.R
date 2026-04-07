testthat::test_that("multi models agree across simulate/probability/likelihood", {
  # Helper to compute manual log-likelihood from response probabilities
  manual_loglik <- function(probs, outcomes) {
    labels <- as.character(outcomes)
    labels[is.na(labels)] <- "NA"
    if (!all(labels %in% names(probs))) {
      missing <- setdiff(unique(labels), names(probs))
      stop("Missing probabilities for outcomes: ", paste(missing, collapse = ", "))
    }
    sum(log(probs[labels]))
  }

  # Model definitions (copied from dev/scripts/vignettes/*.Rmd) --------------
  vignette_multi_outcome <- function() {
    race_spec(n_outcomes = 2L) |>
      add_accumulator("A", "lognormal") |>
      add_accumulator("B", "lognormal") |>
      add_outcome("A", "A") |>
      add_outcome("B", "B") |>
      finalize_model()
  }
  params_vignette_multi_outcome <- c(
    A.m = log(0.30),
    A.s = 0.18,
    B.m = log(0.38),
    B.s = 0.22
  )

  vignette_chained_onset <- function() {
    race_spec() |>
      add_accumulator("A", "lognormal") |>
      add_accumulator("B", "lognormal") |>
      add_accumulator("C", "lognormal", onset = after("B")) |>
      add_outcome("A", "A") |>
      add_outcome("C", "C") |>
      finalize_model()
  }
  params_vignette_chained_onset <- c(
    A.m = log(0.28),
    A.s = 0.14,
    B.m = log(0.1),
    B.s = 0.1,
    C.m = log(0.15),
    C.s = 0.1
  )

  # Bundle target models -------------------------------------------------------
  models <- list(
    vignette_multi_outcome = list(spec = vignette_multi_outcome, params = params_vignette_multi_outcome),
    vignette_chained_onset = list(spec = vignette_chained_onset, params = params_vignette_chained_onset)
  )

  # Test each model -----------------------------------------------------------
  snapshot_res <- list()
  for (name in names(models)) {
    mod <- models[[name]]
    spec_obj <- mod$spec()
    structure <- finalize_model(spec_obj)
    params_vec <- mod$params

    # Simulate with components retained
    params_df <- build_param_matrix(spec_obj, params_vec, n_trials = 500)
    data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)

    # Analytic probabilities (single-trial param set)
    analytic <- response_probabilities(
      structure,
      build_param_matrix(spec_obj, params_vec, n_trials = 1L),
      include_na = TRUE
    )

    # Observed response probabilities from simulation
    emp <- prop.table(table(data_df$R, useNA = "ifany"))
    emp_names <- names(emp)
    emp_names[is.na(emp_names)] <- "NA"
    emp <- as.numeric(emp)
    names(emp) <- emp_names

    # Non-marginalized likelihood with slim param matrix aligned to components
    prepared <- prepare_data(structure, data_df)
    ctx <- make_context(structure)
    params_df_slim <- build_param_matrix(
      spec_obj,
      params_vec,
      trial_df = prepared
    )
    ll <- as.numeric(log_likelihood(ctx, prepared, params_df_slim))

    # Order and round for stable snapshots
    analytic <- round(analytic[order(names(analytic))], 6)
    emp <- round(emp[order(names(emp))], 6)
    ll <- round(ll, 6)

    snapshot_res[[name]] <- list(
      observed = emp,
      analytic = analytic,
      loglik = ll
    )
  }

  testthat::expect_snapshot_value(snapshot_res, style = "json2")
})
