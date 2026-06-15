cens_trunc_simple_structure <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    test_separate_all_parameters() |>
    finalize_model()
}

cens_trunc_simple_params <- c(
  A.m = log(0.30), A.s = 0.16, A.t0 = 0,
  B.m = log(0.42), B.s = 0.20, B.t0 = 0
)

cens_trunc_simple_eval <- function(data_df) {
  structure <- cens_trunc_simple_structure()
  prepared <- prepare_data(structure, data_df)
  params <- build_param_matrix(
    structure$model_spec,
    cens_trunc_simple_params,
    trial_df = prepared
  )
  as.numeric(log_likelihood(make_context(structure), prepared, params))
}

cens_trunc_surv <- function(label, t) {
  if (!is.finite(t)) {
    return(0)
  }
  mean <- cens_trunc_simple_params[[paste0(label, ".m")]]
  sd <- cens_trunc_simple_params[[paste0(label, ".s")]]
  1 - stats::plnorm(t, mean, sd)
}

cens_trunc_density <- function(label, t) {
  mean <- cens_trunc_simple_params[[paste0(label, ".m")]]
  sd <- cens_trunc_simple_params[[paste0(label, ".s")]]
  stats::dlnorm(t, mean, sd)
}

cens_trunc_race_mass <- function(lower, upper) {
  cens_trunc_surv("A", lower) * cens_trunc_surv("B", lower) -
    cens_trunc_surv("A", upper) * cens_trunc_surv("B", upper)
}

cens_trunc_label_density <- function(label, t) {
  other <- if (identical(label, "A")) "B" else "A"
  cens_trunc_density(label, t) * cens_trunc_surv(other, t)
}

cens_trunc_label_mass <- function(label, lower, upper) {
  stats::integrate(
    function(x) vapply(x, function(xx) cens_trunc_label_density(label, xx), numeric(1)),
    lower,
    upper,
    rel.tol = 1e-10
  )$value
}

testthat::test_that("exact RT rows with repeated censoring cutoffs remain exact observations", {
  base <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = "A",
    rt = 0.35,
    stringsAsFactors = FALSE
  ))
  with_cutoffs <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = "A",
    rt = 0.35,
    LC = 0.10,
    UC = 0.80,
    stringsAsFactors = FALSE
  ))

  testthat::expect_equal(with_cutoffs, base, tolerance = 1e-10)
})

testthat::test_that("exact RT rows are normalized by truncation window mass", {
  out <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = "A",
    rt = 0.35,
    LT = 0.20,
    UT = 0.70,
    stringsAsFactors = FALSE
  ))
  expected <- log(
    cens_trunc_label_density("A", 0.35) /
      cens_trunc_race_mass(0.20, 0.70)
  )

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("unknown-label upper censoring integrates race tail", {
  out <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = NA_character_,
    rt = NA_real_,
    UC = 0.50,
    stringsAsFactors = FALSE
  ))
  expected <- log(cens_trunc_race_mass(0.50, Inf))

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("unknown-label lower censoring integrates early race mass", {
  out <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = NA_character_,
    rt = NA_real_,
    LC = 0.25,
    stringsAsFactors = FALSE
  ))
  expected <- log(cens_trunc_race_mass(0, 0.25))

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("label-specific interval censoring integrates defective outcome density", {
  out <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = "A",
    rt = NA_real_,
    UC = 0.20,
    LC = 0.60,
    stringsAsFactors = FALSE
  ))
  expected <- log(cens_trunc_label_mass("A", 0.20, 0.60))

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("missing RT with truncation only is a bounded interval observation", {
  out <- cens_trunc_simple_eval(data.frame(
    trials = 1L,
    R = "A",
    rt = NA_real_,
    LT = 0.20,
    UT = 0.60,
    stringsAsFactors = FALSE
  ))
  expected <- log(
    cens_trunc_label_mass("A", 0.20, 0.60) /
      cens_trunc_race_mass(0.20, 0.60)
  )

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("latent mixture truncation uses ratio of weighted sums", {
  structure <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("Left", "A") |>
    add_outcome("Right", "B") |>
    add_component("left_only", members = "A") |>
    add_component("right_only", members = "B") |>
    set_mixture(mode = "fixed", weights = c(left_only = 0.25, right_only = 0.75)) |>
    test_separate_all_parameters() |>
    finalize_model()
  data_df <- data.frame(
    trials = 1L,
    R = "Left",
    rt = 0.32,
    LT = 0.20,
    UT = 0.70,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_data(structure, data_df)
  params <- c(
    A.m = log(0.30), A.s = 0.16, A.t0 = 0,
    B.m = log(0.42), B.s = 0.20, B.t0 = 0
  )
  out <- as.numeric(log_likelihood(
    make_context(structure),
    prepared,
    build_param_matrix(structure$model_spec, params, trial_df = prepared)
  ))

  mass_a <- stats::plnorm(0.70, params[["A.m"]], params[["A.s"]]) -
    stats::plnorm(0.20, params[["A.m"]], params[["A.s"]])
  mass_b <- stats::plnorm(0.70, params[["B.m"]], params[["B.s"]]) -
    stats::plnorm(0.20, params[["B.m"]], params[["B.s"]])
  expected <- log(
    0.25 * stats::dlnorm(0.32, params[["A.m"]], params[["A.s"]]) /
      (0.25 * mass_a + 0.75 * mass_b)
  )

  testthat::expect_equal(out, expected, tolerance = 2e-3)
})

testthat::test_that("ranked observations reject censoring and truncation columns", {
  structure <- cens_trunc_simple_structure()
  testthat::expect_error(
    prepare_data(
      structure,
      data.frame(
        trials = 1L,
        R = "A",
        rt = 0.30,
        R2 = "B",
        rt2 = 0.50,
        UT = 0.80,
        stringsAsFactors = FALSE
      )
    ),
    "ranked observations do not support censoring or truncation"
  )
})
