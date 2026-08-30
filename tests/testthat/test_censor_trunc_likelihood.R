censor_trunc_race <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    test_separate_all_parameters() |>
    finalize_model()
}

testthat::test_that("response censor codes condition race regions on truncation", {
  model <- censor_trunc_race()
  data <- data.frame(
    trials = 1:4,
    R = c("A", NA, "A", NA),
    rt = c(0.4, NA, NA, NA),
    LT = 0.2,
    UT = 0.8,
    LC = c(0, 0.3, 0, 0.3),
    UC = c(Inf, Inf, 0.6, 0.6),
    missingness = c(NA, 1L, 2L, 3L)
  )
  prepared <- prepare_data(model, data)
  pars <- c(
    A.m = log(0.3), A.s = 0.16, A.t0 = 0,
    B.m = log(0.42), B.s = 0.2, B.t0 = 0
  )
  observed <- as.numeric(log_likelihood(
    make_context(model),
    prepared,
    build_param_matrix(model, pars, trial_df = prepared),
    sum = FALSE
  ))

  survivor <- function(t) {
    plnorm(t, pars[["A.m"]], pars[["A.s"]], lower.tail = FALSE) *
      plnorm(t, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE)
  }
  race_mass <- function(lower, upper) survivor(lower) - survivor(upper)
  a_mass <- integrate(
    function(t) dlnorm(t, pars[["A.m"]], pars[["A.s"]]) *
      plnorm(t, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE),
    0.6,
    0.8,
    rel.tol = 1e-12
  )$value
  normalizer <- race_mass(0.2, 0.8)
  expected <- log(c(
    dlnorm(0.4, pars[["A.m"]], pars[["A.s"]]) *
      plnorm(0.4, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE),
    race_mass(0.2, 0.3),
    a_mass,
    race_mass(0.2, 0.3) + race_mass(0.6, 0.8)
  ) / normalizer)

  testthat::expect_identical(unique(prepared$missingness), c(NA_integer_, 1:3))
  testthat::expect_equal(observed, expected, tolerance = 5e-7)
})

testthat::test_that("known censored responses combine every mapped outcome", {
  model <- race_spec() |>
    add_accumulator("left", "lognormal") |>
    add_accumulator("right", "lognormal") |>
    add_accumulator("timeout", "lognormal") |>
    add_outcome("Left", "left") |>
    add_outcome("Right", "right") |>
    add_outcome("TIMEOUT", "timeout", options = list(
      guess = list(
        labels = c("Left", "Right"),
        weights = c(0.2, 0.8),
        rt_policy = "keep"
      )
    )) |>
    test_separate_all_parameters() |>
    finalize_model()
  data <- prepare_data(model, data.frame(
    trials = 1L,
    R = "Right",
    rt = NA_real_,
    UC = 0.4,
    missingness = 2L
  ))
  pars <- c(
    left.m = log(0.3), left.s = 0.18, left.t0 = 0,
    right.m = log(0.36), right.s = 0.2, right.t0 = 0,
    timeout.m = log(0.45), timeout.s = 0.22, timeout.t0 = 0
  )
  observed <- as.numeric(log_likelihood(
    make_context(model),
    data,
    build_param_matrix(model, pars, trial_df = data)
  ))
  winner_density <- function(t, winner) {
    names <- c("left", "right", "timeout")
    stats::dlnorm(t, pars[[paste0(winner, ".m")]], pars[[paste0(winner, ".s")]]) *
      Reduce(`*`, lapply(setdiff(names, winner), function(other) {
        stats::plnorm(
          t,
          pars[[paste0(other, ".m")]],
          pars[[paste0(other, ".s")]],
          lower.tail = FALSE
        )
      }))
  }
  expected <- integrate(
    function(t) winner_density(t, "right") +
      0.8 * winner_density(t, "timeout"),
    0.4,
    30,
    rel.tol = 1e-12
  )$value

  testthat::expect_equal(observed, log(expected), tolerance = 1e-9)
})

testthat::test_that("unknown censoring and truncation use observable responses", {
  model <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("hidden", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    add_outcome("hidden", "hidden", options = list(
      map_outcome_to = NA_character_
    )) |>
    test_separate_all_parameters() |>
    finalize_model()
  data <- prepare_data(model, data.frame(
    trials = 1:2,
    R = c(NA, "A"),
    rt = c(NA, 0.35),
    LT = c(0, 0.2),
    UT = c(Inf, 0.8),
    LC = c(0.5, 0),
    missingness = c(1L, NA)
  ))
  pars <- c(
    A.m = log(0.3), A.s = 0.18, A.t0 = 0,
    B.m = log(0.4), B.s = 0.2, B.t0 = 0,
    hidden.m = log(0.34), hidden.s = 0.21, hidden.t0 = 0
  )
  observed <- as.numeric(log_likelihood(
    make_context(model),
    data,
    build_param_matrix(model, pars, trial_df = data),
    sum = FALSE
  ))
  observable_density <- function(t) {
    s_a <- plnorm(t, pars[["A.m"]], pars[["A.s"]], lower.tail = FALSE)
    s_b <- plnorm(t, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE)
    s_h <- plnorm(t, pars[["hidden.m"]], pars[["hidden.s"]], lower.tail = FALSE)
    dlnorm(t, pars[["A.m"]], pars[["A.s"]]) * s_b * s_h +
      dlnorm(t, pars[["B.m"]], pars[["B.s"]]) * s_a * s_h
  }
  numerator_1 <- integrate(observable_density, 0, 0.5, rel.tol = 1e-12)$value
  denominator_2 <- integrate(observable_density, 0.2, 0.8, rel.tol = 1e-12)$value
  numerator_2 <- dlnorm(0.35, pars[["A.m"]], pars[["A.s"]]) *
    plnorm(0.35, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE) *
    plnorm(0.35, pars[["hidden.m"]], pars[["hidden.s"]], lower.tail = FALSE)

  testthat::expect_equal(
    observed,
    c(log(numerator_1), log(numerator_2 / denominator_2)),
    tolerance = 1e-8
  )
})

testthat::test_that("only selected-outcome integration uses the EMC time limit", {
  model <- race_spec() |>
    add_accumulator("A", "lba") |>
    add_accumulator("B", "lba") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    test_separate_all_parameters() |>
    finalize_model()
  data <- prepare_data(model, data.frame(
    trials = 1:2,
    R = c("A", NA),
    rt = NA_real_,
    UC = 0.5,
    missingness = 2L
  ))
  pars <- c(
    A.v = 2, A.B = 1, A.A = 0.5, A.sv = 1, A.t0 = 0,
    B.v = 1.5, B.B = 1, B.A = 0.4, B.sv = 1, B.t0 = 0
  )
  observed <- as.numeric(log_likelihood(
    make_context(model),
    data,
    build_param_matrix(model, pars, trial_df = data),
    sum = FALSE
  ))
  expected <- c(
    integrate(
      function(t) {
        AccumulatR:::dist_lba_pdf(t, 2, 1, 0.5, 1) *
          (1 - AccumulatR:::dist_lba_cdf(t, 1.5, 1, 0.4, 1))
      },
      0.5,
      30,
      rel.tol = 1e-11
    )$value,
    (1 - AccumulatR:::dist_lba_cdf(0.5, 2, 1, 0.5, 1)) *
      (1 - AccumulatR:::dist_lba_cdf(0.5, 1.5, 1, 0.4, 1))
  )

  testthat::expect_equal(observed, log(expected), tolerance = 1e-9)
})

testthat::test_that("unknown-response tails retain completed guarded outcomes", {
  model <- race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("D", "lognormal") |>
    add_outcome("A", inhibit("A", by = "B")) |>
    add_outcome("D", "D") |>
    test_separate_all_parameters() |>
    finalize_model()
  data <- prepare_data(model, data.frame(
    trials = 1:2,
    R = NA_character_,
    rt = NA_real_,
    LC = c(0, 0.4),
    UC = c(0.4, Inf),
    missingness = 2:1
  ))
  pars <- c(
    A.m = log(0.3), A.s = 0.18, A.t0 = 0,
    B.m = log(0.36), B.s = 0.2, B.t0 = 0,
    D.m = log(0.5), D.s = 0.22, D.t0 = 0
  )
  observed <- as.numeric(log_likelihood(
    make_context(model),
    data,
    build_param_matrix(model, pars, trial_df = data),
    sum = FALSE
  ))
  a_completed <- integrate(
    function(t) dlnorm(t, pars[["A.m"]], pars[["A.s"]]) *
      plnorm(t, pars[["B.m"]], pars[["B.s"]], lower.tail = FALSE),
    0,
    0.4,
    rel.tol = 1e-12
  )$value
  survival <- plnorm(
    0.4, pars[["D.m"]], pars[["D.s"]], lower.tail = FALSE
  ) * (1 - a_completed)

  testthat::expect_equal(
    observed,
    c(log(survival), log1p(-survival)),
    tolerance = 1e-9
  )
})

testthat::test_that("prepare_data enforces response-level censoring codes", {
  model <- censor_trunc_race()
  row <- data.frame(trials = 1L, R = "A", rt = NA_real_)

  testthat::expect_error(
    prepare_data(model, transform(row, missingness = 0L)),
    "1, 2, 3, or NA"
  )
  testthat::expect_error(
    prepare_data(model, transform(row, rt = 0.2, LC = 0.3, missingness = 1L)),
    "rt = NA"
  )
  testthat::expect_error(
    prepare_data(model, transform(
      row, LT = 0.2, LC = 0.6, UC = 0.4, UT = 0.8, missingness = 3L
    )),
    "LC <= UC"
  )
  testthat::expect_no_error(prepare_data(model, transform(
    row, LT = 0.2, LC = 0.2, UT = 0.8, missingness = 1L
  )))
  testthat::expect_no_error(prepare_data(model, transform(
    row, LT = 0.2, UC = 0.8, UT = 0.8, missingness = 2L
  )))
  testthat::expect_no_error(prepare_data(model, data.frame(
    trials = 1L, R = "A", rt = 0.4, LT = NA, missingness = NA
  )))
  testthat::expect_error(
    prepare_data(model, transform(
      row, LT = 0.2, UT = 0.8, missingness = NA_integer_
    )),
    "finite first-rank"
  )
})
