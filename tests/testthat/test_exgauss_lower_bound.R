testthat::test_that("generic lower-bound helper renormalizes exgauss", {
  dist_code_exgauss <- 3L
  mu <- 0.18
  sigma <- 0.06
  tau <- 0.09
  lower <- 0.0
  x <- c(-0.05, 0.00, 0.14, Inf)

  base_pdf <- AccumulatR:::dist_exgauss_pdf(x, mu, sigma, tau)
  base_cdf <- AccumulatR:::dist_exgauss_cdf(x, mu, sigma, tau)
  cdf0 <- AccumulatR:::dist_exgauss_cdf(0, mu, sigma, tau)[[1]]
  z <- 1 - cdf0

  got_pdf <- AccumulatR:::dist_pdf_lower_bound(
    x, dist_code_exgauss, lower, mu, sigma, tau, 0
  )
  got_cdf <- AccumulatR:::dist_cdf_lower_bound(
    x, dist_code_exgauss, lower, mu, sigma, tau, 0
  )

  expected_pdf <- c(0, base_pdf[2] / z, base_pdf[3] / z, 0)
  expected_cdf <- c(0, 0, (base_cdf[3] - cdf0) / z, 1)

  testthat::expect_equal(got_pdf, expected_pdf, tolerance = 1e-12)
  testthat::expect_equal(got_cdf, expected_cdf, tolerance = 1e-12)
})

testthat::test_that("generic lower-bound helper is identity when family already respects the bound", {
  dist_code_lognormal <- 1L
  meanlog <- log(0.30)
  sdlog <- 0.20
  lower <- 0.0
  x <- c(-0.05, 0.00, 0.14, Inf)

  got_pdf <- AccumulatR:::dist_pdf_lower_bound(
    x, dist_code_lognormal, lower, meanlog, sdlog, 0, 0
  )
  got_cdf <- AccumulatR:::dist_cdf_lower_bound(
    x, dist_code_lognormal, lower, meanlog, sdlog, 0, 0
  )

  testthat::expect_equal(
    got_pdf,
    AccumulatR:::dist_lognormal_pdf(x, meanlog, sdlog),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    got_cdf,
    AccumulatR:::dist_lognormal_cdf(x, meanlog, sdlog),
    tolerance = 1e-12
  )
})

testthat::test_that("direct exgauss likelihood uses truncated renormalized t0 semantics", {
  spec <- race_spec() |>
    add_accumulator("a", "exgauss") |>
    add_outcome("A", "a")

  structure <- finalize_model(spec)
  params <- c(
    a.mu = 0.18,
    a.sigma = 0.06,
    a.tau = 0.09,
    a.q = 0.00,
    a.t0 = 0.20
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  rt_obs <- 0.35
  data_df <- data.frame(
    trial = 1L,
    R = factor("A", levels = "A"),
    rt = rt_obs
  )

  ll <- as.numeric(
    log_likelihood(build_likelihood_context(structure, data_df), data_df, params_df)
  )

  shifted <- rt_obs - params[["a.t0"]]
  cdf0 <- AccumulatR:::dist_exgauss_cdf(
    0, params[["a.mu"]], params[["a.sigma"]], params[["a.tau"]]
  )[[1]]
  z <- 1 - cdf0
  expected <- AccumulatR:::dist_exgauss_pdf(
    shifted, params[["a.mu"]], params[["a.sigma"]], params[["a.tau"]]
  )[[1]] / z

  testthat::expect_equal(ll, log(expected), tolerance = 1e-8)
})
