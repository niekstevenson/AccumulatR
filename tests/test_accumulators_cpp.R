rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_primitives.R")

approx_equal <- function(a, b, tol = 1e-10) {
  isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
}

make_acc <- function(dist, params, onset = 0, q = 0, t0 = NULL) {
  if (is.null(params)) params <- list()
  if (is.null(t0)) {
    if (is.null(params$t0) || length(params$t0) == 0L) {
      params$t0 <- 0
    } else {
      params$t0 <- as.numeric(params$t0)[1]
    }
  } else {
    params$t0 <- as.numeric(t0)[1]
  }
  list(
    dist = dist,
    params = params,
    onset = onset,
    q = q
  )
}

# ------------------------------------------------------------------------------
# Lognormal accumulator
# ------------------------------------------------------------------------------

acc_ln <- make_acc("lognormal", list(meanlog = 0.2, sdlog = 0.5), onset = 0.3, q = 0.1, t0 = 0.12)
t_ln <- 1.5
t_ln_adj <- t_ln - (acc_ln$onset + acc_ln$params$t0)
expected_pdf_ln <- dlnorm(t_ln_adj, meanlog = acc_ln$params$meanlog, sdlog = acc_ln$params$sdlog)
expected_cdf_ln <- plnorm(t_ln_adj, meanlog = acc_ln$params$meanlog, sdlog = acc_ln$params$sdlog)

if (!approx_equal(.acc_density(acc_ln, t_ln), (1 - acc_ln$q) * expected_pdf_ln)) {
  stop("lognormal accumulator density mismatch")
}

if (!approx_equal(.acc_density_success(acc_ln, t_ln), expected_pdf_ln)) {
  stop("lognormal accumulator success density mismatch")
}

if (!approx_equal(.acc_cdf_success(acc_ln, t_ln), expected_cdf_ln)) {
  stop("lognormal accumulator success CDF mismatch")
}

expected_surv_ln <- acc_ln$q + (1 - acc_ln$q) * (1 - expected_cdf_ln)
if (!approx_equal(.acc_survival(acc_ln, t_ln), expected_surv_ln)) {
  stop("lognormal accumulator survival mismatch")
}

vec_t <- acc_ln$onset + acc_ln$params$t0 + c(0.2, 0.7, 1.3)
vec_expected <- (1 - acc_ln$q) * dlnorm(vec_t - (acc_ln$onset + acc_ln$params$t0),
                                       meanlog = acc_ln$params$meanlog,
                                       sdlog = acc_ln$params$sdlog)
if (!approx_equal(.acc_density(acc_ln, vec_t), vec_expected)) {
  stop("lognormal accumulator vector density mismatch")
}

t_before <- max(0, acc_ln$onset + acc_ln$params$t0 - 0.05)
if (!approx_equal(.acc_density(acc_ln, t_before), 0.0)) {
  stop("lognormal accumulator early density should be zero")
}
if (!approx_equal(.acc_density_success(acc_ln, t_before), 0.0)) {
  stop("lognormal success density before onset should be zero")
}
if (!approx_equal(.acc_cdf_success(acc_ln, t_before), 0.0)) {
  stop("lognormal success CDF before onset should be zero")
}
if (!approx_equal(.acc_survival(acc_ln, t_before), 1.0)) {
  stop("lognormal survival before onset should be one")
}

if (!approx_equal(.acc_density(acc_ln, Inf), 0.0)) {
  stop("lognormal density at Inf should be zero")
}
if (!approx_equal(.acc_density_success(acc_ln, Inf), 0.0)) {
  stop("lognormal success density at Inf should be zero")
}
if (!approx_equal(.acc_cdf_success(acc_ln, Inf), 1.0)) {
  stop("lognormal success CDF at Inf should be one")
}
if (!approx_equal(.acc_survival(acc_ln, Inf), 0.0)) {
  stop("lognormal survival at Inf should be zero")
}

acc_fail <- make_acc("lognormal", acc_ln$params, onset = acc_ln$onset, q = 1.0)
if (!approx_equal(.acc_density(acc_fail, t_ln), 0.0)) {
  stop("Accumulator with q=1 should have zero density")
}
if (!approx_equal(.acc_survival(acc_fail, t_ln), 1.0)) {
  stop("Accumulator with q=1 should have survival one for finite t")
}

# ------------------------------------------------------------------------------
# Gamma accumulator
# ------------------------------------------------------------------------------

acc_gam <- make_acc("gamma", list(shape = 2.3, rate = 1.2), onset = 0.1, q = 0.25, t0 = 0.05)
t_gam <- 1.8
t_gam_adj <- t_gam - (acc_gam$onset + acc_gam$params$t0)
expected_pdf_gam <- dgamma(t_gam_adj, shape = acc_gam$params$shape, rate = acc_gam$params$rate)
expected_cdf_gam <- pgamma(t_gam_adj, shape = acc_gam$params$shape, rate = acc_gam$params$rate)

if (!approx_equal(.acc_density(acc_gam, t_gam), (1 - acc_gam$q) * expected_pdf_gam)) {
  stop("gamma accumulator density mismatch")
}
if (!approx_equal(.acc_density_success(acc_gam, t_gam), expected_pdf_gam)) {
  stop("gamma accumulator success density mismatch")
}
if (!approx_equal(.acc_cdf_success(acc_gam, t_gam), expected_cdf_gam)) {
  stop("gamma accumulator success CDF mismatch")
}

expected_surv_gam <- acc_gam$q + (1 - acc_gam$q) * (1 - expected_cdf_gam)
if (!approx_equal(.acc_survival(acc_gam, t_gam), expected_surv_gam)) {
  stop("gamma accumulator survival mismatch")
}

# ------------------------------------------------------------------------------
# Ex-Gaussian accumulator
# ------------------------------------------------------------------------------

acc_ex <- make_acc("exgauss", list(mu = 0.5, sigma = 0.6, tau = 0.9), onset = 0.2, q = 0.05, t0 = 0.08)
t_ex <- 1.7
t_ex_adj <- t_ex - (acc_ex$onset + acc_ex$params$t0)
expected_pdf_ex <- {
  mu <- acc_ex$params$mu
  sigma <- acc_ex$params$sigma
  tau <- acc_ex$params$tau
  (1 / tau) * exp((sigma * sigma) / (2 * tau * tau) - (t_ex_adj - mu) / tau) *
    pnorm((t_ex_adj - mu) / sigma - sigma / tau)
}
expected_cdf_ex <- {
  mu <- acc_ex$params$mu
  sigma <- acc_ex$params$sigma
  tau <- acc_ex$params$tau
  z <- (t_ex_adj - mu) / sigma
  w <- sigma / tau
  pnorm(z) - exp((sigma * sigma) / (2 * tau * tau) - (t_ex_adj - mu) / tau) * pnorm(z - w)
}

if (!approx_equal(.acc_density(acc_ex, t_ex), (1 - acc_ex$q) * expected_pdf_ex)) {
  stop("exgaussian accumulator density mismatch")
}
if (!approx_equal(.acc_density_success(acc_ex, t_ex), expected_pdf_ex)) {
  stop("exgaussian accumulator success density mismatch")
}
if (!approx_equal(.acc_cdf_success(acc_ex, t_ex), expected_cdf_ex)) {
  stop("exgaussian accumulator success CDF mismatch")
}

# ------------------------------------------------------------------------------
# Invalid parameter handling
# ------------------------------------------------------------------------------

acc_bad <- make_acc("gamma", list(shape = -1, rate = 1), onset = 0, q = 0)
bad_val <- .acc_density_success(acc_bad, 1.0)
if (!is.na(bad_val)) {
  stop("Invalid gamma parameters should yield NA density")
}

cat("accumulators_cpp: ok\n")
