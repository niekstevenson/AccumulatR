rm(list = ls())

base::source("R/dist.R")

approx_equal <- function(a, b, tol = 1e-10) {
  diff <- abs(a - b)
  diff <- diff[is.finite(diff)]
  if (length(diff) == 0L) return(TRUE)
  max(diff) <= tol
}

# ------------------------------------------------------------------------------
# Lognormal
# ------------------------------------------------------------------------------

lognormal_par <- list(meanlog = 0.3, sdlog = 0.8)
lognormal_x <- c(0, 0.2, 0.5, 1.0, 2.5, NA_real_, Inf, -Inf)

lognormal_reg <- dist_registry("lognormal")

lognormal_pdf_cpp <- lognormal_reg$d(lognormal_x, lognormal_par)
lognormal_pdf_ref <- dlnorm(lognormal_x, meanlog = lognormal_par$meanlog, sdlog = lognormal_par$sdlog)
if (!approx_equal(lognormal_pdf_cpp, lognormal_pdf_ref)) {
  stop("lognormal pdf mismatch between C++ and reference")
}

lognormal_cdf_cpp <- lognormal_reg$p(lognormal_x, lognormal_par)
lognormal_cdf_ref <- plnorm(lognormal_x, meanlog = lognormal_par$meanlog, sdlog = lognormal_par$sdlog)
if (!approx_equal(lognormal_cdf_cpp, lognormal_cdf_ref)) {
  stop("lognormal cdf mismatch between C++ and reference")
}

set.seed(2025)
lognormal_rng_cpp <- lognormal_reg$r(1000, lognormal_par)
set.seed(2025)
lognormal_rng_ref <- rlnorm(1000, meanlog = lognormal_par$meanlog, sdlog = lognormal_par$sdlog)
if (!isTRUE(all.equal(lognormal_rng_cpp, lognormal_rng_ref, tolerance = 1e-12))) {
  stop("lognormal rng sequence differs from reference")
}

if (!all(is.na(lognormal_reg$d(1.0, list(meanlog = 0, sdlog = -1))))) {
  stop("lognormal pdf should return NA for invalid sdlog")
}

# ------------------------------------------------------------------------------
# Gamma
# ------------------------------------------------------------------------------

gamma_par <- list(shape = 2.3, rate = 1.7)
gamma_x <- seq(0, 5, length.out = 11)
gamma_reg <- dist_registry("gamma")

gamma_pdf_cpp <- gamma_reg$d(gamma_x, gamma_par)
gamma_pdf_ref <- dgamma(gamma_x, shape = gamma_par$shape, rate = gamma_par$rate)
if (!approx_equal(gamma_pdf_cpp, gamma_pdf_ref)) {
  stop("gamma pdf mismatch between C++ and reference")
}

gamma_cdf_cpp <- gamma_reg$p(gamma_x, gamma_par)
gamma_cdf_ref <- pgamma(gamma_x, shape = gamma_par$shape, rate = gamma_par$rate)
if (!approx_equal(gamma_cdf_cpp, gamma_cdf_ref)) {
  stop("gamma cdf mismatch between C++ and reference")
}

set.seed(88)
gamma_rng_cpp <- gamma_reg$r(2000, gamma_par)
set.seed(88)
gamma_rng_ref <- rgamma(2000, shape = gamma_par$shape, rate = gamma_par$rate)
if (!isTRUE(all.equal(gamma_rng_cpp, gamma_rng_ref, tolerance = 1e-12))) {
  stop("gamma rng sequence differs from reference")
}

if (!all(is.na(gamma_reg$d(1.0, list(shape = -1, rate = 1))))) {
  stop("gamma pdf should return NA for invalid shape")
}

# ------------------------------------------------------------------------------
# Ex-Gaussian
# ------------------------------------------------------------------------------

exgauss_pdf_ref <- function(x, par) {
  mu <- par$mu; sigma <- par$sigma; tau <- par$tau
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[[i]]
    if (!is.finite(xi)) {
      out[[i]] <- 0.0
      next
    }
    if (!is.finite(sigma) || !is.finite(tau) || sigma <= 0 || tau <= 0) {
      out[[i]] <- NA_real_
      next
    }
    z <- (xi - mu) / sigma
    w <- sigma / tau
    out[[i]] <- (1.0 / tau) *
      exp((sigma * sigma) / (2.0 * tau * tau) - (xi - mu) / tau) *
      pnorm(z - w)
  }
  out
}

exgauss_cdf_ref <- function(x, par) {
  mu <- par$mu; sigma <- par$sigma; tau <- par$tau
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[[i]]
    if (is.na(xi)) {
      out[[i]] <- NA_real_
      next
    }
    if (!is.finite(xi)) {
      out[[i]] <- if (xi < 0) 0.0 else 1.0
      next
    }
    if (!is.finite(sigma) || !is.finite(tau) || sigma <= 0 || tau <= 0) {
      out[[i]] <- NA_real_
      next
    }
    z <- (xi - mu) / sigma
    w <- sigma / tau
    out[[i]] <- pnorm(z) - exp((sigma * sigma) / (2.0 * tau * tau) - (xi - mu) / tau) * pnorm(z - w)
  }
  out
}

exgauss_par <- list(mu = 0.4, sigma = 0.6, tau = 0.9)
exgauss_x <- seq(-1, 5, length.out = 13)
exgauss_reg <- dist_registry("exgauss")

exgauss_pdf_cpp <- exgauss_reg$d(exgauss_x, exgauss_par)
exgauss_pdf_calc <- exgauss_pdf_ref(exgauss_x, exgauss_par)
if (!approx_equal(exgauss_pdf_cpp, exgauss_pdf_calc)) {
  stop("exgaussian pdf mismatch between C++ and reference")
}

exgauss_cdf_cpp <- exgauss_reg$p(exgauss_x, exgauss_par)
exgauss_cdf_calc <- exgauss_cdf_ref(exgauss_x, exgauss_par)
if (!approx_equal(exgauss_cdf_cpp, exgauss_cdf_calc)) {
  stop("exgaussian cdf mismatch between C++ and reference")
}

set.seed(321)
exgauss_rng_cpp <- exgauss_reg$r(2000, exgauss_par)
set.seed(321)
exgauss_rng_ref <- rnorm(2000, mean = exgauss_par$mu, sd = exgauss_par$sigma) +
  rexp(2000, rate = 1.0 / exgauss_par$tau)
if (!isTRUE(all.equal(exgauss_rng_cpp, exgauss_rng_ref, tolerance = 1e-12))) {
  stop("exgaussian rng sequence differs from reference")
}

if (!all(is.na(exgauss_reg$d(1.0, list(mu = 0, sigma = -1, tau = 1))))) {
  stop("exgaussian pdf should return NA for invalid sigma")
}

cat("distributions_cpp: ok\n")
