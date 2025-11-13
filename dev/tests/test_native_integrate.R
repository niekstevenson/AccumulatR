rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")

gaussian_fn <- function(x) exp(-x * x)
val_native <- .native_integrate(gaussian_fn, 0, 1)
val_stats <- stats::integrate(gaussian_fn, lower = 0, upper = 1, rel.tol = 1e-8)$value
if (abs(val_native - val_stats) > 1e-6) {
  stop(sprintf("Native integrator mismatch for Gaussian integral: %.8f vs %.8f", val_native, val_stats))
}

poly_fn <- function(x) x^3 - 2 * x + 1
val_poly <- .native_integrate(poly_fn, -2, 3)
expected_poly <- stats::integrate(poly_fn, lower = -2, upper = 3)$value
if (abs(val_poly - expected_poly) > 1e-8) {
  stop("Native integrator mismatch for polynomial integral")
}

cat("native_integrate: ok\n")
