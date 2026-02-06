#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rcpp)
})

sourceCpp("dev/exgauss_verify.cpp")

# Reference formula using base R pnorm.
d_exgauss_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(rep(NA_real_, length(x)))
  }
  out <- numeric(length(x))
  is_na <- is.na(x)
  out[is_na] <- NA_real_
  is_finite <- is.finite(x)
  idx <- which(!is_na & is_finite)
  if (length(idx)) {
    xx <- x[idx]
    z <- (xx - mu) / sigma
    exponent <- (sigma^2) / (2 * tau^2) - (xx - mu) / tau
    val <- (1 / tau) * exp(exponent) * pnorm(z - sigma / tau)
    val[!is.finite(val) | val < 0] <- 0
    out[idx] <- val
  }
  out[!is_na & !is_finite] <- 0
  out
}

p_exgauss_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(rep(NA_real_, length(x)))
  }
  out <- numeric(length(x))
  is_na <- is.na(x)
  out[is_na] <- NA_real_
  is_finite <- is.finite(x)
  idx <- which(!is_na & is_finite)
  if (length(idx)) {
    xx <- x[idx]
    z <- (xx - mu) / sigma
    exponent <- (sigma^2) / (2 * tau^2) - (xx - mu) / tau
    val <- pnorm(z) - exp(exponent) * pnorm(z - sigma / tau)
    val[!is.finite(val)] <- NA_real_
    val <- pmin(pmax(val, 0), 1)
    out[idx] <- val
  }
  nonfinite_idx <- !is_na & !is_finite
  if (any(nonfinite_idx)) {
    out[nonfinite_idx & x < 0] <- 0
    out[nonfinite_idx & x > 0] <- 1
  }
  out
}

r_exgauss_ref <- function(n, mu, sigma, tau) {
  if (n <= 0 || !is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(numeric(0))
  }
  rnorm(n, mu, sigma) + rexp(n, rate = 1 / tau)
}

find_backend <- function() {
  pkgs <- c("exgauss", "retimes", "VGAM")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      next
    }
    exports <- getNamespaceExports(pkg)
    if (all(c("dexgauss", "pexgauss", "rexgauss") %in% exports)) {
      return(pkg)
    }
  }
  NULL
}

backend_pkg <- find_backend()
backend <- NULL
if (!is.null(backend_pkg)) {
  ns <- asNamespace(backend_pkg)
  backend <- list(
    pkg = backend_pkg,
    d = get("dexgauss", ns),
    p = get("pexgauss", ns),
    r = get("rexgauss", ns)
  )
}

compare_vec <- function(a, b) {
  na_mismatch <- sum(is.na(a) != is.na(b))
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) {
    return(list(max_abs = NA_real_, max_rel = NA_real_, na_mismatch = na_mismatch))
  }
  abs_err <- abs(a[ok] - b[ok])
  rel_err <- abs_err / pmax(1e-12, abs(b[ok]))
  list(max_abs = max(abs_err), max_rel = max(rel_err), na_mismatch = na_mismatch)
}

make_x <- function(mu, sigma, tau) {
  set.seed(1)
  grid <- seq(mu - 6 * sigma, mu + 10 * tau + 6 * sigma, length.out = 241)
  draws <- rnorm(200, mu, sigma) + rexp(200, rate = 1 / tau)
  specials <- mu + c(-100, -10, -5, -1, 0, 1, 5, 10, 100) * sigma
  specials <- c(specials, mu + c(-10, -5, -1, 0, 1, 5, 10) * tau)
  unique(sort(c(grid, draws, specials, -Inf, Inf, NA_real_)))
}

params <- expand.grid(
  mu = c(-5, 0, 2, 10),
  sigma = c(0.1, 0.5, 2, 5),
  tau = c(0.05, 0.5, 2, 10)
)

summary_rows <- vector("list", nrow(params))

for (i in seq_len(nrow(params))) {
  mu <- params$mu[i]
  sigma <- params$sigma[i]
  tau <- params$tau[i]
  x <- make_x(mu, sigma, tau)

  cpp_pdf <- exgauss_pdf_cpp(x, mu, sigma, tau)
  cpp_cdf <- exgauss_cdf_cpp(x, mu, sigma, tau)
  ref_pdf <- d_exgauss_ref(x, mu, sigma, tau)
  ref_cdf <- p_exgauss_ref(x, mu, sigma, tau)

  pdf_stats <- compare_vec(cpp_pdf, ref_pdf)
  cdf_stats <- compare_vec(cpp_cdf, ref_cdf)

  finite_cdf <- is.finite(cpp_cdf)
  cdf_nonmono <- FALSE
  if (any(finite_cdf)) {
    cdf_vals <- cpp_cdf[finite_cdf]
    if (length(cdf_vals) > 1) {
      cdf_nonmono <- any(diff(cdf_vals) < -1e-12)
    }
  }

  pdf_negative <- sum(is.finite(cpp_pdf) & cpp_pdf < -1e-12)
  cdf_outside <- sum(is.finite(cpp_cdf) & (cpp_cdf < -1e-12 | cpp_cdf > 1 + 1e-12))

  summary_rows[[i]] <- data.frame(
    mu = mu,
    sigma = sigma,
    tau = tau,
    pdf_max_abs = pdf_stats$max_abs,
    pdf_max_rel = pdf_stats$max_rel,
    pdf_na_mismatch = pdf_stats$na_mismatch,
    cdf_max_abs = cdf_stats$max_abs,
    cdf_max_rel = cdf_stats$max_rel,
    cdf_na_mismatch = cdf_stats$na_mismatch,
    pdf_negative = pdf_negative,
    cdf_outside = cdf_outside,
    cdf_nonmono = cdf_nonmono
  )
}

summary_df <- do.call(rbind, summary_rows)

cat("Ex-Gaussian verification vs base R formula\n")
cat(sprintf("Parameter grid: %d combinations\n", nrow(params)))
cat(sprintf("Max pdf abs err: %.3e\n", max(summary_df$pdf_max_abs, na.rm = TRUE)))
cat(sprintf("Max pdf rel err: %.3e\n", max(summary_df$pdf_max_rel, na.rm = TRUE)))
cat(sprintf("Max cdf abs err: %.3e\n", max(summary_df$cdf_max_abs, na.rm = TRUE)))
cat(sprintf("Max cdf rel err: %.3e\n", max(summary_df$cdf_max_rel, na.rm = TRUE)))
cat(sprintf("Any pdf negative: %s\n", any(summary_df$pdf_negative > 0)))
cat(sprintf("Any cdf outside [0,1]: %s\n", any(summary_df$cdf_outside > 0)))
cat(sprintf("Any cdf non-monotone: %s\n", any(summary_df$cdf_nonmono)))
cat(sprintf("Any pdf NA mismatch: %s\n", any(summary_df$pdf_na_mismatch > 0)))
cat(sprintf("Any cdf NA mismatch: %s\n", any(summary_df$cdf_na_mismatch > 0)))

# RNG checks
rng_params <- params[c(1, 10, 25, 40), ]
identical_rng <- logical(nrow(rng_params))
for (i in seq_len(nrow(rng_params))) {
  mu <- rng_params$mu[i]
  sigma <- rng_params$sigma[i]
  tau <- rng_params$tau[i]
  set.seed(123)
  cpp_draws <- exgauss_rng_cpp(1000, mu, sigma, tau)
  set.seed(123)
  ref_draws <- r_exgauss_ref(1000, mu, sigma, tau)
  identical_rng[i] <- isTRUE(all.equal(cpp_draws, ref_draws, tolerance = 0))
}
cat(sprintf("RNG exact match to rnorm+rexp (seeded): %s\n", all(identical_rng)))

# Optional backend comparison (if installed)
if (!is.null(backend)) {
  cat(sprintf("\nFound backend package: %s\n", backend$pkg))
  mu <- 0
  sigma <- 1
  tau <- 2
  x <- make_x(mu, sigma, tau)
  backend_pdf <- backend$d(x, mu, sigma, tau)
  backend_cdf <- backend$p(x, mu, sigma, tau)
  bpdf_stats <- compare_vec(backend_pdf, d_exgauss_ref(x, mu, sigma, tau))
  bcdf_stats <- compare_vec(backend_cdf, p_exgauss_ref(x, mu, sigma, tau))
  cat(sprintf("Backend vs formula pdf max abs err: %.3e\n", bpdf_stats$max_abs))
  cat(sprintf("Backend vs formula cdf max abs err: %.3e\n", bcdf_stats$max_abs))
  if (is.finite(bpdf_stats$max_abs) && bpdf_stats$max_abs > 1e-6) {
    cat("Backend parameterization may differ; review backend docs.\n")
  }
} else {
  cat("\nNo dexgauss/pexgauss/rexgauss backend found among installed packages.\n")
}
