.dist_param_scalar <- function(par, name, default = NULL, required = TRUE) {
  val <- par[[name]]
  if (is.null(val) || length(val) == 0L) {
    if (!isTRUE(required)) {
      return(default)
    }
    stop(sprintf("Missing parameter '%s' for distribution", name), call. = FALSE)
  }
  as.numeric(val)[1]
}

.dist_param_t0 <- function(par) {
  .dist_param_scalar(par, "t0", default = 0, required = FALSE)
}

.dist_zero_before_t0 <- function(values, shifted) {
  if (length(values) == 0L) return(values)
  mask <- !is.na(shifted) & shifted <= 0 & !is.na(values)
  if (any(mask)) values[mask] <- 0
  values
}

.dist_as_count <- function(n) {
  if (length(n) == 0L) return(0L)
  n_val <- as.numeric(n)[1]
  if (is.na(n_val) || n_val < 0) {
    stop("Sample size 'n' must be a non-negative numeric value", call. = FALSE)
  }
  as.integer(n_val)
}

.dist_param_values <- function(par, dist) {
  param_names <- dist_param_names(dist)
  if (length(param_names) == 0L) {
    stop(sprintf("Unknown distribution '%s'", dist), call. = FALSE)
  }
  values <- lapply(param_names, function(name) .dist_param_scalar(par, name))
  names(values) <- param_names
  values
}

.dist_shifted_rng <- function(fun, n, par, dist) {
  t0 <- .dist_param_t0(par)
  args <- .dist_param_values(par, dist)
  do.call(fun, c(list(.dist_as_count(n)), args)) + t0
}

.dist_shifted_eval <- function(fun, x, par, dist) {
  t0 <- .dist_param_t0(par)
  shifted <- x - t0
  args <- .dist_param_values(par, dist)
  out <- do.call(fun, c(list(shifted), args))
  .dist_zero_before_t0(out, shifted)
}

.dist_vapply_numeric <- function(x, fun) {
  vapply(x, fun, numeric(1), USE.NAMES = FALSE)
}

dist_lognormal_pdf <- function(x, m, s) {
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(0)
    stats::dlnorm(xi, meanlog = m, sdlog = s)
  })
}

dist_lognormal_cdf <- function(x, m, s) {
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(if (xi < 0) 0 else 1)
    min(max(stats::plnorm(xi, meanlog = m, sdlog = s), 0), 1)
  })
}

dist_lognormal_rng <- function(n, m, s) {
  n <- .dist_as_count(n)
  if (n == 0L || !is.finite(s) || s <= 0) return(numeric(0))
  exp(stats::rnorm(n, mean = m, sd = s))
}

dist_gamma_pdf <- function(x, shape, rate) {
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    return(rep(NA_real_, length(x)))
  }
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(0)
    stats::dgamma(xi, shape = shape, rate = rate)
  })
}

dist_gamma_cdf <- function(x, shape, rate) {
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    return(rep(NA_real_, length(x)))
  }
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(if (xi < 0) 0 else 1)
    min(max(stats::pgamma(xi, shape = shape, rate = rate), 0), 1)
  })
}

dist_gamma_rng <- function(n, shape, rate) {
  n <- .dist_as_count(n)
  if (n == 0L || !is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    return(numeric(0))
  }
  stats::rgamma(n, shape = shape, rate = rate)
}

dist_exgauss_pdf <- function(x, mu, sigma, tau) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(rep(NA_real_, length(x)))
  }
  inv_tau <- 1 / tau
  sigma_sq <- sigma * sigma
  tau_sq <- tau * tau
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(0)
    z <- (xi - mu) / sigma
    exponent <- sigma_sq / (2 * tau_sq) - (xi - mu) * inv_tau
    tail <- stats::pnorm(z - sigma * inv_tau)
    out <- inv_tau * exp(exponent) * tail
    if (!is.finite(out) || out <= 0) 0 else out
  })
}

dist_exgauss_cdf <- function(x, mu, sigma, tau) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(rep(NA_real_, length(x)))
  }
  inv_tau <- 1 / tau
  sigma_sq <- sigma * sigma
  tau_sq <- tau * tau
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(if (xi < 0) 0 else 1)
    z <- (xi - mu) / sigma
    exponent <- sigma_sq / (2 * tau_sq) - (xi - mu) * inv_tau
    tail <- stats::pnorm(z - sigma * inv_tau)
    base <- stats::pnorm(z)
    min(max(base - exp(exponent) * tail, 0), 1)
  })
}

dist_exgauss_rng <- function(n, mu, sigma, tau) {
  n <- .dist_as_count(n)
  if (n == 0L || !is.finite(sigma) || sigma <= 0 || !is.finite(tau) || tau <= 0) {
    return(numeric(0))
  }
  stats::rnorm(n, mean = mu, sd = sigma) + stats::rexp(n, rate = 1 / tau)
}

.dist_lba_denom <- function(v, sv) {
  denom <- stats::pnorm(v / sv)
  if (!is.finite(denom) || denom < 1e-10) denom <- 1e-10
  denom
}

dist_lba_pdf <- function(x, v, B, A, sv) {
  if (!is.finite(v) || !is.finite(B) || !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(rep(NA_real_, length(x)))
  }
  denom <- .dist_lba_denom(v, sv)
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi) || xi <= 0) return(0)
    if (A > 1e-10) {
      zs <- xi * sv
      cmz <- B - xi * v
      cz <- cmz / zs
      cz_max <- (cmz - A) / zs
      pdf <- (v * (stats::pnorm(cz) - stats::pnorm(cz_max)) +
        sv * (stats::dnorm(cz_max) - stats::dnorm(cz))) / (A * denom)
    } else {
      pdf <- stats::dnorm(B / xi, mean = v, sd = sv) * B / (xi * xi * denom)
    }
    if (!is.finite(pdf) || pdf <= 0) 0 else pdf
  })
}

dist_lba_cdf <- function(x, v, B, A, sv) {
  if (!is.finite(v) || !is.finite(B) || !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(rep(NA_real_, length(x)))
  }
  denom <- .dist_lba_denom(v, sv)
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(if (xi < 0) 0 else 1)
    if (xi <= 0) return(0)
    if (A > 1e-10) {
      zs <- xi * sv
      cmz <- B - xi * v
      xx <- cmz - A
      cz <- cmz / zs
      cz_max <- xx / zs
      cdf <- (1 +
        (zs * (stats::dnorm(cz_max) - stats::dnorm(cz)) +
           xx * stats::pnorm(cz_max) -
           cmz * stats::pnorm(cz)) / A) / denom
    } else {
      cdf <- stats::pnorm(B / xi, mean = v, sd = sv, lower.tail = FALSE) / denom
    }
    min(max(cdf, 0), 1)
  })
}

dist_lba_rng <- function(n, v, B, A, sv) {
  n <- .dist_as_count(n)
  if (n == 0L || !is.finite(v) || !is.finite(B) || !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(numeric(0))
  }
  lower_mass <- stats::pnorm(0, mean = v, sd = sv)
  upper_mass <- 1 - lower_mass
  if (!is.finite(upper_mass) || upper_mass <= 0) {
    return(rep(Inf, n))
  }
  out <- numeric(n)
  for (i in seq_len(n)) {
    start <- B - A * stats::runif(1)
    u <- stats::runif(1)
    q <- min(1 - 1e-15, lower_mass + u * upper_mass)
    drift <- stats::qnorm(q, mean = v, sd = sv)
    dt <- start / drift
    out[i] <- if (is.finite(dt) && dt > 0) dt else Inf
  }
  out
}

.rdm_pigt0_ref <- function(x, k, l) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l)) return(0)
  if (abs(l) < 1e-12) {
    z <- k / sqrt(x)
    return(min(max(2 * (1 - stats::pnorm(z)), 0), 1))
  }
  mu <- k / l
  lambda <- k * k
  p1 <- 1 - stats::pnorm(sqrt(lambda / x) * (1 + x / mu))
  p2 <- 1 - stats::pnorm(sqrt(lambda / x) * (1 - x / mu))
  part <- exp(log(2 * lambda / mu) + log(max(1e-300, p1)))
  min(max(part + p2, 0), 1)
}

.rdm_digt0_ref <- function(x, k, l) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l)) return(0)
  lambda <- k * k
  exponent <- if (l == 0) {
    -0.5 * lambda / x
  } else {
    mu <- k / l
    -(lambda / (2 * x)) * ((x * x) / (mu * mu) - 2 * x / mu + 1)
  }
  exp(exponent + 0.5 * log(lambda) - 0.5 * log(2 * x * x * x * pi))
}

.rdm_pigt_ref <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) || !is.finite(a)) return(0)
  if (a < threshold) return(.rdm_pigt0_ref(x, k, l))
  sqt <- sqrt(x)
  lgt <- log(x)
  if (l < threshold) {
    t5a <- 2 * stats::pnorm((k + a) / sqt) - 1
    t5b <- 2 * stats::pnorm((-k - a) / sqt) - 1
    t6a <- -0.5 * ((k + a) * (k + a) / x - log(2) - log(pi) + lgt) - log(a)
    t6b <- -0.5 * ((k - a) * (k - a) / x - log(2) - log(pi) + lgt) - log(a)
    cdf <- 1 + exp(t6a) - exp(t6b) + ((-k + a) * t5a - (k - a) * t5b) / (2 * a)
  } else {
    t1a <- exp(-0.5 * (k - a - x * l)^2 / x)
    t1b <- exp(-0.5 * (a + k - x * l)^2 / x)
    t1 <- exp(0.5 * (lgt - log(2) - log(pi))) * (t1a - t1b)
    t2a <- exp(2 * l * (k - a) + stats::pnorm(-(k - a + x * l) / sqt, log.p = TRUE))
    t2b <- exp(2 * l * (k + a) + stats::pnorm(-(k + a + x * l) / sqt, log.p = TRUE))
    t2 <- a + (t2b - t2a) / (2 * l)
    t4a <- 2 * stats::pnorm((k + a) / sqt - sqt * l) - 1
    t4b <- 2 * stats::pnorm((k - a) / sqt - sqt * l) - 1
    t4 <- 0.5 * (x * l - a - k + 0.5 / l) * t4a +
      0.5 * (k - a - x * l - 0.5 / l) * t4b
    cdf <- 0.5 * (t4 + t2 + t1) / a
  }
  if (!is.finite(cdf) || cdf < 0) return(0)
  min(max(cdf, 0), 1)
}

.rdm_digt_ref <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) || !is.finite(a)) return(0)
  if (a < threshold) return(.rdm_digt0_ref(x, k, l))
  if (l < threshold) {
    term <- exp(-(k - a)^2 / (2 * x)) - exp(-(k + a)^2 / (2 * x))
    pdf <- exp(-0.5 * (log(2) + log(pi) + log(x)) +
      log(max(1e-300, term)) - log(2) - log(a))
  } else {
    sqt <- sqrt(x)
    t1a <- -(a - k + x * l)^2 / (2 * x)
    t1b <- -(a + k - x * l)^2 / (2 * x)
    t1 <- (1 / sqrt(2)) * (exp(t1a) - exp(t1b)) / (sqrt(pi) * sqt)
    t2a <- 2 * stats::pnorm((-k + a) / sqt + sqt * l) - 1
    t2b <- 2 * stats::pnorm((k + a) / sqt - sqt * l) - 1
    t2 <- exp(log(0.5) + log(l)) * (t2a + t2b)
    pdf <- exp(log(max(1e-300, t1 + t2)) - log(2) - log(a))
  }
  if (!is.finite(pdf) || pdf <= 0) 0 else pdf
}

dist_rdm_pdf <- function(x, v, B, A, s) {
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(0)
    .rdm_digt_ref(xi, B / s + 0.5 * A / s, v / s, 0.5 * A / s)
  })
}

dist_rdm_cdf <- function(x, v, B, A, s) {
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  .dist_vapply_numeric(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (!is.finite(xi)) return(if (xi < 0) 0 else 1)
    .rdm_pigt_ref(xi, B / s + 0.5 * A / s, v / s, 0.5 * A / s)
  })
}

.rdm_rng_single <- function(k, l, tiny = 1e-6) {
  if (!is.finite(k) || !is.finite(l) || l < 0) return(Inf)
  if (l <= tiny) {
    q <- max(1e-15, 1 - stats::runif(1) / 2)
    z <- stats::qnorm(q)
    if (!is.finite(z) || abs(z) < 1e-15) return(Inf)
    return((k * k) / (z * z))
  }
  mu <- k / l
  lambda <- k * k
  if (!is.finite(mu) || !is.finite(lambda) || mu <= 0 || lambda <= 0) return(Inf)
  y <- stats::rnorm(1)^2
  root <- sqrt(4 * mu * lambda * y + mu * mu * y * y)
  x <- mu + (mu * mu * y) / (2 * lambda) - (mu * root) / (2 * lambda)
  if (stats::runif(1) > mu / (mu + x)) {
    x <- (mu * mu) / x
  }
  if (!is.finite(x) || x <= 0) Inf else x
}

dist_rdm_rng <- function(n, v, B, A, s) {
  n <- .dist_as_count(n)
  if (n == 0L || !is.finite(v) || !is.finite(B) || !is.finite(A) || !is.finite(s) || s <= 0) {
    return(numeric(0))
  }
  v_sc <- v / s
  B_sc <- max(0, B / s)
  A_sc <- max(0, A / s)
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- .rdm_rng_single(B_sc + stats::runif(1, 0, A_sc), v_sc)
  }
  out
}

.dist_exgauss_shifted_pdf <- function(x, par) {
  t0 <- .dist_param_t0(par)
  shifted <- x - t0
  args <- .dist_param_values(par, "exgauss")
  lower_cdf <- do.call(dist_exgauss_cdf, c(list(0), args))
  lower_survival <- 1 - lower_cdf
  out <- rep(0, length(shifted))
  ok <- !is.na(shifted) & shifted > 0 & is.finite(lower_survival) & lower_survival > 0
  if (any(ok)) {
    raw <- do.call(dist_exgauss_pdf, c(list(shifted[ok]), args))
    out[ok] <- raw / lower_survival
  }
  out[!is.finite(out) | out < 0] <- 0
  out
}

.dist_exgauss_shifted_cdf <- function(x, par) {
  t0 <- .dist_param_t0(par)
  shifted <- x - t0
  args <- .dist_param_values(par, "exgauss")
  lower_cdf <- do.call(dist_exgauss_cdf, c(list(0), args))
  lower_survival <- 1 - lower_cdf
  out <- rep(0, length(shifted))
  ok <- !is.na(shifted) & shifted > 0 & is.finite(lower_survival) & lower_survival > 0
  if (any(ok)) {
    raw <- do.call(dist_exgauss_cdf, c(list(shifted[ok]), args))
    out[ok] <- (raw - lower_cdf) / lower_survival
  }
  out[!is.finite(out)] <- 0
  pmin(pmax(out, 0), 1)
}

.dist_exgauss_shifted_rng <- function(n, par) {
  n_int <- .dist_as_count(n)
  if (n_int == 0L) {
    return(numeric(0))
  }
  t0 <- .dist_param_t0(par)
  args <- .dist_param_values(par, "exgauss")
  out <- numeric(n_int)
  filled <- 0L
  attempts <- 0L
  while (filled < n_int) {
    attempts <- attempts + 1L
    if (attempts > 10000L) {
      stop("exgauss truncation sampler could not draw support above onset + t0", call. = FALSE)
    }
    draw <- do.call(dist_exgauss_rng, c(list(n_int - filled), args))
    keep <- draw[is.finite(draw) & draw > 0]
    if (length(keep) == 0L) {
      next
    }
    take <- min(length(keep), n_int - filled)
    out[seq.int(filled + 1L, length.out = take)] <- keep[seq_len(take)]
    filled <- filled + take
  }
  out + t0
}

.dist_make_entry <- function(dist, pdf_fun, cdf_fun, rng_fun) {
  force(dist)
  force(pdf_fun)
  force(cdf_fun)
  force(rng_fun)
  if (identical(dist, "exgauss")) {
    return(list(
      r = function(n, par) .dist_exgauss_shifted_rng(n, par),
      d = function(x, par) .dist_exgauss_shifted_pdf(x, par),
      p = function(x, par) .dist_exgauss_shifted_cdf(x, par)
    ))
  }
  list(
    r = function(n, par) .dist_shifted_rng(rng_fun, n, par, dist),
    d = function(x, par) .dist_shifted_eval(pdf_fun, x, par, dist),
    p = function(x, par) .dist_shifted_eval(cdf_fun, x, par, dist)
  )
}

dist_registry <- local({
  reg <- new.env(parent = emptyenv())
  families <- list(
    lognormal = list(pdf = dist_lognormal_pdf, cdf = dist_lognormal_cdf, rng = dist_lognormal_rng),
    gamma = list(pdf = dist_gamma_pdf, cdf = dist_gamma_cdf, rng = dist_gamma_rng),
    exgauss = list(pdf = dist_exgauss_pdf, cdf = dist_exgauss_cdf, rng = dist_exgauss_rng),
    lba = list(pdf = dist_lba_pdf, cdf = dist_lba_cdf, rng = dist_lba_rng),
    rdm = list(pdf = dist_rdm_pdf, cdf = dist_rdm_cdf, rng = dist_rdm_rng)
  )

  for (dist in names(families)) {
    family_fns <- families[[dist]]
    reg[[dist]] <- .dist_make_entry(
      dist = dist,
      pdf_fun = family_fns$pdf,
      cdf_fun = family_fns$cdf,
      rng_fun = family_fns$rng
    )
  }

  function(name) {
    key <- if (is.null(name) || length(name) == 0L) "" else tolower(as.character(name)[1])
    if (!nzchar(key) || !exists(key, envir = reg, inherits = FALSE)) {
      stop(sprintf("Unknown distribution '%s'", name), call. = FALSE)
    }
    get(key, envir = reg, inherits = FALSE)
  }
})
