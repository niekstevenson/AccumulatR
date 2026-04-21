lba_denom_ref <- function(v, sv) {
  denom <- pnorm(v / sv)
  if (!is.finite(denom) || denom < 1e-10) denom <- 1e-10
  denom
}

lba_pdf_ref <- function(x, v, B, A, sv) {
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) || !is.finite(A) ||
      !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom_ref(v, sv)
  if (A > 1e-10) {
    zs <- x * sv
    cmz <- B - x * v
    cz <- cmz / zs
    cz_max <- (cmz - A) / zs
    pdf <- (v * (pnorm(cz) - pnorm(cz_max)) +
      sv * (dnorm(cz_max) - dnorm(cz))) / (A * denom)
  } else {
    pdf <- dnorm(B / x, mean = v, sd = sv) * B / (x * x * denom)
  }
  if (!is.finite(pdf) || pdf <= 0) 0 else pdf
}

lba_cdf_ref <- function(x, v, B, A, sv) {
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) || !is.finite(A) ||
      !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom_ref(v, sv)
  if (A > 1e-10) {
    zs <- x * sv
    cmz <- B - x * v
    xx <- cmz - A
    cz <- cmz / zs
    cz_max <- xx / zs
    cdf <- (1 +
      (zs * (dnorm(cz_max) - dnorm(cz)) +
         xx * pnorm(cz_max) -
         cmz * pnorm(cz)) / A) / denom
  } else {
    cdf <- pnorm(B / x, mean = v, sd = sv, lower.tail = FALSE) / denom
  }
  min(max(cdf, 0), 1)
}

rdm_pigt0_ref <- function(x, k, l) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l)) return(0)
  if (abs(l) < 1e-12) {
    z <- k / sqrt(x)
    return(min(max(2 * (1 - pnorm(z)), 0), 1))
  }
  mu <- k / l
  lambda <- k * k
  p1 <- 1 - pnorm(sqrt(lambda / x) * (1 + x / mu))
  p2 <- 1 - pnorm(sqrt(lambda / x) * (1 - x / mu))
  part <- exp(exp(log(2 * lambda) - log(mu)) + log(max(1e-300, p1)))
  min(max(part + p2, 0), 1)
}

rdm_digt0_ref <- function(x, k, l) {
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

rdm_pigt_ref <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) || !is.finite(a)) return(0)
  if (a < threshold) return(rdm_pigt0_ref(x, k, l))
  sqt <- sqrt(x)
  lgt <- log(x)
  if (l < threshold) {
    t5a <- 2 * pnorm((k + a) / sqt) - 1
    t5b <- 2 * pnorm((-k - a) / sqt) - 1
    t6a <- -0.5 * ((k + a) * (k + a) / x - log(2) - log(pi) + lgt) - log(a)
    t6b <- -0.5 * ((k - a) * (k - a) / x - log(2) - log(pi) + lgt) - log(a)
    cdf <- 1 + exp(t6a) - exp(t6b) + ((-k + a) * t5a - (k - a) * t5b) / (2 * a)
  } else {
    t1a <- exp(-0.5 * (k - a - x * l)^2 / x)
    t1b <- exp(-0.5 * (a + k - x * l)^2 / x)
    t1 <- exp(0.5 * (lgt - log(2) - log(pi))) * (t1a - t1b)
    t2a <- exp(2 * l * (k - a) + pnorm(-(k - a + x * l) / sqt, log.p = TRUE))
    t2b <- exp(2 * l * (k + a) + pnorm(-(k + a + x * l) / sqt, log.p = TRUE))
    t2 <- a + (t2b - t2a) / (2 * l)
    t4a <- 2 * pnorm((k + a) / sqt - sqt * l) - 1
    t4b <- 2 * pnorm((k - a) / sqt - sqt * l) - 1
    t4 <- 0.5 * (x * l - a - k + 0.5 / l) * t4a +
      0.5 * (k - a - x * l - 0.5 / l) * t4b
    cdf <- 0.5 * (t4 + t2 + t1) / a
  }
  if (!is.finite(cdf) || cdf < 0) return(0)
  min(max(cdf, 0), 1)
}

rdm_digt_ref <- function(x, k, l, a, threshold = 1e-10) {
  if (!is.finite(x) || x <= 0 || !is.finite(k) || !is.finite(l) || !is.finite(a)) return(0)
  if (a < threshold) {
    pdf <- rdm_digt0_ref(x, k, l)
    return(if (!is.finite(pdf) || pdf <= 0) 0 else pdf)
  }
  if (l < threshold) {
    term <- exp(-(k - a)^2 / (2 * x)) - exp(-(k + a)^2 / (2 * x))
    pdf <- exp(-0.5 * (log(2) + log(pi) + log(x)) +
      log(max(1e-300, term)) - log(2) - log(a))
  } else {
    sqt <- sqrt(x)
    t1a <- -(a - k + x * l)^2 / (2 * x)
    t1b <- -(a + k - x * l)^2 / (2 * x)
    t1 <- (1 / sqrt(2)) * (exp(t1a) - exp(t1b)) / (sqrt(pi) * sqt)
    t2a <- 2 * pnorm((-k + a) / sqt + sqt * l) - 1
    t2b <- 2 * pnorm((k + a) / sqt - sqt * l) - 1
    t2 <- exp(log(0.5) + log(l)) * (t2a + t2b)
    pdf <- exp(log(max(1e-300, t1 + t2)) - log(2) - log(a))
  }
  if (!is.finite(pdf) || pdf <= 0) 0 else pdf
}

rdm_pdf_ref <- function(x, v, B, A, s) {
  if (!is.finite(s) || s <= 0) return(0)
  rdm_digt_ref(x, B / s + 0.5 * A / s, v / s, 0.5 * A / s)
}

rdm_cdf_ref <- function(x, v, B, A, s) {
  if (!is.finite(s) || s <= 0) return(0)
  rdm_pigt_ref(x, B / s + 0.5 * A / s, v / s, 0.5 * A / s)
}

exgauss_raw_pdf_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || !is.finite(mu) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(tau) || tau <= 0) {
    return(0)
  }
  inv_tau <- 1 / tau
  sigma_sq <- sigma * sigma
  tau_sq <- tau * tau
  z <- (x - mu) / sigma
  exponent <- sigma_sq / (2 * tau_sq) - (x - mu) * inv_tau
  tail <- pnorm(z - sigma * inv_tau)
  out <- inv_tau * exp(exponent) * tail
  if (!is.finite(out) || out <= 0) 0 else out
}

exgauss_raw_cdf_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || !is.finite(mu) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(tau) || tau <= 0) {
    return(0)
  }
  inv_tau <- 1 / tau
  sigma_sq <- sigma * sigma
  tau_sq <- tau * tau
  z <- (x - mu) / sigma
  exponent <- sigma_sq / (2 * tau_sq) - (x - mu) * inv_tau
  tail <- pnorm(z - sigma * inv_tau)
  base <- pnorm(z)
  min(max(base - exp(exponent) * tail, 0), 1)
}

exgauss_trunc_pdf_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || x <= 0) return(0)
  lower_cdf <- exgauss_raw_cdf_ref(0, mu, sigma, tau)
  lower_survival <- 1 - lower_cdf
  if (!is.finite(lower_survival) || lower_survival <= 0) return(0)
  exgauss_raw_pdf_ref(x, mu, sigma, tau) / lower_survival
}

exgauss_trunc_cdf_ref <- function(x, mu, sigma, tau) {
  if (!is.finite(x) || x <= 0) return(0)
  lower_cdf <- exgauss_raw_cdf_ref(0, mu, sigma, tau)
  lower_survival <- 1 - lower_cdf
  if (!is.finite(lower_survival) || lower_survival <= 0) return(0)
  min(max((exgauss_raw_cdf_ref(x, mu, sigma, tau) - lower_cdf) / lower_survival, 0), 1)
}

run_public_loglik <- function(spec, trial_df, params, min_ll = log(1e-10)) {
  structure <- finalize_model(spec)
  context <- make_context(structure)
  prepared <- prepare_data(structure, trial_df)
  params_mat <- build_param_matrix(spec, params, trial_df = prepared)
  as.numeric(log_likelihood(context, prepared, params_mat, min_ll = min_ll))
}

testthat::test_that("direct kernel matches simple two-leaf top-1 formula", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "gamma") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  trial_df <- data.frame(
    trial = c(1L, 2L),
    R = c("A", "B"),
    rt = c(0.30, 0.55),
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.32), a.s = 0.18, a.q = 0.10, a.t0 = 0.02,
    b.shape = 3.0, b.rate = 8.0, b.q = 0.20, b.t0 = 0.05
  )
  out <- run_public_loglik(spec, trial_df, params)

  pa1 <- (1 - params[["a.q"]]) * dlnorm(trial_df$rt[[1]] - params[["a.t0"]], params[["a.m"]], params[["a.s"]])
  sb1 <- 1 - (1 - params[["b.q"]]) * pgamma(trial_df$rt[[1]] - params[["b.t0"]], shape = params[["b.shape"]], rate = params[["b.rate"]])
  pb2 <- (1 - params[["b.q"]]) * dgamma(trial_df$rt[[2]] - params[["b.t0"]], shape = params[["b.shape"]], rate = params[["b.rate"]])
  sa2 <- 1 - (1 - params[["a.q"]]) * plnorm(trial_df$rt[[2]] - params[["a.t0"]], params[["a.m"]], params[["a.s"]])
  expected <- log(c(pa1 * sb1, pb2 * sa2))

  testthat::expect_equal(out, sum(expected), tolerance = 1e-10)
})

testthat::test_that("direct kernel matches pooled top-1 formula", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("c", "lognormal") |>
    add_pool("ab", members = c("a", "b"), k = 1L) |>
    add_outcome("AB", "ab") |>
    add_outcome("C", "c")

  trial_df <- data.frame(
    trial = 1L,
    R = "AB",
    rt = 0.42,
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.35), a.s = 0.16, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.45), b.s = 0.18, b.q = 0.00, b.t0 = 0.00,
    c.m = log(0.55), c.s = 0.20, c.q = 0.00, c.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  t <- trial_df$rt[[1]]
  fa <- dlnorm(t, params[["a.m"]], params[["a.s"]])
  fb <- dlnorm(t, params[["b.m"]], params[["b.s"]])
  fc <- dlnorm(t, params[["c.m"]], params[["c.s"]])
  sa <- 1 - plnorm(t, params[["a.m"]], params[["a.s"]])
  sb <- 1 - plnorm(t, params[["b.m"]], params[["b.s"]])
  sc <- 1 - plnorm(t, params[["c.m"]], params[["c.s"]])
  expected <- log((fa * sb + fb * sa) * sc)

  testthat::expect_equal(out, expected, tolerance = 1e-10)
})

testthat::test_that("direct likelihood handles contiguous trials across components", {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b") |>
    add_component("one", members = c("a")) |>
    add_component("two", members = c("a", "b"))

  trial_df <- data.frame(
    trial = 1:3,
    component = c("one", "one", "two"),
    R = c("A", "A", "B"),
    rt = c(0.30, 0.33, 0.40),
    stringsAsFactors = FALSE
  )
  params <- c(
    a.m = log(0.31), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
    b.m = log(0.44), b.s = 0.18, b.q = 0.00, b.t0 = 0.00
  )
  out <- run_public_loglik(spec, trial_df, params)

  expected <- c(
    log(dlnorm(0.30, params[["a.m"]], params[["a.s"]])),
    log(dlnorm(0.33, params[["a.m"]], params[["a.s"]])),
    log(dlnorm(0.40, params[["b.m"]], params[["b.s"]]) *
          (1 - plnorm(0.40, params[["a.m"]], params[["a.s"]])))
  )

  testthat::expect_equal(out, sum(expected), tolerance = 1e-10)
})

testthat::test_that("direct kernel matches LBA top-1 formula", {
  spec <- race_spec() |>
    add_accumulator("a", "LBA") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  trial_df <- data.frame(
    trial = c(1L, 2L),
    R = c("A", "B"),
    rt = c(0.55, 0.68),
    stringsAsFactors = FALSE
  )
  params <- c(
    a.v = 2.0, a.B = 1.2, a.A = 0.4, a.sv = 0.6, a.q = 0.0, a.t0 = 0.05,
    b.m = log(0.60), b.s = 0.18, b.q = 0.0, b.t0 = 0.02
  )
  out <- run_public_loglik(spec, trial_df, params)

  xa <- trial_df$rt - params[["a.t0"]]
  xb <- trial_df$rt - params[["b.t0"]]
  expected <- log(c(
    lba_pdf_ref(xa[[1]], params[["a.v"]], params[["a.B"]], params[["a.A"]], params[["a.sv"]]) *
      (1 - plnorm(xb[[1]], params[["b.m"]], params[["b.s"]])),
    dlnorm(xb[[2]], params[["b.m"]], params[["b.s"]]) *
      (1 - lba_cdf_ref(xa[[2]], params[["a.v"]], params[["a.B"]], params[["a.A"]], params[["a.sv"]]))
  ))

  testthat::expect_equal(out, sum(expected), tolerance = 1e-10)
})

testthat::test_that("direct kernel matches RDM top-1 formula", {
  spec <- race_spec() |>
    add_accumulator("a", "RDM") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  trial_df <- data.frame(
    trial = c(1L, 2L),
    R = c("A", "B"),
    rt = c(0.52, 0.72),
    stringsAsFactors = FALSE
  )
  params <- c(
    a.v = 1.5, a.B = 1.0, a.A = 0.3, a.s = 1.1, a.q = 0.0, a.t0 = 0.04,
    b.m = log(0.62), b.s = 0.17, b.q = 0.0, b.t0 = 0.01
  )
  out <- run_public_loglik(spec, trial_df, params)

  xa <- trial_df$rt - params[["a.t0"]]
  xb <- trial_df$rt - params[["b.t0"]]
  expected <- log(c(
    rdm_pdf_ref(xa[[1]], params[["a.v"]], params[["a.B"]], params[["a.A"]], params[["a.s"]]) *
      (1 - plnorm(xb[[1]], params[["b.m"]], params[["b.s"]])),
    dlnorm(xb[[2]], params[["b.m"]], params[["b.s"]]) *
      (1 - rdm_cdf_ref(xa[[2]], params[["a.v"]], params[["a.B"]], params[["a.A"]], params[["a.s"]]))
  ))

  testthat::expect_equal(out, sum(expected), tolerance = 1e-10)
})

testthat::test_that("direct kernel truncates exgauss at the elapsed-time lower bound", {
  spec <- race_spec() |>
    add_accumulator("a", "exgauss") |>
    add_accumulator("b", "lognormal") |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")

  trial_df <- data.frame(
    trial = c(1L, 2L),
    R = c("A", "B"),
    rt = c(0.30, 0.58),
    stringsAsFactors = FALSE
  )
  params <- c(
    a.mu = -0.08, a.sigma = 0.10, a.tau = 0.22, a.q = 0.00, a.t0 = 0.07,
    b.m = log(0.52), b.s = 0.16, b.q = 0.00, b.t0 = 0.03
  )
  out <- run_public_loglik(spec, trial_df, params)

  xa <- trial_df$rt - params[["a.t0"]]
  xb <- trial_df$rt - params[["b.t0"]]
  expected <- log(c(
    exgauss_trunc_pdf_ref(xa[[1]], params[["a.mu"]], params[["a.sigma"]], params[["a.tau"]]) *
      (1 - plnorm(xb[[1]], params[["b.m"]], params[["b.s"]])),
    dlnorm(xb[[2]], params[["b.m"]], params[["b.s"]]) *
      (1 - exgauss_trunc_cdf_ref(xa[[2]], params[["a.mu"]], params[["a.sigma"]], params[["a.tau"]]))
  ))

  testthat::expect_equal(out, sum(expected), tolerance = 1e-8)
})

testthat::test_that("identity likelihood supports shared-trigger no-response trials", {
  spec <- race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_outcome("R1", "go1") |>
    add_outcome("R2", "go2") |>
    add_trigger("shared_trigger",
      members = c("go1", "go2"),
      q = 0.10,
      draw = "shared"
    )

  trial_df <- data.frame(
    trial = 1L,
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  params <- c(
    go1.m = log(0.30), go1.s = 0.18, go1.t0 = 0.00,
    go2.m = log(0.32), go2.s = 0.18, go2.t0 = 0.00,
    shared_trigger = 0.10
  )

  out <- run_public_loglik(spec, trial_df, params)

  testthat::expect_equal(out, log(0.10), tolerance = 1e-8)
})

testthat::test_that("identity likelihood supports independent-trigger no-response trials", {
  spec <- race_spec() |>
    add_accumulator("go_left", "lognormal") |>
    add_accumulator("go_right", "lognormal") |>
    add_outcome("Left", "go_left") |>
    add_outcome("Right", "go_right") |>
    add_trigger("q_shared",
      members = c("go_left", "go_right"),
      q = 0.10,
      draw = "independent"
    )

  trial_df <- data.frame(
    trial = 1L,
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  params <- c(
    go_left.m = log(0.30), go_left.s = 0.18, go_left.t0 = 0.00,
    go_right.m = log(0.32), go_right.s = 0.18, go_right.t0 = 0.00,
    q_shared = 0.10
  )

  out <- run_public_loglik(spec, trial_df, params)

  testthat::expect_equal(out, log(0.01), tolerance = 1e-8)
})
