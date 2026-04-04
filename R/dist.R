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

dist_registry <- local({
  reg <- new.env(parent = emptyenv())

  reg$lognormal <- list(
    r = function(n, par) {
      m <- .dist_param_scalar(par, "m")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      dist_lognormal_rng(.dist_as_count(n), m, s) + t0
    },
    d = function(x, par) {
      m <- .dist_param_scalar(par, "m")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_lognormal_pdf(shifted, m, s)
      .dist_zero_before_t0(out, shifted)
    },
    p = function(x, par) {
      m <- .dist_param_scalar(par, "m")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_lognormal_cdf(shifted, m, s)
      .dist_zero_before_t0(out, shifted)
    }
  )

  reg$gamma <- list(
    r = function(n, par) {
      shape <- .dist_param_scalar(par, "shape")
      rate <- .dist_param_scalar(par, "rate")
      t0 <- .dist_param_t0(par)
      dist_gamma_rng(.dist_as_count(n), shape, rate) + t0
    },
    d = function(x, par) {
      shape <- .dist_param_scalar(par, "shape")
      rate <- .dist_param_scalar(par, "rate")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_gamma_pdf(shifted, shape, rate)
      .dist_zero_before_t0(out, shifted)
    },
    p = function(x, par) {
      shape <- .dist_param_scalar(par, "shape")
      rate <- .dist_param_scalar(par, "rate")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_gamma_cdf(shifted, shape, rate)
      .dist_zero_before_t0(out, shifted)
    }
  )

  reg$exgauss <- list(
    r = function(n, par) {
      mu <- .dist_param_scalar(par, "mu")
      sigma <- .dist_param_scalar(par, "sigma")
      tau <- .dist_param_scalar(par, "tau")
      t0 <- .dist_param_t0(par)
      dist_exgauss_rng(.dist_as_count(n), mu, sigma, tau) + t0
    },
    d = function(x, par) {
      mu <- .dist_param_scalar(par, "mu")
      sigma <- .dist_param_scalar(par, "sigma")
      tau <- .dist_param_scalar(par, "tau")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_exgauss_pdf(shifted, mu, sigma, tau)
      .dist_zero_before_t0(out, shifted)
    },
    p = function(x, par) {
      mu <- .dist_param_scalar(par, "mu")
      sigma <- .dist_param_scalar(par, "sigma")
      tau <- .dist_param_scalar(par, "tau")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_exgauss_cdf(shifted, mu, sigma, tau)
      .dist_zero_before_t0(out, shifted)
    }
  )

  reg$lba <- list(
    r = function(n, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      sv <- .dist_param_scalar(par, "sv")
      t0 <- .dist_param_t0(par)
      dist_lba_rng(.dist_as_count(n), v, B, A, sv) + t0
    },
    d = function(x, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      sv <- .dist_param_scalar(par, "sv")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_lba_pdf(shifted, v, B, A, sv)
      .dist_zero_before_t0(out, shifted)
    },
    p = function(x, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      sv <- .dist_param_scalar(par, "sv")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_lba_cdf(shifted, v, B, A, sv)
      .dist_zero_before_t0(out, shifted)
    }
  )

  reg$rdm <- list(
    r = function(n, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      dist_rdm_rng(.dist_as_count(n), v, B, A, s) + t0
    },
    d = function(x, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_rdm_pdf(shifted, v, B, A, s)
      .dist_zero_before_t0(out, shifted)
    },
    p = function(x, par) {
      v <- .dist_param_scalar(par, "v")
      B <- .dist_param_scalar(par, "B")
      A <- .dist_param_scalar(par, "A")
      s <- .dist_param_scalar(par, "s")
      t0 <- .dist_param_t0(par)
      shifted <- x - t0
      out <- dist_rdm_cdf(shifted, v, B, A, s)
      .dist_zero_before_t0(out, shifted)
    }
  )

  function(name) {
    key <- if (is.null(name) || length(name) == 0L) "" else tolower(as.character(name)[1])
    if (!nzchar(key) || !exists(key, envir = reg, inherits = FALSE)) {
      stop(sprintf("Unknown distribution '%s'", name))
    }
    get(key, envir = reg, inherits = FALSE)
  }
})
