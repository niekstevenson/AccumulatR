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

.dist_make_entry <- function(dist, pdf_fun, cdf_fun, rng_fun) {
  force(dist)
  force(pdf_fun)
  force(cdf_fun)
  force(rng_fun)
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
