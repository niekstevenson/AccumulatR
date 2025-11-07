.dist_native_env <- local({
  env <- new.env(parent = baseenv())
  env$loaded <- FALSE
  env
})

.load_dist_native <- function() {
  env <- .dist_native_env
  if (isTRUE(env$loaded)) return(env)
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp package is required to use the native distribution kernels", call. = FALSE)
  }
  cpp_path <- file.path("src", "distributions.cpp")
  if (!file.exists(cpp_path)) {
    stop(sprintf("Native distribution source not found at '%s'", cpp_path), call. = FALSE)
  }
  boost_dirs <- c(
    Sys.getenv("UUBER_BOOST_INCLUDE", unset = NA_character_),
    "/opt/homebrew/opt/boost/include"
  )
  boost_dirs <- boost_dirs[nzchar(boost_dirs)]
  boost_dir <- NULL
  if (length(boost_dirs) > 0) {
    for (bd in boost_dirs) {
      if (!is.null(bd) && length(bd) > 0 && nzchar(bd) && file.exists(bd)) {
        boost_dir <- bd
        break
      }
    }
  }
  old_cppflags <- NULL
  if (!is.null(boost_dir)) {
    old_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = "")
    new_flags <- paste(sprintf("-I%s", boost_dir), old_cppflags)
    Sys.setenv(PKG_CPPFLAGS = new_flags)
    on.exit({
      Sys.setenv(PKG_CPPFLAGS = old_cppflags)
    }, add = TRUE)
  }
  rebuild_flag <- getOption("uuber.native.rebuild", TRUE)
  Rcpp::sourceCpp(cpp_path, env = env, rebuild = isTRUE(rebuild_flag))
  env$loaded <- TRUE
  env
}

.get_dist_native <- function() {
  env <- .dist_native_env
  if (!isTRUE(env$loaded)) {
    env <- .load_dist_native()
  }
  env
}

.lik_native_fn <- function(name) {
  env <- .get_dist_native()
  fn <- env[[name]]
  if (!is.function(fn)) {
    stop(sprintf("Native likelihood routine '%s' is unavailable", name), call. = FALSE)
  }
  fn
}

.dist_param_scalar <- function(par, name) {
  val <- par[[name]]
  if (is.null(val) || length(val) == 0L) {
    stop(sprintf("Missing parameter '%s' for distribution", name), call. = FALSE)
  }
  as.numeric(val)[1]
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
      native <- .get_dist_native()
      native$dist_lognormal_rng(
        .dist_as_count(n),
        .dist_param_scalar(par, "meanlog"),
        .dist_param_scalar(par, "sdlog")
      )
    },
    d = function(x, par) {
      native <- .get_dist_native()
      native$dist_lognormal_pdf(
        x,
        .dist_param_scalar(par, "meanlog"),
        .dist_param_scalar(par, "sdlog")
      )
    },
    p = function(x, par) {
      native <- .get_dist_native()
      native$dist_lognormal_cdf(
        x,
        .dist_param_scalar(par, "meanlog"),
        .dist_param_scalar(par, "sdlog")
      )
    }
  )

  reg$gamma <- list(
    r = function(n, par) {
      native <- .get_dist_native()
      native$dist_gamma_rng(
        .dist_as_count(n),
        .dist_param_scalar(par, "shape"),
        .dist_param_scalar(par, "rate")
      )
    },
    d = function(x, par) {
      native <- .get_dist_native()
      native$dist_gamma_pdf(
        x,
        .dist_param_scalar(par, "shape"),
        .dist_param_scalar(par, "rate")
      )
    },
    p = function(x, par) {
      native <- .get_dist_native()
      native$dist_gamma_cdf(
        x,
        .dist_param_scalar(par, "shape"),
        .dist_param_scalar(par, "rate")
      )
    }
  )

  reg$exgauss <- list(
    r = function(n, par) {
      native <- .get_dist_native()
      native$dist_exgauss_rng(
        .dist_as_count(n),
        .dist_param_scalar(par, "mu"),
        .dist_param_scalar(par, "sigma"),
        .dist_param_scalar(par, "tau")
      )
    },
    d = function(x, par) {
      native <- .get_dist_native()
      native$dist_exgauss_pdf(
        x,
        .dist_param_scalar(par, "mu"),
        .dist_param_scalar(par, "sigma"),
        .dist_param_scalar(par, "tau")
      )
    },
    p = function(x, par) {
      native <- .get_dist_native()
      native$dist_exgauss_cdf(
        x,
        .dist_param_scalar(par, "mu"),
        .dist_param_scalar(par, "sigma"),
        .dist_param_scalar(par, "tau")
      )
    }
  )

  function(name) {
    if (!exists(name, envir = reg, inherits = FALSE)) {
      stop(sprintf("Unknown distribution '%s'", name))
    }
    get(name, envir = reg, inherits = FALSE)
  }
})
