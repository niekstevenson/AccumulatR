dist_registry <- local({
  # Registry of base distributions supported: lognormal, gamma, exgauss
  # Each entry returns list(r, d, p) functions operating on scalar x and param list.
  exgauss_pdf <- function(x, par) {
    mu <- par$mu; sigma <- par$sigma; tau <- par$tau
    if (!is.finite(x)) return(0.0)
    # Handle sigma or tau edge cases
    if (sigma <= 0 || tau <= 0) return(NA_real_)
    z <- (x - mu) / sigma
    w <- sigma / tau
    # f(x) = (1/tau) * exp( (sigma^2)/(2 tau^2) - (x-mu)/tau ) * Phi( z - w )
    return((1.0/tau) * exp( (sigma*sigma)/(2.0*tau*tau) - (x - mu)/tau ) * pnorm(z - w))
  }
  exgauss_cdf <- function(x, par) {
    mu <- par$mu; sigma <- par$sigma; tau <- par$tau
    if (!is.finite(x)) return(NA_real_)
    if (sigma <= 0 || tau <= 0) return(NA_real_)
    z <- (x - mu) / sigma
    w <- sigma / tau
    # F(x) = Phi(z) - exp( (sigma^2)/(2 tau^2) - (x-mu)/tau ) * Phi( z - w )
    return(pnorm(z) - exp( (sigma*sigma)/(2.0*tau*tau) - (x - mu)/tau ) * pnorm(z - w))
  }
  exgauss_rng <- function(n, par) {
    mu <- par$mu; sigma <- par$sigma; tau <- par$tau
    if (sigma <= 0 || tau <= 0) stop("Invalid exgauss parameters")
    rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1.0/tau)
  }

  reg <- new.env(parent = emptyenv())
  reg$lognormal <- list(
    r = function(n, par) rlnorm(n, meanlog = par$meanlog, sdlog = par$sdlog),
    d = function(x, par) dlnorm(x, meanlog = par$meanlog, sdlog = par$sdlog),
    p = function(x, par) plnorm(x, meanlog = par$meanlog, sdlog = par$sdlog)
  )
  reg$gamma <- list(
    r = function(n, par) rgamma(n, shape = par$shape, rate = par$rate),
    d = function(x, par) dgamma(x, shape = par$shape, rate = par$rate),
    p = function(x, par) pgamma(x, shape = par$shape, rate = par$rate)
  )
  reg$exgauss <- list(
    r = function(n, par) exgauss_rng(n, par),
    d = function(x, par) exgauss_pdf(x, par),
    p = function(x, par) exgauss_cdf(x, par)
  )

  # accessor
  function(name) {
    if (!exists(name, envir = reg, inherits = FALSE)) {
      stop(sprintf("Unknown distribution '%s'", name))
    }
    get(name, envir = reg, inherits = FALSE)
  }
})

