suppressPackageStartupMessages({
  library(pkgload)
  loaded <- FALSE
  try({
    load_all(".", compile = FALSE, quiet = TRUE)
    loaded <- TRUE
  }, silent = TRUE)
  if (!loaded) {
    load_all(".", compile = TRUE, quiet = TRUE)
  }
})

check_true <- function(ok, msg) {
  if (!isTRUE(ok)) {
    stop(msg, call. = FALSE)
  }
  cat("PASS:", msg, "\n")
}

check_close <- function(x, y, tol, msg) {
  if (!is.finite(x) || !is.finite(y) || abs(x - y) > tol) {
    stop(sprintf("%s (x=%.12f, y=%.12f, tol=%.3g)", msg, x, y, tol), call. = FALSE)
  }
  cat("PASS:", msg, "\n")
}

check_component_inactive_source <- function() {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("B", "b") |>
    add_component("with_source", members = c("a", "b"), weight = 0.5) |>
    add_component("without_source", members = c("b"), weight = 0.5)
  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.22), a.sdlog = 0.10, a.q = 0.00, a.t0 = 0.00,
    b.meanlog = log(0.20), b.sdlog = 0.10, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  min_ll <- -1e6
  observed_df <- data.frame(
    trial = 1L,
    component = "without_source",
    R = "B",
    rt = 0.70,
    stringsAsFactors = FALSE
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, observed_df),
    params_df,
    min_ll = min_ll
  ))
  check_true(identical(ll, min_ll), "inactive source component forces non-start (min_ll)")
}

check_ranked_observed_source_shortcut <- function() {
  spec <- race_spec(n_outcomes = 2L) |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal", onset = after("a")) |>
    add_outcome("A", "a") |>
    add_outcome("B", "b")
  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.25), a.sdlog = 0.10, a.q = 0.00, a.t0 = 0.00,
    b.meanlog = log(0.20), b.sdlog = 0.10, b.q = 0.00, b.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  obs <- data.frame(
    trial = 1L,
    R = "A",
    rt = 0.40,
    R2 = "B",
    rt2 = 0.65,
    stringsAsFactors = FALSE
  )
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, obs),
    params_df,
    min_ll = -1e9
  ))
  manual <- dlnorm(0.40, meanlog = log(0.25), sdlog = 0.10) *
    dlnorm(0.65 - 0.40, meanlog = log(0.20), sdlog = 0.10)
  check_close(ll, log(manual), tol = 1e-6, msg = "ranked A->B uses exact shifted density shortcut")
}

check_pool_source_finite <- function() {
  spec <- race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_pool("p", c("a", "b"), k = 1L) |>
    add_accumulator("c", "lognormal", onset = after("p", lag = 0.03)) |>
    add_outcome("C", "c")
  structure <- finalize_model(spec)
  params <- c(
    a.meanlog = log(0.18), a.sdlog = 0.12, a.q = 0.20, a.t0 = 0.00,
    b.meanlog = log(0.22), b.sdlog = 0.12, b.q = 0.20, b.t0 = 0.00,
    c.meanlog = log(0.15), c.sdlog = 0.10, c.q = 0.00, c.t0 = 0.00
  )
  params_df <- build_param_matrix(spec, params, n_trials = 1L)
  obs <- data.frame(trial = 1L, R = "C", rt = 0.70, stringsAsFactors = FALSE)
  ll <- as.numeric(log_likelihood(
    build_likelihood_context(structure, obs),
    params_df,
    min_ll = -1e9
  ))
  check_true(is.finite(ll), "pool-source chained onset returns finite log-likelihood")
}

check_component_inactive_source()
check_ranked_observed_source_shortcut()
check_pool_source_finite()

cat("\nAll manual chained-onset checks passed.\n")
