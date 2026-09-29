# Independent unit-exponential oracles for positive-mass ties at shared events.
shared_release_model <- function(outcomes, trigger = FALSE, pool = FALSE) {
  spec <- race_spec()
  for (id in c("G", "A", "B", "C", "D", "X")) {
    spec <- add_accumulator(spec, id, "gamma")
  }
  if (trigger) spec <- add_trigger(spec, "q", members = c("B", "C"))
  if (pool) spec <- add_pool(spec, "P", members = c("G", "A", "B"), k = 2L)
  for (id in names(outcomes)) spec <- add_outcome(spec, id, outcomes[[id]])
  finalize_model(spec)
}

shared_release_density <- function(model, responses, times, q = NULL) {
  data <- expand.grid(rt = times, R = responses, stringsAsFactors = FALSE)
  data$trials <- seq_len(nrow(data))
  params <- build_param_matrix(model, c(shape = 1, rate = 1, t0 = 0, if (!is.null(q)) c(q = q)),
                               n_trials = nrow(data))
  exp(log_likelihood(make_context(model), prepare_data(model, data),
                     params, sum = FALSE))
}

testthat::test_that("overlapping first_of routes count their shared release once", {
  times <- c(0.2, 0.5, 1.5)
  s <- exp(-times)
  # H(t) = integral_0^t PDF_G(u) CDF_B(u) du.
  h <- (1 - s)^2 / 2
  f <- s^2 + 3 * s * h
  survivor <- s * (1 + h)
  for (q in c(0, 0.3, 1)) {
    for (reverse in c(FALSE, TRUE)) {
      children <- list(inhibit("G", by = "B"), all_of("G", "C"))
      if (reverse) children <- rev(children)
      model <- shared_release_model(
        list(R = do.call(first_of, children), D = "D"), trigger = TRUE)
      expected <- c((q * s + (1 - q) * f) * s,
                    s * (q * s + (1 - q) * survivor))
      testthat::expect_equal(shared_release_density(model, c("R", "D"), times, q),
                             expected, tolerance = 2e-5)
    }
  }
})

testthat::test_that("all_of merges prerequisites of children with the same release", {
  times <- c(0.15, 0.4, 1, 2)
  s <- exp(-times)
  for (competitor in c(FALSE, TRUE)) {
    f_r <- if (competitor) s^2 * (1 - s^2) else
      s^2 * (1 - s) + s * (1 - s^2) / 2
    expected <- if (competitor) c(f_r, f_r + 1.5 * s * (1 - s)^2) else f_r
    outputs <- list()
    for (repeated in c(FALSE, TRUE)) {
      expr <- all_of(inhibit("G", by = "B"),
                     if (repeated) all_of("G", "A") else "A")
      outcomes <- list(R = expr)
      if (competitor) outcomes$Q <- all_of("G", "C")
      model <- shared_release_model(outcomes)
      testthat::expect_equal(shared_release_density(model, names(outcomes), times),
                             expected, tolerance = 2e-5)
      if (!competitor) {
        probabilities <- response_probabilities(make_context(model),
          build_param_matrix(model, c(shape = 1, rate = 1, t0 = 0), n_trials = 1))
        testthat::expect_equal(unname(probabilities["R"]), 0.5, tolerance = 2e-5)
      }
      outputs[[length(outputs) + 1L]] <- simulate(model,
        build_param_matrix(model, c(shape = 1, rate = 1, t0 = 0), n_trials = 2000),
        seed = 84231L)
    }
    testthat::expect_identical(outputs[[1]], outputs[[2]])
  }
})

testthat::test_that("redundant guarded routes preserve equal readiness", {
  times <- c(0.2, 0.5, 1.5)
  s <- exp(-times)
  f <- s^2 * (1 - s) + s * (1 - s^2) / 2
  base <- all_of("G", "A")
  for (reverse in c(FALSE, TRUE)) {
    children <- list(base, inhibit(base, by = "B"))
    if (reverse) children <- rev(children)
    model <- shared_release_model(list(R = do.call(first_of, children),
                                       Q = all_of("G", "C")))
    testthat::expect_equal(shared_release_density(model, c("R", "Q"), times),
                           c(f, f), tolerance = 2e-5)
  }
})

testthat::test_that("two guarded conjunctions can finish at their shared event", {
  times <- c(0.2, 0.5, 1.5)
  s <- exp(-times)
  response <- all_of(
    inhibit(all_of("G", "A"), by = "C"),
    inhibit(all_of("G", "B"), by = "D"))
  model <- shared_release_model(list(R = response))
  # Integrating over all strict completion orders gives response mass 2/15.
  expected <- (2 * s^2 + 3 * s^3 - 12 * s^4 + 7 * s^5) / 3
  testthat::expect_equal(shared_release_density(model, "R", times),
                         expected, tolerance = 2e-5)
  model <- shared_release_model(list(R = response, Q = all_of("G", "X")))
  # Enumerating the finishing orders of six exponential runners gives these subdensities.
  f_r <- s^3 - 3 * s^5 + 2 * s^6
  f_q <- 28 * s / 15 - 2 * s^2 + 4 * s^4 / 3 - 2 * s^5 + 4 * s^6 / 5
  testthat::expect_equal(shared_release_density(model, c("R", "Q"), times),
                         c(f_r, f_q), tolerance = 2e-5)
})

testthat::test_that("shared completion does not erase differences in readiness", {
  times <- c(0.2, 0.5, 1.5)
  s <- exp(-times)
  h <- (1 - s)^2 / 2
  for (explicit in c(FALSE, TRUE)) {
    release <- if (explicit) all_of("G", "B", "C") else all_of("G", "C")
    model <- shared_release_model(list(
      R = first_of(inhibit("G", by = "B"), release), Q = all_of("G", "D")))
    at_g <- if (explicit) (1 - s^2) - 2 * (1 - s^3) / 3 else
      (1 - s) * (1 - s^2) / 2
    f_r <- s * (s + at_g) + s^2 * h
    total <- s * (s + (1 - s) * (1 - s^2)) + 2 * s^2 * h
    testthat::expect_equal(shared_release_density(model, c("R", "Q"), times),
                           c(f_r, total - f_r), tolerance = 2e-5)
  }
})

testthat::test_that("reused aggregate pools retain their source dependence", {
  times <- c(0.2, 0.5, 1.5)
  s <- exp(-times)
  # P is second of three: PDF_P = 6 s^2 (1-s).
  f <- 6 * s^3 * (1 - s)
  cdf <- 0.5 - 2 * s^3 + 1.5 * s^4
  model <- shared_release_model(
    list(R = all_of("P", inhibit("P", by = "C")), D = "D"), pool = TRUE)
  testthat::expect_equal(shared_release_density(model, c("R", "D"), times),
                         c(f * s, s * (1 - cdf)), tolerance = 2e-5)
})
