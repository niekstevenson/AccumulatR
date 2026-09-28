testthat::test_that("ranked races remain factored and match their joint density", {
  # Six racers already distinguish the old 192-cell expansion from <= 36 cells,
  # without exhausting memory if that expansion is accidentally reintroduced.
  n <- 6L
  ids <- paste0("R", seq_len(n))
  spec <- race_spec(n_outcomes = 3L)
  for (id in ids) {
    spec <- add_accumulator(spec, id, "lognormal")
    spec <- add_outcome(spec, id, id)
  }
  model <- finalize_model(test_separate_all_parameters(spec))
  means <- log(seq(0.3, 0.9, length.out = n))
  params <- unlist(stats::setNames(lapply(means, function(m) {
    c(m = m, s = 0.4, t0 = 0)
  }), ids))
  data <- data.frame(
    trials = 1:4, R = ids[c(1, 4, 6, 2)], rt = c(0.2, 0.3, 0.25, 0.2),
    R2 = ids[c(4, 6, 2, 1)], rt2 = c(0.3, 0.4, 0.35, 0.3),
    R3 = ids[c(6, 2, 1, 4)], rt3 = c(0.4, 0.5, 0.45, 0.4)
  )
  expected <- vapply(seq_len(nrow(data)), function(i) {
    winners <- match(unlist(data[i, c("R", "R2", "R3")]), ids)
    times <- unlist(data[i, c("rt", "rt2", "rt3")])
    sum(dlnorm(times, means[winners], 0.4, log = TRUE)) +
      sum(plnorm(data$rt3[i], means[-winners], 0.4,
                 lower.tail = FALSE, log.p = TRUE))
  }, numeric(1))
  context <- make_context(model, diagnostics = TRUE)
  metrics <- complexity_metrics(context)$total
  testthat::expect_lte(metrics$symbolic_cells, n * n)
  testthat::expect_lte(metrics$compiled_nodes, 8L * n)
  actual <- log_likelihood(
    context, prepare_data(model, data),
    build_param_matrix(model, params, n_trials = nrow(data)),
    sum = FALSE, min_ll = -1000
  )
  testthat::expect_equal(as.numeric(actual), expected, tolerance = 1e-10)
})
