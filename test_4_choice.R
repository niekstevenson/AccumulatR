#!/usr/bin/env Rscript

# Quick consistency checks for example 16 (four-choice race with shared inhibitors)

suppressWarnings({
  source("examples/new_API.R")
  source("R/model_tables.R")
  source("R/generator_new.R")
  source("R/likelihood_new.R")
})

set.seed(2025)

model_spec <- new_api_examples[[16]]
model_tables <- model_to_tables(model_spec)

# Analytic probabilities
analytic_probs <- observed_response_probabilities(model_tables, include_na = TRUE)
required_labels <- c("A", "B", "C", "D")

cat("Analytic probabilities:\n")
print(round(analytic_probs, 6))

# Basic sanity checks
stopifnot(
  abs(sum(analytic_probs) - 1) < 1e-8,
  all(required_labels %in% names(analytic_probs))
)

min_margin <- 0.02
issues <- character(0)
analytic_margin <- analytic_probs[["B"]] - analytic_probs[["C"]]
if (analytic_margin <= min_margin) {
  issues <- c(issues, sprintf(
    "Analytic probabilities should give B at least %.3f higher than C (B=%.4f, C=%.4f)",
    min_margin, analytic_probs[["B"]], analytic_probs[["C"]]
  ))
}

# Fast simulation for empirical comparison
sim_trials <- 8000
sim_data <- simulate_model(model_tables, n_trials = sim_trials)
sim_counts <- table(sim_data$outcome, useNA = "ifany")
sim_probs <- sim_counts / sum(sim_counts)

cat(sprintf("\nSimulation probabilities (n = %d):\n", sim_trials))
print(round(sim_probs, 6))

# Ensure simulated outcomes cover the required labels
stopifnot(all(required_labels %in% names(sim_probs)))

# B should still outrank C in simulation
sim_margin <- sim_probs[["B"]] - sim_probs[["C"]]
if (sim_margin <= min_margin) {
  issues <- c(issues, sprintf(
    "Simulation suggests B should exceed C by more than %.3f (B=%.4f, C=%.4f)",
    min_margin, sim_probs[["B"]], sim_probs[["C"]]
  ))
}

# Analytic vs simulation should be reasonably close
tolerance <- 0.02
diffs <- abs(as.numeric(sim_probs[required_labels]) - as.numeric(analytic_probs[required_labels]))
names(diffs) <- required_labels

cat("\nAbsolute differences (simulation - analytic):\n")
print(round(diffs, 6))

if (!all(diffs < tolerance)) {
  issues <- c(issues, sprintf(
    "Differences exceed tolerance %.3f (max diff %.3f)",
    tolerance, max(diffs)
  ))
}

if (length(issues) == 0) {
  cat("\nAll checks passed.\n")
} else {
  cat("\nIssues detected:\n")
  for (msg in issues) cat(" - ", msg, "\n", sep = "")
}

# Diagnostic: show unconditional vs conditional survival for B_no_b1 relative to blocker b2
diagnostic_t <- 0.30
prep <- .prepare_model_for_likelihood(model_tables)
S_blocker <- .acc_survival(prep$accumulators[["b2"]], diagnostic_t)
S_pool_uncond <- .pool_survival(prep, "B_no_b1", "__default__", diagnostic_t)

cat(sprintf("\nDiagnostic at t = %.2f:\n", diagnostic_t))
cat(sprintf("  S_blocker(b2 > t)          = %.4f\n", S_blocker))
cat(sprintf("  S_pool_uncond(B_no_b1 > t) = %.4f\n", S_pool_uncond))

detail_trials <- 5000
sim_detail <- simulate_model(model_tables, n_trials = detail_trials, keep_detail = TRUE)
detail_info <- attr(sim_detail, "details")
b2_times <- vapply(detail_info, function(d) as.numeric(d$acc_times[["b2"]]), numeric(1))
b3_times <- vapply(detail_info, function(d) as.numeric(d$acc_times[["b3"]]), numeric(1))

cond_mask <- b2_times > diagnostic_t
cond_denom <- mean(cond_mask)
cond_num <- mean(cond_mask & pmax(b2_times, b3_times) > diagnostic_t)
S_pool_cond_hat <- if (cond_denom > 0) cond_num / cond_denom else NA_real_

cat(sprintf("  Empirical P(B_no_b1 > t | b2 > t) = %.4f\n", S_pool_cond_hat))

if (length(issues) > 0) {
  stop("Test failed – see issues above.")
}
