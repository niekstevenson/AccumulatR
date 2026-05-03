#!/usr/bin/env Rscript

set.seed(20260218L)
if (requireNamespace("devtools", quietly = TRUE)) {
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(AccumulatR))
}

source(file.path("tests", "testthat", "helper-loglik-golden.R"))

fixture_path <- file.path("tests", "testthat", "fixtures", "loglik_golden_v1.rds")
if (!file.exists(fixture_path)) {
  stop("Golden fixture missing: ", fixture_path)
}
golden <- readRDS(fixture_path)
golden_map <- setNames(golden$cases, vapply(golden$cases, `[[`, character(1), "label"))
registry <- loglik_golden_case_registry()

if (!setequal(names(registry), names(golden_map))) {
  stop("Golden fixture case set does not match registry")
}

max_diff <- 0.0
max_tolerance <- 0.0
failed <- FALSE
for (label in names(registry)) {
  gl <- golden_map[[label]]
  values <- evaluate_loglik_golden_case(
    registry[[label]]$build(),
    gl$data_df,
    gl$params
  )
  diff_total <- abs(values$total_loglik - gl$total_loglik)
  diff_trial <- max(abs(values$trial_loglik - gl$trial_loglik))
  case_max <- max(diff_total, diff_trial)
  tolerance <- gl$tolerance %||% golden$tolerance %||% 1e-4
  max_diff <- max(max_diff, case_max)
  max_tolerance <- max(max_tolerance, tolerance)
  failed <- failed || !is.finite(case_max) || case_max > tolerance
  cat(sprintf(
    "%s: total_diff=%.6g trial_max_diff=%.6g tolerance=%.6g\n",
    label, diff_total, diff_trial, tolerance
  ))
}

cat(sprintf("max_abs_diff=%.6g max_tolerance=%.6g\n", max_diff, max_tolerance))
if (failed) {
  stop("Golden check failed")
}
cat("Golden check passed\n")
