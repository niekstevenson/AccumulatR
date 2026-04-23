#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(AccumulatR)
  }
})

source(file.path("tests", "testthat", "helper-loglik-golden.R"))

fixture <- build_loglik_golden_fixture()
out_dir <- file.path("tests", "testthat", "fixtures")
out_path <- file.path(out_dir, "loglik_golden_v1.rds")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(fixture, out_path)
cat("Wrote", out_path, "\n")
for (case in fixture$cases) {
  cat(sprintf(
    "%s: trials=%d total_loglik=%.10f\n",
    case$label,
    length(case$trial_loglik),
    case$total_loglik
  ))
}
