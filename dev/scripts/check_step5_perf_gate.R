#!/usr/bin/env Rscript

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

diff_file <- file.path("dev", "scripts", "scratch_outputs", "benchmark_centralized_diff.csv")
if (!file.exists(diff_file)) {
  stop("Missing benchmark diff file: ", diff_file)
}

diff_df <- utils::read.csv(diff_file, stringsAsFactors = FALSE)
if (nrow(diff_df) == 0L) {
  stop("Benchmark diff is empty: ", diff_file)
}

get_ratio <- function(label, metric = "median_per_eval_sec", suite = "nested") {
  rows <- diff_df[
    diff_df$section == "speed" &
      diff_df$suite == suite &
      diff_df$label == label &
      diff_df$metric == metric &
      diff_df$mode == "ir",
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) {
    stop("Missing or duplicated ratio row for label=", label, ", metric=", metric)
  }
  as.numeric(rows$ratio_current_over_baseline[[1]])
}

targets <- c("shared_gate_pair", "shared_gate_nway", "shared_gate_nway_shared_triggers")
target_ratios <- vapply(targets, get_ratio, numeric(1))

guard_ratio <- get_ratio("depth3_guard_competitor")

if (any(!is.finite(target_ratios))) {
  stop("Non-finite shared-gate performance ratios detected")
}
if (!is.finite(guard_ratio)) {
  stop("Non-finite depth3_guard_competitor ratio detected")
}

if (any(target_ratios > 1.0)) {
  bad <- paste(names(target_ratios)[target_ratios > 1.0], collapse = ", ")
  stop("Shared-gate regression detected for: ", bad)
}

target_median_ratio <- stats::median(target_ratios)
if (!(target_median_ratio <= 0.98)) {
  stop(
    "Target improvement not achieved; median shared-gate ratio=",
    format(target_median_ratio, digits = 6),
    " (need <= 0.98)"
  )
}

if (guard_ratio > 1.05) {
  stop(
    "depth3_guard_competitor regression too large: ratio=",
    format(guard_ratio, digits = 6),
    " (limit 1.05)"
  )
}

cat(
  "Step 5 performance gate passed: median shared-gate ratio=",
  format(target_median_ratio, digits = 6),
  ", depth3_guard_competitor ratio=",
  format(guard_ratio, digits = 6),
  "\n",
  sep = ""
)
