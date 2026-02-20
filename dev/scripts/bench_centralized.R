#!/usr/bin/env Rscript

set.seed(as.integer(Sys.getenv("ACCUMULATR_BENCH_SEED", "20260220")))

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

script_file <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) == 0L) {
    stop("This benchmark script must be run with Rscript --file")
  }
  sub("^--file=", "", file_arg[[1]])
}

repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, "benchmark_centralized.csv")
nested_file <- file.path(out_dir, "bench_nested_integral_baseline.csv")
examples_file <- file.path(out_dir, "bench_examples_models_likelihood.csv")

capture_file_bytes <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "raw", n = file.info(path)$size)
}

restore_file_bytes <- function(path, bytes) {
  if (is.null(bytes)) {
    if (file.exists(path)) {
      file.remove(path)
    }
    return(invisible(NULL))
  }
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(bytes, con)
}

run_script <- function(path) {
  out <- system2("Rscript", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status") %||% 0L
  list(status = as.integer(status), output = out)
}

parse_golden_output <- function(lines) {
  case_rows <- list()
  max_abs_diff <- NA_real_
  tolerance <- NA_real_
  passed <- FALSE

  for (ln in lines) {
    case_m <- regexec(
      "^([^:]+): total_diff=([0-9.eE+-]+) trial_max_diff=([0-9.eE+-]+)$",
      ln
    )
    case_hit <- regmatches(ln, case_m)[[1]]
    if (length(case_hit) == 4L) {
      case_rows[[length(case_rows) + 1L]] <- data.frame(
        section = "correctness",
        suite = "golden",
        label = case_hit[[2]],
        metric = c("total_diff", "trial_max_diff"),
        value = c(as.numeric(case_hit[[3]]), as.numeric(case_hit[[4]])),
        mode = "",
        stringsAsFactors = FALSE
      )
      next
    }

    overall_m <- regexec(
      "^max_abs_diff=([0-9.eE+-]+) tolerance=([0-9.eE+-]+)$",
      ln
    )
    overall_hit <- regmatches(ln, overall_m)[[1]]
    if (length(overall_hit) == 3L) {
      max_abs_diff <- as.numeric(overall_hit[[2]])
      tolerance <- as.numeric(overall_hit[[3]])
      next
    }

    if (identical(trimws(ln), "Golden check passed")) {
      passed <- TRUE
    }
  }

  rows <- do.call(rbind, case_rows)
  overall_rows <- data.frame(
    section = "correctness",
    suite = "golden",
    label = "overall",
    metric = c("max_abs_diff", "tolerance", "passed"),
    value = c(max_abs_diff, tolerance, if (passed) 1.0 else 0.0),
    mode = "",
    stringsAsFactors = FALSE
  )

  if (is.null(rows)) {
    rows <- overall_rows
  } else {
    rows <- rbind(rows, overall_rows)
  }
  rows
}

flatten_speed <- function(df, suite_name) {
  metrics <- c("median_sec", "min_sec", "max_sec", "inner_reps", "median_per_eval_sec")
  out <- do.call(
    rbind,
    lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, , drop = FALSE]
      data.frame(
        section = "speed",
        suite = suite_name,
        label = as.character(row$label[[1]]),
        metric = metrics,
        value = as.numeric(unlist(row[1, metrics], use.names = FALSE)),
        mode = as.character(row$mode[[1]] %||% ""),
        stringsAsFactors = FALSE
      )
    })
  )
  out
}

nested_script <- file.path("dev", "scripts", "bench_nested_integral_baseline.R")
examples_script <- file.path("dev", "scripts", "bench_examples_models_likelihood.R")
golden_script <- file.path("dev", "scripts", "check_loglik_golden.R")

old_nested_bytes <- capture_file_bytes(nested_file)
old_examples_bytes <- capture_file_bytes(examples_file)
on.exit({
  restore_file_bytes(nested_file, old_nested_bytes)
  restore_file_bytes(examples_file, old_examples_bytes)
}, add = TRUE)

cat("Running nested benchmark...\n")
nested_run <- run_script(nested_script)
if (nested_run$status != 0L) {
  cat(paste(nested_run$output, collapse = "\n"), "\n")
  stop("Nested benchmark failed")
}

cat("Running examples benchmark...\n")
examples_run <- run_script(examples_script)
if (examples_run$status != 0L) {
  cat(paste(examples_run$output, collapse = "\n"), "\n")
  stop("Examples benchmark failed")
}

cat("Running golden correctness check...\n")
golden_run <- run_script(golden_script)
if (golden_run$status != 0L) {
  cat(paste(golden_run$output, collapse = "\n"), "\n")
  stop("Golden correctness check failed")
}

if (!file.exists(nested_file)) {
  stop("Missing nested benchmark output: ", nested_file)
}
if (!file.exists(examples_file)) {
  stop("Missing examples benchmark output: ", examples_file)
}

nested_df <- utils::read.csv(nested_file, stringsAsFactors = FALSE)
examples_df <- utils::read.csv(examples_file, stringsAsFactors = FALSE)

speed_rows <- rbind(
  flatten_speed(nested_df, "nested"),
  flatten_speed(examples_df, "examples_models")
)
correctness_rows <- parse_golden_output(golden_run$output)

git_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)[1]),
  error = function(e) ""
)
run_ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

combined <- rbind(speed_rows, correctness_rows)
combined$run_timestamp_utc <- run_ts
combined$git_commit <- git_commit
combined <- combined[, c(
  "run_timestamp_utc", "git_commit",
  "section", "suite", "label", "metric", "value", "mode"
)]

utils::write.csv(combined, out_file, row.names = FALSE)

cat("\nWrote centralized benchmark output to", out_file, "\n")
cat("Rows:", nrow(combined), "\n")
