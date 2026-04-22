#!/usr/bin/env Rscript

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

suppressPackageStartupMessages({
  library(pkgload)
})

load_all(repo_root, quiet = TRUE, helpers = FALSE)
source(file.path("dev", "examples", "new_API.R"))

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_csv <- file.path(out_dir, "benchmark_all_examples_capped_100ms.csv")
helper <- file.path("dev", "scripts", "bench_one_example.R")

n_trials <- 250L
cap_ms <- 100
target_batch_ms <- 400
model_timeout_sec <- 20

labels <- names(new_api_examples)
results <- vector("list", length(labels))

for (i in seq_along(labels)) {
  label <- labels[[i]]
  cat(sprintf("[%d/%d] %s\n", i, length(labels), label))
  flush.console()

  tmp_csv <- tempfile("bench_one_", fileext = ".csv")
  cmd <- c(
    helper,
    label,
    as.character(n_trials),
    as.character(cap_ms),
    as.character(target_batch_ms),
    tmp_csv
  )

  out <- tryCatch(
    system2(
      "Rscript",
      cmd,
      stdout = TRUE,
      stderr = TRUE,
      timeout = model_timeout_sec
    ),
    error = function(e) structure(character(), status = 1L, errmsg = conditionMessage(e))
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }

  if (file.exists(tmp_csv)) {
    row <- utils::read.csv(tmp_csv, stringsAsFactors = FALSE)
  } else {
    timed_out <- any(grepl("timed out", out, ignore.case = TRUE)) ||
      identical(as.integer(status), 124L)
    msg <- if (length(out) > 0L) paste(out, collapse = " | ") else "subprocess failed"
    row <- data.frame(
      label = label,
      per_eval_ms = NA_real_,
      batch_ms = NA_real_,
      inner_reps = NA_integer_,
      capped = if (timed_out) TRUE else NA,
      status = if (timed_out) "timeout" else "error",
      error = msg,
      stringsAsFactors = FALSE
    )
  }

  results[[i]] <- row

  partial <- do.call(rbind, results[seq_len(i)])
  partial$rank <- NA_integer_
  ok <- partial$status == "ok"
  partial$rank[ok] <- rank(partial$per_eval_ms[ok], ties.method = "first")
  partial <- partial[, c(
    "rank", "label", "per_eval_ms", "batch_ms", "inner_reps",
    "capped", "status", "error"
  )]
  partial <- partial[order(is.na(partial$rank), partial$rank, partial$label), ]
  utils::write.csv(partial, out_csv, row.names = FALSE)
}

results <- do.call(rbind, results)
results$rank <- NA_integer_
ok <- results$status == "ok"
results$rank[ok] <- rank(results$per_eval_ms[ok], ties.method = "first")
results <- results[, c(
  "rank", "label", "per_eval_ms", "batch_ms", "inner_reps",
  "capped", "status", "error"
)]
results <- results[order(is.na(results$rank), results$rank, results$label), ]

utils::write.csv(results, out_csv, row.names = FALSE)

cat("Wrote capped full-suite benchmark to", out_csv, "\n")
cat("Rows:", nrow(results), "\n")
