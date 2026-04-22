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

out_dir <- file.path("dev", "scripts", "scratch_outputs")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cpp_file <- file.path("dev", "scripts", "bench_integral_families.cpp")
bin_file <- file.path(out_dir, "bench_integral_families")
csv_file <- file.path(out_dir, "integral_family_benchmark.csv")

bh_include <- system.file("include", package = "BH")
if (!nzchar(bh_include)) {
  stop("BH include path not found")
}

cmd <- c(
  "-std=c++17",
  "-O3",
  "-DNDEBUG",
  paste0("-I", repo_root),
  paste0("-I", bh_include),
  cpp_file,
  "-o",
  bin_file
)

status <- system2("clang++", cmd)
if (!identical(status, 0L)) {
  stop("Failed to compile integral-family benchmark")
}

status <- system2(bin_file, csv_file)
if (!identical(status, 0L)) {
  stop("Integral-family benchmark failed")
}

cat("Wrote integral-family benchmark to", csv_file, "\n")
