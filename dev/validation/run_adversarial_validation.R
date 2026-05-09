args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("dev/validation/run_adversarial_validation.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

source(file.path(repo_root, "dev", "validation", "helpers.R"))
source(file.path(repo_root, "dev", "validation", "cases.R"))

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

suppressPackageStartupMessages({
  library(pkgload)
})
pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
validation_source_rebuild(repo_root)

adversarial_case_names <- c(
  "oracle_repeated_shared_gate_six_way",
  "oracle_deep_composite_blocker",
  "oracle_pool_k2_shared_gate_guard"
)
case_args <- sub(
  "^--case=",
  "",
  grep("^--case=", commandArgs(trailingOnly = TRUE), value = TRUE)
)
if (length(case_args) > 0L) {
  adversarial_case_names <- intersect(adversarial_case_names, case_args)
}
if (length(adversarial_case_names) == 0L) {
  stop("no adversarial validation cases selected", call. = FALSE)
}
cases <- validation_cases(include_adversarial = TRUE)
cases <- cases[adversarial_case_names]

results <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    out <- cases[[case_name]]()
    out$model_name <- case_name
    out
  })
)

results <- results[, c(
  "model_name",
  "check_id",
  "description",
  "engine",
  "manual",
  "abs_diff",
  "tolerance",
  "passed"
)]
row.names(results) <- NULL

print(results, row.names = FALSE)

summary_df <- aggregate(
  passed ~ model_name,
  data = results,
  FUN = all
)
summary_df$n_checks <- as.integer(table(results$model_name)[summary_df$model_name])
summary_df$n_failed <- as.integer(tapply(!results$passed, results$model_name, sum)[summary_df$model_name])

cat("\nAdversarial model summary\n")
print(summary_df, row.names = FALSE)

all_passed <- all(results$passed)
cat(sprintf(
  "\nOverall: %d/%d checks passed across %d adversarial models\n",
  sum(results$passed),
  nrow(results),
  nrow(summary_df)
))

if (!all_passed) {
  failures <- results[!results$passed, c("model_name", "check_id", "abs_diff", "tolerance")]
  cat("\nFailing adversarial checks\n")
  print(failures, row.names = FALSE)
  quit(status = 1L)
}
