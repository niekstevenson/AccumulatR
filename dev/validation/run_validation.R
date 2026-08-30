pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
source("dev/validation/helpers.R")
source("dev/validation/cases.R")

arguments <- commandArgs(trailingOnly = TRUE)
adversarial <- "--adversarial" %in% arguments
standard_cases <- validation_cases()
cases <- if (adversarial) {
  all_cases <- validation_cases(include_adversarial = TRUE)
  all_cases[setdiff(names(all_cases), names(standard_cases))]
} else {
  standard_cases
}
selected <- sub("^--case=", "", grep("^--case=", arguments, value = TRUE))
if (length(selected)) {
  cases <- cases[intersect(names(cases), selected)]
}
if (!length(cases)) {
  stop("no validation cases selected", call. = FALSE)
}

results <- do.call(rbind, lapply(names(cases), function(case_name) {
  result <- cases[[case_name]]()
  result$model_name <- case_name
  result
}))
results <- results[, c(
  "model_name", "check_id", "description", "engine", "manual",
  "abs_diff", "tolerance", "passed"
)]
row.names(results) <- NULL
print(results, row.names = FALSE)

summary <- aggregate(passed ~ model_name, results, all)
summary$n_checks <- as.integer(table(results$model_name)[summary$model_name])
summary$n_failed <- as.integer(tapply(!results$passed, results$model_name, sum)[summary$model_name])
cat("\nModel summary\n")
print(summary, row.names = FALSE)
cat(sprintf(
  "\nOverall: %d/%d checks passed across %d models\n",
  sum(results$passed), nrow(results), nrow(summary)
))

if (!all(results$passed)) {
  cat("\nFailing checks\n")
  print(
    results[!results$passed, c("model_name", "check_id", "abs_diff", "tolerance")],
    row.names = FALSE
  )
  quit(status = 1L)
}
