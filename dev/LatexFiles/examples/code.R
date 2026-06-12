# Standalone replication script for the AccumulatR JSS manuscript.

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  this_file <- normalizePath(file.path("dev", "LatexFiles", "examples", "code.R"))
}
example_dir <- dirname(this_file)

source(file.path(example_dir, "01-basic-workflow.R"), local = TRUE)
source(file.path(example_dir, "02-structured-models.R"), local = TRUE)
source(file.path(example_dir, "03-mixtures-and-fitting.R"), local = TRUE)

sessionInfo()
