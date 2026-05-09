args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("dev/validation/compiler_architecture_acceptance.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

suppressPackageStartupMessages({
  library(pkgload)
})
pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)

model_simple_first_of <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_outcome("A", "A") |>
    add_outcome("B", "B") |>
    finalize_model()
}

model_all_of_competitor <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("C", "lognormal") |>
    add_outcome("AB", all_of("A", "B")) |>
    add_outcome("C", "C") |>
    finalize_model()
}

model_shared_gate_competitors <- function() {
  race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("gate", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("C1", all_of("a", "gate")) |>
    add_outcome("C2", all_of("b", "gate")) |>
    add_outcome("D", "d") |>
    finalize_model()
}

model_guarded_overlap <- function() {
  race_spec() |>
    add_accumulator("a", "lognormal") |>
    add_accumulator("b", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("d", "lognormal") |>
    add_outcome("C1", inhibit("a", by = "stop")) |>
    add_outcome("C2", inhibit("b", by = "stop")) |>
    add_outcome("D", "d") |>
    finalize_model()
}

model_stim_selective_stop2 <- function() {
  race_spec() |>
    add_accumulator("A", "lognormal") |>
    add_accumulator("B", "lognormal") |>
    add_accumulator("S1", "lognormal") |>
    add_accumulator("IS", "lognormal") |>
    add_accumulator("S2", "lognormal") |>
    add_outcome(
      "A",
      first_of(
        inhibit("A", by = "S1"),
        all_of("A", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome(
      "B",
      first_of(
        inhibit("B", by = "S1"),
        all_of("B", "S1", inhibit("IS", by = "S2"))
      )
    ) |>
    add_outcome("STOP", all_of("S1", "S2")) |>
    add_component("go", members = c("A", "B")) |>
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    set_parameters(list(
      m_go = c("A.m", "B.m"),
      s_go = c("A.s", "B.s"),
      t0_go = c("A.t0", "B.t0")
    )) |>
    finalize_model()
}

compiler_metrics <- function(structure) {
  as.list(complexity_metrics(make_context(structure))$total)
}

make_check <- function(model_name, metric, observed, relation, expected, passed) {
  data.frame(
    model_name = model_name,
    metric = metric,
    observed = unname(as.numeric(observed)),
    relation = relation,
    expected = unname(as.numeric(expected)),
    passed = unname(as.logical(passed)),
    stringsAsFactors = FALSE
  )
}

check_equal <- function(model_name, metrics, metric, expected) {
  observed <- metrics[[metric]]
  make_check(model_name, metric, observed, "==", expected, identical(observed, expected))
}

check_at_most <- function(model_name, metrics, metric, expected) {
  observed <- metrics[[metric]]
  make_check(model_name, metric, observed, "<=", expected, observed <= expected)
}

acceptance_cases <- list(
  simple_first_of = list(
    build = model_simple_first_of,
    max_integral_depth = 0L,
    max_integral_nodes = 0L,
    max_symbolic_cells = 4L
  ),
  all_of_competitor = list(
    build = model_all_of_competitor,
    max_integral_depth = 0L,
    max_integral_nodes = 0L,
    max_symbolic_cells = 7L
  ),
  shared_gate_competitors = list(
    build = model_shared_gate_competitors,
    max_integral_depth = 1L,
    max_integral_nodes = 4L,
    max_symbolic_cells = 31L
  ),
  guarded_overlap = list(
    build = model_guarded_overlap,
    max_integral_depth = 2L,
    max_integral_nodes = 14L,
    max_symbolic_cells = 21L
  ),
  stim_selective_stop2 = list(
    build = model_stim_selective_stop2,
    max_integral_depth = 3L,
    max_integral_nodes = 83L,
    max_symbolic_cells = 60L
  )
)

results <- do.call(rbind, lapply(names(acceptance_cases), function(model_name) {
  spec <- acceptance_cases[[model_name]]
  metrics <- compiler_metrics(spec$build())
  checks <- list(
    check_equal(model_name, metrics, "negative_symbolic_cells", 0L),
    check_equal(model_name, metrics, "overlapping_symbolic_cell_pairs", 0L),
    check_equal(model_name, metrics, "generic_integral_kernels", 0L),
    check_at_most(model_name, metrics, "max_integral_depth", spec$max_integral_depth),
    check_at_most(model_name, metrics, "integral_nodes", spec$max_integral_nodes),
    check_at_most(model_name, metrics, "symbolic_cells", spec$max_symbolic_cells),
    make_check(
      model_name,
      "expr_relation_atoms_diagnostic",
      metrics[["expr_relation_atoms"]],
      "reported",
      NA_real_,
      TRUE
    )
  )
  do.call(rbind, checks)
}))

row.names(results) <- NULL
print(results, row.names = FALSE)

summary_df <- aggregate(
  passed ~ model_name,
  data = results[results$relation != "reported", , drop = FALSE],
  FUN = all
)
summary_df$n_checks <- as.integer(table(
  results$model_name[results$relation != "reported"]
)[summary_df$model_name])
summary_df$n_failed <- as.integer(tapply(
  !results$passed[results$relation != "reported"],
  results$model_name[results$relation != "reported"],
  sum
)[summary_df$model_name])

cat("\nArchitecture acceptance summary\n")
print(summary_df, row.names = FALSE)

checked <- results$relation != "reported"
all_passed <- all(results$passed[checked])
cat(sprintf(
  "\nOverall: %d/%d structural checks passed across %d models\n",
  sum(results$passed[checked]),
  sum(checked),
  nrow(summary_df)
))

if (!all_passed) {
  cat("\nFailing structural checks\n")
  print(results[checked & !results$passed, , drop = FALSE], row.names = FALSE)
  quit(status = 1L)
}
