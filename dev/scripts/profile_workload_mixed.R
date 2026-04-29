#!/usr/bin/env Rscript

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) ""
)
repo_root <- normalizePath(file.path(dirname(script_path %||% "."), "..", ".."),
                           mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "dev", "examples"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

load_accumulatr <- function() {
  if (requireNamespace("AccumulatR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(AccumulatR))
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(pkgload))
  load_all(repo_root, quiet = TRUE, helpers = FALSE)
}

load_accumulatr()
source(file.path("dev", "examples", "new_API.R"))

parse_int_env <- function(name, default, min_value = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < min_value) {
    return(as.integer(default))
  }
  value
}

parse_num_env <- function(name, default, min_value = 0) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value <= min_value) {
    return(default)
  }
  value
}

stop_change_model <- function() {
  structure <- race_spec() |>
    add_accumulator("S", "lognormal") |>
    add_accumulator("stop", "lognormal") |>
    add_accumulator("change", "lognormal") |>
    add_outcome("S", inhibit("S", by = "stop")) |>
    add_outcome("X", all_of("change", "stop")) |>
    add_component("go_only", members = "S", weight = .75) |>
    add_component("go_stop", members = c("S", "stop", "change"), weight = .25) |>
    add_trigger("stop_trigger", members = c("stop", "change"), q = 0.05) |>
    set_mixture_options(mode = "fixed") |>
    set_parameters(list(
      m_go = "S.m",
      m_stop = "stop.m",
      m_change = "change.m",
      s_go = "S.s",
      s_stop = "stop.s",
      s_change = "change.s",
      t0_go = "S.t0",
      t0_change = "change.t0",
      q = c("stop.q", "change.q")
    )) |>
    finalize_model()
  pars <- c(
    m_go = log(0.30), s_go = 0.18, t0_go = 0.00,
    m_stop = log(0.22), s_stop = 0.18,
    m_change = log(0.40), s_change = 0.18, t0_change = 0.00,
    q = 0.05, S.q = 0.00
  )
  list(structure = structure, pars = pars)
}

stim_selective_model <- function() {
  structure <- race_spec() |>
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
    add_outcome(
      "STOP",
      all_of("S1", inhibit("S2", by = "IS")),
      options = list(map_outcome_to = NA_character_)
    ) |>
    add_component("go", members = c("A", "B")) |>
    add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
    set_parameters(list(
      m_go = c("A.m", "B.m"),
      s_go = c("A.s", "B.s"),
      t0_go = c("A.t0", "B.t0")
    )) |>
    finalize_model()
  pars <- c(
    m_go = log(0.30), s_go = 0.18, t0_go = 0.05,
    S1.m = log(0.26), S1.s = 0.18, S1.t0 = 0.00,
    IS.m = log(0.35), IS.s = 0.18, IS.t0 = 0.00,
    S2.m = log(0.32), S2.s = 0.18, S2.t0 = 0.00
  )
  list(structure = structure, pars = pars)
}

models <- list(
  example_1_simple = list(
    structure = new_api_examples[["example_1_simple"]],
    pars = new_api_example_params[["example_1_simple"]]
  ),
  example_2_stop_mixture = list(
    structure = new_api_examples[["example_2_stop_mixture"]],
    pars = new_api_example_params[["example_2_stop_mixture"]]
  ),
  example_3_stop_na = list(
    structure = new_api_examples[["example_3_stop_na"]],
    pars = new_api_example_params[["example_3_stop_na"]]
  ),
  example_5_timeout_guess = list(
    structure = new_api_examples[["example_5_timeout_guess"]],
    pars = new_api_example_params[["example_5_timeout_guess"]]
  ),
  example_6_dual_path = list(
    structure = new_api_examples[["example_6_dual_path"]],
    pars = new_api_example_params[["example_6_dual_path"]]
  ),
  example_7_mixture = list(
    structure = new_api_examples[["example_7_mixture"]],
    pars = new_api_example_params[["example_7_mixture"]]
  ),
  example_10_exclusion = list(
    structure = new_api_examples[["example_10_exclusion"]],
    pars = new_api_example_params[["example_10_exclusion"]]
  ),
  example_16_guard_tie_simple = list(
    structure = new_api_examples[["example_15_guard_tie_simple"]],
    pars = new_api_example_params[["example_15_guard_tie_simple"]]
  ),
  example_21_simple_q = list(
    structure = new_api_examples[["example_21_simple_q"]],
    pars = new_api_example_params[["example_21_simple_q"]]
  ),
  example_22_shared_q = list(
    structure = new_api_examples[["example_22_shared_q"]],
    pars = new_api_example_params[["example_22_shared_q"]]
  ),
  example_23_ranked_chain = local({
    structure <- race_spec(n_outcomes = 2L) |>
      add_accumulator("a", "lognormal") |>
      add_accumulator("b", "lognormal", onset = after("a")) |>
      add_outcome("A", "a") |>
      add_outcome("B", "b") |>
      finalize_model()
    pars <- c(a.m = log(0.30), a.s = 0.16, b.m = log(0.22), b.s = 0.16)
    list(structure = structure, pars = pars)
  }),
  stop_change_shared_trigger = stop_change_model(),
  stim_selective_stop = stim_selective_model()
)

balanced_labels <- c(
  "example_1_simple",
  "example_5_timeout_guess",
  "example_21_simple_q",
  "example_22_shared_q",
  "example_2_stop_mixture",
  "example_7_mixture",
  "example_10_exclusion",
  "example_23_ranked_chain",
  "stop_change_shared_trigger"
)

case_env <- Sys.getenv("ACCUMULATR_PROFILE_CASES", "")
case_labels <- if (nzchar(case_env)) {
  trimws(strsplit(case_env, ",", fixed = TRUE)[[1]])
} else {
  balanced_labels
}
case_labels <- case_labels[nzchar(case_labels)]
unknown <- setdiff(case_labels, names(models))
if (length(unknown) > 0L) {
  stop("Unknown profile case(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

n_trials <- parse_int_env("ACCUMULATR_PROFILE_TRIALS", 50L)
workload_seconds <- parse_num_env(
  "ACCUMULATR_PROFILE_WORKLOAD_SECONDS",
  parse_num_env("ACCUMULATR_PROFILE_DURATION", 20)
)

make_case <- function(mod, n_trials) {
  params_df <- build_param_matrix(mod$structure, mod$pars, n_trials = n_trials)
  data_df <- simulate(mod$structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(mod$structure, data_df)
  ctx <- make_context(mod$structure)
  params_slim <- build_param_matrix(mod$structure, mod$pars, trial_df = prepared)
  value <- as.numeric(log_likelihood(ctx, prepared, params_slim))
  if (length(value) != 1L || !is.finite(value)) {
    stop("Non-finite log-likelihood during profile setup", call. = FALSE)
  }
  list(ctx = ctx, prepared = prepared, params = params_slim)
}

cases <- lapply(models[case_labels], make_case, n_trials = n_trials)
counts <- setNames(integer(length(cases)), names(cases))

start_file <- Sys.getenv("ACCUMULATR_PROFILE_START_FILE", "")
end_file <- Sys.getenv("ACCUMULATR_PROFILE_END_FILE", "")
if (nzchar(start_file)) {
  writeLines("start", start_file)
}
on.exit({
  if (nzchar(end_file)) {
    writeLines("end", end_file)
  }
}, add = TRUE)

deadline <- proc.time()[["elapsed"]] + workload_seconds
repeat {
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    value <- as.numeric(log_likelihood(case$ctx, case$prepared, case$params))
    if (length(value) != 1L || !is.finite(value)) {
      stop("Non-finite log-likelihood during profile loop", call. = FALSE)
    }
    counts[[case_name]] <- counts[[case_name]] + 1L
  }
  if (proc.time()[["elapsed"]] >= deadline) {
    break
  }
}

cat("Mixed profile workload completed\n")
cat("Trials per case:", n_trials, "\n")
print(counts)
