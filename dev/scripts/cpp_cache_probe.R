rm(list = ls())
# Cache probe script for native likelihood path.

if (!"AccumulatR" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    stop("Load AccumulatR (or run devtools::load_all()) before sourcing this script.")
  }
}

source("dev/examples/new_API.R")

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

simulate_data_from_params_table <- function(structure, param_table, seed = NULL) {
  sim <- simulate_trials_from_params(structure, param_table, seed = seed, keep_detail = FALSE)
  data.frame(
    trial = sim$trial,
    outcome = sim$outcome,
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
}

extract_native_stats <- function(prep) {
  stats <- AccumulatR:::likelihood_cache_stats(prep)
  native <- stats$native
  r_stats <- stats$r %||% list()
  struct_hits <- r_stats$struct_hits %||% 0
  struct_misses <- r_stats$struct_misses %||% 0
  if (is.null(native)) {
    return(data.frame(
      struct_hits = struct_hits,
      struct_misses = struct_misses,
      na_hits = NA_real_,
      na_misses = NA_real_,
      scratch_hits = NA_real_,
      scratch_misses = NA_real_,
      context_builds = NA_real_,
      context_reuses = NA_real_
    ))
  }
  data.frame(
    struct_hits = struct_hits,
    struct_misses = struct_misses,
    na_hits = native$na_hits %||% 0,
    na_misses = native$na_misses %||% 0,
    scratch_hits = native$scratch_hits %||% 0,
    scratch_misses = native$scratch_misses %||% 0,
    context_builds = native$context_builds %||% 0,
    context_reuses = native$context_reuses %||% 0
  )
}

perturb_core_params <- function(core_params, scale = 0.02, seed = NULL) {
  if (is.null(core_params)) stop("core_params must be supplied")
  if (!is.null(seed)) set.seed(seed)
  vals <- core_params
  if (is.null(names(vals))) {
    stop("core_params must be a named vector or list")
  }
  for (nm in names(vals)) {
    if (!is.numeric(vals[[nm]])) next
    jitter <- stats::runif(length(vals[[nm]]), 1 - scale, 1 + scale)
    vals[[nm]] <- vals[[nm]] * jitter
  }
  vals
}

run_cpp_evals <- function(example_id,
                          n_trials = 200L,
                          reps = 3L,
                          seed = 1234L) {
  model_spec <- new_api_examples[[example_id]]
  if (is.null(model_spec)) stop(sprintf("Unknown example id '%s'", example_id))
  core_params <- new_api_example_params[[example_id]]
  if (is.null(core_params)) stop(sprintf("Missing parameter vector for '%s'", example_id))

  options_stack <- options(
    uuber.use.native.node.eval = TRUE,
    uuber.use.native.param.rows = TRUE
  )
  on.exit(options(options_stack), add = TRUE)

  structure <- build_generator_structure(model_spec)
  param_table <- build_params_df(model_spec, core_params, n_trials = n_trials)
  data_df <- simulate_data_from_params_table(structure, param_table, seed = seed)

  stats_rows <- vector("list", reps)
  for (i in seq_len(reps)) {
    perturbed_core <- perturb_core_params(core_params, seed = seed + i)
    iter_table <- build_params_df(model_spec, perturbed_core, n_trials = n_trials)
    context <- build_likelihood_context(structure, iter_table)
    prep <- context$prep
    AccumulatR:::likelihood_cache_reset_stats(prep)
    invisible(log_likelihood_from_context(
      context,
      data_df = data_df
    ))
    stats_rows[[i]] <- cbind(
      data.frame(example = example_id, iteration = i),
      extract_native_stats(prep)
    )
  }
  do.call(rbind, stats_rows)
}

example_ids <- c(
  "example_2_stop_mixture",
  "example_3_stop_na",
  "example_6_dual_path"
)

results <- do.call(
  rbind,
  lapply(example_ids, run_cpp_evals, n_trials = 400L, reps = 3L, seed = 2025L)
)

print(results)
