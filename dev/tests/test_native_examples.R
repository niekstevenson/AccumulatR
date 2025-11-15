rm(list = ls())
options(uuber.native.rebuild = TRUE)

base::source("examples/new_API.R")
base::source("R/dist.R")
base::source("R/utils.R")
base::source("R/pool_math.R")
base::source("R/model_tables.R")
base::source("R/generator_new.R")
base::source("R/likelihood_cache.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_prep.R")
base::source("R/likelihood_primitives.R")
base::source("R/likelihood_kernels.R")
base::source("R/likelihood_integrate.R")
base::source("R/likelihood_param_interface.R")
base::source("R/super_large_likelihood.R")

build_param_table <- function(structure, n_trials) {
  acc_lookup <- structure$accumulators
  param_table <- data.frame(
    trial = rep(seq_len(n_trials), each = nrow(acc_lookup)),
    accumulator_id = rep(acc_lookup$accumulator_id, times = n_trials),
    accumulator = rep(acc_lookup$accumulator_index, times = n_trials),
    component = ifelse(
      length(structure$components$component_id) == 1L,
      structure$components$component_id,
      NA_character_
    ),
    role = rep(acc_lookup$role, times = n_trials),
    onset = rep(acc_lookup$onset, times = n_trials),
    q = rep(acc_lookup$q, times = n_trials),
    stringsAsFactors = FALSE
  )
  param_names <- unique(unlist(lapply(acc_lookup$params, names)))
  if (length(param_names) > 0L) {
    for (param_name in param_names) {
      base_vals <- vapply(acc_lookup$params, function(p) {
        val <- p[[param_name]]
        if (is.null(val)) NA_real_ else as.numeric(val)
      }, numeric(1))
      param_table[[param_name]] <- rep(base_vals, times = n_trials)
    }
  }
  param_table
}

run_example_check <- function(example_idx, seed = 2025, n_trials = 400L,
                              tol = 5e-3, t0_shift = NULL) {
  set.seed(seed)
  model_spec <- new_api_examples[[example_idx]]
  if (!is.null(t0_shift)) {
    model_spec$accumulators <- lapply(model_spec$accumulators, function(acc) {
      params <- acc$params %||% list()
      params$t0 <- t0_shift
      acc$params <- params
      acc
    })
  }
  model_tables <- model_to_tables(model_spec)
  structure <- build_generator_structure(model_spec)
  sim_data <- simulate_model(model_tables, n_trials = n_trials)
  param_table <- build_param_table(structure, n_trials)
  native_bundle <- likelihood_build_native_bundle(model_spec = model_spec)
  prep <- native_bundle$prep
  trial_plan <- .likelihood_build_trial_plan(structure, param_table, prep)
trial_plan_buffer <- .likelihood_build_trial_plan(structure,
                                                 param_table,
                                                 prep)
  system.time(
    ll_r <- compute_loglik(model_tables, sim_data)
  )

  ll_params <- log_likelihood_from_params(structure,
                                          param_table,
                                          sim_data,
                                          prep = prep,
                                          trial_plan = trial_plan)
  ll_buffer <- log_likelihood_from_params_buffer(structure,
                                                 param_table,
                                                 sim_data,
                                                 prep = prep,
                                                 trial_plan = trial_plan_buffer)
  if (!isTRUE(all.equal(as.numeric(ll_r), as.numeric(ll_params$loglik), tolerance = tol))) {
    stop(sprintf("Example %d: native log-likelihood %.12f differs from baseline %.12f",
                 example_idx, as.numeric(ll_params$loglik), as.numeric(ll_r)))
  }
  if (!isTRUE(all.equal(as.numeric(ll_r), as.numeric(ll_buffer$loglik), tolerance = tol))) {
    stop(sprintf("Example %d (buffer): native log-likelihood %.12f differs from baseline %.12f",
                 example_idx, as.numeric(ll_buffer$loglik), as.numeric(ll_r)))
  }
  invisible(TRUE)
}

debug(run_example_check)
examples_to_check <- c(2, 3, 6)
for (idx in examples_to_check) {
  run_example_check(idx)
  cat(sprintf("Example %d native likelihood matches baseline\n", idx))
}
run_example_check(2, t0_shift = 0.05)
cat("Example 2 native likelihood matches baseline with t0 shift\n")
