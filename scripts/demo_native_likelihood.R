rm(list = ls())

# Load example specs and likelihood machinery ---------------------------------
source("examples/new_API.R")
source("R/dist.R")
source("R/utils.R")
source("R/pool_math.R")
source("R/model_tables.R")
source("R/generator_new.R")
source("R/likelihood_cache.R")
source("R/likelihood_common.R")
source("R/likelihood_prep.R")
source("R/likelihood_primitives.R")
source("R/likelihood_kernels.R")
source("R/likelihood_integrate.R")
source("R/likelihood_param_interface.R")
source("R/super_large_likelihood.R")

set.seed(1234)

# -----------------------------------------------------------------------------#
# 1. Build a model, simulate data, and capture its parameter table
# -----------------------------------------------------------------------------#
model_spec <- new_api_examples[[1]]
model_tables <- model_to_tables(model_spec)
structure <- build_generator_structure(model_spec)
sim_data <- simulate_model(model_tables, n_trials = 10000)
cat(sprintf("Simulated %d trials for example '%s'\n",
            nrow(sim_data),
            names(new_api_examples)[[1]]))

n_trials <- nrow(sim_data)
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

# Shared caches (compiled nodes, pool templates, guard quadrature, proto blob)
# stay in prep$.runtime$cache_bundle. Per-trial variability lives entirely in
# the trial_rows/component_rows data frames, so we never mutate the shared prep.

# -----------------------------------------------------------------------------#
# 2. Prep the likelihood once and export the native proto/context bundle
# -----------------------------------------------------------------------------#
native_bundle <- likelihood_build_native_bundle(model_spec = model_spec)
prep <- native_bundle$prep
proto_blob <- native_bundle$proto
ctx_ptr <- native_bundle$context
cat(sprintf("Proto blob size: %d bytes\n", length(proto_blob)))

# Demonstrate that we can rebuild the context from only the proto
ctx_from_proto <- .lik_native_fn("native_context_from_proto_cpp")(proto_blob)

trial_plan <- .likelihood_build_trial_plan(structure, param_table, prep)

# -----------------------------------------------------------------------------#
# 3. Compare R vs C++ outcome evaluation on a single trial (table-driven)
# -----------------------------------------------------------------------------#
trial <- sim_data[1, , drop = FALSE]
trial_id <- trial$trial
outcome_lbl <- as.character(trial$outcome)
component_lbl <- trial$component
if (length(component_lbl) == 0L || is.na(component_lbl)) {
  component_lbl <- NULL
} else {
  component_lbl <- as.character(component_lbl)
}
component_for_r <- component_lbl %||% "__default__"
rt_val <- as.numeric(trial$rt)
trial_entry <- trial_plan$trials[[as.character(trial_id)]] %||% list()
trial_rows <- trial_entry$rows %||% param_table[param_table$trial == trial_id, , drop = FALSE]
component_entries <- trial_entry$components %||% list()
comp_entry <- component_entries[[component_for_r]] %||% component_entries[["__default__"]] %||% list(rows = NULL)
comp_rows_df <- comp_entry$rows
trial_override_ptr <- trial_entry$override_ptr

# R path: use the standard outcome likelihood helper with the component rows
r_density <- .outcome_likelihood(
  outcome_lbl,
  rt_val,
  prep,
  component_for_r,
  trial_rows = comp_rows_df,
  trial_overrides = trial_override_ptr
)

# Collect competitor node IDs for the native call
competitor_map <- .prep_competitors(prep) %||% list()
competitor_exprs <- competitor_map[[outcome_lbl]] %||% list()
comp_nodes <- lapply(competitor_exprs, function(ex) .expr_lookup_compiled(ex, prep))
comp_ids <- vapply(comp_nodes, function(node) {
  id <- node$id %||% NA_integer_
  as.integer(id)
}, integer(1))
comp_ids <- comp_ids[!is.na(comp_ids)]

expr <- prep$outcomes[[outcome_lbl]]$expr
node_id <- attr(expr, ".lik_id", exact = TRUE)

native_density_fn <- .lik_native_fn("native_density_with_competitors_overrides_cpp")
native_density <- native_density_fn(
  ctx_from_proto,
  as.integer(node_id),
  as.numeric(rt_val),
  component_lbl,
  integer(0),
  integer(0),
  comp_ids,
  trial_override_ptr
)$density

cat("\nOutcome density comparison (single trial)\n")
cat(sprintf("  R   : %.6f\n", as.numeric(r_density)))
cat(sprintf("  C++ : %.6f\n", as.numeric(native_density)))

# -----------------------------------------------------------------------------#
# 4. Compare integrated outcome probabilities (still table-driven)
# -----------------------------------------------------------------------------#
deadline <- .get_component_attr(prep, component_for_r, "deadline") %||% prep$default_deadline
upper_limit <- min(deadline, 2)
r_prob <- .integrate_outcome_probability(
  expr = expr,
  prep = prep,
  component = component_for_r,
  upper_limit = upper_limit,
  competitor_exprs = competitor_exprs,
  trial_rows = comp_rows_df
)

native_prob_fn <- .lik_native_fn("native_outcome_probability_overrides_cpp")
native_prob <- native_prob_fn(
  ctx_from_proto,
  as.integer(node_id),
  as.numeric(upper_limit),
  component_lbl,
  integer(0),
  integer(0),
  comp_ids,
  .integrate_rel_tol(),
  .integrate_abs_tol(),
  12L,
  trial_override_ptr
)

cat("\nOutcome probability up to deadline\n")
cat(sprintf("  Upper limit: %.3f\n", as.numeric(upper_limit)))
cat(sprintf("  R   : %.6f\n", as.numeric(r_prob)))
cat(sprintf("  C++ : %.6f\n", as.numeric(native_prob)))

# -----------------------------------------------------------------------------#
# 5. Baseline: table-driven log-likelihood vs original R likelihood
# -----------------------------------------------------------------------------#
library(profvis)
ll_r <- NA_real_
profvis({
  ll_r <<- compute_loglik(model_tables, sim_data)
}, prof_output = "profiles/ll_tables.out")
df_tables <- summaryRprof("profiles/ll_tables.out")
write.csv(df_tables$by.total, file = "profiles/ll_tables_total.csv")
write.csv(df_tables$by.self, file = "profiles/ll_tables_self.csv")


ll_params <- NA_real_
profvis({
  ll_params <<- log_likelihood_from_params(structure, param_table, sim_data,
                                           prep = prep,
                                           trial_plan = trial_plan)
}, prof_output = "profiles/ll_params.out")
df_params <- summaryRprof("profiles/ll_params.out")
write.csv(df_tables$by.total, file = "profiles/ll_params_total.csv")
write.csv(df_tables$by.self, file = "profiles/ll_params_self.csv")

cat(sprintf("\nFull R (model_tables) log-likelihood : %.4f\n", ll_r))
cat(sprintf("Table/log_likelihood_from_params result: %.4f\n", ll_params$loglik))
