rm(list = ls())
options(uuber.native.rebuild = TRUE)

base::source("R/super_large_likelihood.R")
base::source("R/generator_new.R")
base::source("R/likelihood_old.R")
base::source("R/likelihood_prep.R")
base::source("R/likelihood_param_interface.R")
base::source("examples/new_API.R")

model_spec <- new_api_examples[[3]]
structure <- build_generator_structure(model_spec)
model_tables <- model_to_tables(model_spec)

set.seed(2025)
n_trials <- 200L
acc_lookup <- structure$accumulators

base_table <- data.frame(
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
    base_table[[param_name]] <- rep(base_vals, times = n_trials)
  }
}

data_sim <- simulate_trials_from_params(structure, base_table, seed = 2025)
data_df <- data.frame(
  trial = data_sim$trial,
  outcome = data_sim$outcome,
  rt = data_sim$rt,
  stringsAsFactors = FALSE
)

fmt_rt <- function(x) if (is.na(x)) "NA" else format(x, digits = 22, scientific = FALSE, trim = TRUE)
component_ids <- structure$components$component_id
component_key <- if (length(component_ids) > 0) component_ids[[1]] else "__default__"
combined_keys <- sprintf(
  "%s|%s|%s",
  component_key,
  data_df$outcome,
  vapply(data_df$rt, fmt_rt, character(1))
)
na_mask <- is.na(data_df$rt)
na_keys <- combined_keys[na_mask]
na_unique_count <- length(unique(na_keys))
na_trial_count <- sum(na_mask)

cat(sprintf("na trials %d unique %d\n", na_trial_count, na_unique_count))

count_integrals <- function(expr) {
  expr_call <- substitute(expr)
  integral_counter <<- 0L
  trace(".integrate_outcome_probability",
        tracer = quote({ integral_counter <<- integral_counter + 1L }),
        print = FALSE)
  on.exit(untrace(".integrate_outcome_probability"), add = TRUE)
  invisible(eval.parent(expr_call))
  integral_counter
}

count_old <- count_integrals(compute_loglik_old(model_tables, data_df))
count_new_first <- count_integrals(log_likelihood_from_params(structure, base_table, data_df))
count_new_second <- count_integrals(log_likelihood_from_params(structure, base_table, data_df))

cat(sprintf("old=%d new_first=%d new_second=%d\n", count_old, count_new_first, count_new_second))

if (count_old > na_unique_count) {
  stop(sprintf(
    "Legacy interface integrated more than expected: unique=%d observed=%d",
    na_unique_count, count_old
  ))
}

if (count_new_first != na_unique_count) {
  stop(sprintf(
    "Expected param interface to integrate %d times on first pass, saw %d",
    na_unique_count, count_new_first
  ))
}

native_option_old <- getOption("uuber.use.native.param.rows")
options(uuber.use.native.param.rows = FALSE)
legacy_res <- log_likelihood_from_params(structure, base_table, data_df)
options(uuber.use.native.param.rows = TRUE)
native_res <- log_likelihood_from_params(structure, base_table, data_df)
options(uuber.use.native.param.rows = native_option_old)

if (!isTRUE(all.equal(as.numeric(legacy_res$loglik), as.numeric(native_res$loglik), tolerance = 1e-9))) {
  stop("Native param-row likelihood does not match legacy override path")
}
if (!isTRUE(all.equal(as.numeric(legacy_res$per_trial), as.numeric(native_res$per_trial), tolerance = 1e-9))) {
  stop("Per-trial log-likelihood differs between native and override paths")
}
