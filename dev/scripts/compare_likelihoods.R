# Helper utilities for comparing native vs. R likelihood paths.

if (!"AccumulatR" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    stop("Load the AccumulatR package (or run devtools::load_all()) before using this script.")
  }
}

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

with_native_flags <- function(node_eval, param_rows, expr) {
  old_opts <- options(
    uuber.use.native.node.eval = node_eval,
    uuber.use.native.param.rows = param_rows
  )
  on.exit(options(old_opts), add = TRUE)
  force(expr)
}

compare_response_probabilities <- function(model, component = NULL) {
  resp_r <- with_native_flags(FALSE, FALSE, {
    response_probabilities(model, component = component)
  })

  resp_cpp <- with_native_flags(TRUE, TRUE, {
    response_probabilities(model, component = component)
  })

  diff <- resp_cpp - resp_r
  list(
    R = resp_r,
    native = resp_cpp,
    diff = diff,
    max_abs_diff = max(abs(diff), na.rm = TRUE)
  )
}

compare_log_likelihood <- function(structure,
                                   params_df,
                                   data_df,
                                   component_weights = NULL,
                                   prep = NULL,
                                   trial_plan = NULL) {
  structure_obj <- AccumulatR:::.as_generator_structure(structure)

  r_path <- with_native_flags(FALSE, FALSE, {
    log_likelihood_from_params(
      structure_obj,
      params_df,
      data_df,
      component_weights = component_weights,
      prep = prep,
      trial_plan = trial_plan
    )
  })

  native_env <- new.env(parent = emptyenv())
  cpp_path <- with_native_flags(TRUE, TRUE, {
    prep_base <- prep %||% AccumulatR:::.prepare_model_for_likelihood(structure_obj$model_spec)
    if (is.null(prep_base[[".runtime"]]) || is.null(prep_base$.runtime$cache_bundle)) {
      prep_base <- AccumulatR:::.prep_set_cache_bundle(
        prep_base,
        AccumulatR:::.build_likelihood_cache_bundle(prep_base)
      )
    }
    plan_base <- trial_plan %||% AccumulatR:::.likelihood_build_trial_plan(
      structure_obj,
      params_df,
      prep_base,
      build_override_ptr = TRUE,
      keep_trial_rows = TRUE,
      keep_component_rows = TRUE
    )
    data_df$trial <- data_df$trial
    trial_ids <- unique(data_df$trial)
    data_row_indices <- split(seq_len(nrow(data_df)), as.character(data_df$trial %||% NA))
    native_env$batch_result <- AccumulatR:::.native_loglikelihood_batch(
      structure = structure_obj,
      prep = prep_base,
      plan = plan_base,
      trial_ids = trial_ids,
      data_df = data_df,
      data_row_indices = data_row_indices,
      component_weights = component_weights,
      params_df = params_df,
      use_buffer = FALSE
    )
    log_likelihood_from_params(
      structure_obj,
      params_df,
      data_df,
      component_weights = component_weights,
      prep = prep_base,
      trial_plan = plan_base
    )
  })

  native_used <- !is.null(native_env$batch_result)

  list(
    R = r_path,
    native = cpp_path,
    native_batch_used = native_used,
    loglik_diff = cpp_path$loglik - r_path$loglik,
    per_trial_diff = cpp_path$per_trial - r_path$per_trial
  )
}

bench_likelihood <- function(expr, times = 3, warmup = TRUE) {
  stopifnot(is.expression(expr) || is.call(expr) || is.language(expr))
  if (isTRUE(warmup)) {
    force(eval(expr, envir = parent.frame()))
  }
  timings <- numeric(times)
  for (i in seq_len(times)) {
    start <- proc.time()[["elapsed"]]
    force(eval(expr, envir = parent.frame()))
    timings[[i]] <- proc.time()[["elapsed"]] - start
  }
  timings
}

compare_speed <- function(model = NULL,
                          response_component = NULL,
                          structure = NULL,
                          params_df = NULL,
                          data_df = NULL,
                          component_weights = NULL,
                          prep = NULL,
                          trial_plan = NULL,
                          n = 3) {
  if (!is.null(model)) {
    timings_r <- bench_likelihood(quote(with_native_flags(FALSE, FALSE, {
      response_probabilities(model, component = response_component)
    })), times = n)
    timings_cpp <- bench_likelihood(quote(with_native_flags(TRUE, TRUE, {
      response_probabilities(model, component = response_component)
    })), times = n)
    timings <- data.frame(
      path = c("R", "native"),
      seconds = c(mean(timings_r), mean(timings_cpp))
    )
    timings$relative <- timings$seconds / min(timings$seconds, na.rm = TRUE)
    return(list(
      task = "response_probabilities",
      times = n,
      timings = timings
    ))
  }
  if (!is.null(structure) && !is.null(params_df) && !is.null(data_df)) {
    timings_r <- bench_likelihood(quote(with_native_flags(FALSE, FALSE, {
      log_likelihood_from_params(
        structure,
        params_df,
        data_df,
        component_weights = component_weights,
        prep = prep,
        trial_plan = trial_plan
      )
    })), times = n)
    timings_cpp <- bench_likelihood(quote(with_native_flags(TRUE, TRUE, {
      log_likelihood_from_params(
        structure,
        params_df,
        data_df,
        component_weights = component_weights,
        prep = prep,
        trial_plan = trial_plan
      )
    })), times = n)
    timings <- data.frame(
      path = c("R", "native"),
      seconds = c(mean(timings_r), mean(timings_cpp))
    )
    timings$relative <- timings$seconds / min(timings$seconds, na.rm = TRUE)
    return(list(
      task = "log_likelihood_from_params",
      times = n,
      timings = timings
    ))
  }
  stop("Provide either a model (for response comparisons) or structure + params/data frames (for likelihood benchmarks).")
}

compare_response_suite <- function(model_list,
                                   example_ids = names(model_list),
                                   component = NULL) {
  if (length(example_ids) == 0) stop("No example ids supplied")
  rows <- lapply(example_ids, function(id) {
    model <- model_list[[id]]
    if (is.null(model)) stop(sprintf("Model '%s' not found in supplied list", id))
    cmp <- compare_response_probabilities(model, component = component)
    data.frame(
      model = id,
      max_abs_diff = cmp$max_abs_diff,
      R = I(list(cmp$R)),
      native = I(list(cmp$native)),
      diff = I(list(cmp$diff)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_param_table_from_structure <- function(structure, n_trials = 1000L) {
  acc_lookup <- structure$accumulators
  if (is.null(acc_lookup) || nrow(acc_lookup) == 0L) {
    stop("Structure must contain accumulator definitions")
  }
  comp_ids <- structure$components$component_id
  component_label <- if (length(comp_ids) == 1L) comp_ids[[1]] else NA_character_
  n_acc <- nrow(acc_lookup)
  base <- data.frame(
    trial = rep(seq_len(n_trials), each = n_acc),
    accumulator_id = rep(acc_lookup$accumulator_id, times = n_trials),
    accumulator = rep(acc_lookup$accumulator_index, times = n_trials),
    component = rep(component_label, times = n_trials * n_acc),
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
      base[[param_name]] <- rep(base_vals, times = n_trials)
    }
  }
  base
}

simulate_data_from_params_table <- function(structure,
                                            param_table,
                                            seed = NULL) {
  sim <- simulate_trials_from_params(structure, param_table, seed = seed)
  data.frame(
    trial = sim$trial,
    outcome = sim$outcome,
    rt = sim$rt,
    stringsAsFactors = FALSE
  )
}

run_param_table_benchmark <- function(model_spec,
                                      n_trials = 1000L,
                                      seed = 2025,
                                      bench_reps = 3) {
  structure <- build_generator_structure(model_spec)
  param_table <- build_param_table_from_structure(structure, n_trials = n_trials)
  data_df <- simulate_data_from_params_table(structure, param_table, seed = seed)
  prep_base <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  if (is.null(prep_base[[".runtime"]]) || is.null(prep_base$.runtime$cache_bundle)) {
    prep_base <- AccumulatR:::.prep_set_cache_bundle(
      prep_base,
      AccumulatR:::.build_likelihood_cache_bundle(prep_base)
    )
  }
  plan_base <- AccumulatR:::.likelihood_build_trial_plan(
    structure,
    param_table,
    prep_base,
    build_override_ptr = TRUE,
    keep_trial_rows = TRUE,
    keep_component_rows = TRUE
  )
  loglik_cmp <- compare_log_likelihood(
    structure,
    param_table,
    data_df,
    prep = prep_base,
    trial_plan = plan_base
  )
  speed_cmp <- compare_speed(
    structure = structure,
    params_df = param_table,
    data_df = data_df,
    prep = prep_base,
    trial_plan = plan_base,
    n = bench_reps
  )
  list(
    structure = structure,
    params = param_table,
    data = data_df,
    loglik = loglik_cmp,
    speed = speed_cmp
  )
}

# Example usage (run inside dev/scripts/test.R or the console):
# mdl <- new_api_examples[[1]]
# compare_response_probabilities(mdl)
# structure <- build_generator_structure(mdl)
# params_df <- ...; data_df <- ...
# compare_log_likelihood(structure, params_df, data_df)
