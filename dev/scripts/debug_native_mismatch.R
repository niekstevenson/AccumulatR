# Helper script to diagnose native vs. R likelihood differences for a single example.

if (!"AccumulatR" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    stop("Load the AccumulatR package (or run devtools::load_all()) before using this script.")
  }
}

source("dev/examples/new_API.R")
source("dev/scripts/compare_likelihoods.R")

collect_native_debug <- function(example_id = "example_5_timeout_guess",
                                 n_trials = 200L,
                                 seed = 2025) {
  model_spec <- new_api_examples[[example_id]]
  params_vec <- new_api_example_params[[example_id]]
  if (is.null(model_spec) || is.null(params_vec)) {
    stop(sprintf("Example '%s' not found", example_id))
  }
  structure <- build_generator_structure(model_spec)
  params_df <- build_params_df(model_spec, params_vec, n_trials = n_trials)
  data_df <- simulate_data_from_params_table(structure, params_df, seed = seed)

  prep <- AccumulatR:::.prepare_model_for_likelihood(structure$model_spec)
  ctx <- build_likelihood_context(
    structure = structure,
    params_df = params_df,
    data_df = data_df,
    prep = prep
  )
  plan <- ctx$plan
  trial_ids <- ctx$trial_ids
  data_row_indices <- ctx$data_row_indices
  trial_lookup <- split(ctx$data_df, as.character(ctx$data_df$trial))
  trial_frames <- lapply(trial_ids, function(tid) trial_lookup[[as.character(tid)]])

  r_eval <- with(list(), {
    old_opts <- options(uuber.use.native.node.eval = FALSE, uuber.use.native.param.rows = FALSE)
    on.exit(options(old_opts), add = TRUE)
    log_likelihood_from_params(structure, params_df, data_df, prep = prep, trial_plan = plan)
  })

  native_eval <- AccumulatR:::.native_loglikelihood_batch(
    structure = structure,
    prep = ctx$prep,
    plan = plan,
    trial_ids = trial_ids
  )

  builder <- AccumulatR:::.build_native_trial_entries(
    structure = structure,
    prep = ctx$prep,
    plan = plan,
    trial_ids = trial_ids,
    data_df = ctx$data_df,
    data_row_indices = data_row_indices,
    component_weights = NULL
  )

  manual_native <- numeric(length(trial_ids))
  for (i in seq_along(trial_ids)) {
    trial_key <- as.character(trial_ids[[i]])
    trial_entry <- plan$trials[[trial_key]]
    trial_data <- trial_frames[[i]]
    forced_component <- if (!is.null(trial_data$component)) trial_data$component[[1]] else NULL
    outcome_lbl <- trial_data$outcome[[1]]
    rt_val <- trial_data$rt[[1]]
    comp_plan <- AccumulatR:::.native_component_plan(
      structure = structure,
      trial_rows = trial_entry$rows,
      forced_component = forced_component,
      component_weights = NULL,
      component_ids = structure$components$component_id
    )
    manual_native[[i]] <- log(
      AccumulatR:::.native_trial_mixture_eval(
        prep = ctx$prep,
        outcome_label = outcome_lbl,
        rt_val = rt_val,
        component_plan = comp_plan,
        forced_component = forced_component,
        trial_rows = trial_entry$rows,
        guess_donors = builder$entries[[i]]$guess_donors
      )
    )
  }

  data.frame(
    trial = trial_ids,
    outcome = vapply(trial_frames, function(df) df$outcome[[1]], character(1)),
    rt = vapply(trial_frames, function(df) df$rt[[1]], numeric(1)),
    log_r = r_eval$per_trial,
    log_native_batch = native_eval$per_trial,
    log_native_manual = manual_native,
    diff_batch_manual = native_eval$per_trial - manual_native,
    diff_manual_r = manual_native - r_eval$per_trial,
    diff_batch_r = native_eval$per_trial - r_eval$per_trial,
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0 || identical(parent.frame(), .GlobalEnv)) {
  debug_tbl <- collect_native_debug()
  print(head(debug_tbl, 10))
  cat("\nSummary of batch vs. manual native log-density diff:\n")
  print(summary(debug_tbl$diff_batch_manual))
  cat("\nSummary of manual native vs. R log-density diff:\n")
  print(summary(debug_tbl$diff_manual_r))
}
