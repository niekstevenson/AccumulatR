# Archived helpers that are no longer part of the active package API.
# These were unused in the current C++-first pipeline but kept here
# for reference while the old R paths are retired.

# Expression constructors (previously exported)
expr_event <- function(source_id, k = NULL) {
  list(kind = "event", source = source_id, k = k)
}

expr_and <- function(...) list(kind = "and", args = list(...))
expr_or  <- function(...) list(kind = "or",  args = list(...))
expr_not <- function(expr) list(kind = "not", arg = expr)

# Likelihood/plan helpers no longer referenced
likelihood_reset_cache <- function(prep) {
  if (is.null(prep)) return(prep)
  runtime <- prep[[".runtime"]]
  if (is.null(runtime)) return(prep)
  .prep_set_cache_bundle(prep, .build_likelihood_cache_bundle(prep))
}

.inspect_likelihood_plan <- function(prep, include_cache = FALSE) {
  comp <- .prep_expr_compiled(prep)
  if (is.null(comp)) {
    if (!isTRUE(include_cache)) return(data.frame())
    bundle <- .prep_cache_bundle(prep)
    guard_meta <- if (is.null(bundle) || is.null(bundle$guard_quadrature_meta)) list() else {
      keys <- ls(bundle$guard_quadrature_meta, all.names = TRUE)
      stats::setNames(lapply(keys, function(k) bundle$guard_quadrature_meta[[k]]), keys)
    }
    return(list(nodes = data.frame(), cache = list(
      precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
      pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
      guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE),
      guard_quadrature_orders = guard_meta
    )))
  }
  nodes <- comp$nodes %||% list()
  node_df <- do.call(rbind, lapply(nodes, function(node) {
    data.frame(
      id = node$id,
      kind = node$kind,
      needs_forced = isTRUE(node$needs_forced),
      scenario_sensitive = isTRUE(node$scenario_sensitive),
      sources = paste(node$sources %||% integer(0), collapse = ","),
      args = paste(node$args %||% integer(0), collapse = ","),
      has_fast_density = is.function(node$density_fast_fn),
      has_fast_survival = is.function(node$surv_fast_fn),
      has_fast_cdf = is.function(node$cdf_fast_fn),
      stringsAsFactors = FALSE
    )
  }))
  if (!isTRUE(include_cache)) return(node_df)
  bundle <- .prep_cache_bundle(prep)
  guard_meta <- if (is.null(bundle) || is.null(bundle$guard_quadrature_meta)) list() else {
    keys <- ls(bundle$guard_quadrature_meta, all.names = TRUE)
    stats::setNames(lapply(keys, function(k) bundle$guard_quadrature_meta[[k]]), keys)
  }
  cache_summary <- list(
    precomputed_values = if (is.null(bundle)) character(0) else ls(bundle$precomputed_values, all.names = TRUE),
    pool_templates = if (is.null(bundle)) character(0) else ls(bundle$pool_templates, all.names = TRUE),
    guard_quadrature = if (is.null(bundle)) character(0) else ls(bundle$guard_quadrature, all.names = TRUE),
    guard_quadrature_orders = guard_meta
  )
  list(nodes = node_df, cache = cache_summary)
}

.native_node_batch_eval <- function(prep, node_ids, times, component = NULL,
                                    forced_complete = integer(0),
                                    forced_survive = integer(0)) {
  node_ids <- as.integer(node_ids)
  times <- as.numeric(times)
  if (length(node_ids) == 0L || length(times) == 0L) return(list())
  native_ctx <- .prep_native_context(prep)
  tasks <- lapply(node_ids, function(nid) {
    list(
      node_id = as.integer(nid),
      times = times,
      component = component,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  })
  native_likelihood_eval_cpp(native_ctx, tasks)
}

.state_entry_is_new <- function(entry) {
  FALSE
}

.likelihood_params_by_trial_component <- function(trial_map) {
  if (length(trial_map) == 0) return(list())
  lapply(trial_map, function(df) {
    if (nrow(df) == 0L) return(list())
    comp_vals <- if ("component" %in% names(df)) df$component else rep("__default__", nrow(df))
    if (is.factor(comp_vals)) comp_vals <- as.character(comp_vals)
    comp_vals[is.na(comp_vals)] <- "__default__"
    split(df, comp_vals)
  })
}

.build_native_trial_entries <- function(structure, prep, plan, trial_ids,
                                        data_df, data_row_indices,
                                        component_weights = NULL) {
  trials <- plan$trials %||% list()
  if (length(trials) == 0L) {
    return(NULL)
  }
  metadata <- .likelihood_trial_metadata(structure, prep)
  if (is.null(metadata)) return(NULL)
  outcome_defs <- prep[["outcomes"]] %||% list()
  compiled_nodes <- metadata$compiled_nodes
  compiled_competitors <- metadata$compiled_competitors
  na_source_specs <- metadata$na_source_specs
  guess_target_specs <- metadata$guess_target_specs
  alias_specs <- metadata$alias_specs
  rel_tol <- metadata$rel_tol
  abs_tol <- metadata$abs_tol
  max_depth <- metadata$max_depth

  entries <- vector("list", length(trial_ids))
  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_key <- as.character(tid %||% NA)
    entry <- trials[[trial_key]]
    if (is.null(entry)) return(NULL)
    data_idx <- data_row_indices[[trial_key]]
    if (is.null(data_idx) || length(data_idx) != 1L) return(NULL)
    data_row <- data_df[data_idx, , drop = FALSE]
    outcome_lbl <- data_row$outcome[[1]] %||% NA_character_
    rt_val <- data_row$rt[[1]] %||% NA_real_
    forced_component <- if ("component" %in% names(data_row)) data_row$component[[1]] else NULL
    comp_plan <- .native_component_plan(
      structure = structure,
      trial_rows = entry$rows,
      forced_component = forced_component,
      component_weights = component_weights,
      component_ids = structure$components$component_id
    )
    if (is.na(outcome_lbl) || identical(outcome_lbl, "NA")) {
      guess_key <- "NA"
      guess_donors <- guess_target_specs[[guess_key]] %||% NULL
      entry_na <- list(
        type = "direct",
        trial_rows = entry$rows %||% data.frame(),
        rt = if (is.null(rt_val)) NA_real_ else as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        na_sources = na_source_specs %||% list(),
        component_plan = comp_plan,
        outcome_label = "NA"
      )
      if (!is.null(guess_donors)) entry_na$guess_donors <- guess_donors
      entries[[i]] <- entry_na
      next
    }
    outcome_def <- outcome_defs[[outcome_lbl]]
    if (is.null(outcome_def)) return(NULL)
    alias_sources <- alias_specs[[outcome_lbl]]
    if (!is.null(alias_sources)) {
      entries[[i]] <- list(
        type = "alias_sum",
        trial_rows = entry$rows %||% data.frame(),
        rt = as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        alias_sources = alias_sources,
        component_plan = comp_plan,
        outcome_label = outcome_lbl
      )
      next
    }
    compiled <- compiled_nodes[[outcome_lbl]]
    if (is.null(compiled)) return(NULL)
    comp_ids <- compiled_competitors[[outcome_lbl]]
    if (is.null(comp_ids)) return(NULL)
    if (is.na(rt_val) || !is.finite(rt_val) || rt_val < 0) return(NULL)
    guess_key <- as.character(outcome_lbl)
    guess_donors <- guess_target_specs[[guess_key]] %||% NULL
    entry_fields <- list(
      type = "direct",
      trial_rows = entry$rows %||% data.frame(),
      node_id = as.integer(compiled$id),
      rt = as.numeric(rt_val),
      forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
      competitor_ids = comp_ids,
      component_plan = comp_plan,
      outcome_label = outcome_lbl
    )
    if (!is.null(guess_donors)) entry_fields$guess_donors <- guess_donors
    entries[[i]] <- entry_fields
  }
  component_weights_df <- if (!is.null(component_weights)) as.data.frame(component_weights) else NULL
  list(
    entries = entries,
    component_weights_df = component_weights_df,
    rel_tol = as.numeric(rel_tol),
    abs_tol = as.numeric(abs_tol),
    max_depth = as.integer(max_depth)
  )
}

.native_loglikelihood_batch <- function(structure, prep, plan, trial_ids) {
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) return(NULL)
  cache <- plan$.native_cache %||% NULL
  entries_all <- cache$entries %||% NULL
  entry_names_norm <- cache$entry_names %||% NULL
  eval_cfg <- cache$eval_config %||% NULL
  component_weights_df <- cache$component_weights_df %||% NULL
  if (is.null(entries_all) || is.null(entry_names_norm) || is.null(eval_cfg)) {
    stop("Plan is missing cached native entries; rebuild the likelihood context with data_df", call. = FALSE)
  }
  entry_names_norm[is.na(entry_names_norm)] <- "NA"
  selected_entries <- entries_all
  if (!is.null(trial_ids) && length(trial_ids) > 0L) {
    trial_keys_norm <- as.character(trial_ids %||% NA)
    trial_keys_norm[is.na(trial_keys_norm)] <- "NA"
    idx <- match(trial_keys_norm, entry_names_norm)
    if (any(is.na(idx))) {
      missing_keys <- unique(trial_keys_norm[is.na(idx)])
      stop(sprintf("Missing native entries for: %s", paste(missing_keys, collapse = ", ")), call. = FALSE)
    }
    selected_entries <- entries_all[idx]
    names(selected_entries) <- NULL
  }
  res <- tryCatch(
    native_loglik_from_params_cpp(
      native_ctx,
      structure,
      selected_entries,
      component_weights_df,
      as.numeric(eval_cfg$rel_tol %||% .integrate_rel_tol()),
      as.numeric(eval_cfg$abs_tol %||% .integrate_abs_tol()),
      as.integer(eval_cfg$max_depth %||% 12L)
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  list(
    loglik = as.numeric(res$loglik),
    per_trial = as.numeric(res$per_trial)
  )
}

.likelihood_mixture_likelihood <- function(structure, prep_eval_base, trial_rows,
                                           component_ids, component_weights,
                                           outcome_label, rt_val,
                                           forced_component = NULL,
                                           component_plan_entries = NULL,
                                           use_native_rows = isTRUE(getOption("uuber.use.native.param.rows", TRUE)),
                                           trial_state = NULL) {
  results <- numeric(0)
  plan <- .native_component_plan(
    structure = structure,
    trial_rows = trial_rows,
    forced_component = forced_component,
    component_weights = component_weights,
    component_ids = component_ids
  )
  comps <- plan$components
  weights <- plan$weights
  component_plan_entries <- component_plan_entries %||% list()
  if (use_native_rows) {
    native_mix <- .native_trial_mixture_eval(
      prep = prep_eval_base,
      outcome_label = outcome_label,
      rt_val = rt_val,
      component_plan = plan,
      forced_component = forced_component,
      trial_rows = trial_rows
    )
    if (!is.null(native_mix)) {
      return(as.numeric(native_mix))
    }
  }
  total <- 0.0
  for (idx in seq_along(comps)) {
    comp_id <- comps[[idx]]
    plan_entry <- component_plan_entries[[comp_id]] %||% component_plan_entries[["__default__"]] %||% NULL
    comp_rows_df <- NULL
    if (!is.null(plan_entry)) {
      comp_rows_df <- plan_entry$rows
    }
    if (is.null(comp_rows_df) || nrow(comp_rows_df) == 0L) {
      comp_rows_df <- .likelihood_component_rows(trial_rows, comp_id)
    }
    trial_type_key <- .likelihood_component_label(comp_id)
    cache_key <- .likelihood_outcome_cache_key(trial_type_key, "__structure__", outcome_label, rt_val)
    res <- .likelihood_outcome_cached(prep_eval_base, cache_key, function() {
      .outcome_likelihood(
        outcome_label,
        rt_val,
        prep_eval_base,
        comp_id,
        trial_rows = comp_rows_df,
        state = trial_state
      )
    })
    prep_eval_base <- res$prep
    lik_val <- res$value
    total <- total + weights[[idx]] * as.numeric(lik_val)
  }
  total
}

