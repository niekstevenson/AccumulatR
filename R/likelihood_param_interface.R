.likelihood_component_rows <- function(trial_rows, component) {
  if (is.null(trial_rows) || nrow(trial_rows) == 0L) return(trial_rows)
  if (!"component" %in% names(trial_rows)) return(trial_rows)
  comp_label <- component %||% "__default__"
  mask <- is.na(trial_rows$component)
  if (!identical(comp_label, "__default__")) {
    mask <- mask | trial_rows$component == comp_label
  }
  trial_rows[mask, , drop = FALSE]
}

.likelihood_seed_prep_from_params <- function(prep, params_df) {
  if (is.null(prep) || is.null(params_df) || nrow(params_df) == 0L) {
    return(prep)
  }
  bundle <- .build_trial_param_bundle(prep, params_df)
  acc_overrides <- bundle$acc %||% list()
  if (length(acc_overrides) > 0L) {
    acc_defs <- prep$accumulators %||% list()
    for (acc_id in names(acc_overrides)) {
      override <- acc_overrides[[acc_id]]
      if (is.null(override)) next
      base <- acc_defs[[acc_id]] %||% list(id = acc_id)
      base$params <- override$params
      base$onset <- override$onset
      base$q <- override$q
      base$shared_trigger_id <- override$shared_trigger_id
      acc_defs[[acc_id]] <- base
    }
    prep$accumulators <- acc_defs
  }
  shared_overrides <- bundle$shared %||% list()
  if (length(shared_overrides) > 0L) {
    prep$shared_triggers <- prep$shared_triggers %||% list()
    for (shared_id in names(shared_overrides)) {
      entry <- shared_overrides[[shared_id]]
      prob <- entry$prob
      if (is.null(prob) || length(prob) == 0L) next
      val <- as.numeric(prob[[1]])
      if (!is.finite(val)) next
      if (is.null(prep$shared_triggers[[shared_id]])) {
        prep$shared_triggers[[shared_id]] <- list(
          id = shared_id,
          members = character(0),
          q = val
        )
      } else {
        prep$shared_triggers[[shared_id]]$q <- val
      }
    }
  }
  prep
}

.model_spec_with_params <- function(model_spec, params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) return(model_spec)
  spec_copy <- unserialize(serialize(model_spec, NULL))
  prep_tmp <- .prepare_model_for_likelihood(model_spec)
  bundle <- .build_trial_param_bundle(prep_tmp, params_df)
  acc_overrides <- bundle$acc %||% list()
  if (length(acc_overrides) > 0L) {
    acc_defs <- spec_copy$accumulators %||% list()
    acc_index <- setNames(seq_along(acc_defs), vapply(acc_defs, `[[`, character(1), "id"))
    for (acc_id in names(acc_overrides)) {
      idx <- acc_index[[acc_id]]
      if (is.na(idx)) next
      override <- acc_overrides[[acc_id]]
      acc <- acc_defs[[idx]]
      acc$params <- override$params
      acc$onset <- override$onset
      acc$q <- override$q
      acc_defs[[idx]] <- acc
    }
    spec_copy$accumulators <- acc_defs
  }
  spec_copy
}

.likelihood_build_trial_plan <- function(structure, params_df, prep) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    return(list(order = character(0), trials = list()))
  }
  params_df <- as.data.frame(params_df)
  if (is.null(rownames(params_df))) {
    rownames(params_df) <- as.character(seq_len(nrow(params_df)))
  }
  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  params_df <- .params_hash_attach(params_df)
  structure_hash <- .structure_hash_value(prep)
  trial_rows_map <- .likelihood_params_by_trial(params_df)
  component_ids <- structure$components$component_id
  comp_labels <- component_ids
  if (length(comp_labels) == 0L) comp_labels <- "__default__"
  comp_labels <- unique(c(comp_labels, "__default__"))

  trials <- list()
  trial_order <- character(0)
  for (trial_key in names(trial_rows_map)) {
    trial_df <- as.data.frame(trial_rows_map[[trial_key]])
    trial_id <- trial_df$trial[[1]] %||% NA
    comp_entries <- list()
    for (comp_id in comp_labels) {
      rows_df <- .likelihood_component_rows(trial_df, comp_id)
      hash_val <- NULL
      if (!is.null(rows_df) && nrow(rows_df) > 0L && "params_hash" %in% names(rows_df)) {
        vals <- rows_df$params_hash
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0L) hash_val <- vals[[1]]
      }
      if (is.null(hash_val) || !nzchar(hash_val)) {
        temp_rows <- rows_df
        if (!is.null(temp_rows)) temp_rows$params_hash <- NULL
        hash_val <- .param_rows_hash_value(temp_rows)
      }
      if (is.null(hash_val) || !nzchar(hash_val)) hash_val <- "__base__"
      key <- .likelihood_param_key(structure_hash, comp_id)
      comp_entries[[comp_id]] <- list(rows = rows_df, key = key, hash = hash_val)
    }
    trials[[trial_key]] <- list(
      trial_id = trial_id,
      rows = trial_df,
      components = comp_entries
    )
    trial_order <- c(trial_order, trial_key)
  }
  order_idx <- order(suppressWarnings(as.numeric(trial_order)))
  trial_order <- trial_order[order_idx]
  trials <- trials[trial_order]
  list(order = trial_order, trials = trials)
}

.likelihood_params_by_trial <- function(params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) return(list())
  trials <- params_df$trial %||% seq_len(nrow(params_df))
  if (is.factor(trials)) trials <- as.character(trials)
  split(params_df, as.character(trials))
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

.likelihood_component_weights <- function(structure, trial_id, available_components,
                                          base_weights, component_weights) {
  weights <- base_weights[match(available_components, structure$components$component_id)]
  if (!is.null(component_weights)) {
    rows <- component_weights[
      component_weights$trial == trial_id &
        component_weights$component %in% available_components,
      ,
      drop = FALSE
    ]
    if (nrow(rows) > 0L) {
      weights <- rows$weight
    }
  }
  if (length(weights) != length(available_components) || all(is.na(weights))) {
    weights <- rep(1 / length(available_components), length(available_components))
  } else {
    if (any(!is.finite(weights) | weights < 0)) {
      stop("Component weights must be non-negative finite numbers")
    }
    total <- sum(weights)
    if (!is.finite(total) || total <= 0) {
      weights <- rep(1 / length(available_components), length(available_components))
    } else {
      weights <- weights / total
    }
  }
  weights
}

.likelihood_response_prob_component <- function(prep, outcome_label, component, trial_rows = NULL) {
  use_fastpath <- getOption("uuber.shared_gate_fastpath", default = TRUE)
  if (isTRUE(use_fastpath)) {
    pair <- .find_shared_gate_pair(prep[["outcomes"]])
    if (!is.null(pair) && outcome_label %in% c(pair[['label_x']], pair[['label_y']])) {
      vals <- .shared_gate_pair_probs(prep, component, pair, trial_rows = trial_rows)
      if (identical(outcome_label, pair[['label_x']])) {
        return(vals[[1]])
      } else {
        return(vals[[2]])
      }
    }
  }

  out_def <- prep[["outcomes"]][[outcome_label]]
  if (!is.null(out_def[['options']][['alias_of']])) {
    refs <- as.character(out_def[['options']][['alias_of']])
    vals <- vapply(refs, function(lbl) {
      .outcome_likelihood(
        outcome_label = lbl,
        rt = NA_real_,
        prep = prep,
        component = component,
        trial_rows = trial_rows
      )
    }, numeric(1))
    return(sum(vals))
  }

  base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, component, trial_rows = trial_rows))
  if (!identical(outcome_label, "GUESS")) {
    gp <- .get_component_attr(prep, component, "guess")
    if (!is.null(gp) && !is.null(gp[['weights']])) {
      keep <- gp[['weights']][[outcome_label]] %||% gp[['weights']][[normalize_label(outcome_label)]] %||% 1.0
      base <- base * as.numeric(keep)
    }
  }
  base
}

.component_plan_fallback <- function(structure, component_ids, trial_rows,
                                     forced_component = NULL, component_weights = NULL) {
  comps <- component_ids
  if (!is.null(forced_component) && !is.na(forced_component)) {
    comps <- intersect(component_ids, as.character(forced_component))
  } else if (!is.null(trial_rows) && "component" %in% names(trial_rows)) {
    listed <- unique(trial_rows$component)
    listed <- listed[!is.na(listed)]
    if (length(listed) > 0L) {
      comps <- intersect(component_ids, listed)
    }
  }
  if (length(comps) == 0L) comps <- component_ids
  trial_id <- if (!is.null(trial_rows) && "trial" %in% names(trial_rows) && nrow(trial_rows) > 0L) {
    trial_rows$trial[[1]]
  } else NA_integer_
  weights <- .likelihood_component_weights(
    structure,
    trial_id,
    comps,
    structure$components$weight,
    component_weights
  )
  list(components = comps, weights = weights)
}

.native_component_plan <- function(structure, trial_rows,
                                   forced_component = NULL,
                                   component_weights = NULL,
                                   component_ids = NULL) {
  native_fn <- .lik_native_fn("native_component_plan_exported")
  forced_chr <- forced_component
  if (length(forced_chr) == 0L || is.null(forced_chr) || is.na(forced_chr)) {
    forced_chr <- NULL
  } else {
    forced_chr <- as.character(forced_chr)[[1]]
  }
  trial_df <- if (is.null(trial_rows)) {
    data.frame()
  } else {
    trial_rows
  }
  component_weights_df <- component_weights %||% NULL
  plan <- tryCatch(
    native_fn(
      structure,
      trial_df,
      forced_chr,
      component_weights_df
    ),
    error = function(e) NULL
  )
  if (is.null(plan)) {
    return(.component_plan_fallback(structure, component_ids, trial_rows, forced_component, component_weights))
  }
  comps <- as.character(plan$components %||% character(0))
  weights <- as.numeric(plan$weights %||% numeric(0))
  if (length(comps) == 0L || length(weights) != length(comps)) {
    return(.component_plan_fallback(structure, component_ids, trial_rows, forced_component, component_weights))
  }
  list(components = comps, weights = weights)
}

.likelihood_component_label <- function(component) {
  if (is.null(component) || length(component) == 0L) {
    "__default__"
  } else {
    as.character(component)[[1]]
  }
}

.likelihood_param_key <- function(structure_hash, component) {
  struct_chr <- structure_hash %||% "__structure__"
  comp_chr <- .likelihood_component_label(component)
  paste(struct_chr, comp_chr, sep = "|")
}
.native_trial_mixture_eval <- function(prep, outcome_label, rt_val, component_plan,
                                       forced_component = NULL,
                                       trial_rows = NULL) {
  if (is.null(prep) || is.null(outcome_label) || length(outcome_label) == 0L) return(NULL)
  outcome_chr <- as.character(outcome_label)[[1]]
  if (is.na(outcome_chr) || is.na(rt_val) || !is.finite(rt_val) || rt_val < 0) return(NULL)
  components <- component_plan$components %||% character(0)
  if (length(components) == 0L) return(NULL)
  outcome_defs <- prep[["outcomes"]] %||% list()
  outcome_def <- outcome_defs[[outcome_chr]]
  if (is.null(outcome_def)) return(NULL)
  expr <- outcome_def[["expr"]]
  if (is.null(expr)) return(NULL)
  compiled <- .expr_lookup_compiled(expr, prep)
  if (is.null(compiled)) return(NULL)
  competitor_map <- .prep_competitors(prep) %||% list()
  competitor_exprs <- competitor_map[[outcome_chr]] %||% list()
  if (length(competitor_exprs) > 0L) {
    comp_nodes <- lapply(competitor_exprs, function(ex) .expr_lookup_compiled(ex, prep))
    if (any(vapply(comp_nodes, is.null, logical(1)))) return(NULL)
    comp_ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
    if (any(is.na(comp_ids))) return(NULL)
  } else {
    comp_ids <- integer(0)
  }
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) return(NULL)
  native_fn <- .lik_native_fn("native_trial_mixture_cpp")
  forced_chr <- NULL
  if (!is.null(forced_component) && !is.na(forced_component)) {
    forced_chr <- as.character(forced_component)[[1]]
  }
  trial_df <- if (is.null(trial_rows)) data.frame() else trial_rows
  res <- tryCatch(
    native_fn(
      native_ctx,
      as.integer(compiled$id),
      as.numeric(rt_val),
      as.character(components),
      as.numeric(component_plan$weights %||% numeric(0)),
      forced_chr,
      as.integer(comp_ids),
      trial_df
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(res) || res < 0) return(NULL)
  res
}

.build_native_trial_entries <- function(structure, prep, plan, trial_ids,
                                        data_df, data_row_indices,
                                        component_weights = NULL) {
  trials <- plan$trials %||% list()
  if (length(trials) == 0L) {
    return(NULL)
  }
  structure_hash <- .structure_hash_value(prep)
  outcome_defs <- prep[["outcomes"]] %||% list()
  competitor_map <- .prep_competitors(prep) %||% list()
  compiled_nodes <- lapply(outcome_defs, function(def) {
    expr <- def[["expr"]]
    if (is.null(expr)) return(NULL)
    .expr_lookup_compiled(expr, prep)
  })
  names(compiled_nodes) <- names(outcome_defs)
  compiled_competitors <- lapply(names(outcome_defs), function(lbl) {
    comp_exprs <- competitor_map[[lbl]] %||% list()
    if (length(comp_exprs) == 0) return(integer(0))
    comp_nodes <- lapply(comp_exprs, function(ex) .expr_lookup_compiled(ex, prep))
    if (any(vapply(comp_nodes, is.null, logical(1)))) return(NULL)
    ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
    if (any(is.na(ids))) return(NULL)
    ids
  })
  names(compiled_competitors) <- names(outcome_defs)
  shared_gate_specs <- lapply(names(outcome_defs), function(lbl) {
    info <- .detect_shared_gate(outcome_defs, lbl)
    if (is.null(info)) return(NULL)
    x_lbl <- info$x_id %||% NULL
    y_lbl <- info$y_id %||% NULL
    c_lbl <- info$c_id %||% NULL
    if (is.null(x_lbl) || is.null(y_lbl) || is.null(c_lbl)) return(NULL)
    list(
      x_label = as.character(x_lbl),
      y_label = as.character(y_lbl),
      c_label = as.character(c_lbl)
    )
  })
  names(shared_gate_specs) <- names(outcome_defs)
  na_source_labels <- Filter(function(lbl) {
    def <- outcome_defs[[lbl]]
    map_to <- def[['options']][['map_outcome_to']]
    if (is.null(map_to)) return(FALSE)
    if (is.na(map_to)) return(TRUE)
    identical(map_to, "NA")
  }, names(outcome_defs))
  na_source_specs <- NULL
  if (length(na_source_labels) > 0L) {
    na_source_specs <- lapply(na_source_labels, function(lbl) {
      node <- compiled_nodes[[lbl]]
      comp_ids <- compiled_competitors[[lbl]]
      if (is.null(node) || is.null(comp_ids)) return(NULL)
      list(
        node_id = as.integer(node$id),
        competitor_ids = comp_ids %||% integer(0)
      )
    })
    if (any(vapply(na_source_specs, is.null, logical(1)))) na_source_specs <- NULL
  }
  default_deadline <- prep$default_deadline %||% Inf
  rel_tol <- .integrate_rel_tol()
  abs_tol <- .integrate_abs_tol()
  max_depth <- 12L
  entries <- vector("list", length(trial_ids))
  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_key <- as.character(tid %||% NA)
    entry <- trials[[trial_key]]
    if (is.null(entry)) return(NULL)
    raw_components <- entry$components %||% list()
    component_cache <- lapply(names(raw_components), function(lbl) {
      comp_entry <- raw_components[[lbl]] %||% list()
      key_val <- comp_entry$key %||% .likelihood_param_key(structure_hash, lbl)
      hash_val <- comp_entry$hash %||% "__base__"
      if (is.null(hash_val) || !nzchar(hash_val)) hash_val <- "__base__"
      list(
        component = lbl,
        key = key_val,
        hash = hash_val
      )
    })
    names(component_cache) <- names(raw_components)
    data_idx <- data_row_indices[[trial_key]]
    if (is.null(data_idx) || length(data_idx) != 1L) return(NULL)
    data_row <- data_df[data_idx, , drop = FALSE]
    outcome_lbl <- data_row$outcome[[1]] %||% NA_character_
    rt_val <- data_row$rt[[1]] %||% NA_real_
    forced_component <- if ("component" %in% names(data_row)) data_row$component[[1]] else NULL
    if (is.na(outcome_lbl) || identical(outcome_lbl, "NA")) {
      if (is.null(na_source_specs)) return(NULL)
      entries[[i]] <- list(
        type = "na_map",
        trial_rows = entry$rows %||% data.frame(),
        rt = if (is.null(rt_val)) NA_real_ else as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        na_sources = na_source_specs,
        component_cache = component_cache
      )
      next
    }
    outcome_def <- outcome_defs[[outcome_lbl]]
    if (is.null(outcome_def)) return(NULL)
    alias_refs <- outcome_def[['options']][['alias_of']]
    if (!is.null(alias_refs)) {
      refs <- as.character(alias_refs)
      alias_sources <- lapply(refs, function(ref_lbl) {
        node <- compiled_nodes[[ref_lbl]]
        comp_ids <- compiled_competitors[[ref_lbl]]
        if (is.null(node) || is.null(comp_ids)) return(NULL)
        list(
          node_id = as.integer(node$id),
          competitor_ids = comp_ids %||% integer(0)
        )
      })
      if (any(vapply(alias_sources, is.null, logical(1)))) return(NULL)
      entries[[i]] <- list(
        type = "alias_sum",
        trial_rows = entry$rows %||% data.frame(),
        rt = as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        alias_sources = alias_sources,
        component_cache = component_cache
      )
      next
    }
    compiled <- compiled_nodes[[outcome_lbl]]
    if (is.null(compiled)) return(NULL)
    comp_ids <- compiled_competitors[[outcome_lbl]]
    if (is.null(comp_ids)) return(NULL)
    if (is.na(rt_val) || !is.finite(rt_val) || rt_val < 0) return(NULL)
    entry_fields <- list(
      type = "direct",
      trial_rows = entry$rows %||% data.frame(),
      node_id = as.integer(compiled$id),
      rt = as.numeric(rt_val),
      forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
      competitor_ids = comp_ids,
      component_cache = component_cache
    )
    shared_info <- shared_gate_specs[[outcome_lbl]]
    if (!is.null(shared_info)) {
      entry_fields$shared_gate <- shared_info
    }
    entries[[i]] <- entry_fields
  }
  component_weights_df <- if (!is.null(component_weights)) as.data.frame(component_weights) else NULL
  list(
    entries = entries,
    component_weights_df = component_weights_df,
    default_deadline = as.numeric(default_deadline),
    rel_tol = as.numeric(rel_tol),
    abs_tol = as.numeric(abs_tol),
    max_depth = as.integer(max_depth)
  )
}

.native_loglikelihood_batch <- function(structure, prep, plan, trial_ids,
                                        data_df, data_row_indices,
                                        component_weights = NULL) {
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) return(NULL)
  builder <- .build_native_trial_entries(
    structure = structure,
    prep = prep,
    plan = plan,
    trial_ids = trial_ids,
    data_df = data_df,
    data_row_indices = data_row_indices,
    component_weights = component_weights
  )
  if (is.null(builder)) return(NULL)
  entries <- builder$entries
  component_weights_df <- builder$component_weights_df
  default_deadline <- builder$default_deadline
  rel_tol <- builder$rel_tol
  abs_tol <- builder$abs_tol
  max_depth <- builder$max_depth
  native_fn <- .lik_native_fn("native_loglik_from_params_cpp")
  res <- native_fn(
    native_ctx,
    structure,
    entries,
    component_weights_df,
    as.numeric(default_deadline),
    as.numeric(rel_tol),
    as.numeric(abs_tol),
    as.integer(max_depth)
  )
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
                                           use_native_rows = isTRUE(getOption("uuber.use.native.param.rows", TRUE))) {
  results <- numeric(0)
  structure_hash <- .structure_hash_value(prep_eval_base)
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
    trial_type_key <- if (!is.null(plan_entry)) plan_entry$key else .likelihood_param_key(structure_hash, comp_id)
    params_hash <- NULL
    if (!is.null(plan_entry)) params_hash <- plan_entry$hash
    if (is.null(params_hash) || !nzchar(params_hash)) {
      if (!is.null(comp_rows_df) && nrow(comp_rows_df) > 0L && "params_hash" %in% names(comp_rows_df)) {
        vals <- comp_rows_df$params_hash
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0L) params_hash <- vals[[1]]
      }
    }
    if (is.null(params_hash) || !nzchar(params_hash)) params_hash <- "__base__"
    cache_key <- .likelihood_outcome_cache_key(trial_type_key, params_hash, outcome_label, rt_val)
    res <- .likelihood_outcome_cached(prep_eval_base, cache_key, function() {
      .outcome_likelihood(
        outcome_label,
        rt_val,
        prep_eval_base,
        comp_id,
        trial_rows = comp_rows_df
      )
    })
    prep_eval_base <- res$prep
    lik_val <- res$value
    total <- total + weights[[idx]] * as.numeric(lik_val)
  }
  total
}

log_likelihood_from_params <- function(structure, params_df, data_df,
                                       component_weights = NULL,
                                       prep = NULL,
                                       trial_plan = NULL,
                                       native_bundle = NULL) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  params_df <- .params_hash_attach(params_df)
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain outcome/rt per trial")
  }
  structure <- .as_generator_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
  prep_eval_base <- .likelihood_seed_prep_from_params(prep_eval_base, params_df)
  if (is.null(prep_eval_base[[".runtime"]])) {
    prep_eval_base <- .prep_set_cache_bundle(
      prep_eval_base,
      .build_likelihood_cache_bundle(prep_eval_base)
    )
  }
  comp_ids <- structure$components$component_id

  if (!"trial" %in% names(data_df)) {
    stop("Data frame must include a 'trial' column")
  }
  data_df$trial <- data_df$trial
  if ("component" %in% names(data_df)) {
    data_df$component <- as.character(data_df$component)
  }
  trial_lookup_keys <- as.character(data_df$trial %||% NA)
  data_row_indices <- split(seq_len(nrow(data_df)), trial_lookup_keys)

  plan <- trial_plan %||% .likelihood_build_trial_plan(structure,
                                                      params_df,
                                                      prep_eval_base)
  trial_ids <- unique(data_df$trial)
  per_trial_loglik <- numeric(length(trial_ids))
  use_native_rows <- isTRUE(getOption("uuber.use.native.param.rows", TRUE))
  batch_res <- NULL
  if (use_native_rows) {
    batch_res <- .native_loglikelihood_batch(
      structure = structure,
      prep = prep_eval_base,
      plan = plan,
      trial_ids = trial_ids,
      data_df = data_df,
      data_row_indices = data_row_indices,
      component_weights = component_weights
    )
  }
  if (!is.null(batch_res)) {
    return(list(
      loglik = batch_res$loglik,
      per_trial = batch_res$per_trial
    ))
  }

  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_key <- as.character(tid %||% NA)
    entry <- plan$trials[[trial_key]]
    if (is.null(entry)) {
      fallback_rows <- params_df[params_df$trial == tid, , drop = FALSE]
      entry <- list(
        trial_id = tid,
        rows = fallback_rows,
        components = list()
      )
    }
    trial_rows <- entry$rows
    data_idx <- data_row_indices[[trial_key]]
    if (is.null(data_idx) || length(data_idx) != 1L) {
      stop(sprintf("Expected exactly one data row for trial %s", tid))
    }
    data_row <- data_df[data_idx, , drop = FALSE]
    outcome <- data_row$outcome[[1]] %||% NA_character_
    rt_val <- data_row$rt[[1]] %||% NA_real_
    forced_component <- if ("component" %in% names(data_row)) data_row$component[[1]] else NULL
    component_entries <- entry$components %||% list()

    mixture <- .likelihood_mixture_likelihood(
      structure = structure,
      prep_eval_base = prep_eval_base,
      trial_rows = trial_rows,
      component_ids = comp_ids,
      component_weights = component_weights,
      outcome_label = outcome,
      rt_val = rt_val,
      forced_component = forced_component,
      component_plan_entries = component_entries,
      use_native_rows = use_native_rows
    )
    if (!is.finite(mixture) || mixture <= 0) {
      per_trial_loglik[[i]] <- -Inf
    } else {
      per_trial_loglik[[i]] <- log(mixture)
    }
  }
  list(
    loglik = sum(per_trial_loglik),
    per_trial = per_trial_loglik
  )
}

log_likelihood_from_params_buffer <- function(structure, params_df, data_df,
                                              component_weights = NULL,
                                              prep = NULL,
                                              trial_plan = NULL,
                                              native_bundle = NULL) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain outcome/rt per trial")
  }
  structure <- .as_generator_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
  prep_eval_base <- .likelihood_seed_prep_from_params(prep_eval_base, params_df)
  if (is.null(prep_eval_base[[".runtime"]])) {
    prep_eval_base <- .prep_set_cache_bundle(
      prep_eval_base,
      .build_likelihood_cache_bundle(prep_eval_base)
    )
  }
  if (!"trial" %in% names(data_df)) {
    stop("Data frame must include a 'trial' column")
  }
  data_df$trial <- data_df$trial
  if ("component" %in% names(data_df)) {
    data_df$component <- as.character(data_df$component)
  }
  trial_lookup_keys <- as.character(data_df$trial %||% NA)
  data_row_indices <- split(seq_len(nrow(data_df)), trial_lookup_keys)

  plan <- trial_plan %||% .likelihood_build_trial_plan(
    structure,
    params_df,
    prep_eval_base
  )
  trial_ids <- unique(data_df$trial)
  batch_res <- .native_loglikelihood_batch(
    structure = structure,
    prep = prep_eval_base,
    plan = plan,
    trial_ids = trial_ids,
    data_df = data_df,
    data_row_indices = data_row_indices,
    component_weights = component_weights
  )
  if (is.null(batch_res)) {
    stop("native log-likelihood batch evaluation failed; rebuild the trial plan or inspect params_df")
  }
  list(
    loglik = batch_res$loglik,
    per_trial = batch_res$per_trial
  )
}

build_buffer_trial_entries <- function(structure, params_df, data_df,
                                       component_weights = NULL,
                                       prep = NULL,
                                       trial_plan = NULL,
                                       native_bundle = NULL) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain outcome/rt per trial")
  }
  structure <- .as_generator_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
  prep_eval_base <- .likelihood_seed_prep_from_params(prep_eval_base, params_df)
  if (is.null(prep_eval_base[[".runtime"]])) {
    prep_eval_base <- .prep_set_cache_bundle(
      prep_eval_base,
      .build_likelihood_cache_bundle(prep_eval_base)
    )
  }
  if (!"trial" %in% names(data_df)) {
    stop("Data frame must include a 'trial' column")
  }
  data_df$trial <- data_df$trial
  if ("component" %in% names(data_df)) {
    data_df$component <- as.character(data_df$component)
  }
  trial_lookup_keys <- as.character(data_df$trial %||% NA)
  data_row_indices <- split(seq_len(nrow(data_df)), trial_lookup_keys)

  plan <- trial_plan %||% .likelihood_build_trial_plan(
    structure,
    params_df,
    prep_eval_base
  )
  trial_ids <- unique(data_df$trial)
  builder <- .build_native_trial_entries(
    structure = structure,
    prep = prep_eval_base,
    plan = plan,
    trial_ids = trial_ids,
    data_df = data_df,
    data_row_indices = data_row_indices,
    component_weights = component_weights
  )
  list(
    entries = builder$entries,
    default_deadline = builder$default_deadline,
    rel_tol = builder$rel_tol,
    abs_tol = builder$abs_tol,
    max_depth = builder$max_depth,
    component_weights = builder$component_weights_df,
    prep = prep_eval_base,
    trial_plan = plan,
    data_row_indices = data_row_indices
  )
}

.aggregate_observed_probs <- function(prep, probs, include_na = TRUE) {
  labels <- names(prep[["outcomes"]])
  labels <- Filter(function(lbl) {
    is.null(prep[["outcomes"]][[lbl]][['options']][['alias_of']])
  }, labels)
  obs <- numeric(0)
  na_sum <- 0.0
  for (lbl in labels) {
    out_def <- prep[["outcomes"]][[lbl]]
    map_to <- out_def[['options']][['map_outcome_to']] %||% NULL
    prob <- probs[[lbl]] %||% 0.0
    if (is.null(map_to)) {
      current <- obs[lbl]
      if (length(current) == 0L || is.na(current)) current <- 0.0
      obs[lbl] <- current + prob
    } else {
      map_label <- as.character(map_to)[1]
      if (!nzchar(map_label) || identical(map_label, "NA") || is.na(map_label)) {
        na_sum <- na_sum + prob
      } else {
        current <- obs[map_label]
        if (length(current) == 0L || is.na(current)) current <- 0.0
        obs[map_label] <- current + prob
      }
    }
  }
  if (include_na && (na_sum > 0 || length(obs) == 0L)) {
    obs <- c(obs, "NA" = na_sum)
  }
  if (include_na) {
    total <- sum(obs)
    resid <- 1.0 - total
    if (!is.finite(resid)) resid <- 0.0
    if (resid > .Machine$double.eps) {
      if ("NA" %in% names(obs)) {
        obs[["NA"]] <- obs[["NA"]] + resid
      } else {
        obs <- c(obs, "NA" = resid)
      }
    }
  }
  obs
}

observed_response_probabilities_from_params <- function(structure, params_df,
                                                        component_weights = NULL,
                                                        include_na = TRUE) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  params_df <- .params_hash_attach(params_df)
  structure <- .as_generator_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  model_spec <- .model_spec_with_params(structure$model_spec, params_df)
  prep_eval_base <- .prepare_model_for_likelihood(model_spec)
  comp_ids <- structure$components$component_id

  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }

  trial_ids <- sort(unique(params_df$trial))
  accum <- list()
  weights_per_trial <- numeric(length(trial_ids))

  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_rows <- params_df[params_df$trial == tid, , drop = FALSE]
    comps <- comp_ids
    if ("component" %in% names(trial_rows)) {
      listed <- unique(trial_rows$component)
      listed <- listed[!is.na(listed)]
      if (length(listed) > 0L) comps <- intersect(comp_ids, listed)
    }
    if (length(comps) == 0L) comps <- comp_ids
    weights <- .likelihood_component_weights(
      structure,
      tid,
      comps,
      structure$components$weight,
      component_weights
    )
    trial_probs <- numeric(0)
    for (idx in seq_along(comps)) {
      comp_id <- comps[[idx]]
      comp_rows <- .likelihood_component_rows(trial_rows, comp_id)
      labels <- names(prep_eval_base$outcomes)
      base_probs <- setNames(vapply(labels, function(lbl) {
        .likelihood_response_prob_component(
          prep_eval_base,
          lbl,
          comp_id,
          trial_rows = comp_rows
        )
      }, numeric(1)), labels)
      weighted <- weights[[idx]] * base_probs
      if (length(trial_probs) == 0L) {
        trial_probs <- weighted
      } else {
        for (lbl in names(weighted)) {
          current <- trial_probs[[lbl]]
          if (is.null(current) || is.na(current)) current <- 0
          trial_probs[[lbl]] <- current + weighted[[lbl]]
        }
      }
    }
    agg <- .aggregate_observed_probs(prep_eval_base, trial_probs, include_na = include_na)
    accum[[i]] <- agg
    weights_per_trial[[i]] <- 1.0
  }
  total_weight <- sum(weights_per_trial)
  combined_labels <- Reduce(union, lapply(accum, names))
  result <- setNames(numeric(length(combined_labels)), combined_labels)
  for (i in seq_along(accum)) {
    vec <- accum[[i]]
    w <- weights_per_trial[[i]] / total_weight
    for (lbl in names(vec)) {
      result[[lbl]] <- result[[lbl]] + w * vec[[lbl]]
    }
  }
  result
}
