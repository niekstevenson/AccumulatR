.likelihood_apply_overrides <- function(prep, acc_overrides, shared_overrides) {
  trial_prep <- prep
  modified <- FALSE
  if (!is.null(acc_overrides) && length(acc_overrides) > 0) {
    for (acc_id in names(acc_overrides)) {
      override <- acc_overrides[[acc_id]]
      if (is.null(override)) next
      trial_prep$accumulators[[acc_id]] <- override
      modified <- TRUE
    }
  }
  if (!is.null(shared_overrides) && length(shared_overrides) > 0) {
    if (is.null(trial_prep$shared_triggers)) trial_prep$shared_triggers <- list()
    for (shared_id in names(shared_overrides)) {
      entry <- shared_overrides[[shared_id]]
      if (is.null(entry)) next
      prob <- entry$prob
      if (length(prob) == 0L || is.null(prob) || !is.finite(prob)) next
      if (is.null(trial_prep$shared_triggers[[shared_id]])) {
        trial_prep$shared_triggers[[shared_id]] <- list(
          id = shared_id,
          members = character(0),
          q = as.numeric(prob)
        )
      } else {
        trial_prep$shared_triggers[[shared_id]]$q <- as.numeric(prob)
      }
      modified <- TRUE
    }
  }

  runtime <- trial_prep[[".runtime"]] %||% list()
  base_bundle <- runtime$cache_bundle
  if (modified) {
    if (!is.null(base_bundle)) {
      runtime$cache_bundle <- .likelihood_cache_bundle_clone(base_bundle)
    } else {
      runtime$cache_bundle <- .build_likelihood_cache_bundle(trial_prep)
    }
    trial_prep[[".runtime"]] <- runtime
  }

  trial_prep
}

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

.likelihood_build_trial_plan <- function(structure, params_df, prep,
                                         build_override_ptr = TRUE,
                                         keep_trial_rows = TRUE,
                                         keep_component_rows = TRUE) {
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
  trial_rows_map <- .likelihood_params_by_trial(params_df)
  native_ctx <- .prep_native_context(prep)
  builder <- .lik_native_fn("native_build_trial_overrides_cpp")
  component_ids <- structure$components$component_id
  comp_labels <- component_ids
  if (length(comp_labels) == 0L) comp_labels <- "__default__"
  comp_labels <- unique(c(comp_labels, "__default__"))

  trials <- list()
  trial_order <- character(0)
  for (trial_key in names(trial_rows_map)) {
    trial_df <- as.data.frame(trial_rows_map[[trial_key]])
    trial_id <- trial_df$trial[[1]] %||% NA
    override_ptr <- NULL
    if (build_override_ptr && nrow(trial_df) > 0L) {
      override_ptr <- builder(native_ctx, trial_df)
    }
    trial_row_index <- suppressWarnings(as.integer(rownames(trial_df)))
    if (anyNA(trial_row_index)) {
      trial_row_index <- seq_len(nrow(trial_df)) - 1L
    } else {
      trial_row_index <- trial_row_index - 1L
    }
    comp_entries <- list()
    for (comp_id in comp_labels) {
      rows_df <- NULL
      if (keep_component_rows) {
        rows_df <- .likelihood_component_rows(trial_df, comp_id)
      }
      key <- if (!is.null(rows_df) && nrow(rows_df) > 0L) {
        .likelihood_override_key(comp_id, rows_df)
      } else {
        .likelihood_base_override_key(comp_id)
      }
      comp_row_index <- if (!is.null(rows_df) && nrow(rows_df) > 0L) {
        idx <- suppressWarnings(as.integer(rownames(rows_df)))
        if (anyNA(idx)) {
          seq_len(nrow(rows_df)) - 1L
        } else {
          idx - 1L
        }
      } else {
        integer(0)
      }
      comp_entries[[comp_id]] <- list(rows = if (keep_component_rows) rows_df else NULL,
                                      key = key,
                                      row_index = comp_row_index)
    }
    trials[[trial_key]] <- list(
      trial_id = trial_id,
      rows = if (keep_trial_rows) trial_df else NULL,
      row_index = trial_row_index,
      components = comp_entries,
      override_ptr = override_ptr
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

.likelihood_response_prob_component <- function(prep, outcome_label, component) {
  use_fastpath <- getOption("uuber.shared_gate_fastpath", default = TRUE)
  if (isTRUE(use_fastpath)) {
    pair <- .find_shared_gate_pair(prep[["outcomes"]])
    if (!is.null(pair) && outcome_label %in% c(pair[['label_x']], pair[['label_y']])) {
      vals <- .shared_gate_pair_probs(prep, component, pair)
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
    vals <- vapply(refs, function(lbl) as.numeric(.outcome_likelihood(lbl, NA_real_, prep, component)), numeric(1))
    return(sum(vals))
  }

  base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, component))
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

.likelihood_base_override_key <- function(component) {
  paste(.likelihood_component_label(component), "__base__", sep = "|")
}

.likelihood_canonicalise_params <- function(entry) {
  if (is.null(entry) || !is.list(entry)) return(entry)
  if (is.null(names(entry))) return(entry)
  entry[order(names(entry))]
}

.likelihood_override_key <- function(component, rows) {
  base_key <- .likelihood_base_override_key(component)
  if (is.null(rows) || nrow(rows) == 0L) return(base_key)

  key_rows <- rows
  drop_cols <- c(
    "trial", "component", "accumulator_index", "acc_idx",
    "type", "role", "outcome", "rt", "condition",
    "component_weight"
  )
  keep_cols <- setdiff(names(key_rows), drop_cols)
  if (length(keep_cols) == 0L) return(base_key)
  key_rows <- key_rows[, keep_cols, drop = FALSE]

  for (col in names(key_rows)) {
    if (is.factor(key_rows[[col]])) {
      key_rows[[col]] <- as.character(key_rows[[col]])
    }
  }

  if ("params" %in% names(key_rows)) {
    key_rows$params <- lapply(key_rows$params, .likelihood_canonicalise_params)
  }

  order_cols <- intersect(
    c("accumulator_id", "accumulator", "accumulator_index", "acc_idx"),
    names(key_rows)
  )
  if (length(order_cols) > 0L) {
    ord_args <- lapply(order_cols, function(col) key_rows[[col]])
    ord <- do.call(order, ord_args)
    key_rows <- key_rows[ord, , drop = FALSE]
  }
  rownames(key_rows) <- NULL

  raw_bytes <- serialize(key_rows, connection = NULL, ascii = TRUE)
  raw_hex <- paste(as.character(raw_bytes), collapse = "")
  paste(.likelihood_component_label(component), raw_hex, sep = "|")
}

.likelihood_component_override_bundle <- function(structure, trial_rows, component, component_rows = NULL) {
  base_key <- .likelihood_base_override_key(component)
  if (is.null(trial_rows) || nrow(trial_rows) == 0L) {
    return(list(acc = list(), shared = list(), key = base_key))
  }

  effective_rows <- component_rows
  if (is.null(effective_rows)) {
    effective_rows <- trial_rows
    if ("component" %in% names(effective_rows)) {
      effective_rows <- effective_rows[
        is.na(effective_rows$component) | effective_rows$component == component,
        ,
        drop = FALSE
      ]
    }
  }
  if (nrow(effective_rows) == 0L) {
    return(list(acc = list(), shared = list(), key = base_key))
  }

  overrides <- .build_trial_overrides(structure, effective_rows)
  acc_overrides <- overrides$acc %||% list()
  if (length(acc_overrides) > 0L) {
    acc_overrides <- Filter(function(x) !is.null(x), acc_overrides)
  }
  shared_overrides <- overrides$shared %||% list()
  if (length(shared_overrides) > 0L) {
    shared_overrides <- Filter(function(x) !is.null(x), shared_overrides)
  }
  if (length(acc_overrides) == 0L && length(shared_overrides) == 0L) {
    return(list(acc = list(), shared = list(), key = base_key))
  }

  key <- .likelihood_override_key(component, effective_rows)
  list(acc = acc_overrides, shared = shared_overrides, key = key)
}

.likelihood_fetch_component_prep <- function(structure, prep_base, trial_rows, component,
                                             cache_env = NULL, override_bundle = NULL) {
  override_bundle <- override_bundle %||% .likelihood_component_override_bundle(structure, trial_rows, component)
  key <- override_bundle$key
  base_key <- .likelihood_base_override_key(component)

  if (!is.null(cache_env) && is.environment(cache_env)) {
    cached <- cache_env[[key]]
    if (!is.null(cached)) return(cached)
  }

  if (identical(key, base_key)) {
    prep_val <- prep_base
  } else {
    prep_val <- .likelihood_apply_overrides(prep_base, override_bundle$acc, override_bundle$shared)
  }

  if (!is.null(cache_env) && is.environment(cache_env) && !identical(key, base_key)) {
    cache_env[[key]] <- prep_val
  }
  prep_val
}

.likelihood_outcome_cache_key <- function(bundle_key, outcome_label, rt_val) {
  outcome_chr <- if (is.null(outcome_label) || length(outcome_label) == 0L || is.na(outcome_label[[1]])) {
    "NA"
  } else {
    as.character(outcome_label)[[1]]
  }
  rt_key <- .eval_state_time_key(rt_val)
  paste(bundle_key, outcome_chr, rt_key, sep = "|")
}

.native_trial_mixture_eval <- function(prep, outcome_label, rt_val, component_plan,
                                       forced_component = NULL,
                                       trial_override_ptr = NULL) {
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
  res <- tryCatch(
    native_fn(
      native_ctx,
      as.integer(compiled$id),
      as.numeric(rt_val),
      as.character(components),
      as.numeric(component_plan$weights %||% numeric(0)),
      forced_chr,
      as.integer(comp_ids),
      trial_override_ptr %||% NULL
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(res) || res < 0) return(NULL)
  res
}

.native_loglikelihood_batch <- function(structure, prep, plan, trial_ids,
                                        data_df, data_row_indices,
                                        component_weights = NULL,
                                        params_df = NULL,
                                        use_buffer = FALSE) {
  if (use_buffer && is.null(params_df)) {
    stop("params_df must be supplied when use_buffer = TRUE")
  }
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) return(NULL)
  trials <- plan$trials %||% list()
  if (length(trials) == 0L) return(NULL)
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
    if (is.null(entry)) {
      if (use_buffer) stop(sprintf("Buffer plan missing trial entry for key '%s'", trial_key))
      return(NULL)
    }
    row_index <- entry$row_index %||% integer(0)
    if (use_buffer && length(row_index) == 0L && is.null(entry$override_ptr)) {
      stop("Buffer trial plan is missing row indices; rebuild with keep_trial_rows = FALSE / keep_component_rows = FALSE")
    }
    data_idx <- data_row_indices[[trial_key]]
    if (is.null(data_idx) || length(data_idx) != 1L) {
      if (use_buffer) stop(sprintf("Expected exactly one data row for trial %s", trial_key))
      return(NULL)
    }
    data_row <- data_df[data_idx, , drop = FALSE]
    outcome_lbl <- data_row$outcome[[1]] %||% NA_character_
    rt_val <- data_row$rt[[1]] %||% NA_real_
    forced_component <- if ("component" %in% names(data_row)) data_row$component[[1]] else NULL
    if (is.na(outcome_lbl) || identical(outcome_lbl, "NA")) {
      if (is.null(na_source_specs)) {
        if (use_buffer) stop("NA map encountered but na_source_specs unavailable")
        return(NULL)
      }
      entries[[i]] <- list(
        type = "na_map",
        trial_rows = entry$rows %||% data.frame(),
        rt = if (is.null(rt_val)) NA_real_ else as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        na_sources = na_source_specs,
        override_ptr = entry$override_ptr,
        row_index = row_index
      )
      next
    }
    outcome_def <- outcome_defs[[outcome_lbl]]
    if (is.null(outcome_def)) {
      if (use_buffer) stop(sprintf("Outcome '%s' not found in prep", outcome_lbl))
      return(NULL)
    }
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
      if (any(vapply(alias_sources, is.null, logical(1)))) {
        if (use_buffer) stop("Alias sources missing compiled nodes")
        return(NULL)
      }
      entries[[i]] <- list(
        type = "alias_sum",
        trial_rows = entry$rows %||% data.frame(),
        rt = as.numeric(rt_val),
        forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
        alias_sources = alias_sources,
        override_ptr = entry$override_ptr,
        row_index = row_index
      )
      next
    }
    compiled <- compiled_nodes[[outcome_lbl]]
    if (is.null(compiled)) {
      if (use_buffer) stop(sprintf("Outcome '%s' missing compiled node", outcome_lbl))
      return(NULL)
    }
    comp_ids <- compiled_competitors[[outcome_lbl]]
    if (is.null(comp_ids)) {
      if (use_buffer) stop(sprintf("Outcome '%s' missing competitor ids", outcome_lbl))
      return(NULL)
    }
    if (is.na(rt_val) || !is.finite(rt_val) || rt_val < 0) {
      if (use_buffer) stop("Encountered invalid RT in buffer path")
      return(NULL)
    }
    entry_fields <- list(
      type = "direct",
      trial_rows = entry$rows %||% data.frame(),
      node_id = as.integer(compiled$id),
      rt = as.numeric(rt_val),
      forced_component = if (!is.null(forced_component) && !is.na(forced_component)) as.character(forced_component) else NULL,
      competitor_ids = comp_ids,
      override_ptr = entry$override_ptr,
      row_index = row_index
    )
    shared_info <- shared_gate_specs[[outcome_lbl]]
    if (!is.null(shared_info)) {
      entry_fields$shared_gate <- shared_info
    }
    entries[[i]] <- entry_fields
  }
  component_weights_df <- if (!is.null(component_weights)) as.data.frame(component_weights) else NULL
  if (use_buffer) {
    if (is.null(params_df)) {
      stop("params_df must be supplied when use_buffer = TRUE")
    }
    params_df_export <- as.data.frame(params_df)
    native_fn <- .lik_native_fn("native_loglik_from_buffer_cpp")
    res <- native_fn(
      native_ctx,
      structure,
      entries,
      params_df_export,
      component_weights_df,
      as.numeric(default_deadline),
      as.numeric(rel_tol),
      as.numeric(abs_tol),
      as.integer(max_depth)
    )
  } else {
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
  }
  list(
    loglik = as.numeric(res$loglik),
    per_trial = as.numeric(res$per_trial)
  )
}

.likelihood_mixture_likelihood <- function(structure, prep_eval_base, trial_rows,
                                           component_ids, component_weights,
                                           outcome_label, rt_val,
                                           forced_component = NULL,
                                           prep_cache = NULL,
                                           component_rows_map = NULL,
                                           component_plan_entries = NULL,
                                           trial_override_ptr = NULL,
                                           use_native_rows = isTRUE(getOption("uuber.use.native.param.rows", TRUE))) {
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
      trial_override_ptr = trial_override_ptr
    )
    if (!is.null(native_mix)) {
      return(as.numeric(native_mix))
    }
  }
  total <- 0.0
  for (idx in seq_along(comps)) {
    comp_id <- comps[[idx]]
    comp_rows <- NULL
    if (!is.null(component_rows_map) && length(component_rows_map) > 0L) {
      comp_rows <- component_rows_map[[comp_id]] %||% component_rows_map[["__default__"]] %||% NULL
    }
    comp_rows_df <- if (!is.null(comp_rows)) as.data.frame(comp_rows) else NULL
    plan_entry <- component_plan_entries[[comp_id]] %||% component_plan_entries[["__default__"]] %||% NULL
    if (!is.null(plan_entry)) {
      comp_rows_df <- plan_entry$rows
      if (is.null(comp_rows)) {
        comp_rows <- plan_entry$rows
      }
    }
    use_trial_rows <- use_native_rows &&
      !is.null(comp_rows_df) &&
      nrow(comp_rows_df) > 0L

    if (use_trial_rows) {
      override_key <- if (!is.null(plan_entry)) plan_entry$key else .likelihood_override_key(comp_id, comp_rows_df)
      cache_key <- .likelihood_outcome_cache_key(override_key, outcome_label, rt_val)
      res <- .likelihood_outcome_cached(prep_eval_base, cache_key, function() {
        .outcome_likelihood(
          outcome_label,
          rt_val,
          prep_eval_base,
          comp_id,
          trial_rows = comp_rows_df,
          trial_overrides = trial_override_ptr
        )
      })
      prep_eval_base <- res$prep
      lik_val <- res$value
    } else {
      override_bundle <- .likelihood_component_override_bundle(structure, trial_rows, comp_id, component_rows = comp_rows)
      trial_prep <- .likelihood_fetch_component_prep(
        structure,
        prep_eval_base,
        trial_rows,
        comp_id,
        cache_env = prep_cache,
        override_bundle = override_bundle
      )
      cache_key <- .likelihood_outcome_cache_key(override_bundle$key, outcome_label, rt_val)
      res <- .likelihood_outcome_cached(trial_prep, cache_key, function() {
        .outcome_likelihood(outcome_label, rt_val, trial_prep, comp_id)
      })
      trial_prep <- res$prep
      if (identical(override_bundle$key, .likelihood_base_override_key(comp_id))) {
        prep_eval_base <- trial_prep
      }
      lik_val <- res$value
    }
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
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain outcome/rt per trial")
  }
  structure <- .as_generator_structure(structure)
  if (!exists(".build_trial_overrides", mode = "function")) {
    stop("likelihood param interface requires generator_new.R to be sourced")
  }
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
  if (is.null(prep_eval_base[[".runtime"]])) {
    prep_eval_base <- .prep_set_cache_bundle(
      prep_eval_base,
      .build_likelihood_cache_bundle(prep_eval_base)
    )
  }
  comp_ids <- structure$components$component_id

  prep_cache <- new.env(parent = emptyenv(), hash = TRUE)

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
                                                      prep_eval_base,
                                                      build_override_ptr = TRUE,
                                                      keep_trial_rows = TRUE,
                                                      keep_component_rows = TRUE)
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
        components = list(),
        override_ptr = NULL
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
    overrides_ptr <- entry$override_ptr

    mixture <- .likelihood_mixture_likelihood(
      structure = structure,
      prep_eval_base = prep_eval_base,
      trial_rows = trial_rows,
      component_ids = comp_ids,
      component_weights = component_weights,
      outcome_label = outcome,
      rt_val = rt_val,
      forced_component = forced_component,
      prep_cache = prep_cache,
      component_rows_map = NULL,
      component_plan_entries = component_entries,
      trial_override_ptr = overrides_ptr,
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
  if (!exists(".build_trial_overrides", mode = "function")) {
    stop("likelihood param interface requires generator_new.R to be sourced")
  }
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
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
    prep_eval_base,
    build_override_ptr = FALSE,
    keep_trial_rows = FALSE,
    keep_component_rows = FALSE
  )
  trial_ids <- unique(data_df$trial)
  batch_res <- .native_loglikelihood_batch(
    structure = structure,
    prep = prep_eval_base,
    plan = plan,
    trial_ids = trial_ids,
    data_df = data_df,
    data_row_indices = data_row_indices,
    component_weights = component_weights,
    params_df = params_df,
    use_buffer = TRUE
  )
  if (is.null(batch_res)) {
    stop("native_loglik_from_buffer_cpp returned NULL; rebuild the trial plan or inspect params_df")
  }
  list(
    loglik = batch_res$loglik,
    per_trial = batch_res$per_trial
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
    } else if (is.na(map_to)) {
      na_sum <- na_sum + prob
    } else {
      obs_lbl <- as.character(map_to)
      current <- obs[obs_lbl]
      if (length(current) == 0L || is.na(current)) current <- 0.0
      obs[obs_lbl] <- current + prob
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
  structure <- .as_generator_structure(structure)
  if (!exists(".build_trial_overrides", mode = "function")) {
    stop("likelihood param interface requires generator_new.R to be sourced")
  }
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  prep_eval_base <- .prepare_model_for_likelihood(structure$model_spec)
  comp_ids <- structure$components$component_id
  prep_cache <- new.env(parent = emptyenv(), hash = TRUE)

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
      trial_prep <- .likelihood_fetch_component_prep(structure, prep_eval_base, trial_rows, comp_id, cache_env = prep_cache)
      labels <- names(trial_prep$outcomes)
      base_probs <- setNames(vapply(labels, function(lbl) {
        .likelihood_response_prob_component(trial_prep, lbl, comp_id)
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
