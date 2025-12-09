.likelihood_context_structure <- function(structure, params_df, data_df,
                                          component_weights = NULL,
                                          prep = NULL,
                                          trial_plan = NULL,
                                          native_bundle = NULL) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row", call. = FALSE)
  }
  if (is.null(data_df) || nrow(data_df) == 0L) {
    stop("Data frame must contain outcome/rt per trial", call. = FALSE)
  }
  structure <- .as_generator_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("generator structure must include model_spec; rebuild with build_generator_structure")
  }
  data_df <- as.data.frame(data_df)
  required_cols <- c("trial", "outcome", "rt")
  missing_cols <- setdiff(required_cols, names(data_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Data frame must include columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  data_df$trial <- data_df$trial
  data_df$outcome <- as.character(data_df$outcome)
  data_df$rt <- as.numeric(data_df$rt)
  if ("component" %in% names(data_df)) {
    data_df$component <- as.character(data_df$component)
  }
  trial_lookup_keys <- as.character(data_df$trial %||% NA)
  data_row_indices <- split(seq_len(nrow(data_df)), trial_lookup_keys)
  trial_ids <- unique(data_df$trial)
  if (!is.null(native_bundle)) {
    prep <- native_bundle$prep %||% prep
  }
  prep_eval_base <- prep %||% .prepare_model_for_likelihood(structure$model_spec)
  if (is.null(prep_eval_base[[".runtime"]]) || is.null(prep_eval_base$.runtime$cache_bundle)) {
    prep_eval_base <- .prep_set_cache_bundle(
      prep_eval_base,
      .build_likelihood_cache_bundle(prep_eval_base)
    )
  }
  plan <- trial_plan %||% .likelihood_build_trial_plan(structure, params_df, prep_eval_base)
  if (is.null(plan$.native_cache)) {
    plan$.native_cache <- new.env(parent = emptyenv())
  }
  native_ctx <- .prep_native_context(prep_eval_base)
  if (!inherits(native_ctx, "externalptr")) {
    native_ctx <- NULL
  }
  plan_cache <- plan$.native_cache
  plan_cache$entries <- NULL
  plan_cache$entry_names <- NULL
  plan_cache$eval_config <- NULL
  plan_cache$component_weights_df <- NULL
  component_weights_df <- if (!is.null(component_weights)) as.data.frame(component_weights) else NULL
  if (!is.null(native_ctx)) {
    metadata <- .likelihood_trial_metadata(structure, prep_eval_base)
    if (!is.null(metadata)) {
      plan_keys_raw <- as.character((plan$order %||% names(plan$trials)) %||% character(0))
      if (length(plan_keys_raw) == 0L) {
        stop("Plan contains no trials", call. = FALSE)
      }
      data_trial_keys <- as.character(data_df$trial %||% NA)
      entries_all <- tryCatch(
        native_plan_entries_cpp(
          native_ctx,
          structure,
          plan,
          plan_keys_raw,
          data_trial_keys,
          data_df,
          component_weights_df,
          metadata$na_source_specs %||% list(),
          metadata$guess_target_specs %||% list(),
          metadata$alias_specs %||% list()
        ),
        error = function(e) NULL
      )
      if (!is.null(entries_all)) {
        entry_names <- names(entries_all)
        if (is.null(entry_names)) {
          entry_names <- plan_keys_raw
        }
        entry_names_norm <- entry_names
        entry_names_norm[is.na(entry_names_norm)] <- "NA"
        plan_cache$entries <- entries_all
        plan_cache$entry_names <- entry_names_norm
        plan_cache$eval_config <- list(
          rel_tol = as.numeric(metadata$rel_tol),
          abs_tol = as.numeric(metadata$abs_tol),
          max_depth = as.integer(metadata$max_depth)
        )
        plan_cache$component_weights_df <- component_weights_df
      }
    }
  }
  structure(list(
    structure = structure,
    params_df = params_df,
    prep = prep_eval_base,
    plan = plan,
    component_weights = component_weights,
    native_ctx = native_ctx,
    data_df = data_df,
    data_row_indices = data_row_indices,
    trial_ids = trial_ids
  ), class = "likelihood_context")
}

build_likelihood_context <- function(structure, params_df, data_df,
                                     component_weights = NULL,
                                     prep = NULL,
                                     trial_plan = NULL,
                                     native_bundle = NULL) {
  .likelihood_context_structure(
    structure = structure,
    params_df = params_df,
    data_df = data_df,
    component_weights = component_weights,
    prep = prep,
    trial_plan = trial_plan,
    native_bundle = native_bundle
  )
}

.validate_likelihood_context <- function(context) {
  if (inherits(context, "likelihood_context")) return(context)
  stop("likelihood context must be created via build_likelihood_context()", call. = FALSE)
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

.model_spec_with_params <- function(model_spec, params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) return(model_spec)
  spec_copy <- unserialize(serialize(model_spec, NULL))
  prep_tmp <- .prepare_model_for_likelihood(model_spec)
  param_state <- .generator_param_state_from_rows(prep_tmp, params_df)
  acc_overrides <- param_state$accs %||% list()
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
  shared_overrides <- param_state$shared %||% list()
  if (length(shared_overrides) > 0L && length(spec_copy$groups) > 0L) {
    for (g in seq_along(spec_copy$groups)) {
      grp <- spec_copy$groups[[g]] %||% list()
      trig <- grp$attrs$shared_trigger %||% NULL
      if (is.null(trig)) next
      trig_id <- trig$id %||% grp$id
      if (is.null(trig_id) || !nzchar(trig_id)) next
      override <- shared_overrides[[trig_id]] %||% NULL
      if (is.null(override) || is.null(override$prob)) next
      grp$attrs$shared_trigger$q <- override$prob
      spec_copy$groups[[g]] <- grp
    }
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
      comp_entries[[comp_id]] <- list(rows = rows_df)
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
  list(
    order = trial_order,
    trials = trials,
    .native_cache = new.env(parent = emptyenv())
  )
}

.likelihood_trial_metadata <- function(structure, prep) {
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
      def <- outcome_defs[[lbl]] %||% list()
      comp_ids <- compiled_competitors[[lbl]]
      if (is.null(node) || is.null(comp_ids)) return(NULL)
      list(
        node_id = as.integer(node$id),
        competitor_ids = comp_ids %||% integer(0),
        source_label = lbl
      )
    })
    if (any(vapply(na_source_specs, is.null, logical(1)))) na_source_specs <- NULL
  }
  guess_target_specs <- list()
  for (lbl in names(outcome_defs)) {
    def <- outcome_defs[[lbl]]
    opts <- def[['options']] %||% list()
    guess_opts <- opts[['guess']]
    if (is.null(guess_opts)) next
    donor_node <- compiled_nodes[[lbl]]
    if (is.null(donor_node)) next
    donor_comp <- compiled_competitors[[lbl]] %||% integer(0)
    labels <- guess_opts[['labels']] %||% character(0)
    weights <- guess_opts[['weights']] %||% numeric(0)
    if (!length(labels) || length(labels) != length(weights)) next
    rt_policy <- guess_opts[['rt_policy']] %||% "keep"
    for (j in seq_along(labels)) {
      tgt <- labels[[j]]
      tgt_key <- if (is.na(tgt)) "NA" else as.character(tgt)
      release_val <- as.numeric(weights[[j]])
      release_val <- if (is.finite(release_val)) release_val else 0
      if (release_val < 0) release_val <- 0
      if (release_val > 1) release_val <- 1
      donor_rec <- list(
        node_id = as.integer(donor_node$id),
        competitor_ids = donor_comp %||% integer(0),
        source_label = lbl,
        release = release_val,
        rt_policy = rt_policy
      )
      guess_target_specs[[tgt_key]] <- c(guess_target_specs[[tgt_key]] %||% list(), list(donor_rec))
    }
  }
  alias_specs <- lapply(names(outcome_defs), function(lbl) {
    def <- outcome_defs[[lbl]]
    alias_refs <- def[['options']][['alias_of']]
    if (is.null(alias_refs)) return(NULL)
    refs <- as.character(alias_refs)
    alias_sources <- lapply(refs, function(ref_lbl) {
      node <- compiled_nodes[[ref_lbl]]
      comp_ids <- compiled_competitors[[ref_lbl]]
      if (is.null(node) || is.null(comp_ids)) return(NULL)
      list(
        node_id = as.integer(node$id),
        competitor_ids = comp_ids %||% integer(0),
        source_label = ref_lbl
      )
    })
    if (any(vapply(alias_sources, is.null, logical(1)))) return(NULL)
    alias_sources
  })
  names(alias_specs) <- names(outcome_defs)
  rel_tol <- .integrate_rel_tol()
  abs_tol <- .integrate_abs_tol()
  max_depth <- 12L
  list(
    structure_hash = structure_hash,
    compiled_nodes = compiled_nodes,
    compiled_competitors = compiled_competitors,
    na_source_specs = na_source_specs,
    guess_target_specs = guess_target_specs,
    alias_specs = alias_specs,
    rel_tol = as.numeric(rel_tol),
    abs_tol = as.numeric(abs_tol),
    max_depth = as.integer(max_depth)
  )
}

.params_df_to_matrix <- function(df) {
  if (!"trial" %in% names(df)) {
    stop("Parameter table must include 'trial'", call. = FALSE)
  }

  n <- nrow(df)
  if (is.null(n) || n == 0L) {
    stop("Parameter table must contain rows", call. = FALSE)
  }

  # Core columns
  trial <- as.numeric(df$trial)
  acc_idx <- if ("accumulator_index" %in% names(df)) {
    as.numeric(df$accumulator_index)
  } else if ("accumulator" %in% names(df)) {
    as.numeric(df$accumulator)
  } else {
    stop("Parameter table must include 'accumulator_index' (or numeric 'accumulator')", call. = FALSE)
  }
  onset <- as.numeric(df$onset %||% 0)
  q <- as.numeric(df$q %||% 0)
  t0 <- as.numeric(df$t0 %||% 0)
  shared_id <- if ("shared_trigger_id" %in% names(df)) as.character(df$shared_trigger_id) else rep(NA_character_, n)
  shared_q <- if ("shared_trigger_q" %in% names(df)) as.numeric(df$shared_trigger_q) else rep(NA_real_, n)

  # Helper: per-row first-non-NA across candidate columns
  pick_slot <- function(candidates) {
    out <- rep(NA_real_, n)
    for (nm in candidates) {
      if (is.null(nm) || !nm %in% names(df)) next
      col <- as.numeric(df[[nm]])
      update_mask <- is.na(out) & !is.na(col)
      if (any(update_mask)) {
        out[update_mask] <- col[update_mask]
      }
      if (all(!is.na(out))) break
    }
    out
  }

  # Slot priorities (t0 is always the first parameter column)
  p1 <- pick_slot(c("p1", "meanlog", "shape", "mu"))
  p2 <- pick_slot(c("p2", "sdlog", "rate", "sigma"))
  p3 <- pick_slot(c("p3", "tau"))

  max_params <- 0L
  if (!all(is.na(p1))) max_params <- max(max_params, 1L)
  if (!all(is.na(p2))) max_params <- max(max_params, 2L)
  if (!all(is.na(p3))) max_params <- max(max_params, 3L)

  cols <- list(
    trial = trial,
    accumulator_index = acc_idx,
    onset = onset,
    q = q,
    t0 = t0,
    shared_trigger_q = shared_q
  )
  if (max_params >= 1L) cols$p1 <- p1
  if (max_params >= 2L) cols$p2 <- p2
  if (max_params >= 3L) cols$p3 <- p3

  mat <- do.call(cbind, cols)
  storage.mode(mat) <- "double"
  attr(mat, "trial") <- trial
  colnames(mat) <- names(cols)
  mat
}

.likelihood_params_by_trial <- function(params_df) {
  if (is.null(params_df) || nrow(params_df) == 0L) return(list())
  trials <- params_df$trial %||% seq_len(nrow(params_df))
  if (is.factor(trials)) trials <- as.character(trials)
  split(params_df, as.character(trials))
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

.likelihood_response_prob_component <- function(prep, outcome_label, component,
                                                trial_rows = NULL,
                                                trial_state = NULL) {
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

  # Prefer native probability (upper = Inf) so shared triggers and guards follow the
  # same path as the main likelihood. Fall back to the pure-R integrator if needed.
  base <- NA_real_
  native_ctx <- .prep_native_context(prep)
  if (inherits(native_ctx, "externalptr")) {
    compiled <- .expr_lookup_compiled(out_def[['expr']], prep)
    comp_exprs <- (.prep_competitors(prep) %||% list())[[outcome_label]] %||% list()
    comp_ids <- integer(0)
    if (length(comp_exprs) > 0L) {
      comp_nodes <- lapply(comp_exprs, function(ex) .expr_lookup_compiled(ex, prep))
      if (!any(vapply(comp_nodes, is.null, logical(1)))) {
        comp_ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
      } else {
        comp_ids <- integer(0)
      }
    }
    trial_rows_df <- if (is.null(trial_rows)) data.frame() else as.data.frame(trial_rows)
    if (!is.null(compiled) && length(comp_ids) > 0L) {
      prob_native <- tryCatch(
        native_outcome_probability_params_cpp(
          native_ctx,
          as.integer(compiled$id),
          Inf,
          component,
          integer(0),
          integer(0),
          as.integer(comp_ids),
          .integrate_rel_tol(),
          .integrate_abs_tol(),
          12L,
          trial_rows_df
        ),
        error = function(e) NA_real_
      )
      if (is.finite(prob_native) && prob_native >= 0) {
        base <- as.numeric(prob_native)
      }
    }
  }
  if (is.na(base)) {
    base <- as.numeric(.outcome_likelihood(
      outcome_label,
      NA_real_,
      prep,
      component,
      trial_rows = trial_rows,
      state = trial_state
    ))
  }
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
    native_component_plan_exported(
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
.native_trial_mixture_eval <- function(prep, outcome_label, rt_val, component_plan,
                                       forced_component = NULL,
                                       trial_rows = NULL,
                                       guess_donors = NULL) {
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
  forced_chr <- NULL
  if (!is.null(forced_component) && !is.na(forced_component)) {
    forced_chr <- as.character(forced_component)[[1]]
  }
  trial_df <- if (is.null(trial_rows)) data.frame() else trial_rows
  res <- tryCatch(
    native_trial_mixture_cpp(
      native_ctx,
      as.integer(compiled$id),
      as.numeric(rt_val),
      as.character(components),
      as.numeric(component_plan$weights %||% numeric(0)),
      forced_chr,
      as.integer(comp_ids),
      trial_df,
      guess_donors %||% list()
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(res) || res < 0) return(NULL)
  res
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
    trial_state <- .eval_state_create()
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
          trial_rows = comp_rows,
          trial_state = trial_state
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


native_loglikelihood_param_repeat <- function(likelihood_context, params_list) {
  ctx <- .validate_likelihood_context(likelihood_context)
  cache <- ctx$plan$.native_cache %||% NULL
  entries <- cache$entries %||% NULL
  eval_cfg <- cache$eval_config %||% NULL
  native_ctx <- ctx$native_ctx
  if (is.null(entries) || is.null(eval_cfg) || !inherits(native_ctx, "externalptr")) {
    # Fallback: rebuild context from stored data/params in the likelihood_context
    if (!is.null(ctx$data_df) && !is.null(ctx$params_df)) {
      rebuilt <- build_likelihood_context(
        structure = ctx$structure,
        params_df = ctx$params_df,
        data_df = ctx$data_df,
        component_weights = ctx$component_weights,
        prep = ctx$prep
      )
      ctx <- rebuilt
      cache <- ctx$plan$.native_cache %||% cache
      entries <- cache$entries %||% entries
      eval_cfg <- cache$eval_config %||% eval_cfg
      native_ctx <- ctx$native_ctx
    }
    if (is.null(entries) || is.null(eval_cfg) || !inherits(native_ctx, "externalptr")) {
      stop("Cached native entries are required; build the context with data_df first", call. = FALSE)
    }
  }
  params_list <- lapply(params_list, function(df) .params_df_to_matrix(as.data.frame(df)))
  native_loglik_param_repeat_cpp(
    native_ctx,
    ctx$structure,
    entries,
    cache$component_weights_df %||% NULL,
    params_list,
    as.numeric(eval_cfg$rel_tol %||% .integrate_rel_tol()),
    as.numeric(eval_cfg$abs_tol %||% .integrate_abs_tol()),
    as.integer(eval_cfg$max_depth %||% 12L)
  )
}
