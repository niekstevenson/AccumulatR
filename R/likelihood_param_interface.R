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
  structure <- .as_model_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("model structure must include model_spec; rebuild with finalize_model")
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
  component_weights_df <- .likelihood_component_weights_df(
    structure = structure,
    params_df = params_df,
    component_weights = component_weights
  )
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

#' Build a likelihood context. Needed for efficient C++ evaluation
#'
#' @param structure Generator structure
#' @param params_df Parameter data frame
#' @param data_df Data frame with observed outcome/rt
#' @param component_weights Optional component weights
#' @param prep Optional prep bundle
#' @param trial_plan Optional trial plan
#' @param native_bundle Optional native bundle
#' @return likelihood_context object
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_params_df(
#'   spec,
#'   c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' build_likelihood_context(structure, params_df, data_df)
#' @export
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
                                          base_weights, component_weights,
                                          trial_rows = NULL) {
  mode <- structure$components$mode[[1]] %||% "fixed"
  comp_table <- structure$components
  if (!identical(mode, "sample")) {
    weights <- base_weights[match(available_components, comp_table$component_id)]
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
    return(weights)
  }
  ref_id <- comp_table$reference[[1]] %||% available_components[[1]]
  comp_attrs <- comp_table$attrs
  weights <- numeric(length(available_components))
  names(weights) <- available_components
  sum_nonref <- 0
  for (i in seq_along(available_components)) {
    cid <- available_components[[i]]
    if (identical(cid, ref_id)) next
    attrs <- comp_attrs[[cid]] %||% list()
    wp <- attrs$weight_param
    if (is.null(wp) || !nzchar(wp)) {
      stop(sprintf("Component '%s' missing weight_param for sampled mixture", cid))
    }
    val <- NA_real_
    if (!is.null(trial_rows) && wp %in% names(trial_rows)) {
      vals <- trial_rows[[wp]]
      if (!is.null(vals) && length(vals) > 0) val <- as.numeric(vals[[1]])
    }
    if (is.na(val) && !is.null(component_weights) && wp %in% names(component_weights)) {
      rows <- component_weights[
        component_weights$trial == trial_id &
          component_weights$component %in% available_components,
        ,
        drop = FALSE
      ]
      if (nrow(rows) > 0L && wp %in% names(rows)) {
        val <- as.numeric(rows[[wp]][[1]])
      }
    }
    if (is.na(val)) {
      idx <- which(comp_table$component_id == cid)
      if (length(idx) == 1) val <- comp_table$weight[[idx]]
    }
    if (is.na(val) || !is.finite(val) || val < 0 || val > 1) {
      stop(sprintf("Mixture weight '%s' must be a probability in [0,1]", wp))
    }
    weights[[i]] <- val
    sum_nonref <- sum_nonref + val
  }
  if (ref_id %in% available_components) {
    ref_idx <- match(ref_id, available_components)
    ref_weight <- 1 - sum_nonref
    if (ref_weight < -1e-8) stop("Mixture weights sum to >1; check non-reference weight params")
    if (ref_weight < 0) ref_weight <- 0
    weights[[ref_idx]] <- ref_weight
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    weights <- rep(1 / length(available_components), length(available_components))
  } else {
    weights <- weights / total
  }
  weights
}

.likelihood_component_weights_df <- function(structure, params_df, component_weights = NULL) {
  if (!is.null(component_weights)) {
    return(as.data.frame(component_weights))
  }
  comp_table <- structure$components %||% NULL
  if (is.null(comp_table) || nrow(comp_table) == 0L) return(NULL)
  mode <- comp_table$mode[[1]] %||% "fixed"
  if (!identical(mode, "sample")) return(NULL)
  params_df <- as.data.frame(params_df)
  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  trial_rows_map <- .likelihood_params_by_trial(params_df)
  comp_ids <- comp_table$component_id
  weight_rows <- list()
  for (trial_key in names(trial_rows_map)) {
    trial_df <- as.data.frame(trial_rows_map[[trial_key]])
    tid <- trial_df$trial[[1]] %||% NA
    comps <- comp_ids
    if ("component" %in% names(trial_df)) {
      listed <- unique(trial_df$component)
      listed <- listed[!is.na(listed)]
      if (length(listed) > 0L) {
        comps <- intersect(comps, listed)
      }
    }
    if (length(comps) == 0L) comps <- comp_ids
    weights <- .likelihood_component_weights(
      structure,
      tid,
      comps,
      comp_table$weight,
      NULL,
      trial_rows = trial_df
    )
    weight_rows[[length(weight_rows) + 1L]] <- data.frame(
      trial = tid,
      component = comps,
      weight = weights,
      stringsAsFactors = FALSE
    )
  }
  if (length(weight_rows) == 0L) return(NULL)
  do.call(rbind, weight_rows)
}

.likelihood_response_prob_component <- function(prep, outcome_label, component,
                                                trial_rows = NULL,
                                                trial_state = NULL) {
  out_def <- prep[["outcomes"]][[outcome_label]]
  if (!is.null(out_def[['options']][['alias_of']])) {
    refs <- as.character(out_def[['options']][['alias_of']])
    vals <- vapply(refs, function(lbl) {
      .likelihood_response_prob_component(
        prep = prep,
        outcome_label = lbl,
        component = component,
        trial_rows = trial_rows,
        trial_state = trial_state
      )
    }, numeric(1))
    return(sum(vals, na.rm = TRUE))
  }

  # Native-only probability (upper = Inf) so shared triggers and guards follow the
  # same path as the main likelihood.
  native_ctx <- .prep_native_context(prep)
  if (!inherits(native_ctx, "externalptr")) {
    stop("Native context required for response probabilities; R fallback removed", call. = FALSE)
  }

  compiled <- .expr_lookup_compiled(out_def[['expr']], prep)
  if (is.null(compiled)) {
    stop(sprintf("No compiled node for outcome '%s'", outcome_label), call. = FALSE)
  }

  comp_exprs <- (.prep_competitors(prep) %||% list())[[outcome_label]] %||% list()
  comp_ids <- integer(0)
  if (length(comp_exprs) > 0L) {
    comp_nodes <- lapply(comp_exprs, function(ex) .expr_lookup_compiled(ex, prep))
    if (!any(vapply(comp_nodes, is.null, logical(1)))) {
      comp_ids <- vapply(comp_nodes, function(node) as.integer(node$id %||% NA_integer_), integer(1))
    }
  }

  trial_rows_df <- if (is.null(trial_rows)) data.frame() else as.data.frame(trial_rows)
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
  if (!is.finite(prob_native) || prob_native < 0) {
    stop(sprintf("Native outcome probability failed for outcome '%s'", outcome_label), call. = FALSE)
  }
  base <- as.numeric(prob_native)

  if (!identical(outcome_label, "GUESS")) {
    gp <- .get_component_attr(prep, component, "guess")
    if (!is.null(gp) && !is.null(gp[['weights']])) {
      keep <- gp[['weights']][[outcome_label]] %||% gp[['weights']][[normalize_label(outcome_label)]] %||% 1.0
      base <- base * as.numeric(keep)
    }
  }
  base
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

#' @rdname response_probabilities
#' @export
response_probabilities.default <- function(structure, ...) {
  stop("response_probabilities() expects a model_structure", call. = FALSE)
}

#' Analytic outcome probabilities for a parameter set
#'
#' @param structure Model structure
#' @param params_df Parameter data frame (one or more trials)
#' @param component_weights Optional component weights
#' @param include_na Whether to include NA outcome mass
#' @param ... Unused; for S3 compatibility
#' @return Named numeric vector of probabilities
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_params_df(
#'   spec,
#'   c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 1
#' )
#' response_probabilities(structure, params_df)
#' @export
response_probabilities <- function(structure,
                                   params_df,
                                   component_weights = NULL,
                                   include_na = TRUE,
                                   ...) {
  UseMethod("response_probabilities")
}

#' @rdname response_probabilities
#' @export
response_probabilities.model_structure <- function(structure,
                                                   params_df,
                                                   component_weights = NULL,
                                                   include_na = TRUE,
                                                   ...) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  structure <- .as_model_structure(structure)
  if (is.null(structure$model_spec)) {
    stop("model structure must include model_spec; rebuild with finalize_model")
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
      component_weights,
      trial_rows = trial_rows
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


#' Native log-likelihood evaluation over multiple parameter sets
#'
#' @param likelihood_context Context built by build_likelihood_context
#' @param parameters Parameter data frame, or list of parameter data frames
#' @param ... Unused; for S3 compatibility
#' @return Numeric vector of log-likelihood values
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params_df <- build_params_df(
#'   spec,
#'   c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
#'   n_trials = 2
#' )
#' data_df <- simulate(structure, params_df, seed = 1)
#' ctx <- build_likelihood_context(structure, params_df, data_df)
#' log_likelihood(ctx, list(params_df))
#' @export
log_likelihood <- function(likelihood_context, parameters, ...) {
  UseMethod("log_likelihood")
}

#' @rdname log_likelihood
#' @export
log_likelihood.likelihood_context <- function(likelihood_context, parameters, ...) {
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
  if (is.data.frame(parameters) || is.matrix(parameters)) {
    parameters <- list(as.data.frame(parameters))
  } else {
    parameters <- lapply(parameters, as.data.frame)
  }

  comp_mode <- ctx$structure$components$mode[[1]] %||% "fixed"
  if (identical(comp_mode, "sample")) {
    if (is.null(native_ctx) || !inherits(native_ctx, "externalptr")) {
      stop("Native context required for log_likelihood", call. = FALSE)
    }
    plan_keys_raw <- as.character((ctx$plan$order %||% names(ctx$plan$trials)) %||% character(0))
    if (length(plan_keys_raw) == 0L) {
      stop("Plan contains no trials", call. = FALSE)
    }
    data_trial_keys <- as.character(ctx$data_df$trial %||% NA)
    metadata <- .likelihood_trial_metadata(ctx$structure, ctx$prep)
    rel_tol <- as.numeric(metadata$rel_tol %||% .integrate_rel_tol())
    abs_tol <- as.numeric(metadata$abs_tol %||% .integrate_abs_tol())
    max_depth <- as.integer(metadata$max_depth %||% 12L)
    na_specs <- metadata$na_source_specs %||% list()
    guess_specs <- metadata$guess_target_specs %||% list()
    alias_specs <- metadata$alias_specs %||% list()
    out <- numeric(length(parameters))
    for (i in seq_along(parameters)) {
      df <- parameters[[i]]
      comp_weights_df <- .likelihood_component_weights_df(
        structure = ctx$structure,
        params_df = df,
        component_weights = ctx$component_weights
      )
      params_mat <- .params_df_to_matrix(df)
      entries_i <- native_plan_entries_cpp(
        native_ctx,
        ctx$structure,
        ctx$plan,
        plan_keys_raw,
        data_trial_keys,
        ctx$data_df,
        comp_weights_df,
        na_specs,
        guess_specs,
        alias_specs
      )
      res <- cpp_loglik_multiple(
        native_ctx,
        ctx$structure,
        entries_i,
        comp_weights_df,
        list(params_mat),
        rel_tol,
        abs_tol,
        max_depth
      )
      out[[i]] <- as.numeric(res[[1]])
    }
    return(out)
  }

  params_list <- lapply(parameters, function(df) .params_df_to_matrix(as.data.frame(df)))
  cpp_loglik_multiple(
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

#' @rdname log_likelihood
#' @export
log_likelihood.default <- function(likelihood_context, ...) {
  stop("log_likelihood() expects a likelihood_context", call. = FALSE)
}
