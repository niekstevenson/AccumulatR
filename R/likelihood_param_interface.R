.likelihood_apply_overrides <- function(prep, acc_overrides, shared_overrides) {
  trial_prep <- prep
  if (!is.null(acc_overrides) && length(acc_overrides) > 0) {
    for (acc_id in names(acc_overrides)) {
      override <- acc_overrides[[acc_id]]
      if (is.null(override)) next
      trial_prep$accumulators[[acc_id]] <- override
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
    }
  }
  trial_prep
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

.likelihood_component_override_bundle <- function(structure, trial_rows, component) {
  base_key <- .likelihood_base_override_key(component)
  if (is.null(trial_rows) || nrow(trial_rows) == 0L) {
    return(list(acc = list(), shared = list(), key = base_key))
  }

  effective_rows <- trial_rows
  if ("component" %in% names(effective_rows)) {
    effective_rows <- effective_rows[
      is.na(effective_rows$component) | effective_rows$component == component,
      ,
      drop = FALSE
    ]
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
                                             cache_env = NULL, bundle = NULL) {
  bundle <- bundle %||% .likelihood_component_override_bundle(structure, trial_rows, component)
  key <- bundle$key
  base_key <- .likelihood_base_override_key(component)

  if (!is.null(cache_env) && is.environment(cache_env)) {
    cached <- cache_env[[key]]
    if (!is.null(cached)) return(cached)
  }

  if (identical(key, base_key)) {
    prep_val <- prep_base
  } else {
    prep_val <- .likelihood_apply_overrides(prep_base, bundle$acc, bundle$shared)
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

.likelihood_mixture_likelihood <- function(structure, prep_eval_base, trial_rows,
                                           component_ids, component_weights,
                                           outcome_label, rt_val,
                                           forced_component = NULL,
                                           prep_cache = NULL,
                                           lik_cache = NULL) {
  results <- numeric(0)
  comps <- component_ids
  if (!is.null(forced_component) && !is.na(forced_component)) {
    comps <- intersect(component_ids, as.character(forced_component))
  } else if ("component" %in% names(trial_rows)) {
    listed <- unique(trial_rows$component)
    listed <- listed[!is.na(listed)]
    if (length(listed) > 0L) {
      comps <- intersect(component_ids, listed)
    }
  }
  if (length(comps) == 0L) comps <- component_ids
  weights <- .likelihood_component_weights(
    structure,
    if ("trial" %in% names(trial_rows) && nrow(trial_rows) > 0L) trial_rows$trial[[1]] else NA_integer_,
    comps,
    structure$components$weight,
    component_weights
  )
  total <- 0.0
  for (idx in seq_along(comps)) {
    comp_id <- comps[[idx]]
    bundle <- .likelihood_component_override_bundle(structure, trial_rows, comp_id)
    trial_prep <- .likelihood_fetch_component_prep(
      structure,
      prep_eval_base,
      trial_rows,
      comp_id,
      cache_env = prep_cache,
      bundle = bundle
    )
    cache_key <- .likelihood_outcome_cache_key(bundle$key, outcome_label, rt_val)
    lik_val <- NULL
    if (!is.null(lik_cache) && is.environment(lik_cache)) {
      cached <- lik_cache[[cache_key]]
      if (!is.null(cached)) {
        lik_val <- cached
      }
    }
    if (is.null(lik_val)) {
      lik_val <- .outcome_likelihood(outcome_label, rt_val, trial_prep, comp_id)
      if (!is.null(lik_cache) && is.environment(lik_cache)) {
        lik_cache[[cache_key]] <- lik_val
      }
    }
    total <- total + weights[[idx]] * as.numeric(lik_val)
  }
  total
}

log_likelihood_from_params <- function(structure, params_df, data_df,
                                       component_weights = NULL) {
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
  prep_eval_base <- .prepare_model_for_likelihood(structure$model_spec)
  comp_ids <- structure$components$component_id

  prep_cache <- new.env(parent = emptyenv(), hash = TRUE)
  cache_across_trials <- isTRUE(getOption("uuber.param_cache_across_trials", TRUE))
  likelihood_cache <- if (cache_across_trials) new.env(parent = emptyenv(), hash = TRUE) else NULL

  if (!"trial" %in% names(params_df)) params_df$trial <- 1L
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  if (!"trial" %in% names(data_df)) {
    stop("Data frame must include a 'trial' column")
  }
  data_df$trial <- data_df$trial
  if ("component" %in% names(data_df)) {
    data_df$component <- as.character(data_df$component)
  }

  trial_ids <- sort(unique(params_df$trial))
  per_trial_loglik <- numeric(length(trial_ids))

  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_rows <- params_df[params_df$trial == tid, , drop = FALSE]
    data_row <- data_df[data_df$trial == tid, , drop = FALSE]
    if (nrow(data_row) != 1L) {
      stop(sprintf("Expected exactly one data row for trial %s", tid))
    }
    outcome <- data_row$outcome[[1]] %||% NA_character_
    rt_val <- data_row$rt[[1]] %||% NA_real_
    forced_component <- if ("component" %in% names(data_row)) data_row$component[[1]] else NULL
    lik_cache <- if (cache_across_trials && !is.null(likelihood_cache)) {
      likelihood_cache
    } else {
      new.env(parent = emptyenv(), hash = TRUE)
    }
    mixture <- .likelihood_mixture_likelihood(
      structure,
      prep_eval_base,
      trial_rows,
      comp_ids,
      component_weights,
      outcome,
      rt_val,
      forced_component = forced_component,
      prep_cache = prep_cache,
      lik_cache = lik_cache
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
