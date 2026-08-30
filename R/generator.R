# ---- Component activity helpers ----------------------------------------------

.acc_active_in_component <- function(acc_def, component) {
  comps <- acc_def$components
  !length(comps) || component %in% comps
}

.outcome_allowed_in_component <- function(options, component) {
  is.null(options$component) || component %in% options$component
}

.component_readout_count <- function(prep, component) {
  obs <- prep$observation
  base_readout <- obs$global_n_outcomes
  comp_override <- obs$component_n_outcomes[[component]]
  if (!is.null(comp_override)) {
    base_readout <- comp_override[[1L]]
  }
  min(base_readout, obs$n_outcomes)
}

# ---- Sampling primitives ------------------------------------------------------

.shared_trigger_fail <- function(ctx, trigger_id) {
  cached <- ctx$shared_trigger_state[[trigger_id]]
  if (!is.null(cached)) return(cached)
  fail <- stats::runif(1) < ctx$trial_shared_triggers[[trigger_id]]
  ctx$shared_trigger_state[[trigger_id]] <- fail
  fail
}

.resolve_effective_onset <- function(ctx, acc_id, acc_def) {
  cached <- ctx$onset_cache[[acc_id]]
  if (!is.null(cached)) {
    return(as.numeric(cached[[1]]))
  }
  spec <- acc_def$onset_spec

  if (identical(spec$kind, "absolute")) {
    onset <- spec$value
  } else {
    source_time <- if (identical(spec$source_kind, "pool")) {
      .resolve_pool(ctx, spec$source)$time
    } else {
      .get_acc_time(ctx, spec$source)
    }
    onset <- if (is.finite(source_time)) source_time + spec$lag else Inf
  }
  ctx$onset_cache[[acc_id]] <- onset
  onset
}

.sample_accumulator <- function(acc_id, acc_def, ctx) {
  onset <- .resolve_effective_onset(ctx, acc_id, acc_def)
  if (!is.finite(onset)) {
    return(Inf)
  }
  shared_id <- acc_def$shared_trigger_id
  if (!is.null(shared_id)) {
    if (.shared_trigger_fail(ctx, shared_id)) {
      return(Inf)
    }
  } else {
    success <- stats::runif(1) < (1 - acc_def$q)
    if (!success) return(Inf)
  }
  draw <- dist_registry(acc_def$dist)$r(1L, acc_def$params)
  onset + as.numeric(draw[[1]])
}

# ---- Expression evaluation ----------------------------------------------------

.get_acc_time <- function(ctx, acc_id) {
  acc_defs <- ctx$model[["accumulators"]]
  if (!.acc_active_in_component(acc_defs[[acc_id]], ctx$component)) {
    ctx$acc_times[[acc_id]] <- Inf
    return(Inf)
  }
  acc_def <- ctx$trial_accs[[acc_id]]
  if (is.null(ctx$acc_times[[acc_id]])) {
    ctx$acc_times[[acc_id]] <- .sample_accumulator(acc_id, acc_def, ctx)
  }
  ctx$acc_times[[acc_id]]
}

.members_for_component <- function(ctx, member_ids) {
  if (length(member_ids) == 0) return(character(0))
  out <- character(0)
  acc_defs <- ctx$model[["accumulators"]]
  pool_defs <- ctx$model[["pools"]]
  for (m in member_ids) {
    if (!is.null(acc_defs[[m]])) {
      if (.acc_active_in_component(acc_defs[[m]], ctx$component)) {
        out <- c(out, m)
      }
    } else {
      out <- c(out, m)
    }
  }
  out
}

.resolve_pool <- function(ctx, pool_id) {
  if (!is.null(ctx$pool_cache[[pool_id]])) return(ctx$pool_cache[[pool_id]])
  pool_defs <- ctx$model[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  members <- .members_for_component(ctx, pool_def$members)
  member_times <- list()
  cores <- numeric(0)
  acc_defs <- ctx$model[["accumulators"]]
  for (m in members) {
    if (!is.null(acc_defs[[m]])) {
      t_m <- .get_acc_time(ctx, m)
      member_times[[m]] <- t_m
      cores[m] <- t_m
    } else {
      sub <- .resolve_pool(ctx, m)
      member_times[[m]] <- sub$time
      if (length(sub$core) > 0) cores <- c(cores, sub$core)
    }
  }

  finite_vals <- as.numeric(member_times)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  k <- pool_def$k
  if (length(finite_vals) < k) {
    res <- list(time = Inf, core = cores)
  } else {
    sorted <- sort(finite_vals, partial = k)
    res <- list(time = sorted[[k]], core = cores)
  }
  ctx$pool_cache[[pool_id]] <- res
  res
}

.event_key <- function(ev) {
  as.character(ev$source)
}

.resolve_event <- function(ctx, ev) {
  key <- .event_key(ev)
  if (!is.null(ctx$event_cache[[key]])) return(ctx$event_cache[[key]])
  source_id <- ev$source
  if (!is.null(ctx$model[["pools"]][[source_id]])) {
    result <- .resolve_pool(ctx, source_id)
  } else if (!is.null(ctx$model[["accumulators"]][[source_id]])) {
    acc_def <- ctx$model[["accumulators"]][[source_id]]
    if (!.acc_active_in_component(acc_def, ctx$component)) {
      result <- list(time = Inf, core = numeric(0))
    } else {
      t_acc <- .get_acc_time(ctx, source_id)
      result <- list(time = t_acc, core = setNames(t_acc, source_id))
    }
  }
  ctx$event_cache[[key]] <- result
  result
}

.eval_expr <- function(expr, ctx) {
  kind <- expr$kind
  if (identical(kind, "event")) {
    ev <- list(source = expr$source)
    return(.resolve_event(ctx, ev))
  }
  if (identical(kind, "and")) {
    inputs <- lapply(expr$args, .eval_expr, ctx = ctx)
    times <- vapply(inputs, function(x) x$time, numeric(1))
    cores <- unlist(lapply(inputs, function(x) x$core), use.names = TRUE)
    if (any(!is.finite(times))) {
      return(list(time = Inf, core = cores))
    }
    return(list(time = max(times), core = cores))
  }
  if (identical(kind, "or")) {
    inputs <- lapply(expr$args, .eval_expr, ctx = ctx)
    times <- vapply(inputs, function(x) x$time, numeric(1))
    finite_idx <- which(is.finite(times))
    if (length(finite_idx) == 0) return(list(time = Inf, core = numeric(0)))
    times <- times[finite_idx]
    inputs <- inputs[finite_idx]
    tmin <- min(times)
    tied <- which(times == tmin)
    chosen <- tied[[1]]
    if (length(tied) > 1) {
      tie_vals <- vapply(tied, function(idx) {
        core <- inputs[[idx]]$core
        if (length(core) == 0) inputs[[idx]]$time else min(core)
      }, numeric(1))
      chosen <- tied[which.min(tie_vals)]
    }
    return(list(time = times[[chosen]], core = inputs[[chosen]]$core))
  }
  if (identical(kind, "not")) {
    inner <- .eval_expr(expr$arg, ctx)
    if (is.finite(inner$time)) {
      return(list(time = Inf, core = numeric(0)))
    } else {
      return(list(time = 0, core = numeric(0)))
    }
  }
  if (identical(kind, "guard")) {
    ref <- .eval_expr(expr$reference, ctx)
    if (!is.finite(ref$time)) return(list(time = Inf, core = ref$core))
    blocker <- .eval_expr(expr$blocker, ctx)
    if (is.finite(blocker$time) && blocker$time < ref$time) {
      return(list(time = Inf, core = ref$core))
    } else {
      return(list(time = ref$time, core = ref$core))
    }
  }
  stop(sprintf("Unsupported expression kind '%s'", kind))
}

.evaluate_outcomes <- function(ctx) {
  outs <- ctx$model[["outcomes"]]
  res <- vector("list", length(outs))
  labels <- names(outs)
  names(res) <- labels
  for (i in seq_along(outs)) {
    def <- outs[[i]]
    options <- def$options
    if (!.outcome_allowed_in_component(options, ctx$component)) {
      res[[i]] <- list(time = Inf, core = numeric(0), options = options)
      next
    }
    expr <- def[["expr"]]
    eval <- .eval_expr(expr, ctx)
    res[[i]] <- list(
      time = eval$time,
      core = eval$core,
      options = options
    )
  }
  res
}

# ---- Trial simulation ---------------------------------------------------------

.simulate_trial <- function(prep, component, trial_accs,
                            trial_shared_triggers, n_outcomes,
                            keep_detail = FALSE) {
  ctx <- list(
    model = prep,
    component = component,
    acc_times = new.env(parent = emptyenv()),
    onset_cache = new.env(parent = emptyenv()),
    pool_cache = new.env(parent = emptyenv()),
    event_cache = new.env(parent = emptyenv()),
    trial_accs = trial_accs,
    trial_shared_triggers = trial_shared_triggers,
    shared_trigger_state = new.env(parent = emptyenv())
  )

  outcomes <- .evaluate_outcomes(ctx)
  outcome_labels <- names(outcomes)
  cand_labels <- character(0)
  cand_times <- numeric(0)
  cand_core <- list()
  cand_options <- list()

  for (i in seq_along(outcomes)) {
    entry <- outcomes[[i]]
    time <- entry$time
    if (!is.finite(time)) next
    cand_labels <- c(cand_labels, outcome_labels[[i]])
    cand_times <- c(cand_times, time)
    cand_core[[length(cand_core) + 1L]] <- entry$core
    cand_options[[length(cand_options) + 1L]] <- entry$options
  }

  if (length(cand_times) == 0L) {
    return(list(
      outcomes = rep(NA_character_, n_outcomes),
      rts = rep(NA_real_, n_outcomes),
      detail = NULL
    ))
  }

  if (n_outcomes > 1L) {
    order_idx <- order(cand_times, seq_along(cand_times))
    keep <- min(length(order_idx), n_outcomes)
    ranked_labels <- rep(NA_character_, n_outcomes)
    ranked_times <- rep(NA_real_, n_outcomes)
    ranked_labels[seq_len(keep)] <- cand_labels[order_idx[seq_len(keep)]]
    ranked_times[seq_len(keep)] <- cand_times[order_idx[seq_len(keep)]]

    detail <- NULL
    if (isTRUE(keep_detail)) {
      detail <- list(
        component = component,
        acc_times = as.list(ctx$acc_times),
        pool_times = as.list(ctx$pool_cache),
        event_times = as.list(ctx$event_cache),
        outcome_candidates = data.frame(label = cand_labels, time = cand_times, stringsAsFactors = FALSE),
        ranked_outcomes = data.frame(
          rank = seq_len(n_outcomes),
          label = ranked_labels,
          time = ranked_times,
          stringsAsFactors = FALSE
        )
      )
    }

    return(list(
      outcomes = ranked_labels,
      rts = ranked_times,
      detail = detail
    ))
  }

  tmin <- min(cand_times)
  tied <- which(cand_times == tmin)
  chosen_idx <- tied[[1]]
  if (length(tied) > 1) {
    tie_scores <- vapply(tied, function(i) {
      core <- cand_core[[i]]
      if (length(core) == 0) cand_times[[i]] else min(core)
    }, numeric(1))
    chosen_idx <- tied[which.min(tie_scores)]
  }

  chosen_label <- cand_labels[[chosen_idx]]
  chosen_time <- cand_times[[chosen_idx]]
  options <- cand_options[[chosen_idx]]

  # Outcome-level guess policies
  if (!is.null(options$guess)) {
    gp <- options$guess
    draw <- sample(gp$labels, size = 1L, prob = gp$weights)
    if (!is.null(gp$rt_policy) && identical(gp$rt_policy, "na")) {
      chosen_time <- NA_real_
    }
    chosen_label <- draw
  }

  if (!is.null(options$map_outcome_to)) {
    # If mapping to NA, also drop the response time
    if (is.na(options$map_outcome_to)) {
      chosen_label <- NA_character_
      chosen_time <- NA_real_
    } else {
      chosen_label <- options$map_outcome_to
    }
  }

  detail <- NULL
  if (isTRUE(keep_detail)) {
    detail <- list(
      component = component,
      acc_times = as.list(ctx$acc_times),
      pool_times = as.list(ctx$pool_cache),
      event_times = as.list(ctx$event_cache),
      outcome_candidates = data.frame(label = cand_labels, time = cand_times, stringsAsFactors = FALSE)
    )
  }

  list(
    outcomes = chosen_label,
    rts = chosen_time,
    detail = detail
  )
}

# ---- Simulation helpers ------------------------------------------------------

.normalize_trial_df <- function(trial_df, n_trials, acc_ids, comp_ids) {
  if (is.null(trial_df)) {
    return(list(
      component = rep.int(NA_character_, n_trials),
      onset = NULL
    ))
  }

  trial_df <- as.data.frame(trial_df)
  if (!all(c("trials", "racer") %in% names(trial_df))) {
    stop("trial_df must contain trials and racer columns", call. = FALSE)
  }
  n_acc <- length(acc_ids)
  expected_trials <- rep.int(seq_len(n_trials), rep.int(n_acc, n_trials))
  if (!is.numeric(trial_df$trials) ||
      length(trial_df$trials) != length(expected_trials) ||
      any(!is.finite(trial_df$trials)) ||
      any(trial_df$trials != expected_trials) ||
      !identical(as.character(trial_df$racer), rep(acc_ids, times = n_trials))) {
    stop(
      "trial_df must contain one complete trials-by-accumulator block in model order",
      call. = FALSE
    )
  }

  component <- rep.int(NA_character_, n_trials)
  if ("component" %in% names(trial_df)) {
    values <- as.character(trial_df$component)
    component <- values[seq.int(1L, length(values), by = n_acc)]
    if (!identical(values, rep.int(component, rep.int(n_acc, n_trials)))) {
      stop("trial_df$component must be constant within each trial", call. = FALSE)
    }
    unknown <- unique(component[!is.na(component) & !component %in% comp_ids])
    if (length(unknown)) {
      stop("Unknown trial_df component: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  onset <- NULL
  if ("onset" %in% names(trial_df)) {
    if (!is.numeric(trial_df$onset) || any(!is.na(trial_df$onset) & !is.finite(trial_df$onset))) {
      stop("trial_df$onset must be numeric", call. = FALSE)
    }
    onset <- matrix(
      trial_df$onset,
      nrow = n_trials,
      byrow = TRUE,
      dimnames = list(NULL, acc_ids)
    )
  }

  list(component = component, onset = onset)
}

# ---- Public API ----------------------------------------------------------------

#' Simulate behavioral data from a model
#'
#' @param structure Finalized model structure.
#' @param params_df Trial-level parameter values.
#' @param trial_df Optional trial-by-accumulator conditioning table containing
#'   complete `trials`/`racer` blocks in model order. It may supply per-trial
#'   `component` and per-accumulator `onset` values.
#' @param seed Optional random-number seed.
#' @param keep_detail If `TRUE`, keep additional simulation detail.
#' @param keep_component Whether to keep the chosen mixture component in the
#'   output when the model has multiple components. If `NULL`, fixed mixtures
#'   keep the component label and sampled mixtures drop it.
#' @return A data frame of simulated behavioral data. For standard models this
#'   includes `trials`, `R`, and `rt`. If `n_outcomes > 1`, additional ordered
#'   response columns such as `R2`/`rt2` are included.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params <- c(m = 0, s = 0.1)
#' df <- build_param_matrix(structure, params, n_trials = 3)
#' simulate(structure, df, seed = 123)
#' @export
simulate <- function(structure,
                     params_df,
                     trial_df = NULL,
                     seed = NULL,
                     keep_detail = FALSE,
                     keep_component = NULL) {
  prep <- structure$prep
  acc_defs <- prep$accumulators
  acc_ids <- names(acc_defs)
  n_acc <- length(acc_ids)
  params_mat <- params_df
  n_trials <- nrow(params_mat) %/% n_acc
  p_cols <- grep("^p[0-9]+$", colnames(params_mat), value = TRUE)
  p_cols <- p_cols[order(as.integer(sub("^p", "", p_cols)))]

  components <- prep$components
  comp_ids <- components$ids
  mix_mode <- components$mode
  component_weight_names <- vapply(comp_ids, function(component) {
    components$attrs[[component]]$weight_param %||% NA_character_
  }, character(1))

  trial_info <- .normalize_trial_df(trial_df, n_trials, acc_ids, comp_ids)
  component_vec <- trial_info$component
  onset_mat <- trial_info$onset

  row_map <- lapply(seq_len(n_trials), function(trial) {
    setNames((trial - 1L) * n_acc + seq_len(n_acc), acc_ids)
  })
  if (!is.null(seed)) set.seed(seed)
  trial_ids <- seq_len(n_trials)
  n_readout <- prep$observation$n_outcomes
  outcomes <- matrix(NA_character_, nrow = length(trial_ids), ncol = n_readout)
  rts <- matrix(NA_real_, nrow = length(trial_ids), ncol = n_readout)
  comp_record <- rep(NA_character_, length(trial_ids))
  details <- if (keep_detail) vector("list", length(trial_ids)) else NULL

  acc_idx_map <- setNames(seq_along(acc_ids), acc_ids)

  for (i in seq_along(trial_ids)) {
    forced_component <- component_vec[[i]]
    if (!is.na(forced_component)) {
      chosen_component <- forced_component
    } else if (length(comp_ids) == 1L) {
      chosen_component <- comp_ids[[1]]
    } else {
      comp_weights <- components$weights
      if (identical(mix_mode, "sample")) {
        sampled <- !is.na(component_weight_names)
        comp_weights[sampled] <- params_mat[
          row_map[[i]][[1L]],
          component_weight_names[sampled]
        ]
        comp_weights[!sampled] <- 1 - sum(comp_weights[sampled])
      }
      chosen_component <- sample(comp_ids, size = 1L, prob = comp_weights)
    }
    comp_record[[i]] <- chosen_component

    trial_accs <- list()
    shared_map <- list()
    trial_rows <- row_map[[i]]

    for (acc_id in names(trial_rows)) {
      row_idx <- trial_rows[[acc_id]]
      row_vals <- params_mat[row_idx, ]
      dist_params <- dist_param_names(acc_defs[[acc_id]]$dist)
      dist_list <- setNames(as.list(row_vals[p_cols][seq_along(dist_params)]), dist_params)
      dist_list$t0 <- row_vals[["t0"]]
      onset_override <- NA_real_
      if (!is.null(onset_mat)) {
        onset_override <- as.numeric(onset_mat[i, acc_idx_map[[acc_id]]])
      }
      onset_spec <- acc_defs[[acc_id]]$onset_spec
      if (identical(onset_spec$kind, "absolute")) {
        onset_spec$value <- if (is.finite(onset_override)) {
          onset_override
        } else {
          onset_spec$value
        }
      } else {
        lag_val <- onset_spec$lag
        if (is.finite(onset_override)) {
          lag_val <- lag_val + onset_override
        }
        onset_spec <- list(
          kind = "after",
          source = onset_spec$source,
          source_kind = onset_spec$source_kind,
          lag = lag_val
        )
      }
      acc_entry <- list(
        dist = acc_defs[[acc_id]]$dist,
        params = dist_list,
        onset_spec = onset_spec,
        q = row_vals[["q"]],
        shared_trigger_id = acc_defs[[acc_id]]$shared_trigger_id
      )
      trial_accs[[acc_id]] <- acc_entry
      stid <- acc_entry$shared_trigger_id
      if (!is.null(stid)) {
        if (is.null(shared_map[[stid]])) {
          shared_map[[stid]] <- row_vals[["q"]]
        }
      }
    }

    result <- .simulate_trial(
      prep,
      chosen_component,
      keep_detail = keep_detail,
      trial_accs = trial_accs,
      trial_shared_triggers = shared_map,
      n_outcomes = .component_readout_count(prep, chosen_component)
    )
    out_vals <- result$outcomes
    rt_vals <- result$rts
    if (length(out_vals) > 0L) outcomes[i, seq_along(out_vals)] <- as.character(out_vals)
    if (length(rt_vals) > 0L) rts[i, seq_along(rt_vals)] <- as.numeric(rt_vals)
    if (keep_detail) details[[i]] <- result$detail
  }

  out_df <- data.frame(
    trials = trial_ids,
    R = outcomes[, 1L],
    rt = rts[, 1L],
    stringsAsFactors = FALSE
  )
  if (n_readout > 1L) {
    for (rank_idx in 2:n_readout) {
      out_df[[paste0("R", rank_idx)]] <- outcomes[, rank_idx]
      out_df[[paste0("rt", rank_idx)]] <- rts[, rank_idx]
    }
  }
  keep_component <- keep_component %||% if (identical(mix_mode, "sample")) FALSE else TRUE
  if (keep_component && (length(comp_ids) > 1L)) out_df$component <- comp_record
  if (keep_detail) attr(out_df, "details") <- details
  out_df
}
