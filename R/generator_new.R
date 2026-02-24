# Interpreter-driven race simulator for the boolean-expression DSL defined in new_API.R
# Depends on core distribution and sampling helpers from the existing codebase.

# ---- Component activity helpers ----------------------------------------------

.acc_active_in_component <- function(acc_def, component) {
  comps <- acc_def$components
  if (length(comps) == 0 || is.null(component) || identical(component, "__default__")) return(TRUE)
  component %in% comps
}

.get_component_attr <- function(prep, component, name) {
  comp_info <- prep[["components"]]
  if (is.null(comp_info)) return(NULL)
  attrs <- comp_info$attrs[[component]]
  if (is.null(attrs)) return(NULL)
  attrs[[name]]
}

.outcome_allowed_in_component <- function(options, component) {
  comps <- options$component %||% NULL
  if (is.null(comps) || length(comps) == 0L) return(TRUE)
  comp_label <- component %||% "__default__"
  if (is.na(comp_label) || !nzchar(comp_label)) return(FALSE)
  comp_label %in% as.character(comps)
}

.component_readout_count <- function(prep, component) {
  obs <- prep$observation %||% list()
  max_readout <- obs$n_outcomes %||% 1L
  max_readout <- as.integer(max_readout)
  if (is.na(max_readout) || max_readout < 1L) {
    max_readout <- 1L
  }
  base_readout <- obs$global_n_outcomes %||% max_readout
  base_readout <- as.integer(base_readout)
  if (is.na(base_readout) || base_readout < 1L) {
    base_readout <- 1L
  }
  comp_label <- component %||% "__default__"
  overrides <- obs$component_n_outcomes %||% list()
  comp_override <- overrides[[comp_label]] %||% NULL
  if (!is.null(comp_override)) {
    base_readout <- as.integer(comp_override[[1]])
    if (is.na(base_readout) || base_readout < 1L) {
      base_readout <- 1L
    }
  }
  min(base_readout, max_readout)
}

# ---- Sampling primitives ------------------------------------------------------

.shared_trigger_fail <- function(ctx, trigger_id) {
  if (is.null(trigger_id) || is.na(trigger_id) || trigger_id == "") return(FALSE)
  cached <- ctx$shared_trigger_state[[trigger_id]]
  if (!is.null(cached)) return(isTRUE(cached$fail))
  base_info <- ctx$model[["shared_triggers"]][[trigger_id]] %||% list()
  prob <- base_info$q %||% 0
  override <- ctx$trial_shared_triggers[[trigger_id]] %||% list()
  if (!is.null(override$prob) && !is.na(override$prob)) {
    prob <- as.numeric(override$prob)
  }
  if (!is.numeric(prob) || length(prob) != 1L || prob < 0 || prob > 1) {
    stop(sprintf("Shared trigger '%s' requires a probability between 0 and 1", trigger_id))
  }
  fail <- stats::runif(1) < prob
  ctx$shared_trigger_state[[trigger_id]] <- list(fail = fail, prob = prob)
  fail
}

.resolve_effective_onset <- function(ctx, acc_id, acc_def) {
  cached <- ctx$onset_cache[[acc_id]]
  if (!is.null(cached)) {
    return(as.numeric(cached[[1]]))
  }
  if (isTRUE(ctx$onset_inflight[[acc_id]])) {
    stop(sprintf("Chained onset cycle encountered while resolving '%s'", acc_id))
  }
  ctx$onset_inflight[[acc_id]] <- TRUE
  on.exit(rm(list = acc_id, envir = ctx$onset_inflight), add = TRUE)

  spec <- acc_def$onset_spec %||% list(kind = "absolute", value = acc_def$onset %||% 0)
  kind <- spec$kind %||% "absolute"
  onset <- NA_real_

  if (identical(kind, "absolute")) {
    onset <- as.numeric(spec$value %||% acc_def$onset %||% 0)[1]
  } else if (identical(kind, "after")) {
    source <- spec$source %||% NULL
    source_kind <- spec$source_kind %||% NULL
    lag <- as.numeric(spec$lag %||% 0)[1]
    if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
      stop(sprintf("Accumulator '%s' has invalid chained onset source", acc_id))
    }
    if (!is.finite(lag) || lag < 0) {
      stop(sprintf("Accumulator '%s' has invalid chained onset lag", acc_id))
    }

    source_time <- if (identical(source_kind, "pool")) {
      .resolve_pool(ctx, source)$time
    } else if (identical(source_kind, "accumulator")) {
      .get_acc_time(ctx, source)
    } else if (!is.null(ctx$model[["pools"]][[source]])) {
      .resolve_pool(ctx, source)$time
    } else if (!is.null(ctx$model[["accumulators"]][[source]])) {
      .get_acc_time(ctx, source)
    } else {
      stop(sprintf("Accumulator '%s' references unknown chained onset source '%s'", acc_id, source))
    }

    onset <- if (is.finite(source_time)) source_time + lag else Inf
  } else {
    stop(sprintf("Accumulator '%s' has unsupported onset kind '%s'", acc_id, kind))
  }

  if (is.na(onset)) {
    stop(sprintf("Accumulator '%s' produced NA onset", acc_id))
  }
  if (!is.finite(onset) && !is.infinite(onset)) {
    stop(sprintf("Accumulator '%s' produced non-finite onset", acc_id))
  }
  ctx$onset_cache[[acc_id]] <- onset
  onset
}

.sample_accumulator <- function(acc_id, acc_def, ctx) {
  onset <- .resolve_effective_onset(ctx, acc_id, acc_def)
  if (!is.finite(onset)) {
    return(Inf)
  }
  shared_id <- acc_def$shared_trigger_id %||% NULL
  if (!is.null(shared_id) && !is.na(shared_id) && shared_id != "") {
    if (.shared_trigger_fail(ctx, shared_id)) {
      return(Inf)
    }
  } else {
    success <- stats::runif(1) < (1 - acc_def$q)
    if (!success) return(Inf)
  }
  reg <- dist_registry(acc_def$dist)
  if (is.null(reg) || is.null(reg$r)) stop(sprintf("No sampler registered for distribution '%s'", acc_def$dist))
  draw <- reg$r(1L, acc_def$params)
  onset + as.numeric(draw[[1]])
}

# ---- Expression evaluation ----------------------------------------------------

.get_acc_time <- function(ctx, acc_id) {
  acc_defs <- ctx$model[["accumulators"]]
  if (!acc_id %in% names(acc_defs)) stop(sprintf("Unknown accumulator '%s'", acc_id))
  if (!.acc_active_in_component(acc_defs[[acc_id]], ctx$component)) {
    ctx$acc_times[[acc_id]] <- Inf
    return(Inf)
  }
  acc_def <- ctx$trial_accs[[acc_id]] %||% acc_defs[[acc_id]]
  if (is.null(acc_def)) stop(sprintf("Accumulator '%s' is not available in this trial context", acc_id))
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
    } else if (!is.null(pool_defs[[m]])) {
      out <- c(out, m)
    }
  }
  out
}

.resolve_pool <- function(ctx, pool_id) {
  if (!is.null(ctx$pool_cache[[pool_id]])) return(ctx$pool_cache[[pool_id]])
  pool_defs <- ctx$model[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  members <- .members_for_component(ctx, pool_def$members)
  member_times <- list()
  cores <- numeric(0)
  acc_defs <- ctx$model[["accumulators"]]
  for (m in members) {
    if (!is.null(acc_defs[[m]])) {
      t_m <- .get_acc_time(ctx, m)
      member_times[[m]] <- t_m
      cores[m] <- t_m
    } else if (!is.null(pool_defs[[m]])) {
      sub <- .resolve_pool(ctx, m)
      member_times[[m]] <- sub$time
      if (length(sub$core) > 0) cores <- c(cores, sub$core)
    } else {
      stop(sprintf("Pool '%s' references unknown member '%s'", pool_id, m))
    }
  }

  finite_vals <- as.numeric(member_times)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  k <- pool_def$k
  if (length(finite_vals) < k || k < 1) {
    res <- list(time = Inf, core = cores)
  } else {
    sorted <- sort(finite_vals, partial = k)
    res <- list(time = sorted[[k]], core = cores)
  }
  ctx$pool_cache[[pool_id]] <- res
  res
}

.event_key <- function(ev) {
  paste(ev$source, ev$k %||% "", sep = "::")
}

.resolve_event <- function(ctx, ev) {
  key <- .event_key(ev)
  if (!is.null(ctx$event_cache[[key]])) return(ctx$event_cache[[key]])
  source_id <- ev$source
  if (identical(source_id, "__DEADLINE__")) {
    result <- list(time = Inf, core = numeric(0))
  } else if (identical(source_id, "__GUESS__")) {
    result <- list(time = Inf, core = numeric(0))
  } else if (!is.null(ctx$model[["pools"]][[source_id]])) {
    result <- .resolve_pool(ctx, source_id)
  } else if (!is.null(ctx$model[["accumulators"]][[source_id]])) {
    acc_def <- ctx$model[["accumulators"]][[source_id]]
    if (!.acc_active_in_component(acc_def, ctx$component)) {
      result <- list(time = Inf, core = numeric(0))
    } else {
      t_acc <- .get_acc_time(ctx, source_id)
      result <- list(time = t_acc, core = setNames(t_acc, source_id))
    }
  } else {
    result <- list(time = Inf, core = numeric(0))
  }
  ctx$event_cache[[key]] <- result
  result
}

.eval_expr <- function(expr, ctx) {
  kind <- expr$kind
  if (identical(kind, "event")) {
    ev <- list(source = expr$source, k = expr$k %||% NULL)
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
      list(time = Inf, core = numeric(0))
    } else {
      list(time = 0, core = numeric(0))
    }
  }
  if (identical(kind, "guard")) {
    ref <- .eval_expr(expr$reference, ctx)
    if (!is.finite(ref$time)) return(list(time = Inf, core = ref$core))
    blocker <- .eval_expr(expr$blocker, ctx)
    unless_list <- expr$unless %||% list()
    guard_blocked <- FALSE
    if (length(unless_list) > 0) {
      for (unl in unless_list) {
        unl_eval <- .eval_expr(unl, ctx)
        if (is.finite(unl_eval$time) && unl_eval$time <= blocker$time) {
          guard_blocked <- TRUE
          break
        }
      }
    }
    if (!guard_blocked && is.finite(blocker$time) && blocker$time < ref$time) {
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
    def <- outs[[i]] %||% list()
    options <- def[["options"]] %||% list()
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

.resolve_param_row_acc_id <- function(row, acc_ids) {
  pick_id <- function(val) {
    if (is.null(val) || length(val) == 0L) return(NULL)
    v <- val[[1]]
    if (is.null(v) || (length(v) == 1L && is.na(v))) return(NULL)
    v
  }
  raw <- if ("accumulator" %in% names(row)) pick_id(row[["accumulator"]]) else NULL
  if (!is.null(raw)) {
    if (is.numeric(raw)) {
      idx <- as.integer(raw)
      if (!is.na(idx) && idx >= 1L && idx <= length(acc_ids)) {
        return(acc_ids[[idx]])
      }
    }
  }
  stop("Parameter rows must include a valid numeric 'accumulator' column")
}

.param_override_columns <- function(params_rows) {
  base_cols <- c(
    "trial", "component", "accumulator",
    "type",
    "outcome", "rt", "params", "condition",
    "component_weight"
  )
  setdiff(names(params_rows), c(base_cols, "onset", "q", "shared_trigger_id"))
}

.generator_component_rows <- function(trial_rows, component) {
  if (is.null(trial_rows) || nrow(trial_rows) == 0L) return(trial_rows)
  if (!"component" %in% names(trial_rows)) return(trial_rows)
  comp_label <- component %||% "__default__"
  mask <- is.na(trial_rows$component)
  if (!identical(comp_label, "__default__")) {
    mask <- mask | trial_rows$component == comp_label
  }
  trial_rows[mask, , drop = FALSE]
}

.component_weights_for_trial <- function(structure, trial_id, available_components,
                                         trial_rows = NULL) {
  comp_table <- structure$components
  base_weights <- comp_table$weight[match(available_components, comp_table$component_id)]
  mode <- comp_table$mode[[1]] %||% "fixed"
  if (!identical(mode, "sample")) {
    if (length(base_weights) != length(available_components) || all(is.na(base_weights))) {
      base_weights <- rep(1 / length(available_components), length(available_components))
    } else {
      if (any(!is.finite(base_weights) | base_weights < 0)) {
        stop("Component weights must be non-negative finite numbers")
      }
      total <- sum(base_weights)
      if (!is.finite(total) || total <= 0) {
        base_weights <- rep(1 / length(available_components), length(available_components))
      } else {
        base_weights <- base_weights / total
      }
    }
    return(base_weights)
  }
  # mode == "sample": use per-component weight_param values (probabilities) for non-reference,
  # reference weight is residual.
  comp_attrs <- comp_table$attrs
  ref_id <- comp_table$reference[[1]] %||% available_components[[1]]
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
    if (is.na(val)) {
      # fallback to base weight if provided
      idx <- which(comp_table$component_id == cid)
      if (length(idx) == 1) val <- comp_table$weight[[idx]]
    }
    if (is.na(val) || !is.finite(val) || val < 0 || val > 1) {
      stop(sprintf("Mixture weight '%s' must be a probability in [0,1]", wp))
    }
    weights[[i]] <- val
    sum_nonref <- sum_nonref + val
  }
  # reference weight
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

.generator_param_state_from_rows <- function(prep, params_rows) {
  if (is.null(params_rows) || nrow(params_rows) == 0L) {
    return(list(accs = list(), shared = list()))
  }
  acc_defs <- prep$accumulators %||% list()
  if (length(acc_defs) == 0L) return(list(accs = list(), shared = list()))
  acc_ids <- names(acc_defs)
  param_cols <- .param_override_columns(params_rows)
  acc_entries <- list()
  shared_entries <- list()

  # Shared trigger parameter overrides (group-level param names)
  if (!is.null(prep$shared_triggers) && length(prep$shared_triggers) > 0) {
    for (trig in prep$shared_triggers) {
      param_name <- trig$param %||% NULL
      if (is.null(param_name) || !nzchar(param_name)) next
      if (!param_name %in% names(params_rows)) next
      vals <- params_rows[[param_name]]
      if (is.null(vals) || length(vals) == 0L) next
      val <- vals[[1]]
      if (length(val) == 1L && !is.na(val)) {
        if (val < 0 || val > 1) {
          stop(sprintf("Shared trigger '%s' probability must be in [0,1]", param_name))
        }
        prob <- as.numeric(val)
        shared_entries[[trig$id]] <- list(prob = prob)
      }
    }
  }

  for (i in seq_len(nrow(params_rows))) {
    row <- params_rows[i, , drop = FALSE]
    entry <- .acc_override_entry(row, acc_defs, acc_ids, param_cols)
    acc_entries[[entry$acc_id]] <- entry$acc
    if (!is.null(entry$shared_id) && !is.na(entry$shared_id) && nzchar(entry$shared_id)) {
      shared_entries[[entry$shared_id]] <- entry$shared
    }
  }
  list(accs = acc_entries, shared = shared_entries)
}

.acc_override_entry <- function(row, acc_defs, acc_ids, param_cols) {
  acc_id <- .resolve_param_row_acc_id(row, acc_ids)
  base_def <- acc_defs[[acc_id]]
  if (is.null(base_def)) {
    stop(sprintf("No accumulator definition found for '%s'", acc_id))
  }
  override_params <- base_def$params %||% list()
  if ("params" %in% names(row)) {
    entry <- row$params[[1]]
    if (!is.null(entry)) {
      if (!is.list(entry)) stop("Column 'params' must contain lists of parameter values")
      override_params <- modifyList(override_params, entry, keep.null = TRUE)
    }
  }
  if (length(param_cols) > 0L) {
    for (col in param_cols) {
      vals <- row[[col]]
      if (is.null(vals) || length(vals) == 0L) next
      val <- vals[[1]]
      if (length(val) == 1L && is.na(val)) next
      override_params[[col]] <- val
    }
  }
  override_params <- .ensure_acc_param_t0(override_params)
  onset_val <- row$onset
  if (!is.null(onset_val) && length(onset_val) > 0L && !is.na(onset_val[[1]])) {
    onset_val <- as.numeric(onset_val[[1]])
  } else {
    onset_val <- base_def$onset %||% 0
  }
  # Distinguish between gate probability (shared trigger) and per-acc q
  has_shared <- !is.null(base_def$shared_trigger_id) && !is.na(base_def$shared_trigger_id)
  q_raw <- row$q
  q_override <- NULL
  if (!is.null(q_raw) && length(q_raw) > 0L && !is.na(q_raw[[1]])) {
    tmp <- as.numeric(q_raw[[1]])
    if (tmp < 0 || tmp > 1) {
      stop("q must be a probability in [0,1]")
    }
    if (has_shared && identical(tmp, base_def$q %||% 0)) {
      q_override <- NULL
    } else {
      q_override <- tmp
    }
  }
  base_shared_prob <- if (!is.null(base_def$shared_trigger_q)) as.numeric(base_def$shared_trigger_q) else NULL
  base_q_prob <- if (!is.null(base_def$q)) as.numeric(base_def$q) else NULL
  gate_prob <- q_override %||% base_shared_prob %||% base_q_prob %||% 0
  q_val <- if (has_shared) 0 else (q_override %||% base_q_prob %||% 0)
  shared_id <- NULL
  if ("shared_trigger_id" %in% names(row)) {
    shared_raw <- row$shared_trigger_id
    if (!is.null(shared_raw) && length(shared_raw) > 0L) {
      val <- shared_raw[[1]]
      if (!is.null(val) && !(length(val) == 1L && is.na(val))) {
        shared_id <- as.character(val)
      }
    }
  }
  base_shared <- base_def$shared_trigger_id %||% NA_character_
  trig_id <- shared_id %||% base_shared %||% NA_character_
  acc_entry <- list(
    dist = base_def$dist,
    params = override_params,
    onset = onset_val,
    q = q_val,
    shared_trigger_id = trig_id,
    components = base_def$components %||% character(0)
  )
  shared_entry <- NULL
  if (!is.null(trig_id) && !is.na(trig_id) && nzchar(trig_id)) {
    shared_entry <- list(prob = as.numeric(gate_prob))
  }
  list(acc_id = acc_id, acc = acc_entry, shared_id = trig_id, shared = shared_entry)
}

# ---- Trial simulation ---------------------------------------------------------

.simulate_trial <- function(prep, component, keep_detail = FALSE,
                            trial_accs = NULL,
                            trial_shared_triggers = NULL,
                            n_outcomes = NULL) {
  comp_info <- prep[["components"]]
  component_attr <- comp_info$attrs[[component]]
  if (is.null(component_attr)) component_attr <- list()
  ctx <- list(
    model = prep,
    component = component,
    acc_times = new.env(parent = emptyenv()),
    onset_cache = new.env(parent = emptyenv()),
    onset_inflight = new.env(parent = emptyenv()),
    pool_cache = new.env(parent = emptyenv()),
    event_cache = new.env(parent = emptyenv()),
    trial_accs = trial_accs %||% list(),
    trial_shared_triggers = trial_shared_triggers %||% list(),
    shared_trigger_state = new.env(parent = emptyenv())
  )

  outcomes <- .evaluate_outcomes(ctx)
  if (is.null(n_outcomes)) {
    n_outcomes <- .component_readout_count(prep, component)
  } else {
    n_outcomes <- as.integer(n_outcomes)
    if (is.na(n_outcomes) || n_outcomes < 1L) {
      n_outcomes <- 1L
    }
  }
  outcome_labels <- names(outcomes)
  cand_labels <- character(0)
  cand_times <- numeric(0)
  cand_core <- list()
  cand_options <- list()

  for (i in seq_along(outcomes)) {
    entry <- outcomes[[i]]
    if (is.null(entry)) next
    time <- entry$time
    if (!is.finite(time)) next
    cand_labels <- c(cand_labels, outcome_labels[[i]])
    cand_times <- c(cand_times, time)
    cand_core[[length(cand_core) + 1L]] <- entry$core
    cand_options[[length(cand_options) + 1L]] <- entry$options
  }

  if (length(cand_times) == 0L) {
    result <- list(
      outcome = NA_character_,
      rt = NA_real_,
      outcomes = rep(NA_character_, n_outcomes),
      rts = rep(NA_real_, n_outcomes),
      detail = NULL
    )
    return(result)
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
      outcome = ranked_labels[[1]],
      rt = ranked_times[[1]],
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
    labels <- gp$labels
    weights <- gp$weights
    if (is.null(labels) || is.null(weights)) stop("Guess policy requires labels and weights")
    draw <- sample(labels, size = 1L, prob = weights)
    if (!is.null(gp$rt_policy) && identical(gp$rt_policy, "na")) {
      chosen_time <- NA_real_
    }
    chosen_label <- draw
  }

  # Component-level guess overrides
  if (!is.null(component_attr$guess)) {
    guess <- component_attr$guess
    w_keep <- guess$weights[[chosen_label]] %||% 1
    keep <- stats::runif(1) <= w_keep
    if (!keep) {
      chosen_label <- guess$outcome %||% chosen_label
    }
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
    outcome = chosen_label,
    rt = chosen_time,
    outcomes = chosen_label,
    rts = chosen_time,
    detail = detail
  )
}

# ---- Simulation helpers ------------------------------------------------------

.normalize_trial_df <- function(trial_df, acc_ids, comp_ids) {
  if (is.null(trial_df)) {
    return(list(
      trial_df = NULL,
      trial_has_acc = FALSE,
      trial_has_component = FALSE,
      trial_has_onset = FALSE
    ))
  }

  trial_df <- as.data.frame(trial_df)
  trial_has_acc <- "accumulator" %in% names(trial_df)
  if (!"trial" %in% names(trial_df)) {
    if (trial_has_acc) stop("trial_df with accumulator must include a trial column")
    trial_df$trial <- seq_len(nrow(trial_df))
  }
  if (is.factor(trial_df$trial)) trial_df$trial <- as.character(trial_df$trial)
  trial_df$trial <- suppressWarnings(as.integer(trial_df$trial))
  if (any(is.na(trial_df$trial)) || any(trial_df$trial < 1L)) {
    stop("trial_df$trial must contain positive integers")
  }

  if (trial_has_acc) {
    acc_raw <- trial_df$accumulator
    if (is.factor(acc_raw)) acc_raw <- as.character(acc_raw)
    if (is.numeric(acc_raw)) {
      acc_idx <- suppressWarnings(as.integer(acc_raw))
      if (any(is.na(acc_idx)) || any(acc_idx < 1L | acc_idx > length(acc_ids))) {
        stop("trial_df$accumulator numeric values must be 1..n_acc")
      }
      trial_df$accumulator <- acc_ids[acc_idx]
    } else {
      acc_ids_vec <- as.character(acc_raw)
      acc_ids_vec[!is.na(acc_ids_vec) & !nzchar(acc_ids_vec)] <- NA_character_
      if (any(is.na(acc_ids_vec))) stop("trial_df$accumulator must include valid accumulator ids")
      bad_vals <- unique(acc_ids_vec[!acc_ids_vec %in% acc_ids])
      if (length(bad_vals) > 0L) stop("trial_df accumulator values must match model accumulators: ",
                                      paste(bad_vals, collapse = ", "))
      trial_df$accumulator <- acc_ids_vec
    }
    pair_key <- paste(trial_df$trial, trial_df$accumulator, sep = "::")
    if (any(duplicated(pair_key))) {
      stop("trial_df with accumulator must have at most one row per trial/accumulator")
    }
  }

  trial_has_component <- FALSE
  if ("component" %in% names(trial_df)) {
    comp_col <- trial_df$component
    if (is.factor(comp_col)) comp_col <- as.character(comp_col)
    comp_col <- as.character(comp_col)
    comp_col[!is.na(comp_col) & !nzchar(comp_col)] <- NA_character_
    bad_vals <- unique(comp_col[!is.na(comp_col) & !comp_col %in% comp_ids])
    if (length(bad_vals) > 0L) stop("component values must match model components: ", paste(bad_vals, collapse = ", "))
    trial_df$component <- comp_col
    trial_has_component <- TRUE
  }

  trial_has_onset <- FALSE
  if ("onset" %in% names(trial_df)) {
    onset_col <- trial_df$onset
    if (is.factor(onset_col)) onset_col <- as.character(onset_col)
    if (!is.numeric(onset_col)) {
      coerced <- suppressWarnings(as.numeric(onset_col))
      if (any(!is.na(onset_col) & is.na(coerced))) stop("trial_df$onset must be numeric")
      onset_col <- coerced
    }
    trial_df$onset <- as.numeric(onset_col)
    trial_has_onset <- TRUE
  }

  list(
    trial_df = trial_df,
    trial_has_acc = trial_has_acc,
    trial_has_component = trial_has_component,
    trial_has_onset = trial_has_onset
  )
}

# ---- Public API ----------------------------------------------------------------

#' Simulate trials from parameter rows
#'
#' @param structure Model structure
#' @param params_df Parameter data frame
#' @param trial_df Optional trial-level or expanded data frame used to condition
#'   simulation. When `trial_df` includes an `accumulator` column, `onset` and
#'   `component` values are matched by trial and accumulator. If `layout = "long"`
#'   and `trial_df` rows align to `params_df`, those rows are treated as the
#'   accumulator mapping for each trial. Otherwise, values apply to all accumulators
#'   for a trial.
#' @param seed Optional RNG seed
#' @param keep_detail Whether to retain per-trial detail
#' @param keep_component Whether to include the chosen component in the output (if multiple components exist).
#'   If `NULL` (default), fixed mixtures keep the component; sampled mixtures drop it.
#' @param layout Parameter layout. `"rectangular"` expects `n_trials * n_acc` rows ordered
#'   by trial then accumulator; `"long"` expects either a `trial_df` with `accumulator`
#'   rows aligned to `params_df`, or a `trial_df` with `component` to infer active
#'   accumulators. `"auto"` chooses based on row counts and component structure.
#' @param ... Unused; for S3 compatibility
#' @return Data frame of simulated outcomes/RTs (with optional detail). When
#'   `race_spec(n_outcomes = k)` has `k > 1`, additional ordered readout columns
#'   are appended as `R2`/`rt2`, `R3`/`rt3`, ..., `Rk`/`rtk`.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params <- c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0)
#' df <- build_param_matrix(spec, params, n_trials = 3)
#' simulate(structure, df, seed = 123)
#' @export
simulate <- function(structure,
                     params_df,
                     trial_df = NULL,
                     seed = NULL,
                     keep_detail = FALSE,
                     keep_component = NULL,
                     layout = c("auto", "rectangular", "long"),
                     ...) {
  UseMethod("simulate")
}

#' @rdname simulate
#' @export
simulate.model_structure <- function(structure,
                                     params_df,
                                     trial_df = NULL,
                                     seed = NULL,
                                     keep_detail = FALSE,
                                     keep_component = NULL,
                                     layout = c("auto", "rectangular", "long"),
                                     ...) {
  layout <- match.arg(layout)

  if (is.null(params_df)) stop("Parameter matrix must be provided")
  structure <- .as_model_structure(structure)
  prep <- structure$prep
  acc_defs <- prep$accumulators %||% list()
  acc_ids <- names(acc_defs)
  n_acc <- length(acc_ids)
  if (n_acc == 0L) stop("No accumulators defined in model structure")

  params_mat <- if (is.matrix(params_df)) params_df else as.matrix(params_df)
  if (!is.numeric(params_mat)) stop("Parameter matrix must be numeric")
  if (nrow(params_mat) == 0L) stop("Parameter matrix must have rows")
  if (is.null(colnames(params_mat))) stop("Parameter matrix must have column names")
  required_cols <- c("q", "w", "t0")
  missing_cols <- setdiff(required_cols, colnames(params_mat))
  if (length(missing_cols) > 0L) stop("Parameter matrix missing columns: ", paste(missing_cols, collapse = ", "))
  p_cols <- grep("^p[0-9]+$", colnames(params_mat), value = TRUE)
  if (length(p_cols) == 0L) stop("Parameter matrix must include at least p1")
  p_cols <- p_cols[order(suppressWarnings(as.integer(sub("^p", "", p_cols))))]
  p_cols <- p_cols[seq_len(min(length(p_cols), 8L))]

  comp_table <- structure$components
  comp_ids <- comp_table$component_id %||% "__default__"
  mix_mode <- comp_table$mode[[1]] %||% "fixed"

  trial_info <- .normalize_trial_df(trial_df, acc_ids, comp_ids)
  trial_df <- trial_info$trial_df
  trial_has_acc <- trial_info$trial_has_acc
  trial_has_component <- trial_info$trial_has_component
  trial_has_onset <- trial_info$trial_has_onset
  trial_ids <- if (!is.null(trial_df)) sort(unique(trial_df$trial)) else integer(0)

  layout_mode <- layout
  if (identical(layout, "auto")) {
    if (nrow(params_mat) %% n_acc != 0L) {
      layout_mode <- "long"
    } else if (trial_has_component) {
      comp_levels <- unique(trial_df$component)
      comp_levels <- comp_levels[!is.na(comp_levels)]
      if (length(comp_levels) > 0L) {
        acc_counts <- vapply(comp_levels, function(comp_label) {
          sum(vapply(acc_defs, function(acc_def) {
            .acc_active_in_component(acc_def, comp_label)
          }, logical(1)))
        }, integer(1))
        if (any(acc_counts != n_acc)) layout_mode <- "long"
      }
    }
    if (identical(layout_mode, "auto")) layout_mode <- "rectangular"
  }

  mapping_mode <- "implicit"
  if (identical(layout_mode, "long") && identical(layout, "long") &&
      trial_has_acc && !is.null(trial_df) && nrow(params_mat) == nrow(trial_df)) {
    mapping_mode <- "explicit"
  }

  if (identical(layout_mode, "rectangular")) {
    n_trials <- nrow(params_mat) / n_acc
    if (!is.finite(n_trials) || n_trials != floor(n_trials)) {
      stop(sprintf("Parameter rows (%d) not divisible by number of accumulators (%d); ",
                   nrow(params_mat), n_acc),
           "use layout = \"long\" for non-rectangular matrices")
    }
    n_trials <- as.integer(n_trials)
    if (!is.null(trial_df)) {
      if (any(trial_df$trial > n_trials)) {
        stop("trial_df$trial values must be between 1 and n_trials")
      }
    }
  } else if (identical(mapping_mode, "explicit")) {
    if (is.null(trial_df) || !trial_has_acc) {
      stop("layout = \"long\" with explicit mapping requires trial_df with accumulator")
    }
    if (nrow(params_mat) != nrow(trial_df)) stop("For explicit long layout, params_df rows must match trial_df rows")
    if (length(trial_ids) == 0L) stop("trial_df must include trial rows")
    n_trials <- max(trial_ids)
  } else {
    if (is.null(trial_df) || !trial_has_component) stop("Non-rectangular layout requires trial_df with component when accumulator mapping is absent")
    n_trials <- length(trial_ids)
    if (n_trials == 0L) stop("trial_df must include trial rows")
  }

  component_vec <- NULL
  if (!is.null(trial_df) && trial_has_component) {
    component_vec <- rep(NA_character_, n_trials)
    comp_by_trial <- split(trial_df$component, trial_df$trial)
    for (tid in names(comp_by_trial)) {
      vals <- unique(comp_by_trial[[tid]])
      vals <- vals[!is.na(vals)]
      if (length(vals) > 1L) stop(sprintf("Multiple component values for trial %s", tid))
      if (length(vals) == 1L) component_vec[[as.integer(tid)]] <- vals[[1]]
    }
    if (identical(layout_mode, "long") && identical(mapping_mode, "implicit") &&
        any(is.na(component_vec))) {
      stop("component values must be provided for each trial when accumulator mapping is absent")
    }
  }

  comp_active_map <- NULL
  if (identical(layout_mode, "long") && identical(mapping_mode, "implicit")) {
    comp_levels <- unique(component_vec)
    comp_active_map <- lapply(comp_levels, function(comp_label) {
      vapply(acc_defs, function(acc_def) {
        .acc_active_in_component(acc_def, comp_label)
      }, logical(1))
    })
    names(comp_active_map) <- comp_levels
  }

  onset_mat <- NULL
  if (!is.null(trial_df) && trial_has_onset) {
    onset_mat <- matrix(NA_real_, nrow = n_trials, ncol = n_acc, dimnames = list(NULL, acc_ids))
    if (trial_has_acc) {
      acc_idx <- match(trial_df$accumulator, acc_ids)
      for (i in seq_len(nrow(trial_df))) {
        t <- trial_df$trial[[i]]
        a <- acc_idx[[i]]
        onset_val <- trial_df$onset[[i]]
        if (!is.finite(onset_val)) next
        onset_mat[t, a] <- onset_val
      }
    } else {
      for (i in seq_len(nrow(trial_df))) {
        t <- trial_df$trial[[i]]
        onset_val <- trial_df$onset[[i]]
        if (!is.finite(onset_val)) next
        onset_mat[t, ] <- onset_val
      }
    }
  }

  comp_leader_acc <- setNames(rep(NA_character_, length(comp_ids)), comp_ids)
  for (acc_id in acc_ids) {
    comps <- acc_defs[[acc_id]]$components %||% character(0)
    for (c_id in comps) {
      if (is.na(comp_leader_acc[[c_id]])) comp_leader_acc[[c_id]] <- acc_id
    }
  }

  row_map <- vector("list", n_trials)
  if (identical(layout_mode, "rectangular")) {
    for (i in seq_len(n_trials)) {
      start <- (i - 1L) * n_acc
      row_map[[i]] <- setNames(start + seq_len(n_acc), acc_ids)
    }
  } else if (identical(mapping_mode, "explicit")) {
    rows_by_trial <- split(seq_len(nrow(trial_df)), trial_df$trial)
    for (i in seq_len(n_trials)) {
      rows <- rows_by_trial[[as.character(i)]]
      if (is.null(rows) || length(rows) == 0L) stop(sprintf("No parameter rows mapped to trial %d", i))
      accs <- trial_df$accumulator[rows]
      mapped <- setNames(rows, accs)
      mapped <- mapped[acc_ids[acc_ids %in% names(mapped)]]
      row_map[[i]] <- mapped
    }
  } else {
    row_cursor <- 1L
    for (i in seq_len(n_trials)) {
      active_mask <- comp_active_map[[component_vec[[i]]]]
      accs <- acc_ids[active_mask]
      if (length(accs) == 0L) stop("No accumulator rows matched the provided component labels")
      rows <- row_cursor + seq_len(length(accs)) - 1L
      row_map[[i]] <- setNames(rows, accs)
      row_cursor <- row_cursor + length(accs)
    }
    if (row_cursor != (nrow(params_mat) + 1L)) {
      stop("Parameter rows did not align with component assignments")
    }
  }

  if (!is.null(seed)) set.seed(seed)
  trial_ids <- seq_len(n_trials)
  n_readout <- prep$observation$n_outcomes %||% 1L
  n_readout <- as.integer(n_readout)
  if (is.na(n_readout) || n_readout < 1L) {
    n_readout <- 1L
  }
  outcomes <- matrix(NA_character_, nrow = length(trial_ids), ncol = n_readout)
  rts <- matrix(NA_real_, nrow = length(trial_ids), ncol = n_readout)
  comp_record <- rep(NA_character_, length(trial_ids))
  details <- if (keep_detail) vector("list", length(trial_ids)) else NULL

  acc_idx_map <- setNames(seq_along(acc_ids), acc_ids)

  for (i in seq_along(trial_ids)) {
    forced_component <- if (!is.null(component_vec)) component_vec[[i]] else NA_character_
    if (!is.na(forced_component)) {
      chosen_component <- forced_component
    } else if (length(comp_ids) == 1L) {
      chosen_component <- comp_ids[[1]]
    } else {
      comp_weights <- numeric(length(comp_ids))
      trial_rows <- row_map[[i]]
      for (ci in seq_along(comp_ids)) {
        leader_acc <- comp_leader_acc[[comp_ids[[ci]]]]
        if (!is.na(leader_acc) && !is.null(trial_rows[[leader_acc]])) {
          row_idx <- trial_rows[[leader_acc]]
          w_val <- params_mat[row_idx, "w"]
          comp_weights[[ci]] <- if (is.finite(w_val) && w_val >= 0) w_val else 0
        } else {
          comp_weights[[ci]] <- 1
        }
      }
      if (sum(comp_weights) <= 0 || any(!is.finite(comp_weights))) {
        comp_weights <- rep(1 / length(comp_weights), length(comp_weights))
      } else {
        comp_weights <- comp_weights / sum(comp_weights)
      }
      chosen_component <- sample(comp_ids, size = 1L, prob = comp_weights)
    }
    comp_record[[i]] <- chosen_component

    trial_accs <- list()
    shared_map <- list()
    trial_rows <- row_map[[i]]
    if (length(trial_rows) == 0L) stop(sprintf("No parameter rows mapped to trial %d", i))

    # preserve accumulator order where possible
    trial_rows <- trial_rows[acc_ids[acc_ids %in% names(trial_rows)]]

    for (acc_id in names(trial_rows)) {
      row_idx <- trial_rows[[acc_id]]
      row_vals <- params_mat[row_idx, ]
      dist_params <- dist_param_names(acc_defs[[acc_id]]$dist)
      if (length(dist_params) > length(p_cols)) {
        stop(
          sprintf(
            "Accumulator '%s' requires %d distribution params, but matrix has only %d p-slots",
            acc_id,
            length(dist_params),
            length(p_cols)
          )
        )
      }
      dist_list <- setNames(as.list(row_vals[p_cols][seq_along(dist_params)]), dist_params)
      dist_list$t0 <- row_vals[["t0"]]
      onset_override <- NA_real_
      if (!is.null(onset_mat)) {
        onset_override <- as.numeric(onset_mat[i, acc_idx_map[[acc_id]]])
      }
      onset_spec <- acc_defs[[acc_id]]$onset_spec %||% list(
        kind = "absolute",
        value = acc_defs[[acc_id]]$onset %||% 0
      )
      onset_kind <- onset_spec$kind %||% "absolute"
      onset_val <- acc_defs[[acc_id]]$onset %||% 0
      if (identical(onset_kind, "absolute")) {
        onset_val <- as.numeric(onset_spec$value %||% onset_val)[1]
        if (is.finite(onset_override)) {
          onset_val <- onset_override
        }
        onset_spec <- list(kind = "absolute", value = onset_val)
      } else if (identical(onset_kind, "after")) {
        lag_val <- as.numeric(onset_spec$lag %||% 0)[1]
        if (is.finite(onset_override)) {
          lag_val <- lag_val + onset_override
        }
        onset_spec <- list(
          kind = "after",
          source = onset_spec$source,
          source_kind = onset_spec$source_kind,
          lag = lag_val
        )
        onset_val <- 0
      } else {
        stop(sprintf("Unsupported onset kind '%s' for accumulator '%s'", onset_kind, acc_id))
      }
      acc_entry <- list(
        dist = acc_defs[[acc_id]]$dist,
        params = dist_list,
        onset = onset_val,
        onset_spec = onset_spec,
        q = row_vals[["q"]] %||% 0,
        shared_trigger_id = acc_defs[[acc_id]]$shared_trigger_id %||% NA_character_,
        components = acc_defs[[acc_id]]$components %||% character(0)
      )
      trial_accs[[acc_id]] <- acc_entry
      stid <- acc_entry$shared_trigger_id
      if (!is.null(stid) && !is.na(stid) && nzchar(stid)) {
        if (is.null(shared_map[[stid]])) {
          shared_map[[stid]] <- list(prob = row_vals[["q"]] %||% 0)
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
    out_vals <- result$outcomes %||% result$outcome
    rt_vals <- result$rts %||% result$rt
    if (length(out_vals) > n_readout) out_vals <- out_vals[seq_len(n_readout)]
    if (length(rt_vals) > n_readout) rt_vals <- rt_vals[seq_len(n_readout)]
    if (length(out_vals) > 0L) outcomes[i, seq_along(out_vals)] <- as.character(out_vals)
    if (length(rt_vals) > 0L) rts[i, seq_along(rt_vals)] <- as.numeric(rt_vals)
    if (keep_detail) details[[i]] <- result$detail
  }

  out_df <- data.frame(
    trial = trial_ids,
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

#' @rdname simulate
#' @export
simulate.default <- function(structure, ...) {
  stop("simulate() expects a model_structure", call. = FALSE)
}
