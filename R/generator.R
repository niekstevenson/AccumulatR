# ---- Component activity helpers ----------------------------------------------

.acc_active_in_component <- function(acc_def, component) {
  comps <- acc_def$components
  if (length(comps) == 0 || is.null(component) || identical(component, "__default__")) return(TRUE)
  component %in% comps
}

.outcome_allowed_in_component <- function(options, component) {
  comps <- options$component %||% NULL
  if (is.null(comps) || length(comps) == 0L) return(TRUE)
  comp_label <- component %||% "__default__"
  if (is.na(comp_label) || !nzchar(comp_label)) return(FALSE)
  comp_label %in% as.character(comps)
}

.component_readout_count <- function(prep, component) {
  obs <- prep$observation
  base_readout <- obs$global_n_outcomes
  comp_label <- component %||% "__default__"
  comp_override <- obs$component_n_outcomes[[comp_label]]
  if (!is.null(comp_override)) {
    base_readout <- comp_override[[1L]]
  }
  min(base_readout, obs$n_outcomes)
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
  shared_id <- acc_def$shared_trigger_id %||% NULL
  if (!is.null(shared_id) && !is.na(shared_id) && shared_id != "") {
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
  trial_has_acc <- "racer" %in% names(trial_df)
  if (!"trials" %in% names(trial_df)) {
    if (trial_has_acc) stop("trial_df with racer must include a trials column")
    trial_df$trials <- seq_len(nrow(trial_df))
  }
  if (is.factor(trial_df$trials)) trial_df$trials <- as.character(trial_df$trials)
  trial_df$trials <- suppressWarnings(as.integer(trial_df$trials))
  if (any(is.na(trial_df$trials)) || any(trial_df$trials < 1L)) {
    stop("trial_df$trials must contain positive integers")
  }

  if (trial_has_acc) {
    acc_raw <- trial_df$racer
    if (is.factor(acc_raw)) acc_raw <- as.character(acc_raw)
    if (is.numeric(acc_raw)) {
      acc_idx <- suppressWarnings(as.integer(acc_raw))
      if (any(is.na(acc_idx)) || any(acc_idx < 1L | acc_idx > length(acc_ids))) {
        stop("trial_df$racer numeric values must be 1..n_acc")
      }
      trial_df$racer <- acc_ids[acc_idx]
    } else {
      acc_ids_vec <- as.character(acc_raw)
      acc_ids_vec[!is.na(acc_ids_vec) & !nzchar(acc_ids_vec)] <- NA_character_
      if (any(is.na(acc_ids_vec))) stop("trial_df$racer must include valid racer ids")
      bad_vals <- unique(acc_ids_vec[!acc_ids_vec %in% acc_ids])
      if (length(bad_vals) > 0L) stop("trial_df racer values must match model racers: ",
                                      paste(bad_vals, collapse = ", "))
      trial_df$racer <- acc_ids_vec
    }
    pair_key <- paste(trial_df$trials, trial_df$racer, sep = "::")
    if (any(duplicated(pair_key))) {
      stop("trial_df with racer must have at most one row per trials/racer")
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

#' Simulate behavioral data from a model
#'
#' @param structure Finalized model structure.
#' @param params_df Trial-level parameter values.
#' @param trial_df Optional data frame used to condition the simulation. When it
#'   includes a `racer` column, `onset` and `component` values are matched
#'   by `trials` and racer. Otherwise, values apply at the trial level.
#' @param seed Optional random-number seed.
#' @param keep_detail If `TRUE`, keep additional simulation detail.
#' @param keep_component Whether to keep the chosen mixture component in the
#'   output when the model has multiple components. If `NULL`, fixed mixtures
#'   keep the component label and sampled mixtures drop it.
#' @param layout Parameter layout. `"rectangular"` expects rows ordered by trial
#'   and racer. `"long"` expects rows aligned to `trial_df` or inferable
#'   from `component`. `"auto"` chooses automatically.
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
                     keep_component = NULL,
                     layout = c("auto", "rectangular", "long")) {
  layout <- match.arg(layout)

  if (is.null(params_df)) stop("Parameter matrix must be provided")
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
  trial_ids <- if (!is.null(trial_df)) sort(unique(trial_df$trials)) else integer(0)

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
      if (any(trial_df$trials > n_trials)) {
        stop("trial_df$trials values must be between 1 and n_trials")
      }
    }
  } else if (identical(mapping_mode, "explicit")) {
    if (is.null(trial_df) || !trial_has_acc) {
      stop("layout = \"long\" with explicit mapping requires trial_df with racer")
    }
    if (nrow(params_mat) != nrow(trial_df)) stop("For explicit long layout, params_df rows must match trial_df rows")
    if (length(trial_ids) == 0L) stop("trial_df must include trial rows")
    n_trials <- max(trial_ids)
  } else {
    if (is.null(trial_df) || !trial_has_component) stop("Non-rectangular layout requires trial_df with component when racer mapping is absent")
    n_trials <- length(trial_ids)
    if (n_trials == 0L) stop("trial_df must include trial rows")
  }

  component_vec <- NULL
  if (!is.null(trial_df) && trial_has_component) {
    component_vec <- rep(NA_character_, n_trials)
    comp_by_trial <- split(trial_df$component, trial_df$trials)
    for (tid in names(comp_by_trial)) {
      vals <- unique(comp_by_trial[[tid]])
      vals <- vals[!is.na(vals)]
      if (length(vals) > 1L) stop(sprintf("Multiple component values for trial %s", tid))
      if (length(vals) == 1L) component_vec[[as.integer(tid)]] <- vals[[1]]
    }
    if (identical(layout_mode, "long") && identical(mapping_mode, "implicit") &&
        any(is.na(component_vec))) {
      stop("component values must be provided for each trial when racer mapping is absent")
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
      acc_idx <- match(trial_df$racer, acc_ids)
      for (i in seq_len(nrow(trial_df))) {
        t <- trial_df$trials[[i]]
        a <- acc_idx[[i]]
        onset_val <- trial_df$onset[[i]]
        if (!is.finite(onset_val)) next
        onset_mat[t, a] <- onset_val
      }
    } else {
      for (i in seq_len(nrow(trial_df))) {
        t <- trial_df$trials[[i]]
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
    rows_by_trial <- split(seq_len(nrow(trial_df)), trial_df$trials)
    for (i in seq_len(n_trials)) {
      rows <- rows_by_trial[[as.character(i)]]
      if (is.null(rows) || length(rows) == 0L) stop(sprintf("No parameter rows mapped to trial %d", i))
      accs <- trial_df$racer[rows]
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
  n_readout <- prep$observation$n_outcomes
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
        row_idx <- trial_rows[[leader_acc]]
        comp_weights[[ci]] <- params_mat[row_idx, "w"]
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
      dist_list <- setNames(as.list(row_vals[p_cols][seq_along(dist_params)]), dist_params)
      dist_list$t0 <- row_vals[["t0"]]
      onset_override <- NA_real_
      if (!is.null(onset_mat)) {
        onset_override <- as.numeric(onset_mat[i, acc_idx_map[[acc_id]]])
      }
      onset_spec <- acc_defs[[acc_id]]$onset_spec
      onset_val <- acc_defs[[acc_id]]$onset
      if (identical(onset_spec$kind, "absolute")) {
        onset_val <- onset_spec$value
        if (is.finite(onset_override)) {
          onset_val <- onset_override
        }
        onset_spec <- list(kind = "absolute", value = onset_val)
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
        onset_val <- 0
      }
      acc_entry <- list(
        dist = acc_defs[[acc_id]]$dist,
        params = dist_list,
        onset = onset_val,
        onset_spec = onset_spec,
        q = row_vals[["q"]],
        shared_trigger_id = acc_defs[[acc_id]]$shared_trigger_id %||% NA_character_,
        components = acc_defs[[acc_id]]$components %||% character(0)
      )
      trial_accs[[acc_id]] <- acc_entry
      stid <- acc_entry$shared_trigger_id
      if (!is.null(stid) && !is.na(stid) && nzchar(stid)) {
        if (is.null(shared_map[[stid]])) {
          shared_map[[stid]] <- list(prob = row_vals[["q"]])
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
