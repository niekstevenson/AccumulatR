# Interpreter-driven race simulator for the boolean-expression DSL defined in new_API.R
# Depends on core distribution and sampling helpers from the existing codebase.

source("R/dist.R")
source("R/utils.R")
source("R/model_tables.R")
`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

# ---- Model preparation helpers ------------------------------------------------

.normalize_model <- function(model) {
  if (inherits(model, "race_spec")) {
    if (exists("build_model", mode = "function")) {
      return(build_model(model))
    }
    return(structure(list(
      accumulators = unname(model$accumulators),
      pools = unname(model$pools),
      outcomes = unname(model$outcomes),
      groups = unname(model$groups),
      metadata = model$metadata %||% list()
    ), class = "race_model_spec"))
  }
  if (inherits(model, "race_model_spec")) return(model)
  model
}

.prepare_acc_defs <- function(model) {
  accs <- model$accumulators
  if (length(accs) == 0) stop("Model must define at least one accumulator")

  # Initialize accumulator definitions keyed by id
  defs <- setNames(vector("list", length(accs)), vapply(accs, `[[`, character(1), "id"))
  for (acc in accs) {
    defs[[acc$id]] <- list(
      id = acc$id,
      dist = acc$dist,
      onset = acc$onset %||% 0,
      q = acc$q %||% 0,
      params = acc$params %||% list(),
      components = character(0)
    )
  }

  # Apply group-level attributes (shared params, components)
  if (!is.null(model$groups) && length(model$groups) > 0) {
    for (grp in model$groups) {
      members <- grp$members
      attrs <- grp$attrs
      if (is.null(members) || length(members) == 0) next
      if (!is.null(attrs$shared_params)) {
        for (m in members) {
          if (!is.null(defs[[m]])) {
            defs[[m]]$params <- modifyList(attrs$shared_params, defs[[m]]$params, keep.null = TRUE)
          }
        }
      }
      if (!is.null(attrs$component)) {
        comp_id <- attrs$component
        for (m in members) {
          if (!is.null(defs[[m]])) {
            defs[[m]]$components <- unique(c(defs[[m]]$components, comp_id))
          }
        }
      }
    }
  }

  # Final sanity: ensure each accumulator has params list ready
  for (acc_id in names(defs)) {
    if (is.null(defs[[acc_id]]$params)) defs[[acc_id]]$params <- list()
  }
  defs
}

.extract_components <- function(model) {
  mixture <- model$metadata$mixture
  if (is.null(mixture) || is.null(mixture$components) || length(mixture$components) == 0) {
    list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list(deadline = model$metadata$deadline %||% Inf, guess = NULL)),
      has_weight_param = FALSE
    )
  } else {
    comps <- mixture$components
    ids <- vapply(comps, `[[`, character(1), "id")
    weights <- vapply(comps, function(cmp) cmp$weight %||% 1, numeric(1))
    if (all(is.na(weights))) weights <- rep(1, length(weights))
    weights <- weights / sum(weights)
    attrs <- setNames(vector("list", length(ids)), ids)
    has_wparam <- logical(length(ids))
    for (i in seq_along(ids)) {
      wp <- comps[[i]]$attrs$weight_param %||% NULL
      has_wparam[[i]] <- !is.null(wp)
      attrs[[ids[[i]]]] <- list(
        deadline = comps[[i]]$attrs$deadline %||% model$metadata$deadline %||% Inf,
        guess = comps[[i]]$attrs$guess %||% NULL,
        weight_param = wp
      )
    }
    list(ids = ids, weights = weights, attrs = attrs, has_weight_param = has_wparam)
  }
}

.prepare_pool_defs <- function(model) {
  pools <- model$pools
  if (is.null(pools)) pools <- list()
  defs <- setNames(vector("list", length(pools)), vapply(pools, `[[`, character(1), "id"))
  for (pl in pools) {
    rule <- pl$rule
    k <- rule$k
    if (is.infinite(k)) k <- length(pl$members)
    defs[[pl$id]] <- list(
      id = pl$id,
      members = pl$members,
      k = k,
      weights = rule$weights %||% NULL,
      tags = pl$tags %||% list()
    )
  }
  defs
}

.prepare_outcomes <- function(model) {
  outs <- model$outcomes
  if (is.null(outs) || length(outs) == 0) stop("Model must define outcomes")
  outs <- lapply(outs, function(out) {
    if (is.null(out$label)) stop("Outcome missing label")
    if (is.null(out$expr)) stop(sprintf("Outcome '%s' missing expr", out$label))
    out$options <- out$options %||% list()
    out
  })
  setNames(outs, vapply(outs, `[[`, character(1), "label"))
}

prepare_model <- function(model) {
  model <- .normalize_model(model)
  acc_defs <- .prepare_acc_defs(model)
  pool_defs <- .prepare_pool_defs(model)
  outcome_defs <- .prepare_outcomes(model)
  component_defs <- .extract_components(model)
  list(
    accumulators = acc_defs,
    pools = pool_defs,
    outcomes = outcome_defs,
    components = component_defs,
    default_deadline = model$metadata$deadline %||% Inf,
    special_outcomes = model$metadata$special_outcomes %||% list()
  )
}

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

# ---- Sampling primitives ------------------------------------------------------

.sample_accumulator <- function(acc_def) {
  success <- stats::runif(1) < (1 - acc_def$q)
  if (!success) return(Inf)
  reg <- dist_registry(acc_def$dist)
  if (is.null(reg) || is.null(reg$r)) stop(sprintf("No sampler registered for distribution '%s'", acc_def$dist))
  draw <- reg$r(1L, acc_def$params)
  acc_def$onset + as.numeric(draw[[1]])
}

# ---- Expression evaluation ----------------------------------------------------

.get_acc_time <- function(ctx, acc_id) {
  acc_defs <- ctx$model[["accumulators"]]
  if (!acc_id %in% names(acc_defs)) stop(sprintf("Unknown accumulator '%s'", acc_id))
  if (is.null(ctx$acc_times[[acc_id]])) {
    ctx$acc_times[[acc_id]] <- .sample_accumulator(acc_defs[[acc_id]])
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
    deadline <- .get_component_attr(ctx$model, ctx$component, "deadline")
    if (is.null(deadline) || !is.finite(deadline)) deadline <- ctx$model[["default_deadline"]]
    result <- list(time = deadline %||% Inf, core = numeric(0))
  } else if (identical(source_id, "__GUESS__")) {
    deadline <- .get_component_attr(ctx$model, ctx$component, "deadline")
    result <- list(time = deadline %||% Inf, core = numeric(0))
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
  for (label in labels) {
    expr <- outs[[label]][["expr"]]
    eval <- .eval_expr(expr, ctx)
    res[[label]] <- list(
      time = eval$time,
      core = eval$core,
      options = outs[[label]][["options"]] %||% list()
    )
  }
  res
}

# ---- Trial simulation ---------------------------------------------------------

.simulate_trial <- function(prep, component, keep_detail = FALSE) {
  comp_info <- prep[["components"]]
  component_attr <- comp_info$attrs[[component]]
  if (is.null(component_attr)) component_attr <- list()
  ctx <- list(
    model = prep,
    component = component,
    acc_times = new.env(parent = emptyenv()),
    pool_cache = new.env(parent = emptyenv()),
    event_cache = new.env(parent = emptyenv())
  )

  outcomes <- .evaluate_outcomes(ctx)
  cand_labels <- character(0)
  cand_times <- numeric(0)
  cand_core <- list()
  cand_options <- list()

  for (lbl in names(outcomes)) {
    time <- outcomes[[lbl]]$time
    if (!is.finite(time)) next
    cand_labels <- c(cand_labels, lbl)
    cand_times <- c(cand_times, time)
    cand_core[[lbl]] <- outcomes[[lbl]]$core
    cand_options[[lbl]] <- outcomes[[lbl]]$options
  }

  # Include global deadline if not already represented
  default_deadline <- component_attr$deadline %||% prep[["default_deadline"]]
  if (is.finite(default_deadline) && !("NR_DEADLINE" %in% cand_labels)) {
    special_deadline_label <- prep[["special_outcomes"]]$deadline %||% "NR_DEADLINE"
    cand_labels <- c(cand_labels, special_deadline_label)
    cand_times <- c(cand_times, default_deadline)
    cand_core[[tail(cand_labels, 1)]] <- numeric(0)
    cand_options[[tail(cand_labels, 1)]] <- list()
  }

  if (length(cand_times) == 0) {
    result <- list(outcome = NA_character_, rt = NA_real_, detail = NULL)
    return(result)
  }

  tmin <- min(cand_times)
  tied <- which(cand_times == tmin)
  chosen_idx <- tied[[1]]
  if (length(tied) > 1) {
    tie_scores <- vapply(tied, function(i) {
      core <- cand_core[[cand_labels[[i]]]]
      if (length(core) == 0) cand_times[[i]] else min(core)
    }, numeric(1))
    chosen_idx <- tied[which.min(tie_scores)]
  }

  chosen_label <- cand_labels[[chosen_idx]]
  chosen_time <- cand_times[[chosen_idx]]
  options <- cand_options[[chosen_label]]

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

  list(outcome = chosen_label, rt = chosen_time, detail = detail)
}

# ---- Public API ----------------------------------------------------------------

simulate_model <- function(model, n_trials, seed = NULL, keep_detail = FALSE) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  if (!is.null(seed)) set.seed(seed)
  prep <- prepare_model(model)
  components <- prep[["components"]]
  comp_ids <- components$ids
  weights <- components$weights

  trial_components <- if (length(comp_ids) == 1L) {
    rep(comp_ids, n_trials)
  } else {
    sample(comp_ids, size = n_trials, replace = TRUE, prob = weights)
  }

  outcomes <- character(n_trials)
  rts <- numeric(n_trials)
  details <- if (keep_detail) vector("list", n_trials) else NULL

  for (i in seq_len(n_trials)) {
    comp <- trial_components[[i]]
    res <- .simulate_trial(prep, comp, keep_detail = keep_detail)
    outcomes[[i]] <- res$outcome
    rts[[i]] <- res$rt
    if (keep_detail) details[[i]] <- res$detail
  }

  out_df <- data.frame(
    trial = seq_len(n_trials),
    outcome = outcomes,
    rt = rts,
    stringsAsFactors = FALSE
  )
  if (length(comp_ids) > 1L) {
    out_df$component <- trial_components
  }
  if (keep_detail) attr(out_df, "details") <- details
  out_df
}
