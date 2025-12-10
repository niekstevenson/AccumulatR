# Interpreter-driven race simulator for the boolean-expression DSL defined in new_API.R
# Depends on core distribution and sampling helpers from the existing codebase.
`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.ensure_acc_param_t0 <- function(params) {
  if (is.null(params) || length(params) == 0L) {
    params <- list()
  }
  if (is.null(params$t0) || length(params$t0) == 0L) {
    params$t0 <- 0
  } else {
    params$t0 <- as.numeric(params$t0)[1]
  }
  params
}

# ---- Model preparation helpers ------------------------------------------------

.normalize_model <- function(model) {
  if (inherits(model, "race_spec")) {
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
      params = .ensure_acc_param_t0(acc$params %||% list()),
      components = character(0),
      shared_trigger_id = NULL,
      shared_trigger_q = NULL
    )
  }

  shared_triggers <- list()

  # Apply group-level attributes (shared params, components)
  if (!is.null(model$groups) && length(model$groups) > 0) {
    for (grp in model$groups) {
      members <- grp$members
      attrs <- grp$attrs
      if (is.null(members) || length(members) == 0) next
      if (!is.null(attrs$shared_params)) {
        for (m in members) {
          if (!is.null(defs[[m]])) {
            # Shared parameters should override member defaults so they stay linked
            defs[[m]]$params <- modifyList(defs[[m]]$params, attrs$shared_params, keep.null = TRUE)
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
      shared_trigger <- attrs$shared_trigger %||% NULL
      if (!is.null(shared_trigger)) {
        trig_id <- shared_trigger$id %||% grp$id
        if (is.null(trig_id) || trig_id == "") {
          stop(sprintf("Group '%s' shared_trigger must provide a non-empty id", grp$id %||% "<unnamed>"))
        }
        default_q <- shared_trigger$q %||% NA_real_
        if (is.na(default_q)) {
          q_vals <- vapply(members, function(m) {
            acc_def <- defs[[m]]
            if (is.null(acc_def)) return(NA_real_)
            acc_def$q %||% NA_real_
          }, numeric(1))
          q_vals <- q_vals[!is.na(q_vals)]
          if (length(q_vals) == 0) {
            default_q <- 0
          } else if (length(unique(q_vals)) == 1) {
            default_q <- q_vals[[1]]
          } else {
            stop(sprintf(
              "Group '%s' shared_trigger requires a single q value; found %s",
              trig_id,
              paste(format(unique(q_vals)), collapse = ", ")
            ))
          }
        }
        param_name <- shared_trigger$param %||% NULL
        shared_triggers[[trig_id]] <- list(
          id = trig_id,
          group_id = grp$id %||% trig_id,
          members = members,
          q = as.numeric(default_q),
          param = param_name
        )
        for (m in members) {
          if (!is.null(defs[[m]])) {
            defs[[m]]$shared_trigger_id <- trig_id
            defs[[m]]$shared_trigger_q <- as.numeric(default_q)
            # Keep q carrying the gate probability for downstream likelihood/native use.
            defs[[m]]$q <- as.numeric(default_q)
          }
        }
      }
    }
  }

  # Final sanity: ensure each accumulator has params list ready
  for (acc_id in names(defs)) {
    defs[[acc_id]]$params <- .ensure_acc_param_t0(defs[[acc_id]]$params)
  }
  list(acc = defs, shared_triggers = shared_triggers)
}

.extract_components <- function(model) {
  mixture <- model$metadata$mixture
  if (is.null(mixture) || is.null(mixture$components) || length(mixture$components) == 0) {
    list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list(guess = NULL)),
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
  acc_prep <- .prepare_acc_defs(model)
  acc_defs <- acc_prep$acc
  pool_defs <- .prepare_pool_defs(model)
  outcome_defs <- .prepare_outcomes(model)
  component_defs <- .extract_components(model)
  list(
    accumulators = acc_defs,
    pools = pool_defs,
    outcomes = outcome_defs,
    components = component_defs,
    special_outcomes = model$metadata$special_outcomes %||% list(),
    shared_triggers = acc_prep$shared_triggers %||% list()
  )
}

# ---- Structure + parameter helpers ------------------------------------------

.build_component_table <- function(comp_defs) {
  comp_ids <- comp_defs$ids
  comp_tbl <- data.frame(
    component_id = comp_ids,
    weight = comp_defs$weights,
    stringsAsFactors = FALSE
  )
  comp_tbl$has_weight_param <- comp_defs$has_weight_param
  comp_tbl$attrs <- I(comp_defs$attrs[match(comp_ids, names(comp_defs$attrs))])
  comp_tbl
}

.build_accumulator_template <- function(acc_defs) {
  acc_ids <- names(acc_defs)
  if (length(acc_ids) == 0L) {
    return(data.frame(
      accumulator_id = character(0),
      accumulator_index = integer(0),
      dist = character(0),
      onset = numeric(0),
      q = numeric(0),
      role = character(0),
      shared_trigger_id = character(0),
      shared_trigger_q = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  acc_df <- data.frame(
    accumulator_id = acc_ids,
    accumulator_index = seq_along(acc_ids),
    dist = vapply(acc_defs, function(acc) acc$dist %||% NA_character_, character(1)),
    onset = vapply(acc_defs, function(acc) acc$onset %||% 0, numeric(1)),
    q = vapply(acc_defs, function(acc) acc$q %||% 0, numeric(1)),
    role = rep("std", length(acc_ids)),
    shared_trigger_id = vapply(acc_defs, function(acc) acc$shared_trigger_id %||% NA_character_, character(1)),
    shared_trigger_q = vapply(acc_defs, function(acc) acc$shared_trigger_q %||% NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )
  acc_df$params <- I(lapply(acc_defs, function(acc) .ensure_acc_param_t0(acc$params %||% list())))
  acc_df$components <- I(lapply(acc_defs, function(acc) acc$components %||% character(0)))
  acc_df
}

#' Finalize a model for simulation/likelihood
#'
#' @param model Model specification
#' @return A model_structure list
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' finalize_model(spec)
#' @export
finalize_model <- function(model) {
  unwrap_model_spec <- function(x) {
    seen <- 0L
    repeat {
      if (inherits(x, "generator_structure")) {
        x <- x$model_spec %||% stop("generator_structure missing model_spec")
      } else if (!inherits(x, "race_model_spec") && !is.null(x$model_spec)) {
        x <- x$model_spec
      } else {
        break
      }
      seen <- seen + 1L
      if (seen > 20L) stop("generator_structure nesting too deep")
    }
    x
  }

  # Unwrap any generator_structures and normalize to a plain race_model_spec
  model_norm <- unwrap_model_spec(model)
  model_norm <- .normalize_model(model_norm)
  model_norm <- unwrap_model_spec(model_norm)
  if (!inherits(model_norm, "race_model_spec")) {
    stop("finalize_model requires a race_model_spec")
  }

  # Deep copy to prevent accidental aliasing back into the enclosing structure.
  model_norm <- unserialize(serialize(model_norm, NULL))
  prep <- prepare_model(model_norm)
  structure <- list(
    model_spec = model_norm,
    prep = prep,
    accumulators = .build_accumulator_template(prep$accumulators),
    components = .build_component_table(prep$components),
    special_outcomes = prep$special_outcomes,
    shared_triggers = prep$shared_triggers %||% list()
  )
  class(structure) <- c("model_structure", "generator_structure", class(structure))
  structure
}

.as_model_structure <- function(x) {
  if (inherits(x, "generator_structure")) return(x)
  if (inherits(x, "model_structure")) return(x)
  if (is.list(x) && !is.null(x$prep) && !is.null(x$accumulators) && !is.null(x$components)) {
    class(x) <- unique(c("model_structure", "generator_structure", class(x)))
    return(x)
  }
  finalize_model(x)
}

# Backward-compatible alias for older code paths
.as_generator_structure <- .as_model_structure

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

.sample_accumulator <- function(acc_def, ctx) {
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
  acc_def$onset + as.numeric(draw[[1]])
}

# ---- Expression evaluation ----------------------------------------------------

.get_acc_time <- function(ctx, acc_id) {
  acc_defs <- ctx$model[["accumulators"]]
  if (!acc_id %in% names(acc_defs)) stop(sprintf("Unknown accumulator '%s'", acc_id))
  acc_def <- ctx$trial_accs[[acc_id]] %||% acc_defs[[acc_id]]
  if (is.null(acc_def)) stop(sprintf("Accumulator '%s' is not available in this trial context", acc_id))
  if (is.null(ctx$acc_times[[acc_id]])) {
    ctx$acc_times[[acc_id]] <- .sample_accumulator(acc_def, ctx)
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

.resolve_param_row_acc_id <- function(row, acc_ids) {
  pick_id <- function(val) {
    if (is.null(val) || length(val) == 0L) return(NULL)
    v <- val[[1]]
    if (is.null(v) || (length(v) == 1L && is.na(v))) return(NULL)
    v
  }
  raw <- NULL
  if ("accumulator_id" %in% names(row)) {
    raw <- pick_id(row[["accumulator_id"]])
    if (!is.null(raw)) {
      id <- as.character(raw)
      if (nzchar(id)) return(id)
    }
  }
  if ("accumulator" %in% names(row)) {
    raw <- pick_id(row[["accumulator"]])
  } else if ("accumulator_index" %in% names(row)) {
    raw <- pick_id(row[["accumulator_index"]])
  } else if ("acc_idx" %in% names(row)) {
    raw <- pick_id(row[["acc_idx"]])
  }
  if (!is.null(raw)) {
    if (is.numeric(raw)) {
      idx <- as.integer(raw)
      if (!is.na(idx) && idx >= 1L && idx <= length(acc_ids)) {
        return(acc_ids[[idx]])
      }
    } else {
      id <- as.character(raw)
      if (nzchar(id)) return(id)
    }
  }
  stop("Parameter rows must include a valid 'accumulator_id' or index column")
}

.param_override_columns <- function(params_rows) {
  base_cols <- c(
    "trial", "component", "accumulator", "accumulator_id",
    "accumulator_index", "acc_idx", "type", "role",
    "outcome", "rt", "params", "condition",
    "component_weight", "params_hash"
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

.generator_component_hash <- function(rows_df) {
  if (is.null(rows_df) || nrow(rows_df) == 0L) return("__base__")
  if ("params_hash" %in% names(rows_df)) {
    vals <- rows_df$params_hash
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) > 0L) return(vals[[1]])
  }
  tmp <- rows_df
  if (!is.null(tmp)) tmp$params_hash <- NULL
  val <- .param_rows_hash_value(tmp)
  if (is.null(val) || !nzchar(val)) "__base__" else val
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
        prob <- if (val < 0 || val > 1) plogis(val) else as.numeric(val)
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
    if (has_shared && identical(tmp, base_def$q %||% 0)) {
      q_override <- NULL
    } else {
      q_override <- tmp
    }
  }
  # q values are supplied on the logit scale; convert to probability when used
  q_prob_override <- if (!is.null(q_override)) plogis(q_override) else NULL
  base_shared_prob <- if (!is.null(base_def$shared_trigger_q)) plogis(base_def$shared_trigger_q) else NULL
  base_q_prob <- if (!is.null(base_def$q)) plogis(base_def$q) else NULL
  gate_prob <- q_prob_override %||% base_shared_prob %||% base_q_prob %||% 0
  q_val <- if (has_shared) 0 else (q_prob_override %||% base_q_prob %||% 0)
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
                            trial_shared_triggers = NULL) {
  comp_info <- prep[["components"]]
  component_attr <- comp_info$attrs[[component]]
  if (is.null(component_attr)) component_attr <- list()
  ctx <- list(
    model = prep,
    component = component,
    acc_times = new.env(parent = emptyenv()),
    pool_cache = new.env(parent = emptyenv()),
    event_cache = new.env(parent = emptyenv()),
    trial_accs = trial_accs %||% list(),
    trial_shared_triggers = trial_shared_triggers %||% list(),
    shared_trigger_state = new.env(parent = emptyenv())
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

#' Simulate trials from parameter rows
#'
#' @param structure Model structure
#' @param params_df Parameter data frame
#' @param component_weights Optional component weight overrides
#' @param seed Optional RNG seed
#' @param keep_detail Whether to retain per-trial detail
#' @param ... Unused; for S3 compatibility
#' @return Data frame of simulated outcomes/RTs (with optional detail)
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' structure <- finalize_model(spec)
#' params <- c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0)
#' df <- build_params_df(spec, params, n_trials = 3)
#' simulate(structure, df, seed = 123)
#' @export
simulate <- function(structure,
                     params_df,
                     component_weights = NULL,
                     seed = NULL,
                     keep_detail = FALSE,
                     ...) {
  UseMethod("simulate")
}

#' @rdname simulate
#' @export
simulate.model_structure <- function(structure,
                                     params_df,
                                     component_weights = NULL,
                                     seed = NULL,
                                     keep_detail = FALSE,
                                     ...) {
  if (is.null(params_df) || nrow(params_df) == 0L) {
    stop("Parameter data frame must contain at least one row")
  }
  structure <- .as_model_structure(structure)
  prep <- structure$prep
  comp_table <- structure$components
  comp_ids <- comp_table$component_id
  params_df <- .params_hash_attach(params_df)
  if (!"trial" %in% names(params_df)) {
    params_df$trial <- 1L
  }
  params_df$trial <- params_df$trial
  if ("component" %in% names(params_df)) {
    params_df$component <- as.character(params_df$component)
  }
  if (!is.null(component_weights)) {
    if (!all(c("trial", "component", "weight") %in% names(component_weights))) {
      stop("component_weights must have columns 'trial', 'component', and 'weight'")
    }
    component_weights$component <- as.character(component_weights$component)
  }
  if (!is.null(seed)) set.seed(seed)
  trial_ids <- unique(params_df$trial)
  trial_ids <- trial_ids[order(trial_ids)]
  outcomes <- character(length(trial_ids))
  rts <- numeric(length(trial_ids))
  comp_record <- rep(NA_character_, length(trial_ids))
  details <- if (keep_detail) vector("list", length(trial_ids)) else NULL

  param_state_cache <- new.env(parent = emptyenv())
  get_component_state <- function(component_label, rows_df) {
    comp_key <- component_label %||% "__default__"
    hash_val <- .generator_component_hash(rows_df)
    cache_key <- paste(comp_key, hash_val, sep = "::")
    entry <- param_state_cache[[cache_key]]
    if (is.null(entry)) {
      entry <- .generator_param_state_from_rows(prep, rows_df)
      param_state_cache[[cache_key]] <- entry
    }
    entry
  }

  for (i in seq_along(trial_ids)) {
    tid <- trial_ids[[i]]
    trial_rows <- params_df[params_df$trial == tid, , drop = FALSE]
    available_components <- comp_ids
    if ("component" %in% names(trial_rows)) {
      trial_comp <- unique(trial_rows$component)
      trial_comp <- trial_comp[!is.na(trial_comp)]
      if (length(trial_comp) > 0L) {
        available_components <- intersect(comp_ids, trial_comp)
      }
    }
    if (length(available_components) == 0L) available_components <- comp_ids
    weights <- comp_table$weight[match(available_components, comp_table$component_id)]
    if (!is.null(component_weights)) {
      w_rows <- component_weights[
        component_weights$trial == tid &
          component_weights$component %in% available_components,
        ,
        drop = FALSE
      ]
      if (nrow(w_rows) > 0L) {
        weights <- w_rows$weight
      }
    }
    if (length(weights) != length(available_components) || all(is.na(weights))) {
      weights <- rep(1 / length(available_components), length(available_components))
    } else {
      if (any(!is.finite(weights) | weights < 0)) {
        stop("Component weights must be non-negative finite numbers")
      }
      total_w <- sum(weights)
      if (!is.finite(total_w) || total_w <= 0) {
        weights <- rep(1 / length(available_components), length(available_components))
      } else {
        weights <- weights / total_w
      }
    }
    chosen_component <- if (length(available_components) == 1L) {
      available_components[[1]]
    } else {
      sample(available_components, size = 1L, prob = weights)
    }
    comp_record[[i]] <- chosen_component
    component_rows <- .generator_component_rows(trial_rows, chosen_component)
    param_state <- get_component_state(chosen_component, component_rows)
    result <- .simulate_trial(
      prep,
      chosen_component,
      keep_detail = keep_detail,
      trial_accs = param_state$accs,
      trial_shared_triggers = param_state$shared
    )
    outcomes[[i]] <- result$outcome
    rts[[i]] <- result$rt
    if (keep_detail) details[[i]] <- result$detail
  }

  out_df <- data.frame(
    trial = trial_ids,
    outcome = outcomes,
    rt = rts,
    stringsAsFactors = FALSE
  )
  if (length(comp_ids) > 1L || "component" %in% names(params_df)) {
    out_df$component <- comp_record
  }
  if (keep_detail) attr(out_df, "details") <- details
  out_df
}

#' @rdname simulate
#' @export
simulate.default <- function(structure, ...) {
  stop("simulate() expects a model_structure", call. = FALSE)
}
