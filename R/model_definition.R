# ------------------------------------------------------------------------------
# Observation helpers
# ------------------------------------------------------------------------------

.validate_n_outcomes <- function(n_outcomes) {
  if (length(n_outcomes) != 1L || is.null(n_outcomes) || is.na(n_outcomes)) {
    stop("n_outcomes must be a single non-missing value")
  }
  if (!is.numeric(n_outcomes) || !is.finite(n_outcomes)) {
    stop("n_outcomes must be a finite numeric value")
  }
  if (!isTRUE(all.equal(n_outcomes, round(n_outcomes)))) {
    stop("n_outcomes must be an integer value")
  }
  out <- as.integer(round(n_outcomes))
  if (out < 1L) {
    stop("n_outcomes must be >= 1")
  }
  out
}

.normalize_onset_spec <- function(onset, acc_ids, pool_ids, outcome_labels) {
  if (is.numeric(onset) && length(onset) == 1L && is.finite(onset)) {
    return(list(kind = "absolute", value = as.numeric(onset)))
  }
  if (!is.list(onset) || !identical(onset$kind, "after")) {
    stop("onset must be a finite number or after(source, lag)")
  }
  source <- onset$source
  lag <- onset$lag
  if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
    stop("after() requires one non-empty source id")
  }
  if (!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag < 0) {
    stop("after() lag must be one finite non-negative number")
  }
  source_kind <- if (source %in% acc_ids) {
    "accumulator"
  } else if (source %in% pool_ids) {
    "pool"
  } else if (source %in% outcome_labels) {
    stop("Onset sources must be accumulator or pool ids, not outcome labels")
  } else {
    stop("Onset source '", source, "' is not a declared accumulator or pool")
  }
  list(kind = "after", source = source, source_kind = source_kind, lag = lag)
}

.expand_pool_accumulator_dependencies <- function(pool_id, pool_defs, acc_ids, stack = character(0)) {
  if (pool_id %in% stack) {
    cycle <- c(stack, pool_id)
    stop("Pool dependency cycle detected while resolving onsets: ", paste(cycle, collapse = " -> "))
  }
  members <- pool_defs[[pool_id]]$members
  deps <- character(0)
  for (member in members) {
    if (member %in% acc_ids) {
      deps <- c(deps, member)
    } else {
      deps <- c(
        deps,
        .expand_pool_accumulator_dependencies(member, pool_defs, acc_ids, c(stack, pool_id))
      )
    }
  }
  unique(deps)
}

.find_onset_cycle_path <- function(nodes, adjacency) {
  state <- setNames(integer(length(nodes)), nodes)
  stack <- character(0)
  found <- character(0)

  visit <- function(node) {
    s <- state[[node]]
    if (s == 1L) {
      idx <- match(node, stack)
      if (is.na(idx)) {
        found <<- c(node, node)
      } else {
        found <<- c(stack[idx:length(stack)], node)
      }
      return(TRUE)
    }
    if (s == 2L) {
      return(FALSE)
    }
    state[[node]] <<- 1L
    stack <<- c(stack, node)
    children <- adjacency[[node]]
    for (child in children) {
      if (visit(child)) {
        return(TRUE)
      }
    }
    stack <<- stack[-length(stack)]
    state[[node]] <<- 2L
    FALSE
  }

  for (node in nodes) {
    if (state[[node]] == 0L && visit(node)) {
      break
    }
  }
  found
}

.validate_onset_dependencies <- function(acc_defs, pool_defs) {
  acc_ids <- names(acc_defs)
  dependencies <- setNames(vector("list", length(acc_ids)), acc_ids)
  for (acc_id in acc_ids) {
    spec <- acc_defs[[acc_id]]$onset_spec
    dependencies[[acc_id]] <- if (!identical(spec$kind, "after")) {
      character(0)
    } else if (identical(spec$source_kind, "accumulator")) {
      spec$source
    } else {
      .expand_pool_accumulator_dependencies(spec$source, pool_defs, acc_ids)
    }
  }
  cycle <- .find_onset_cycle_path(acc_ids, dependencies)
  if (length(cycle)) {
    stop("Chained onset dependency cycle detected: ", paste(cycle, collapse = " -> "))
  }
  invisible(NULL)
}

.deterministic_source_signature <- function(source, acc_ids, pool_defs) {
  if (source %in% acc_ids) {
    return(paste0("acc:", source))
  }
  pool <- pool_defs[[source]]
  members <- pool$members
  k <- pool$k
  if (k == 1L && length(members) == 1L) {
    return(.deterministic_source_signature(members[[1L]], acc_ids, pool_defs))
  }
  paste0("pool:", source)
}

.outcome_is_direct_event <- function(outcome_def, acc_ids, pool_ids) {
  identical(outcome_def$expr$kind, "event") &&
    outcome_def$expr$source %in% c(acc_ids, pool_ids)
}

# Validator for multi-outcome readout declarations.
.validate_multi_outcome_dsl <- function(model_or_prep) {
  obs <- model_or_prep$observation
  n_outcomes <- obs$n_outcomes
  if (n_outcomes <= 1L) {
    return(invisible(NULL))
  }
  outcomes <- model_or_prep$outcomes
  n_defined <- length(outcomes)
  if (n_outcomes > n_defined) {
    stop(sprintf(
      "n_outcomes (%d) cannot exceed number of declared outcomes (%d)",
      n_outcomes, n_defined
    ))
  }

  acc_ids <- names(model_or_prep$accumulators)
  pool_ids <- names(model_or_prep$pools)
  issues <- character(0)
  direct_sources <- character(0)
  direct_labels <- character(0)
  for (i in seq_along(outcomes)) {
    out <- outcomes[[i]]
    label <- out$label
    opts <- out$options
    is_direct <- .outcome_is_direct_event(out, acc_ids, pool_ids)
    if (!is_direct) {
      issues <- c(issues, sprintf("%s (must be a direct event outcome)", label))
    } else {
      direct_sources <- c(direct_sources, out$expr$source)
      direct_labels <- c(direct_labels, label)
    }
    if (!is.null(opts$guess)) {
      issues <- c(issues, sprintf("%s (guess option not supported)", label))
    }
    if (!is.null(opts$map_outcome_to)) {
      issues <- c(issues, sprintf("%s (map_outcome_to option not supported)", label))
    }
  }
  if (length(issues) > 0L) {
    stop(
      "n_outcomes > 1 currently supports only direct event outcomes ",
      "with no guess/map_outcome_to options. Invalid outcomes: ",
      paste(unique(issues), collapse = ", ")
    )
  }

  pool_defs <- model_or_prep$pools
  signatures <- vapply(direct_sources, function(source) {
    .deterministic_source_signature(source, acc_ids, pool_defs)
  }, character(1))
  overlap_sig <- unique(signatures[duplicated(signatures) | duplicated(signatures, fromLast = TRUE)])
  if (length(overlap_sig) > 0L) {
    groups <- vapply(overlap_sig, function(sig) {
      lbls <- direct_labels[signatures == sig]
      paste(lbls, collapse = " = ")
    }, character(1))
    stop(
      "n_outcomes > 1 rejects deterministic overlapping outcomes. Overlap groups: ",
      paste(groups, collapse = "; ")
    )
  }
  invisible(NULL)
}

.build_outcome_component_maps <- function(outcomes, component_ids) {
  allowed <- observed <- setNames(
    rep(list(character(0)), length(component_ids)),
    component_ids
  )
  labels <- names(outcomes)
  for (i in seq_along(outcomes)) {
    options <- outcomes[[i]]$options
    components <- options$component %||% component_ids
    observed_labels <- labels[[i]]
    if (!is.null(options$guess)) {
      observed_labels <- options$guess$labels
    }
    if (!is.null(options$map_outcome_to)) {
      observed_labels <- options$map_outcome_to
    }
    observed_labels <- as.character(observed_labels[!is.na(observed_labels)])
    for (component in components) {
      allowed[[component]] <- c(allowed[[component]], labels[[i]])
      observed[[component]] <- c(observed[[component]], observed_labels)
    }
  }
  list(
    allowed = lapply(allowed, unique),
    observed = lapply(observed, unique)
  )
}

# ------------------------------------------------------------------------------
# Expression parsing utilities
# ------------------------------------------------------------------------------

.is_expr_node <- function(x) is.list(x) && length(x) > 0 && !is.null(x$kind)

.expr_from_value <- function(val) {
  if (.is_expr_node(val)) {
    return(val)
  }
  if (is.character(val) && length(val) == 1L) {
    return(list(kind = "event", source = val))
  }
  if (is.symbol(val)) {
    return(list(kind = "event", source = as.character(val)))
  }
  if (is.call(val)) {
    return(.parse_expr_call(val))
  }
  stop("Cannot interpret expression component of type '", typeof(val), "'")
}

.parse_expr_call <- function(call) {
  op <- as.character(call[[1]])
  if (op == "(") {
    return(.expr_from_value(call[[2]]))
  }
  if (op %in% c("&", "&&")) {
    parts <- lapply(as.list(call)[-1], .expr_from_value)
    return(list(kind = "and", args = parts))
  }
  if (op %in% c("|", "||")) {
    parts <- lapply(as.list(call)[-1], .expr_from_value)
    return(list(kind = "or", args = parts))
  }
  if (identical(op, "!")) {
    return(list(kind = "not", arg = .expr_from_value(call[[2]])))
  }
  stop(sprintf("Unsupported token '%s' in expression", op))
}

.build_expr <- function(expr) {
  if (.is_expr_node(expr)) {
    return(expr)
  }
  .expr_from_value(expr)
}

#' Turn a response rule into an internal expression
#'
#' Use this when you want to write a response rule programmatically rather than
#' through the helper functions such as `all_of()` or `inhibit()`.
#'
#' @param expr Expression or symbol describing an event or blocking rule.
#' @return An expression object used inside model specifications.
#' @examples
#' build_outcome_expr(quote(A & !B))
#' @export
build_outcome_expr <- function(expr) {
  .build_expr(expr)
}

# ------------------------------------------------------------------------------
# Public DSL helpers
# ------------------------------------------------------------------------------

#' Define a response that is blocked by another process
#'
#' @param reference Response rule or accumulator label to be blocked.
#' @param by Blocking process or expression.
#' @return A guarded expression object.
#' @examples
#' inhibit("A", "B")
#' @export
inhibit <- function(reference, by) {
  list(
    kind = "guard",
    blocker = .build_expr(by),
    reference = .build_expr(reference)
  )
}

#' Define a response that occurs when the first listed process finishes
#'
#' @param ... Accumulator labels or expression objects to combine with OR.
#' @return An expression object.
#' @examples
#' first_of("A", "B")
#' @export
first_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("first_of() requires at least one argument")
  list(kind = "or", args = lapply(args, .expr_from_value))
}

#' Define a response that requires several processes to finish
#'
#' @param ... Accumulator labels or expression objects to combine with AND.
#' @return An expression object.
#' @examples
#' all_of("A", "B")
#' @export
all_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("all_of() requires at least one argument")
  list(kind = "and", args = lapply(args, .expr_from_value))
}

#' Define the absence of an event
#'
#' @param expr Accumulator label or expression to negate.
#' @return An expression object.
#' @export
#' @examples
#' none_of("A")
none_of <- function(expr) {
  list(kind = "not", arg = .build_expr(expr))
}

#' Start one accumulator after another process finishes
#'
#' This is useful for staged or contingent architectures, where one process can
#' only begin after an earlier accumulator or pool has finished.
#'
#' @param source Accumulator or pool label that must finish first.
#' @param lag Optional non-negative delay added after `source` finishes.
#' @return A chained-onset specification for `add_accumulator(onset = ...)`.
#' @examples
#' after("A")
#' after("pool1", lag = 0.05)
#' @export
after <- function(source, lag = 0) {
  list(kind = "after", source = source, lag = lag)
}

# ------------------------------------------------------------------------------
# Model builder
# ------------------------------------------------------------------------------

#' Start a race-model specification
#'
#' @param n_outcomes Number of ordered observed responses to retain per trial.
#'   Use `1` for standard choice/RT data, `2` when you also observe the second
#'   finishing response, and so on.
#' @return A `race_spec` object.
#' @examples
#' race_spec()
#' @export
race_spec <- function(n_outcomes = 1L) {
  structure(list(
    accumulators = list(),
    pools = list(),
    outcomes = list(),
    triggers = list(),
    parameters = .normalize_parameter_spec(),
    components = list(),
    mixture_options = list(),
    observation = list(n_outcomes = .validate_n_outcomes(n_outcomes))
  ), class = "race_spec")
}

.validate_race_spec_input <- function(spec, fn_name) {
  if (inherits(spec, "race_spec")) {
    return(spec)
  }
  stop(
    sprintf(
      "%s() expects a model specification created with race_spec().",
      fn_name
    ),
    call. = FALSE
  )
}

#' Add an accumulator to a model
#'
#' @param spec A `race_spec` object.
#' @param id Label for the accumulator.
#' @param dist Distribution family used for that accumulator.
#' @param onset Start time for the accumulator. This can be a fixed numeric
#'   onset or a chained onset created with `after()`.
#' @return The updated `race_spec`.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' @export
add_accumulator <- function(spec, id, dist, onset = 0) {
  spec <- .validate_race_spec_input(spec, "add_accumulator")
  spec$accumulators[[length(spec$accumulators) + 1L]] <- list(
    id = id,
    dist = dist,
    onset = onset
  )
  spec
}

#' Pool several accumulators under a shared label
#'
#' Pools let you talk about several accumulators as one source when defining
#' observed responses.
#'
#' @param spec A `race_spec` object.
#' @param id Label for the pool.
#' @param members Accumulator labels included in the pool.
#' @param k Threshold for a `k`-of-`n` pool rule.
#' @return The updated `race_spec`.
#' @examples
#' spec <- race_spec()
#' spec <- add_pool(spec, "P1", members = c("A", "B"), k = 1L)
#' @export
add_pool <- function(spec, id, members, k = 1L) {
  spec <- .validate_race_spec_input(spec, "add_pool")
  if (missing(members) || length(members) == 0) {
    stop("Pool must define at least one member")
  }
  spec$pools[[length(spec$pools) + 1L]] <- list(
    id = id,
    members = as.character(members),
    rule = list(kind = "k_of_n", k = k)
  )
  spec
}

#' Define an observed response
#'
#' @param spec A `race_spec` object.
#' @param label Response label that should appear in the behavioral data.
#' @param expr Rule describing when that response is observed.
#' @param options Optional response settings.
#' @return The updated `race_spec`.
#' @examples
#' spec <- race_spec()
#' spec <- add_outcome(spec, "A_win", "A")
#' @export
add_outcome <- function(spec, label, expr, options = list()) {
  spec <- .validate_race_spec_input(spec, "add_outcome")
  options <- options %||% list()
  if (!is.list(options)) {
    stop("Outcome options must be a list", call. = FALSE)
  }
  unknown_options <- setdiff(names(options), c("component", "map_outcome_to", "guess"))
  if (length(unknown_options) > 0L) {
    stop(
      sprintf("Unknown outcome option(s): %s", paste(unknown_options, collapse = ", ")),
      call. = FALSE
    )
  }
  spec$outcomes[[length(spec$outcomes) + 1L]] <- list(
    label = label,
    expr = build_outcome_expr(expr),
    options = options
  )
  spec
}

#' Define a mixture component
#'
#' Components are useful when trials can come from qualitatively different
#' processing modes, such as fast versus slow processing. Component declarations
#' define membership only; component probabilities are configured with
#' `set_mixture()`.
#'
#' @param spec A `race_spec` object.
#' @param id Component label.
#' @param members Accumulator labels that belong to this component.
#' @param n_outcomes Optional component-specific override for the number of
#'   observed ordered responses.
#' @return The updated `race_spec`.
#' @export
add_component <- function(spec, id, members, n_outcomes = NULL) {
  spec <- .validate_race_spec_input(spec, "add_component")
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Component must specify members")
  }
  comp_attrs <- list()
  if (!is.null(n_outcomes)) {
    comp_attrs$n_outcomes <- .validate_n_outcomes(n_outcomes)
  }
  spec$components[[length(spec$components) + 1L]] <- list(
    id = id,
    members = as.character(members),
    attrs = comp_attrs
  )
  spec
}

#' Add a shared absence trigger
#'
#' A trigger is a named absence-probability parameter. All members in one
#' trigger call share the same absence draw. Use separate trigger calls for
#' independent absence draws.
#'
#' @param spec A `race_spec` object.
#' @param name Trigger parameter name.
#' @param members Accumulator labels controlled by the shared absence draw.
#' @return The updated `race_spec`.
#' @export
add_trigger <- function(spec, name, members) {
  spec <- .validate_race_spec_input(spec, "add_trigger")
  if (missing(name) || is.null(name) || length(name) != 1L || !nzchar(name)) {
    stop("Trigger must have a non-empty name", call. = FALSE)
  }
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Trigger must specify members")
  }
  spec$triggers[[length(spec$triggers) + 1L]] <- list(
    id = as.character(name),
    members = as.character(members)
  )
  spec
}

#' Control parameter grouping and names
#'
#' Parameters are grouped by compatible type by default. For example, two
#' lognormal accumulators expose one `m`, one `s`, and one `t0` parameter unless
#' you ask for specific parameters to be separate. Triggers expose their trigger name
#' directly, and sampled mixtures expose automatic `p.<component>` parameters.
#'
#' @param spec A `race_spec` object.
#' @param separate Named list. Each name is a grouped public parameter, and each
#'   value is one or more accumulator ids to split from that group. Use `TRUE`
#'   to split every member of a group.
#' @param share Named list mapping a new public name to default public names or
#'   internal parameter names that should share one value.
#' @param rename Named character vector mapping current public names to new
#'   public names.
#' @return The updated `race_spec`.
#' @examples
#' spec <- race_spec() |>
#'   add_accumulator("go", "lognormal") |>
#'   add_accumulator("stop", "lognormal") |>
#'   add_outcome("go", "go") |>
#'   add_outcome("stop", "stop") |>
#'   set_parameters(
#'     separate = list(m = c("go", "stop")),
#'     rename = c(s = "spread", t0 = "onset")
#'   )
#'
#' par_names(spec)
#' @export
set_parameters <- function(spec, separate = NULL, share = NULL, rename = NULL) {
  spec <- .validate_race_spec_input(spec, "set_parameters")
  spec$parameters <- .normalize_parameter_spec(separate = separate, share = share, rename = rename)
  spec
}

#' Control how mixture components are combined
#'
#' Fixed mixtures use known component probabilities. Sampled mixtures expose
#' automatic `p.<component>` parameters for every non-reference component; the
#' reference component receives the residual probability.
#'
#' @param spec A `race_spec` object.
#' @param mode Mixture mode. Fixed mixtures use known component probabilities;
#'   sampled mixtures estimate probabilities for all non-reference components.
#' @param weights Named numeric component probabilities for fixed mixtures. If
#'   `NULL`, fixed mixtures use uniform component probabilities.
#' @param reference Reference component for sampled mixtures. Its probability is
#'   the residual probability after non-reference component probabilities.
#' @return The updated `race_spec`.
#' @export
set_mixture <- function(spec, mode = c("fixed", "sample"), weights = NULL, reference = NULL) {
  spec <- .validate_race_spec_input(spec, "set_mixture")
  mode <- match.arg(mode)
  if (identical(mode, "sample") && !is.null(weights)) {
    stop("Sampled mixtures use automatic p.<component> parameters, not fixed weights", call. = FALSE)
  }
  if (!is.null(weights)) {
    if (!is.numeric(weights) || is.null(names(weights)) || any(!nzchar(names(weights)))) {
      stop("Fixed mixture weights must be a named numeric vector", call. = FALSE)
    }
    if (any(!is.finite(weights) | weights < 0)) {
      stop("Fixed mixture weights must be non-negative finite probabilities", call. = FALSE)
    }
    if (!isTRUE(all.equal(sum(weights), 1, tolerance = 1e-8))) {
      stop("Fixed mixture weights must sum to 1", call. = FALSE)
    }
  }
  spec$mixture_options <- list(
    mode = mode,
    weights = weights,
    reference = reference
  )
  spec
}

# ----------------------------------------------------------------------
# Model normalization/finalization (shared by simulation and likelihood)
# ----------------------------------------------------------------------

.model_ids <- function(items, field, type) {
  ids <- vapply(items, function(item) item[[field]] %||% "", character(1))
  if (any(!nzchar(ids))) {
    stop(type, " ids must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop(type, " ids must be unique", call. = FALSE)
  }
  ids
}

.expression_sources <- function(expr) {
  switch(
    expr$kind,
    event = expr$source,
    and = ,
    or = unlist(lapply(expr$args, .expression_sources), use.names = FALSE),
    not = .expression_sources(expr$arg),
    guard = c(
      .expression_sources(expr$reference),
      .expression_sources(expr$blocker)
    ),
    stop("Unsupported expression kind '", expr$kind, "'", call. = FALSE)
  )
}

.validate_model_spec <- function(model) {
  if (!length(model$accumulators) || !length(model$outcomes)) {
    stop("Model must define accumulators and outcomes", call. = FALSE)
  }
  acc_ids <- .model_ids(model$accumulators, "id", "Accumulator")
  pool_ids <- .model_ids(model$pools, "id", "Pool")
  component_ids <- .model_ids(model$components, "id", "Component")
  .model_ids(model$triggers, "id", "Trigger")
  .validate_n_outcomes(model$observation$n_outcomes)
  if (length(intersect(acc_ids, pool_ids))) {
    stop("Accumulator and pool ids must be distinct", call. = FALSE)
  }
  invisible(lapply(model$accumulators, function(acc) dist_registry(acc$dist)))

  pool_defs <- setNames(model$pools, pool_ids)
  for (pool in model$pools) {
    members <- as.character(pool$members)
    if (any(!members %in% c(acc_ids, pool_ids))) {
      stop("Pool '", pool$id, "' references an unknown member", call. = FALSE)
    }
    k <- pool$rule$k
    if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
        k != as.integer(k) || k < 1L || k > length(members)) {
      stop("Pool '", pool$id, "' requires an integer k between 1 and its member count", call. = FALSE)
    }
    .expand_pool_accumulator_dependencies(pool$id, pool_defs, acc_ids)
  }

  for (component in model$components) {
    if (any(!component$members %in% acc_ids)) {
      stop("Component '", component$id, "' references an unknown accumulator", call. = FALSE)
    }
  }
  trigger_members <- unlist(lapply(model$triggers, `[[`, "members"), use.names = FALSE)
  if (any(!trigger_members %in% acc_ids)) {
    stop("Triggers may reference only declared accumulators", call. = FALSE)
  }
  if (anyDuplicated(trigger_members)) {
    stop("An accumulator may belong to only one trigger", call. = FALSE)
  }

  known_components <- if (length(component_ids)) component_ids else "__default__"
  outcome_labels <- vapply(
    model$outcomes, function(outcome) outcome$label %||% "", character(1)
  )
  if (any(!nzchar(outcome_labels))) {
    stop("Outcome labels must be non-empty", call. = FALSE)
  }
  for (label in unique(outcome_labels[duplicated(outcome_labels)])) {
    used_components <- character(0)
    for (index in which(outcome_labels == label)) {
      components <- model$outcomes[[index]]$options$component %||% known_components
      if (length(intersect(used_components, components))) {
        stop(
          "Outcome '", label,
          "' has multiple definitions active in the same component",
          call. = FALSE
        )
      }
      used_components <- c(used_components, components)
    }
  }
  for (outcome in model$outcomes) {
    unknown_sources <- setdiff(
      unique(.expression_sources(outcome$expr)), c(acc_ids, pool_ids)
    )
    if (length(unknown_sources)) {
      stop("Outcome '", outcome$label, "' references unknown source(s): ",
           paste(unknown_sources, collapse = ", "), call. = FALSE)
    }
    options <- outcome$options
    components <- options$component
    if (!is.null(components) &&
        (!is.character(components) || !length(components) || anyNA(components) ||
         any(!nzchar(components)) || anyDuplicated(components) ||
         any(!components %in% known_components))) {
      stop(
        "Outcome '", outcome$label,
        "' must name one or more unique declared components",
        call. = FALSE
      )
    }
    if (!is.null(options$map_outcome_to) &&
        !is.na(options$map_outcome_to) &&
        !options$map_outcome_to %in% outcome_labels) {
      stop("Outcome '", outcome$label, "' maps to an unknown label", call. = FALSE)
    }
    if (!is.null(options$guess)) {
      guess <- options$guess
      valid <- is.character(guess$labels) && is.numeric(guess$weights) &&
        length(guess$labels) == length(guess$weights) &&
        length(guess$labels) > 0L &&
        all(guess$labels %in% outcome_labels) &&
        all(is.finite(guess$weights) & guess$weights >= 0) &&
        isTRUE(all.equal(sum(guess$weights), 1, tolerance = 1e-8)) &&
        (guess$rt_policy %||% "keep") %in% c("keep", "na")
      if (!valid) {
        stop("Outcome '", outcome$label, "' has an invalid guess policy", call. = FALSE)
      }
    }
  }
  invisible(model)
}

.prepare_acc_defs <- function(model) {
  accs <- model$accumulators
  acc_ids <- vapply(accs, `[[`, character(1), "id")
  pool_ids <- vapply(model$pools, `[[`, character(1), "id")
  outcome_labels <- vapply(model$outcomes, `[[`, character(1), "label")

  defs <- setNames(vector("list", length(accs)), acc_ids)
  for (acc in accs) {
    onset_spec <- .normalize_onset_spec(
      onset = acc$onset,
      acc_ids = acc_ids,
      pool_ids = pool_ids,
      outcome_labels = outcome_labels
    )
    onset_value <- if (identical(onset_spec$kind, "absolute")) onset_spec$value else 0
    defs[[acc$id]] <- list(
      id = acc$id,
      dist = acc$dist,
      onset = onset_value,
      onset_spec = onset_spec,
      components = character(0),
      shared_trigger_id = NULL
    )
  }

  for (component in model$components) {
    for (member in component$members) {
      defs[[member]]$components <- c(defs[[member]]$components, component$id)
    }
  }

  shared_triggers <- setNames(
    model$triggers,
    vapply(model$triggers, `[[`, character(1), "id")
  )
  for (trigger in model$triggers) {
    for (member in trigger$members) {
      defs[[member]]$shared_trigger_id <- trigger$id
    }
  }
  list(acc = defs, shared_triggers = shared_triggers)
}

.prepare_components <- function(model) {
  comps <- model$components
  mix_opts <- model$mixture_options
  mode <- mix_opts$mode %||% "fixed"
  reference <- mix_opts$reference %||% NA_character_
  if (length(comps) == 0) {
    return(list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list()),
      mode = "fixed",
      reference = "__default__"
    ))
  }

  ids <- vapply(comps, `[[`, character(1), "id")
  attrs <- setNames(vector("list", length(ids)), ids)

  if (identical(mode, "sample")) {
    if (is.na(reference) || !nzchar(reference)) {
      reference <- ids[[length(ids)]]
    } else if (!reference %in% ids) {
      stop("mixture reference '", reference, "' must match a component id")
    }
    weights <- stats::setNames(rep(1 / length(ids), length(ids)), ids)
  } else {
    weights <- mix_opts$weights %||% NULL
    if (is.null(weights)) {
      weights <- stats::setNames(rep(1 / length(ids), length(ids)), ids)
    } else {
      extra <- setdiff(names(weights), ids)
      missing <- setdiff(ids, names(weights))
      if (length(extra) > 0L || length(missing) > 0L) {
        stop("Fixed mixture weights must be named for every component and no unknown components", call. = FALSE)
      }
      weights <- weights[ids]
    }
    if (is.na(reference) || !nzchar(reference)) {
      reference <- ids[[1]]
    }
  }

  for (i in seq_along(ids)) {
    cmp_id <- ids[[i]]
    attrs_cmp <- comps[[i]]$attrs
    if (identical(mode, "sample") && !identical(cmp_id, reference)) {
      attrs_cmp$weight_param <- paste0("p.", cmp_id)
    }
    attrs[[cmp_id]] <- attrs_cmp
  }

  list(
    ids = ids,
    weights = as.numeric(weights),
    attrs = attrs,
    mode = mode,
    reference = reference
  )
}

.prepare_pool_defs <- function(model) {
  pools <- model$pools
  setNames(lapply(pools, function(pool) {
    list(id = pool$id, members = pool$members, k = pool$rule$k)
  }), vapply(pools, `[[`, character(1), "id"))
}

.prepare_model <- function(model) {
  .validate_model_spec(model)
  acc_prep <- .prepare_acc_defs(model)
  acc_defs <- acc_prep$acc
  pool_defs <- .prepare_pool_defs(model)
  outcome_defs <- setNames(
    model$outcomes,
    vapply(model$outcomes, `[[`, character(1), "label")
  )
  component_defs <- .prepare_components(model)
  observation <- model$observation
  component_overrides <- vapply(
    component_defs$attrs,
    function(attributes) attributes$n_outcomes %||% NA_integer_,
    integer(1)
  )
  max_component_outcomes <- max(c(1L, component_overrides), na.rm = TRUE)
  .validate_onset_dependencies(acc_defs, pool_defs)
  observation$global_n_outcomes <- observation$n_outcomes
  observation$component_n_outcomes <- as.list(component_overrides[!is.na(component_overrides)])
  observation$n_outcomes <- max(observation$n_outcomes, max_component_outcomes)
  prep <- list(
    accumulators = acc_defs,
    pools = pool_defs,
    outcomes = outcome_defs,
    components = component_defs,
    observation = observation,
    shared_triggers = acc_prep$shared_triggers,
    parameter_lookup = .parameter_name_lookup(model)
  )
  outcome_maps <- .build_outcome_component_maps(
    outcome_defs,
    component_defs$ids
  )
  prep$outcomes_by_component <- outcome_maps$allowed
  prep$observed_outcomes_by_component <- outcome_maps$observed
  .validate_multi_outcome_dsl(prep)
  prep
}

#' Compile a model for simulation and fitting
#'
#' This converts a human-readable model specification into the finalized object
#' used by `simulate()`, `prepare_data()`, `make_context()`, and related functions.
#'
#' @param model Model specification.
#' @return A `model_structure` object.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' finalize_model(spec)
#' @export
finalize_model <- function(model) {
  model <- .validate_race_spec_input(model, "finalize_model")
  structure <- list(
    model_spec = model,
    prep = .prepare_model(model)
  )
  class(structure) <- c("model_structure", class(structure))
  structure
}

# ------------------------------------------------------------------------------
# Parameter utilities
# ------------------------------------------------------------------------------

dist_param_names <- function(dist) {
  dist_registry(dist)$params
}

.parameter_character <- function(value, what) {
  if (is.null(value)) {
    return(character(0))
  }
  if (is.character(value)) {
    return(as.character(value))
  }
  stop(sprintf("%s must be a character vector", what), call. = FALSE)
}

.normalize_parameter_spec <- function(separate = NULL, share = NULL, rename = NULL) {
  if (is.null(separate)) {
    separate <- list()
  }
  if (!is.list(separate)) {
    stop("separate must be a named list", call. = FALSE)
  }
  separate_names <- names(separate)
  if (length(separate) > 0L && (is.null(separate_names) || any(!nzchar(separate_names)))) {
    stop("separate must be a named list", call. = FALSE)
  }
  if (anyDuplicated(separate_names)) {
    stop("separate names must be unique", call. = FALSE)
  }
  separate <- lapply(names(separate), function(nm) {
    value <- separate[[nm]]
    if (isTRUE(value)) {
      return(TRUE)
    }
    members <- .parameter_character(value, sprintf("separate[['%s']]", nm))
    members <- unique(members[nzchar(members)])
    if (length(members) == 0L) {
      stop(sprintf("separate parameter '%s' must name at least one accumulator or use TRUE", nm), call. = FALSE)
    }
    members
  }) |>
    stats::setNames(separate_names)

  if (is.null(share)) {
    share <- list()
  }
  if (!is.list(share)) {
    stop("share must be a named list", call. = FALSE)
  }
  share_names <- names(share)
  if (length(share) > 0L && (is.null(share_names) || any(!nzchar(share_names)))) {
    stop("share must be a named list", call. = FALSE)
  }
  if (anyDuplicated(share_names)) {
    stop("share names must be unique", call. = FALSE)
  }
  share <- lapply(names(share), function(nm) {
    targets <- .parameter_character(share[[nm]], sprintf("share[['%s']]", nm))
    targets <- targets[nzchar(targets)]
    if (length(targets) == 0L) {
      stop(sprintf("share parameter '%s' must reference at least one parameter", nm), call. = FALSE)
    }
    unname(targets)
  }) |>
    stats::setNames(share_names)

  if (is.null(rename)) {
    rename <- character(0)
  }
  if (!is.character(rename)) {
    stop("rename must be a named character vector like c(old = 'new')", call. = FALSE)
  }
  if (length(rename) > 0L && (is.null(names(rename)) || any(!nzchar(names(rename))))) {
    stop("rename must be a named character vector like c(old = 'new')", call. = FALSE)
  }
  rename_names <- names(rename)
  rename <- as.character(rename)
  names(rename) <- rename_names
  rename <- rename[nzchar(rename)]
  if (anyDuplicated(names(rename))) {
    stop("rename source names must be unique", call. = FALSE)
  }
  if (anyDuplicated(unname(rename))) {
    stop("rename target names must be unique", call. = FALSE)
  }

  list(
    separate = separate,
    share = share,
    rename = rename
  )
}

.mixture_weight_parameter_names <- function(spec) {
  comps <- spec$components
  if (length(comps) == 0L) {
    return(character(0))
  }
  mix <- spec$mixture_options
  if (!identical(mix$mode %||% "fixed", "sample")) {
    return(character(0))
  }
  ids <- vapply(comps, `[[`, character(1), "id")
  reference <- mix$reference %||% ids[[length(ids)]]
  setdiff(paste0("p.", ids), paste0("p.", reference))
}

.parameter_base_table <- function(spec) {
  rows <- list()
  append_row <- function(internal, public, kind, acc = NA_character_, dist = NA_character_, param = NA_character_, slot = NA_integer_) {
    rows[[length(rows) + 1L]] <<- list(
      internal = internal,
      public = public,
      kind = kind,
      racer = acc,
      dist = dist,
      param = param,
      slot = slot
    )
  }

  accs <- spec$accumulators %||% list()
  for (acc in accs) {
    acc_id <- acc$id
    dist <- tolower(acc$dist)
    params <- dist_param_names(dist)
    for (slot in seq_along(params)) {
      param <- params[[slot]]
      append_row(
        internal = paste(acc_id, param, sep = "."),
        public = param,
        kind = "distribution",
        acc = acc_id,
        dist = dist,
        param = param,
        slot = slot
      )
    }
    append_row(
      internal = paste(acc_id, "t0", sep = "."),
      public = "t0",
      kind = "t0",
      acc = acc_id,
      dist = dist,
      param = "t0",
      slot = NA_integer_
    )
  }

  trigger_ids <- vapply(spec$triggers %||% list(), function(trig) {
    trig$id %||% NA_character_
  }, character(1))
  trigger_ids <- trigger_ids[!is.na(trigger_ids) & nzchar(trigger_ids)]
  for (trigger_id in trigger_ids) {
    append_row(
      internal = trigger_id,
      public = trigger_id,
      kind = "trigger",
      param = trigger_id
    )
  }

  weight_params <- .mixture_weight_parameter_names(spec)
  for (weight_param in weight_params) {
    append_row(
      internal = weight_param,
      public = weight_param,
      kind = "mixture",
      param = weight_param
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(
      internal = character(0),
      public = character(0),
      kind = character(0),
      racer = character(0),
      dist = character(0),
      param = character(0),
      slot = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  out$slot <- as.integer(out$slot)
  collision <- intersect(
    unique(out$public[out$kind %in% c("distribution", "t0")]),
    c(trigger_ids, weight_params)
  )
  if (length(collision) > 0L) {
    stop(
      sprintf(
        "Parameter names conflict with trigger or mixture names: %s",
        paste(collision, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out
}

.default_parameter_lookup <- function(spec) {
  params <- .parameter_base_table(spec)
  if (nrow(params) == 0L) {
    return(setNames(character(0), character(0)))
  }

  distribution <- params$kind == "distribution"
  if (any(distribution)) {
    dist_rows <- params[distribution, , drop = FALSE]
    split_names <- names(which(vapply(split(dist_rows, dist_rows$param), function(group) {
      length(unique(group$slot)) > 1L
    }, logical(1))))
    for (param in split_names) {
      idx <- distribution & params$param == param
      params$public[idx] <- paste(params$dist[idx], params$param[idx], sep = ".")
    }
  }

  lookup <- setNames(params$public, params$internal)
  duplicate_public <- unique(unname(lookup)[duplicated(unname(lookup))])
  for (public in duplicate_public) {
    internals <- names(lookup)[lookup == public]
    kinds <- unique(params$kind[match(internals, params$internal)])
    if (length(kinds) > 1L || "trigger" %in% kinds || "mixture" %in% kinds) {
      lookup[internals] <- internals
    }
  }
  lookup
}

.resolve_parameter_targets <- function(targets, lookup, what) {
  out <- character(0)
  for (target in targets) {
    if (target %in% names(lookup)) {
      out <- c(out, target)
      next
    }
    if (target %in% unname(lookup)) {
      out <- c(out, names(lookup)[lookup == target])
      next
    }
    stop(sprintf("%s references unknown parameter '%s'", what, target), call. = FALSE)
  }
  unique(out)
}

.resolve_separate_targets <- function(group_name, members, lookup) {
  group_targets <- names(lookup)[lookup == group_name]
  if (length(group_targets) == 0L) {
    stop(sprintf("separate references unknown grouped parameter '%s'", group_name), call. = FALSE)
  }
  if (isTRUE(members)) {
    return(group_targets)
  }

  out <- character(0)
  for (member in members) {
    internal_name <- paste(member, group_name, sep = ".")
    if (internal_name %in% group_targets) {
      out <- c(out, internal_name)
      next
    }
    if (member %in% group_targets) {
      out <- c(out, member)
      next
    }
    stop(
      sprintf(
        "separate[['%s']] references '%s', but '%s' is not in that parameter group",
        group_name,
        member,
        internal_name
      ),
      call. = FALSE
    )
  }
  unique(out)
}

.parameter_name_lookup <- function(spec) {
  lookup <- .default_parameter_lookup(spec)
  mappings <- spec$parameters

  separate <- mappings$separate
  for (group_name in names(separate)) {
    targets <- .resolve_separate_targets(group_name, separate[[group_name]], lookup)
    lookup[targets] <- targets
  }

  seen_targets <- character(0)
  for (external_name in names(mappings$share)) {
    targets <- .resolve_parameter_targets(mappings$share[[external_name]], lookup, sprintf("share[['%s']]", external_name))
    duplicated_targets <- intersect(seen_targets, targets)
    if (length(duplicated_targets) > 0L) {
      stop(
        sprintf(
          "share targets cannot be assigned more than once: %s",
          paste(duplicated_targets, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    lookup[targets] <- external_name
    seen_targets <- c(seen_targets, targets)
  }

  rename <- mappings$rename
  pre_rename_lookup <- lookup
  for (old_name in names(rename)) {
    if (!old_name %in% unname(lookup)) {
      stop(sprintf("rename references unknown public parameter '%s'", old_name), call. = FALSE)
    }
    lookup[lookup == old_name] <- unname(rename[[old_name]])
  }
  final_groups <- split(names(lookup), unname(lookup))
  for (public_name in names(final_groups)) {
    previous_names <- unique(unname(pre_rename_lookup[final_groups[[public_name]]]))
    if (length(previous_names) > 1L) {
      stop(
        sprintf(
          "rename target '%s' would merge existing parameters (%s); use share for intentional sharing",
          public_name,
          paste(previous_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  lookup
}

.expand_parameter_values <- function(lookup, param_values) {
  required <- unique(unname(lookup))
  expanded <- setNames(numeric(length(lookup)), names(lookup))
  missing <- character(0)

  for (external_name in required) {
    targets <- names(lookup)[lookup == external_name]
    if (external_name %in% names(param_values)) {
      expanded[targets] <- as.numeric(param_values[[external_name]])[1]
      next
    }
    suffixes <- sub("^.*\\.", "", targets)
    defaultable <- length(suffixes) > 0L && all(suffixes %in% "t0")
    if (defaultable) {
      expanded[targets] <- 0
      next
    }
    missing <- c(missing, external_name)
  }

  if (length(missing) > 0L) {
    stop("Missing parameter values for: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(names(param_values), required)
  if (length(unknown) > 0L) {
    stop("Unknown parameter values: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  expanded
}

#' List the free parameters implied by a model
#'
#' @param model A `race_spec` or finalized `model_structure` object.
#' @return A character vector of parameter names.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' par_names(spec)
#' @export
par_names <- function(model) {
  lookup <- if (inherits(model, "model_structure")) {
    model$prep$parameter_lookup
  } else {
    .parameter_name_lookup(.validate_race_spec_input(model, "par_names"))
  }
  unique(unname(lookup))
}

#' Create trial-level parameter values
#'
#' This expands a named parameter vector into the trial-by-trial format expected
#' by `simulate()` and `log_likelihood()`.
#'
#' @param model Finalized model structure.
#' @param param_values Named numeric vector of parameter values.
#' @param n_trials Number of trials to generate.
#' @return A numeric parameter matrix with one row per trial/accumulator pair.
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' vals <- c(m = 0, s = 0.1)
#' build_param_matrix(finalize_model(spec), vals, n_trials = 2)
#' @export
build_param_matrix <- function(model,
                               param_values,
                               n_trials = 1L) {
  if (!inherits(model, "model_structure")) {
    stop("build_param_matrix() expects a finalized model.", call. = FALSE)
  }
  prep <- model$prep
  accs <- prep$accumulators
  if (is.null(param_values) || length(param_values) == 0L) {
    stop("Named parameter values are required")
  }
  if (is.null(names(param_values)) || any(!nzchar(names(param_values)))) {
    stop("Parameter values must be a named vector")
  }
  raw_names <- names(param_values)
  param_values <- as.numeric(param_values)
  names(param_values) <- raw_names
  param_values <- .expand_parameter_values(prep$parameter_lookup, param_values)

  acc_ids <- names(accs)
  dist_param_list <- lapply(accs, function(acc) dist_param_names(acc$dist))
  max_p <- max(vapply(dist_param_list, length, integer(1)))
  p_col_names <- paste0("p", seq_len(max_p))
  weight_params <- vapply(
    prep$components$attrs,
    function(attributes) attributes$weight_param %||% NA_character_,
    character(1)
  )
  weight_params <- unname(weight_params[!is.na(weight_params)])
  col_names <- c("q", "t0", p_col_names, weight_params)
  if (length(weight_params)) {
    weights <- unname(param_values[weight_params])
    if (any(!is.finite(weights) | weights < 0) || sum(weights) > 1) {
      stop("Sampled mixture weights must define a probability simplex", call. = FALSE)
    }
  }
  trigger_params <- names(prep$shared_triggers)
  if (length(trigger_params)) {
    triggers <- unname(param_values[trigger_params])
    if (any(!is.finite(triggers) | triggers < 0 | triggers > 1)) {
      stop("Trigger probabilities must lie in [0, 1]", call. = FALSE)
    }
  }

  base_mat <- matrix(NA_real_, nrow = length(accs), ncol = length(col_names))
  colnames(base_mat) <- col_names
  for (i in seq_along(accs)) {
    acc_id <- acc_ids[[i]]
    dist_params <- dist_param_list[[i]]
    p_vals <- numeric(max_p)
    for (j in seq_along(dist_params)) {
      p_vals[[j]] <- param_values[[paste0(acc_id, ".", dist_params[[j]])]]
    }
    trigger <- accs[[i]]$shared_trigger_id
    q <- if (is.null(trigger)) {
      0
    } else {
      param_values[[trigger]]
    }
    base_mat[i, ] <- c(
      q,
      param_values[[paste0(acc_id, ".t0")]],
      p_vals,
      unname(param_values[weight_params])
    )
  }

  if (length(n_trials) != 1L || !is.numeric(n_trials) ||
      !is.finite(n_trials) || n_trials < 1L || n_trials != as.integer(n_trials)) {
    stop("n_trials must be a positive integer", call. = FALSE)
  }
  n_trials <- as.integer(n_trials)

  params <- base_mat[rep(seq_along(accs), times = n_trials), , drop = FALSE]
  class(params) <- c("accumulatr_parameters", "matrix", "array")
  params
}

.expand_accumulator_rows <- function(structure, data) {
  acc_ids <- names(structure$prep$accumulators)
  df <- as.data.frame(data)
  df$trials <- seq_len(nrow(df))
  out <- df[rep(seq_len(nrow(df)), each = length(acc_ids)), , drop = FALSE]
  out$racer <- rep(acc_ids, times = nrow(df))
  rownames(out) <- NULL
  out
}
