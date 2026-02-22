`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

# ------------------------------------------------------------------------------
# Observation metadata helpers
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

.normalize_observation_metadata <- function(metadata = list(), n_outcomes = NULL) {
  metadata <- metadata %||% list()
  obs <- metadata$observation %||% list()
  if (!is.list(obs)) {
    stop("metadata$observation must be a list")
  }
  if (!is.null(obs$mode) && !identical(obs$mode, "top_k")) {
    stop("Only metadata$observation$mode = 'top_k' is currently supported")
  }
  if (!is.null(n_outcomes)) {
    obs$n_outcomes <- n_outcomes
  }
  n_out <- .validate_n_outcomes(obs$n_outcomes %||% 1L)
  metadata$observation <- list(mode = "top_k", n_outcomes = n_out)
  metadata
}

.extract_observation_spec <- function(x) {
  obs <- NULL
  if (is.list(x) && !is.null(x$observation)) {
    obs <- x$observation
  }
  if (is.null(obs) && is.list(x) && !is.null(x$metadata)) {
    obs <- x$metadata$observation
  }
  if (is.null(obs)) {
    return(list(mode = "top_k", n_outcomes = 1L))
  }
  if (!is.list(obs)) {
    stop("observation specification must be a list")
  }
  mode <- obs$mode %||% "top_k"
  if (!identical(mode, "top_k")) {
    stop("Only observation mode 'top_k' is currently supported")
  }
  list(
    mode = "top_k",
    n_outcomes = .validate_n_outcomes(obs$n_outcomes %||% 1L)
  )
}

.is_after_onset <- function(x) {
  is.list(x) && identical(x$kind %||% NULL, "after")
}

.validate_after_literal <- function(x) {
  if (!.is_after_onset(x)) {
    stop("after() onset must be a list with kind = 'after'")
  }
  source <- x$source %||% NULL
  if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
    stop("after(source, ...) requires a single non-empty character source id")
  }
  lag <- x$lag %||% 0
  if (!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag < 0) {
    stop("after(..., lag =) must be a single finite numeric value >= 0")
  }
  list(kind = "after", source = as.character(source), lag = as.numeric(lag))
}

.normalize_onset_literal <- function(onset) {
  if (is.numeric(onset) && length(onset) == 1L && is.finite(onset)) {
    return(as.numeric(onset))
  }
  if (is.character(onset)) {
    stop("Character onset shorthand is not supported; use after(\"id\", lag = ...)")
  }
  if (.is_after_onset(onset)) {
    after_spec <- .validate_after_literal(onset)
    return(structure(after_spec, class = c("race_onset_after", "list")))
  }
  stop("onset must be a single finite numeric value or after(source, lag = ...)")
}

.as_named_pool_defs <- function(pools) {
  pools <- pools %||% list()
  if (length(pools) == 0L) {
    return(list())
  }
  nms <- names(pools)
  if (!is.null(nms) && length(nms) == length(pools) && all(nzchar(nms))) {
    return(pools)
  }
  ids <- vapply(pools, function(pool) pool$id %||% "", character(1))
  if (any(!nzchar(ids))) {
    stop("All pools must define non-empty ids")
  }
  setNames(pools, ids)
}

.pool_threshold <- function(pool_def) {
  members <- pool_def$members %||% character(0)
  k <- pool_def$k %||% pool_def$rule$k %||% length(members)
  if (is.infinite(k)) {
    k <- length(members)
  }
  as.integer(k)
}

.resolve_onset_source_kind <- function(source, acc_ids, pool_ids, outcome_labels) {
  if (source %in% acc_ids) {
    return("accumulator")
  }
  if (source %in% pool_ids) {
    return("pool")
  }
  if (source %in% outcome_labels) {
    stop(sprintf(
      "Onset source '%s' references an outcome label; onset sources must be accumulator or pool ids",
      source
    ))
  }
  stop(sprintf(
    "Unknown onset source '%s'; onset sources must be declared accumulators or pools",
    source
  ))
}

.normalize_onset_spec <- function(onset, acc_ids, pool_ids, outcome_labels) {
  if (is.numeric(onset) && length(onset) == 1L && is.finite(onset)) {
    return(list(kind = "absolute", value = as.numeric(onset)))
  }
  if (is.character(onset)) {
    stop("Character onset shorthand is not supported; use after(\"id\", lag = ...)")
  }
  if (.is_after_onset(onset)) {
    after_spec <- .validate_after_literal(onset)
    source_kind <- .resolve_onset_source_kind(
      source = after_spec$source,
      acc_ids = acc_ids,
      pool_ids = pool_ids,
      outcome_labels = outcome_labels
    )
    return(list(
      kind = "after",
      source = after_spec$source,
      source_kind = source_kind,
      lag = after_spec$lag
    ))
  }
  stop("onset must be a single finite numeric value or after(source, lag = ...)")
}

.expand_pool_accumulator_dependencies <- function(pool_id, pool_defs, acc_ids, stack = character(0)) {
  if (!pool_id %in% names(pool_defs)) {
    stop(sprintf("Unknown pool '%s' while resolving onset dependencies", pool_id))
  }
  if (pool_id %in% stack) {
    cycle <- c(stack, pool_id)
    stop("Pool dependency cycle detected while resolving onsets: ", paste(cycle, collapse = " -> "))
  }
  pool <- pool_defs[[pool_id]]
  members <- as.character(pool$members %||% character(0))
  deps <- character(0)
  for (member in members) {
    if (member %in% acc_ids) {
      deps <- c(deps, member)
    } else if (member %in% names(pool_defs)) {
      deps <- c(
        deps,
        .expand_pool_accumulator_dependencies(member, pool_defs, acc_ids, c(stack, pool_id))
      )
    } else {
      stop(sprintf("Pool '%s' references unknown member '%s'", pool_id, member))
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
    children <- unique(adjacency[[node]] %||% character(0))
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

.build_onset_dependency_metadata <- function(acc_defs, pool_defs) {
  acc_ids <- names(acc_defs %||% list())
  pool_defs <- .as_named_pool_defs(pool_defs %||% list())
  onset_specs <- setNames(vector("list", length(acc_ids)), acc_ids)
  dependencies <- setNames(vector("list", length(acc_ids)), acc_ids)
  onset_sources <- setNames(vector("list", length(acc_ids)), acc_ids)

  for (acc_id in acc_ids) {
    acc_def <- acc_defs[[acc_id]] %||% list()
    spec <- acc_def$onset_spec %||% list(kind = "absolute", value = acc_def$onset %||% 0)
    onset_specs[[acc_id]] <- spec
    if (identical(spec$kind, "after")) {
      src_kind <- spec$source_kind %||% NULL
      src <- spec$source %||% NULL
      if (identical(src_kind, "accumulator")) {
        preds <- src
      } else if (identical(src_kind, "pool")) {
        preds <- .expand_pool_accumulator_dependencies(src, pool_defs, acc_ids)
      } else {
        stop(sprintf("Unsupported onset source kind '%s' for accumulator '%s'", src_kind, acc_id))
      }
      dependencies[[acc_id]] <- unique(as.character(preds))
      onset_sources[[acc_id]] <- list(
        source = src,
        source_kind = src_kind,
        lag = spec$lag %||% 0
      )
      if (acc_id %in% dependencies[[acc_id]]) {
        stop("Chained onset dependency cycle detected: ", acc_id, " -> ", acc_id)
      }
    } else {
      dependencies[[acc_id]] <- character(0)
      onset_sources[[acc_id]] <- NULL
    }
  }

  adjacency <- setNames(vector("list", length(acc_ids)), acc_ids)
  incoming <- setNames(integer(length(acc_ids)), acc_ids)
  for (target in acc_ids) {
    preds <- dependencies[[target]] %||% character(0)
    if (length(preds) == 0L) {
      next
    }
    for (pred in preds) {
      if (!pred %in% acc_ids) {
        stop(sprintf("Unknown accumulator dependency '%s' for onset of '%s'", pred, target))
      }
      adjacency[[pred]] <- c(adjacency[[pred]], target)
      incoming[[target]] <- incoming[[target]] + 1L
    }
  }
  for (node in acc_ids) {
    adjacency[[node]] <- unique(adjacency[[node]] %||% character(0))
  }

  queue <- acc_ids[incoming == 0L]
  topo <- character(0)
  while (length(queue) > 0L) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    topo <- c(topo, node)
    children <- adjacency[[node]] %||% character(0)
    for (child in children) {
      incoming[[child]] <- incoming[[child]] - 1L
      if (incoming[[child]] == 0L) {
        queue <- c(queue, child)
      }
    }
  }
  if (length(topo) != length(acc_ids)) {
    cycle <- .find_onset_cycle_path(acc_ids, adjacency)
    cycle_text <- if (length(cycle) > 0L) paste(cycle, collapse = " -> ") else "unknown cycle"
    stop("Chained onset dependency cycle detected: ", cycle_text)
  }

  has_dependencies <- any(vapply(onset_specs, function(spec) identical(spec$kind, "after"), logical(1)))
  list(
    onset_specs = onset_specs,
    onset_dependencies = dependencies,
    onset_sources = onset_sources,
    onset_topology = topo,
    onset_has_dependencies = has_dependencies
  )
}

.deterministic_source_signature <- function(source, acc_ids, pool_defs, stack = character(0)) {
  if (source %in% acc_ids) {
    return(paste0("acc:", source))
  }
  if (!source %in% names(pool_defs)) {
    return(paste0("id:", source))
  }
  if (source %in% stack) {
    return(paste0("pool_cycle:", paste(c(stack, source), collapse = "->")))
  }
  pool <- pool_defs[[source]]
  members <- as.character(pool$members %||% character(0))
  k <- .pool_threshold(pool)
  if (k == 1L && length(members) == 1L) {
    return(.deterministic_source_signature(members[[1L]], acc_ids, pool_defs, c(stack, source)))
  }
  paste0("pool:", source)
}

.outcome_is_direct_event <- function(outcome_def, acc_ids, pool_ids) {
  expr <- outcome_def$expr %||% list()
  if (!is.list(expr) || !identical(expr$kind, "event")) {
    return(FALSE)
  }
  source <- expr$source %||% NULL
  if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
    return(FALSE)
  }
  if (source %in% c("__GUESS__", "__DEADLINE__")) {
    return(FALSE)
  }
  source %in% c(acc_ids, pool_ids)
}

# v1 validator for multi-outcome readout declarations.
.validate_multi_outcome_dsl <- function(model_or_prep) {
  obs <- .extract_observation_spec(model_or_prep)
  n_outcomes <- obs$n_outcomes
  if (n_outcomes <= 1L) {
    return(invisible(NULL))
  }
  outcomes <- model_or_prep$outcomes %||% list()
  n_defined <- length(outcomes)
  if (n_defined == 0L) {
    stop("n_outcomes > 1 requires declared outcomes")
  }
  if (n_outcomes > n_defined) {
    stop(sprintf(
      "n_outcomes (%d) cannot exceed number of declared outcomes (%d)",
      n_outcomes, n_defined
    ))
  }

  acc_ids <- names(model_or_prep$accumulators %||% list())
  pool_ids <- names(model_or_prep$pools %||% list())
  issues <- character(0)
  direct_sources <- character(0)
  direct_labels <- character(0)
  for (i in seq_along(outcomes)) {
    out <- outcomes[[i]] %||% list()
    label <- out$label %||% names(outcomes)[i] %||% sprintf("outcome_%d", i)
    opts <- out$options %||% list()
    is_direct <- .outcome_is_direct_event(out, acc_ids, pool_ids)
    if (!is_direct) {
      issues <- c(issues, sprintf("%s (must be a direct event outcome)", label))
    } else {
      direct_sources <- c(direct_sources, out$expr$source %||% "")
      direct_labels <- c(direct_labels, label)
    }
    if (!is.null(opts$guess)) {
      issues <- c(issues, sprintf("%s (guess option not supported)", label))
    }
    if (!is.null(opts$alias_of)) {
      issues <- c(issues, sprintf("%s (alias_of option not supported)", label))
    }
    if (!is.null(opts$map_outcome_to)) {
      issues <- c(issues, sprintf("%s (map_outcome_to option not supported)", label))
    }
  }
  if (length(issues) > 0L) {
    stop(
      "n_outcomes > 1 currently supports only direct event outcomes ",
      "with no guess/alias_of/map_outcome_to options. Invalid outcomes: ",
      paste(unique(issues), collapse = ", ")
    )
  }

  pool_defs <- .as_named_pool_defs(model_or_prep$pools %||% list())
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

# ------------------------------------------------------------------------------
# Expression parsing utilities
# ------------------------------------------------------------------------------

.is_expr_node <- function(x) is.list(x) && length(x) > 0 && !is.null(x$kind)

.expr_from_symbol <- function(sym) {
  list(kind = "event", source = as.character(sym), k = NULL)
}

.expr_label <- function(node) {
  if (is.symbol(node)) {
    return(as.character(node))
  }
  if (is.character(node) && length(node) == 1L) {
    return(node)
  }
  stop("Unable to convert expression label of type '", typeof(node), "'")
}

.expr_from_value <- function(val) {
  if (.is_expr_node(val)) {
    return(val)
  }
  if (inherits(val, "formula")) {
    return(.build_expr(val[[length(val)]]))
  }
  if (is.character(val) && length(val) == 1L) {
    token <- trimws(val)
    has_logic <- grepl("[&|!()]", token)
    looks_call <- grepl("\\w+\\s*\\(", token)
    if (!has_logic && !looks_call) {
      return(list(kind = "event", source = token, k = NULL))
    }
    parsed <- tryCatch(str2lang(token), error = function(e) NULL)
    if (is.null(parsed)) {
      return(list(kind = "event", source = token, k = NULL))
    }
    return(.build_expr(parsed))
  }
  if (is.symbol(val)) {
    return(.expr_from_symbol(val))
  }
  if (is.call(val)) {
    return(.parse_expr_call(val))
  }
  stop("Cannot interpret expression component of type '", typeof(val), "'")
}

.parse_event_call <- function(call) {
  args <- as.list(call)[-1]
  nm <- names(call)[-1]
  id_arg <- NULL
  k_arg <- NULL

  if (length(args) >= 1L) {
    arg_name <- if (length(nm) >= 1L) nm[[1L]] else ""
    if (is.null(arg_name) || arg_name == "" || identical(arg_name, "id")) {
      id_arg <- args[[1L]]
    }
  }
  if (length(args) >= 2L) {
    arg_name <- if (length(nm) >= 2L) nm[[2L]] else ""
    if (is.null(arg_name) || arg_name == "" || identical(arg_name, "k")) {
      k_arg <- args[[2L]]
    }
  }
  if ("id" %in% nm) id_arg <- args[[which(nm == "id")[1L]]]
  if ("k" %in% nm) k_arg <- args[[which(nm == "k")[1L]]]

  if (is.null(id_arg)) stop("event() requires an id argument")
  source <- .expr_label(id_arg)
  k <- NULL
  if (!is.null(k_arg)) {
    if (!is.numeric(k_arg) || length(k_arg) != 1L) {
      stop("event(k =) must be a single numeric value")
    }
    k <- as.integer(k_arg)
  }
  list(kind = "event", source = source, k = k)
}



.parse_guard_call <- function(call) {
  args <- as.list(call)[-1]
  nm <- names(call)[-1]

  blocker <- reference <- NULL
  if (length(args) >= 1L && (is.null(nm[[1L]]) || nm[[1L]] == "")) blocker <- args[[1L]]
  if (length(args) >= 2L && (is.null(nm[[2L]]) || nm[[2L]] == "")) reference <- args[[2L]]
  if ("blocker" %in% nm) blocker <- args[[which(nm == "blocker")[1L]]]
  if ("reference" %in% nm) reference <- args[[which(nm == "reference")[1L]]]
  if (is.null(blocker) || is.null(reference)) {
    stop("guard() requires both blocker and reference arguments")
  }
  list(
    kind = "guard",
    blocker = .expr_from_value(blocker),
    reference = .expr_from_value(reference)
  )
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
  if (op %in% c("!", "~")) {
    return(list(kind = "not", arg = .expr_from_value(call[[2]])))
  }
  if (identical(op, "event")) {
    return(.parse_event_call(call))
  }
  if (identical(op, "guard")) {
    return(.parse_guard_call(call))
  }
  stop(sprintf("Unsupported token '%s' in expression", op))
}

.build_expr <- function(expr) {
  if (.is_expr_node(expr)) {
    return(expr)
  }
  .expr_from_value(expr)
}

#' Normalize an outcome expression
#'
#' @param expr Expression or symbol describing an event/guard
#' @return Expression list used in model specs
#' @examples
#' build_outcome_expr(quote(A & !B))
#' @export
build_outcome_expr <- function(expr) {
  .build_expr(expr)
}

# ------------------------------------------------------------------------------
# Public DSL helpers
# ------------------------------------------------------------------------------

#' Define an inhibitory relationship
#'
#' @param reference Outcome expression or label to inhibit
#' @param by Blocking expression or label
#' @return Guard expression list
#' @examples
#' inhibit("A", "B")
#' @export
inhibit <- function(reference, by) {
  ref_expr <- if (is.character(reference) && length(reference) == 1L) {
    list(kind = "event", source = reference, k = NULL)
  } else {
    .build_expr(reference)
  }
  blocker_expr <- if (is.character(by) && length(by) == 1L) {
    list(kind = "event", source = by, k = NULL)
  } else {
    .build_expr(by)
  }
  list(kind = "guard", blocker = blocker_expr, reference = ref_expr)
}

#' First-of expression helper
#'
#' @param ... Expressions or labels to combine with OR
#' @return Expression list
#' @examples
#' first_of("A", "B")
#' @export
first_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("first_of() requires at least one argument")
  if (length(args) == 1 && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }
  list(kind = "or", args = lapply(args, .expr_from_value))
}

#' All-of expression helper
#'
#' @param ... Expressions or labels to combine with AND
#' @return Expression list
#' @examples
#' all_of("A", "B")
#' @export
all_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("all_of() requires at least one argument")
  if (length(args) == 1 && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }
  list(kind = "and", args = lapply(args, .expr_from_value))
}

#' None-of expression helper
#'
#' @param expr Expression or label to negate
#' @return Expression list
#' @export
#' @aliases exclude
#' @examples
#' none_of("A")
none_of <- function(expr) {
  inner <- if (is.character(expr) && length(expr) == 1L) {
    list(kind = "event", source = expr, k = NULL)
  } else {
    .build_expr(expr)
  }
  list(kind = "not", arg = inner)
}

#' @rdname none_of
#' @export
exclude <- none_of

#' Define a chained onset relative to another source
#'
#' @param source Accumulator or pool id that must finish first.
#' @param lag Non-negative delay added after source completion.
#' @return A chained onset specification for `add_accumulator(onset = ...)`.
#' @examples
#' after("A")
#' after("pool1", lag = 0.05)
#' @export
after <- function(source, lag = 0) {
  if (!is.character(source) || length(source) != 1L || !nzchar(source)) {
    stop("after(source, ...) requires a single non-empty character source id")
  }
  if (!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag < 0) {
    stop("after(..., lag =) must be a single finite numeric value >= 0")
  }
  structure(
    list(kind = "after", source = as.character(source), lag = as.numeric(lag)),
    class = c("race_onset_after", "list")
  )
}

# ------------------------------------------------------------------------------
# Model builder
# ------------------------------------------------------------------------------

#' Create an empty race specification
#'
#' @param n_outcomes Number of observed outcomes to retain per trial.
#' @return A race_spec object
#' @examples
#' race_spec()
#' @export
race_spec <- function(n_outcomes = 1L) {
  metadata <- .normalize_observation_metadata(list(), n_outcomes = n_outcomes)
  structure(list(
    accumulators = list(),
    pools = list(),
    outcomes = list(),
    triggers = list(),
    shared_params = list(),
    components = list(),
    mixture_options = list(),
    metadata = metadata
  ), class = "race_spec")
}

.collect_params <- function(params, dots) {
  if (!is.null(params) && length(dots) > 0) {
    stop("Provide either params=list(...) or named arguments, not both")
  }
  if (!is.null(params)) {
    if (!is.list(params)) stop("params must be a list")
    return(params)
  }
  dots
}

.ensure_acc_param_t0 <- function(params) {
  if (is.null(params) || length(params) == 0L) params <- list()
  if (is.null(params$t0) || length(params$t0) == 0L) {
    params$t0 <- 0
  } else {
    params$t0 <- as.numeric(params$t0)[1]
  }
  params
}

#' Add an accumulator definition
#'
#' @param spec race_spec object
#' @param id Accumulator id
#' @param dist Distribution name
#' @param onset Onset shift; either numeric or `after(source, lag = ...)`.
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' @export
add_accumulator <- function(spec, id, dist, onset = 0) {
  stopifnot(inherits(spec, "race_spec"))
  onset <- .normalize_onset_literal(onset)
  spec$accumulators[[length(spec$accumulators) + 1L]] <- list(
    id = id,
    dist = dist,
    onset = onset
  )
  spec
}

#' Add a pool definition
#'
#' @param spec race_spec object
#' @param id Pool id
#' @param members Member ids
#' @param k Threshold parameter
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_pool(spec, "P1", members = c("A", "B"), k = 1L)
#' @export
add_pool <- function(spec, id, members, k = 1L) {
  stopifnot(inherits(spec, "race_spec"))
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

#' Add an outcome definition
#'
#' @param spec race_spec object
#' @param label Outcome label
#' @param expr Outcome expression
#' @param options Optional list of options
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_outcome(spec, "A_win", "A")
#' @export
add_outcome <- function(spec, label, expr, options = list()) {
  stopifnot(inherits(spec, "race_spec"))
  spec$outcomes[[length(spec$outcomes) + 1L]] <- list(
    label = label,
    expr = build_outcome_expr(expr),
    options = options %||% list()
  )
  spec
}

#' Add a group of members
#'
#' @param spec race_spec object
#' @param id Group id
#' Define a component (unified definition and membership)
#'
#' @param spec race_spec object
#' @param id Component id (e.g. "fast")
#' @param members Character vector of accumulator IDs belonging to this component.
#' @param weight Optional fixed weight (0-1).
#' @param weight_param Optional parameter name for weight.
#' @param n_outcomes Optional component-level override for observed outcomes.
#' @param attrs Additional attributes.
#' @return Updated race_spec
#' @export
add_component <- function(spec, id, members, weight = NULL, weight_param = NULL, n_outcomes = NULL, attrs = list()) {
  stopifnot(inherits(spec, "race_spec"))
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Component must specify members")
  }
  comp_attrs <- attrs %||% list()
  if (!is.list(comp_attrs)) {
    stop("Component attrs must be a list")
  }
  if (!is.null(comp_attrs$n_outcomes)) {
    comp_attrs$n_outcomes <- .validate_n_outcomes(comp_attrs$n_outcomes)
  }
  if (!is.null(n_outcomes)) {
    comp_attrs$n_outcomes <- .validate_n_outcomes(n_outcomes)
  }
  spec$components[[length(spec$components) + 1L]] <- list(
    id = id,
    members = as.character(members),
    weight = weight,
    weight_param = weight_param,
    attrs = comp_attrs
  )
  spec
}

#' Add a shared trigger (gate)
#'
#' @param spec race_spec object
#' @param id Trigger id (optional; defaults to first member name if missing)
#' @param members Character vector of accumulator IDs sharing this trigger.
#' @param q Probability (0-1) for trigger failure.
#' @param param Optional parameter name supplying q.
#' @param draw Either "shared" (one draw for all members) or "independent" (same q, independent draws).
#' @return Updated race_spec
#' @export
add_trigger <- function(spec, id, members, q = NULL, param = NULL, draw = c("shared", "independent")) {
  stopifnot(inherits(spec, "race_spec"))
  draw <- match.arg(draw)
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Trigger must specify members")
  }
  if (is.null(id)) {
    id <- members[[1]]
  }
  spec$triggers[[length(spec$triggers) + 1L]] <- list(
    id = id,
    members = as.character(members),
    q = q,
    param = param,
    draw = draw
  )
  if (!is.null(q)) {
    for (m in members) {
      for (i in seq_along(spec$accumulators)) {
        if (identical(spec$accumulators[[i]]$id, m)) {
          spec$accumulators[[i]]$q <- q
          break
        }
      }
    }
  }
  spec
}

#' Add shared parameters for a set of accumulators
#'
#' @param spec race_spec object
#' @param members Character vector of accumulator IDs.
#' @param ... Named parameters to share (e.g., meanlog = "log(0.5)", or explicit value).
#' @return Updated race_spec
#' @export
add_shared_params <- function(spec, members, ...) {
  stopifnot(inherits(spec, "race_spec"))
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Shared params must specify members")
  }
  params <- list(...)
  if (length(params) == 0) {
    warning("add_shared_params called with no parameters")
    return(spec)
  }
  spec$shared_params[[length(spec$shared_params) + 1L]] <- list(
    members = as.character(members),
    params = params
  )
  spec
}

#' Set mixture options (global)
#'
#' @param spec race_spec object
#' @param mode Mixture mode (e.g. "sample")
#' @param reference Component ID to use as reference (base)
#' @return Updated race_spec
#' @export
set_mixture_options <- function(spec, mode = "sample", reference = NULL) {
  stopifnot(inherits(spec, "race_spec"))
  spec$mixture_options <- list(
    mode = mode,
    reference = reference
  )
  spec
}

#' Set metadata on a specification
#'
#' @param spec race_spec object
#' @param ... Named metadata entries
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- set_metadata(spec, rel_tol = 1e-4)
#' @export
set_metadata <- function(spec, ...) {
  stopifnot(inherits(spec, "race_spec"))
  updates <- list(...)
  for (nm in names(updates)) {
    spec$metadata[[nm]] <- updates[[nm]]
  }
  spec
}



#' Build an expression guard node
#'
#' @param blocker Blocker expression or id
#' @param reference Reference expression or id
#' @return Guard expression list
#' @examples
#' expr_guard("B", "A")
#' @export
expr_guard <- function(blocker, reference) {
  list(kind = "guard", blocker = blocker, reference = reference)
}

clone_obj <- function(x) {
  if (is.environment(x)) stop("Cannot clone environments for race models")
  if (is.list(x)) {
    out <- lapply(x, clone_obj)
    attrs <- attributes(x)
    if (!is.null(attrs)) attributes(out) <- lapply(attrs, clone_obj)
    return(out)
  }
  x
}

outcome_def <- function(label, expr, options = list()) {
  list(
    label = label,
    expr = build_outcome_expr(expr),
    options = options %||% list()
  )
}



race_model <- function(accumulators, pools = list(), outcomes, triggers = list(), shared_params = list(), components = list(), mixture_options = list(), metadata = list()) {
  if (inherits(accumulators, "race_spec")) {
    if (!missing(pools) || !missing(outcomes) || !missing(triggers) || !missing(shared_params) || !missing(components) || !missing(mixture_options) || !missing(metadata)) {
      stop("When passing a race_spec, provide it alone and call race_model(spec)")
    }
    finalized <- finalize_model(accumulators)
    return(finalized$model_spec)
  }
  if (inherits(accumulators, "model_structure") || inherits(accumulators, "generator_structure")) {
    spec <- accumulators$model_spec %||% NULL
    if (is.null(spec)) stop("model_structure is missing model_spec")
    return(spec)
  }
  if (inherits(accumulators, "race_model_spec")) {
    accumulators$metadata <- .normalize_observation_metadata(accumulators$metadata %||% list())
    return(accumulators)
  }
  if (missing(accumulators) || missing(outcomes)) {
    stop("accumulators and outcomes must be supplied")
  }
  structure(list(
    accumulators = lapply(accumulators, clone_obj),
    pools = lapply(pools %||% list(), clone_obj),
    outcomes = lapply(outcomes, function(out) {
      out$expr <- build_outcome_expr(out$expr)
      out$options <- out$options %||% list()
      out
    }),
    triggers = lapply(triggers %||% list(), clone_obj),
    shared_params = lapply(shared_params %||% list(), clone_obj),
    components = lapply(components %||% list(), clone_obj),
    mixture_options = mixture_options %||% list(),
    metadata = .normalize_observation_metadata(metadata %||% list())
  ), class = "race_model_spec")
}

# ----------------------------------------------------------------------
# Model normalization/finalization (shared by simulation and likelihood)
# ----------------------------------------------------------------------

.normalize_model <- function(model) {
  if (inherits(model, "race_spec")) {
    metadata <- .normalize_observation_metadata(model$metadata %||% list())
    return(structure(list(
      accumulators = unname(model$accumulators),
      pools = unname(model$pools),
      outcomes = unname(model$outcomes),
      triggers = unname(model$triggers),
      shared_params = unname(model$shared_params),
      components = unname(model$components),
      mixture_options = model$mixture_options %||% list(),
      metadata = metadata
    ), class = "race_model_spec"))
  }
  if (inherits(model, "race_model_spec")) {
    model$metadata <- .normalize_observation_metadata(model$metadata %||% list())
    return(model)
  }
  model
}

.prepare_acc_defs <- function(model) {
  accs <- model$accumulators
  if (length(accs) == 0) stop("Model must define at least one accumulator")

  acc_ids <- vapply(accs, `[[`, character(1), "id")
  pool_ids <- vapply(model$pools %||% list(), function(pool) pool$id %||% "", character(1))
  pool_ids <- pool_ids[nzchar(pool_ids)]
  outcome_labels <- vapply(model$outcomes %||% list(), function(out) out$label %||% "", character(1))
  outcome_labels <- outcome_labels[nzchar(outcome_labels)]

  defs <- setNames(vector("list", length(accs)), acc_ids)
  for (acc in accs) {
    onset_spec <- .normalize_onset_spec(
      onset = acc$onset %||% 0,
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
      q = acc$q %||% 0,
      params = .ensure_acc_param_t0(acc$params %||% list()),
      components = character(0),
      shared_trigger_id = NULL,
      shared_trigger_q = NULL
    )
  }

  # Process shared parameters
  if (length(model$shared_params) > 0) {
    for (sp in model$shared_params) {
      members <- sp$members
      par_updates <- sp$params
      if (length(members) == 0 || length(par_updates) == 0) next
      for (m in members) {
        if (!is.null(defs[[m]])) {
          defs[[m]]$params <- modifyList(defs[[m]]$params, par_updates, keep.null = TRUE)
        }
      }
    }
  }

  # Process components (membership)
  if (length(model$components) > 0) {
    for (cmp in model$components) {
      comp_id <- cmp$id
      members <- cmp$members
      if (length(members) == 0) next
      for (m in members) {
        if (!is.null(defs[[m]])) {
          defs[[m]]$components <- unique(c(defs[[m]]$components, comp_id))
        }
      }
    }
  }

  # Process triggers
  shared_triggers <- list()
  if (length(model$triggers) > 0) {
    for (trig in model$triggers) {
      members <- trig$members
      if (length(members) == 0) next
      trig_id <- trig$id
      draw_mode <- trig$draw %||% "shared"
      default_q <- trig$q %||% NA_real_
      param_name <- trig$param %||% NULL

      # Validate default_q if NA
      if (is.na(default_q)) {
        q_vals <- vapply(members, function(m) {
          acc_def <- defs[[m]]
          if (is.null(acc_def)) {
            return(NA_real_)
          }
          acc_def$q %||% NA_real_
        }, numeric(1))
        q_vals <- q_vals[!is.na(q_vals)]
        if (length(q_vals) == 0) {
          default_q <- 0
        } else if (length(unique(q_vals)) == 1) {
          default_q <- q_vals[[1]]
        } else {
          stop(sprintf("Trigger '%s' requires a single q value; found disparate values", trig_id))
        }
      }

      if (identical(draw_mode, "shared")) {
        shared_triggers[[trig_id]] <- list(
          id = trig_id,
          group_id = trig_id, # for backwards compat internally, just use trig_id
          members = members,
          q = as.numeric(default_q),
          param = param_name,
          draw = draw_mode
        )
        for (m in members) {
          if (!is.null(defs[[m]])) {
            defs[[m]]$shared_trigger_id <- trig_id
            defs[[m]]$shared_trigger_q <- as.numeric(default_q)
            defs[[m]]$q <- as.numeric(default_q)
          }
        }
      } else if (identical(draw_mode, "independent")) {
        for (m in members) {
          if (!is.null(defs[[m]])) {
            defs[[m]]$q <- as.numeric(default_q)
          }
        }
      }
    }
  }

  for (acc_id in names(defs)) {
    defs[[acc_id]]$params <- .ensure_acc_param_t0(defs[[acc_id]]$params)
  }
  list(acc = defs, shared_triggers = shared_triggers)
}

.extract_components <- function(model) {
  # Components are now stored directly in model$components
  comps <- model$components %||% list()

  # Mixture options stored in model$mixture_options
  mix_opts <- model$mixture_options %||% list()
  mode <- mix_opts$mode %||% "fixed"
  reference <- mix_opts$reference %||% NA_character_
  if (length(comps) == 0) {
    return(list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list(guess = NULL, weight_param = NULL)),
      has_weight_param = FALSE,
      mode = mode,
      reference = "__default__"
    ))
  }
  ids <- vapply(comps, `[[`, character(1), "id")
  attrs <- setNames(vector("list", length(ids)), ids)
  has_wparam <- logical(length(ids))
  weights <- vapply(comps, function(cmp) {
    w <- cmp$weight %||% 1
    w <- as.numeric(w)[1]
    if (!is.finite(w) || w < 0) NA_real_ else w
  }, numeric(1))
  weight_params <- vapply(comps, function(cmp) {
    attrs_cmp <- cmp$attrs %||% list()
    cmp$weight_param %||% attrs_cmp$weight_param %||% NA_character_
  }, character(1))

  if (identical(mode, "sample")) {
    if (is.na(reference) || !nzchar(reference)) {
      # choose first component without a weight_param, else first component
      idx_no_param <- which(is.na(weight_params) | !nzchar(weight_params))
      ref_idx <- if (length(idx_no_param) > 0) idx_no_param[[length(idx_no_param)]] else 1L
      reference <- ids[[ref_idx]]
      message(sprintf("mixture reference not provided; using '%s' as reference component", reference))
    } else if (!reference %in% ids) {
      stop("mixture reference '", reference, "' must match a component id")
    }
    for (i in seq_along(ids)) {
      cmp_id <- ids[[i]]
      attrs_cmp <- comps[[i]]$attrs %||% list()
      if (!is.list(attrs_cmp)) {
        stop(sprintf("Component '%s' attrs must be a list", cmp_id))
      }
      wp <- comps[[i]]$weight_param %||% attrs_cmp$weight_param %||% NULL
      if (!identical(cmp_id, reference) && (is.null(wp) || !nzchar(wp))) {
        stop(sprintf("Component '%s' must define weight_param in sampled mixtures (reference is '%s')", cmp_id, reference))
      }
      has_wparam[[i]] <- !is.null(wp) && nzchar(wp)
      n_out <- attrs_cmp$n_outcomes %||% NULL
      if (!is.null(n_out)) {
        n_out <- .validate_n_outcomes(n_out)
      }
      attrs_cmp$guess <- attrs_cmp$guess %||% NULL
      attrs_cmp$weight_param <- wp
      attrs_cmp$n_outcomes <- n_out
      attrs[[cmp_id]] <- attrs_cmp
    }
  } else {
    # fixed mode: weights only used for defaults
    if (all(is.na(weights))) {
      weights <- rep(1, length(weights))
    } else {
      weights[is.na(weights)] <- 0
      if (sum(weights) <= 0) {
        weights <- rep(1, length(weights))
      }
    }
    weights <- weights / sum(weights)
    for (i in seq_along(ids)) {
      cmp_id <- ids[[i]]
      attrs_cmp <- comps[[i]]$attrs %||% list()
      if (!is.list(attrs_cmp)) {
        stop(sprintf("Component '%s' attrs must be a list", cmp_id))
      }
      wp <- comps[[i]]$weight_param %||% attrs_cmp$weight_param %||% NULL
      has_wparam[[i]] <- !is.null(wp) && nzchar(wp)
      n_out <- attrs_cmp$n_outcomes %||% NULL
      if (!is.null(n_out)) {
        n_out <- .validate_n_outcomes(n_out)
      }
      attrs_cmp$guess <- attrs_cmp$guess %||% NULL
      attrs_cmp$weight_param <- wp
      attrs_cmp$n_outcomes <- n_out
      attrs[[cmp_id]] <- attrs_cmp
    }
    if (is.na(reference) || !nzchar(reference)) reference <- ids[[1]]
  }

  list(
    ids = ids,
    weights = weights,
    attrs = attrs,
    has_weight_param = has_wparam,
    mode = mode,
    reference = reference
  )
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
      k = k
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
  observation <- .extract_observation_spec(model)
  component_attrs <- component_defs$attrs %||% list()
  component_overrides <- integer(0)
  if (length(component_attrs) > 0L) {
    component_overrides <- vapply(component_attrs, function(attr) {
      if (!is.list(attr) || is.null(attr$n_outcomes)) {
        return(NA_integer_)
      }
      .validate_n_outcomes(attr$n_outcomes)
    }, integer(1))
  }
  max_component_outcomes <- 1L
  if (length(component_overrides) > 0L && any(!is.na(component_overrides))) {
    max_component_outcomes <- max(component_overrides, na.rm = TRUE)
  }
  onset_metadata <- .build_onset_dependency_metadata(acc_defs, pool_defs)
  observation$global_n_outcomes <- observation$n_outcomes
  observation$component_n_outcomes <- as.list(component_overrides[!is.na(component_overrides)])
  observation$n_outcomes <- max(observation$n_outcomes, max_component_outcomes)
  prep <- list(
    accumulators = acc_defs,
    pools = pool_defs,
    outcomes = outcome_defs,
    components = component_defs,
    observation = observation,
    special_outcomes = model$metadata$special_outcomes %||% list(),
    shared_triggers = acc_prep$shared_triggers %||% list(),
    onset_specs = onset_metadata$onset_specs,
    onset_dependencies = onset_metadata$onset_dependencies,
    onset_sources = onset_metadata$onset_sources,
    onset_topology = onset_metadata$onset_topology,
    onset_has_dependencies = onset_metadata$onset_has_dependencies
  )
  .validate_multi_outcome_dsl(prep)
  prep
}

.build_component_table <- function(comp_defs) {
  comp_ids <- comp_defs$ids
  comp_tbl <- data.frame(
    component_id = comp_ids,
    weight = comp_defs$weights,
    stringsAsFactors = FALSE
  )
  comp_tbl$has_weight_param <- comp_defs$has_weight_param
  comp_tbl$attrs <- I(comp_defs$attrs[match(comp_ids, names(comp_defs$attrs))])
  comp_tbl$mode <- comp_defs$mode %||% "fixed"
  comp_tbl$reference <- comp_defs$reference %||% comp_ids[[1]]
  comp_tbl
}

.build_accumulator_template <- function(acc_defs) {
  acc_ids <- names(acc_defs)
  if (length(acc_ids) == 0L) {
    return(data.frame(
      dist = character(0),
      onset = numeric(0),
      q = numeric(0),
      shared_trigger_id = character(0),
      shared_trigger_q = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  acc_df <- data.frame(
    dist = vapply(acc_defs, function(acc) acc$dist %||% NA_character_, character(1)),
    onset = vapply(acc_defs, function(acc) acc$onset %||% 0, numeric(1)),
    q = vapply(acc_defs, function(acc) acc$q %||% 0, numeric(1)),
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
#' spec <- add_accumulator(spec, "A", "lognormal")
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

  model_norm <- unwrap_model_spec(model)
  model_norm <- .normalize_model(model_norm)
  model_norm <- unwrap_model_spec(model_norm)
  if (!inherits(model_norm, "race_model_spec")) {
    stop("finalize_model requires a race_model_spec")
  }

  model_norm <- unserialize(serialize(model_norm, NULL))
  model_norm <- .normalize_model(model_norm)
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
  if (inherits(x, "generator_structure")) {
    return(x)
  }
  if (inherits(x, "model_structure")) {
    return(x)
  }
  if (is.list(x) && !is.null(x$prep) && !is.null(x$accumulators) && !is.null(x$components)) {
    class(x) <- unique(c("model_structure", "generator_structure", class(x)))
    return(x)
  }
  finalize_model(x)
}

# ------------------------------------------------------------------------------
# Parameter utilities
# ------------------------------------------------------------------------------

dist_param_names <- function(dist) {
  dist <- tolower(dist)
  if (dist == "lognormal") {
    return(c("meanlog", "sdlog"))
  }
  if (dist == "gamma") {
    return(c("shape", "rate"))
  }
  if (dist == "exgauss") {
    return(c("mu", "sigma", "tau"))
  }
  character(0)
}

#' List sampled parameter names for a model
#'
#' @param model race_spec or race_model_spec
#' @return Character vector of parameter names
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' sampled_pars(spec)
#' @export
sampled_pars <- function(model) {
  spec <- race_model(model)
  params <- character(0)
  for (acc in spec$accumulators) {
    base <- dist_param_names(acc$dist)
    # onset is a fixed accumulator attribute; only q and t0 are estimated
    acc_params <- c(base, "q", "t0")
    params <- c(params, paste0(acc$id, ".", acc_params))
  }
  for (grp in spec$groups %||% list()) {
    shared <- grp$attrs$shared_params
    if (is.null(shared)) next
    shared_names <- unlist(shared, use.names = FALSE)
    shared_names <- as.character(shared_names)
    params <- c(params, paste0(grp$id, ".", shared_names))
  }
  comps <- spec$metadata$mixture$components %||% list()
  for (comp in comps) {
    weight_param <- comp$attrs$weight_param %||% NULL
    if (!is.null(weight_param) && nzchar(weight_param)) {
      params <- c(params, weight_param)
    }
  }
  unique(params)
}

param_table <- function(model, ...) {
  values <- c(...)
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("Named parameter values are required")
  }
  params <- sampled_pars(model)
  missing <- setdiff(params, names(values))
  if (length(missing) > 0L) {
    stop("Missing parameter values for: ", paste(missing, collapse = ", "))
  }
  data.frame(
    parameter = params,
    value = as.numeric(values[params]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build a parameter data frame for likelihood/simulation
#'
#' @param model Model definition
#' @param param_values Named numeric vector of parameter values
#' @param n_trials Number of trial rows to generate
#' @param component Optional component label(s)
#' @param layout Optional layout for slimmed down construction
#' @return Data frame of parameters per trial
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal")
#' spec <- add_outcome(spec, "A_win", "A")
#' vals <- c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0)
#' build_param_matrix(spec, vals, n_trials = 2)
#' @export
build_param_matrix <- function(model,
                               param_values,
                               n_trials = 1L,
                               component = NULL,
                               layout = NULL) {
  spec <- race_model(model)
  accs <- spec$accumulators %||% list()
  if (length(accs) == 0L) {
    stop("Model must define accumulators before building parameter matrices")
  }
  if (is.null(param_values) || length(param_values) == 0L) {
    stop("Named parameter values are required")
  }
  if (is.null(names(param_values)) || any(!nzchar(names(param_values)))) {
    stop("Parameter values must be a named vector")
  }
  if (!is.numeric(n_trials) || length(n_trials) != 1L || n_trials < 1) {
    stop("n_trials must be a positive integer")
  }
  n_trials <- as.integer(n_trials)

  param_names <- names(param_values)
  param_values <- as.numeric(param_values)
  names(param_values) <- param_names

  acc_ids <- vapply(accs, `[[`, character(1), "id")
  acc_dists <- vapply(accs, `[[`, character(1), "dist")
  dist_param_list <- lapply(acc_dists, dist_param_names)
  max_p <- max(vapply(dist_param_list, length, integer(1)), 0L)
  max_p <- max(3L, max_p)
  col_names <- c("q", "w", "t0", paste0("p", seq_len(max_p)))

  # Identify mixture weight parameters to include as columns
  mix <- spec$mixture_options %||% spec$metadata$mixture %||% list()
  comp_defs <- spec$components %||% mix$components %||% list()
  weight_params <- unique(vapply(comp_defs, function(c) {
    p <- c$weight_param %||% c$attrs$weight_param %||% NA_character_
    if (!is.na(p) && nzchar(p)) p else NA_character_
  }, character(1)))
  weight_params <- weight_params[!is.na(weight_params)]
  col_names <- c(col_names, weight_params)

  # Shared parameter groups (R-side only): expand group-level values to member accumulators.
  shared_param_values <- list()
  # Legacy groups
  for (grp in spec$groups %||% list()) {
    attrs <- grp$attrs %||% list()
    shared <- attrs$shared_params %||% NULL
    if (is.null(shared)) next
    shared <- unique(as.character(unlist(shared)))
    shared <- shared[nzchar(shared)]
    if (length(shared) == 0L) next
    members <- grp$members %||% character(0)
    if (length(members) == 0L) next
    for (par_name in shared) {
      group_nm <- paste0(grp$id, ".", par_name)
      if (!group_nm %in% names(param_values)) next
      value <- as.numeric(param_values[[group_nm]])
      for (m in members) {
        if (!m %in% acc_ids) next
        leaf_nm <- paste0(m, ".", par_name)
        if (leaf_nm %in% names(param_values)) next
        if (!is.null(shared_param_values[[leaf_nm]]) &&
          !isTRUE(all.equal(shared_param_values[[leaf_nm]], value))) {
          stop(sprintf("Conflicting shared parameter values for '%s'", leaf_nm))
        }
        shared_param_values[[leaf_nm]] <- value
      }
    }
  }
  # New shared_params
  for (sp in spec$shared_params %||% list()) {
    members <- sp$members
    params <- sp$params

    # Handle unnamed arguments (interpreted as param_name = param_name)
    param_names <- names(params)
    if (is.null(param_names)) param_names <- rep("", length(params))

    for (i in seq_along(params)) {
      pname <- param_names[i]
      pvalue <- params[[i]]

      # Determine accumulator side name and shared side name
      if (nzchar(pname)) {
        # named: meanlog = "shared_m" -> acc sets meanlog to shared_m
        target_param <- pname
        shared_name <- pvalue
      } else {
        # unnamed: "meanlog" -> acc sets meanlog to meanlog
        if (!is.character(pvalue) || length(pvalue) != 1) {
          stop("Unnamed shared param arguments must be parameter name strings")
        }
        target_param <- pvalue
        shared_name <- pvalue
      }

      # Now propagate constraint:
      # For each member, member.target_param MUST equal member.target_param (if shared)
      # OR: we create a mapping where 'member.target_param' looks up 'shared_name' from inputs?

      # Logic: The input `param_values` (passed to build_param_matrix) contains the RAW values.
      # If we are sharing, we expect the USER to provide `shared_name` in the input params (e.g. "meanlog")
      # And we must ensure that for all members, we use that value.

      # Check if shared_name is in param_values
      if (shared_name %in% names(param_values)) {
        val <- as.numeric(param_values[[shared_name]])

        for (m in members) {
          # Check if this accumulator even needs this param?
          # We blindly assign/overwrite.
          # But wait, build_param_matrix constructs the matrix row from `base_params` + overrides.
          # We need to inject 'val' into the correct column for this accumulator.

          fullname <- paste0(m, ".", target_param)

          # Check consistency if already set?
          if (fullname %in% names(shared_param_values) &&
            !isTRUE(all.equal(shared_param_values[[fullname]], val))) {
            stop(sprintf("Conflicting shared parameter values for '%s'", fullname))
          }
          shared_param_values[[fullname]] <- val
        }
      }
    }
  }

  if (length(shared_param_values) > 0L) {
    expanded <- unlist(shared_param_values, use.names = TRUE)
    param_values <- c(param_values, expanded[setdiff(names(expanded), names(param_values))])
  }

  # Map accumulators to components from group/component declarations.
  comp_attr <- list()
  for (grp in spec$groups %||% list()) {
    attrs <- grp$attrs %||% list()
    comp_val <- attrs$component %||% NULL
    if (is.null(comp_val)) next
    members <- grp$members %||% character(0)
    for (m in members) {
      comp_attr[[m]] <- comp_val
    }
  }
  for (cmp in spec$components %||% list()) {
    cid <- cmp$id
    members <- cmp$members %||% character(0)
    for (m in members) {
      comp_attr[[m]] <- cid
    }
  }

  # Mixture component weights (mapped to components, then to member rows)
  mix <- spec$mixture_options %||% spec$metadata$mixture %||% list()
  comp_defs <- spec$components %||% mix$components %||% list()
  comp_ids <- vapply(comp_defs, `[[`, character(1), "id")
  comp_weights <- setNames(rep(NA_real_, length(comp_defs)), comp_ids)
  comp_mode <- mix$mode %||% "fixed"
  comp_ref <- mix$reference %||% if (length(comp_ids) > 0) comp_ids[[1]] else NA_character_
  comp_index <- setNames(seq_along(comp_ids), comp_ids)
  if (length(comp_defs) > 0) {
    for (i in seq_along(comp_defs)) {
      comp <- comp_defs[[i]]
      attrs <- comp$attrs %||% list()
      wp <- comp$weight_param %||% attrs$weight_param %||% NULL
      if (!is.null(wp) && wp %in% names(param_values)) {
        comp_weights[[i]] <- as.numeric(param_values[[wp]])
      } else if (!is.null(comp$weight)) {
        comp_weights[[i]] <- as.numeric(comp$weight)
      }
    }
    if (identical(comp_mode, "sample")) {
      non_ref_ids <- setdiff(comp_ids, comp_ref)
      non_ref_sum <- sum(comp_weights[non_ref_ids], na.rm = TRUE)
      if (is.na(comp_weights[[comp_ref]])) {
        ref_val <- 1 - non_ref_sum
        if (!is.finite(ref_val) || ref_val < 0) ref_val <- 0
        comp_weights[[comp_ref]] <- ref_val
      }
    }
    if (any(is.na(comp_weights))) {
      if (length(comp_weights) > 0) {
        comp_weights[is.na(comp_weights)] <- 1 / length(comp_weights)
      }
    }
    total_w <- sum(comp_weights)
    if (is.finite(total_w) && total_w > 0) {
      comp_weights <- comp_weights / total_w
    }
  }

  # Shared trigger groups (q may be shared or independent draws)
  acc_trigger_map <- setNames(vector("list", length(acc_ids)), acc_ids)
  all_trigger_defs <- c(
    lapply(spec$groups %||% list(), function(g) {
      t <- g$attrs$shared_trigger
      if (!is.null(t)) {
        list(
          id = t$id %||% g$id,
          param = t$param,
          q = t$q,
          draw = t$draw,
          members = g$members
        )
      } else {
        NULL
      }
    }),
    lapply(spec$triggers %||% list(), function(t) {
      list(
        id = t$id,
        param = t$param,
        q = t$q,
        draw = t$draw,
        members = t$members
      )
    })
  )
  all_trigger_defs <- Filter(Negate(is.null), all_trigger_defs)

  for (trig in all_trigger_defs) {
    trig_id <- trig$id %||% trig$param %||% NA_character_
    if (is.null(trig_id) || !nzchar(trig_id)) next
    draw_mode <- trig$draw %||% "shared"
    q_val <- NA_real_
    if (!is.null(trig$param) && trig$param %in% names(param_values)) {
      q_val <- as.numeric(param_values[[trig$param]])
    } else if (!is.null(trig$q)) {
      q_val <- as.numeric(trig$q)
    }
    if (!is.finite(q_val) || is.na(q_val)) q_val <- 0

    members <- grp$members %||% character(0)
    if (length(members) == 0L) next
    for (m in members) {
      if (!m %in% acc_ids) next
      existing <- acc_trigger_map[[m]]
      if (!is.null(existing)) {
        stop(sprintf(
          "Accumulator '%s' is assigned to multiple trigger groups ('%s' and '%s')",
          m, existing$id %||% "<unknown>", trig_id
        ))
      }
      acc_trigger_map[[m]] <- list(
        id = trig_id,
        draw = draw_mode,
        q = q_val,
        param = trig$param %||% NULL
      )
    }
  }

  # Build one row per accumulator
  base_mat <- matrix(NA_real_, nrow = length(accs), ncol = length(col_names))
  colnames(base_mat) <- col_names
  colnames(base_mat) <- col_names
  for (i in seq_along(accs)) {
    acc <- accs[[i]]
    acc_id <- acc_ids[[i]]
    dist <- acc$dist
    dist_params <- dist_param_names(dist)
    needed <- length(dist_params)
    if (needed > max_p) stop("Unexpected dist param count > max_p")
    # Distribution parameters
    p_vals <- numeric(max_p)
    for (j in seq_along(dist_params)) {
      nm <- paste0(acc_id, ".", dist_params[[j]])
      if (nm %in% names(param_values)) {
        p_vals[[j]] <- as.numeric(param_values[[nm]])
      } else if (!is.null(acc$params[[dist_params[[j]]]])) {
        p_vals[[j]] <- as.numeric(acc$params[[dist_params[[j]]]])
      } else {
        stop(sprintf("Missing required parameter '%s' for accumulator '%s'", dist_params[[j]], acc_id))
      }
    }
    # q
    q_nm <- paste0(acc_id, ".q")
    has_explicit_q <- q_nm %in% names(param_values)
    q_val <- if (has_explicit_q) {
      as.numeric(param_values[[q_nm]])
    } else if (!is.null(acc$q)) {
      as.numeric(acc$q)
    } else {
      NA_real_
    }
    trig_info <- acc_trigger_map[[acc_id]] %||% NULL
    if (!is.null(trig_info)) {
      if (identical(trig_info$draw, "shared")) {
        if (has_explicit_q && !isTRUE(all.equal(q_val, trig_info$q))) {
          stop(sprintf(
            "Accumulator '%s' provides '%s' but belongs to shared trigger '%s' with q=%.6f",
            acc_id, q_nm, trig_info$id, trig_info$q
          ))
        }
        q_val <- trig_info$q
      } else if (identical(trig_info$draw, "independent")) {
        if (!has_explicit_q) q_val <- trig_info$q
      } else {
        stop(sprintf("Unknown trigger draw mode '%s'", trig_info$draw))
      }
    }
    if (!is.finite(q_val) || is.na(q_val)) q_val <- 0
    # t0
    t0_nm <- paste0(acc_id, ".t0")
    t0_val <- if (t0_nm %in% names(param_values)) {
      as.numeric(param_values[[t0_nm]])
    } else if (!is.null(acc$params$t0)) {
      as.numeric(acc$params$t0)
    } else {
      0
    }
    # w
    w_nm <- paste0(acc_id, ".w")
    w_val <- if (w_nm %in% names(param_values)) {
      as.numeric(param_values[[w_nm]])
    } else {
      comp_id <- comp_attr[[acc_id]] %||% NA_character_
      if (!is.null(comp_id) && !is.na(comp_id) && comp_id %in% names(comp_weights)) {
        as.numeric(comp_weights[[comp_id]])
      } else {
        1
      }
    }
    row_vals <- c(q_val, w_val, t0_val, p_vals)
    for (wp in weight_params) {
      val <- if (wp %in% names(param_values)) as.numeric(param_values[[wp]]) else NA_real_
      row_vals <- c(row_vals, val)
    }
    base_mat[i, ] <- row_vals
  }

  # If a pre-built layout is provided (static mapping stored in context), honor it.
  if (!is.null(layout)) {
    row_trial <- layout$row_trial %||% layout$trial
    row_acc <- layout$row_acc %||% layout$acc
    if (is.null(row_trial) || is.null(row_acc)) {
      stop("layout must include row_trial and row_acc")
    }
    if (length(row_trial) != length(row_acc)) {
      stop("layout row_trial/row_acc lengths must match")
    }
    rows <- lapply(seq_along(row_trial), function(idx) {
      t <- as.integer(row_trial[[idx]])
      a <- as.integer(row_acc[[idx]])
      if (is.na(t) || is.na(a) || a < 1L || a > length(accs)) {
        stop("layout indices out of range")
      }
      base_mat[a, , drop = FALSE]
    })
    params_mat <- do.call(rbind, rows)
    colnames(params_mat) <- col_names
    return(params_mat)
  }

  # If component labels are provided per trial, build only the rows for accumulators
  # that participate in that component; otherwise, fall back to the full rectangular
  # layout (all accumulators per trial).
  if (!is.null(component)) {
    if (length(component) != n_trials) {
      stop("component vector must have length n_trials when provided")
    }
    rows <- list()
    for (t in seq_len(n_trials)) {
      comp_lbl <- as.character(component[[t]])
      for (i in seq_along(accs)) {
        acc_comp <- accs[[i]]$components %||% character(0)
        if (length(acc_comp) > 0 && !comp_lbl %in% acc_comp) {
          next
        }
        rows[[length(rows) + 1L]] <- base_mat[i, , drop = FALSE]
      }
    }
    if (length(rows) == 0L) {
      stop("No accumulator rows matched the provided component labels")
    }
    params_mat <- do.call(rbind, rows)
    colnames(params_mat) <- col_names
    params_mat
  } else {
    # Replicate for n_trials (pure positional, all accumulators)
    params_mat <- base_mat[rep(seq_len(nrow(base_mat)), times = n_trials), , drop = FALSE]
    colnames(params_mat) <- col_names
    params_mat
  }
}

#' Nest data by accumulator
#'
#' Expand a trial-level data frame so each trial is repeated for every accumulator
#' in a finalized model, adding `trial` (if missing) and `accumulator` columns.
#'
#' @param model_spec Finalized model (from `finalize_model()`)
#' @param data Data frame with one row per trial
#' @return Data frame with rows repeated per accumulator and an `accumulator` column
#' @export
nest_accumulators <- function(model_spec, data) {
  structure <- .as_model_structure(model_spec)
  acc_defs <- structure$prep$accumulators %||% list()
  acc_ids <- names(acc_defs)
  if (length(acc_ids) == 0L) {
    stop("Model must define accumulators")
  }
  df <- as.data.frame(data)
  if (!"trial" %in% names(df)) {
    df$trial <- seq_len(nrow(df))
  }

  use_component <- "component" %in% names(df)
  comp_ids <- structure$components$component_id %||% character(0)
  acc_by_comp <- NULL
  if (use_component) {
    comp_col <- df$component
    comp_levels <- if (is.factor(comp_col)) levels(comp_col) else NULL
    data_comps <- unique(as.character(stats::na.omit(comp_col)))
    if (length(data_comps) > 0L && !all(data_comps %in% comp_ids)) {
      stop("component column values must match model components")
    }
    acc_by_comp <- lapply(comp_ids, function(cmp) {
      keep <- vapply(acc_defs, function(a) {
        comps <- a$components %||% character(0)
        length(comps) == 0L || cmp %in% comps
      }, logical(1))
      acc_ids[keep]
    })
    names(acc_by_comp) <- comp_ids
  }

  n_trials <- nrow(df)
  accs_per_trial <- vector("list", n_trials)
  if (use_component) {
    for (t in seq_len(n_trials)) {
      comp_val <- comp_col[t]
      if (is.na(comp_val) || !nzchar(as.character(comp_val))) {
        accs_per_trial[[t]] <- acc_ids
      } else {
        accs_per_trial[[t]] <- acc_by_comp[[as.character(comp_val)]] %||% acc_ids
      }
    }
  } else {
    for (t in seq_len(n_trials)) accs_per_trial[[t]] <- acc_ids
  }

  total_rows <- sum(vapply(accs_per_trial, length, integer(1)))
  out <- vector("list", length(df) + 1L)
  names(out) <- c(names(df), "accumulator")
  idx <- 1L
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.factor(col)) {
      out[[nm]] <- character(total_rows)
    } else if (is.logical(col)) {
      out[[nm]] <- logical(total_rows)
    } else if (is.integer(col)) {
      out[[nm]] <- integer(total_rows)
    } else if (is.numeric(col)) {
      out[[nm]] <- numeric(total_rows)
    } else {
      out[[nm]] <- vector(mode(col), total_rows)
    }
  }
  out$accumulator <- character(total_rows)

  for (t in seq_len(n_trials)) {
    accs_t <- accs_per_trial[[t]]
    len <- length(accs_t)
    if (len == 0L) next
    rng <- idx:(idx + len - 1L)
    for (nm in names(df)) {
      col <- df[[nm]]
      val <- col[t]
      if (is.factor(col)) {
        out[[nm]][rng] <- as.character(val)
      } else {
        out[[nm]][rng] <- val
      }
    }
    out$accumulator[rng] <- accs_t
    idx <- idx + len
  }
  out_df <- as.data.frame(out, stringsAsFactors = FALSE)
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.factor(col)) {
      out_df[[nm]] <- factor(out_df[[nm]], levels = levels(col), ordered = is.ordered(col))
    }
  }
  if (use_component && is.factor(df$component)) {
    out_df$component <- factor(out_df$component, levels = levels(df$component), ordered = is.ordered(df$component))
  }
  out_df
}
