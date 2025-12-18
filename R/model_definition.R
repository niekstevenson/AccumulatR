`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

# ------------------------------------------------------------------------------
# Expression parsing utilities
# ------------------------------------------------------------------------------

.is_expr_node <- function(x) is.list(x) && length(x) > 0 && !is.null(x$kind)

.expr_from_symbol <- function(sym) {
  list(kind = "event", source = as.character(sym), k = NULL)
}

.expr_label <- function(node) {
  if (is.symbol(node)) return(as.character(node))
  if (is.character(node) && length(node) == 1L) return(node)
  stop("Unable to convert expression label of type '", typeof(node),"'")
}

.expr_from_value <- function(val) {
  if (.is_expr_node(val)) return(val)
  if (inherits(val, "formula")) return(.build_expr(val[[length(val)]]))
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
  if (is.symbol(val)) return(.expr_from_symbol(val))
  if (is.call(val)) return(.parse_expr_call(val))
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

.coerce_unless_list <- function(x) {
  if (is.null(x)) return(list())
  if (.is_expr_node(x)) return(list(x))
  if (is.symbol(x)) return(list(.expr_from_symbol(x)))
  if (inherits(x, "formula")) return(list(.expr_from_value(x)))
  if (is.atomic(x)) return(lapply(as.list(x), .expr_from_value))
  if (is.list(x) && length(x) > 0L) return(lapply(x, .expr_from_value))
  if (is.call(x)) return(list(.expr_from_value(x)))
  stop("Cannot interpret guard 'unless' expression")
}

.parse_guard_call <- function(call) {
  args <- as.list(call)[-1]
  nm <- names(call)[-1]

  blocker <- reference <- unless <- NULL
  if (length(args) >= 1L && (is.null(nm[[1L]]) || nm[[1L]] == "")) blocker <- args[[1L]]
  if (length(args) >= 2L && (is.null(nm[[2L]]) || nm[[2L]] == "")) reference <- args[[2L]]
  if ("blocker" %in% nm) blocker <- args[[which(nm == "blocker")[1L]]]
  if ("reference" %in% nm) reference <- args[[which(nm == "reference")[1L]]]
  if ("unless" %in% nm) unless <- args[[which(nm == "unless")[1L]]]
  if (is.null(blocker) || is.null(reference)) {
    stop("guard() requires both blocker and reference arguments")
  }
  list(
    kind = "guard",
    blocker = .expr_from_value(blocker),
    reference = .expr_from_value(reference),
    unless = .coerce_unless_list(unless)
  )
}

.parse_expr_call <- function(call) {
  op <- as.character(call[[1]])
  if (op == "(") return(.expr_from_value(call[[2]]))
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
  if (identical(op, "event")) return(.parse_event_call(call))
  if (identical(op, "guard")) return(.parse_guard_call(call))
  stop(sprintf("Unsupported token '%s' in expression", op))
}

.build_expr <- function(expr) {
  if (.is_expr_node(expr)) return(expr)
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
#' @param unless Optional list of unless expressions
#' @return Guard expression list
#' @examples
#' inhibit("A", "B")
#' @export
inhibit <- function(reference, by, unless = NULL) {
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
  unless_list <- if (is.null(unless)) list() else .coerce_unless_list(unless)
  list(kind = "guard", blocker = blocker_expr, reference = ref_expr, unless = unless_list)
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

# ------------------------------------------------------------------------------
# Model builder
# ------------------------------------------------------------------------------

#' Create an empty race specification
#'
#' @return A race_spec object
#' @examples
#' race_spec()
#' @export
race_spec <- function() {
  structure(list(
    accumulators = list(),
    pools = list(),
    outcomes = list(),
    groups = list(),
    metadata = list()
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
#' @param onset Onset shift
#' @param q Guess probability (0-1)
#' @param tags Optional tags list
#' @param params Optional parameter list
#' @param ... Additional named parameters
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' @export
add_accumulator <- function(spec, id, dist, onset = 0, q = 0,
                            tags = list(), params = NULL, ...) {
  stopifnot(inherits(spec, "race_spec"))
  param_list <- .collect_params(params, list(...))
  param_list <- .ensure_acc_param_t0(param_list)
  spec$accumulators[[length(spec$accumulators) + 1L]] <- list(
    id = id,
    dist = dist,
    params = param_list,
    onset = onset,
    q = q,
    tags = tags %||% list()
  )
  spec
}

#' Add a pool definition
#'
#' @param spec race_spec object
#' @param id Pool id
#' @param members Member ids
#' @param k Threshold parameter
#' @param weights Optional weights
#' @param tags Optional tags
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_pool(spec, "P1", members = c("A", "B"), k = 1L)
#' @export
add_pool <- function(spec, id, members, k = 1L, weights = NULL, tags = list()) {
  stopifnot(inherits(spec, "race_spec"))
  if (missing(members) || length(members) == 0) {
    stop("Pool must define at least one member")
  }
  spec$pools[[length(spec$pools) + 1L]] <- list(
    id = id,
    members = as.character(members),
    rule = list(kind = "k_of_n", k = k, weights = weights),
    tags = tags %||% list()
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
#' @param members Member ids
#' @param attrs Attributes list
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_group(spec, "G1", members = c("A", "B"),
#'   attrs = list(shared_params = list(q = 0.2)))
#' @export
add_group <- function(spec, id, members, attrs = list()) {
  stopifnot(inherits(spec, "race_spec"))
  if (missing(members) || is.null(members) || length(members) == 0) {
    stop("Group must specify members")
  }
  spec$groups[[length(spec$groups) + 1L]] <- list(
    id = id,
    members = as.character(members),
    attrs = attrs %||% list()
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

# Trigger/helpers for clearer API ---------------------------------------------------

#' Define a joint trigger gate shared across members
#' 
#' @param id Trigger id (optional; defaults to group id)
#' @param q Probability (0-1) for trigger failure
#' @param param Optional parameter name supplying q
#' @param draw Either "shared" (one draw for all members) or "independent" (same q, independent draws)
#' @return Trigger descriptor list
#' @examples
#' trigger(q = 0.1, param = "gate_q", draw = "shared")
#' @export
trigger <- function(id = NULL, q = NULL, param = NULL, draw = c("shared", "independent")) {
  draw <- match.arg(draw)
  list(shared_trigger = list(id = id, q = q, param = param, draw = draw))
}

#' Build an expression guard node
#'
#' @param blocker Blocker expression or id
#' @param reference Reference expression or id
#' @param unless Optional list of unless expressions
#' @return Guard expression list
#' @examples
#' expr_guard("B", "A")
#' @export
expr_guard <- function(blocker, reference, unless = list()) {
  list(kind = "guard", blocker = blocker, reference = reference, unless = unless)
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

#' Define a component entry
#'
#' @param id Component id
#' @param weight Optional weight
#' @param attrs Optional attributes list
#' @return Component list
#' @examples
#' component("fast", weight = 0.7)
#' @export
component <- function(id, weight = NULL, attrs = list()) {
  list(id = id, weight = weight, attrs = attrs %||% list())
}

race_model <- function(accumulators, pools = list(), outcomes, groups = list(), metadata = list()) {
  if (inherits(accumulators, "race_spec")) {
    if (!missing(pools) || !missing(outcomes) || !missing(groups) || !missing(metadata)) {
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
  if (inherits(accumulators, "race_model_spec")) return(accumulators)
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
    groups = lapply(groups %||% list(), clone_obj),
    metadata = metadata %||% list()
  ), class = "race_model_spec")
}

# ----------------------------------------------------------------------
# Model normalization/finalization (shared by simulation and likelihood)
# ----------------------------------------------------------------------

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
  if (!is.null(model$groups) && length(model$groups) > 0) {
    for (grp in model$groups) {
      members <- grp$members
      attrs <- grp$attrs
      if (is.null(members) || length(members) == 0) next
      trig <- attrs$shared_trigger %||% NULL
      if (!is.null(trig)) {
        draw_mode <- trig$draw %||% "shared"
        trig_id <- trig$id %||% grp$id
        if (is.null(trig_id) || trig_id == "") {
          stop(sprintf("Group '%s' trigger must provide a non-empty id", grp$id %||% "<unnamed>"))
        }
        default_q <- trig$q %||% NA_real_
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
              "Group '%s' trigger requires a single q value; found %s",
              trig_id,
              paste(format(unique(q_vals)), collapse = ", ")
            ))
          }
        }
        param_name <- trig$param %||% NULL
        if (identical(draw_mode, "shared")) {
          shared_triggers[[trig_id]] <- list(
            id = trig_id,
            group_id = grp$id %||% trig_id,
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
        } else {
          stop(sprintf("Unknown trigger draw mode '%s'", draw_mode))
        }
      }
      if (!is.null(attrs$shared_params)) {
        for (m in members) {
          if (!is.null(defs[[m]])) {
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
    }
  }

  for (acc_id in names(defs)) {
    defs[[acc_id]]$params <- .ensure_acc_param_t0(defs[[acc_id]]$params)
  }
  list(acc = defs, shared_triggers = shared_triggers)
}

.extract_components <- function(model) {
  mixture <- model$metadata$mixture %||% list()
  mode <- mixture$mode %||% "fixed"
  comps <- mixture$components %||% list()
  if (length(comps) == 0) {
    return(list(
      ids = "__default__",
      weights = 1,
      attrs = list(`__default__` = list(guess = NULL)),
      has_weight_param = FALSE,
      mode = mode,
      reference = "__default__"
    ))
  }
  ids <- vapply(comps, `[[`, character(1), "id")
  attrs <- setNames(vector("list", length(ids)), ids)
  has_wparam <- logical(length(ids))
  weights <- vapply(comps, function(cmp) cmp$weight %||% 1, numeric(1))
  weight_params <- vapply(comps, function(cmp) cmp$attrs$weight_param %||% NA_character_, character(1))

  reference <- mixture$reference %||% NA_character_
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
      wp <- comps[[i]]$attrs$weight_param %||% NULL
      if (!identical(cmp_id, reference) && (is.null(wp) || !nzchar(wp))) {
        stop(sprintf("Component '%s' must define weight_param in sampled mixtures (reference is '%s')", cmp_id, reference))
      }
      has_wparam[[i]] <- !is.null(wp) && nzchar(wp)
      attrs[[cmp_id]] <- list(
        guess = comps[[i]]$attrs$guess %||% NULL,
        weight_param = wp
      )
    }
  } else {
    # fixed mode: weights only used for defaults
    if (all(is.na(weights))) weights <- rep(1, length(weights))
    weights <- weights / sum(weights)
    for (i in seq_along(ids)) {
      wp <- comps[[i]]$attrs$weight_param %||% NULL
      has_wparam[[i]] <- !is.null(wp) && nzchar(wp)
      attrs[[ids[[i]]]] <- list(
        guess = comps[[i]]$attrs$guess %||% NULL,
        weight_param = wp
      )
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

  model_norm <- unwrap_model_spec(model)
  model_norm <- .normalize_model(model_norm)
  model_norm <- unwrap_model_spec(model_norm)
  if (!inherits(model_norm, "race_model_spec")) {
    stop("finalize_model requires a race_model_spec")
  }

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

# ------------------------------------------------------------------------------
# Parameter utilities
# ------------------------------------------------------------------------------

dist_param_names <- function(dist) {
  dist <- tolower(dist)
  if (dist == "lognormal") return(c("meanlog", "sdlog"))
  if (dist == "gamma") return(c("shape", "rate"))
  if (dist == "exgauss") return(c("mu", "sigma", "tau"))
  character(0)
}

#' List sampled parameter names for a model
#'
#' @param model race_spec or race_model_spec
#' @return Character vector of parameter names
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
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
#' @return Data frame of parameters per trial
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' spec <- add_outcome(spec, "A_win", "A")
#' vals <- c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0)
#' build_params_df(spec, vals, n_trials = 2)
#' @export
build_params_df <- function(model,
                            param_values,
                            n_trials = 1L,
                            component = NULL) {
  build_param_matrix(model, param_values, n_trials = n_trials)
}

#' Build a positional parameter matrix (q, w, t0, p1..pK)
#'
#' @param model Model definition (race_spec or race_model_spec)
#' @param param_values Named numeric vector of core parameters
#' @param n_trials Number of trials to replicate
#' @return Numeric matrix with columns `q`, `w`, `t0`, `p1`, `p2`, `p3`
#'   and attributes `acc_ids`, `n_trials`, `layout`
#' @export
build_param_matrix <- function(model,
                               param_values,
                               n_trials = 1L) {
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

  # Shared parameter groups (R-side only): expand group-level values to member accumulators.
  shared_param_values <- list()
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
  if (length(shared_param_values) > 0L) {
    expanded <- unlist(shared_param_values, use.names = TRUE)
    param_values <- c(param_values, expanded[setdiff(names(expanded), names(param_values))])
  }

  # Map accumulators to components (if any) via group attrs$component
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

  # Mixture component weights (mapped to components, then to member rows)
  mix <- spec$metadata$mixture %||% list()
  comp_defs <- mix$components %||% list()
  comp_ids <- vapply(comp_defs, `[[`, character(1), "id")
  comp_weights <- setNames(rep(NA_real_, length(comp_defs)), comp_ids)
  comp_mode <- mix$mode %||% "fixed"
  comp_ref <- mix$reference %||% if (length(comp_ids) > 0) comp_ids[[1]] else NA_character_
  if (length(comp_defs) > 0) {
    for (i in seq_along(comp_defs)) {
      comp <- comp_defs[[i]]
      attrs <- comp$attrs %||% list()
      wp <- attrs$weight_param %||% NULL
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
  for (grp in spec$groups %||% list()) {
    trig <- grp$attrs$shared_trigger %||% NULL
    if (is.null(trig)) next
    trig_id <- trig$id %||% grp$id %||% trig$param %||% NA_character_
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
    base_mat[i, ] <- c(q_val, w_val, t0_val, p_vals)
  }

  # Replicate for n_trials (pure positional)
  params_mat <- base_mat[rep(seq_len(nrow(base_mat)), times = n_trials), , drop = FALSE]
  layout <- data.frame(
    row = seq_len(nrow(params_mat)),
    trial = rep(seq_len(n_trials), each = nrow(base_mat)),
    accumulator = rep(acc_ids, times = n_trials),
    stringsAsFactors = FALSE
  )
  attr(params_mat, "acc_ids") <- acc_ids
  attr(params_mat, "n_trials") <- n_trials
  attr(params_mat, "layout") <- layout
  class(params_mat) <- c("param_matrix", class(params_mat))
  params_mat
}
