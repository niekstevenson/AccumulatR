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
#' @param q Guess probability (logit allowed)
#' @param tags Optional tags list
#' @param params Optional parameter list
#' @param ... Additional named parameters
#' @return Updated race_spec
#' @examples
#' spec <- race_spec()
#' spec <- add_accumulator(spec, "A", "lognormal",
#'   params = list(meanlog = 0, sdlog = 0.1))
#' @export
add_accumulator <- function(spec, id, dist, onset = 0, q = -Inf,
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
#' @param id Trigger id
#' @param q Optional probability value
#' @param param Optional parameter name supplying q
#' @return Trigger descriptor list
#' @examples
#' joint_trigger("gate1", q = 0.2)
#' @export
joint_trigger <- function(id, q = NULL, param = NULL) {
  if (missing(id) || is.null(id) || !nzchar(id)) {
    stop("joint_trigger requires a non-empty id")
  }
  list(shared_trigger = list(id = id, q = q, param = param))
}

#' Shared trigger probability placeholder
#'
#' @return Placeholder used in model specs
#' @examples
#' shared_q()
#' @export
shared_q <- function() {
  list(shared_params = list("q"))
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
                            component = NA_character_) {
  spec <- race_model(model)
  accs <- spec[["accumulators"]] %||% list()
  if (length(accs) == 0) {
    stop("Model must define accumulators before building parameter tables")
  }
  prep_defaults <- .prepare_model_for_likelihood(spec)
  prep_accs <- prep_defaults$accumulators %||% list()
  if (is.null(param_values) || length(param_values) == 0) {
    stop("Named parameter values are required")
  }
  if (is.null(names(param_values)) || any(!nzchar(names(param_values)))) {
    stop("Parameter values must be a named vector")
  }
  if (!is.numeric(n_trials) || length(n_trials) != 1L || n_trials < 1) {
    stop("n_trials must be a positive integer")
  }
  n_trials <- as.integer(n_trials)

  value_names <- names(param_values)
  param_values <- as.numeric(param_values)
  names(param_values) <- value_names

  acc_ids <- vapply(accs, `[[`, character(1), "id")
  acc_dists <- vapply(accs, `[[`, character(1), "dist")

  per_acc_params <- setNames(
    lapply(acc_dists, function(dist) c(dist_param_names(dist), "onset", "q", "t0")),
    acc_ids
  )
  all_required_names <- unlist(mapply(function(acc, fields) {
    paste0(acc, ".", fields)
  }, acc_ids, per_acc_params, SIMPLIFY = FALSE))
  optional_suffixes <- c("onset", "q", "t0")
  suffix_union <- unique(unlist(per_acc_params))
  shared_map <- setNames(vector("list", length(acc_ids)), acc_ids)
  spec_groups <- spec[["groups"]] %||% list()
  if (length(spec_groups) > 0) {
    for (grp in spec_groups) {
      grp_attrs <- grp[["attrs"]] %||% list()
      sp <- grp_attrs[["shared_params"]] %||% NULL
      members <- grp[["members"]] %||% character(0)
      if (is.null(sp) || length(members) == 0) next
      fields <- names(sp)
      if (is.null(fields) || any(!nzchar(fields))) {
        fields <- as.character(unlist(sp, use.names = FALSE))
      }
      grp_id <- grp[["id"]] %||% NA_character_
      for (member in members) {
        if (is.null(shared_map[[member]])) shared_map[[member]] <- list()
        for (field in fields) {
          shared_map[[member]][[field]] <- paste0(grp_id, ".", field)
        }
      }
    }
  }
  missing <- character(0)
  for (acc_id in acc_ids) {
    required_suffixes <- setdiff(per_acc_params[[acc_id]], optional_suffixes)
    for (suffix in required_suffixes) {
      nm <- paste0(acc_id, ".", suffix)
      shared_nm <- shared_map[[acc_id]][[suffix]] %||% NA_character_
      if (nm %in% names(param_values)) next
      if (!is.na(shared_nm) && shared_nm %in% names(param_values)) next
      missing <- c(missing, nm)
    }
  }
  if (length(missing) > 0L) {
    stop("Missing parameter values for: ", paste(missing, collapse = ", "))
  }
  optional_default <- function(spec_acc, prep_acc, suffix) {
    acc <- prep_acc %||% spec_acc %||% list()
    if (identical(suffix, "onset")) {
      return(as.numeric(acc$onset %||% 0))
    }
    if (identical(suffix, "q")) {
      val <- acc$q %||% acc$shared_trigger_q %||% spec_acc$q %||% 0
      return(as.numeric(val %||% 0))
    }
    if (identical(suffix, "t0")) {
      params <- acc$params %||% spec_acc$params %||% list()
      return(as.numeric(params[[suffix]] %||% 0))
    }
    0
  }

  base_rows <- lapply(seq_along(accs), function(idx) {
    acc <- accs[[idx]]
    acc_id <- acc_ids[[idx]]
    prep_acc <- prep_accs[[acc_id]] %||% list()
    role_value <- (acc[["tags"]] %||% list())[["role"]] %||% "std"
    row_vals <- setNames(as.list(rep(NA_real_, length(suffix_union))), suffix_union)
    for (suffix in per_acc_params[[acc_id]]) {
      nm <- paste0(acc_id, ".", suffix)
      if (nm %in% names(param_values)) {
        row_vals[[suffix]] <- as.numeric(param_values[[nm]])
      } else if (!is.null(shared_map[[acc_id]][[suffix]])) {
        shared_nm <- shared_map[[acc_id]][[suffix]]
        if (!shared_nm %in% names(param_values)) {
          stop(sprintf("Missing shared parameter value '%s' for member '%s'", shared_nm, acc_id))
        }
        row_vals[[suffix]] <- as.numeric(param_values[[shared_nm]])
      } else if (suffix %in% optional_suffixes) {
        row_vals[[suffix]] <- optional_default(acc, prep_acc, suffix)
      } else {
        stop(sprintf("Missing parameter value for '%s.%s'", acc_id, suffix))
      }
    }
    data.frame(
      trial = 1L,
      accumulator_id = acc_id,
      accumulator = idx,
      accumulator_index = idx,
      component = if (is.null(component)) NA_character_ else as.character(component),
      role = role_value,
      row_vals,
      stringsAsFactors = FALSE
    )
  })
  base_df <- do.call(rbind, base_rows)
  extra_names <- setdiff(names(param_values), all_required_names)
  if (length(extra_names) > 0L) {
    for (nm in extra_names) {
      base_df[[nm]] <- as.numeric(param_values[[nm]])
    }
  }
  # Propagate shared trigger parameters (group-level) into per-acc rows:
  # - set per-acc q to 0 (success-path) and attach shared_trigger_id so the
  #   native layer applies the joint gate once using shared_trigger_q.
  # - if the param value is on the logit scale, store raw value; native will
  #   convert as needed.
  if (!is.null(spec$groups) && length(spec$groups) > 0) {
    if (!"shared_trigger_id" %in% names(base_df)) {
      base_df$shared_trigger_id <- NA_character_
    }
    if (!"shared_trigger_q" %in% names(base_df)) {
      base_df$shared_trigger_q <- NA_real_
    }
    for (grp in spec$groups) {
      trig <- grp$attrs$shared_trigger %||% NULL
      if (is.null(trig)) next
      param_nm <- trig$param %||% NULL
      members <- grp$members %||% character(0)
      if (length(members) == 0) next
      mask <- base_df$accumulator_id %in% members
      if (!any(mask)) next
      base_df$shared_trigger_id[mask] <- grp$id %||% trig$id %||% param_nm
      if (!is.null(param_nm) && nzchar(param_nm) && param_nm %in% names(param_values)) {
        trig_q <- as.numeric(param_values[[param_nm]])
        base_df$shared_trigger_q[mask] <- trig_q
      } else if (!is.null(trig$q)) {
        base_df$shared_trigger_q[mask] <- as.numeric(trig$q)
      }
      # Joint gate: per-acc q is handled by the gate; avoid per-acc q overrides.
      # Leave q as NA for these rows so downstream logic does not treat it as an override.
      base_df$q[mask] <- NA_real_
      # Keep q on logit scale for inspection; add passthrough column
      if (!"shared_trigger_q_logit" %in% names(base_df)) {
        base_df$shared_trigger_q_logit <- NA_real_
      }
      base_df$shared_trigger_q_logit[mask] <- base_df$shared_trigger_q[mask]
    }
  }
  base_n <- nrow(base_df)
  if (n_trials > 1L) {
    base_df <- base_df[rep(seq_len(base_n), times = n_trials), , drop = FALSE]
  }
  base_df$trial <- rep(seq_len(n_trials), each = base_n)
  rownames(base_df) <- NULL
  base_df <- .params_hash_attach(base_df)
  base_df
}
