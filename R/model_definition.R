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

build_outcome_expr <- function(expr) {
  .build_expr(expr)
}

# ------------------------------------------------------------------------------
# Public DSL helpers
# ------------------------------------------------------------------------------

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

first_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("first_of() requires at least one argument")
  if (length(args) == 1 && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }
  list(kind = "or", args = lapply(args, .expr_from_value))
}

all_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("all_of() requires at least one argument")
  if (length(args) == 1 && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }
  list(kind = "and", args = lapply(args, .expr_from_value))
}

none_of <- function(expr) {
  inner <- if (is.character(expr) && length(expr) == 1L) {
    list(kind = "event", source = expr, k = NULL)
  } else {
    .build_expr(expr)
  }
  list(kind = "not", arg = inner)
}

exclude <- none_of

# ------------------------------------------------------------------------------
# Model builder
# ------------------------------------------------------------------------------

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

add_outcome <- function(spec, label, expr, options = list()) {
  stopifnot(inherits(spec, "race_spec"))
  spec$outcomes[[length(spec$outcomes) + 1L]] <- list(
    label = label,
    expr = build_outcome_expr(expr),
    options = options %||% list()
  )
  spec
}

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

set_metadata <- function(spec, ...) {
  stopifnot(inherits(spec, "race_spec"))
  updates <- list(...)
  for (nm in names(updates)) {
    spec$metadata[[nm]] <- updates[[nm]]
  }
  spec
}

build_model <- function(spec) {
  stopifnot(inherits(spec, "race_spec"))
  structure(list(
    accumulators = unname(spec$accumulators),
    pools = unname(spec$pools),
    outcomes = unname(spec$outcomes),
    groups = unname(spec$groups),
    metadata = spec$metadata
  ), class = "race_model_spec")
}

# ------------------------------------------------------------------------------
# Low-level constructors
# ------------------------------------------------------------------------------

acc <- function(id, dist, params, onset = 0, q = 0, tags = list()) {
  params <- .ensure_acc_param_t0(params %||% list())
  list(id = id, dist = dist, params = params, onset = onset, q = q, tags = tags)
}

rule_k_of_n <- function(k = 1L, weights = NULL) {
  list(kind = "k_of_n", k = k, weights = weights)
}

pool <- function(id, members, rule = rule_k_of_n(1L), tags = list()) {
  list(id = id, members = members, rule = rule, tags = tags)
}

expr_event <- function(source_id, k = NULL) {
  list(kind = "event", source = source_id, k = k)
}

expr_and <- function(...) list(kind = "and", args = list(...))
expr_or  <- function(...) list(kind = "or",  args = list(...))
expr_not <- function(expr) list(kind = "not", arg = expr)

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

group <- function(id, members, attrs = list()) {
  list(id = id, members = members, attrs = attrs %||% list())
}

component <- function(id, weight = NULL, attrs = list()) {
  list(id = id, weight = weight, attrs = attrs %||% list())
}

race_model <- function(accumulators, pools = list(), outcomes, groups = list(), metadata = list()) {
  if (inherits(accumulators, "race_spec")) {
    if (!missing(pools) || !missing(outcomes) || !missing(groups) || !missing(metadata)) {
      stop("When passing a race_spec, provide it alone and call race_model(spec)")
    }
    return(build_model(accumulators))
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

sampled_pars <- function(model) {
  spec <- race_model(model)
  params <- character(0)
  for (acc in spec$accumulators) {
    base <- dist_param_names(acc$dist)
    acc_params <- c(base, "onset", "q", "t0")
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

populate_model_defaults <- function(model, param_values) {
  spec <- race_model(model)
  accs <- spec[["accumulators"]] %||% list()
  groups <- spec[["groups"]] %||% list()
  if (missing(param_values) || length(param_values) == 0) return(spec)
  if (is.null(names(param_values)) || any(!nzchar(names(param_values)))) {
    stop("populate_model_defaults() requires named parameter values")
  }

  acc_ids <- vapply(accs, `[[`, character(1), "id")
  acc_lookup <- setNames(seq_along(acc_ids), acc_ids)
  group_ids <- vapply(groups, `[[`, character(1), "id")
  group_lookup <- setNames(seq_along(group_ids), group_ids)

  split_name <- function(x) {
    idx <- regexpr("\\.[^.]+$", x)
    if (idx <= 0) return(NULL)
    list(
      target = substr(x, 1, idx - 1),
      param = substr(x, idx + 1, nchar(x))
    )
  }

  coerce_shared_params <- function(shared) {
    if (is.null(shared)) return(list())
    if (!is.list(shared)) shared <- as.list(shared)
    nm <- names(shared)
    out <- list()
    for (i in seq_along(shared)) {
      key <- nm[[i]]
      val <- shared[[i]]
      if (!is.null(key) && nzchar(key)) {
        out[[key]] <- val
      } else if (is.character(val) && length(val) == 1L && nzchar(val)) {
        out[[val]] <- out[[val]] %||% NULL
      }
    }
    out
  }

  for (nm in names(param_values)) {
    pieces <- split_name(nm)
    if (is.null(pieces)) next
    target <- pieces$target
    param <- pieces$param
    value <- param_values[[nm]]
    if (target %in% names(acc_lookup)) {
      idx <- acc_lookup[[target]]
      acc <- accs[[idx]]
      if (identical(param, "onset")) {
        acc$onset <- value
      } else if (identical(param, "q")) {
        acc$q <- value
      } else {
        acc$params <- acc$params %||% list()
        acc$params[[param]] <- value
      }
      accs[[idx]] <- acc
    } else if (target %in% names(group_lookup)) {
      idx <- group_lookup[[target]]
      grp <- groups[[idx]]
      grp_attrs <- grp$attrs %||% list()
      shared <- coerce_shared_params(grp_attrs$shared_params)
      shared[[param]] <- value
      grp_attrs$shared_params <- shared
      grp$attrs <- grp_attrs
      groups[[idx]] <- grp
    }
  }
  spec[["accumulators"]] <- accs
  spec[["groups"]] <- groups
  spec
}
