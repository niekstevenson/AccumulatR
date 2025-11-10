# Simplified declarative API for building race-model specifications
# -----------------------------------------------------------------
# The DSL below keeps things explicit: accumulators feed pools, pools are
# combined through boolean expressions, and outcomes attach labels to those
# expressions. Everything is ordinary R lists so the simulator can walk the
# structure without surprises.

# --- Expression parsing --------------------------------------------------------

# Internal helper to recognise already-built expression nodes
.is_expr_node <- function(x) is.list(x) && length(x) > 0 && !is.null(x$kind)

.expr_from_symbol <- function(sym) {
  list(kind = "event", source = as.character(sym), k = NULL)
}

.expr_from_value <- function(val) {
  if (.is_expr_node(val)) return(val)
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

.coerce_unless_list <- function(x) {
  if (is.null(x)) return(list())
  if (.is_expr_node(x)) return(list(x))
  if (is.symbol(x)) return(list(.expr_from_symbol(x)))
  if (inherits(x, "formula")) return(list(.expr_from_value(x)))
  if (is.atomic(x)) {
    return(lapply(as.list(x), .expr_from_value))
  }
  if (is.list(x) && length(x) > 0L) {
    return(lapply(x, .expr_from_value))
  }
  if (is.call(x)) {
    return(list(.expr_from_value(x)))
  }
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

.expr_label <- function(node) {
  if (is.symbol(node)) return(as.character(node))
  if (is.character(node) && length(node) == 1L) return(node)
  stop("Unable to convert expression label of type '", typeof(node), "'")
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
  if (identical(op, "event")) {
    return(.parse_event_call(call))
  }
  if (identical(op, "guard")) {
    return(.parse_guard_call(call))
  }
  stop(sprintf("Unsupported token '%s' in expression", op))
}

.build_expr <- function(expr) {
  if (.is_expr_node(expr)) return(expr)
  .expr_from_value(expr)
}

build_outcome_expr <- function(expr) {
  .build_expr(expr)
}

# --- Intuitive expression builders --------------------------------------------

#' Inhibit a reference accumulator/pool by a blocker
#'
#' @param reference The accumulator/pool that should fire (if not blocked)
#' @param by The blocker that prevents reference from firing if it finishes first
#' @param unless Optional protector(s) that disable the blocker if they finish before it
#' @return Expression node representing the guard logic
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

  unless_list <- if (is.null(unless)) {
    list()
  } else {
    .coerce_unless_list(unless)
  }

  list(
    kind = "guard",
    blocker = blocker_expr,
    reference = ref_expr,
    unless = unless_list
  )
}

#' First to finish wins (OR logic)
#'
#' @param ... Accumulators/pools/expressions to race
#' @return Expression node representing OR logic
first_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("first_of() requires at least one argument")

  # Handle case where a single list is passed
  if (length(args) == 1L && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }

  expr_list <- lapply(args, function(arg) {
    if (is.character(arg) && length(arg) == 1L) {
      list(kind = "event", source = arg, k = NULL)
    } else {
      .build_expr(arg)
    }
  })

  if (length(expr_list) == 1L) return(expr_list[[1]])
  list(kind = "or", args = expr_list)
}

#' All must finish (AND logic)
#'
#' @param ... Accumulators/pools/expressions that must all finish
#' @return Expression node representing AND logic
all_of <- function(...) {
  args <- list(...)
  if (length(args) == 0) stop("all_of() requires at least one argument")

  # Handle case where a single list is passed
  if (length(args) == 1L && is.list(args[[1]]) && !.is_expr_node(args[[1]])) {
    args <- args[[1]]
  }

  expr_list <- lapply(args, function(arg) {
    if (is.character(arg) && length(arg) == 1L) {
      list(kind = "event", source = arg, k = NULL)
    } else {
      .build_expr(arg)
    }
  })

  if (length(expr_list) == 1L) return(expr_list[[1]])
  list(kind = "and", args = expr_list)
}

#' None must finish (NOT logic)
#'
#' @param expr Accumulator/pool/expression that must not finish
#' @return Expression node representing NOT logic
none_of <- function(expr) {
  if (is.character(expr) && length(expr) == 1L) {
    inner <- list(kind = "event", source = expr, k = NULL)
  } else {
    inner <- .build_expr(expr)
  }
  list(kind = "not", arg = inner)
}

# Alias for none_of
exclude <- none_of

# --- Builder helpers -----------------------------------------------------------

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

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

# --- Core constructors --------------------------------------------------------

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

# Expressions compose pools/accumulators
expr_event <- function(source_id, k = NULL) {
  list(kind = "event", source = source_id, k = k)
}

expr_and <- function(...) list(kind = "and", args = list(...))
expr_or  <- function(...) list(kind = "or",  args = list(...))
expr_not <- function(expr) list(kind = "not", arg = expr)

# Guard: suppress reference unless protectors fire before blocker
expr_guard <- function(blocker, reference, unless = list()) {
  list(kind = "guard", blocker = blocker, reference = reference, unless = unless)
}

clone_obj <- function(x) {
  if (is.environment(x)) stop("Cannot clone environments for race models")
  if (is.list(x)) {
    out <- lapply(x, clone_obj)
    attrs <- attributes(x)
    if (!is.null(attrs)) {
      attributes(out) <- lapply(attrs, clone_obj)
    }
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
  if (inherits(accumulators, "race_model_spec")) {
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
    groups = lapply(groups %||% list(), clone_obj),
    metadata = metadata %||% list()
  ), class = "race_model_spec")
}

# --- Example library ----------------------------------------------------------

# Example 1 – simple two-response race built with the helper pipeline
example_1_simple <- race_spec() |>
  add_accumulator("go1", "lognormal", meanlog = log(0.30), sdlog = 0.18) |>
  add_accumulator("go2", "lognormal", meanlog = log(0.32), sdlog = 0.18) |>
  add_pool("R1", "go1") |>
  add_pool("R2", "go2") |>
  add_outcome("R1", "R1") |>
  add_outcome("R2", "R2") |>
  build_model()

# Example 2 – stop/go mixture with gating and inhibition
example_2_stop_mixture <- race_spec() |>
  add_accumulator("go1", "lognormal", meanlog = log(0.35), sdlog = 0.2) |>
  add_accumulator("stop", "exgauss", onset = 0.20,
                  params = list(mu = 0.1, sigma = 0.04, tau = 0.1)) |>
  add_accumulator("go2", "lognormal", meanlog = log(0.60), sdlog = 0.18,
                  onset = 0.20) |>
  add_pool("GO1", "go1") |>
  add_pool("STOP", "stop") |>
  add_pool("GO2", "go2") |>
  add_outcome("R1", inhibit("GO1", by = "STOP"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("R2", all_of("GO2", "STOP"), options = list(component = "go_stop")) |>
  add_group("component:go_only", members = c("go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:go_stop", members = c("go1", "stop", "go2"),
            attrs = list(component = "go_stop")) |>
  set_metadata(mixture = list(
    components = list(component("go_only", weight = 0.5),
                       component("go_stop", weight = 0.5))
  )) |>
  build_model()

# Example 3 – stop outcome mapped to NA
example_3_stop_na <- race_spec() |>
  add_accumulator("go_left", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("go_right", "lognormal", params = list(meanlog = log(0.32), sdlog = 0.20)) |>
  add_accumulator("stop", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.15), sdlog = 0.18)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("STOP", "stop") |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("STOP", "STOP", options = list(map_outcome_to = NA_character_)) |>
  build_model()

# Example 4 – two accumulators vs one
example_4_two_on_one <- race_spec() |>
  add_accumulator("R1_A", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.15)) |>
  add_accumulator("R1_B", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("R2", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.18)) |>
  add_pool("R1", c("R1_A", "R1_B")) |>
  add_pool("R2", "R2") |>
  add_outcome("R1", "R1") |>
  add_outcome("R2", "R2") |>
  build_model()

# Example 5 – timeout with weighted guess
example_5_timeout_guess <- race_spec() |>
  add_accumulator("go_left", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.18)) |>
  add_accumulator("go_right", "lognormal", params = list(meanlog = log(0.325), sdlog = 0.18)) |>
  add_accumulator("timeout", "lognormal", onset = 0.05,
                  params = list(meanlog = log(0.25), sdlog = 0.10)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("TO", "timeout") |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("TIMEOUT", "TO", options = list(
    guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
  )) |>
  build_model()

# Example 6 – dual path (A & C) OR (B & C)
example_6_dual_path <- race_spec() |>
  add_accumulator("acc_taskA", "lognormal", params = list(meanlog = log(0.28), sdlog = 0.18)) |>
  add_accumulator("acc_taskB", "lognormal", params = list(meanlog = log(0.32), sdlog = 0.18)) |>
  add_accumulator("acc_gateC", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.18)) |>
  add_pool("TaskA", "acc_taskA") |>
  add_pool("TaskB", "acc_taskB") |>
  add_pool("GateC", "acc_gateC") |>
  add_outcome("Outcome_via_A", all_of("TaskA", "GateC")) |>
  add_outcome("Outcome_via_B", all_of("TaskB", "GateC")) |>
  build_model()

# Example 7 – fast/slow mixture
example_7_mixture <- race_spec() |>
  add_accumulator("target_fast", "lognormal",
                  params = list(meanlog = log(0.25), sdlog = 0.15)) |>
  add_accumulator("target_slow", "lognormal",
                  params = list(meanlog = log(0.45), sdlog = 0.20)) |>
  add_accumulator("competitor", "lognormal",
                  params = list(meanlog = log(0.35), sdlog = 0.18)) |>
  add_pool("TARGET", c("target_fast", "target_slow")) |>
  add_pool("COMP", "competitor") |>
  add_outcome("R1", "TARGET") |>
  add_outcome("R2", "COMP") |>
  add_group("component:fast", members = c("target_fast", "competitor"),
            attrs = list(component = "fast")) |>
  add_group("component:slow", members = c("target_slow", "competitor"),
            attrs = list(component = "slow")) |>
  set_metadata(mixture = list(
    components = list(
      component("fast", weight = 0.2, attrs = list(weight_param = "logit_fast")),
      component("slow", weight = 0.8, attrs = list(weight_param = "logit_slow"))
    )
  )) |>
  build_model()

# Example 8 – shared parameter group
example_8_shared_params <- race_spec() |>
  add_accumulator("go_left", "lognormal", params = list(sdlog = 0.18)) |>
  add_accumulator("go_right", "lognormal", params = list(sdlog = 0.24)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_group("par:shared_mean", members = c("go_left", "go_right"),
            attrs = list(shared_params = list(meanlog = log(0.32)))) |>
  build_model()

# Example 9 – advanced k-of-n logic
example_9_advanced_k <- race_spec() |>
  add_accumulator("a1", "lognormal", params = list()) |>
  add_accumulator("a2", "lognormal", params = list()) |>
  add_accumulator("a3", "lognormal", params = list()) |>
  add_accumulator("b1_1", "lognormal", params = list()) |>
  add_accumulator("b1_2", "lognormal", params = list(meanlog = log(0.50), sdlog = 0.25)) |>
  add_accumulator("b2_1", "lognormal", params = list(meanlog = log(0.48), sdlog = 0.25)) |>
  add_accumulator("b2_2", "lognormal", params = list()) |>
  add_pool("A", c("a1", "a2", "a3"), k = 2L) |>
  add_pool("B1", c("b1_1", "b1_2")) |>
  add_pool("B2", c("b2_1", "b2_2")) |>
  add_outcome("A", "A") |>
  add_outcome("B", all_of("B1", "B2")) |>
  add_group("par:A_shared", members = c("a2", "a3"),
            attrs = list(shared_params = list(meanlog = log(0.575), sdlog = 0.25))) |>
  add_group("par:Cross_shared", members = c("a1", "b1_1", "b2_2"),
            attrs = list(shared_params = list(meanlog = log(0.50), sdlog = 0.25))) |>
  build_model()

# Example 10 – exclusion via guard
example_10_exclusion <- race_spec() |>
  add_accumulator("R1_acc", "lognormal", params = list(meanlog = log(0.35), sdlog = 0.18)) |>
  add_accumulator("R2_acc", "lognormal", params = list(meanlog = log(0.45), sdlog = 0.18)) |>
  add_accumulator("X_acc", "lognormal", params = list(meanlog = log(0.35), sdlog = 0.18)) |>
  add_pool("R1", "R1_acc") |>
  add_pool("R2", "R2_acc") |>
  add_pool("X", "X_acc") |>
  add_outcome("R1", inhibit("R1", by = "X")) |>
  add_outcome("R2", "R2") |>
  build_model()

# Example 11 – censoring pool and deadline
example_11_censor_deadline <- race_spec() |>
  add_accumulator("go_left", "lognormal", params = list(meanlog = log(0.28), sdlog = 0.18)) |>
  add_accumulator("go_right", "lognormal", params = list(meanlog = log(0.34), sdlog = 0.18)) |>
  add_accumulator("censor_watch", "lognormal", params = list(meanlog = log(0.40), sdlog = 0.12)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("CENSOR", "censor_watch") |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("NR_CENSOR", "CENSOR", options = list(class = "censor")) |>
  add_outcome("NR_DEADLINE", "__DEADLINE__", options = list(class = "deadline")) |>
  set_metadata(deadline = 0.55,
               special_outcomes = list(censor = "NR_CENSOR", deadline = "NR_DEADLINE")) |>
  build_model()

# Example 12 – inhibitor with protector
example_12_inhibitor_with_protector <- race_spec() |>
  # Like example 2 (go/stop) but with a safety signal that disables STOP
  add_accumulator("go1", "lognormal", params = list(meanlog = log(0.35), sdlog = 0.2)) |>
  add_accumulator("stop", "exgauss", onset = 0.20,
                  params = list(mu = 0.1, sigma = 0.04, tau = 0.1)) |>
  add_accumulator("go2", "lognormal", meanlog = log(0.60), sdlog = 0.18,
                  onset = 0.20) |>
  add_accumulator("safety", "lognormal", onset = 0.125,
                  params = list(meanlog = log(0.22), sdlog = 0.13)) |>
  add_pool("GO1", "go1") |>
  add_pool("STOP", "stop") |>
  add_pool("GO2", "go2") |>
  add_pool("SAFETY", "safety") |>
  add_outcome("R1", inhibit("GO1", by = "STOP", unless = "SAFETY"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("R2", all_of("GO2", "STOP"), options = list(component = "go_stop")) |>
  add_group("component:go_only", members = c("go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:go_stop", members = c("go1", "stop", "go2", "safety"),
            attrs = list(component = "go_stop")) |>
  set_metadata(mixture = list(
    components = list(component("go_only", weight = 0.5),
                       component("go_stop", weight = 0.5))
  )) |>
  build_model()

# Example 13 – nested pools
example_13_nested_pools <- race_spec() |>
  add_accumulator("a1", "lognormal", params = list(meanlog = log(0.27), sdlog = 0.15)) |>
  add_accumulator("a2", "lognormal", params = list(meanlog = log(0.29), sdlog = 0.15)) |>
  add_accumulator("a3", "lognormal", params = list(meanlog = log(0.33), sdlog = 0.16)) |>
  add_accumulator("b1", "lognormal", params = list(meanlog = log(0.31), sdlog = 0.17)) |>
  add_accumulator("b2", "lognormal", params = list(meanlog = log(0.36), sdlog = 0.19)) |>
  add_pool("A_line", c("a1", "a2")) |>
  add_pool("A_team", c("A_line", "a3"), k = 2L) |>
  add_pool("B_team", c("b1", "b2")) |>
  add_outcome("TeamA", "A_team") |>
  add_outcome("TeamB", "B_team") |>
  build_model()

# Example 14 – weighted pool
example_14_weighted_pool <- race_spec() |>
  add_accumulator("w1", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.18)) |>
  add_accumulator("w2", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.18)) |>
  add_accumulator("w3", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  # Compete two items in WEIGHTED against an external competitor w3
  add_pool("WEIGHTED", c("w1", "w2"), k = 1L, weights = c(0.6, 0.4)) |>
  add_pool("COMP", "w3") |>
  add_outcome("WeightedChoice", "WEIGHTED") |>
  add_outcome("Competitor", "COMP") |>
  build_model()

# Example 15 – mixture with component metadata
example_15_component_metadata <- race_spec() |>
  # Separate shared accumulators per component to enhance identifiability
  add_accumulator("acc_shared_fast", "lognormal",
                  params = list(meanlog = log(0.32), sdlog = 0.10)) |>
  add_accumulator("acc_shared_slow", "lognormal",
                  params = list(meanlog = log(0.60), sdlog = 0.12)) |>
  add_accumulator("acc_fast", "lognormal",
                  params = list(meanlog = log(0.24), sdlog = 0.10)) |>
  add_accumulator("acc_slow", "lognormal",
                  params = list(meanlog = log(0.50), sdlog = 0.15)) |>
  add_pool("FAST", c("acc_fast", "acc_shared_fast")) |>
  add_pool("SLOW", c("acc_slow", "acc_shared_slow")) |>
  add_outcome("Response", first_of("FAST", "SLOW")) |>
  add_outcome("GUESS", "__GUESS__") |>
  add_group("component:fast", members = c("acc_fast", "acc_shared_fast"),
            attrs = list(component = "fast")) |>
  add_group("component:slow", members = c("acc_slow", "acc_shared_slow"),
            attrs = list(component = "slow")) |>
  set_metadata(mixture = list(
    components = list(
      component("fast", weight = 0.6, attrs = list(
        deadline = 0.45,
        guess = list(outcome = "GUESS", weights = c(Response = 0.7))
      )),
      component("slow", weight = 0.4, attrs = list(deadline = 0.75))
    )
  )) |>
  build_model()

# Guard tie – shared gate with stop control
example_16_guard_tie_simple <- race_spec() |>
  add_accumulator("go_fast", "lognormal", params = list(meanlog = log(0.28), sdlog = 0.18)) |>
  add_accumulator("go_slow", "lognormal", params = list(meanlog = log(0.34), sdlog = 0.18)) |>
  add_accumulator("gate_shared", "lognormal", params = list(meanlog = log(0.30), sdlog = 0.16)) |>
  add_accumulator("stop_control", "lognormal", onset = 0.05,
                  params = list(meanlog = log(0.27), sdlog = 0.15)) |>
  add_pool("FAST", "go_fast") |>
  add_pool("SLOW", "go_slow") |>
  add_pool("GATE", "gate_shared") |>
  add_pool("STOP", "stop_control") |>
  add_outcome("Fast", inhibit(all_of("FAST", "GATE"), by = "STOP")) |>
  add_outcome("Slow", all_of("SLOW", "GATE")) |>
  build_model()

example_17_k_of_n_inhibitors <- race_spec() |>
  add_accumulator("a1", "lognormal") |>
  add_accumulator("a2", "lognormal") |>
  add_accumulator("a3", "lognormal") |>
  # Define A as a k-of-n pool with 2 members
  add_pool("A", c("a1", "a2", "a3"), k = 2L) |>
  add_accumulator("b1", "lognormal") |>
  add_accumulator("b2", "lognormal") |>
  add_accumulator("b3", "lognormal") |>
  # Define B as a k-of-n pool with 2 members
  add_pool("B", c("b1", "b2", "b3"), k = 2L) |>
  add_pool("B_other_k1", c("b2", "b3"), k = 1L) |>
  add_pool("B_no_b1", c("b2", "b3"), k = 2L) |>
  add_accumulator("c1", "lognormal") |>
  add_accumulator("c2", "lognormal") |>
  add_accumulator("c3", "lognormal") |>
  # Define C as a k-of-n pool with 2 members
  add_pool("C", c("c1", "c2", "c3"), k = 2L) |>
  add_pool("C_k1_c2c3", c("c2", "c3"), k = 1L) |>
  add_pool("C_k1_c1c3", c("c1", "c3"), k = 1L) |>
  add_accumulator("d1", "lognormal") |>
  add_accumulator("d2", "lognormal") |>
  add_accumulator("d3", "lognormal") |>
  # Define D as a k-of-n pool with 2 members
  add_pool("D", c("d1", "d2", "d3"), k = 2L) |>
  add_outcome("A", "A") |>
  add_outcome("B", first_of(
    inhibit(all_of("b1", "B_other_k1"), by = "a1"),
    "B_no_b1"
  )) |>
  add_outcome("C", first_of(
    inhibit(all_of("c1", "C_k1_c2c3"), by = "b2"),
    inhibit(all_of("c2", "C_k1_c1c3"), by = "b2")
  )) |>
  add_outcome("D", inhibit("D", by = "c3")) |>
  add_group("par:shared_A", members = c("a1", "a2", "a3"),
            attrs = list(shared_params = list(meanlog = log(0.30), sdlog = 0.15))) |>
  add_group("par:shared_B", members = c("b1", "b2", "b3"),
            attrs = list(shared_params = list(meanlog = log(0.30), sdlog = 0.15))) |>
  add_group("par:shared_C", members = c("c1", "c2", "c3"),
            attrs = list(shared_params = list(meanlog = log(0.30), sdlog = 0.15))) |>
  add_group("par:shared_D", members = c("d1", "d2", "d3"),
            attrs = list(shared_params = list(meanlog = log(0.30), sdlog = 0.15))) |>
  build_model()

example_18_shared_triggers <- race_spec() |>
  add_accumulator("go_left", "lognormal",
                  params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("go_right", "lognormal",
                  params = list(meanlog = log(0.32), sdlog = 0.20)) |>
  add_accumulator("stop_fast", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.16), sdlog = 0.15)) |>
  add_accumulator("stop_slow", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.18), sdlog = 0.15)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("STOP", c("stop_fast", "stop_slow")) |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("STOP", "STOP", options = list(map_outcome_to = NA_character_)) |>
  add_group("stop_fast_component", members = "stop_fast", attrs = list(component = "fast")) |>
  add_group("stop_slow_component", members = "stop_slow", attrs = list(component = "slow")) |>
  add_group("shared_stop_trigger",
            members = c("stop_fast", "stop_slow"),
            attrs = list(shared_trigger = list(id = "stop_gate", q = 0.25))) |>
  set_metadata(mixture = list(
    components = list(
      component("fast", weight = 0.25),
      component("slow", weight = 0.75)
    )
  )) |>
  build_model()


new_api_examples <- list(
  example_1_simple = example_1_simple,
  example_2_stop_mixture = example_2_stop_mixture,
  example_3_stop_na = example_3_stop_na,
  example_4_two_on_one = example_4_two_on_one,
  example_5_timeout_guess = example_5_timeout_guess,
  example_6_dual_path = example_6_dual_path,
  example_7_mixture = example_7_mixture,
  example_8_shared_params = example_8_shared_params,
  example_9_advanced_k = example_9_advanced_k,
  example_10_exclusion = example_10_exclusion,
  example_11_censor_deadline = example_11_censor_deadline,
  example_12_inhibitor_with_protector = example_12_inhibitor_with_protector,
  example_13_nested_pools = example_13_nested_pools,
  example_14_weighted_pool = example_14_weighted_pool,
  example_15_component_metadata = example_15_component_metadata,
  example_16_guard_tie_simple = example_16_guard_tie_simple,
  example_17_k_of_n_inhibitors = example_17_k_of_n_inhibitors,
  example_18_shared_triggers = example_18_shared_triggers
)
