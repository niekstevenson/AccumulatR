`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.integrate_rel_tol <- function() getOption("uuber.integrate.rel.tol", 1e-5)
.integrate_abs_tol <- function() getOption("uuber.integrate.abs.tol", 1e-6)

.expr_is_event <- function(expr) {
  !is.null(expr) && !is.null(expr[["kind"]]) && identical(expr[["kind"]], "event")
}

.collect_guards <- function(expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  if (identical(kind, "guard")) return(list(expr))
  if (kind %in% c("and", "or")) {
    return(unlist(lapply(expr[["args"]], .collect_guards), recursive = FALSE))
  }
  if (identical(kind, "not")) {
    return(.collect_guards(expr[["arg"]]))
  }
  list()
}

.guards_conflict <- function(primary_guard, other_expr) {
  primary_blocker <- primary_guard[["blocker"]]
  primary_reference <- primary_guard[["reference"]]
  if (!.expr_is_event(primary_blocker) || !.expr_is_event(primary_reference)) return(FALSE)
  if (length(primary_guard[["unless"]] %||% list()) > 0) return(FALSE)
  guards <- .collect_guards(other_expr)
  if (length(guards) == 0) return(FALSE)
  blk_src <- primary_blocker[["source"]]
  ref_src <- primary_reference[["source"]]
  for (g in guards) {
    if (length(g[["unless"]] %||% list()) > 0) next
    g_blocker <- g[["blocker"]]
    g_reference <- g[["reference"]]
    if (!.expr_is_event(g_blocker) || !.expr_is_event(g_reference)) next
    if (identical(blk_src, g_reference[["source"]]) &&
        identical(ref_src, g_blocker[["source"]])) {
      return(TRUE)
    }
  }
  FALSE
}

.expr_signature <- function(expr) {
  cached <- attr(expr, ".lik_signature", exact = TRUE)
  if (!is.null(cached)) return(cached)
  if (is.null(expr) || is.null(expr[["kind"]])) return("null")
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]] %||% ""
    kval <- expr[["k"]] %||% ""
    return(sprintf("event:%s:%s", as.character(src), as.character(kval)))
  }
  if (kind %in% c("and", "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0L) {
      return(sprintf("%s:empty", kind))
    }
    parts <- vapply(args, .expr_signature, character(1))
    return(sprintf("%s[%s]", kind, paste(parts, collapse = ",")))
  }
  if (identical(kind, "guard")) {
    blk_sig <- .expr_signature(expr[["blocker"]])
    ref_sig <- .expr_signature(expr[["reference"]])
    unl <- expr[["unless"]] %||% list()
    unl_sig <- if (length(unl) == 0L) "" else paste(vapply(unl, .expr_signature, character(1)), collapse = ";")
    return(sprintf("guard[%s|%s|%s]", blk_sig, ref_sig, unl_sig))
  }
  if (identical(kind, "not")) {
    return(sprintf("not(%s)", .expr_signature(expr[["arg"]])))
  }
  sprintf("other:%s", kind)
}

.compile_once <- function(fn) {
  if (!is.function(fn)) return(NULL)
  if (!isTRUE(getOption("uuber.compile.likelihood", TRUE))) return(fn)
  if (!requireNamespace("compiler", quietly = TRUE)) return(fn)
  compiler::cmpfun(fn)
}

.make_scenario_record <- function(prep, weight, forced_complete, forced_survive) {
  if (!is.finite(weight) || weight <= 0) return(NULL)
  list(
    weight = as.numeric(weight),
    forced_complete = .coerce_forced_ids(prep, forced_complete),
    forced_survive = .coerce_forced_ids(prep, forced_survive)
  )
}

.prep_runtime_get <- function(prep, key, default = NULL) {
  runtime <- prep[[".runtime"]]
  if (is.null(runtime)) return(default)
  val <- runtime[[key]]
  if (is.null(val)) default else val
}

.prep_expr_compiled <- function(prep) .prep_runtime_get(prep, "expr_compiled")
.prep_label_cache <- function(prep) .prep_runtime_get(prep, "label_cache")
.prep_competitors <- function(prep) .prep_runtime_get(prep, "competitor_map")
.prep_id_index <- function(prep) .prep_runtime_get(prep, "id_index")
