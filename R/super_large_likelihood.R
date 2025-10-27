# Likelihood system for the new expression-based API (new_API.R)
# Matches the generator logic in generator_new.R

source("R/dist.R")
source("R/utils.R")
source("R/pool_math.R")
source("R/model_tables.R")

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.integrate_rel_tol <- function() getOption("uuber.integrate.rel.tol", 1e-5)
.integrate_abs_tol <- function() getOption("uuber.integrate.abs.tol", 1e-6)

.time_key <- function(x) {
  vapply(x, function(xx) {
    if (is.na(xx)) return("NA")
    if (!is.finite(xx)) {
      if (xx > 0) return("Inf")
      if (xx < 0) return("-Inf")
      return("NA")
    }
    sprintf("%.15g", xx)
  }, character(1), USE.NAMES = FALSE)
}

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
  if (is.null(expr) || is.null(expr[["kind"]])) return("null")
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]] %||% ""
    kval <- expr[["k"]] %||% ""
    return(sprintf("event:%s:%s", as.character(src), as.character(kval)))
  }
  if (kind %in% c("and", "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0L) return(sprintf("%s:empty", kind))
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

.likelihood_cache_env <- function(cache = NULL) {
  if (is.environment(cache)) return(cache)
  new.env(parent = emptyenv(), hash = TRUE)
}

.forced_key <- function(ids) {
  if (length(ids) == 0L) return(".")
  paste(as.integer(ids), collapse = ",")
}

.likelihood_cache_key <- function(prefix, expr_id, component, t,
                                  forced_complete, forced_survive) {
  if (is.null(expr_id) || is.na(expr_id)) return(NULL)
  comp_key <- component %||% "__default__"
  t_key <- .time_key(t)[[1]]
  fc_key <- .forced_key(forced_complete)
  fs_key <- .forced_key(forced_survive)
  paste(prefix, expr_id, comp_key, t_key, fc_key, fs_key, sep = "|")
}

.cache_lookup <- function(cache, key) {
  if (is.null(cache) || is.null(key) || !is.environment(cache)) return(NULL)
  cache[[key]]
}

.cache_store <- function(cache, key, value) {
  if (is.null(cache) || is.null(key) || !is.environment(cache)) return(value)
  cache[[key]] <- value
  value
}

# ==============================================================================
# PART 1: Model Preparation
# ==============================================================================
# Reuse the preparation logic from generator_new.R to ensure consistency

.prepare_model_for_likelihood <- function(model) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  # Use the same preparation as the generator
  if (!exists("prepare_model", mode = "function")) {
    stop("prepare_model function not found - source generator_new.R first")
  }
  prep <- prepare_model(model)
  acc_ids <- names(prep[["accumulators"]] %||% list())
  pool_ids <- names(prep[["pools"]] %||% list())
  all_ids <- unique(c(acc_ids, pool_ids))
  prep[[".id_index"]] <- setNames(seq_along(all_ids), all_ids)
  if (length(all_ids) > 0) {
    idx_map <- prep[[".id_index"]]
    rev_labels <- rep(NA_character_, length(idx_map))
    rev_labels[as.integer(idx_map)] <- names(idx_map)
    prep[[".id_labels"]] <- rev_labels
  } else {
    prep[[".id_labels"]] <- character(0)
  }
  prep[[".label_cache"]] <- new.env(parent = emptyenv(), hash = TRUE)
  prep <- .precompile_likelihood_expressions(prep)
  prep[[".competitors"]] <- .prepare_competitor_map(prep)
  prep
}

.precompile_likelihood_expressions <- function(prep) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  if (length(outcome_defs) == 0L) {
    prep[[".expr_compiled"]] <- NULL
    return(prep)
  }
  sig_env <- new.env(parent = emptyenv(), hash = TRUE)
  nodes_env <- new.env(parent = emptyenv())
  nodes_env$nodes <- list()
  
  compile_expr <- function(expr) {
    if (is.null(expr) || is.null(expr[["kind"]])) return(expr)
    sig <- .expr_signature(expr)
    existing_id <- sig_env[[sig]]
    if (!is.null(existing_id)) {
      attr(expr, ".lik_id") <- as.integer(existing_id)
      return(expr)
    }
    kind <- expr[["kind"]]
    child_info <- list()
    if (kind %in% c("and", "or")) {
      args <- expr[["args"]] %||% list()
      if (length(args) > 0L) {
        expr[["args"]] <- lapply(args, compile_expr)
        child_info$args <- vapply(expr[["args"]], function(arg) {
          attr(arg, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      } else {
        child_info$args <- integer(0)
      }
    } else if (identical(kind, "guard")) {
      expr[["reference"]] <- compile_expr(expr[["reference"]])
      expr[["blocker"]] <- compile_expr(expr[["blocker"]])
      unless_list <- expr[["unless"]] %||% list()
      if (length(unless_list) > 0L) {
        expr[["unless"]] <- lapply(unless_list, compile_expr)
        child_info$unless <- vapply(expr[["unless"]], function(unl) {
          attr(unl, ".lik_id", exact = TRUE) %||% NA_integer_
        }, integer(1))
      } else {
        expr[["unless"]] <- list()
        child_info$unless <- integer(0)
      }
      child_info$reference <- attr(expr[["reference"]], ".lik_id", exact = TRUE) %||% NA_integer_
      child_info$blocker <- attr(expr[["blocker"]], ".lik_id", exact = TRUE) %||% NA_integer_
    } else if (identical(kind, "not")) {
      expr[["arg"]] <- compile_expr(expr[["arg"]])
      child_info$arg <- attr(expr[["arg"]], ".lik_id", exact = TRUE) %||% NA_integer_
    }
    node <- .build_compiled_node(expr, child_info, prep, nodes_env)
    node_id <- length(nodes_env$nodes) + 1L
    node$id <- node_id
    nodes_env$nodes[[node_id]] <- node
    sig_env[[sig]] <- node_id
    attr(expr, ".lik_id") <- node_id
    expr
  }
  
  for (lbl in names(outcome_defs)) {
    out_expr <- outcome_defs[[lbl]][["expr"]]
    outcome_defs[[lbl]][["expr"]] <- compile_expr(out_expr)
  }
  
  prep[["outcomes"]] <- outcome_defs
  prep[[".expr_compiled"]] <- list(
    nodes = nodes_env$nodes,
    signatures = sig_env
  )
  prep
}

.expr_lookup_compiled <- function(expr, prep) {
  comp <- prep[[".expr_compiled"]]
  if (is.null(comp)) return(NULL)
  node_id <- attr(expr, ".lik_id", exact = TRUE)
  if (is.null(node_id) || is.na(node_id)) {
    sig <- .expr_signature(expr)
    node_id <- comp[["signatures"]][[sig]]
  }
  if (is.null(node_id) || is.na(node_id)) return(NULL)
  node <- comp[["nodes"]][[as.integer(node_id)]]
  if (is.null(node)) return(NULL)
  node
}

.build_compiled_node <- function(expr, child_info, prep, nodes_env) {
  kind <- expr[["kind"]]
  node <- list(
    kind = kind,
    sources = .expr_sources(expr, prep),
    expr = expr
  )
  if (identical(kind, "event")) {
    node$cdf_fn <- .compile_once(.make_event_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_event_survival_fn(expr, prep))
    node$scenario_fn <- .compile_once(.make_event_scenario_fn(expr, prep))
    return(node)
  }
  if (kind %in% c("and", "or")) {
    child_ids <- as.integer(child_info$args %||% integer(0))
    node$args <- child_ids
    if (length(child_ids) == 0L || any(is.na(child_ids))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$scenario_fn <- NULL
      return(node)
    }
    child_nodes <- lapply(child_ids, function(id) nodes_env$nodes[[id]])
    if (any(vapply(child_nodes, is.null, logical(1)))) {
      node$cdf_fn <- NULL
      node$surv_fn <- NULL
      node$scenario_fn <- NULL
      return(node)
    }
    if (identical(kind, "and")) {
      have_child_cdf <- all(vapply(child_nodes, function(n) is.function(n$cdf_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      node$cdf_fn <- if (have_child_cdf) .compile_once(.make_and_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (is.function(node$cdf_fn)) .compile_once(.make_and_survival_fn(node$cdf_fn)) else NULL
      node$scenario_fn <- if (have_child_cdf && have_child_scen) .compile_once(.make_and_scenario_fn(child_nodes, prep)) else NULL
    } else {
      have_child_surv <- all(vapply(child_nodes, function(n) is.function(n$surv_fn), logical(1)))
      have_child_scen <- all(vapply(child_nodes, function(n) is.function(n$scenario_fn), logical(1)))
      node$cdf_fn <- if (have_child_surv) .compile_once(.make_or_cdf_fn(child_nodes, prep)) else NULL
      node$surv_fn <- if (have_child_surv) .compile_once(.make_or_survival_fn(child_nodes, prep)) else NULL
      node$scenario_fn <- if (have_child_scen) .compile_once(.make_or_scenario_fn(child_nodes, prep)) else NULL
    }
    return(node)
  }
  if (identical(kind, "guard")) {
    node$reference_id <- child_info$reference %||% NA_integer_
    node$blocker_id <- child_info$blocker %||% NA_integer_
    node$unless_ids <- as.integer(child_info$unless %||% integer(0))
    node$cdf_fn <- .compile_once(.make_guard_cdf_fn(expr, prep))
    node$surv_fn <- .compile_once(.make_guard_survival_fn(expr, prep))
    node$scenario_fn <- NULL
    return(node)
  }
  node$cdf_fn <- NULL
  node$surv_fn <- NULL
  node$scenario_fn <- NULL
  node
}

.make_event_cdf_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    .event_cdf_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  }
}

.make_event_survival_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    .event_survival_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  }
}

.make_event_scenario_fn <- function(expr, prep) {
  source_id <- expr[["source"]]
  source_idx <- .label_to_id(prep, source_id)
  source_ids <- if (is.na(source_idx)) integer(0) else source_idx
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    weight <- .event_density_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
    if (!is.finite(weight) || weight <= 0) return(list())
    pool_scenarios <- attr(weight, "scenarios")
    if (!is.null(pool_scenarios) && length(pool_scenarios) > 0L) {
      out <- list()
      for (psc in pool_scenarios) {
        if (is.null(psc)) next
        w <- psc$weight
        if (!is.finite(w) || w <= 0) next
        fcomp <- .forced_union(prep, forced_complete, c(source_ids, psc$forced_complete))
        fsurv <- .forced_union(prep, forced_survive, psc$forced_survive)
        sc <- .make_scenario_record(prep, w, fcomp, fsurv)
        if (!is.null(sc)) out[[length(out) + 1L]] <- sc
      }
      return(out)
    }
    sc <- .make_scenario_record(
      prep,
      weight,
      .forced_union(prep, forced_complete, source_ids),
      forced_survive
    )
    if (is.null(sc)) list() else list(sc)
  }
}

.make_and_cdf_fn <- function(child_nodes, prep) {
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fn")
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    if (length(child_cdf) == 0L) return(0.0)
    prod_val <- 1.0
    for (fn in child_cdf) {
      val <- fn(t, component, forced_complete, forced_survive)
      prod_val <- prod_val * val
      if (prod_val == 0) break
    }
    prod_val
  }
}

.make_and_survival_fn <- function(cdf_fn) {
  function(t, component, forced_complete, forced_survive) {
    val <- cdf_fn(t, component, forced_complete, forced_survive)
    if (!is.finite(val)) val <- 0.0
    out <- 1.0 - val
    if (out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    out
  }
}

.make_or_survival_fn <- function(child_nodes, prep) {
  child_surv <- lapply(child_nodes, `[[`, "surv_fn")
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    if (length(child_surv) == 0L) return(1.0)
    prod_val <- 1.0
    for (fn in child_surv) {
      val <- fn(t, component, forced_complete, forced_survive)
      prod_val <- prod_val * val
      if (prod_val == 0) break
    }
    prod_val
  }
}

.make_or_cdf_fn <- function(child_nodes, prep) {
  surv_fn <- .make_or_survival_fn(child_nodes, prep)
  function(t, component, forced_complete, forced_survive) {
    val <- surv_fn(t, component, forced_complete, forced_survive)
    out <- 1.0 - val
    if (!is.finite(out)) out <- 0.0
    if (out < 0) out <- 0.0
    if (out > 1) out <- 1.0
    out
  }
}

.make_and_scenario_fn <- function(child_nodes, prep) {
  child_scen <- lapply(child_nodes, `[[`, "scenario_fn")
  child_cdf <- lapply(child_nodes, `[[`, "cdf_fn")
  child_sources <- lapply(child_nodes, function(node) node$sources %||% integer(0))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    out <- list()
    n <- length(child_scen)
    if (n == 0L) return(out)
    for (i in seq_len(n)) {
      si <- child_scen[[i]](t, component, forced_complete, forced_survive)
      if (length(si) == 0L) next
      others <- seq_len(n)
      if (length(others) > 0L) others <- others[others != i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        ok <- TRUE
        if (length(others) > 0L) {
          for (j in others) {
            Fj <- child_cdf[[j]](t, component, fcomp, fsurv)
            if (!is.finite(Fj) || Fj <= 0) {
              ok <- FALSE
              break
            }
            weight <- weight * Fj
            if (!is.finite(weight) || weight <= 0) {
              ok <- FALSE
              break
            }
            fcomp <- .forced_union(prep, fcomp, child_sources[[j]])
          }
        }
        if (!ok) next
        sc_out <- .make_scenario_record(prep, weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    out
  }
}

.make_or_scenario_fn <- function(child_nodes, prep) {
  child_scen <- lapply(child_nodes, `[[`, "scenario_fn")
  child_sources <- lapply(child_nodes, function(node) node$sources %||% integer(0))
  function(t, component, forced_complete, forced_survive) {
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    out <- list()
    n <- length(child_scen)
    if (n == 0L) return(out)
    for (i in seq_len(n)) {
      si <- child_scen[[i]](t, component, forced_complete, forced_survive)
      if (length(si) == 0L) next
      others <- seq_len(n)
      if (length(others) > 0L) others <- others[others != i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        valid <- TRUE
        if (length(others) > 0L) {
          for (j in others) {
            req_j <- child_sources[[j]]
            if (length(req_j) == 0L) next
            if (all(req_j %in% fcomp)) {
              valid <- FALSE
              break
            }
            witness <- req_j[!(req_j %in% fcomp)]
            if (length(witness) >= 1L) {
              fsurv <- .forced_union(prep, fsurv, witness[[1]])
            } else {
              valid <- FALSE
              break
            }
          }
        }
        if (!valid) next
        sc_out <- .make_scenario_record(prep, weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    out
  }
}

.make_guard_cdf_fn <- function(expr, prep) {
  function(t, component, forced_complete, forced_survive) {
    if (t <= 0) return(0.0)
    if (!is.finite(t)) return(1.0)
    forced_complete <- .coerce_forced_ids(prep, forced_complete)
    forced_survive <- .coerce_forced_ids(prep, forced_survive)
    local_cache <- .likelihood_cache_env()
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        .eval_expr_likelihood(
          expr, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = local_cache
        )
      }, numeric(1))
    }
    val <- tryCatch(
      stats::integrate(
        dens_fun,
        lower = 0,
        upper = t,
        rel.tol = .integrate_rel_tol(),
        abs.tol = .integrate_abs_tol(),
        stop.on.error = FALSE
      )[["value"]],
      error = function(e) 0.0
    )
    out <- as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
}

.make_guard_survival_fn <- function(expr, prep) {
  cdf_fn <- .make_guard_cdf_fn(expr, prep)
  function(t, component, forced_complete, forced_survive) {
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    val <- cdf_fn(t, component, forced_complete, forced_survive)
    out <- 1.0 - val
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
}

.label_to_id <- function(prep, label) {
  if (is.null(label) || length(label) == 0L) return(NA_integer_)
  idx_map <- prep[[".id_index"]]
  if (is.null(idx_map)) return(NA_integer_)
  val <- idx_map[[as.character(label)]]
  if (is.null(val) || is.na(val)) return(NA_integer_)
  as.integer(val)
}

.labels_to_ids <- function(prep, labels) {
  if (is.null(labels) || length(labels) == 0L) return(integer(0))
  idx_map <- prep[[".id_index"]]
  if (is.null(idx_map)) return(integer(0))
  chr <- as.character(labels)
  if (length(chr) == 0L) return(integer(0))
  cache_env <- prep[[".label_cache"]]
  cache_key <- NULL
  if (!is.null(cache_env)) {
    cache_key <- paste(chr, collapse = "\r")
    cached <- cache_env[[cache_key]]
    if (!is.null(cached)) return(cached)
  }
  vals <- idx_map[chr]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(integer(0))
  ids <- as.integer(vals)
  if (length(ids) > 1L) ids <- sort(unique(ids))
  if (!is.null(cache_env) && nzchar(cache_key %||% "")) {
    cache_env[[cache_key]] <- ids
  }
  ids
}

.ids_to_labels <- function(prep, ids) {
  if (is.null(ids) || length(ids) == 0L) return(character(0))
  labels <- prep[[".id_labels"]]
  if (is.null(labels) || length(labels) == 0L) return(character(0))
  ids <- as.integer(ids)
  out <- labels[ids]
  out[!is.na(out)]
}

.coerce_forced_ids <- function(prep, values) {
  if (is.null(values) || length(values) == 0L) return(integer(0))
  if (is.integer(values)) {
    if (length(values) > 1L) return(sort(unique(values)))
    return(values)
  }
  .labels_to_ids(prep, values)
}

.forced_union <- function(prep, base_set, additions) {
  base_ids <- .coerce_forced_ids(prep, base_set)
  add_ids <- .coerce_forced_ids(prep, additions)
  if (length(add_ids) == 0L) return(base_ids)
  if (length(base_ids) == 0L) return(add_ids)
  out <- c(base_ids, add_ids)
  out <- as.integer(out)
  if (length(out) > 1L) out <- sort(unique(out))
  out
}

# ==============================================================================
# PART 2: Distribution Primitives
# ==============================================================================

# Get density function for an accumulator at time t
.acc_density <- function(acc_def, t) {
  if (!is.finite(t) || t < 0) return(0.0)
  if (t < acc_def[['onset']]) return(0.0)
  
  # Account for q (probability of non-completion)
  success_prob <- 1 - acc_def[['q']]
  if (success_prob <= 0) return(0.0)
  
  reg <- dist_registry(acc_def[['dist']])
  if (is.null(reg) || is.null(reg[['d']])) {
    stop(sprintf("No density function for distribution '%s'", acc_def[['dist']]))
  }
  
  # Density of the underlying distribution at (t - onset)
  t_adj <- t - acc_def[['onset']]
  dens <- reg[['d']](t_adj, acc_def[['params']])
  
  # Scale by success probability
  success_prob * dens
}

# Get survival function for an accumulator at time t
.acc_survival <- function(acc_def, t) {
  if (!is.finite(t)) return(0.0)
  if (t < 0) return(1.0)
  if (t < acc_def[['onset']]) return(1.0)
  
  # Survival = q + (1-q) * S(t - onset)
  reg <- dist_registry(acc_def[['dist']])
  if (is.null(reg) || is.null(reg[['p']])) {
    stop(sprintf("No CDF function for distribution '%s'", acc_def[['dist']]))
  }
  
  t_adj <- t - acc_def[['onset']]
  surv_underlying <- 1 - reg[['p']](t_adj, acc_def[['params']])
  
  acc_def[['q']] + (1 - acc_def[['q']]) * surv_underlying
}

# ==============================================================================
# PART 3: Pool Likelihood (Scenario-based)
# ==============================================================================

# Internal: return active members for a pool in this component
.pool_active_members <- function(prep, pool_id, component) {
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  members <- pool_def[["members"]]
  active <- character(0)
  for (m in members) {
    if (!is.null(acc_defs[[m]])) {
      comps <- acc_defs[[m]][["components"]]
      if (length(comps) == 0 || is.null(component) ||
          identical(component, "__default__") || component %in% comps) {
        active <- c(active, m)
      }
    } else if (!is.null(pool_defs[[m]])) {
      active <- c(active, m)
    }
  }
  active
}

.pool_density_fast_value <- function(prep, pool_id, component, t) {
  if (is.na(t) || !is.finite(t) || t < 0) return(0.0)
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) return(0.0)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(0.0)
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k < 1L || k > n) return(0.0)
  
  total_density <- 0.0
  for (i in seq_len(n)) {
    mid <- members[[i]]
    dens_i <- if (!is.null(acc_defs[[mid]])) {
      .acc_density(acc_defs[[mid]], t)
    } else if (!is.null(pool_defs[[mid]])) {
      .pool_density_fast_value(prep, mid, component, t)
    } else {
      0.0
    }
    if (!is.finite(dens_i) || dens_i <= 0) next
    others <- members[-i]
    if (length(others) == 0L) {
      if (k == 1L) total_density <- total_density + dens_i
      next
    }
    Svec <- Fvec <- numeric(length(others))
    for (j in seq_along(others)) {
      oid <- others[[j]]
      surv_j <- if (!is.null(acc_defs[[oid]])) {
        .acc_survival(acc_defs[[oid]], t)
      } else if (!is.null(pool_defs[[oid]])) {
        .pool_survival_fast_value(prep, oid, component, t)
      } else {
        1.0
      }
      Svec[[j]] <- surv_j
      Fvec[[j]] <- 1.0 - surv_j
    }
    coeffs <- pool_coeffs(Svec, Fvec)
    if (k <= length(coeffs)) {
      total_density <- total_density + dens_i * coeffs[[k]]
    }
  }
  total_density
}

.pool_survival_fast_value <- function(prep, pool_id, component, t) {
  if (is.na(t)) return(1.0)
  pool_defs <- prep[["pools"]]
  acc_defs <- prep[["accumulators"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) return(1.0)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(1.0)
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k > n) return(1.0)
  if (k < 1L) return(0.0)
  
  if (!is.finite(t)) {
    if (t < 0) return(1.0)
    t <- Inf
  }
  
  Svec <- numeric(n)
  Fvec <- numeric(n)
  for (i in seq_len(n)) {
    mid <- members[[i]]
    surv <- if (!is.null(acc_defs[[mid]])) {
      .acc_survival(acc_defs[[mid]], t)
    } else if (!is.null(pool_defs[[mid]])) {
      .pool_survival_fast_value(prep, mid, component, t)
    } else {
      1.0
    }
    Svec[[i]] <- surv
    Fvec[[i]] <- 1.0 - surv
  }
  coeffs <- pool_coeffs(Svec, Fvec)
  upto <- min(length(coeffs), k)
  if (upto <= 0L) return(0.0)
  sum(coeffs[seq_len(upto)])
}

.build_pool_templates <- function(pool_id, members, member_ids, pool_idx, k) {
  n <- length(members)
  if (n == 0L || k < 1L || k > n) {
    return(list(templates = list(), finisher_map = vector("list", n)))
  }
  need <- k - 1L
  finisher_map <- vector("list", n)
  templates <- list()
  template_count <- 0L
  for (idx in seq_len(n)) {
    others <- seq_len(n)
    others <- others[others != idx]
    if (need > length(others)) {
      finisher_map[[idx]] <- integer(0)
      next
    }
    combos <- if (need <= 0L) {
      list(integer(0))
    } else if (need == length(others)) {
      list(others)
    } else if (need == 1L) {
      lapply(others, function(j) j)
    } else {
      utils::combn(others, need, simplify = FALSE)
    }
    idx_entries <- integer(length(combos))
    for (j in seq_along(combos)) {
      combo_idx <- combos[[j]]
      survivors <- if (length(combo_idx) == 0L) others else {
        remainder <- others[!(others %in% combo_idx)]
        if (length(remainder) == 0L) integer(0) else remainder
      }
      template_count <- template_count + 1L
      finisher_ids <- member_ids[idx]
      if (length(combo_idx) > 0L) finisher_ids <- c(finisher_ids, member_ids[combo_idx])
      if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
      finisher_ids <- finisher_ids[!is.na(finisher_ids)]
      templates[[template_count]] <- list(
        finisher_idx = idx,
        complete_idx = if (length(combo_idx) > 0L) combo_idx else integer(0),
        survivor_idx = survivors,
        forced_complete_ids = finisher_ids,
        forced_survive_ids = if (length(survivors) > 0L) member_ids[survivors] else integer(0)
      )
      idx_entries[[j]] <- template_count
    }
    finisher_map[[idx]] <- idx_entries
  }
  list(templates = templates, finisher_map = finisher_map)
}

.event_survival_at <- function(prep, id, component, t,
                               forced_complete = integer(0),
                               forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(1.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    return(.pool_survival(prep, id, component, t,
                          forced_complete = forced_complete,
                          forced_survive = forced_survive))
  }
  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) return(.acc_survival(acc_def, t))
  0.0
}

.event_cdf_at <- function(prep, id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  1.0 - .event_survival_at(prep, id, component, t,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive)
}

.event_density_at <- function(prep, id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0)) {
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(0.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    return(.pool_density(prep, id, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive))
  }
  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) return(.acc_density(acc_def, t))
  0.0
}

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  if (!is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  k <- as.integer(pool_def[["k"]] %||% 1L)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L || k < 1L || k > n) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  member_ids <- .labels_to_ids(prep, members)
  pool_idx <- .label_to_id(prep, pool_id)
  
  idx_seq <- seq_len(n)
  dens_vec <- numeric(n)
  cdf_vec <- numeric(n)
  surv_vec <- numeric(n)
  for (i in idx_seq) {
    mid <- members[[i]]
    dens_vec[[i]] <- .event_density_at(
      prep, mid, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    cdf_vec[[i]] <- .event_cdf_at(
      prep, mid, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    surv_vec[[i]] <- 1.0 - cdf_vec[[i]]
  }
  
  scenarios <- list()
  add_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return()
    scenarios[[length(scenarios) + 1L]] <<- list(
      weight = as.numeric(weight),
      forced_complete = .coerce_forced_ids(prep, fcomp),
      forced_survive = .coerce_forced_ids(prep, fsurv)
    )
  }
  
  fast_mode <- (length(forced_complete) == 0L && length(forced_survive) == 0L)
  fast_done <- FALSE
  if (fast_mode) {
    if (k == 1L) {
      if (n == 1L) {
        weight <- dens_vec[[1]]
        if (is.finite(weight) && weight > 0) {
          finisher_ids <- member_ids[1]
          if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
          finisher_ids <- finisher_ids[!is.na(finisher_ids)]
          add_scenario(weight, finisher_ids, integer(0))
          fast_done <- TRUE
        }
      } else {
        for (idx in idx_seq) {
          dens_mid <- dens_vec[[idx]]
          if (!is.finite(dens_mid) || dens_mid <= 0) next
          survivors_idx <- idx_seq[idx_seq != idx]
          weight <- dens_mid
          if (length(survivors_idx) > 0L) {
            weight <- weight * prod(surv_vec[survivors_idx])
          }
          if (!is.finite(weight) || weight <= 0) next
          finisher_ids <- member_ids[idx]
          if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
          finisher_ids <- finisher_ids[!is.na(finisher_ids)]
          survivor_ids <- if (length(survivors_idx) > 0L) member_ids[survivors_idx] else integer(0)
          add_scenario(weight, finisher_ids, survivor_ids)
          fast_done <- TRUE
        }
      }
    } else if (k == n) {
      for (idx in idx_seq) {
        dens_mid <- dens_vec[[idx]]
        if (!is.finite(dens_mid) || dens_mid <= 0) next
        others_idx <- idx_seq[idx_seq != idx]
        weight <- dens_mid
        if (length(others_idx) > 0L) {
          weight <- weight * prod(cdf_vec[others_idx])
        }
        if (!is.finite(weight) || weight <= 0) next
        finisher_ids <- member_ids
        if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
        finisher_ids <- finisher_ids[!is.na(finisher_ids)]
        add_scenario(weight, finisher_ids, integer(0))
        fast_done <- TRUE
      }
    } else {
      template_info <- .build_pool_templates(pool_id, members, member_ids, pool_idx, k)
      if (!is.null(template_info) && length(template_info$templates) > 0L) {
        templates <- template_info$templates
        finisher_map <- template_info$finisher_map
        for (idx in idx_seq) {
          dens_mid <- dens_vec[[idx]]
          if (!is.finite(dens_mid) || dens_mid <= 0) next
          template_ids <- finisher_map[[idx]]
          if (length(template_ids) == 0L) next
          for (tid in template_ids) {
            tmpl <- templates[[tid]]
            weight <- dens_mid
            if (length(tmpl$complete_idx) > 0L) {
              weight <- weight * prod(cdf_vec[tmpl$complete_idx])
            }
            if (length(tmpl$survivor_idx) > 0L) {
              weight <- weight * prod(surv_vec[tmpl$survivor_idx])
            }
            if (!is.finite(weight) || weight <= 0) next
            add_scenario(weight, tmpl$forced_complete_ids, tmpl$forced_survive_ids)
            fast_done <- TRUE
          }
        }
      }
    }
  }
  
  if (!fast_mode || !fast_done) {
    others_idx_list <- lapply(idx_seq, function(i) idx_seq[idx_seq != i])
    forced_complete_flags <- if (length(forced_complete) > 0L) member_ids %in% forced_complete else rep(FALSE, n)
    forced_survive_flags <- if (length(forced_survive) > 0L) member_ids %in% forced_survive else rep(FALSE, n)
    free_mask_base <- !(forced_complete_flags | forced_survive_flags)
    
    for (idx in idx_seq) {
      if (forced_survive_flags[[idx]]) next
      dens_mid <- dens_vec[[idx]]
      if (!is.finite(dens_mid) || dens_mid <= 0) next
      
      need <- k - 1L
      others_idx <- others_idx_list[[idx]]
      if (need > length(others_idx)) next
      
      fc_idx <- if (length(others_idx) > 0L) others_idx[forced_complete_flags[others_idx]] else integer(0)
      fs_idx <- if (length(others_idx) > 0L) others_idx[forced_survive_flags[others_idx]] else integer(0)
      forced_count <- length(fc_idx)
      if (forced_count > need) next
      remaining_need <- need - forced_count
      
      free_idx <- if (length(others_idx) > 0L) others_idx[free_mask_base[others_idx]] else integer(0)
      if (remaining_need > length(free_idx)) next
      
      base_factor <- dens_mid
      if (length(fc_idx) > 0L) base_factor <- base_factor * prod(cdf_vec[fc_idx])
      if (length(fs_idx) > 0L) base_factor <- base_factor * prod(surv_vec[fs_idx])
      if (!is.finite(base_factor) || base_factor <= 0) next
      
      finisher_ids <- member_ids[idx]
      if (length(fc_idx) > 0L) finisher_ids <- c(finisher_ids, member_ids[fc_idx])
      if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
      forced_complete_ids <- .forced_union(prep, forced_complete, finisher_ids)
      base_survive_ids <- .forced_union(prep, forced_survive, member_ids[fs_idx])
      
      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
          )
        }
        next
      }
      
      combos <- if (remaining_need == length(free_idx)) {
        list(free_idx)
      } else if (remaining_need == 1L) {
        lapply(free_idx, function(i) i)
      } else if (remaining_need <= 0L) {
        list(integer(0))
      } else {
        utils::combn(free_idx, remaining_need, simplify = FALSE)
      }
      
      if (length(combos) == 0L) next
      
      for (combo_idx in combos) {
        if (length(combo_idx) != remaining_need) next
        
        combo_weight <- base_factor
        if (length(combo_idx) > 0L) combo_weight <- combo_weight * prod(cdf_vec[combo_idx])
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        if (length(free_idx) > 0L) {
          drop_pos <- match(combo_idx, free_idx, nomatch = 0L)
          drop_pos <- drop_pos[drop_pos > 0L]
          survivors_idx <- if (length(drop_pos) > 0L) free_idx[-drop_pos] else free_idx
        } else {
          survivors_idx <- integer(0)
        }
        
        if (length(survivors_idx) > 0L) {
          combo_weight <- combo_weight * prod(surv_vec[survivors_idx])
        }
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        add_scenario(
          combo_weight,
          .forced_union(prep, forced_complete_ids, member_ids[combo_idx]),
          .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
        )
      }
    }
  }
  
  if (!fast_done || !fast_mode) {
    others_idx_list <- lapply(idx_seq, function(i) idx_seq[idx_seq != i])
    forced_complete_flags <- if (length(forced_complete) > 0L) member_ids %in% forced_complete else rep(FALSE, n)
    forced_survive_flags <- if (length(forced_survive) > 0L) member_ids %in% forced_survive else rep(FALSE, n)
    free_mask_base <- !(forced_complete_flags | forced_survive_flags)
    
    for (idx in idx_seq) {
      if (forced_survive_flags[[idx]]) next
      dens_mid <- dens_vec[[idx]]
      if (!is.finite(dens_mid) || dens_mid <= 0) next
      
      need <- k - 1L
      others_idx <- others_idx_list[[idx]]
      if (need > length(others_idx)) next
      
      fc_idx <- if (length(others_idx) > 0L) others_idx[forced_complete_flags[others_idx]] else integer(0)
      fs_idx <- if (length(others_idx) > 0L) others_idx[forced_survive_flags[others_idx]] else integer(0)
      forced_count <- length(fc_idx)
      if (forced_count > need) next
      remaining_need <- need - forced_count
      
      free_idx <- if (length(others_idx) > 0L) others_idx[free_mask_base[others_idx]] else integer(0)
      if (remaining_need > length(free_idx)) next
      
      base_factor <- dens_mid
      if (length(fc_idx) > 0L) base_factor <- base_factor * prod(cdf_vec[fc_idx])
      if (length(fs_idx) > 0L) base_factor <- base_factor * prod(surv_vec[fs_idx])
      if (!is.finite(base_factor) || base_factor <= 0) next
      
      finisher_ids <- member_ids[idx]
      if (length(fc_idx) > 0L) finisher_ids <- c(finisher_ids, member_ids[fc_idx])
      if (!is.na(pool_idx)) finisher_ids <- c(finisher_ids, pool_idx)
      forced_complete_ids <- .forced_union(prep, forced_complete, finisher_ids)
      base_survive_ids <- .forced_union(prep, forced_survive, member_ids[fs_idx])
      
      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
          )
        }
        next
      }
      
      combos <- if (remaining_need == length(free_idx)) {
        list(free_idx)
      } else if (remaining_need == 1L) {
        lapply(free_idx, function(i) i)
      } else if (remaining_need <= 0L) {
        list(integer(0))
      } else {
        utils::combn(free_idx, remaining_need, simplify = FALSE)
      }
      
      if (length(combos) == 0L) next
      
      for (combo_idx in combos) {
        if (length(combo_idx) != remaining_need) next
        
        combo_weight <- base_factor
        if (length(combo_idx) > 0L) combo_weight <- combo_weight * prod(cdf_vec[combo_idx])
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        if (length(free_idx) > 0L) {
          drop_pos <- match(combo_idx, free_idx, nomatch = 0L)
          drop_pos <- drop_pos[drop_pos > 0L]
          survivors_idx <- if (length(drop_pos) > 0L) free_idx[-drop_pos] else free_idx
        } else {
          survivors_idx <- integer(0)
        }
        
        if (length(survivors_idx) > 0L) {
          combo_weight <- combo_weight * prod(surv_vec[survivors_idx])
        }
        if (!is.finite(combo_weight) || combo_weight <= 0) next
        
        add_scenario(
          combo_weight,
          .forced_union(prep, forced_complete_ids, member_ids[combo_idx]),
          .forced_union(prep, base_survive_ids, member_ids[survivors_idx])
        )
      }
    }
  }
  
  total_density <- if (length(scenarios) == 0L) 0.0 else sum(vapply(scenarios, `[[`, numeric(1), "weight"))
  if (length(scenarios) > 0L) {
    attr(total_density, "scenarios") <- scenarios
  } else {
    attr(total_density, "scenarios") <- list()
  }
  total_density
}

.pool_survival <- function(prep, pool_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)) {
  if (!is.finite(t)) return(0.0)
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  k <- as.integer(pool_def[["k"]] %||% 1L)
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) return(1.0)
  if (k > n) return(1.0)
  if (k < 1L) return(0.0)
  
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  
  Fvec <- numeric(n)
  for (i in seq_len(n)) {
    id <- members[[i]]
    Fvec[[i]] <- .event_cdf_at(prep, id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive)
  }
  Svec <- 1.0 - Fvec
  coeffs <- pool_coeffs(Svec, Fvec)
  upto <- min(length(coeffs), k)
  sum(coeffs[seq_len(upto)])
}

# ==============================================================================
# PART 4: Expression Evaluation (Likelihood)
# ==============================================================================

# Scenario helpers for generic expressions at time t.
# Returns a list of scenario records:
#   list(weight = numeric(1),
#        forced_complete = integer(),
#        forced_survive  = integer())
.expr_scenarios_at <- function(expr, t, prep, component,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               cache = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- .likelihood_cache_key("scen", node_id, component, t,
                                     forced_complete, forced_survive)
  cached <- .cache_lookup(cache, cache_key)
  if (!is.null(cached)) return(cached)
  if (!is.null(compiled) && is.function(compiled$scenario_fn)) {
    res <- compiled$scenario_fn(t, component, forced_complete, forced_survive)
  } else {
    res <- .expr_scenarios_at_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) .cache_store(cache, cache_key, res) else res
}

.expr_scenarios_at_slow <- function(expr, t, prep, component,
                                    forced_complete = integer(0),
                                    forced_survive = integer(0),
                                    cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  
  make_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return(NULL)
    list(
      weight = as.numeric(weight),
      forced_complete = .coerce_forced_ids(prep, fcomp),
      forced_survive = .coerce_forced_ids(prep, fsurv)
    )
  }
  
  if (identical(kind, "event")) {
    source_id <- expr[["source"]]
    source_idx <- .label_to_id(prep, source_id)
    source_ids <- if (is.na(source_idx)) integer(0) else source_idx
    weight <- .event_density_at(
      prep, source_id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive)
    if (!is.finite(weight) || weight <= 0) return(list())
    pool_scenarios <- attr(weight, "scenarios")
    if (!is.null(pool_scenarios) && length(pool_scenarios) > 0) {
      out <- list()
      for (psc in pool_scenarios) {
        if (is.null(psc)) next
        w <- psc$weight
        if (!is.finite(w) || w <= 0) next
        fcomp <- .forced_union(prep, forced_complete, c(source_ids, psc$forced_complete))
        fsurv <- .forced_union(prep, forced_survive, psc$forced_survive)
        sc <- make_scenario(w, fcomp, fsurv)
        if (!is.null(sc)) out[[length(out) + 1L]] <- sc
      }
      return(out)
    }
    sc <- make_scenario(
      weight,
      fcomp = .forced_union(prep, forced_complete, source_ids),
      fsurv = forced_survive
    )
    return(if (is.null(sc)) list() else list(sc))
  }
  
  if (identical(kind, "and")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(list())
    out <- list()
    for (i in seq_along(args)) {
      si <- .expr_scenarios_at(
        args[[i]], t, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        cache = cache)
      if (length(si) == 0) next
      others <- args[-i]
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        ok <- TRUE
        for (aj in others) {
          Fj <- .eval_expr_cdf_cond(
            aj, t, prep, component,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            cache = cache)
          if (!is.finite(Fj) || Fj <= 0) {
            ok <- FALSE
            break
          }
          weight <- weight * Fj
          if (!is.finite(weight) || weight <= 0) {
            ok <- FALSE
            break
          }
          fcomp <- .forced_union(prep, fcomp, .expr_sources(aj, prep))
        }
        if (!ok) next
        sc_out <- make_scenario(weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    return(out)
  }
  
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(list())
    out <- list()
    for (i in seq_along(args)) {
      si <- .expr_scenarios_at(
        args[[i]], t, prep, component,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        cache = cache)
      if (length(si) == 0) next
      others <- args[-i]
      others_req <- lapply(others, function(aj) .expr_sources(aj, prep))
      for (sc in si) {
        if (is.null(sc) || sc$weight <= 0) next
        weight <- sc$weight
        fcomp <- sc$forced_complete
        fsurv <- sc$forced_survive
        valid <- TRUE
        for (idx in seq_along(others)) {
          req_j <- others_req[[idx]]
          if (length(req_j) == 0L) next
          if (all(req_j %in% fcomp)) {
            valid <- FALSE
            break
          }
          witness <- req_j[!(req_j %in% fcomp)]
          if (length(witness) >= 1L) {
            fsurv <- .forced_union(prep, fsurv, witness[[1]])
          } else {
            valid <- FALSE
            break
          }
        }
        if (!valid) next
        sc_out <- make_scenario(weight, fcomp, fsurv)
        if (!is.null(sc_out)) out[[length(out) + 1L]] <- sc_out
      }
    }
    return(out)
  }
  
  if (identical(kind, "guard")) {
    val <- .eval_expr_likelihood(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
    scen_attr <- attr(val, "scenarios")
    if (is.list(scen_attr) && length(scen_attr) > 0) {
      out <- list()
      for (sc in scen_attr) {
        if (is.null(sc)) next
        w <- sc$weight
        if (!is.finite(w) || w <= 0) next
        fc <- .forced_union(prep, forced_complete, sc$forced_complete)
        fs <- .forced_union(prep, forced_survive, sc$forced_survive)
        out[[length(out) + 1L]] <- list(
          weight = as.numeric(w),
          forced_complete = .coerce_forced_ids(prep, fc),
          forced_survive = .coerce_forced_ids(prep, fs)
        )
      }
      return(out)
    }
    if (!is.finite(val) || val <= 0) return(list())
    sc <- make_scenario(
      val,
      fcomp = .forced_union(prep, forced_complete, .expr_sources(expr, prep)),
      fsurv = forced_survive
    )
    return(if (is.null(sc)) list() else list(sc))
  }
  
  weight <- .eval_expr_likelihood(
    expr, t, prep, component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    cache = cache
  )
  if (!is.finite(weight) || weight <= 0) return(list())
  sc <- make_scenario(
    weight,
    fcomp = .forced_union(prep, forced_complete, .expr_sources(expr, prep)),
    fsurv = forced_survive
  )
  if (is.null(sc)) list() else list(sc)
}

.eval_expr_cdf_cond <- function(expr, t, prep, component,
                                forced_complete = integer(0),
                                forced_survive = integer(0),
                                cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- .likelihood_cache_key("cdf", node_id, component, t,
                                     forced_complete, forced_survive)
  cached <- .cache_lookup(cache, cache_key)
  if (!is.null(cached)) return(cached)
  res <- if (!is.null(compiled) && is.function(compiled$cdf_fn)) {
    compiled$cdf_fn(t, component, forced_complete, forced_survive)
  } else {
    .eval_expr_cdf_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) .cache_store(cache, cache_key, res) else res
}

.eval_expr_cdf_cond_slow <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    return(.event_cdf_at(prep, src, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive))
  }
  if (identical(kind, "and")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(0.0)
    prod_val <- 1.0
    for (a in args) {
      Fa <- .eval_expr_cdf_cond(a, t, prep, component,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                cache = cache)
      prod_val <- prod_val * Fa
      if (prod_val == 0) break
    }
    return(prod_val)
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(0.0)
    prod_surv <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache)
      prod_surv <- prod_surv * Sa
      if (prod_surv == 0) break
    }
    return(1.0 - prod_surv)
  }
  if (identical(kind, "guard")) {
    if (t <= 0) return(0.0)
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        .eval_expr_likelihood(
          expr, ui, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = cache
        )
      }, numeric(1))
    }
    val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                     rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                     stop.on.error = FALSE)$value,
                    error = function(e) 0.0)
    out <- as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    return(max(0.0, min(1.0, out)))
  }
  return(.eval_expr_cdf(expr, t, prep, component))
}

.eval_expr_survival_cond <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  compiled <- .expr_lookup_compiled(expr, prep)
  node_id <- if (!is.null(compiled)) compiled$id else attr(expr, ".lik_id", exact = TRUE)
  cache_key <- .likelihood_cache_key("surv", node_id, component, t,
                                     forced_complete, forced_survive)
  cached <- .cache_lookup(cache, cache_key)
  if (!is.null(cached)) return(cached)
  res <- if (!is.null(compiled) && is.function(compiled$surv_fn)) {
    compiled$surv_fn(t, component, forced_complete, forced_survive)
  } else {
    .eval_expr_survival_cond_slow(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache
    )
  }
  if (!is.null(cache_key)) .cache_store(cache, cache_key, res) else res
}

.eval_expr_survival_cond_slow <- function(expr, t, prep, component,
                                          forced_complete = integer(0),
                                          forced_survive = integer(0),
                                          cache = NULL) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  kind <- expr[["kind"]]
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    return(.event_survival_at(prep, src, component, t,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive))
  }
  if (identical(kind, "and")) {
    return(1.0 - .eval_expr_cdf_cond(expr, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache))
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(1.0)
    prod_val <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive,
                                     cache = cache)
      prod_val <- prod_val * Sa
      if (prod_val == 0) break
    }
    return(prod_val)
  }
  if (identical(kind, "guard")) {
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .eval_expr_likelihood(
          expr, uu, prep, component,
          forced_complete = forced_complete,
          forced_survive = forced_survive,
          cache = cache
        )
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                stop.on.error = FALSE)[["value"]],
               error = function(e) 0.0)
    }
    val <- compute_integral()
    out <- 1.0 - as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    return(max(0.0, min(1.0, out)))
  }
  return(.eval_expr_survival(expr, t, prep, component))
}

.expr_fastpath_supported <- function(expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(FALSE)
  kind <- expr[["kind"]]
  if (identical(kind, "event")) return(TRUE)
  if (kind %in% c("and", "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(FALSE)
    return(all(vapply(args, .expr_fastpath_supported, logical(1))))
  }
  FALSE
}

.expr_cdf_fast <- function(expr, t, prep, component) {
  if (is.na(t)) return(0.0)
  if (is.infinite(t) && t < 0) return(0.0)
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    if (is.null(src)) return(0.0)
    if (identical(src, "__DEADLINE__") || identical(src, "__GUESS__")) {
      return(if (t >= 0) 1.0 else 0.0)
    }
    if (!is.null(prep[["pools"]][[src]])) {
      surv <- .pool_survival_fast_value(prep, src, component, t)
      return(1.0 - surv)
    }
    return(.event_cdf_at(prep, src, component, t,
                         forced_complete = integer(0),
                         forced_survive = integer(0)))
  }
  args <- expr[["args"]] %||% list()
  if (length(args) == 0) return(0.0)
  if (identical(kind, "and")) {
    prod_val <- 1.0
    for (a in args) {
      Fa <- .expr_cdf_fast(a, t, prep, component)
      prod_val <- prod_val * Fa
      if (!is.finite(prod_val) || prod_val == 0) break
    }
    return(max(0.0, min(1.0, prod_val)))
  }
  if (identical(kind, "or")) {
    surv_prod <- 1.0
    for (a in args) {
      Sa <- .expr_survival_fast(a, t, prep, component)
      surv_prod <- surv_prod * Sa
      if (!is.finite(surv_prod) || surv_prod == 0) break
    }
    val <- 1.0 - surv_prod
    if (!is.finite(val)) val <- 0.0
    return(max(0.0, min(1.0, val)))
  }
  0.0
}

.expr_survival_fast <- function(expr, t, prep, component) {
  if (is.na(t)) return(1.0)
  if (is.infinite(t) && t < 0) return(1.0)
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    if (is.null(src)) return(1.0)
    if (identical(src, "__DEADLINE__") || identical(src, "__GUESS__")) {
      return(if (t < 0) 1.0 else 0.0)
    }
    if (!is.null(prep[["pools"]][[src]])) {
      return(.pool_survival_fast_value(prep, src, component, t))
    }
    return(.event_survival_at(prep, src, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0)))
  }
  args <- expr[["args"]] %||% list()
  if (length(args) == 0) return(1.0)
  if (identical(kind, "and")) {
    val <- 1.0 - .expr_cdf_fast(expr, t, prep, component)
    if (!is.finite(val)) val <- 0.0
    return(max(0.0, min(1.0, val)))
  }
  if (identical(kind, "or")) {
    prod_val <- 1.0
    for (a in args) {
      Sa <- .expr_survival_fast(a, t, prep, component)
      prod_val <- prod_val * Sa
      if (!is.finite(prod_val) || prod_val == 0) break
    }
    return(max(0.0, min(1.0, prod_val)))
  }
  1.0
}

.expr_density_fast <- function(expr, t, prep, component) {
  if (is.na(t) || !is.finite(t) || t < 0) return(0.0)
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    if (is.null(src)) return(0.0)
    if (identical(src, "__DEADLINE__") || identical(src, "__GUESS__")) return(0.0)
    if (!is.null(prep[["pools"]][[src]])) {
      return(.pool_density_fast_value(prep, src, component, t))
    }
    return(.event_density_at(prep, src, component, t,
                             forced_complete = integer(0),
                             forced_survive = integer(0)))
  }
  args <- expr[["args"]] %||% list()
  if (length(args) == 0) return(0.0)
  if (identical(kind, "and")) {
    total <- 0.0
    for (i in seq_along(args)) {
      dens_i <- .expr_density_fast(args[[i]], t, prep, component)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_cdf <- 1.0
      for (j in seq_along(args)) {
        if (i == j) next
        Fj <- .expr_cdf_fast(args[[j]], t, prep, component)
        prod_cdf <- prod_cdf * Fj
        if (!is.finite(prod_cdf) || prod_cdf == 0) break
      }
      if (!is.finite(prod_cdf) || prod_cdf == 0) next
      total <- total + dens_i * prod_cdf
    }
    if (!is.finite(total)) total <- 0.0
    return(max(0.0, total))
  }
  if (identical(kind, "or")) {
    total <- 0.0
    for (i in seq_along(args)) {
      dens_i <- .expr_density_fast(args[[i]], t, prep, component)
      if (!is.finite(dens_i) || dens_i <= 0) next
      prod_surv <- 1.0
      for (j in seq_along(args)) {
        if (i == j) next
        Sj <- .expr_survival_fast(args[[j]], t, prep, component)
        prod_surv <- prod_surv * Sj
        if (!is.finite(prod_surv) || prod_surv == 0) break
      }
      if (!is.finite(prod_surv) || prod_surv == 0) next
      total <- total + dens_i * prod_surv
    }
    if (!is.finite(total)) total <- 0.0
    return(max(0.0, total))
  }
  0.0
}

.scenario_density_with_competitors <- function(expr, t, prep, component,
                                               competitor_exprs = list(),
                                               cache = NULL) {
  if (is.na(t) || !is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  cache <- .likelihood_cache_env(cache)
  simple_event <- function(e) {
    if (is.null(e) || is.null(e[['kind']]) || !identical(e[['kind']], "event")) return(FALSE)
    src <- e[['source']]
    !is.null(prep[["accumulators"]][[src]]) && is.null(prep[["pools"]][[src]])
  }
  expr_is_simple <- .expr_fastpath_supported(expr)
  if (expr_is_simple) {
    dens_fast <- .expr_density_fast(expr, t, prep, component)
    if (!is.finite(dens_fast) || dens_fast <= 0) {
      return(max(0.0, dens_fast))
    }
    if (length(competitor_exprs) > 0) {
      if (simple_event(expr) && all(vapply(competitor_exprs, simple_event, logical(1)))) {
        surv_fast <- .compute_survival_product(expr, competitor_exprs, prep, component, t,
                                               cache = cache)
      } else {
        surv_fast <- 1.0
        for (comp_expr in competitor_exprs) {
          if (.expr_fastpath_supported(comp_expr)) {
            Sv <- .expr_survival_fast(comp_expr, t, prep, component)
          } else {
            Sv <- .eval_expr_survival_cond(
              comp_expr, t, prep, component,
              forced_complete = integer(0),
              forced_survive = integer(0),
              cache = cache
            )
          }
          if (!is.finite(Sv) || Sv <= 0) {
            surv_fast <- 0.0
            break
          }
          surv_fast <- surv_fast * Sv
        }
      }
      dens_fast <- dens_fast * surv_fast
    }
    if (!is.finite(dens_fast) || dens_fast <= 0) return(0.0)
    return(dens_fast)
  }
  scenarios <- .expr_scenarios_at(
    expr, t, prep, component,
    forced_complete = integer(0),
    forced_survive = integer(0),
    cache = cache)
  if (length(scenarios) == 0) {
    return(0.0)
  }
  total <- 0.0
  for (sc in scenarios) {
    if (is.null(sc) || sc$weight <= 0) next
    weight <- sc$weight
    if (length(competitor_exprs) > 0) {
      surv_prod <- .compute_survival_product(
        outcome_expr = expr,
        competitor_exprs = competitor_exprs,
        prep = prep,
        component = component,
        t = t,
        cache = cache
      )
      weight <- weight * surv_prod
    }
    total <- total + weight
  }
  if (length(scenarios) > 0) attr(total, "scenarios") <- scenarios
  total
}

# ============================================================================== 
# Fast-path helpers for competitor survival products (ported from legacy API)
# ==============================================================================

.cluster_competitors <- function(exprs, prep) {
  n <- length(exprs)
  if (n == 0) return(list())
  src_sets <- lapply(exprs, .expr_sources, prep)
  visited <- rep(FALSE, n)
  clusters <- list()
  for (i in seq_len(n)) {
    if (visited[[i]]) next
    stack <- i
    visited[[i]] <- TRUE
    members <- integer(0)
    while (length(stack) > 0) {
      idx <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      members <- c(members, idx)
      for (j in seq_len(n)) {
        if (visited[[j]]) next
        if (length(intersect(src_sets[[idx]], src_sets[[j]])) > 0) {
          visited[[j]] <- TRUE
          stack <- c(stack, j)
        }
      }
    }
    clusters[[length(clusters) + 1L]] <- list(indices = members, exprs = exprs[members], sources = src_sets[members])
  }
  clusters
}

.get_event_source <- function(node) {
  if (is.null(node)) return(NULL)
  if (!is.null(node$kind) && identical(node$kind, "event")) {
    return(node$source %||% NULL)
  }
  NULL
}

.extract_guard_plus <- function(expr) {
  if (is.null(expr) || is.null(expr$kind)) return(NULL)
  if (identical(expr$kind, "guard")) {
    ref_id <- .get_event_source(expr$reference)
    blk_id <- .get_event_source(expr$blocker)
    if (is.null(ref_id) || is.null(blk_id)) return(NULL)
    return(list(guard = expr, ref = ref_id, block = blk_id, extras = character(0)))
  }
  if (identical(expr$kind, "and")) {
    args <- expr$args %||% list()
    guard_idx <- which(vapply(args, function(arg) !is.null(arg$kind) && identical(arg$kind, "guard"), logical(1)))
    if (length(guard_idx) != 1) return(NULL)
    guard <- args[[guard_idx[[1]]]]
    ref_id <- .get_event_source(guard$reference)
    blk_id <- .get_event_source(guard$blocker)
    if (is.null(ref_id) || is.null(blk_id)) return(NULL)
    extras <- args[-guard_idx[[1]]]
    extra_sources <- character(0)
    for (extra in extras) {
      src <- .get_event_source(extra)
      if (is.null(src)) return(NULL)
      extra_sources <- c(extra_sources, src)
    }
    return(list(guard = guard, ref = ref_id, block = blk_id, extras = unique(extra_sources)))
  }
  NULL
}

.event_density_fn <- function(source_id, prep, component) {
  event_expr <- list(kind = "event", source = source_id, k = NULL)
  function(t) vapply(
    t,
    function(tt) .eval_expr_likelihood(event_expr, tt, prep, component),
    numeric(1)
  )
}

.event_survival_fn <- function(source_id, prep, component) {
  event_expr <- list(kind = "event", source = source_id, k = NULL)
  function(t) vapply(
    t,
    function(tt) .eval_expr_survival(event_expr, tt, prep, component),
    numeric(1)
  )
}

.guard_pair_joint_survival <- function(t, prep, component, ref_a, block_a, extras) {
  if (t <= 0) return(1.0)
  f_a <- .event_density_fn(ref_a, prep, component)
  S_a <- .event_survival_fn(ref_a, prep, component)
  f_b <- .event_density_fn(block_a, prep, component)
  S_b <- .event_survival_fn(block_a, prep, component)
  S_extra_t <- 1.0
  if (length(extras) > 0) {
    S_extra_t <- prod(vapply(extras, function(src) .event_survival_fn(src, prep, component)(t), numeric(1)))
  }
  rel_tol <- .integrate_rel_tol()
  abs_tol <- .integrate_abs_tol()
  term_case2 <- if (is.finite(t)) {
    stats::integrate(function(v) f_a(v) * S_b(v), lower = t, upper = Inf,
                     rel.tol = rel_tol, abs.tol = abs_tol,
                     stop.on.error = FALSE)$value
  } else 0.0
  term_s1_gt_t <- stats::integrate(function(u) f_b(u) * S_a(u), lower = t, upper = Inf,
                                   rel.tol = rel_tol, abs.tol = abs_tol,
                                   stop.on.error = FALSE)$value
  if (!is.finite(term_s1_gt_t)) term_s1_gt_t <- 0.0
  term_s1_le_t <- stats::integrate(function(u) f_b(u) * S_a(u), lower = 0, upper = t,
                                   rel.tol = rel_tol, abs.tol = abs_tol,
                                   stop.on.error = FALSE)$value
  if (!is.finite(term_s1_le_t)) term_s1_le_t <- 0.0
  val <- term_case2 + term_s1_gt_t + S_extra_t * term_s1_le_t
  if (!is.finite(val)) val <- 0.0
  if (val < 0) val <- 0.0
  if (val > 1) val <- 1.0
  val
}

.compute_survival_product <- function(outcome_expr, competitor_exprs, prep, component, t,
                                      cache = NULL) {
  if (length(competitor_exprs) == 0) return(1.0)
  clusters <- .cluster_competitors(competitor_exprs, prep)
  prod_val <- 1.0
  for (cluster in clusters) {
    exprs <- cluster$exprs
    cluster_val <- NULL
    if (length(exprs) == 2) {
      info1 <- .extract_guard_plus(exprs[[1]])
      info2 <- .extract_guard_plus(exprs[[2]])
      if (!is.null(info1) && !is.null(info2)) {
        if (identical(info1$ref, info2$block) && identical(info1$block, info2$ref)) {
          guard_go1 <- info1
          guard_stop <- info2
        } else if (identical(info2$ref, info1$block) && identical(info2$block, info1$ref)) {
          guard_go1 <- info2
          guard_stop <- info1
        } else {
          guard_go1 <- NULL
          guard_stop <- NULL
        }
        if (!is.null(guard_go1) && length(guard_stop$extras) == length(unique(guard_stop$extras))) {
          cluster_val <- .guard_pair_joint_survival(t, prep, component,
                                                    ref_a = guard_go1$ref,
                                                    block_a = guard_go1$block,
                                                    extras = guard_stop$extras)
        }
      }
    }
    if (is.null(cluster_val)) {
      cluster_val <- prod(vapply(exprs, function(ce) {
        .eval_expr_survival_cond(
          ce, t, prep, component,
          forced_complete = integer(0),
          forced_survive = integer(0),
          cache = cache
        )
      }, numeric(1)))
    }
    prod_val <- prod_val * cluster_val
  }
  prod_val
}

# Get CDF (cumulative probability) for an expression
.eval_expr_cdf <- function(expr, t, prep, component) {
  kind <- expr[['kind']]
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(if (t >= 0) 1.0 else 0.0)
    }
    
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_cdf_at(prep, source_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)))
    }
    
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = instant completion
      }
      return(.event_cdf_at(prep, source_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0)))
    }
    
    return(1.0)
  }
  
  if (identical(kind, "and")) {
    # CDF of max: product of CDFs
    args <- expr[['args']]
    prod_cdf <- 1.0
    for (arg in args) {
      cdf_arg <- .eval_expr_cdf(arg, t, prep, component)
      prod_cdf <- prod_cdf * cdf_arg
    }
    return(prod_cdf)
  }
  
  if (identical(kind, "or")) {
    # CDF of min: 1 - product of survivals
    args <- expr[['args']]
    prod_surv <- 1.0
    for (arg in args) {
      surv_arg <- .eval_expr_survival(arg, t, prep, component)
      prod_surv <- prod_surv * surv_arg
    }
    return(1 - prod_surv)
  }
  
  if (identical(kind, "guard")) {
    # Guard: integrate density to obtain CDF
    if (t <= 0) return(0.0)
    if (!is.finite(t)) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(ui) {
        if (!is.finite(ui) || ui < 0) return(0.0)
        .eval_expr_likelihood(expr, ui, prep, component)
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(
        stats::integrate(dens_fun, lower = 0, upper = t,
                         rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                         stop.on.error = FALSE)$value,
        error = function(e) 0.0
      )
    }
    val <- compute_integral()
    cdf <- as.numeric(val)
    if (!is.finite(cdf)) cdf <- 0.0
    if (cdf < 0) cdf <- 0.0
    if (cdf > 1) cdf <- 1.0
    return(cdf)
  }
  
  return(NA_real_)
}

# Get survival function for an expression
.eval_expr_survival <- function(expr, t, prep, component) {
  kind <- expr[['kind']]
  
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(if (t < 0) 1.0 else 0.0)
    }
    
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_survival_at(prep, source_id, component, t,
                                forced_complete = integer(0),
                                forced_survive = integer(0)))
    }
    
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = always survive (never fires)
      }
      return(.event_survival_at(prep, source_id, component, t,
                                forced_complete = integer(0),
                                forced_survive = integer(0)))
    }
    
    return(1.0)
  }
  
  if (identical(kind, "and")) {
    # Survival of max: at least one must be > t
    # S_max(t) = 1 - F_max(t) = 1 - Π F_i(t)
    return(1 - .eval_expr_cdf(expr, t, prep, component))
  }
  
  if (identical(kind, "or")) {
    # Survival of min: all must be > t
    # S_min(t) = Π S_i(t)
    args <- expr[['args']]
    prod_surv <- 1.0
    for (arg in args) {
      surv_arg <- .eval_expr_survival(arg, t, prep, component)
      prod_surv <- prod_surv * surv_arg
    }
    return(prod_surv)
  }
  
  if (identical(kind, "guard")) {
    # Survival of guard = 1 - ∫_0^t f_guard(u) du
    if (!is.finite(t)) return(0.0)
    if (t <= 0) return(1.0)
    dens_fun <- function(u) {
      vapply(u, function(uu) {
        if (!is.finite(uu) || uu < 0) return(0.0)
        .eval_expr_likelihood(expr, uu, prep, component)
      }, numeric(1))
    }
    compute_integral <- function() {
      tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                stop.on.error = FALSE)[['value']],
               error = function(e) 0.0)
    }
    val <- compute_integral()
    s <- 1.0 - as.numeric(val)
    if (!is.finite(s)) s <- 0.0
    if (s < 0) s <- 0.0
    if (s > 1) s <- 1.0
    return(s)
  }
  
  return(0.0)
}

# Compute likelihood (density or probability) for an expression
.eval_expr_likelihood <- function(expr, t, prep, component,
                                  forced_complete = integer(0),
                                  forced_survive = integer(0),
                                  cache = NULL) {
  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)
  kind <- expr[['kind']]
  
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    
    # Handle special sources
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      # These don't contribute density directly
      return(0.0)
    }
    
    # Check if it's a pool
    if (!is.null(prep[["pools"]][[source_id]])) {
      return(.event_density_at(prep, source_id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive))
    }
    
    # Check if it's an accumulator
    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      # Check component membership
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(0.0)  # Not active in this component
      }
      return(.event_density_at(prep, source_id, component, t,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive))
    }
    
    return(0.0)
  }
  
  if (identical(kind, "and")) {
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    scenarios <- .expr_scenarios_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache)
    if (length(scenarios) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    total_density <- sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total_density, "scenarios") <- scenarios
    return(total_density)
  }
  
  if (identical(kind, "or")) {
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    scenarios <- .expr_scenarios_at(
      expr, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache)
    if (length(scenarios) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    total_density <- sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total_density, "scenarios") <- scenarios
    return(total_density)
  }
  
  if (identical(kind, "not")) {
    # NOT typically doesn't have density, only affects probability
    # If we reach here, it means NOT is used in an outcome expression
    # This is unusual but could mean "outcome never happens"
    return(0.0)
  }
  
  if (identical(kind, "guard")) {
    reference <- expr[["reference"]]
    blocker <- expr[["blocker"]]
    unless_list <- expr[["unless"]] %||% list()
    
    if (!is.finite(t) || t < 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    
    ref_scen <- .expr_scenarios_at(
      reference, t, prep, component,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      cache = cache)
    if (length(ref_scen) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    
    blocker_sources <- integer(0)
    if (!is.null(blocker)) {
      blocker_sources <- .expr_sources(blocker, prep)
      blocker_label <- .label_to_id(prep, blocker[["source"]])
      if (!is.na(blocker_label)) {
        blocker_sources <- .forced_union(prep, blocker_sources, blocker_label)
      }
    }
    
    compute_S_eff <- function(fc_ctx, fs_ctx) {
      if (is.null(blocker)) return(1.0)
      if (length(unless_list) == 0) {
        return(.eval_expr_survival_cond(
          blocker, t, prep, component,
          forced_complete = fc_ctx,
          forced_survive = fs_ctx,
          cache = cache
        ))
      }
      f_blocker <- function(tt) {
        .eval_expr_likelihood(
          blocker, tt, prep, component,
          forced_complete = fc_ctx,
          forced_survive = fs_ctx,
          cache = cache
        )
      }
      S_prot_prod <- function(tt) {
        if (length(unless_list) == 0) return(rep(1.0, length(tt)))
        vapply(tt, function(ui) {
          vals <- vapply(unless_list, function(unl) {
            .eval_expr_survival_cond(
              unl, ui, prep, component,
              forced_complete = fc_ctx,
              forced_survive = fs_ctx,
              cache = cache
            )
          }, numeric(1))
          prod(vals)
        }, numeric(1))
      }
      tryCatch({
        if (!is.finite(t) || t <= 0) return(1.0)
        val <- stats::integrate(
          function(u) {
            vapply(u, function(ui) f_blocker(ui) * S_prot_prod(ui), numeric(1))
          },
          lower = 0,
          upper = t,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]]
        s <- 1.0 - as.numeric(val)
        if (!is.finite(s)) 0.0 else max(0.0, min(1.0, s))
      }, error = function(e) 1.0)
    }
    
    scenarios <- list()
    for (sc in ref_scen) {
      if (is.null(sc)) next
      fc_ctx <- .coerce_forced_ids(prep, sc$forced_complete)
      fs_ctx <- .coerce_forced_ids(prep, sc$forced_survive)
      S_eff <- compute_S_eff(fc_ctx, fs_ctx)
      if (!is.finite(S_eff) || S_eff <= 0) next
      weight <- sc$weight * S_eff
      if (!is.finite(weight) || weight <= 0) next
      fcomp <- fc_ctx
      fsurv <- .forced_union(prep, fs_ctx, blocker_sources)
      scenarios[[length(scenarios) + 1L]] <- list(
        weight = as.numeric(weight),
        forced_complete = fcomp,
        forced_survive = .coerce_forced_ids(prep, fsurv)
      )
    }
    total <- sum(vapply(scenarios, `[[`, numeric(1), "weight"))
    attr(total, "scenarios") <- scenarios
    return(total)
  }
  
  stop(sprintf("Unsupported expression kind '%s'", kind))
}

# ==============================================================================
# PART 5: Outcome Likelihood
# ==============================================================================

# Integrand for scenario-conditioned probability (density already encodes competitors)
.integrand_outcome_density <- function(t, expr, prep, component, competitor_exprs,
                                       cache = NULL) {
  cache <- .likelihood_cache_env(cache)
  vapply(t, function(tt) {
    if (!is.finite(tt) || tt < 0) return(0.0)
    .scenario_density_with_competitors(expr, tt, prep, component,
                                       competitor_exprs = competitor_exprs,
                                       cache = cache)
  }, numeric(1))
}

# Scenario-conditioned probability of an expression over [0, upper_limit].
# competitor_exprs retained for backward compatibility but ignored.
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf,
                                           competitor_exprs = NULL,
                                           cache = NULL) {
  if (!is.finite(upper_limit)) upper_limit <- Inf
  cache <- .likelihood_cache_env(cache)
  integrand <- function(t) {
    val <- .integrand_outcome_density(
      t,
      expr = expr,
      prep = prep,
      component = component,
      competitor_exprs = competitor_exprs %||% list(),
      cache = cache
    )
    val
  }
  res <- tryCatch(
    stats::integrate(
      integrand,
      lower = 0,
      upper = upper_limit,
      rel.tol = .integrate_rel_tol(),
      abs.tol = .integrate_abs_tol(),
      stop.on.error = FALSE
    ),
    error = function(e) list(value = 0.0)
  )
  val <- as.numeric(res[["value"]])
  if (!is.finite(val)) val <- 0.0
  if (val < 0) val <- 0.0
  if (val > 1) val <- 1.0
  val
}

.expr_sources <- function(expr, prep) {
  gather_labels <- function(ex) {
    if (is.null(ex) || is.null(ex[['kind']])) return(character(0))
    kind <- ex[['kind']]
    if (identical(kind, "event")) {
      source_id <- ex[['source']]
      if (!is.null(prep[["pools"]][[source_id]])) {
        members <- prep[["pools"]][[source_id]][['members']] %||% character(0)
        unique(unlist(lapply(members, function(mid) {
          if (!is.null(prep[["accumulators"]][[mid]])) return(mid)
          if (!is.null(prep[["pools"]][[mid]])) {
            return(gather_labels(list(kind = "event", source = mid)))
          }
          character(0)
        })))
      } else {
        source_id
      }
    } else if (identical(kind, "and") || identical(kind, "or")) {
      unique(unlist(lapply(ex[['args']], gather_labels)))
    } else if (identical(kind, "guard")) {
      ref_src <- gather_labels(ex[['reference']])
      blk_src <- gather_labels(ex[['blocker']])
      unless_src <- unique(unlist(lapply(ex[['unless']] %||% list(), gather_labels)))
      unique(c(ref_src, blk_src, unless_src))
    } else if (identical(kind, "not")) {
      gather_labels(ex[['arg']])
    } else {
      character(0)
    }
  }
  .labels_to_ids(prep, gather_labels(expr))
}

.build_competitor_exprs <- function(prep, target_label, target_expr) {
  outcome_defs <- prep[["outcomes"]]
  comp_labels <- setdiff(names(outcome_defs), target_label)
  if (length(comp_labels) == 0) return(list())
  target_opts <- outcome_defs[[target_label]][['options']] %||% list()
  target_map <- target_opts[['map_outcome_to']] %||% target_label
  comp_labels <- Filter(function(lbl) {
    comp_opts <- outcome_defs[[lbl]][['options']] %||% list()
    if (!is.null(comp_opts[['alias_of']])) return(FALSE)
    expr_lbl <- outcome_defs[[lbl]][['expr']]
    if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
      src <- expr_lbl[['source']]
      if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
    }
    comp_target <- comp_opts[['map_outcome_to']] %||% lbl
    if (identical(as.character(comp_target), as.character(target_map))) return(FALSE)
    TRUE
  }, comp_labels)
  if (length(comp_labels) == 0) return(list())
  comps <- lapply(comp_labels, function(lbl) outcome_defs[[lbl]][['expr']])
  target_sources <- .expr_sources(target_expr, prep)
  if (length(comps) > 0) {
    comps <- lapply(comps, function(comp_expr) {
      if (!is.null(comp_expr[['kind']]) && identical(comp_expr[['kind']], "guard")) {
        blocker_sources <- .expr_sources(comp_expr[['blocker']], prep)
        if (length(target_sources) > 0 && all(target_sources %in% blocker_sources)) {
          return(comp_expr[['reference']])
        }
      }
      comp_expr
    })
  }
  if (!is.null(target_expr[['kind']]) && identical(target_expr[['kind']], "guard") && length(comps) > 0) {
    blocker_sources <- .expr_sources(target_expr[['blocker']], prep)
    keep_idx <- vapply(comps, function(comp_expr) {
      comp_sources <- .expr_sources(comp_expr, prep)
      length(comp_sources) == 0 || any(!comp_sources %in% blocker_sources)
    }, logical(1))
    comps <- comps[keep_idx]
    comp_labels <- comp_labels[keep_idx]
  }
  if (length(comps) > 0 && isTRUE(getOption("uuber.filter_guard_competitors", TRUE))) {
    primary_guards <- .collect_guards(target_expr)
    if (length(primary_guards) > 0) {
      guard_keep <- vapply(seq_along(comps), function(i) {
        comp_expr <- comps[[i]]
        !any(vapply(primary_guards, function(pg) .guards_conflict(pg, comp_expr), logical(1)))
      }, logical(1))
      comps <- comps[guard_keep]
      comp_labels <- comp_labels[guard_keep]
    }
  }
  comps
}

.prepare_competitor_map <- function(prep) {
  outcome_defs <- prep[["outcomes"]] %||% list()
  if (length(outcome_defs) == 0) return(list())
  comps <- lapply(names(outcome_defs), function(lbl) {
    .build_competitor_exprs(prep, lbl, outcome_defs[[lbl]][['expr']])
  })
  names(comps) <- names(outcome_defs)
  comps
}

# Detect shared-gate pattern (A & C) vs (B & C) across outcome definitions
.detect_shared_gate <- function(outcome_defs, target_label) {
  out_def <- outcome_defs[[target_label]]
  if (is.null(out_def) || is.null(out_def[['expr']])) return(NULL)
  expr <- out_def[['expr']]
  if (is.null(expr[['kind']]) || !identical(expr[['kind']], "and") || length(expr[['args']]) != 2) return(NULL)
  get_evt_id <- function(node) {
    if (!is.null(node) && !is.null(node[['kind']]) && identical(node[['kind']], "event")) return(node[['source']])
    NULL
  }
  x_id <- get_evt_id(expr[['args']][[1]])
  c_id <- get_evt_id(expr[['args']][[2]])
  if (is.null(x_id) || is.null(c_id)) {
    x_id <- get_evt_id(expr[['args']][[2]])
    c_id <- get_evt_id(expr[['args']][[1]])
  }
  if (is.null(x_id) || is.null(c_id)) return(NULL)
  # Find another outcome with same gate c_id
  for (lbl in names(outcome_defs)) {
    if (identical(lbl, target_label)) next
    other_def <- outcome_defs[[lbl]]
    if (!is.null(other_def[['options']][['alias_of']])) next
    other <- other_def[['expr']]
    if (is.null(other) || is.null(other[['kind']]) || !identical(other[['kind']], "and") || length(other[['args']]) != 2) next
    a1 <- get_evt_id(other[['args']][[1]])
    a2 <- get_evt_id(other[['args']][[2]])
    ids_other <- c(a1, a2)
    if (c_id %in% ids_other) {
      y_id <- setdiff(ids_other, c_id)
      if (length(y_id) == 1 && !identical(y_id[[1]], x_id)) {
        return(list(x_id = x_id, y_id = y_id[[1]], c_id = c_id))
      }
    }
  }
  NULL
}

# Find a single shared-gate pair among non-alias outcomes, if present
.find_shared_gate_pair <- function(outcome_defs) {
  labels <- names(outcome_defs)
  core_labels <- Filter(function(lbl) is.null(outcome_defs[[lbl]][['options']][['alias_of']]), labels)
  # collect all ANDs with two events
  events <- lapply(core_labels, function(lbl) {
    ex <- outcome_defs[[lbl]][['expr']]
    if (!is.null(ex[['kind']]) && identical(ex[['kind']], "and") && length(ex[['args']]) == 2) {
      a1 <- if (!is.null(ex[['args']][[1]][['kind']]) && ex[['args']][[1]][['kind']] == "event") ex[['args']][[1]][['source']] else NULL
      a2 <- if (!is.null(ex[['args']][[2]][['kind']]) && ex[['args']][[2]][['kind']] == "event") ex[['args']][[2]][['source']] else NULL
      if (!is.null(a1) && !is.null(a2)) return(list(label = lbl, ids = c(a1, a2)))
    }
    NULL
  })
  events <- events[!vapply(events, is.null, logical(1))]
  if (length(events) < 2) return(NULL)
  # try to find two with exactly one shared id (the gate)
  for (i in seq_along(events)) {
    for (j in seq_along(events)) {
      if (j <= i) next
      ids_i <- events[[i]][['ids']]
      ids_j <- events[[j]][['ids']]
      shared <- intersect(ids_i, ids_j)
      if (length(shared) == 1) {
        c_id <- shared[[1]]
        x_id <- setdiff(ids_i, c_id)[[1]]
        y_id <- setdiff(ids_j, c_id)[[1]]
        return(list(
          label_x = events[[i]][['label']],
          label_y = events[[j]][['label']],
          x_id = x_id,
          y_id = y_id,
          c_id = c_id
        ))
      }
    }
  }
  NULL
}

# Vectorized integrands for shared-gate pair; returns list(prob_x, prob_y)
.shared_gate_pair_probs <- function(prep, component, pair_info,
                                    rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol()) {
  get_f <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti < 0) return(0.0)
        if (!is.null(prep[["pools"]][[id]])) {
          return(.pool_density(
            prep, id, component, ti,
            forced_complete = integer(0),
            forced_survive = integer(0)
          ))
        }
        if (!is.null(prep[["accumulators"]][[id]])) {
          return(.event_density_at(
            prep, id, component, ti,
            forced_complete = integer(0),
            forced_survive = integer(0)
          ))
        }
        0.0
      }, numeric(1))
    }
  }
  
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        if (!is.null(prep[["pools"]][[id]])) {
          return(.pool_survival(
            prep, id, component, ti,
            forced_complete = integer(0),
            forced_survive = integer(0)
          ))
        }
        if (!is.null(prep[["accumulators"]][[id]])) {
          return(.event_survival_at(
            prep, id, component, ti,
            forced_complete = integer(0),
            forced_survive = integer(0)
          ))
        }
        0.0
      }, numeric(1))
    }
  }
  
  get_F <- function(id) {
    Sfun <- get_S(id)
    function(tt) 1.0 - Sfun(tt)
  }
  
  fX <- get_f(pair_info[['x_id']])
  fY <- get_f(pair_info[['y_id']])
  fC <- get_f(pair_info[['c_id']])
  FX <- get_F(pair_info[['x_id']])
  FY <- get_F(pair_info[['y_id']])
  FC <- get_F(pair_info[['c_id']])
  SX <- get_S(pair_info[['x_id']])
  SY <- get_S(pair_info[['y_id']])
  
  palloc_x_single <- function(tt) {
    if (!is.finite(tt) || tt <= 0) return(0.0)
    denom <- FX(tt) * FY(tt)
    if (!is.finite(denom) || denom <= 0) return(0.0)
    integral_val <- tryCatch(
      stats::integrate(
        function(u) fX(u) * FY(u),
        lower = 0,
        upper = tt,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        stop.on.error = FALSE
      )[["value"]],
      error = function(e) 0.0
    )
    out <- 1.0 - (as.numeric(integral_val) / as.numeric(denom))
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
  
  palloc_x <- function(t) {
    vapply(t, palloc_x_single, numeric(1))
  }
  
  integrand_x <- function(t) {
    vapply(t, function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      term1 <- fX(tt) * FC(tt) * SY(tt)
      termCother <- fC(tt) * FX(tt) * SY(tt)
      den <- FX(tt) * FY(tt)
      tie <- if (is.finite(den) && den > 0) {
        fC(tt) * FX(tt) * FY(tt) * palloc_x(tt)
      } else {
        0.0
      }
      term1 + termCother + tie
    }, numeric(1))
  }
  
  integrand_y <- function(t) {
    vapply(t, function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      term1 <- fY(tt) * FC(tt) * SX(tt)
      termCother <- fC(tt) * FY(tt) * SX(tt)
      den <- FX(tt) * FY(tt)
      tie <- if (is.finite(den) && den > 0) {
        fC(tt) * FX(tt) * FY(tt) * (1.0 - palloc_x(tt))
      } else {
        0.0
      }
      term1 + termCother + tie
    }, numeric(1))
  }
  
  px <- tryCatch(
    stats::integrate(
      integrand_x,
      lower = 0,
      upper = Inf,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      stop.on.error = FALSE
    )[["value"]],
    error = function(e) 0.0
  )
  py <- tryCatch(
    stats::integrate(
      integrand_y,
      lower = 0,
      upper = Inf,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      stop.on.error = FALSE
    )[["value"]],
    error = function(e) 0.0
  )
  
  px <- as.numeric(px)
  py <- as.numeric(py)
  if (!is.finite(px) || px < 0) px <- 0.0
  if (!is.finite(py) || py < 0) py <- 0.0
  c(px, py)
}

# Compute likelihood for a single trial/outcome
.outcome_likelihood <- function(outcome_label, rt, prep, component) {
  outcome_defs <- prep[["outcomes"]]
  cache <- .likelihood_cache_env()
  
  # FIRST: Check for GUESS outcomes (can be both defined and from guess policy)
  # GUESS outcomes arise from component-level guess policies
  if (identical(outcome_label, "GUESS")) {
    # Get the guess policy for this component
    guess_policy <- .get_component_attr(prep, component, "guess")
    
    if (!is.null(guess_policy) && !is.null(guess_policy[['weights']])) {
      # The guess policy has:
      # - outcome: the label for guess outcomes
      # - weights: named vector of keep probabilities
      
      guess_weights <- guess_policy[['weights']]
      
      # Probability of GUESS outcome is:
      # Sum over outcomes: P(outcome fires) * P(guess | outcome)
      # P(guess | outcome) = 1 - weight[outcome]
      
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep[["default_deadline"]]
      
      if (is.finite(deadline)) {
        if (is.na(rt)) {
          total_prob <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            expr <- outcome_defs[[label]][['expr']]
            comp_exprs_guess <- prep[['.competitors']][[label]] %||% list()
            prob_outcome <- .integrate_outcome_probability(expr, prep, component, deadline,
                                                           competitor_exprs = comp_exprs_guess,
                                                           cache = cache)
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_prob <- total_prob + prob_outcome * guess_prob
          }
          return(total_prob)
        } else {
          total_dens <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            expr <- outcome_defs[[label]][['expr']]
            comp_exprs_guess <- prep[['.competitors']][[label]] %||% list()
            dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                         competitor_exprs = comp_exprs_guess,
                                                         cache = cache)
            if (dens_r == 0) next
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_dens <- total_dens + guess_prob * dens_r
          }
          return(total_dens)
        }
      }
    }
    # If no guess policy or deadline, fall through to regular handling
  }
  
  # Find the outcome
  outcome_def <- outcome_defs[[outcome_label]]
  if (is.null(outcome_def)) {
    # Check if this is a special outcome
    if (is.na(outcome_label) || identical(outcome_label, "NA")) {
      # NA outcome - this could be from map_outcome_to = NA
      # Need to find which outcome maps to NA
      for (label in names(outcome_defs)) {
        out_def <- outcome_defs[[label]]
        if (!is.null(out_def[['options']][['map_outcome_to']])) {
          mapped <- out_def[['options']][['map_outcome_to']]
          if (is.na(mapped) || identical(mapped, "NA")) {
            # Found the outcome that maps to NA
            expr <- out_def[['expr']]
            if (is.na(rt)) {
              # Integrate race probability for the mapped source outcome
              deadline <- .get_component_attr(prep, component, "deadline")
              deadline <- deadline %||% prep[["default_deadline"]]
              comp_labels_map <- setdiff(names(outcome_defs), label)
              if (length(comp_labels_map) > 0) {
                comp_labels_map <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels_map)
              }
              comp_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              return(.integrate_outcome_probability(expr, prep, component, deadline,
                                                    competitor_exprs = comp_exprs_map,
                                                    cache = cache))
            } else {
              # Race density at rt for the mapped source outcome
              comp_labels_map <- setdiff(names(outcome_defs), label)
              if (length(comp_labels_map) > 0) {
                comp_labels_map <- Filter(function(lbl) {
                  expr_lbl <- outcome_defs[[lbl]][['expr']]
                  if (!is.null(expr_lbl[['kind']]) && identical(expr_lbl[['kind']], "event")) {
                    src <- expr_lbl[['source']]
                    if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                  }
                  TRUE
                }, comp_labels_map)
              }
              comp_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]][['expr']]) else list()
              dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                           competitor_exprs = comp_exprs_map,
                                                           cache = cache)
              return(as.numeric(dens_r))
            }
          }
        }
      }
      # Couldn't find the NA mapping
      warning(sprintf("Unknown outcome '%s' and no NA mapping found", outcome_label))
      return(0.0)
    }
    
    # Check for special deadline outcomes
    if (identical(outcome_label, "NR_DEADLINE")) {
      # Deadline reached - compute probability that no outcome occurred by deadline
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep[["default_deadline"]]
      
      if (!is.finite(deadline)) return(0.0)
      
      # P(no outcome by deadline) = Π_i S_i(deadline)
      # where i ranges over all possible outcomes
      prob_none <- 1.0
      for (label in names(outcome_defs)) {
        out_def <- outcome_defs[[label]]
        expr <- out_def[['expr']]
        # Survival of this outcome at deadline
        surv <- .eval_expr_survival(expr, deadline, prep, component)
        prob_none <- prob_none * surv
      }
      return(prob_none)
    }
    
    warning(sprintf("Unknown outcome '%s'", outcome_label))
    return(0.0)
  }
  
  expr <- outcome_def[['expr']]
  options <- outcome_def[['options']] %||% list()
  
  # Check if outcome is mapped to something else
  if (!is.null(options[['map_outcome_to']])) {
    # This outcome is mapped, but we shouldn't reach here if data has the mapped label
    # Just proceed with the expression
  }
  
  # If this outcome has an outcome-level guess policy, it never appears as-is.
  # Its mass is redistributed to target labels via donors handled below.
  if (!is.null(options[['guess']])) {
    return(0.0)
  }
  
  competitor_exprs <- prep[['.competitors']][[outcome_label]] %||% list()
  
  # Build donor (guess) contributions: outcomes that guess into this label
  donors <- list()
  for (dlbl in setdiff(names(outcome_defs), outcome_label)) {
    dopt <- outcome_defs[[dlbl]][['options']][['guess']] %||% NULL
    if (is.null(dopt)) next
    glabels <- dopt[['labels']] %||% NULL
    gweights <- dopt[['weights']] %||% NULL
    if (is.null(glabels) || is.null(gweights)) next
    idx <- which(as.character(glabels) == outcome_label)
    if (length(idx) == 0) next
    w <- as.numeric(gweights[[idx[[1]]]])
    rt_policy <- dopt[['rt_policy']] %||% "keep"
    donors[[length(donors) + 1L]] <- list(
      expr = outcome_defs[[dlbl]][['expr']],
      weight = w,
      rt_policy = rt_policy,
      competitors = prep[['.competitors']][[dlbl]] %||% list()
    )
  }
  
  # Detect shared-gate AND pattern: (X & C) vs (Y & C)
  shared_gate_info <- .detect_shared_gate(outcome_defs, outcome_label)
  
  # Helper closures to get density/CDF for pools or accumulators
  get_f <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti < 0) return(0.0)
        .event_density_at(prep, id, component, ti,
                          forced_complete = integer(0),
                          forced_survive = integer(0))
      }, numeric(1))
    }
  }
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        .event_survival_at(prep, id, component, ti,
                           forced_complete = integer(0),
                           forced_survive = integer(0))
      }, numeric(1))
    }
  }
  get_F <- function(id) {
    Sfun <- get_S(id)
    function(tt) { 1.0 - Sfun(tt) }
  }
  
  shared_palloc <- NULL
  if (!is.null(shared_gate_info)) {
    fX <- get_f(shared_gate_info[['x_id']])
    fY <- get_f(shared_gate_info[['y_id']])
    fC <- get_f(shared_gate_info[['c_id']])
    FX <- get_F(shared_gate_info[['x_id']])
    FY <- get_F(shared_gate_info[['y_id']])
    FC <- get_F(shared_gate_info[['c_id']])
    SY <- get_S(shared_gate_info[['y_id']])
    shared_palloc <- function(limit) {
      if (!is.finite(limit) || limit <= 0) return(0.0)
      denom <- FX(limit) * FY(limit)
      if (!is.finite(denom) || denom <= 0) return(0.0)
      integral_val <- tryCatch(
        stats::integrate(
          function(u) fX(u) * FY(u),
          lower = 0,
          upper = limit,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]],
        error = function(e) 0.0
      )
      out <- 1.0 - (as.numeric(integral_val) / as.numeric(denom))
      if (!is.finite(out)) out <- 0.0
      max(0.0, min(1.0, out))
    }
  }
  
  # Handle NA/infinite RT by integration of race density
  if (is.na(rt) || !is.finite(rt) || rt < 0) {
    deadline <- .get_component_attr(prep, component, "deadline")
    deadline <- deadline %||% prep[["default_deadline"]]
    if (!is.null(shared_gate_info)) {
      integrand <- function(t) {
        vapply(t, function(tt) {
          term1 <- fX(tt) * FC(tt) * SY(tt)
          denom <- FX(tt) * FY(tt)
          termCother <- fC(tt) * FX(tt) * SY(tt)
          term2 <- if (!is.null(shared_palloc) && is.finite(denom) && denom > 0) {
            fC(tt) * FX(tt) * FY(tt) * shared_palloc(tt)
          } else {
            0.0
          }
          term1 + termCother + term2
        }, numeric(1))
      }
      keep_prob <- 1.0
      gp_comp <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp) && !is.null(gp_comp[['weights']])) {
        kp <- gp_comp[['weights']][[outcome_label]] %||% gp_comp[['weights']][[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob <- as.numeric(kp)
      }
      upper_lim <- if (is.finite(deadline)) deadline else Inf
      res <- tryCatch(
        stats::integrate(
          integrand,
          lower = 0,
          upper = upper_lim,
          rel.tol = .integrate_rel_tol(),
          abs.tol = .integrate_abs_tol(),
          stop.on.error = FALSE
        )[["value"]],
        error = function(e) 0.0
      )
      res <- as.numeric(res)
      if (!is.finite(res) || res < 0) res <- 0.0
      return(keep_prob * res)
    }
    integrand <- function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        base <- .scenario_density_with_competitors(expr, tt, prep, component,
                                                   competitor_exprs = competitor_exprs,
                                                   cache = cache)
        if (base == 0) return(0.0)
        add <- 0.0
        if (length(donors) > 0) {
          for (d in donors) {
            dens_d <- .scenario_density_with_competitors(d[['expr']], tt, prep, component,
                                                         competitor_exprs = d[['competitors']] %||% list(),
                                                         cache = cache)
            if (dens_d == 0) next
            add <- add + d[['weight']] * dens_d
          }
        }
        base + add
      }, numeric(1))
    }
    upper_lim <- if (is.finite(deadline)) deadline else Inf
    res <- tryCatch(
      stats::integrate(
        integrand,
        lower = 0,
        upper = upper_lim,
        rel.tol = .integrate_rel_tol(),
        abs.tol = .integrate_abs_tol(),
        stop.on.error = FALSE
      )[["value"]],
      error = function(e) 0.0
    )
    res <- as.numeric(res)
    if (!is.finite(res) || res < 0) res <- 0.0
    return(res)
  }
  
  # Race density at observed RT
  if (!is.null(shared_gate_info)) {
    term1 <- fX(rt) * FC(rt) * SY(rt)
    denom <- FX(rt) * FY(rt)
    term2 <- 0.0
    termCother <- fC(rt) * FX(rt) * SY(rt)
    if (!is.null(shared_palloc) && is.finite(denom) && denom > 0) {
      palloc <- shared_palloc(rt)
      term2 <- fC(rt) * FX(rt) * FY(rt) * palloc
    }
    base_val <- term1 + termCother + term2
  } else {
    dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                 competitor_exprs = competitor_exprs,
                                                 cache = cache)
    if (dens_r == 0) return(0.0)
    base_val <- dens_r
    # Include donor mass (e.g., timeout guesses) that keeps the observed RT
    if (length(donors) > 0) {
      donor_add <- 0.0
      for (d in donors) {
        if (identical(d[['rt_policy']], "na")) next
        dens_d <- .scenario_density_with_competitors(d[['expr']], rt, prep, component,
                                                     competitor_exprs = d[['competitors']] %||% list(),
                                                     cache = cache)
        if (dens_d == 0) next
        donor_add <- donor_add + d[['weight']] * dens_d
      }
      base_val <- base_val + donor_add
    }
  }
  # Apply component-level keep probability to base density for observed labels
  keep_prob_rt <- 1.0
  gp_comp_rt <- .get_component_attr(prep, component, "guess")
  if (!is.null(gp_comp_rt) && !is.null(gp_comp_rt[['weights']])) {
    kp <- gp_comp_rt[['weights']][[outcome_label]] %||% gp_comp_rt[['weights']][[normalize_label(outcome_label)]] %||% NULL
    if (!is.null(kp)) keep_prob_rt <- as.numeric(kp)
  }
  base_val <- base_val * keep_prob_rt
  # Component-level guess mass is handled via the explicit GUESS outcome
  base_val
}

# ==============================================================================
# PART 6: Top-level Likelihood Function
# ==============================================================================

#' Compute log-likelihood for data given a model
#'
#' @param model Model specification from new_API.R (race_spec or race_model_spec)
#' @param data Data frame with columns: outcome, rt, (optionally: component)
#' @return Log-likelihood value with per-trial attributes
compute_loglik <- function(model, data) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  
  if (!"outcome" %in% names(data) || !"rt" %in% names(data)) {
    stop("data must contain 'outcome' and 'rt' columns")
  }
  
  prep <- .prepare_model_for_likelihood(model)
  
  comp_info <- prep[["components"]]
  comp_ids <- comp_info[['ids']] %||% list()
  weights <- comp_info[['weights']] %||% numeric(0)
  has_weight_param <- isTRUE(any(comp_info[['has_weight_param']] %||% FALSE))
  is_mixture <- length(comp_ids) > 1
  default_component <- if (length(comp_ids) > 0) comp_ids[[1]] else "__default__"
  
  n_rows <- nrow(data)
  if (n_rows == 0) {
    total_ll <- 0.0
    attr(total_ll, "contributions") <- numeric(0)
    attr(total_ll, "log_contributions") <- numeric(0)
    attr(total_ll, "keys") <- character(0)
    return(total_ll)
  }
  
  has_component_col <- "component" %in% names(data)
  raw_keys <- character(n_rows)
  lik_values <- numeric(n_rows)
  cache_env <- new.env(parent = emptyenv())
  
  for (i in seq_len(n_rows)) {
    outcome <- as.character(data[['outcome']][[i]])
    rt_val <- as.numeric(data[['rt']][[i]])
    rt_key <- .time_key(rt_val)
    component_key <- default_component
    if (is_mixture) {
      if (has_weight_param || !has_component_col) {
        component_key <- "__mixture__"
      } else {
        raw_component <- data[['component']][[i]]
        if (is.na(raw_component)) {
          component_key <- default_component
        } else {
          component_chr <- as.character(raw_component)
          component_key <- if (nzchar(component_chr)) component_chr else default_component
        }
      }
    }
    key <- paste(component_key, outcome, rt_key, sep = "|")
    
    cached_val <- cache_env[[key]]
    if (!is.null(cached_val)) {
      lik <- cached_val
    } else {
      if (is_mixture) {
        if (identical(component_key, "__mixture__")) {
          base_weights <- if (length(comp_ids) > 0 && length(weights) == length(comp_ids)) {
            weights
          } else if (length(comp_ids) > 0) {
            rep(1 / length(comp_ids), length(comp_ids))
          } else {
            numeric(0)
          }
          lik <- 0.0
          if (length(comp_ids) > 0) {
            for (j in seq_along(comp_ids)) {
              lik <- lik + base_weights[[j]] * .outcome_likelihood(outcome, rt_val, prep, comp_ids[[j]])
            }
          }
        } else {
          lik <- .outcome_likelihood(outcome, rt_val, prep, component_key)
        }
      } else {
        lik <- .outcome_likelihood(outcome, rt_val, prep, component_key)
      }
      
      cache_env[[key]] <- lik
    }
    
    raw_keys[[i]] <- key
    lik_values[[i]] <- lik
  }
  log_lik_values <- ifelse(is.finite(lik_values) & lik_values > 0, log(lik_values), -Inf)
  
  total_ll <- sum(log_lik_values)
  attr(total_ll, "contributions") <- lik_values
  attr(total_ll, "log_contributions") <- log_lik_values
  attr(total_ll, "keys") <- raw_keys
  
  total_ll
}

compute_loglik_from_tables <- function(model_tables, data) {
  compute_loglik(model_tables, data)
}

# ==============================================================================
# PART 7b: Response probability helpers (analytic, no simulation)
# ==============================================================================

#' Compute probability of an outcome label under the model
#' If the model has multiple components and no specific component is supplied,
#' returns the mixture-weighted probability.
#'
#' @param model Race model spec (race_spec or race_model_spec)
#' @param outcome_label Outcome label to evaluate (character)
#' @param component Optional component id (character) for mixture models
#' @return Probability (numeric in [0,1])
response_probability <- function(model, outcome_label, component = NULL) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  comp_ids <- prep[["components"]][['ids']]
  weights <- prep[["components"]][['weights']]
  
  # Helper for a single component
  .prob_for_component <- function(comp_id) {
    # Shared-gate pair fast-path (can be disabled globally)
    use_fastpath <- getOption("uuber.shared_gate_fastpath", default = TRUE)
    if (isTRUE(use_fastpath)) {
      pair <- .find_shared_gate_pair(prep[["outcomes"]])
      if (!is.null(pair) && outcome_label %in% c(pair[['label_x']], pair[['label_y']])) {
        vals <- .shared_gate_pair_probs(prep, comp_id, pair)
        if (identical(outcome_label, pair[['label_x']])) return(vals[[1]]) else return(vals[[2]])
      }
    }
    # Alias-of outcomes: sum of referenced labels
    out_def <- prep[["outcomes"]][[outcome_label]]
    if (!is.null(out_def[['options']][['alias_of']])) {
      refs <- out_def[['options']][['alias_of']]
      refs <- as.character(refs)
      vals <- vapply(refs, function(lbl) as.numeric(.outcome_likelihood(lbl, NA_real_, prep, comp_id)), numeric(1))
      return(sum(vals))
    }
    base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, comp_id))
    # Apply component-level keep for non-special labels (complement goes to GUESS)
    if (!identical(outcome_label, "GUESS")) {
      gp <- .get_component_attr(prep, comp_id, "guess")
      if (!is.null(gp) && !is.null(gp[['weights']])) {
        keep <- gp[['weights']][[outcome_label]] %||% gp[['weights']][[normalize_label(outcome_label)]] %||% 1.0
        base <- base * as.numeric(keep)
      }
    }
    base
  }
  
  if (!is.null(component)) {
    return(.prob_for_component(component))
  }
  
  if (length(comp_ids) <= 1L) {
    comp <- if (length(comp_ids) == 1L) comp_ids[[1]] else "__default__"
    return(.prob_for_component(comp))
  }
  
  # Mixture: average with weights
  vals <- vapply(seq_along(comp_ids), function(i) .prob_for_component(comp_ids[[i]]), numeric(1))
  sum(weights * vals)
}

#' Compute probabilities for all defined outcomes in the model
#' If component is NULL, returns mixture-weighted probabilities.
#'
#' @param model Race model spec
#' @param component Optional component id
#' @return Named numeric vector of probabilities
response_probabilities <- function(model, component = NULL) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  labels <- names(prep[["outcomes"]])
  vals <- vapply(labels, function(lbl) response_probability(model, lbl, component = component), numeric(1))
  stats::setNames(vals, labels)
}

#' observed_response_probabilities: Aggregate to observable outcomes
#'
#' For models where some outcomes are mapped (e.g., map_outcome_to = NA),
#' aggregate those mapped probabilities into an NA bucket (if include_na=TRUE)
#' and return only observable labels.
#'
#' @param model Race model spec
#' @param component Optional component id
#' @param include_na If TRUE, include an 'NA' entry with summed probability of
#'                   all outcomes mapped to NA. If FALSE, exclude mapped
#'                   outcomes entirely.
#' @return Named numeric vector of probabilities for observable outcomes (and NA)
observed_response_probabilities <- function(model, component = NULL, include_na = TRUE) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  prep <- .prepare_model_for_likelihood(model)
  labels <- names(prep[["outcomes"]])
  # Exclude alias-of outcomes from observed set
  labels <- Filter(function(lbl) {
    is.null(prep[["outcomes"]][[lbl]][['options']][['alias_of']])
  }, labels)
  # Compute base probabilities for all labels
  base <- vapply(labels, function(lbl) response_probability(model, lbl, component = component), numeric(1))
  names(base) <- labels
  # Partition into observable vs mapped-to-NA
  obs <- numeric(0)
  na_sum <- 0.0
  for (lbl in labels) {
    out_def <- prep[["outcomes"]][[lbl]]
    map_to <- out_def[['options']][['map_outcome_to']] %||% NULL
    if (is.null(map_to)) {
      # observable label
      obs[lbl] <- base[[lbl]]
    } else if (is.na(map_to)) {
      # mapped to NA
      na_sum <- na_sum + base[[lbl]]
    } else {
      # mapped to another label: aggregate into that label
      obs_lbl <- as.character(map_to)
      obs[obs_lbl] <- (obs[obs_lbl] %||% 0.0) + base[[lbl]]
    }
  }
  if (include_na && na_sum > 0) {
    obs <- c(obs, "NA" = na_sum)
  }
  if (include_na) {
    total <- sum(obs)
    resid <- 1.0 - total
    if (!is.finite(resid)) resid <- 0.0
    if (resid > .Machine$double.eps) {
      if ("NA" %in% names(obs)) {
        obs[["NA"]] <- obs[["NA"]] + resid
      } else {
        obs <- c(obs, "NA" = resid)
      }
    }
  }
  obs
}
