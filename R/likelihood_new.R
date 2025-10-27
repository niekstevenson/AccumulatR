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
  prep[[".competitors"]] <- .prepare_competitor_map(prep)
  if (length(all_ids) > 0) {
    idx_map <- prep[[".id_index"]]
    rev_labels <- rep(NA_character_, length(idx_map))
    rev_labels[as.integer(idx_map)] <- names(idx_map)
    prep[[".id_labels"]] <- rev_labels
  } else {
    prep[[".id_labels"]] <- character(0)
  }
  prep[[".label_cache"]] <- new.env(parent = emptyenv(), hash = TRUE)
  outcomes <- prep[["outcomes"]] %||% list()
  if (length(outcomes) > 0) {
    for (lbl in names(outcomes)) {
      outcomes[[lbl]][['expr']] <- outcomes[[lbl]][['expr']]
    }
    prep[["outcomes"]] <- outcomes
  }
  comp_map <- prep[[".competitors"]] %||% list()
  if (length(comp_map) > 0) {
    for (lbl in names(comp_map)) {
      exprs <- comp_map[[lbl]] %||% list()
      if (length(exprs) > 0) {
        comp_map[[lbl]] <- exprs
      }
    }
    prep[[".competitors"]] <- comp_map
  }
  prep[[".expr_registry"]] <- .build_expr_registry(prep)
  # Attach stable numeric IDs directly onto expr objects to avoid signature drift
  .attach_expr_ids <- function(prep, registry) {
    outcomes <- prep[["outcomes"]] %||% list()
    if (length(outcomes) > 0) {
      for (lbl in names(outcomes)) {
        ex <- outcomes[[lbl]][['expr']]
        if (!is.null(ex)) {
          id <- .expr_registry_lookup_id(registry, ex)
          if (!is.na(id)) {
            ex[["._id"]] <- as.integer(id)
            outcomes[[lbl]][['expr']] <- ex
          }
        }
      }
      prep[["outcomes"]] <- outcomes
    }
    comp_map <- prep[[".competitors"]] %||% list()
    if (length(comp_map) > 0) {
      for (lbl in names(comp_map)) {
        exprs <- comp_map[[lbl]] %||% list()
        if (length(exprs) == 0) next
        for (i in seq_along(exprs)) {
          ex <- exprs[[i]]
          if (is.null(ex)) next
          id <- .expr_registry_lookup_id(registry, ex)
          if (!is.na(id)) {
            ex[["._id"]] <- as.integer(id)
            exprs[[i]] <- ex
          }
        }
        comp_map[[lbl]] <- exprs
      }
      prep[[".competitors"]] <- comp_map
    }
    prep
  }
  prep <- .attach_expr_ids(prep, prep[[".expr_registry"]])
  # Precompute outcome expression IDs for O(1) lookup
  reg <- prep[[".expr_registry"]]
  outcomes <- prep[["outcomes"]] %||% list()
  if (!is.null(reg) && length(outcomes) > 0) {
    expr_ids <- vapply(outcomes, function(od) {
      .expr_registry_lookup_id(reg, od[["expr"]])
    }, integer(1))
    prep[[".outcome_expr_ids"]] <- stats::setNames(as.integer(expr_ids), names(outcomes))
  } else {
    prep[[".outcome_expr_ids"]] <- integer(0)
  }
  prep
}

.label_to_id <- function(prep, label) {
  if (is.null(label) || length(label) == 0L) return(NA_integer_)
  idx_map <- prep[[".id_index"]]
  if (is.null(idx_map)) return(NA_integer_)
  val <- idx_map[[as.character(label)]]
  if (is.null(val) || is.na(val)) return(NA_integer_)
  as.integer(val)
}

.expr_signature <- function(expr) {
  if (is.null(expr)) return("NULL")
  rawToChar(serialize(expr, NULL, ascii = TRUE))
}

.expr_registry_lookup_id <- function(registry, expr) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(NA_integer_)
  # Prefer attached numeric id when present
  if (!is.null(expr[["._id"]])) {
    id <- as.integer(expr[["._id"]])
    if (is.finite(id) && id > 0) return(id)
  }
  idx_env <- registry[["index"]]
  if (is.null(idx_env)) return(NA_integer_)
  sig <- .expr_signature(expr)
  val <- idx_env[[sig]]
  if (is.null(val)) return(NA_integer_)
  as.integer(val)
}

.build_expr_registry <- function(prep) {
  nodes <- list()
  index <- new.env(parent = emptyenv(), hash = TRUE)

  register_expr <- function(expr) {
    if (is.null(expr) || is.null(expr[["kind"]])) return(NA_integer_)
    sig <- .expr_signature(expr)
    existing <- index[[sig]]
    if (!is.null(existing)) return(as.integer(existing))

    node_id <- length(nodes) + 1L
    kind <- expr[["kind"]] %||% NA_character_
    # Reserve ID and index EARLY to avoid child registration colliding with this ID
    # Place a minimal placeholder node which we will fill in after registering children
    placeholder <- list(
      id = as.integer(node_id),
      kind = kind,
      expr = expr,
      child_ids = integer(0),
      sources = .expr_sources(expr, prep),
      guard_ids = integer(0),
      finish_ids = integer(0)
    )
    nodes[[node_id]] <<- placeholder
    index[[sig]] <<- as.integer(node_id)

    node <- placeholder

    if (identical(kind, "and") || identical(kind, "or")) {
      args <- expr[["args"]] %||% list()
      if (length(args) > 0) {
        child_ids <- vapply(args, register_expr, integer(1))
        node[["child_ids"]] <- as.integer(child_ids)
        node[["args"]] <- args
        child_complete <- lapply(child_ids, function(cid) nodes[[cid]][["complete_ids"]] %||% integer(0))
        child_guard <- lapply(child_ids, function(cid) nodes[[cid]][["guard_ids"]] %||% integer(0))
        child_finish <- lapply(child_ids, function(cid) nodes[[cid]][["finish_ids"]] %||% integer(0))
        comb <- unique(unlist(child_complete, use.names = FALSE))
        guard_comb <- unique(unlist(child_guard, use.names = FALSE))
        finish_comb <- unique(unlist(child_finish, use.names = FALSE))
        node[["complete_ids"]] <- if (length(comb) > 0) as.integer(sort(comb)) else integer(0)
        node[["guard_ids"]] <- if (length(guard_comb) > 0) as.integer(sort(guard_comb)) else integer(0)
        node[["finish_ids"]] <- if (length(finish_comb) > 0) as.integer(sort(finish_comb)) else integer(0)
      } else {
        node[["child_ids"]] <- integer(0)
        node[["args"]] <- list()
        node[["complete_ids"]] <- integer(0)
        node[["guard_ids"]] <- integer(0)
        node[["finish_ids"]] <- integer(0)
      }
    } else if (identical(kind, "guard")) {
      ref_id <- register_expr(expr[["reference"]])
      blocker_expr <- expr[["blocker"]] %||% NULL
      blocker_id <- if (!is.null(blocker_expr)) register_expr(blocker_expr) else NA_integer_
      unless_list <- expr[["unless"]] %||% list()
      unless_ids <- if (length(unless_list) > 0) vapply(unless_list, register_expr, integer(1)) else integer(0)
      node[["reference_id"]] <- if (is.na(ref_id)) NA_integer_ else as.integer(ref_id)
      node[["blocker_id"]] <- if (is.na(blocker_id)) NA_integer_ else as.integer(blocker_id)
      node[["unless_ids"]] <- as.integer(unless_ids)
      child_ids <- c(ref_id, blocker_id, unless_ids)
      child_ids <- child_ids[!is.na(child_ids)]
      node[["child_ids"]] <- if (length(child_ids) > 0) as.integer(child_ids) else integer(0)
      node[["args"]] <- list()
      ref_complete <- if (!is.na(ref_id)) nodes[[ref_id]][["complete_ids"]] %||% integer(0) else integer(0)
      ref_guard <- if (!is.na(ref_id)) nodes[[ref_id]][["guard_ids"]] %||% integer(0) else integer(0)
      ref_finish <- if (!is.na(ref_id)) nodes[[ref_id]][["finish_ids"]] %||% integer(0) else integer(0)
      blocker_complete <- if (!is.na(blocker_id)) nodes[[blocker_id]][["complete_ids"]] %||% integer(0) else integer(0)
      blocker_guard <- if (!is.na(blocker_id)) nodes[[blocker_id]][["guard_ids"]] %||% integer(0) else integer(0)
      blocker_finish <- if (!is.na(blocker_id)) nodes[[blocker_id]][["finish_ids"]] %||% integer(0) else integer(0)
      unless_complete <- if (length(unless_ids) > 0) unique(unlist(lapply(unless_ids, function(uid) nodes[[uid]][["complete_ids"]] %||% integer(0)), use.names = FALSE)) else integer(0)
      unless_guard <- if (length(unless_ids) > 0) unique(unlist(lapply(unless_ids, function(uid) nodes[[uid]][["guard_ids"]] %||% integer(0)), use.names = FALSE)) else integer(0)
      node[["complete_ids"]] <- if (length(ref_complete) > 0) as.integer(sort(unique(ref_complete))) else integer(0)
      node[["guard_ids"]] <- if (length(ref_guard) > 0) as.integer(sort(unique(c(ref_guard, unless_guard)))) else if (length(unless_guard) > 0) as.integer(sort(unique(unless_guard))) else integer(0)
      node[["finish_ids"]] <- if (length(ref_finish) > 0) as.integer(sort(unique(ref_finish))) else integer(0)
      node[["blocker_complete_ids"]] <- if (length(blocker_complete) > 0) as.integer(sort(unique(blocker_complete))) else integer(0)
      node[["blocker_guard_ids"]] <- if (length(blocker_guard) > 0) as.integer(sort(unique(blocker_guard))) else integer(0)
      node[["blocker_finish_ids"]] <- if (length(blocker_finish) > 0) as.integer(sort(unique(blocker_finish))) else integer(0)
      node[["unless_complete_ids"]] <- if (length(unless_complete) > 0) as.integer(sort(unique(unless_complete))) else integer(0)
      node[["unless_guard_ids"]] <- if (length(unless_guard) > 0) as.integer(sort(unique(unless_guard))) else integer(0)
    } else {
      node[["args"]] <- list()
      node[["complete_ids"]] <- node[["complete_ids"]] %||% integer(0)
      node[["guard_ids"]] <- node[["guard_ids"]] %||% integer(0)
      node[["finish_ids"]] <- node[["finish_ids"]] %||% integer(0)
    }

    if (identical(kind, "event")) {
      node[["source_label"]] <- expr[["source"]] %||% NA_character_
      node[["source_idx"]] <- .label_to_id(prep, node[["source_label"]])
      source_ids <- node[["sources"]] %||% integer(0)
      node[["complete_ids"]] <- if (length(source_ids) > 0) as.integer(sort(unique(source_ids))) else integer(0)
      guard_ids <- source_ids
      finish_ids <- source_ids
      src_idx <- node[["source_idx"]] %||% NA_integer_
      if (!is.na(src_idx)) {
        guard_ids <- c(guard_ids, as.integer(src_idx))
        finish_ids <- c(finish_ids, as.integer(src_idx))
      }
      node[["guard_ids"]] <- if (length(guard_ids) > 0) as.integer(sort(unique(guard_ids))) else integer(0)
      node[["finish_ids"]] <- if (length(finish_ids) > 0) as.integer(sort(unique(finish_ids))) else integer(0)
    }

    if (is.null(node[["complete_ids"]])) node[["complete_ids"]] <- integer(0)
    if (is.null(node[["guard_ids"]])) node[["guard_ids"]] <- integer(0)
    if (is.null(node[["finish_ids"]])) node[["finish_ids"]] <- integer(0)

    nodes[[node_id]] <<- node
    as.integer(node_id)
  }

  # Register all expressions referenced in outcomes and competitor maps
  outcomes <- prep[["outcomes"]] %||% list()
  for (expr in lapply(outcomes, `[[`, "expr")) {
    register_expr(expr)
  }
  comp_map <- prep[[".competitors"]] %||% list()
  if (length(comp_map) > 0) {
    for (expr_list in comp_map) {
      if (length(expr_list) == 0) next
      for (expr in expr_list) {
        register_expr(expr)
      }
    }
  }

  registry <- list(
    nodes = nodes,
    index = index
  )

  registry <- .compile_expr_registry(prep, registry)
  registry
}

.compile_expr_registry <- function(prep, registry) {
  nodes <- registry[["nodes"]] %||% list()
  n <- length(nodes)
  if (n == 0) {
    registry[["scenario_fns"]] <- list()
    registry[["cdf_cond_fns"]] <- list()
    registry[["survival_cond_fns"]] <- list()
    return(registry)
  }

  scenario_memo <- vector("list", n)
  cdf_cond_memo <- vector("list", n)
  surv_cond_memo <- vector("list", n)

  build_scenario_fn <- NULL
  build_cdf_cond_fn <- NULL
  build_surv_cond_fn <- NULL

  build_scenario_fn <- function(id) {
    if (is.null(id) || is.na(id) || id <= 0 || id > n) return(function(...) list())
    cached <- scenario_memo[[id]]
    if (!is.null(cached)) return(cached)
    node <- nodes[[id]]
    kind <- node[["kind"]] %||% NA_character_

    # Re-entrancy guard to break potential cycles during construction
    scenario_memo[[id]] <<- function(...) list()

    fn <- switch(
      kind,
      "event" = {
        source_label <- node[["source_label"]] %||% NA_character_
        source_idx <- node[["source_idx"]] %||% NA_integer_
        complete_ids <- node[["complete_ids"]] %||% integer(0)
        finish_ids <- node[["finish_ids"]] %||% integer(0)
        complete_ids <- .forced_union(prep, complete_ids, if (!is.na(source_idx)) as.integer(source_idx) else integer(0))
        finish_ids <- .forced_union(prep, finish_ids, if (!is.na(source_idx)) as.integer(source_idx) else integer(0))
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (!is.null(ctx)) {
            key <- .ctx_key("sc", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$scenario_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          forced_complete <- .coerce_forced_ids(prep, forced_complete)
          forced_survive <- .coerce_forced_ids(prep, forced_survive)
          weight <- .event_density_at(
            prep, source_label, component, t,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            ctx = ctx
          )
          if (!is.finite(weight) || weight <= 0) return(list())
          pool_scenarios <- attr(weight, "scenarios")
          if (!is.null(pool_scenarios) && length(pool_scenarios) > 0) {
            out <- list()
            for (psc in pool_scenarios) {
              if (is.null(psc)) next
              w <- psc$weight
              if (!is.finite(w) || w <= 0) next
              fcomp <- .forced_union(prep, forced_complete, .forced_union(prep, complete_ids, psc$forced_complete %||% integer(0)))
              fsurv <- .forced_union(prep, forced_survive, psc$forced_survive %||% integer(0))
              fins <- .forced_union(prep, finish_ids, psc$finish_ids %||% integer(0))
              sc <- list(
                weight = as.numeric(w),
                forced_complete = .coerce_forced_ids(prep, fcomp),
                forced_survive = .coerce_forced_ids(prep, fsurv),
                finish_ids = .coerce_forced_ids(prep, fins)
              )
              out[[length(out) + 1L]] <- sc
            }
            if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- out
            return(out)
          }
          sc <- list(
            weight = as.numeric(weight),
            forced_complete = .coerce_forced_ids(prep, .forced_union(prep, forced_complete, complete_ids)),
            forced_survive = .coerce_forced_ids(prep, forced_survive),
            finish_ids = .coerce_forced_ids(prep, finish_ids)
          )
          out <- list(sc)
          if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- out
          out
        }
      },
      "and" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) list()
        } else {
          child_scenario <- lapply(child_ids, build_scenario_fn)
          child_cdf <- lapply(child_ids, build_cdf_cond_fn)
          child_sources <- lapply(child_ids, function(cid) nodes[[cid]][["sources"]] %||% character(0))
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("sc", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$scenario_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            forced_complete <- .coerce_forced_ids(prep, forced_complete)
            forced_survive <- .coerce_forced_ids(prep, forced_survive)
            out <- list()
            idx_all <- seq_along(child_ids)
            for (i in idx_all) {
              si <- child_scenario[[i]](t, component, ctx,
                                        forced_complete = forced_complete,
                                        forced_survive = forced_survive)
              if (length(si) == 0) next
              others_idx <- idx_all[idx_all != i]
              for (sc in si) {
                if (is.null(sc) || sc$weight <= 0) next
                weight <- sc$weight
                fcomp <- sc$forced_complete
                fsurv <- sc$forced_survive
                finish_ids <- sc$finish_ids %||% integer(0)
                ok <- TRUE
                if (length(others_idx) > 0) {
                  for (j in others_idx) {
                    Fj <- child_cdf[[j]](
                      t = t,
                      component = component,
                      ctx = ctx,
                      forced_complete = fcomp,
                      forced_survive = fsurv
                    )
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
                out[[length(out) + 1L]] <- list(
                  weight = as.numeric(weight),
                  forced_complete = .coerce_forced_ids(prep, fcomp),
                  forced_survive = .coerce_forced_ids(prep, fsurv),
                  finish_ids = .coerce_forced_ids(prep, finish_ids)
                )
              }
            }
            if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- out
            out
          }
        }
      },
      "or" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) list()
        } else {
          child_scenario <- lapply(child_ids, build_scenario_fn)
          child_sources <- lapply(child_ids, function(cid) nodes[[cid]][["sources"]] %||% character(0))
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("sc", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$scenario_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            forced_complete <- .coerce_forced_ids(prep, forced_complete)
            forced_survive <- .coerce_forced_ids(prep, forced_survive)
            out <- list()
            idx_all <- seq_along(child_ids)
            for (i in idx_all) {
              si <- child_scenario[[i]](t, component, ctx,
                                        forced_complete = forced_complete,
                                        forced_survive = forced_survive)
              if (length(si) == 0) next
              others_idx <- idx_all[idx_all != i]
              req_list <- child_sources[others_idx]
             for (sc in si) {
                if (is.null(sc) || sc$weight <= 0) next
                weight <- sc$weight
                fcomp <- sc$forced_complete
                fsurv <- sc$forced_survive
                finish_ids <- sc$finish_ids %||% integer(0)
                valid <- TRUE
                if (length(others_idx) > 0) {
                  for (k in seq_along(others_idx)) {
                    req_j <- req_list[[k]]
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
                out[[length(out) + 1L]] <- list(
                  weight = as.numeric(weight),
                  forced_complete = .coerce_forced_ids(prep, fcomp),
                  forced_survive = .coerce_forced_ids(prep, fsurv),
                  finish_ids = .coerce_forced_ids(prep, finish_ids)
                )
              }
            }
            if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- out
            out
          }
        }
      },
      "guard" = {
        ref_id <- node[["reference_id"]] %||% NA_integer_
        blocker_id <- node[["blocker_id"]] %||% NA_integer_
        unless_ids <- node[["unless_ids"]] %||% integer(0)
        blocker_sources <- if (!is.na(blocker_id)) (nodes[[blocker_id]][["sources"]] %||% integer(0)) else integer(0)
        self_complete <- node[["complete_ids"]] %||% integer(0)
        self_finish <- node[["finish_ids"]] %||% integer(0)
        ref_scen_fn <- if (!is.na(ref_id)) build_scenario_fn(ref_id) else function(...) list()
        surv_blocker_fn <- if (!is.na(blocker_id)) build_surv_cond_fn(blocker_id) else NULL
        surv_unless_fns <- if (length(unless_ids) > 0) lapply(unless_ids, build_surv_cond_fn) else list()
        blocker_density_at <- if (!is.na(blocker_id)) {
          bl_scen_fn <- build_scenario_fn(blocker_id)
          function(u, component, ctx, forced_complete, forced_survive) {
            scs <- bl_scen_fn(u, component, ctx,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive)
            if (length(scs) == 0) return(0.0)
            sum(vapply(scs, `[[`, numeric(1), "weight"))
          }
        } else {
          function(u, component, ctx, forced_complete, forced_survive) 0.0
        }
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (!is.null(ctx)) {
            key <- .ctx_key("sc", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$scenario_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          forced_complete <- .coerce_forced_ids(prep, forced_complete)
          forced_survive <- .coerce_forced_ids(prep, forced_survive)
          ref_scen <- ref_scen_fn(t, component, ctx,
                                   forced_complete = forced_complete,
                                   forced_survive = forced_survive)
          if (length(ref_scen) == 0) {
            if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- list()
            return(list())
          }
          S_eff <- (function() {
            if (is.na(blocker_id)) return(1.0)
            if (length(unless_ids) == 0) {
              return(surv_blocker_fn(t, component, ctx,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive))
            }
            if (!is.finite(t) || t <= 0) return(1.0)
            integrand <- function(u) {
              bd <- blocker_density_at(u, component, ctx,
                                       forced_complete = forced_complete,
                                       forced_survive = forced_survive)
              if (!is.finite(bd) || bd <= 0) return(0.0)
              prods <- 1.0
              if (length(surv_unless_fns) > 0) {
                for (fn in surv_unless_fns) {
                  sv <- fn(u, component, ctx,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive)
                  prods <- prods * sv
                  if (!is.finite(prods) || prods == 0) break
                }
              }
              bd * prods
            }
            val <- tryCatch(stats::integrate(integrand, lower = 0, upper = t,
                                             rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                             stop.on.error = FALSE)[["value"]],
                            error = function(e) 0.0)
            s <- 1.0 - as.numeric(val)
            if (!is.finite(s)) s <- 0.0
            if (s < 0) s <- 0.0
            if (s > 1) s <- 1.0
            s
          })()
          scenarios <- list()
          for (sc in ref_scen) {
            if (is.null(sc)) next
            w <- sc$weight * S_eff
            if (!is.finite(w) || w <= 0) next
            fcomp <- .forced_union(prep, sc$forced_complete, self_complete)
            fsurv <- .forced_union(prep, sc$forced_survive, blocker_sources)
            finish_ids <- .forced_union(prep, sc$finish_ids %||% integer(0), self_finish)
            scenarios[[length(scenarios) + 1L]] <- list(
              weight = as.numeric(w),
              forced_complete = .coerce_forced_ids(prep, fcomp),
              forced_survive = .coerce_forced_ids(prep, fsurv),
              finish_ids = .coerce_forced_ids(prep, finish_ids)
            )
          }
          if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("sc", id, component, t, forced_complete, forced_survive)]] <- scenarios
          scenarios
        }
      },
      {
        function(...) list()
      }
    )

    scenario_memo[[id]] <<- fn
    fn
  }

  build_cdf_cond_fn <- function(id) {
    if (is.null(id) || is.na(id) || id <= 0 || id > n) return(function(...) 0.0)
    cached <- cdf_cond_memo[[id]]
    if (!is.null(cached)) return(cached)
    node <- nodes[[id]]
    kind <- node[["kind"]] %||% NA_character_

    # Re-entrancy guard
    cdf_cond_memo[[id]] <<- function(...) 0.0

    fn <- switch(
      kind,
      "event" = {
        src <- node[["source_label"]] %||% NA_character_
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (!is.null(ctx)) {
            key <- .ctx_key("cdf", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$cdf_cond_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          forced_complete <- .coerce_forced_ids(prep, forced_complete)
          forced_survive <- .coerce_forced_ids(prep, forced_survive)
          val <- .event_cdf_at(
            prep, src, component, t,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            ctx = ctx
          )
          if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdf", id, component, t, forced_complete, forced_survive)]] <- val
          val
        }
      },
      "and" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) 0.0
        } else {
          child_cdf <- lapply(child_ids, build_cdf_cond_fn)
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("cdf", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$cdf_cond_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            forced_complete <- .coerce_forced_ids(prep, forced_complete)
            forced_survive <- .coerce_forced_ids(prep, forced_survive)
            prod_val <- 1.0
            for (fn_child in child_cdf) {
              val <- fn_child(t, component, ctx,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive)
              prod_val <- prod_val * val
              if (prod_val == 0) break
            }
            if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdf", id, component, t, forced_complete, forced_survive)]] <- prod_val
            prod_val
          }
        }
      },
      "or" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) 0.0
        } else {
          child_surv <- lapply(child_ids, build_surv_cond_fn)
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("cdf", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$cdf_cond_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            forced_complete <- .coerce_forced_ids(prep, forced_complete)
            forced_survive <- .coerce_forced_ids(prep, forced_survive)
            prod_surv <- 1.0
            for (fn_child in child_surv) {
              val <- fn_child(t, component, ctx,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive)
              prod_surv <- prod_surv * val
              if (prod_surv == 0) break
            }
            val <- 1.0 - prod_surv
            if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdf", id, component, t, forced_complete, forced_survive)]] <- val
            val
          }
        }
      },
      "guard" = {
        scen_self <- NULL
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (is.null(scen_self)) scen_self <<- build_scenario_fn(id)
          if (!is.null(ctx)) {
            key <- .ctx_key("cdf", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$cdf_cond_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          if (t <= 0) {
            if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdf", id, component, t, forced_complete, forced_survive)]] <- 0.0
            return(0.0)
          }
          dens_fun <- function(u) {
            scs <- scen_self(u, component, ctx,
                             forced_complete = forced_complete,
                             forced_survive = forced_survive)
            if (length(scs) == 0) return(0.0)
            sum(vapply(scs, `[[`, numeric(1), "weight"))
          }
          val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                           rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                           stop.on.error = FALSE)[["value"]],
                          error = function(e) 0.0)
          out <- as.numeric(val)
          if (!is.finite(out)) out <- 0.0
          out <- max(0.0, min(1.0, out))
          if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdf", id, component, t, forced_complete, forced_survive)]] <- out
          out
        }
      },
      {
        function(...) 0.0
      }
    )

    cdf_cond_memo[[id]] <<- fn
    fn
  }

  build_surv_cond_fn <- function(id) {
    if (is.null(id) || is.na(id) || id <= 0 || id > n) return(function(...) 1.0)
    cached <- surv_cond_memo[[id]]
    if (!is.null(cached)) return(cached)
    node <- nodes[[id]]
    kind <- node[["kind"]] %||% NA_character_

    # Re-entrancy guard
    surv_cond_memo[[id]] <<- function(...) 1.0

    fn <- switch(
      kind,
      "event" = {
        src <- node[["source_label"]] %||% NA_character_
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (!is.null(ctx)) {
            key <- .ctx_key("surv", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$surv_cond_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          forced_complete <- .coerce_forced_ids(prep, forced_complete)
          forced_survive <- .coerce_forced_ids(prep, forced_survive)
          val <- .event_survival_at(
            prep, src, component, t,
            forced_complete = forced_complete,
            forced_survive = forced_survive,
            ctx = ctx
          )
          if (!is.null(ctx)) ctx$surv_cond_memo[[.ctx_key("surv", id, component, t, forced_complete, forced_survive)]] <- val
          val
        }
      },
      "and" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) 1.0
        } else {
          self_cdf <- build_cdf_cond_fn(id)
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("surv", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$surv_cond_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            val <- 1.0 - self_cdf(
              t, component, ctx,
              forced_complete = forced_complete,
              forced_survive = forced_survive
            )
            if (!is.null(ctx)) ctx$surv_cond_memo[[.ctx_key("surv", id, component, t, forced_complete, forced_survive)]] <- val
            val
          }
        }
      },
      "or" = {
        child_ids <- node[["child_ids"]] %||% integer(0)
        if (length(child_ids) == 0) {
          function(...) 1.0
        } else {
          child_surv <- lapply(child_ids, build_surv_cond_fn)
          function(t, component,
                   ctx,
                   forced_complete = integer(0),
                   forced_survive = integer(0)) {
            if (!is.null(ctx)) {
              key <- .ctx_key("surv", id, component, t, forced_complete, forced_survive)
              cached_val <- ctx$surv_cond_memo[[key]]
              if (!is.null(cached_val)) return(cached_val)
            }
            forced_complete <- .coerce_forced_ids(prep, forced_complete)
            forced_survive <- .coerce_forced_ids(prep, forced_survive)
            prod_val <- 1.0
            for (fn_child in child_surv) {
              val <- fn_child(t, component, ctx,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive)
              prod_val <- prod_val * val
              if (!is.finite(prod_val) || prod_val == 0) break
            }
            if (!is.null(ctx)) ctx$surv_cond_memo[[.ctx_key("surv", id, component, t, forced_complete, forced_survive)]] <- prod_val
            prod_val
          }
        }
      },
      "guard" = {
        self_cdf <- build_cdf_cond_fn(id)
        function(t, component,
                 ctx,
                 forced_complete = integer(0),
                 forced_survive = integer(0)) {
          if (!is.null(ctx)) {
            key <- .ctx_key("surv", id, component, t, forced_complete, forced_survive)
            cached_val <- ctx$surv_cond_memo[[key]]
            if (!is.null(cached_val)) return(cached_val)
          }
          val <- 1.0 - self_cdf(t, component, ctx,
                                 forced_complete = forced_complete,
                                 forced_survive = forced_survive)
          if (!is.null(ctx)) ctx$surv_cond_memo[[.ctx_key("surv", id, component, t, forced_complete, forced_survive)]] <- val
          val
        }
      },
      {
        function(...) 1.0
      }
    )

    surv_cond_memo[[id]] <<- fn
    fn
  }

  scenario_fns <- lapply(seq_len(n), build_scenario_fn)
  cdf_cond_fns <- lapply(seq_len(n), build_cdf_cond_fn)
  survival_cond_fns <- lapply(seq_len(n), build_surv_cond_fn)

  registry[["scenario_fns"]] <- scenario_fns
  registry[["cdf_cond_fns"]] <- cdf_cond_fns
  registry[["survival_cond_fns"]] <- survival_cond_fns
  registry
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

# ============================================================================
# PART 2a: Evaluation Context (per-trial) and cache keys
# ============================================================================

# Create a per-trial evaluation context that holds memo tables only for
# the duration of a single trial.
make_eval_ctx <- function(prep) {
  list(
    scenario_memo = new.env(parent = emptyenv(), hash = TRUE),
    cdf_cond_memo = new.env(parent = emptyenv(), hash = TRUE),
    surv_cond_memo = new.env(parent = emptyenv(), hash = TRUE),
    event_density_memo = new.env(parent = emptyenv(), hash = TRUE),
    event_cdf_memo = new.env(parent = emptyenv(), hash = TRUE),
    event_surv_memo = new.env(parent = emptyenv(), hash = TRUE),
    pool_density_memo = new.env(parent = emptyenv(), hash = TRUE),
    pool_surv_memo = new.env(parent = emptyenv(), hash = TRUE)
  )
}

# Lightweight diagnostic helpers for registry/ID mapping and cache pressure
.diagnose_registry_mapping <- function(prep) {
  reg <- prep[[".expr_registry"]]
  outcomes <- prep[["outcomes"]] %||% list()
  if (is.null(reg) || length(outcomes) == 0) return(invisible(NULL))
  bad <- character(0)
  for (lbl in names(outcomes)) {
    ex <- outcomes[[lbl]][["expr"]]
    id <- .expr_registry_lookup_id(reg, ex)
    if (is.na(id) || id <= 0 || id > length(reg[["nodes"]] %||% list())) {
      bad <- c(bad, lbl)
    }
  }
  if (length(bad) > 0) {
    warning(sprintf("Expr->node mapping failed for: %s", paste(bad, collapse = ", ")))
  }
  invisible(NULL)
}

.ctx_env_size <- function(env) {
  if (is.null(env) || !is.environment(env)) return(0L)
  length(ls(env, all.names = TRUE))
}

.ctx_pressure <- function(ctx) {
  if (is.null(ctx)) return(0L)
  .ctx_env_size(ctx$scenario_memo) +
    .ctx_env_size(ctx$cdf_cond_memo) +
    .ctx_env_size(ctx$surv_cond_memo) +
    .ctx_env_size(ctx$event_density_memo) +
    .ctx_env_size(ctx$event_cdf_memo) +
    .ctx_env_size(ctx$event_surv_memo) +
    .ctx_env_size(ctx$pool_density_memo) +
    .ctx_env_size(ctx$pool_surv_memo)
}

.ids_key <- function(ids) {
  if (is.null(ids) || length(ids) == 0L) return("")
  ids <- as.integer(ids)
  if (length(ids) > 1L) ids <- sort(unique(ids))
  paste(ids, collapse = ",")
}

.ctx_key <- function(prefix, node_id_or_label, component, t, forced_complete, forced_survive) {
  t_key <- .time_key(t)
  fc_key <- .ids_key(forced_complete)
  fs_key <- .ids_key(forced_survive)
  paste(prefix, as.character(node_id_or_label), as.character(component %||% "__default__"), t_key, fc_key, fs_key, sep = "|")
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
      finish_ids <- c(member_ids[idx], if (!is.na(pool_idx)) pool_idx else integer(0))
      templates[[template_count]] <- list(
        finisher_idx = idx,
        complete_idx = if (length(combo_idx) > 0L) combo_idx else integer(0),
        survivor_idx = survivors,
        forced_complete_ids = finisher_ids,
        forced_survive_ids = if (length(survivors) > 0L) member_ids[survivors] else integer(0),
        finish_ids = finish_ids
      )
      idx_entries[[j]] <- template_count
    }
    finisher_map[[idx]] <- idx_entries
  }
  list(templates = templates, finisher_map = finisher_map)
}

.event_survival_at <- function(prep, id, component, t,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               ctx = NULL) {
  if (!is.null(ctx)) {
    key <- .ctx_key("evS", id, component, t, forced_complete, forced_survive)
    cached <- ctx$event_surv_memo[[key]]
    if (!is.null(cached)) return(cached)
  }
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(1.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    val <- .pool_survival(prep, id, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive)
    if (!is.null(ctx)) ctx$event_surv_memo[[key]] <- val
    return(val)
  }
  acc_def <- prep[["accumulators"]][[id]]
  val <- if (!is.null(acc_def)) .acc_survival(acc_def, t) else 0.0
  if (!is.null(ctx)) ctx$event_surv_memo[[key]] <- val
  val
}

.event_cdf_at <- function(prep, id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0),
                          ctx = NULL) {
  if (!is.null(ctx)) {
    key <- .ctx_key("evF", id, component, t, forced_complete, forced_survive)
    cached <- ctx$event_cdf_memo[[key]]
    if (!is.null(cached)) return(cached)
  }
  val <- 1.0 - .event_survival_at(prep, id, component, t,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive,
                           ctx = ctx)
  if (!is.null(ctx)) ctx$event_cdf_memo[[.ctx_key("evF", id, component, t, forced_complete, forced_survive)]] <- val
  val
}

.event_density_at <- function(prep, id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0),
                              ctx = NULL) {
  if (!is.null(ctx)) {
    key <- .ctx_key("evD", id, component, t, forced_complete, forced_survive)
    cached <- ctx$event_density_memo[[key]]
    if (!is.null(cached)) return(cached)
  }
  id_idx <- .label_to_id(prep, id)
  if (!is.na(id_idx)) {
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) return(0.0)
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) return(0.0)
  }
  is_pool <- !is.null(prep[["pools"]][[id]])
  if (is_pool) {
    val <- .pool_density(prep, id, component, t,
                         forced_complete = forced_complete,
                         forced_survive = forced_survive)
    if (!is.null(ctx)) ctx$event_density_memo[[.ctx_key("evD", id, component, t, forced_complete, forced_survive)]] <- val
    return(val)
  }
  acc_def <- prep[["accumulators"]][[id]]
  val <- if (!is.null(acc_def)) .acc_density(acc_def, t) else 0.0
  if (!is.null(ctx)) ctx$event_density_memo[[.ctx_key("evD", id, component, t, forced_complete, forced_survive)]] <- val
  val
}

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0),
                          ctx = NULL) {
  if (!is.null(ctx)) {
    key <- .ctx_key("poolD", pool_id, component, t, forced_complete, forced_survive)
    cached <- ctx$pool_density_memo[[key]]
    if (!is.null(cached)) return(cached)
  }
  if (!is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    if (!is.null(ctx)) ctx$pool_density_memo[[.ctx_key("poolD", pool_id, component, t, forced_complete, forced_survive)]] <- val
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
      forced_survive = forced_survive,
      ctx = ctx)
    cdf_vec[[i]] <- .event_cdf_at(
      prep, mid, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      ctx = ctx)
    surv_vec[[i]] <- 1.0 - cdf_vec[[i]]
  }

  scenarios <- list()
  add_scenario <- function(weight, fcomp, fsurv, finish_ids = integer(0)) {
    if (!is.finite(weight) || weight <= 0) return()
    scenarios[[length(scenarios) + 1L]] <<- list(
      weight = as.numeric(weight),
      forced_complete = .coerce_forced_ids(prep, fcomp),
      forced_survive = .coerce_forced_ids(prep, fsurv),
      finish_ids = .coerce_forced_ids(prep, finish_ids)
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
        add_scenario(weight, finisher_ids, integer(0), finisher_ids)
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
        finish_ids <- c(member_ids[idx], if (!is.na(pool_idx)) pool_idx else integer(0))
        add_scenario(weight, finisher_ids, survivor_ids, finish_ids)
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
      finish_ids <- c(member_ids[idx], if (!is.na(pool_idx)) pool_idx else integer(0))
      add_scenario(weight, finisher_ids, integer(0), finish_ids)
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
        add_scenario(weight, tmpl$forced_complete_ids, tmpl$forced_survive_ids, tmpl$finish_ids %||% integer(0))
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
  if (!is.null(ctx)) ctx$pool_density_memo[[.ctx_key("poolD", pool_id, component, t, forced_complete, forced_survive)]] <- total_density
  total_density
}

.pool_survival <- function(prep, pool_id, component, t,
                           forced_complete = integer(0),
                           forced_survive = integer(0),
                           ctx = NULL) {
  if (!is.null(ctx)) {
    key <- .ctx_key("poolS", pool_id, component, t, forced_complete, forced_survive)
    cached <- ctx$pool_surv_memo[[key]]
    if (!is.null(cached)) return(cached)
  }
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
                               forced_survive = forced_survive,
                               ctx = ctx)
  }
  Svec <- 1.0 - Fvec
  coeffs <- pool_coeffs(Svec, Fvec)
  upto <- min(length(coeffs), k)
  val <- sum(coeffs[seq_len(upto)])
  if (!is.null(ctx)) ctx$pool_surv_memo[[.ctx_key("poolS", pool_id, component, t, forced_complete, forced_survive)]] <- val
  val
}

# ==============================================================================
# PART 4: Expression Evaluation (Likelihood)
# ==============================================================================
 
.scenario_density_with_competitors <- function(expr, t, prep, component,
                                               competitor_exprs = list(),
                                               ctx = NULL) {
  if (is.na(t) || !is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  scenarios <- .expr_scenarios_at(
    expr, t, prep, component,
    forced_complete = integer(0),
    forced_survive = integer(0),
    ctx = ctx)
  if (length(scenarios) == 0) {
    return(0.0)
  }
  registry <- prep[[".expr_registry"]]
  nodes <- registry[["nodes"]] %||% list()
  comp_meta <- NULL
  if (length(competitor_exprs) > 0) {
    comp_meta <- lapply(competitor_exprs, function(expr_i) {
      meta <- list(expr = expr_i,
                   kind = NA_character_,
                   finish_ids = integer(0),
                   guard_blocker_finish = integer(0))
      if (!is.null(registry)) {
        expr_id <- .expr_registry_lookup_id(registry, expr_i)
        if (!is.na(expr_id) && expr_id > 0 && expr_id <= length(nodes)) {
          node <- nodes[[expr_id]]
          meta$kind <- node[["kind"]] %||% NA_character_
          meta$finish_ids <- node[["finish_ids"]] %||% integer(0)
          if (identical(meta$kind, "guard")) {
            meta$guard_blocker_finish <- node[["blocker_finish_ids"]] %||% integer(0)
          }
        }
      }
      meta
    })
  }

  total <- 0.0
  for (sc in scenarios) {
    if (is.null(sc) || sc$weight <= 0) next
    weight <- sc$weight
    finish_ids <- .coerce_forced_ids(prep, sc$finish_ids %||% integer(0))
    forced_complete <- sc$forced_complete
    forced_survive <- sc$forced_survive %||% integer(0)
    if (!is.null(comp_meta) && length(comp_meta) > 0) {
      surv_prod <- 1.0
      for (meta in comp_meta) {
        expr_i <- meta$expr
        kind_i <- meta$kind
        guard_blocker_finish <- .coerce_forced_ids(prep, meta$guard_blocker_finish)
        fsurv_adj <- forced_survive
        surv_val <- 1.0
        if (identical(kind_i, "guard")) {
          if (length(guard_blocker_finish) == 0 || length(finish_ids) == 0 ||
              !any(finish_ids %in% guard_blocker_finish)) {
            surv_val <- .eval_expr_survival_cond(
              expr_i, t, prep, component,
              forced_complete = forced_complete,
              forced_survive = fsurv_adj,
              ctx = ctx)
          }
        } else {
          comp_finish_ids <- .coerce_forced_ids(prep, meta$finish_ids)
          if (length(comp_finish_ids) > 0 && length(finish_ids) > 0) {
            finish_conflict <- finish_ids[finish_ids %in% comp_finish_ids]
            if (length(finish_conflict) > 0) {
              fsurv_adj <- .forced_union(prep, fsurv_adj, finish_conflict)
            }
          }
          surv_val <- .eval_expr_survival_cond(
            expr_i, t, prep, component,
            forced_complete = forced_complete,
            forced_survive = fsurv_adj,
            ctx = ctx)
        }
        if (!is.finite(surv_val) || surv_val <= 0) {
          surv_prod <- 0.0
          break
        }
        surv_prod <- surv_prod * surv_val
      }
      weight <- weight * surv_prod
    }
    total <- total + weight
  }
  if (length(scenarios) > 0) attr(total, "scenarios") <- scenarios
  total
}

 

# ==============================================================================
# PART 5: Outcome Likelihood
# ==============================================================================

# Integrand for scenario-conditioned probability (density already encodes competitors)
.integrand_outcome_density <- function(t, expr, prep, component, competitor_exprs, ctx) {
  # To prevent unbounded memo growth during quadrature, use a fresh
  # per-t evaluation context unless explicitly disabled.
  use_fresh_ctx <- getOption("uuber.integrate.fresh_ctx_per_t", default = TRUE)
  vapply(t, function(tt) {
    if (!is.finite(tt) || tt < 0) return(0.0)
    local_ctx <- if (isTRUE(use_fresh_ctx)) make_eval_ctx(prep) else ctx
    .scenario_density_with_competitors(expr, tt, prep, component,
                                       competitor_exprs = competitor_exprs,
                                       ctx = local_ctx)
  }, numeric(1))
}

# Scenario-conditioned probability of an expression over [0, upper_limit].
# competitor_exprs retained for backward compatibility but ignored.
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf, competitor_exprs = NULL, ctx = NULL) {
  if (!is.finite(upper_limit)) upper_limit <- Inf
  res <- tryCatch(
    stats::integrate(
      .integrand_outcome_density,
      lower = 0,
      upper = upper_limit,
      expr = expr,
      prep = prep,
      component = component,
      competitor_exprs = competitor_exprs %||% list(),
      ctx = if (isTRUE(getOption("uuber.integrate.fresh_ctx_per_t", default = TRUE))) NULL else ctx,
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
        # Match legacy: if ALL target sources are within the guard's blocker set,
        # then use the guard's reference as the effective competitor.
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
  if (length(comps) > 0) {
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

# Compute likelihood for a single trial/outcome
.outcome_likelihood <- function(outcome_label, rt, prep, component, ctx = NULL) {
  outcome_defs <- prep[["outcomes"]]

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
                                                           ctx = ctx)
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
                                                         ctx = ctx)
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
              comp_exprs_map <- prep[['.competitors']][[label]] %||% list()
              return(.integrate_outcome_probability(expr, prep, component, deadline,
                                                    competitor_exprs = comp_exprs_map,
                                                    ctx = ctx))
            } else {
              # Race density at rt for the mapped source outcome
              comp_exprs_map <- prep[['.competitors']][[label]] %||% list()
              dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                           competitor_exprs = comp_exprs_map,
                                                           ctx = ctx)
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
        surv <- .eval_expr_survival_cond(expr, deadline, prep, component, ctx = ctx)
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

  # Handle NA/infinite RT by integration of race density
  if (is.na(rt) || !is.finite(rt) || rt < 0) {
    deadline <- .get_component_attr(prep, component, "deadline")
    deadline <- deadline %||% prep[["default_deadline"]]
    integrand <- function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        base <- .scenario_density_with_competitors(expr, tt, prep, component,
                                                   competitor_exprs = competitor_exprs,
                                                   ctx = ctx)
        if (base == 0) return(0.0)
        add <- 0.0
        if (length(donors) > 0) {
          for (d in donors) {
            dens_d <- .scenario_density_with_competitors(d[['expr']], tt, prep, component,
                                                          competitor_exprs = d[['competitors']] %||% list(),
                                                          ctx = ctx)
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
  dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                               competitor_exprs = competitor_exprs,
                                               ctx = ctx)
  if (dens_r == 0) return(0.0)
  base_val <- dens_r
  # Include donor mass (e.g., timeout guesses) that keeps the observed RT
  if (length(donors) > 0) {
    donor_add <- 0.0
    for (d in donors) {
      if (identical(d[['rt_policy']], "na")) next
      dens_d <- .scenario_density_with_competitors(d[['expr']], rt, prep, component,
                                                   competitor_exprs = d[['competitors']] %||% list(),
                                                   ctx = ctx)
      if (dens_d == 0) next
      donor_add <- donor_add + d[['weight']] * dens_d
    }
    base_val <- base_val + donor_add
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
  # Sanity check: ensure expressions are mapped to compiled node IDs
  .diagnose_registry_mapping(prep)

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
      ctx <- make_eval_ctx(prep)
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
              lik <- lik + base_weights[[j]] * .outcome_likelihood(outcome, rt_val, prep, comp_ids[[j]], ctx = ctx)
            }
          }
        } else {
          lik <- .outcome_likelihood(outcome, rt_val, prep, component_key, ctx = ctx)
        }
      } else {
        lik <- .outcome_likelihood(outcome, rt_val, prep, component_key, ctx = ctx)
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
    # Alias-of outcomes: sum of referenced labels
    out_def <- prep[["outcomes"]][[outcome_label]]
    if (!is.null(out_def[['options']][['alias_of']])) {
      refs <- out_def[['options']][['alias_of']]
      refs <- as.character(refs)
      vals <- vapply(refs, function(lbl) as.numeric(.outcome_likelihood(lbl, NA_real_, prep, comp_id, ctx = make_eval_ctx(prep))), numeric(1))
      return(sum(vals))
    }
    base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, comp_id, ctx = make_eval_ctx(prep)))
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
.expr_scenarios_at <- function(expr, t, prep, component,
                               forced_complete = integer(0),
                               forced_survive = integer(0),
                               ctx = NULL) {
  registry <- prep[[".expr_registry"]]
  if (!is.null(registry)) {
    expr_id <- .expr_registry_lookup_id(registry, expr)
    if (!is.na(expr_id)) {
      if (!is.null(ctx)) {
        key <- .ctx_key("scId", expr_id, component, t, forced_complete, forced_survive)
        cached <- ctx$scenario_memo[[key]]
        if (!is.null(cached)) return(cached)
      }
      scen_fns <- registry[["scenario_fns"]]
      if (!is.null(scen_fns) && length(scen_fns) >= expr_id) {
        fn <- scen_fns[[expr_id]]
        if (is.function(fn)) {
          res <- fn(
            t = t,
            component = component,
            ctx = ctx,
            forced_complete = forced_complete,
            forced_survive = forced_survive
          )
          if (!is.null(ctx)) ctx$scenario_memo[[.ctx_key("scId", expr_id, component, t, forced_complete, forced_survive)]] <- res
          return(res)
        }
      }
    }
  }
  list()
}

.eval_expr_cdf_cond <- function(expr, t, prep, component,
                                forced_complete = integer(0),
                                forced_survive = integer(0),
                                ctx = NULL) {
  registry <- prep[[".expr_registry"]]
  if (!is.null(registry)) {
    expr_id <- .expr_registry_lookup_id(registry, expr)
    if (!is.na(expr_id)) {
      if (!is.null(ctx)) {
        key <- .ctx_key("cdfId", expr_id, component, t, forced_complete, forced_survive)
        cached <- ctx$cdf_cond_memo[[key]]
        if (!is.null(cached)) return(cached)
      }
      cdf_fns <- registry[["cdf_cond_fns"]]
      if (!is.null(cdf_fns) && length(cdf_fns) >= expr_id) {
        fn <- cdf_fns[[expr_id]]
        if (is.function(fn)) {
          val <- fn(
            t = t,
            component = component,
            ctx = ctx,
            forced_complete = forced_complete,
            forced_survive = forced_survive
          )
          if (!is.null(ctx)) ctx$cdf_cond_memo[[.ctx_key("cdfId", expr_id, component, t, forced_complete, forced_survive)]] <- val
          return(val)
        }
      }
    }
  }
  0.0
}

.eval_expr_survival_cond <- function(expr, t, prep, component,
                                     forced_complete = integer(0),
                                     forced_survive = integer(0),
                                     ctx = NULL) {
  registry <- prep[[".expr_registry"]]
  if (!is.null(registry)) {
    expr_id <- .expr_registry_lookup_id(registry, expr)
    if (!is.na(expr_id)) {
      if (!is.null(ctx)) {
        key <- .ctx_key("survId", expr_id, component, t, forced_complete, forced_survive)
        cached <- ctx$surv_cond_memo[[key]]
        if (!is.null(cached)) return(cached)
      }
      surv_fns <- registry[["survival_cond_fns"]]
      if (!is.null(surv_fns) && length(surv_fns) >= expr_id) {
        fn <- surv_fns[[expr_id]]
        if (is.function(fn)) {
          val <- fn(
            t = t,
            component = component,
            ctx = ctx,
            forced_complete = forced_complete,
            forced_survive = forced_survive
          )
          if (!is.null(ctx)) ctx$surv_cond_memo[[.ctx_key("survId", expr_id, component, t, forced_complete, forced_survive)]] <- val
          return(val)
        }
      }
    }
  }
  1.0
}
