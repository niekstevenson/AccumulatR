# Likelihood system for the new expression-based API (new_API.R)
# Matches the generator logic in generator_new.R

source("R/dist.R")
source("R/utils.R")
source("R/pool_math.R")
source("R/model_tables.R")

`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.interp1_builder <- function(x, y) {
  n <- length(x)
  function(t) {
    if (!is.finite(t)) return(NA_real_)
    if (t <= x[1]) return(y[1])
    if (t >= x[n]) return(y[n])
    lo <- 1L
    hi <- n
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      if (x[mid] <= t) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    x0 <- x[lo]; x1 <- x[hi]
    y0 <- y[lo]; y1 <- y[hi]
    if (identical(x1, x0)) return(y0)
    y0 + (y1 - y0) * (t - x0) / (x1 - x0)
  }
}

.normalize_forced_sets <- function(forced_complete = character(0),
                                   forced_survive = character(0)) {
  forced_complete <- unique(as.character(forced_complete %||% character(0)))
  forced_survive <- unique(as.character(forced_survive %||% character(0)))
  forced_survive <- setdiff(forced_survive, forced_complete)
  list(complete = forced_complete, survive = forced_survive)
}

.pool_active_members <- function(prep, pool_id, component) {
  pool_defs <- prep$pools
  acc_defs <- prep$accumulators
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) return(character(0))
  members <- pool_def$members
  active_members <- character(0)
  for (m in members) {
    if (!is.null(acc_defs[[m]])) {
      acc_def <- acc_defs[[m]]
      comps <- acc_def$components
      if (length(comps) == 0 || is.null(component) ||
          identical(component, "__default__") || component %in% comps) {
        active_members <- c(active_members, m)
      }
    } else if (!is.null(pool_defs[[m]])) {
      active_members <- c(active_members, m)
    }
  }
  unique(active_members)
}

.expr_leaf_events <- function(expr) {
  if (is.null(expr)) return(character(0))
  kind <- expr$kind %||% "event"
  if (identical(kind, "event")) {
    src <- expr$source
    if (is.null(src) || is.na(src)) return(character(0))
    return(as.character(src))
  }
  if (kind %in% c("and", "or")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(character(0))
    vals <- unlist(lapply(args, .expr_leaf_events), use.names = FALSE)
    return(unique(as.character(vals)))
  }
  if (identical(kind, "guard")) {
    ref_ids <- .expr_leaf_events(expr$reference)
    block_ids <- .expr_leaf_events(expr$blocker)
    unless_list <- expr$unless %||% list()
    unless_ids <- character(0)
    if (length(unless_list) > 0) {
      unless_ids <- unlist(lapply(unless_list, .expr_leaf_events), use.names = FALSE)
    }
    return(unique(as.character(c(ref_ids, block_ids, unless_ids))))
  }
  if (identical(kind, "not")) {
    return(.expr_leaf_events(expr$arg))
  }
  character(0)
}

.collect_guard_forced_sets <- function(expr) {
  if (is.null(expr) || is.null(expr$kind) || !identical(expr$kind, "guard")) {
    return(list(complete = character(0), survive = character(0)))
  }
  ref_ids <- .expr_leaf_events(expr$reference)
  block_ids <- .expr_leaf_events(expr$blocker)
  unless_list <- expr$unless %||% list()
  unless_ids <- character(0)
  if (length(unless_list) > 0) {
    unless_ids <- unlist(lapply(unless_list, .expr_leaf_events), use.names = FALSE)
  }
  .normalize_forced_sets(ref_ids, c(block_ids, unless_ids))
}

.collect_forced_sets <- function(expr) {
  if (is.null(expr) || is.null(expr$kind)) {
    return(list(complete = character(0), survive = character(0)))
  }
  if (identical(expr$kind, "guard")) {
    return(.collect_guard_forced_sets(expr))
  }
  list(complete = character(0), survive = character(0))
}

# ==============================================================================
# Scenario bookkeeping helpers
# ==============================================================================

#' Construct a scenario object capturing conditional completions/survivals.
#'
#' @param weight Numeric weight/density contribution for the scenario.
#' @param finishers Leaf ids finishing exactly at time t in this scenario.
#' @param self_complete Leaf ids forced complete inside the focal branch.
#' @param self_survive Leaf ids forced to survive inside the focal branch.
#' @param comp_complete Leaf ids forced complete when evaluating competitors.
#' @param comp_survive Leaf ids forced to survive when evaluating competitors.
.scenario_make <- function(weight = 1.0,
                           finishers = character(0),
                           self_complete = character(0),
                           self_survive = character(0),
                           comp_complete = self_complete,
                           comp_survive = self_survive) {
  weight <- as.numeric(weight %||% 0.0)
  if (!is.finite(weight) || weight < 0) weight <- 0.0
  finishers <- unique(as.character(finishers %||% character(0)))
  self_norm <- .normalize_forced_sets(self_complete, self_survive)
  comp_norm <- .normalize_forced_sets(comp_complete, comp_survive)
  list(
    weight = weight,
    finishers = finishers,
    self_complete = self_norm$complete,
    self_survive = self_norm$survive,
    comp_complete = comp_norm$complete,
    comp_survive = comp_norm$survive
  )
}

.scenario_is_zero <- function(scn) {
  is.null(scn) || !is.list(scn) || length(scn) == 0 ||
    !is.finite(scn$weight) || scn$weight <= 0
}

.scenario_merge <- function(a, b) {
  if (.scenario_is_zero(a) || .scenario_is_zero(b)) return(NULL)
  weight <- as.numeric(a$weight) * as.numeric(b$weight)
  if (!is.finite(weight) || weight <= 0) return(NULL)
  finishers <- unique(c(a$finishers, b$finishers))
  self_complete <- union(a$self_complete, b$self_complete)
  self_survive <- union(a$self_survive, b$self_survive)
  comp_complete <- union(a$comp_complete, b$comp_complete)
  comp_survive <- union(a$comp_survive, b$comp_survive)
  .scenario_make(weight = weight,
                 finishers = finishers,
                 self_complete = self_complete,
                 self_survive = self_survive,
                 comp_complete = comp_complete,
                 comp_survive = comp_survive)
}

.scenario_cartesian <- function(scenario_lists) {
  if (length(scenario_lists) == 0) {
    return(list(.scenario_make()))
  }
  result <- list(.scenario_make())
  for (lst in scenario_lists) {
    next_result <- list()
    if (length(lst) == 0) return(list())
    for (scn_a in result) {
      for (scn_b in lst) {
        merged <- .scenario_merge(scn_a, scn_b)
        if (!is.null(merged)) {
          next_result[[length(next_result) + 1L]] <- merged
        }
      }
    }
    result <- next_result
    if (length(result) == 0) break
  }
  result
}

.scenario_key <- function(scn) {
  paste(
    paste(sort(scn$finishers), collapse = "|"),
    paste(sort(scn$self_complete), collapse = "|"),
    paste(sort(scn$self_survive), collapse = "|"),
    paste(sort(scn$comp_complete), collapse = "|"),
    paste(sort(scn$comp_survive), collapse = "|"),
    sep = ";;"
  )
}

.scenario_consolidate <- function(scenarios) {
  if (length(scenarios) == 0) return(list())
  agg <- list()
  for (scn in scenarios) {
    if (.scenario_is_zero(scn)) next
    key <- .scenario_key(scn)
    if (is.null(agg[[key]])) {
      agg[[key]] <- scn
    } else {
      agg[[key]]$weight <- agg[[key]]$weight + scn$weight
    }
  }
  Filter(function(s) !.scenario_is_zero(s), agg)
}

.scenario_apply_forced <- function(scenarios,
                                   forced_complete = character(0),
                                   forced_survive = character(0),
                                   role = c("self", "comp", "both")) {
  role <- match.arg(role)
  if (length(scenarios) == 0) return(list())
  forced_norm <- .normalize_forced_sets(forced_complete, forced_survive)
  res <- lapply(scenarios, function(scn) {
    if (.scenario_is_zero(scn)) return(NULL)
    if (role %in% c("self", "both")) {
      scn$self_complete <- union(scn$self_complete, forced_norm$complete)
      scn$self_survive <- union(scn$self_survive, forced_norm$survive)
    }
    if (role %in% c("comp", "both")) {
      scn$comp_complete <- union(scn$comp_complete, forced_norm$complete)
      scn$comp_survive <- union(scn$comp_survive, forced_norm$survive)
    }
    self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
    comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
    scn$self_complete <- self_norm$complete
    scn$self_survive <- self_norm$survive
    scn$comp_complete <- comp_norm$complete
    scn$comp_survive <- comp_norm$survive
    scn
  })
  Filter(Negate(is.null), res)
}

.scenario_baseline <- function(forced_complete = character(0),
                               forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  .scenario_make(
    weight = 1.0,
    finishers = character(0),
    self_complete = normalized$complete,
    self_survive = normalized$survive,
    comp_complete = normalized$complete,
    comp_survive = normalized$survive
  )
}

.assign_expr_ids <- function(expr, counter_env) {
  if (is.null(expr) || !is.list(expr)) return(expr)
  if (is.null(expr$`.__id`)) {
    counter_env$counter <- counter_env$counter + 1L
    expr$`.__id` <- paste0("expr#", counter_env$counter)
  }
  kind <- expr$kind %||% "event"
  if (kind %in% c("and", "or")) {
    expr$args <- lapply(expr$args %||% list(), .assign_expr_ids, counter_env = counter_env)
  } else if (identical(kind, "guard")) {
    expr$reference <- .assign_expr_ids(expr$reference, counter_env)
    expr$blocker <- .assign_expr_ids(expr$blocker, counter_env)
    unless_list <- expr$unless %||% list()
    if (length(unless_list) > 0) {
      expr$unless <- lapply(unless_list, .assign_expr_ids, counter_env = counter_env)
    }
  } else if (identical(kind, "not")) {
    expr$arg <- .assign_expr_ids(expr$arg, counter_env)
  }
  expr
}

.make_eval_key <- function(expr, component, t,
                           forced_complete = character(0),
                           forced_survive = character(0),
                           metric = "VAL") {
  expr_id <- expr$`.__id` %||% rawToChar(serialize(expr, NULL, ascii = TRUE))
  comp_key <- component %||% "__default__"
  t_key <- if (is.na(t)) "NA" else format(as.numeric(t), digits = 12, scientific = FALSE, trim = TRUE)
  fc_key <- if (length(forced_complete) == 0) "-" else paste(sort(unique(forced_complete)), collapse = ",")
  fs_key <- if (length(forced_survive) == 0) "-" else paste(sort(unique(forced_survive)), collapse = ",")
  paste(metric, expr_id, comp_key, t_key, fc_key, fs_key, sep = "|")
}

.expr_id <- function(expr) {
  expr$`.__id` %||% rawToChar(serialize(expr, NULL, ascii = TRUE))
}

.make_integral_key <- function(expr, component,
                               forced_complete = character(0),
                               forced_survive = character(0),
                               forced_complete_self = forced_complete,
                               forced_survive_self = forced_survive,
                               competitor_exprs = list(),
                               metric = "INT") {
  expr_key <- .expr_id(expr)
  comp_key <- component %||% "__default__"
  fc_base <- if (length(forced_complete) == 0) "-" else paste(sort(unique(forced_complete)), collapse = ",")
  fs_base <- if (length(forced_survive) == 0) "-" else paste(sort(unique(forced_survive)), collapse = ",")
  fc_self <- if (length(forced_complete_self) == 0) "-" else paste(sort(unique(forced_complete_self)), collapse = ",")
  fs_self <- if (length(forced_survive_self) == 0) "-" else paste(sort(unique(forced_survive_self)), collapse = ",")
  comp_ids <- if (length(competitor_exprs) == 0) {
    "-"
  } else {
    ids <- vapply(competitor_exprs, .expr_id, character(1))
    paste(sort(unique(ids)), collapse = ",")
  }
  paste(metric, expr_key, comp_key, fc_base, fs_base, fc_self, fs_self, comp_ids, sep = "|")
}

._cached_integral <- function(prep, key, integrand_scalar, upper,
                              rel.tol = 1e-5, abs.tol = 1e-6) {
  if (!is.finite(upper) || upper <= 0) return(0.0)
  integral_cache <- prep$.__integral_cache
  if (is.null(integral_cache)) {
    integral_cache <- new.env(parent = emptyenv())
    prep$.__integral_cache <- integral_cache
  }
  entry <- integral_cache[[key]]
  integrand_vec <- function(u) vapply(u, integrand_scalar, numeric(1))
  if (is.null(entry)) {
    entry <- list(
      times = numeric(0),
      values = numeric(0),
      tol = 1e-9,
      integrand_vec = integrand_vec
    )
  } else if (is.null(entry$integrand_vec)) {
    entry$integrand_vec <- integrand_vec
  }
  times <- entry$times
  values <- entry$values
  tol <- entry$tol %||% 1e-9

  match_idx <- which(abs(times - upper) <= tol)
  if (length(match_idx) > 0) {
    return(values[[match_idx[[1]]]])
  }

  lower_t <- 0.0
  lower_val <- 0.0
  if (length(times) > 0) {
    idx_less <- which(times < (upper - tol))
    if (length(idx_less) > 0) {
      idx <- idx_less[[which.max(times[idx_less])]]
      lower_t <- times[[idx]]
      lower_val <- values[[idx]]
    }
  }

  val <- lower_val
  if (abs(upper - lower_t) > tol) {
    delta <- tryCatch(
      stats::integrate(
        entry$integrand_vec,
        lower = lower_t,
        upper = upper,
        rel.tol = rel.tol,
        abs.tol = abs.tol,
        stop.on.error = FALSE
      )$value,
      error = function(e) 0.0
    )
    val <- lower_val + as.numeric(delta)
  }

  if (!is.finite(val) || val < 0) val <- 0.0

  if (length(times) == 0 || all(abs(times - upper) > tol)) {
    entry$times <- c(times, upper)
    entry$values <- c(values, val)
    if (length(entry$times) > 1) {
      ord <- order(entry$times)
      entry$times <- entry$times[ord]
      entry$values <- entry$values[ord]
    }
  } else {
    idx <- which.min(abs(entry$times - upper))
    entry$times[[idx]] <- upper
    entry$values[[idx]] <- val
  }

  integral_cache[[key]] <- entry
  val
}

.combination_list <- function(items, size) {
  items <- unique(as.character(items %||% character(0)))
  if (size < 0) return(list())
  if (size == 0) return(list(character(0)))
  if (length(items) < size) return(list())
  utils::combn(items, size, simplify = FALSE)
}

.scenario_extend_with_member <- function(base_scenarios, member_id, mode,
                                         prep, component, t) {
  if (length(base_scenarios) == 0) return(list())
  mode <- match.arg(mode, c("finish", "complete", "survive"))
  acc_defs <- prep$accumulators
  pool_defs <- prep$pools
  next_scenarios <- list()
  for (scn in base_scenarios) {
    if (.scenario_is_zero(scn)) next
    forced_complete <- scn$self_complete
    forced_survive <- scn$self_survive
    candidates <- list()
    if (!is.null(acc_defs[[member_id]])) {
      acc_def <- acc_defs[[member_id]]
      if (mode == "finish") {
        candidates <- .acc_finish_scenarios(member_id, acc_def, t,
                                            forced_complete = forced_complete,
                                            forced_survive = forced_survive)
      } else if (mode == "complete") {
        candidates <- .acc_complete_scenarios(member_id, acc_def, t,
                                              forced_complete = forced_complete,
                                              forced_survive = forced_survive)
      } else {
        candidates <- .acc_survive_scenarios(member_id, acc_def, t,
                                             forced_complete = forced_complete,
                                             forced_survive = forced_survive)
      }
    } else if (!is.null(pool_defs[[member_id]])) {
      if (mode == "finish") {
        candidates <- .pool_density_scenarios(prep, member_id, component, t,
                                              forced_complete = forced_complete,
                                              forced_survive = forced_survive)
      } else if (mode == "complete") {
        candidates <- .pool_complete_scenarios(prep, member_id, component, t,
                                               forced_complete = forced_complete,
                                               forced_survive = forced_survive)
      } else {
        candidates <- .pool_survival_scenarios(prep, member_id, component, t,
                                               forced_complete = forced_complete,
                                               forced_survive = forced_survive)
      }
    } else {
      next
    }
    if (length(candidates) == 0) next
    for (cand in candidates) {
      merged <- .scenario_merge(scn, cand)
      if (!is.null(merged)) {
        next_scenarios[[length(next_scenarios) + 1L]] <- merged
      }
    }
  }
  .scenario_consolidate(next_scenarios)
}

.scenario_extend_with_expr <- function(base_scenarios, expr, mode,
                                       prep, component, t) {
  if (length(base_scenarios) == 0) return(list())
  mode <- match.arg(mode, c("finish", "complete", "survive"))
  next_scenarios <- list()
  for (scn in base_scenarios) {
    if (.scenario_is_zero(scn)) next
    forced_complete <- scn$self_complete
    forced_survive <- scn$self_survive
    candidates <- switch(
      mode,
      finish = .expr_finish_scenarios(expr, t, prep, component,
                                      forced_complete = forced_complete,
                                      forced_survive = forced_survive),
      complete = .expr_complete_scenarios(expr, t, prep, component,
                                          forced_complete = forced_complete,
                                          forced_survive = forced_survive),
      survive = .expr_survive_scenarios(expr, t, prep, component,
                                        forced_complete = forced_complete,
                                        forced_survive = forced_survive)
    )
    if (length(candidates) == 0) next
    for (cand in candidates) {
      merged <- .scenario_merge(scn, cand)
      if (!is.null(merged)) {
        next_scenarios[[length(next_scenarios) + 1L]] <- merged
      }
    }
  }
  .scenario_consolidate(next_scenarios)
}

.acc_finish_scenarios <- function(acc_id, acc_def, t,
                                  forced_complete = character(0),
                                  forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  if (acc_id %in% normalized$complete) return(list())
  if (acc_id %in% normalized$survive) return(list())
  dens <- .acc_density(acc_def, t)
  if (!is.finite(dens) || dens <= 0) return(list())
  list(.scenario_make(
    weight = dens,
    finishers = acc_id,
    self_complete = union(normalized$complete, acc_id),
    self_survive = normalized$survive,
    comp_complete = union(normalized$complete, acc_id),
    comp_survive = normalized$survive
  ))
}

.acc_complete_scenarios <- function(acc_id, acc_def, t,
                                    forced_complete = character(0),
                                    forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  if (acc_id %in% normalized$survive) return(list())
  if (acc_id %in% normalized$complete) {
    return(list(.scenario_make(
      weight = 1.0,
      finishers = character(0),
      self_complete = normalized$complete,
      self_survive = normalized$survive,
      comp_complete = normalized$complete,
      comp_survive = normalized$survive
    )))
  }
  cdf_val <- 1 - .acc_survival(acc_def, t)
  if (!is.finite(cdf_val) || cdf_val <= 0) return(list())
  list(.scenario_make(
    weight = cdf_val,
    finishers = character(0),
    self_complete = union(normalized$complete, acc_id),
    self_survive = normalized$survive,
    comp_complete = union(normalized$complete, acc_id),
    comp_survive = normalized$survive
  ))
}

.acc_survive_scenarios <- function(acc_id, acc_def, t,
                                   forced_complete = character(0),
                                   forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  if (acc_id %in% normalized$complete) return(list())
  if (acc_id %in% normalized$survive) {
    return(list(.scenario_make(
      weight = 1.0,
      finishers = character(0),
      self_complete = normalized$complete,
      self_survive = normalized$survive,
      comp_complete = normalized$complete,
      comp_survive = normalized$survive
    )))
  }
  surv_val <- .acc_survival(acc_def, t)
  if (!is.finite(surv_val) || surv_val <= 0) return(list())
  list(.scenario_make(
    weight = surv_val,
    finishers = character(0),
    self_complete = normalized$complete,
    self_survive = union(normalized$survive, acc_id),
    comp_complete = normalized$complete,
    comp_survive = union(normalized$survive, acc_id)
  ))
}

.pool_density_scenarios <- function(prep, pool_id, component, t,
                                    forced_complete = character(0),
                                    forced_survive = character(0)) {
  if (!is.finite(t) || t < 0) return(list())
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive

  pool_defs <- prep$pools
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) {
    stop(sprintf("Unknown pool '%s'", pool_id))
  }

  active_members <- .pool_active_members(prep, pool_id, component)
  if (length(active_members) == 0) return(list())

  k <- pool_def$k %||% 1L
  forced_complete_members <- intersect(active_members, forced_complete)
  forced_survive_members <- setdiff(intersect(active_members, forced_survive),
                                    forced_complete_members)
  free_members <- setdiff(active_members, union(forced_complete_members,
                                                forced_survive_members))
  effective_k <- k - length(forced_complete_members)

  if (effective_k <= 0) return(list())
  if (length(free_members) == 0) return(list())
  if (effective_k > length(free_members)) return(list())

  base <- list(.scenario_baseline(forced_complete, forced_survive))
  scenario_accum <- list()

  for (cand in free_members) {
    current <- .scenario_extend_with_member(base, cand, "finish",
                                            prep, component, t)
    if (length(current) == 0) next

    required_pre <- effective_k - 1L
    other_free <- setdiff(free_members, cand)
    combos <- .combination_list(other_free, required_pre)
    if (length(combos) == 0 && required_pre > 0) next
    if (length(combos) == 0) combos <- list(character(0))

    for (combo in combos) {
      scn_combo <- current
      if (length(combo) > 0) {
        for (mid in combo) {
          scn_combo <- .scenario_extend_with_member(scn_combo, mid, "complete",
                                                    prep, component, t)
          if (length(scn_combo) == 0) break
        }
      }
      if (length(scn_combo) == 0) next

      survivors <- setdiff(active_members,
                           union(union(forced_complete_members, combo), cand))
      if (length(survivors) > 0) {
        for (mid in survivors) {
          scn_combo <- .scenario_extend_with_member(scn_combo, mid, "survive",
                                                    prep, component, t)
          if (length(scn_combo) == 0) break
        }
      }
      if (length(scn_combo) == 0) next

      scn_combo <- lapply(scn_combo, function(scn) {
        scn$finishers <- unique(c(scn$finishers, pool_id))
        scn$self_complete <- union(scn$self_complete, pool_id)
        scn$comp_complete <- union(scn$comp_complete, pool_id)
        self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
        comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
        scn$self_complete <- self_norm$complete
        scn$self_survive <- self_norm$survive
        scn$comp_complete <- comp_norm$complete
        scn$comp_survive <- comp_norm$survive
        scn
      })
      scenario_accum <- c(scenario_accum, scn_combo)
    }
  }

  .scenario_consolidate(scenario_accum)
}

.pool_complete_scenarios <- function(prep, pool_id, component, t,
                                     forced_complete = character(0),
                                     forced_survive = character(0)) {
  if (!is.finite(t) || t < 0) return(list())
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive

  pool_defs <- prep$pools
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) {
    stop(sprintf("Unknown pool '%s'", pool_id))
  }

  active_members <- .pool_active_members(prep, pool_id, component)
  if (length(active_members) == 0) return(list())

  k <- pool_def$k %||% 1L
  forced_complete_members <- intersect(active_members, forced_complete)
  forced_survive_members <- setdiff(intersect(active_members, forced_survive),
                                    forced_complete_members)
  free_members <- setdiff(active_members, union(forced_complete_members,
                                                forced_survive_members))
  effective_k <- k - length(forced_complete_members)

  if (effective_k <= 0) {
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$self_complete <- union(scn$self_complete, pool_id)
    scn$comp_complete <- union(scn$comp_complete, pool_id)
    scn_norm_self <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
    scn_norm_comp <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
    scn$self_complete <- scn_norm_self$complete
    scn$self_survive <- scn_norm_self$survive
    scn$comp_complete <- scn_norm_comp$complete
    scn$comp_survive <- scn_norm_comp$survive
    scn$weight <- 1.0
    return(list(scn))
  }

  if (length(free_members) == 0) return(list())
  if (effective_k > length(free_members)) return(list())

  base <- list(.scenario_baseline(forced_complete, forced_survive))
  scenario_accum <- list()

  for (m in seq.int(effective_k, length(free_members))) {
    combos <- .combination_list(free_members, m)
    if (length(combos) == 0) next
    for (combo in combos) {
      scn_combo <- base
      for (mid in combo) {
        scn_combo <- .scenario_extend_with_member(scn_combo, mid, "complete",
                                                  prep, component, t)
        if (length(scn_combo) == 0) break
      }
      if (length(scn_combo) == 0) next

      survivors <- setdiff(active_members,
                           union(forced_complete_members, combo))
      if (length(survivors) > 0) {
        for (mid in survivors) {
          scn_combo <- .scenario_extend_with_member(scn_combo, mid, "survive",
                                                    prep, component, t)
          if (length(scn_combo) == 0) break
        }
      }
      if (length(scn_combo) == 0) next

      scn_combo <- lapply(scn_combo, function(scn) {
        scn$self_complete <- union(scn$self_complete, pool_id)
        scn$comp_complete <- union(scn$comp_complete, pool_id)
        self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
        comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
        scn$self_complete <- self_norm$complete
        scn$self_survive <- self_norm$survive
        scn$comp_complete <- comp_norm$complete
        scn$comp_survive <- comp_norm$survive
        scn
      })
      scenario_accum <- c(scenario_accum, scn_combo)
    }
  }

  .scenario_consolidate(scenario_accum)
}

.pool_survival_scenarios <- function(prep, pool_id, component, t,
                                     forced_complete = character(0),
                                     forced_survive = character(0)) {
  if (!is.finite(t) || t < 0) {
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$self_survive <- union(scn$self_survive, pool_id)
    scn$comp_survive <- union(scn$comp_survive, pool_id)
    self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
    comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
    scn$self_complete <- self_norm$complete
    scn$self_survive <- self_norm$survive
    scn$comp_complete <- comp_norm$complete
    scn$comp_survive <- comp_norm$survive
    scn$weight <- 1.0
    return(list(scn))
  }

  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive

  pool_defs <- prep$pools
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) {
    stop(sprintf("Unknown pool '%s'", pool_id))
  }

  active_members <- .pool_active_members(prep, pool_id, component)
  if (length(active_members) == 0) {
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$self_survive <- union(scn$self_survive, pool_id)
    scn$comp_survive <- union(scn$comp_survive, pool_id)
    self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
    comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
    scn$self_complete <- self_norm$complete
    scn$self_survive <- self_norm$survive
    scn$comp_complete <- comp_norm$complete
    scn$comp_survive <- comp_norm$survive
    scn$weight <- 1.0
    return(list(scn))
  }

  k <- pool_def$k %||% 1L
  forced_complete_members <- intersect(active_members, forced_complete)
  forced_survive_members <- setdiff(intersect(active_members, forced_survive),
                                    forced_complete_members)
  free_members <- setdiff(active_members, union(forced_complete_members,
                                                forced_survive_members))
  effective_k <- k - length(forced_complete_members)

  if (effective_k <= 0) return(list())
  if (effective_k > length(free_members)) {
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$self_survive <- union(scn$self_survive, pool_id)
    scn$comp_survive <- union(scn$comp_survive, pool_id)
    self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
    comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
    scn$self_complete <- self_norm$complete
    scn$self_survive <- self_norm$survive
    scn$comp_complete <- comp_norm$complete
    scn$comp_survive <- comp_norm$survive
    scn$weight <- 1.0
    return(list(scn))
  }

  base <- list(.scenario_baseline(forced_complete, forced_survive))
  scenario_accum <- list()
  max_finish <- min(effective_k - 1L, length(free_members))

  for (m in 0:max_finish) {
    combos <- .combination_list(free_members, m)
    if (length(combos) == 0) next
    for (combo in combos) {
      scn_combo <- base
      if (length(combo) > 0) {
        for (mid in combo) {
          scn_combo <- .scenario_extend_with_member(scn_combo, mid, "complete",
                                                    prep, component, t)
          if (length(scn_combo) == 0) break
        }
      }
      if (length(scn_combo) == 0) next

      survivors <- setdiff(active_members,
                           union(forced_complete_members, combo))
      if (length(survivors) > 0) {
        for (mid in survivors) {
          scn_combo <- .scenario_extend_with_member(scn_combo, mid, "survive",
                                                    prep, component, t)
          if (length(scn_combo) == 0) break
        }
      }
      if (length(scn_combo) == 0) next

      scn_combo <- lapply(scn_combo, function(scn) {
        scn$self_survive <- union(scn$self_survive, pool_id)
        scn$comp_survive <- union(scn$comp_survive, pool_id)
        self_norm <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
        comp_norm <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
        scn$self_complete <- self_norm$complete
        scn$self_survive <- self_norm$survive
        scn$comp_complete <- comp_norm$complete
        scn$comp_survive <- comp_norm$survive
        scn
      })
      scenario_accum <- c(scenario_accum, scn_combo)
    }
  }

  .scenario_consolidate(scenario_accum)
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

  counter_env <- new.env(parent = emptyenv())
  counter_env$counter <- 0L
  prep$outcomes <- lapply(prep$outcomes, function(out_def) {
    out_def$expr <- .assign_expr_ids(out_def$expr, counter_env)
    out_def
  })
  prep$`.__expr_id_counter` <- counter_env$counter
  prep$`.__guard_cache` <- new.env(parent = emptyenv())
  prep$`.__memo` <- new.env(parent = emptyenv())
  prep$`.__integral_cache` <- new.env(parent = emptyenv())

  prep
}

# ==============================================================================
# PART 2: Distribution Primitives
# ==============================================================================

# Get density function for an accumulator at time t
.acc_density <- function(acc_def, t) {
  if (!is.finite(t) || t < 0) return(0.0)
  if (t < acc_def$onset) return(0.0)
  
  # Account for q (probability of non-completion)
  success_prob <- 1 - acc_def$q
  if (success_prob <= 0) return(0.0)
  
  reg <- dist_registry(acc_def$dist)
  if (is.null(reg) || is.null(reg$d)) {
    stop(sprintf("No density function for distribution '%s'", acc_def$dist))
  }
  
  # Density of the underlying distribution at (t - onset)
  t_adj <- t - acc_def$onset
  dens <- reg$d(t_adj, acc_def$params)
  
  # Scale by success probability
  success_prob * dens
}

# Get survival function for an accumulator at time t
.acc_survival <- function(acc_def, t) {
  if (!is.finite(t)) return(0.0)
  if (t < 0) return(1.0)
  if (t < acc_def$onset) return(1.0)
  
  # Survival = q + (1-q) * S(t - onset)
  reg <- dist_registry(acc_def$dist)
  if (is.null(reg) || is.null(reg$p)) {
    stop(sprintf("No CDF function for distribution '%s'", acc_def$dist))
  }
  
  t_adj <- t - acc_def$onset
  surv_underlying <- 1 - reg$p(t_adj, acc_def$params)
  
  acc_def$q + (1 - acc_def$q) * surv_underlying
}

# ==============================================================================
# PART 3: Pool Likelihood
# ==============================================================================

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = character(0),
                          forced_survive = character(0)) {
  scenarios <- .pool_density_scenarios(
    prep, pool_id, component, t,
    forced_complete = forced_complete,
    forced_survive = forced_survive
  )
  if (length(scenarios) == 0) return(0.0)
  sum(vapply(scenarios, `[[`, numeric(1), "weight"))
}

.pool_survival <- function(prep, pool_id, component, t,
                           forced_complete = character(0),
                           forced_survive = character(0)) {
  scenarios <- .pool_survival_scenarios(
    prep, pool_id, component, t,
    forced_complete = forced_complete,
    forced_survive = forced_survive
  )
  if (length(scenarios) == 0) return(0.0)
  sum(vapply(scenarios, `[[`, numeric(1), "weight"))
}

# ==============================================================================
# PART 4: Expression Evaluation (Likelihood)
# ==============================================================================

.expr_finish_scenarios <- function(expr, t, prep, component,
                                   forced_complete = character(0),
                                   forced_survive = character(0)) {
  if (is.null(expr)) return(list())
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive
  kind <- expr$kind %||% "event"

  memo_env <- prep$.__memo
  cache_key <- NULL
  memo_store <- function(val) val
  if (!is.null(memo_env) && is.environment(memo_env)) {
    cache_key <- .make_eval_key(expr, component, t,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                metric = "FIN")
    cached <- memo_env[[cache_key]]
    if (!is.null(cached)) return(cached)
    memo_store <- function(val) {
      memo_env[[cache_key]] <- val
      val
    }
  }

  if (identical(kind, "event")) {
    source_id <- expr$source
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(memo_store(list()))
    }
    if (!is.null(prep$pools[[source_id]])) {
      return(memo_store(.pool_density_scenarios(prep, source_id, component, t,
                                                forced_complete = forced_complete,
                                                forced_survive = forced_survive)))
    }
    acc_def <- prep$accumulators[[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def$components
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(memo_store(list()))
      }
      return(memo_store(.acc_finish_scenarios(source_id, acc_def, t,
                                              forced_complete = forced_complete,
                                              forced_survive = forced_survive)))
    }
    return(memo_store(list()))
  }

  if (identical(kind, "and")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(memo_store(list()))
    scenario_accum <- list()
    for (i in seq_along(args)) {
      fin_i <- .expr_finish_scenarios(args[[i]], t, prep, component,
                                      forced_complete = forced_complete,
                                      forced_survive = forced_survive)
      if (length(fin_i) == 0) next
      scn_i <- fin_i
      for (j in seq_along(args)) {
        if (i == j) next
        scn_i <- .scenario_extend_with_expr(scn_i, args[[j]], "complete",
                                            prep, component, t)
        if (length(scn_i) == 0) break
      }
      scenario_accum <- c(scenario_accum, scn_i)
    }
    return(memo_store(.scenario_consolidate(scenario_accum)))
  }

  if (identical(kind, "or")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(memo_store(list()))
    scenario_accum <- list()
    for (i in seq_along(args)) {
      fin_i <- .expr_finish_scenarios(args[[i]], t, prep, component,
                                      forced_complete = forced_complete,
                                      forced_survive = forced_survive)
      if (length(fin_i) == 0) next
      scn_i <- fin_i
      for (j in seq_along(args)) {
        if (i == j) next
        scn_i <- .scenario_extend_with_expr(scn_i, args[[j]], "survive",
                                            prep, component, t)
        if (length(scn_i) == 0) break
      }
      scenario_accum <- c(scenario_accum, scn_i)
    }
    return(memo_store(.scenario_consolidate(scenario_accum)))
  }

  if (identical(kind, "guard")) {
    reference <- expr$reference
    blocker <- expr$blocker
    unless_list <- expr$unless %||% list()
    protect_leaves <- if (length(unless_list) > 0) unique(unlist(lapply(unless_list, .expr_leaf_events), use.names = FALSE)) else character(0)
    blocker_leaves <- .expr_leaf_events(blocker)

    ref_scenarios <- .expr_finish_scenarios(reference, t, prep, component,
                                            forced_complete = forced_complete,
                                            forced_survive = forced_survive)
    if (length(ref_scenarios) == 0) return(memo_store(list()))

    scenario_accum <- list()

    guard_success_weight <- function(scn) {
      forced_self_complete <- scn$self_complete
      forced_self_survive <- scn$self_survive

      if (length(unless_list) == 0) {
        return(.eval_expr_survival(
          blocker, t, prep, component,
          forced_complete = forced_self_complete,
          forced_survive = forced_self_survive
        ))
      }

      if (length(protect_leaves) > 0 &&
          any(protect_leaves %in% forced_self_complete)) {
        return(1.0)
      }

      guard_cache <- prep$.__guard_cache
      if (is.null(guard_cache)) {
        guard_cache <- new.env(parent = emptyenv())
        prep$.__guard_cache <- guard_cache
      }
      guard_id <- expr$`.__id` %||% rawToChar(serialize(expr, NULL, ascii = TRUE))
      cache_key <- paste(
        guard_id,
        paste(sort(unique(forced_self_complete)), collapse = ","),
        paste(sort(unique(forced_self_survive)), collapse = ","),
        sep = "|"
      )
      entry <- guard_cache[[cache_key]]
      if (is.null(entry)) {
        integrand_scalar <- function(u) {
          if (!is.finite(u) || u < 0) return(0.0)
          blocker_scenarios <- .expr_finish_scenarios(
            blocker, u, prep, component,
            forced_complete = forced_self_complete,
            forced_survive = forced_self_survive
          )
          if (length(blocker_scenarios) == 0) return(0.0)
          total <- 0.0
          for (bs in blocker_scenarios) {
            dens <- bs$weight
            if (!is.finite(dens) || dens <= 0) next
            forced_complete_current <- union(forced_self_complete, bs$self_complete)
            forced_survive_current <- union(forced_self_survive, bs$self_survive)
            surv_prod <- 1.0
            for (prot in unless_list) {
              surv_val <- .eval_expr_survival(
                prot, u, prep, component,
                forced_complete = forced_complete_current,
                forced_survive = forced_survive_current
              )
              surv_prod <- surv_prod * surv_val
              if (surv_prod == 0) break
            }
            total <- total + dens * surv_prod
          }
          total
        }
        entry <- list(
          cache = new.env(parent = emptyenv()),
          integrand = integrand_scalar,
          integrand_vec = function(u) vapply(u, integrand_scalar, numeric(1)),
          # Track cumulative guard integrals so repeated calls only integrate forward
          times = numeric(0),
          values = numeric(0),
          tol = 1e-10
        )
        guard_cache[[cache_key]] <- entry
      }

      calc_fail <- function(tt) {
        if (!is.finite(tt) || tt <= 0) return(0.0)
        key_t <- format(as.numeric(tt), digits = 12, scientific = FALSE, trim = TRUE)
        cached <- entry$cache[[key_t]]
        if (!is.null(cached)) return(cached)

        tol <- entry$tol %||% 1e-10
        times <- entry$times
        values <- entry$values

        if (length(times) > 0) {
          idx_exact <- which(abs(times - tt) <= tol)
          if (length(idx_exact) > 0) {
            val_exact <- values[[idx_exact[[1]]]]
            entry$cache[[key_t]] <- val_exact
            guard_cache[[cache_key]] <- entry
            return(val_exact)
          }
          idx_less <- which(times < (tt - tol))
          if (length(idx_less) > 0) {
            base_idx <- idx_less[[which.max(times[idx_less])]]
            lower_t <- times[[base_idx]]
            lower_val <- values[[base_idx]]
          } else {
            lower_t <- 0.0
            lower_val <- 0.0
          }
        } else {
          lower_t <- 0.0
          lower_val <- 0.0
        }

        if (abs(tt - lower_t) <= tol) {
          new_val <- lower_val
        } else {
          delta <- tryCatch(
            stats::integrate(
              entry$integrand_vec,
              lower = lower_t,
              upper = tt,
              rel.tol = 1e-5,
              abs.tol = 1e-6,
              stop.on.error = FALSE
            )$value,
            error = function(e) 0.0
          )
          delta <- as.numeric(delta)
          if (!is.finite(delta)) delta <- 0.0
          new_val <- lower_val + delta
        }

        new_val <- min(max(as.numeric(new_val), 0.0), 1.0)
        entry$cache[[key_t]] <- new_val
        if (length(times) == 0 || all(abs(times - tt) > tol)) {
          entry$times <- c(times, tt)
          entry$values <- c(values, new_val)
          if (length(entry$times) > 1) {
            ord <- order(entry$times)
            entry$times <- entry$times[ord]
            entry$values <- entry$values[ord]
          }
        } else if (length(times) > 0) {
          idx_match <- which.min(abs(times - tt))
          entry$times[[idx_match]] <- times[[idx_match]]
          entry$values[[idx_match]] <- new_val
        }
        guard_cache[[cache_key]] <- entry
        new_val
      }

      fail_prob <- calc_fail(t)
      S_eff <- 1.0 - fail_prob
      max(0.0, min(1.0, S_eff))
    }

    for (scn in ref_scenarios) {
      S_eff <- guard_success_weight(scn)
      if (S_eff <= 0) next

      scn_new <- scn
      scn_new$weight <- scn$weight * S_eff
      scenario_accum[[length(scenario_accum) + 1L]] <- scn_new
    }

    return(memo_store(.scenario_consolidate(scenario_accum)))
  }

  if (identical(kind, "not")) {
    return(memo_store(list()))
  }

  memo_store(list())
}

.expr_complete_scenarios <- function(expr, t, prep, component,
                                     forced_complete = character(0),
                                     forced_survive = character(0)) {
  if (is.null(expr)) return(list())
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive
  kind <- expr$kind %||% "event"

  memo_env <- prep$.__memo
  cache_key <- NULL
  memo_store <- function(val) val
  if (!is.null(memo_env) && is.environment(memo_env)) {
    cache_key <- .make_eval_key(expr, component, t,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                metric = "COM")
    cached <- memo_env[[cache_key]]
    if (!is.null(cached)) return(cached)
    memo_store <- function(val) {
      memo_env[[cache_key]] <- val
      val
    }
  }

  if (identical(kind, "event")) {
    source_id <- expr$source
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      scn <- .scenario_baseline(forced_complete, forced_survive)
      scn$self_complete <- union(scn$self_complete, source_id)
      scn$comp_complete <- union(scn$comp_complete, source_id)
      norm_self <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
      norm_comp <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
      scn$self_complete <- norm_self$complete
      scn$self_survive <- norm_self$survive
      scn$comp_complete <- norm_comp$complete
      scn$comp_survive <- norm_comp$survive
      return(memo_store(list(scn)))
    }
    if (!is.null(prep$pools[[source_id]])) {
      return(memo_store(.pool_complete_scenarios(prep, source_id, component, t,
                                                 forced_complete = forced_complete,
                                                 forced_survive = forced_survive)))
    }
    acc_def <- prep$accumulators[[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def$components
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        scn <- .scenario_baseline(forced_complete, forced_survive)
        scn$self_complete <- union(scn$self_complete, source_id)
        scn$comp_complete <- union(scn$comp_complete, source_id)
        norm_self <- .normalize_forced_sets(scn$self_complete, scn$self_survive)
        norm_comp <- .normalize_forced_sets(scn$comp_complete, scn$comp_survive)
        scn$self_complete <- norm_self$complete
        scn$self_survive <- norm_self$survive
        scn$comp_complete <- norm_comp$complete
        scn$comp_survive <- norm_comp$survive
        return(memo_store(list(scn)))
      }
      return(memo_store(.acc_complete_scenarios(source_id, acc_def, t,
                                                forced_complete = forced_complete,
                                                forced_survive = forced_survive)))
    }
    return(memo_store(list(.scenario_baseline(forced_complete, forced_survive))))
  }

  if (identical(kind, "and")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(memo_store(list(.scenario_baseline(forced_complete, forced_survive))))
    scenario_lists <- lapply(args, function(arg) {
      .expr_complete_scenarios(arg, t, prep, component,
                               forced_complete = forced_complete,
                               forced_survive = forced_survive)
    })
    return(memo_store(.scenario_cartesian(scenario_lists)))
  }

  if (identical(kind, "or")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(memo_store(list()))
    scenario_accum <- list()
    for (arg in args) {
      arg_comp <- .expr_complete_scenarios(arg, t, prep, component,
                                           forced_complete = forced_complete,
                                           forced_survive = forced_survive)
      scenario_accum <- c(scenario_accum, arg_comp)
    }
    return(memo_store(.scenario_consolidate(scenario_accum)))
  }

  if (identical(kind, "guard")) {
    if (!is.finite(t) || t <= 0) return(memo_store(list()))
    forced_complete_self <- forced_complete
    forced_survive_self <- forced_survive

    probe_times <- unique(Filter(function(x) is.finite(x) && x > 0,
                                 c(t, t / 2, t / 4, 1e-3, t + 1e-3)))
    if (length(probe_times) == 0) probe_times <- c(t)

    template_env <- new.env(parent = emptyenv())
    keys <- character(0)
    for (pt in probe_times) {
      scns_pt <- .expr_finish_scenarios(expr, pt, prep, component,
                                        forced_complete = forced_complete_self,
                                        forced_survive = forced_survive_self)
      if (length(scns_pt) == 0) next
      for (scn in scns_pt) {
        key <- .scenario_key(scn)
        if (!(key %in% keys)) keys <- c(keys, key)
        if (is.null(template_env[[key]])) template_env[[key]] <- scn
      }
    }

    if (length(keys) == 0) return(memo_store(list()))

    results <- list()
    for (key in keys) {
      integrand <- (function(target_key) {
        function(u) {
          vapply(u, function(ui) {
            if (!is.finite(ui) || ui < 0) return(0.0)
            scns_ui <- .expr_finish_scenarios(expr, ui, prep, component,
                                              forced_complete = forced_complete_self,
                                              forced_survive = forced_survive_self)
            if (length(scns_ui) == 0) return(0.0)
            total <- 0.0
            for (scn in scns_ui) {
              scn_key <- .scenario_key(scn)
              if (is.null(template_env[[scn_key]])) template_env[[scn_key]] <- scn
              if (identical(scn_key, target_key)) {
                total <- total + scn$weight
              }
            }
            total
          }, numeric(1))
        }
      })(key)

      weight <- tryCatch(
        stats::integrate(integrand, lower = 0, upper = t,
                         rel.tol = 1e-5, abs.tol = 1e-6,
                         stop.on.error = FALSE)$value,
        error = function(e) 0.0
      )
      weight <- as.numeric(weight)
      if (!is.finite(weight) || weight <= 0) next
      scn_template <- template_env[[key]]
      if (is.null(scn_template)) next
      scn_template$weight <- weight
      results[[length(results) + 1L]] <- scn_template
    }

    return(memo_store(.scenario_consolidate(results)))
  }

  if (identical(kind, "not")) {
    return(memo_store(list()))
  }

  memo_store(list())
}

.expr_survive_scenarios <- function(expr, t, prep, component,
                                    forced_complete = character(0),
                                    forced_survive = character(0)) {
  if (is.null(expr)) return(list())
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive
  kind <- expr$kind %||% "event"

  memo_env <- prep$.__memo
  cache_key <- NULL
  memo_store <- function(val) val
  if (!is.null(memo_env) && is.environment(memo_env)) {
    cache_key <- .make_eval_key(expr, component, t,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                metric = "SURVSCN")
    cached <- memo_env[[cache_key]]
    if (!is.null(cached)) return(cached)
    memo_store <- function(val) {
      memo_env[[cache_key]] <- val
      val
    }
  }

  if (identical(kind, "event")) {
    source_id <- expr$source
    if (identical(source_id, "__DEADLINE__") || identical(source_id, "__GUESS__")) {
      return(memo_store(list()))
    }
    if (!is.null(prep$pools[[source_id]])) {
      return(memo_store(.pool_survival_scenarios(prep, source_id, component, t,
                                                 forced_complete = forced_complete,
                                                 forced_survive = forced_survive)))
    }
    acc_def <- prep$accumulators[[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def$components
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(memo_store(list()))
      }
      return(memo_store(.acc_survive_scenarios(source_id, acc_def, t,
                                               forced_complete = forced_complete,
                                               forced_survive = forced_survive)))
    }
    return(memo_store(list()))
  }

  if (identical(kind, "and")) {
    complete_scenarios <- .expr_complete_scenarios(expr, t, prep, component,
                                                   forced_complete = forced_complete,
                                                   forced_survive = forced_survive)
    prob_complete <- sum(vapply(complete_scenarios, `[[`, numeric(1), "weight"))
    prob_complete <- min(max(prob_complete, 0.0), 1.0)
    surv_weight <- 1.0 - prob_complete
    if (surv_weight <= 0) return(memo_store(list()))
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$weight <- surv_weight
    return(memo_store(list(scn)))
  }

  if (identical(kind, "or")) {
    args <- expr$args %||% list()
    if (length(args) == 0) return(memo_store(list(.scenario_baseline(forced_complete, forced_survive))))
    scenario_lists <- lapply(args, function(arg) {
      .expr_survive_scenarios(arg, t, prep, component,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive)
    })
    return(memo_store(.scenario_cartesian(scenario_lists)))
  }

  if (identical(kind, "guard")) {
    if (!is.finite(t) || t <= 0) {
      scn <- .scenario_baseline(forced_complete, forced_survive)
      scn$weight <- 1.0
      return(memo_store(list(scn)))
    }
    complete_scenarios <- .expr_complete_scenarios(expr, t, prep, component,
                                                  forced_complete = forced_complete,
                                                  forced_survive = forced_survive)
    prob_finish <- sum(vapply(complete_scenarios, `[[`, numeric(1), "weight"))
    prob_finish <- min(max(prob_finish, 0.0), 1.0)
    surv_weight <- 1.0 - prob_finish
    if (surv_weight <= 0) return(memo_store(list()))
    scn <- .scenario_baseline(forced_complete, forced_survive)
    scn$weight <- surv_weight
    return(memo_store(list(scn)))
  }

  if (identical(kind, "not")) {
    scn <- .scenario_baseline(forced_complete, forced_survive)
    return(memo_store(list(scn)))
  }

  memo_store(list())
}

# Get CDF (cumulative probability) for an expression
.eval_expr_cdf <- function(expr, t, prep, component,
                           forced_complete = character(0),
                           forced_survive = character(0)) {
  complete_scenarios <- .expr_complete_scenarios(expr, t, prep, component,
                                                 forced_complete = forced_complete,
                                                 forced_survive = forced_survive)
  prob <- sum(vapply(complete_scenarios, `[[`, numeric(1), "weight"))
  prob <- min(max(prob, 0.0), 1.0)
  prob
}

# Get survival function for an expression
.eval_expr_survival <- function(expr, t, prep, component,
                                forced_complete = character(0),
                                forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive

  memo_env <- prep$.__memo
  cache_key <- NULL
  if (!is.null(memo_env) && is.environment(memo_env)) {
    cache_key <- .make_eval_key(expr, component, t,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                metric = "SURVVAL")
    cached <- memo_env[[cache_key]]
    if (!is.null(cached)) return(cached)
  }

  survival_scenarios <- .expr_survive_scenarios(expr, t, prep, component,
                                                forced_complete = forced_complete,
                                                forced_survive = forced_survive)
  prob <- sum(vapply(survival_scenarios, `[[`, numeric(1), "weight"))
  prob <- min(max(prob, 0.0), 1.0)
  if (!is.null(cache_key) && !is.null(memo_env) && is.environment(memo_env)) {
    memo_env[[cache_key]] <- prob
  }
  prob
}

# Compute likelihood (density or probability) for an expression
.eval_expr_likelihood <- function(expr, t, prep, component,
                                  forced_complete = character(0),
                                  forced_survive = character(0)) {
  normalized <- .normalize_forced_sets(forced_complete, forced_survive)
  forced_complete <- normalized$complete
  forced_survive <- normalized$survive

  memo_env <- prep$.__memo
  cache_key <- NULL
  if (!is.null(memo_env) && is.environment(memo_env)) {
    cache_key <- .make_eval_key(expr, component, t,
                                forced_complete = forced_complete,
                                forced_survive = forced_survive,
                                metric = "LIK")
    cached <- memo_env[[cache_key]]
    if (!is.null(cached)) return(cached)
  }

  finish_scenarios <- .expr_finish_scenarios(expr, t, prep, component,
                                             forced_complete = forced_complete,
                                             forced_survive = forced_survive)
  dens <- sum(vapply(finish_scenarios, `[[`, numeric(1), "weight"))
  dens <- max(0.0, dens)
  if (!is.null(cache_key) && !is.null(memo_env) && is.environment(memo_env)) {
    memo_env[[cache_key]] <- dens
  }
  dens
}

# ==============================================================================
# PART 5: Outcome Likelihood
# ==============================================================================

# Integrand function for probability integration
.integrand_outcome_density <- function(t, expr, prep, component,
                                       competitor_exprs = NULL,
                                       forced_complete = character(0),
                                       forced_survive = character(0),
                                       forced_complete_self = forced_complete,
                                       forced_survive_self = forced_survive) {
  vapply(t, function(tt) {
    if (!is.finite(tt) || tt < 0) return(0.0)
    scenarios <- .expr_finish_scenarios(expr, tt, prep, component,
                                        forced_complete = forced_complete_self,
                                        forced_survive = forced_survive_self)
    if (length(scenarios) == 0) return(0.0)
    total <- 0.0
    for (scn in scenarios) {
      dens <- scn$weight
      if (!is.finite(dens) || dens <= 0) next
      if (!is.null(competitor_exprs) && length(competitor_exprs) > 0) {
        comp_complete <- union(forced_complete, scn$comp_complete)
        comp_survive <- union(forced_survive, scn$comp_survive)
        surv_prod <- 1.0
        for (cexpr in competitor_exprs) {
          surv_val <- .eval_expr_survival(cexpr, tt, prep, component,
                                          forced_complete = comp_complete,
                                          forced_survive = comp_survive)
          surv_prod <- surv_prod * surv_val
          if (surv_prod == 0) break
        }
        dens <- dens * surv_prod
      }
      total <- total + dens
    }
    total
  }, numeric(1))
}

# Compute probability by integrating race density
# If competitor_exprs is provided, integrates:
#   f_r(t) * prod_j S_j(t)  over t in [0, upper]
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf, competitor_exprs = NULL,
                                           forced_complete = character(0),
                                           forced_survive = character(0),
                                           forced_complete_self = forced_complete,
                                           forced_survive_self = forced_survive) {
  if (!is.finite(upper_limit)) {
    integrand_vec <- function(u) {
      .integrand_outcome_density(
        u,
        expr = expr,
        prep = prep,
        component = component,
        competitor_exprs = competitor_exprs,
        forced_complete = forced_complete,
        forced_survive = forced_survive,
        forced_complete_self = forced_complete_self,
        forced_survive_self = forced_survive_self
      )
    }
    result <- tryCatch({
      stats::integrate(
        integrand_vec,
        lower = 0,
        upper = upper_limit,
        rel.tol = 1e-5,
        abs.tol = 1e-6,
        stop.on.error = FALSE
      )
    }, error = function(e) {
      list(value = 0.0, message = as.character(e))
    })
    val <- as.numeric(result$value)
    if (!is.finite(val)) val <- 0.0
    if (val < 0) val <- 0.0
    if (val > 1) val <- 1.0
    return(val)
  }

  integrand_scalar <- function(u) {
    val <- .integrand_outcome_density(
      u,
      expr = expr,
      prep = prep,
      component = component,
      competitor_exprs = competitor_exprs,
      forced_complete = forced_complete,
      forced_survive = forced_survive,
      forced_complete_self = forced_complete_self,
      forced_survive_self = forced_survive_self
    )
    if (length(val) == 0) 0.0 else val
  }

  key <- .make_integral_key(
    expr = expr,
    component = component,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    forced_complete_self = forced_complete_self,
    forced_survive_self = forced_survive_self,
    competitor_exprs = competitor_exprs
  )

  val <- .cached_integral(
    prep = prep,
    key = key,
    integrand_scalar = integrand_scalar,
    upper = upper_limit,
    rel.tol = 1e-5,
    abs.tol = 1e-6
  )
  if (val > 1) val <- 1.0
  val
}

# Detect shared-gate pattern (A & C) vs (B & C) across outcome definitions
# Find a single shared-gate pair among non-alias outcomes, if present
# Vectorized integrands for shared-gate pair; returns list(prob_x, prob_y)
.expr_to_branches <- function(expr) {
  if (is.null(expr) || is.null(expr$kind)) return(list())
  if (identical(expr$kind, "event")) {
    return(list(list(events = as.character(expr$source), guards = list())))
  }
  if (identical(expr$kind, "and")) {
    a <- .expr_to_branches(expr$args[[1]])
    b <- .expr_to_branches(expr$args[[2]])
    res <- list()
    if (length(a) == 0 || length(b) == 0) return(res)
    for (ba in a) {
      for (bb in b) {
        res[[length(res) + 1L]] <- list(
          events = unique(c(ba$events, bb$events)),
          guards = c(ba$guards %||% list(), bb$guards %||% list())
        )
      }
    }
    return(res)
  }
  if (identical(expr$kind, "or")) {
    res <- list()
    for (arg in expr$args) {
      brs <- .expr_to_branches(arg)
      if (length(brs) > 0) res <- c(res, brs)
    }
    return(res)
  }
  if (identical(expr$kind, "guard")) {
    inner <- .expr_to_branches(expr$reference)
    out <- list()
    for (br in inner) {
      br$guards <- c(br$guards, list(list(blocker = expr$blocker, unless = expr$unless %||% list())))
      out[[length(out) + 1L]] <- br
    }
    return(out)
  }
  # Unsupported in branch extraction
  list()
}

# Create a simple canonical key for a guard: only support event blocker and event protectors.
.guard_key <- function(guard) {
  blk <- guard$blocker
  if (is.null(blk) || is.null(blk$kind) || !identical(blk$kind, "event")) return(NULL)
  blk_id <- as.character(blk$source)
  unl_ids <- character(0)
  if (!is.null(guard$unless) && length(guard$unless) > 0) {
    for (un in guard$unless) {
      if (!is.null(un$kind) && identical(un$kind, "event")) unl_ids <- c(unl_ids, as.character(un$source))
      else return(NULL)
    }
    unl_ids <- sort(unique(unl_ids))
  }
  paste0(blk_id, "|unless:", paste(unl_ids, collapse = ","))
}

# Create a uniform guard descriptor for a branch: either NULL (no guards) or a single canonical guard key
# Require at most one unique guard per branch; otherwise return NA (unsupported)
.branch_uniform_guard_key <- function(branch) {
  if (length(branch$guards) == 0) return(NA_character_)
  keys <- vapply(branch$guards, .guard_key, character(1))
  if (any(is.na(keys) | is.null(keys) | keys == "")) return(NULL)
  keys <- unique(keys)
  if (length(keys) > 1) return(NULL)
  keys[[1]]
}

# Build shared families across outcomes. Each family:
# list(labels = c(label1, ...), G_ids = character(),
#      X_map = named list(label -> character vector of private ids),
#      guard = list(blocker_expr, unless_list) or NULL,
#      key = canonical string)
.build_shared_families <- function(prep, max_K = 4L, max_G = 3L) {
  outdefs <- prep$outcomes
  labels_all <- names(outdefs)
  # exclude alias-of and special sources
  core_labels <- Filter(function(lbl) {
    opt <- outdefs[[lbl]]$options %||% list()
    if (!is.null(opt$alias_of)) return(FALSE)
    ex <- outdefs[[lbl]]$expr
    if (!is.null(ex$kind) && identical(ex$kind, "event")) {
      src <- ex$source
      if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
    }
    TRUE
  }, labels_all)
  branches <- lapply(core_labels, function(lbl) .expr_to_branches(outdefs[[lbl]]$expr))
  names(branches) <- core_labels
  families <- list()
  family_index <- setNames(vector("list", length(core_labels)), core_labels)
  seen_keys <- character(0)

  # pairwise seed
  for (i_idx in seq_along(core_labels)) {
    li <- core_labels[[i_idx]]
    for (j_idx in seq_along(core_labels)) {
      if (j_idx <= i_idx) next
      lj <- core_labels[[j_idx]]
      for (bi in branches[[li]]) {
        gi_key <- .branch_uniform_guard_key(bi)
        if (is.null(gi_key)) next
        for (bj in branches[[lj]]) {
          gj_key <- .branch_uniform_guard_key(bj)
          if (is.null(gj_key)) next
          # require same guard key (or both NA)
          if (!identical(gi_key, gj_key)) next
          Gcand <- intersect(bi$events, bj$events)
          if (length(Gcand) == 0) next
          if (length(Gcand) > max_G) {
            stop(sprintf("Shared gate set size %d exceeds supported max %d (labels: %s,%s)", length(Gcand), max_G, li, lj))
          }
          # Build candidate label set that contain Gcand
          members <- list()
          members[[li]] <- bi
          members[[lj]] <- bj
          for (lk in core_labels) {
            if (lk %in% c(li, lj)) next
            # find first branch with same guard key and containing Gcand
            found <- NULL
            for (bk in branches[[lk]]) {
              gk_key <- .branch_uniform_guard_key(bk)
              if (!identical(gk_key, gi_key)) next
              if (all(Gcand %in% bk$events)) { found <- bk; break }
            }
            if (!is.null(found)) members[[lk]] <- found
          }
          lbls <- names(members)
          if (length(lbls) < 2) next
          if (length(lbls) > max_K) {
            stop(sprintf("Shared family size %d exceeds supported max %d (seed: %s,%s)", length(lbls), max_K, li, lj))
          }
          # Derive private sets and validate disjointness and non-empty
          Xmap <- list()
          ok <- TRUE
          priv_all <- character(0)
          for (nm in lbls) {
            Xi <- setdiff(members[[nm]]$events, Gcand)
            if (length(Xi) == 0) { ok <- FALSE; break }
            # enforce disjointness
            if (length(intersect(priv_all, Xi)) > 0) { ok <- FALSE; break }
            priv_all <- union(priv_all, Xi)
            Xmap[[nm]] <- Xi
          }
          if (!ok) next
          # Canonical key to avoid duplicates
          fam_key <- paste0(paste(sort(lbls), collapse = "+"), "|G:", paste(sort(Gcand), collapse = ","), "|guard:", as.character(gi_key))
          if (fam_key %in% seen_keys) next
          seen_keys <- c(seen_keys, fam_key)
          # Capture per-label guards (can differ per branch)
          guard_per_label <- lapply(lbls, function(nm) members[[nm]]$guards %||% list())
          names(guard_per_label) <- lbls
          fam <- list(labels = lbls, G_ids = Gcand, X_map = Xmap, guard_per_label = guard_per_label, key = fam_key)
          idx <- length(families) + 1L
          families[[idx]] <- fam
          for (nm in lbls) {
            family_index[[nm]] <- c(family_index[[nm]], idx)
          }
        }
      }
    }
  }
  # Ambiguity check: no label may belong to more than one family
  for (nm in names(family_index)) {
    ids <- family_index[[nm]]
    if (length(ids) > 1) {
      stop(sprintf("Outcome '%s' participates in multiple overlapping shared-gate families; ambiguous and not implemented", nm))
    }
  }
  list(families = families, index = family_index)
}

.make_event_expr <- function(id) {
  list(kind = "event", source = id)
}

.make_and_expr_from_ids <- function(ids) {
  ids <- as.character(ids)
  if (length(ids) == 1) return(.make_event_expr(ids[[1]]))
  ex <- .make_event_expr(ids[[1]])
  if (length(ids) == 1) return(ex)
  for (i in seq(2, length(ids))) {
    ex <- list(kind = "and", args = list(ex, .make_event_expr(ids[[i]])))
  }
  ex
}

.guard_eff_survival_fn <- function(guard_obj, prep, component) {
  if (is.null(guard_obj)) return(function(tt) rep(1.0, length(tt)))
  blk <- guard_obj$blocker
  unl <- guard_obj$unless %||% list()
  f_block <- function(u) { .eval_expr_likelihood(blk, u, prep, component) }
  S_prot <- function(u) {
    if (length(unl) == 0) return(rep(1.0, length(u)))
    vapply(u, function(ui) {
      prod(vapply(unl, function(un) .eval_expr_survival(un, ui, prep, component), numeric(1)))
    }, numeric(1))
  }
  function(tt) {
    vapply(tt, function(ti) {
      if (!is.finite(ti) || ti <= 0) return(1.0)
      val <- tryCatch(stats::integrate(function(u) {
        vapply(u, function(ui) f_block(ui) * S_prot(ui), numeric(1))
      }, lower = 0, upper = ti, rel.tol = 1e-5, abs.tol = 1e-6, stop.on.error = FALSE)$value, error = function(e) 0.0)
      s <- 1.0 - as.numeric(val)
      if (!is.finite(s)) 0.0 else max(0.0, min(1.0, s))
    }, numeric(1))
  }
}

# Combine multiple guard entries (list of {blocker, unless}) into a single survival function
.guard_eff_survival_fn_multi <- function(guard_list, prep, component) {
  if (is.null(guard_list) || length(guard_list) == 0) return(function(tt) rep(1.0, length(tt)))
  fns <- lapply(guard_list, function(g) .guard_eff_survival_fn(g, prep, component))
  function(tt) {
    vapply(tt, function(ti) {
      prod(vapply(fns, function(fn) fn(ti), numeric(1)))
    }, numeric(1))
  }
}

# Build integrand for a specific family member label
.family_integrand_for_label <- function(prep, component, family, label, competitor_exprs = list()) {
  labels <- family$labels
  G_ids <- as.character(family$G_ids)
  X_map <- family$X_map
  if (!(label %in% labels)) stop("label not in family")
  # Closures for gate members
  get_f <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti < 0) return(0.0)
        if (!is.null(prep$pools[[id]])) return(.pool_density(prep, id, component, ti))
        if (!is.null(prep$accumulators[[id]])) return(.acc_density(prep$accumulators[[id]], ti))
        0.0
      }, numeric(1))
    }
  }
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        if (!is.null(prep$pools[[id]])) return(.pool_survival(prep, id, component, ti))
        if (!is.null(prep$accumulators[[id]])) return(.acc_survival(prep$accumulators[[id]], ti))
        0.0
      }, numeric(1))
    }
  }
  get_F <- function(id) { Sfun <- get_S(id); function(tt) { 1.0 - Sfun(tt) } }
  # Private set closures per label using AND expressions
  X_exprs <- lapply(labels, function(l) .make_and_expr_from_ids(X_map[[l]]))
  names(X_exprs) <- labels
  fX <- function(l) { function(tt) { .eval_expr_likelihood(X_exprs[[l]], tt, prep, component) } }
  FX <- function(l) { function(tt) { .eval_expr_cdf(X_exprs[[l]], tt, prep, component) } }
  SX <- function(l) { function(tt) { .eval_expr_survival(X_exprs[[l]], tt, prep, component) } }
  # Guarded CDFs and survivals for private sets
  FXg <- list()
  SXg <- list()
  for (nm in labels) {
    FXg[[nm]] <- (function(nm_inner) {
      function(tt) {
        vapply(tt, function(ti) {
          if (!is.finite(ti) || ti <= 0) return(0.0)
          integ <- tryCatch(stats::integrate(function(u) {
            vapply(u, function(ui) {
              dens <- fX(nm_inner)(ui)
              seff <- S_eff_label[[nm_inner]](ui)
              if (!is.finite(dens) || dens <= 0 || !is.finite(seff) || seff <= 0) 0.0 else dens * seff
            }, numeric(1))
          }, lower = 0, upper = ti, subdivisions = 1200L, rel.tol = 1e-6, abs.tol = 1e-8, stop.on.error = FALSE)$value, error = function(e) 0.0)
          val <- as.numeric(integ)
          if (!is.finite(val) || val < 0) 0.0 else if (val > 1) 1.0 else val
        }, numeric(1))
      }
    })(nm)
    SXg[[nm]] <- (function(nm_inner) {
      function(tt) { 1.0 - FXg[[nm_inner]](tt) }
    })(nm)
  }
  # Gate derived closures
  Fg_list <- lapply(G_ids, get_F)
  fg_list <- lapply(G_ids, get_f)
  # Build per-label guard survival (labels may differ in inhibition structure)
  guard_per_label <- family$guard_per_label %||% setNames(vector("list", length(labels)), labels)
  S_eff_label <- lapply(labels, function(nm) .guard_eff_survival_fn_multi(guard_per_label[[nm]] %||% list(), prep, component))
  names(S_eff_label) <- labels
  # helper: product survival of competitors (not in family)
  comp_surv <- function(tt) {
    if (length(competitor_exprs) == 0) return(rep(1.0, length(tt)))
    vapply(tt, function(ti) {
      prod(vapply(competitor_exprs, function(ce) .eval_expr_survival(ce, ti, prep, component), numeric(1)))
    }, numeric(1))
  }
  # Numerator_i(t) via inclusion-exclusion (exact for any K):
  # numerator_i(t) = Σ_{S⊆others} (-1)^{|S|} [ Π_{k∈others\S} F_k(t) ] × ∫_0^t f_i(u) Π_{j∈S} F_j(u) du
  numerator_i <- function(i_label) {
    fXi <- fX(i_label)
    others <- setdiff(labels, i_label)
    FX_others <- lapply(others, function(l) FX(l))
    idxs <- seq_along(others)
    subsets <- list(integer(0))
    if (length(idxs) > 0) {
      for (m in 1:length(idxs)) {
        cmb <- utils::combn(idxs, m, simplify = FALSE)
        subsets <- c(subsets, cmb)
      }
    }
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti <= 0) return(0.0)
        total <- 0.0
        for (S in subsets) {
          s_size <- length(S)
          sign <- if ((s_size %% 2) == 0) 1.0 else -1.0
          # Product over k not in S of F_k(t)
          if (length(idxs) > 0) {
            comp_idx <- setdiff(idxs, S)
            prod_F_t <- 1.0
            for (k in comp_idx) { prod_F_t <- prod_F_t * FX_others[[k]](ti) }
          } else {
            prod_F_t <- 1.0
          }
          if (s_size == 0) {
            # ∫ f_i(u) du = F_i(t)
            Fi_t <- FX(i_label)(ti)
            total <- total + sign * prod_F_t * Fi_t
            next
          }
          integrand <- function(u) {
            u <- as.numeric(u)
            dens_vals <- fXi(u)
            if (length(dens_vals) == 1 && length(u) > 1) dens_vals <- rep(dens_vals, length(u))
            out <- numeric(length(u))
            for (idx in seq_along(u)) {
              dval <- dens_vals[[idx]]
              if (!is.finite(dval) || dval <= 0) { out[[idx]] <- 0.0; next }
              prod_Fu <- 1.0
              for (k in S) { prod_Fu <- prod_Fu * FX_others[[k]](u[[idx]]) }
              out[[idx]] <- dval * prod_Fu
            }
            out
          }
          val <- tryCatch(stats::integrate(integrand, lower = 0, upper = ti,
                                           subdivisions = 2000L, rel.tol = 1e-9, abs.tol = 1e-10, stop.on.error = FALSE)$value,
                          error = function(e) 0.0)
          total <- total + sign * prod_F_t * as.numeric(val)
        }
        if (!is.finite(total) || total < 0) total <- 0.0
        total
      }, numeric(1))
    }
  }
  # Gate-last factor Σ_g f_g(t) Π_{g'≠g} F_{g'}(t)
  gate_last_sum <- function(tt) {
    vapply(tt, function(ti) {
      s <- 0.0
      for (g_idx in seq_along(G_ids)) {
        fg <- fg_list[[g_idx]]; prodF <- 1.0
        for (h_idx in seq_along(G_ids)) {
          if (h_idx == g_idx) next
          prodF <- prodF * Fg_list[[h_idx]](ti)
        }
        s <- s + fg(ti) * prodF
      }
      s
    }, numeric(1))
  }
  # F_maxG(t)
  FmaxG <- function(tt) {
    vapply(tt, function(ti) { prod(vapply(Fg_list, function(Ffun) Ffun(ti), numeric(1))) }, numeric(1))
  }
  if (length(labels) == 2L) {
    # Pair-case using the general guarded formulation explicitly
    other <- setdiff(labels, label)[[1]]
    fXi <- fX(label)
    FXi_g <- FXg[[label]]
    SXo_g <- SXg[[other]]
    FYo_g <- FXg[[other]]
    # Guard survival for this label at time t
    S_i_fn <- S_eff_label[[label]]
    # helper for ∫ fXi(u) * S_i(u) * FYo_g(u) du
    int_guarded_fXi_FYo <- function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti <= 0) return(0.0)
        tryCatch(stats::integrate(function(u) {
          vapply(u, function(ui) {
            val_f <- fXi(ui)
            if (!is.finite(val_f) || val_f <= 0) return(0.0)
            s_i_u <- S_i_fn(ui)
            fyo_u <- FYo_g(ui)
            if (!is.finite(s_i_u) || s_i_u <= 0 || !is.finite(fyo_u) || fyo_u < 0) 0.0 else (val_f * s_i_u * fyo_u)
          }, numeric(1))
        }, lower = 0, upper = ti, rel.tol = 1e-6, abs.tol = 1e-8, stop.on.error = FALSE)$value, error = function(e) 0.0)
      }, numeric(1))
    }
    function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        # Per-label guard survival at t
        S_i <- S_i_fn(tt)
        # Xi last (guard-conditioned via S_i at t, others via SXo_g)
        term_X_last <- fXi(tt) * FmaxG(tt)
        term_X_last <- term_X_last * SXo_g(tt)
        term_X_last <- term_X_last * S_i
        # Gate last, non-tie: use guarded CDF FXi_g
        term_gate_other <- gate_last_sum(tt) * FXi_g(tt)
        term_gate_other <- term_gate_other * SXo_g(tt)
        term_gate_other <- term_gate_other * S_i
        # Gate last tie with proportional allocation using guarded quantities
        denom <- FXi_g(tt) * FYo_g(tt)
        palloc <- 0.0
        if (is.finite(denom) && denom > 0) {
          num_int <- int_guarded_fXi_FYo(tt)
          palloc <- 1.0 - (as.numeric(num_int) / as.numeric(denom))
          if (!is.finite(palloc)) palloc <- 0.0
          if (palloc < 0) palloc <- 0.0
          if (palloc > 1) palloc <- 1.0
        }
        term_gate_tie <- gate_last_sum(tt) * denom * palloc
        # Multiply tie term by guard survival at t for this label
        term_gate_tie <- term_gate_tie * S_i
        base <- term_X_last + term_gate_other + term_gate_tie
        base <- base * comp_surv(tt)
        base
      }, numeric(1))
    }
  } else {
    # General K, |G| form
    lab <- label
    fXi <- fX(lab); FXi <- FXg[[lab]]; SXi <- SX(lab)
    SX_others <- lapply(setdiff(labels, lab), function(l) SX(l))
    num_fun_lab <- numerator_i(lab)
    # Also compute numerators for other labels to allow normalization
    num_funs_all <- setNames(lapply(labels, numerator_i), labels)
    FX_all_prod <- function(tt) {
      vapply(tt, function(ti) {
        prod(vapply(labels, function(l) FX(l)(ti), numeric(1)))
      }, numeric(1))
    }
    # Dynamic inclusion-exclusion guarded numerator at time t over active others
    compute_numerator_for_i <- function(ti, i_label, active_others) {
      if (length(active_others) == 0) {
        return(FXg[[i_label]](ti))
      }
      idxs <- seq_along(active_others)
      total <- 0.0
      # iterate over all subsets S of active_others
      subsets <- list(integer(0))
      if (length(idxs) > 0) {
        for (m in 1:length(idxs)) {
          cmb <- utils::combn(idxs, m, simplify = FALSE)
          subsets <- c(subsets, cmb)
        }
      }
      for (S in subsets) {
        s_size <- length(S)
        sign <- if ((s_size %% 2) == 0) 1.0 else -1.0
        # product over k not in S of F_k(t)
        if (length(idxs) > 0) {
          comp_idx <- setdiff(idxs, S)
          prod_F_t <- 1.0
          for (k in comp_idx) { prod_F_t <- prod_F_t * FXg[[active_others[[k]]]](ti) }
        } else {
          prod_F_t <- 1.0
        }
        if (s_size == 0) {
          Fi_t <- FXg[[i_label]](ti)
          total <- total + sign * prod_F_t * Fi_t
          next
        }
        integrand <- function(u) {
          u <- as.numeric(u)
          dens_vals <- fX(i_label)(u)
          if (length(dens_vals) == 1 && length(u) > 1) dens_vals <- rep(dens_vals, length(u))
          out <- numeric(length(u))
          for (idx in seq_along(u)) {
            dval <- dens_vals[[idx]]
            if (!is.finite(dval) || dval <= 0) { out[[idx]] <- 0.0; next }
            seff_i <- S_eff_label[[i_label]](u[[idx]])
            prod_Fu <- 1.0
            for (k in S) { prod_Fu <- prod_Fu * FXg[[active_others[[k]]]](u[[idx]]) }
            out[[idx]] <- (dval * seff_i) * prod_Fu
          }
          out
        }
        val <- tryCatch(stats::integrate(integrand, lower = 0, upper = ti,
                                         subdivisions = 2000L, rel.tol = 1e-9, abs.tol = 1e-10, stop.on.error = FALSE)$value,
                        error = function(e) 0.0)
        total <- total + sign * as.numeric(prod_F_t * val)
      }
      if (!is.finite(total) || total < 0) total <- 0.0
      total
    }

    function(t) {
      vapply(t, function(tt) {
        if (!is.finite(tt) || tt < 0) return(0.0)
        # Per-label guard survival at t
        S_i <- S_eff_label[[lab]](tt)
        # Active labels at t
        active_labels <- labels[vapply(labels, function(nm) S_eff_label[[nm]](tt) > 0, logical(1))]
        active_others <- setdiff(active_labels, lab)
        # Xi-last
        term_X_last <- fXi(tt) * FmaxG(tt)
        if (length(active_others) > 0) {
          for (nm in active_others) { term_X_last <- term_X_last * SXg[[nm]](tt) }
        }
        term_X_last <- term_X_last * S_i
        # Gate-last, non-tie
        term_gate_other <- gate_last_sum(tt) * FXi(tt)
        if (length(active_others) > 0) {
          for (nm in active_others) { term_gate_other <- term_gate_other * SXg[[nm]](tt) }
        }
        term_gate_other <- term_gate_other * S_i
        # Gate-last tie share: total tie mass × normalized guard-weighted numerators
        # Guard-weighted numerator for this label
        num_i_eff <- if (lab %in% active_labels) S_i * compute_numerator_for_i(tt, lab, active_others) else 0.0
        # Sum over labels of guard-weighted numerators
        sum_num_eff <- 0.0
        for (nm in active_labels) {
          S_nm <- S_eff_label[[nm]](tt)
          sum_num_eff <- sum_num_eff + (S_nm * compute_numerator_for_i(tt, nm, setdiff(active_labels, nm)))
        }
        eligible <- vapply(labels, function(nm) S_eff_label[[nm]](tt) > 0, logical(1))
        elig_labels <- labels[eligible]
        # Eligible-only FX product
        FX_elig_prod <- if (length(active_labels) > 0) {
          prod(vapply(active_labels, function(nm) FXg[[nm]](tt), numeric(1)))
        } else 0.0
        # Guard-weighted numerator sum over eligible labels
        if (length(active_labels) == 0) {
          weight_lab <- 0.0
        } else if (!is.finite(sum_num_eff) || sum_num_eff <= 0) {
          # Fallback: equal split among eligible
          n_elig <- length(active_labels)
          weight_lab <- if (lab %in% active_labels) 1.0 / n_elig else 0.0
        } else {
          weight_lab <- if (lab %in% active_labels) as.numeric(num_i_eff / sum_num_eff) else 0.0
        }
        tie_mass <- gate_last_sum(tt) * FX_elig_prod
        term_gate_tie <- tie_mass * weight_lab
        # Apply per-label guard survival at t to tie allocation as required
        term_gate_tie <- term_gate_tie * S_i
        base <- term_X_last + term_gate_other + term_gate_tie
        base <- base * comp_surv(tt)
        base
      }, numeric(1))
    }
  }
}

.ensure_shared_families <- function(prep) {
  if (!is.null(prep$shared_families) && !is.null(prep$shared_family_index)) return(prep)
  res <- .build_shared_families(prep)
  prep$shared_families <- res$families
  prep$shared_family_index <- res$index
  prep
}

# Detect OR of two guarded ANDs sharing a common blocker and a common k=1 gate member
# Pattern: or( guard(and(X, K1), by=G, unless=?), guard(and(Y, K2), by=G, unless=?) )
# where K1 and K2 are pools whose member sets overlap in exactly one accumulator Z.
# Returns list(x_id, y_id, z_id, blocker_expr, unless_list) or NULL if not matched.
.detect_or_guard_shared_k1 <- function(expr, prep) {
  if (is.null(expr) || is.null(expr$kind) || !identical(expr$kind, "or")) return(NULL)
  if (length(expr$args) != 2) return(NULL)
  a <- expr$args[[1]]; b <- expr$args[[2]]
  if (is.null(a$kind) || !identical(a$kind, "guard")) return(NULL)
  if (is.null(b$kind) || !identical(b$kind, "guard")) return(NULL)
  # Both guards must share the same blocker (structure equality)
  # We'll be conservative and require source id equality when blocker is an event.
  blocker_a <- a$blocker; blocker_b <- b$blocker
  # quick structural check: if both blocker are events, their source must match
  same_blocker <- FALSE
  if (!is.null(blocker_a$kind) && !is.null(blocker_b$kind)) {
    if (identical(blocker_a$kind, "event") && identical(blocker_b$kind, "event")) {
      same_blocker <- identical(blocker_a$source, blocker_b$source)
    } else {
      # For non-event blockers, fall back to a loose check on deparsed form
      same_blocker <- TRUE
    }
  }
  if (!same_blocker) return(NULL)
  ref_a <- a$reference; ref_b <- b$reference
  if (is.null(ref_a$kind) || !identical(ref_a$kind, "and")) return(NULL)
  if (is.null(ref_b$kind) || !identical(ref_b$kind, "and")) return(NULL)
  if (length(ref_a$args) != 2 || length(ref_b$args) != 2) return(NULL)
  # Helper: identify (primary,event) vs (gate,pool-or-event)
  pick_xy <- function(and_node) {
    e1 <- and_node$args[[1]]; e2 <- and_node$args[[2]]
    get_evt_id <- function(node) {
      if (!is.null(node$kind) && identical(node$kind, "event")) return(node$source)
      NULL
    }
    id1 <- get_evt_id(e1); id2 <- get_evt_id(e2)
    # Prefer choosing accumulator id as primary (X/Y)
    if (!is.null(id1) && !is.null(prep$accumulators[[id1]])) {
      primary <- id1; other_id <- id2
    } else if (!is.null(id2) && !is.null(prep$accumulators[[id2]])) {
      primary <- id2; other_id <- id1
    } else {
      # Fallback: take first as primary if neither is a direct accumulator
      primary <- id1 %||% id2
      other_id <- if (!identical(primary, id1)) id1 else id2
    }
    list(primary = primary, other = other_id)
  }
  xa <- pick_xy(ref_a); xb <- pick_xy(ref_b)
  if (is.null(xa$primary) || is.null(xb$primary)) return(NULL)
  # Flatten pool members of the "gate" side to underlying accumulator ids
  flatten_pool_accs <- function(id) {
    accs <- character(0)
    if (is.null(id)) return(accs)
    if (!is.null(prep$accumulators[[id]])) return(id)
    if (!is.null(prep$pools[[id]])) {
      members <- prep$pools[[id]]$members
      for (m in members) {
        accs <- union(accs, flatten_pool_accs(m))
      }
    }
    accs
  }
  set_a <- flatten_pool_accs(xa$other)
  set_b <- flatten_pool_accs(xb$other)
  if (length(set_a) == 0 || length(set_b) == 0) return(NULL)
  shared <- intersect(set_a, set_b)
  if (length(shared) != 1) return(NULL)
  list(
    x_id = xa$primary,
    y_id = xb$primary,
    z_id = shared[[1]],
    blocker_expr = blocker_a,
    unless_list = a$unless %||% list()
  )
}

# Build integrand for OR of guarded ANDs sharing a single common Z in their k=1 gates
# Density at time t for the combined outcome (union of the two guarded paths):
#   [ fX(t) * FZ(t) * SY(t) + fY(t) * FZ(t) * SX(t) + fZ(t) * FX(t) * FY(t) ]
#   * S_eff_blocker(t)
# Optionally multiplied by survival of competitor outcomes at t (if provided).
.or_guard_shared_k1_integrand <- function(prep, component, info, competitor_exprs = list(),
                                   rel.tol = 1e-5, abs.tol = 1e-6) {
  get_f <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti) || ti < 0) return(0.0)
        if (!is.null(prep$pools[[id]])) return(.pool_density(prep, id, component, ti))
        if (!is.null(prep$accumulators[[id]])) return(.acc_density(prep$accumulators[[id]], ti))
        0.0
      }, numeric(1))
    }
  }
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        if (!is.null(prep$pools[[id]])) return(.pool_survival(prep, id, component, ti))
        if (!is.null(prep$accumulators[[id]])) return(.acc_survival(prep$accumulators[[id]], ti))
        0.0
      }, numeric(1))
    }
  }
  get_F <- function(id) {
    Sfun <- get_S(id)
    function(tt) { 1.0 - Sfun(tt) }
  }
  # Blocker effective survival with optional protectors (unless)
  S_eff_blocker <- function(tt) {
    vapply(tt, function(ti) {
      if (!is.finite(ti) || ti <= 0) return(1.0)
      f_block <- function(u) .eval_expr_likelihood(info$blocker_expr, u, prep, component)
      S_prot <- function(u) {
        if (length(info$unless_list) == 0) return(1.0)
        prod(vapply(info$unless_list, function(un) .eval_expr_survival(un, u, prep, component), numeric(1)))
      }
      val <- tryCatch(stats::integrate(function(u) {
        vapply(u, function(ui) f_block(ui) * S_prot(ui), numeric(1))
      }, lower = 0, upper = ti, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = FALSE)$value,
      error = function(e) 0.0)
      s <- 1.0 - as.numeric(val)
      if (!is.finite(s)) 0.0 else max(0.0, min(1.0, s))
    }, numeric(1))
  }
  # Closures
  fX <- get_f(info$x_id); SX <- get_S(info$x_id); FX <- get_F(info$x_id)
  fY <- get_f(info$y_id); SY <- get_S(info$y_id); FY <- get_F(info$y_id)
  fZ <- get_f(info$z_id); FZ <- get_F(info$z_id)

  # Guard-conditioned CDFs for X and Y (integrate density times blocker survival up to t)
  FXg <- function(tt) {
    vapply(tt, function(ti) {
      if (!is.finite(ti) || ti <= 0) return(0.0)
      val <- tryCatch(stats::integrate(function(u) {
        vapply(u, function(ui) {
          dens <- fX(ui)
          seff <- S_eff_blocker(ui)
          if (!is.finite(dens) || dens <= 0 || !is.finite(seff) || seff <= 0) 0.0 else dens * seff
        }, numeric(1))
      }, lower = 0, upper = ti, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = FALSE)$value, error = function(e) 0.0)
      as.numeric(val)
    }, numeric(1))
  }
  FYg <- function(tt) {
    vapply(tt, function(ti) {
      if (!is.finite(ti) || ti <= 0) return(0.0)
      val <- tryCatch(stats::integrate(function(u) {
        vapply(u, function(ui) {
          dens <- fY(ui)
          seff <- S_eff_blocker(ui)
          if (!is.finite(dens) || dens <= 0 || !is.finite(seff) || seff <= 0) 0.0 else dens * seff
        }, numeric(1))
      }, lower = 0, upper = ti, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = FALSE)$value, error = function(e) 0.0)
      as.numeric(val)
    }, numeric(1))
  }

  integrand <- function(t) {
    vapply(t, function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      # Xi last (guard on blocker at tt)
      term_x_last <- fX(tt) * FZ(tt) * SY(tt)
      term_y_last <- fY(tt) * FZ(tt) * SX(tt)
      # Guard-conditioned tie allocation
      denom <- FXg(tt) * FYg(tt)
      palloc_x <- 0.0
      if (is.finite(denom) && denom > 0) {
        num_int <- tryCatch(stats::integrate(function(u) {
          vapply(u, function(ui) fX(ui) * FYg(ui), numeric(1))
        }, lower = 0, upper = tt, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = FALSE)$value, error = function(e) 0.0)
        palloc_x <- 1.0 - (as.numeric(num_int) / as.numeric(denom))
        if (!is.finite(palloc_x)) palloc_x <- 0.0
        if (palloc_x < 0) palloc_x <- 0.0
        if (palloc_x > 1) palloc_x <- 1.0
      }
      term_tie <- fZ(tt) * FXg(tt) * FYg(tt) * palloc_x
      base <- term_x_last + term_y_last + term_tie
      if (length(competitor_exprs) > 0) {
        comp_surv <- 1.0
        for (ce in competitor_exprs) {
          comp_surv <- comp_surv * .eval_expr_survival(ce, tt, prep, component)
          if (comp_surv == 0) break
        }
        base <- base * comp_surv
      }
      base * S_eff_blocker(tt)
    }, numeric(1))
  }
  integrand
}

# Compute likelihood for a single trial/outcome
.outcome_likelihood <- function(outcome_label, rt, prep, component) {
  outcome_defs <- prep$outcomes
  
  # FIRST: Check for GUESS outcomes (can be both defined and from guess policy)
  # GUESS outcomes arise from component-level guess policies
  if (identical(outcome_label, "GUESS")) {
    # Get the guess policy for this component
    guess_policy <- .get_component_attr(prep, component, "guess")
    
    if (!is.null(guess_policy) && !is.null(guess_policy$weights)) {
      # The guess policy has:
      # - outcome: the label for guess outcomes
      # - weights: named vector of keep probabilities
      
      guess_weights <- guess_policy$weights
      
      # Probability of GUESS outcome is:
      # Sum over outcomes: P(outcome fires) * P(guess | outcome)
      # P(guess | outcome) = 1 - weight[outcome]
      
      deadline <- .get_component_attr(prep, component, "deadline")
      deadline <- deadline %||% prep$default_deadline
      
      if (is.finite(deadline)) {
        # GUESS probability/density comes from diverting mass from base labels
        if (is.na(rt)) {
          total_prob <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            out_def <- outcome_defs[[label]]
            expr <- out_def$expr
            comp_labels <- setdiff(names(outcome_defs), label)
            if (length(comp_labels) > 0) {
              comp_labels <- Filter(function(lbl) {
                expr_lbl <- outcome_defs[[lbl]]$expr
                if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
                  src <- expr_lbl$source
                  if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                }
                TRUE
              }, comp_labels)
            }
            competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]]$expr) else list()
            prob_outcome <- .integrate_outcome_probability(expr, prep, component, deadline, competitor_exprs = competitor_exprs)
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_prob <- total_prob + prob_outcome * guess_prob
          }
          return(total_prob)
        } else {
          total_dens <- 0.0
          for (label in names(guess_weights)) {
            if (!label %in% names(outcome_defs)) next
            out_def <- outcome_defs[[label]]
            expr <- out_def$expr
            # density of label at rt times survival of all competitors excluding this label
            dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
            if (dens_r == 0) next
            comp_labels <- setdiff(names(outcome_defs), label)
            if (length(comp_labels) > 0) {
              comp_labels <- Filter(function(lbl) {
                expr_lbl <- outcome_defs[[lbl]]$expr
                if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
                  src <- expr_lbl$source
                  if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
                }
                TRUE
              }, comp_labels)
            }
            competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]]$expr) else list()
            surv_prod <- if (length(competitor_exprs) > 0) prod(vapply(competitor_exprs, function(ce) .eval_expr_survival(ce, rt, prep, component), numeric(1))) else 1.0
            keep_prob <- as.numeric(guess_weights[[label]])
            guess_prob <- 1.0 - keep_prob
            total_dens <- total_dens + guess_prob * dens_r * surv_prod
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
        if (!is.null(out_def$options$map_outcome_to)) {
          mapped <- out_def$options$map_outcome_to
          if (is.na(mapped) || identical(mapped, "NA")) {
            # Found the outcome that maps to NA
            expr <- out_def$expr
            if (is.na(rt)) {
              # Integrate race probability for the mapped source outcome
              deadline <- .get_component_attr(prep, component, "deadline")
              deadline <- deadline %||% prep$default_deadline
              # competitors exclude the mapped source label
              comp_labels_map <- setdiff(names(outcome_defs), label)
              competitor_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]]$expr) else list()
              return(.integrate_outcome_probability(expr, prep, component, deadline, competitor_exprs = competitor_exprs_map))
            } else {
              # Race density at rt for the mapped source outcome
              dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
              if (dens_r == 0) return(0.0)
              comp_labels_map <- setdiff(names(outcome_defs), label)
              competitor_exprs_map <- if (length(comp_labels_map) > 0) lapply(comp_labels_map, function(lbl) outcome_defs[[lbl]]$expr) else list()
              surv_comp <- if (length(competitor_exprs_map) > 0) prod(vapply(competitor_exprs_map, function(ce) .eval_expr_survival(ce, rt, prep, component), numeric(1))) else 1.0
              return(dens_r * surv_comp)
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
      deadline <- deadline %||% prep$default_deadline
      
      if (!is.finite(deadline)) return(0.0)
      
      # P(no outcome by deadline) = Π_i S_i(deadline)
      # where i ranges over all possible outcomes
      prob_none <- 1.0
      for (label in names(outcome_defs)) {
        out_def <- outcome_defs[[label]]
        expr <- out_def$expr
        # Survival of this outcome at deadline
        surv <- .eval_expr_survival(expr, deadline, prep, component)
        prob_none <- prob_none * surv
      }
      return(prob_none)
    }
    
    warning(sprintf("Unknown outcome '%s'", outcome_label))
    return(0.0)
  }
  
  expr <- outcome_def$expr
  options <- outcome_def$options %||% list()
  
  # Check if outcome is mapped to something else
  if (!is.null(options$map_outcome_to)) {
    # This outcome is mapped, but we shouldn't reach here if data has the mapped label
    # Just proceed with the expression
  }

  # If this outcome has an outcome-level guess policy, it never appears as-is.
  # Its mass is redistributed to target labels via donors handled below.
  if (!is.null(options$guess)) {
    return(0.0)
  }
  
  # Build competitor exprs for this outcome (exclude alias-of, special outcomes, and siblings mapping to the same top-level label)
  comp_labels <- setdiff(names(outcome_defs), outcome_label)
  if (length(comp_labels) > 0) {
    # Determine this outcome's top-level target (itself or map_outcome_to)
    this_target <- options$map_outcome_to %||% outcome_label
    comp_labels <- Filter(function(lbl) {
      # skip alias-of
      if (!is.null(outcome_defs[[lbl]]$options$alias_of)) return(FALSE)
      # skip special sources
      expr_lbl <- outcome_defs[[lbl]]$expr
      if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
        src <- expr_lbl$source
        if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
      }
      # skip siblings that map to the same top-level target label
      comp_target <- outcome_defs[[lbl]]$options$map_outcome_to %||% lbl
      if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
      TRUE
    }, comp_labels)
  }
  competitor_exprs <- if (length(comp_labels) > 0) lapply(comp_labels, function(lbl) outcome_defs[[lbl]]$expr) else list()

  # Build donor (guess) contributions: outcomes that guess into this label
  donors <- list()
  if (length(comp_labels) > 0) {
    for (dlbl in comp_labels) {
      dopt <- outcome_defs[[dlbl]]$options$guess %||% NULL
      if (is.null(dopt)) next
      glabels <- dopt$labels %||% NULL
      gweights <- dopt$weights %||% NULL
      if (is.null(glabels) || is.null(gweights)) next
      idx <- which(as.character(glabels) == outcome_label)
      if (length(idx) == 0) next
      w <- as.numeric(gweights[[idx[[1]]]])
      rt_policy <- dopt$rt_policy %||% "keep"
      donors[[length(donors) + 1L]] <- list(expr = outcome_defs[[dlbl]]$expr, weight = w, rt_policy = rt_policy)
    }
  }

  # No special per-example fast paths; use general family detection

  # Helper closures to get density/CDF for pools or accumulators
  get_f <- function(id) {
    function(tt) {
      if (!is.finite(tt) || tt < 0) return(0.0)
      if (!is.null(prep$pools[[id]])) return(.pool_density(prep, id, component, tt))
      if (!is.null(prep$accumulators[[id]])) return(.acc_density(prep$accumulators[[id]], tt))
      0.0
    }
  }
  get_S <- function(id) {
    function(tt) {
      if (!is.finite(tt)) return(0.0)
      if (tt < 0) return(1.0)
      if (!is.null(prep$pools[[id]])) return(.pool_survival(prep, id, component, tt))
      if (!is.null(prep$accumulators[[id]])) return(.acc_survival(prep$accumulators[[id]], tt))
      0.0
    }
  }
  get_F <- function(id) {
    Sfun <- get_S(id)
    function(tt) { 1.0 - Sfun(tt) }
  }

  # Handle NA/infinite RT by integration of race density
  if (is.na(rt) || !is.finite(rt) || rt < 0) {
    deadline <- .get_component_attr(prep, component, "deadline")
    deadline <- deadline %||% prep$default_deadline
    # Prefer general shared-family fast path (exact) when available
    prep <- .ensure_shared_families(prep)
    fam_ids <- prep$shared_family_index[[outcome_label]] %||% integer(0)
    if (length(fam_ids) == 1) {
      fam <- prep$shared_families[[fam_ids[[1]]]]
      # Build competitor exprs excluding family members
      comp_labels_ex <- setdiff(names(outcome_defs), union(outcome_label, fam$labels))
      if (length(comp_labels_ex) > 0) {
        # Exclude alias-of, special, and siblings mapping to same top-level target
        this_target <- options$map_outcome_to %||% outcome_label
        comp_labels_ex <- Filter(function(lbl) {
          if (!is.null(outcome_defs[[lbl]]$options$alias_of)) return(FALSE)
          expr_lbl <- outcome_defs[[lbl]]$expr
          if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
            src <- expr_lbl$source
            if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
          }
          comp_target <- outcome_defs[[lbl]]$options$map_outcome_to %||% lbl
          if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
          TRUE
        }, comp_labels_ex)
      }
      comp_exprs_ex <- if (length(comp_labels_ex) > 0) lapply(comp_labels_ex, function(lbl) outcome_defs[[lbl]]$expr) else list()
      # Pair-case with any shared gate set: P(label) = P(X_label < X_other)
      if (length(fam$labels) == 2L) {
        labels_pair <- fam$labels
        other_label <- setdiff(labels_pair, outcome_label)[[1]]
        # Private expressions (AND of private members)
        X_i_expr <- .make_and_expr_from_ids(fam$X_map[[outcome_label]])
        X_o_expr <- .make_and_expr_from_ids(fam$X_map[[other_label]])
        # Closures for density and survival of private sets
        fXi <- function(u) { .eval_expr_likelihood(X_i_expr, u, prep, component) }
        SXo <- function(u) { .eval_expr_survival(X_o_expr, u, prep, component) }
        # Integrand for P(X_i < X_other) with outside-family competitors
        integrand <- function(t) {
          vapply(t, function(tt) {
            if (!is.finite(tt) || tt < 0) return(0.0)
            base <- fXi(tt) * SXo(tt)
            if (base == 0) return(0.0)
            comp_surv <- if (length(comp_exprs_ex) == 0) 1.0 else prod(vapply(comp_exprs_ex, function(ce) .eval_expr_survival(ce, tt, prep, component), numeric(1)))
            base * comp_surv
          }, numeric(1))
        }
        keep_prob <- 1.0
        gp_comp <- .get_component_attr(prep, component, "guess")
        if (!is.null(gp_comp) && !is.null(gp_comp$weights)) {
          kp <- gp_comp$weights[[outcome_label]] %||% gp_comp$weights[[normalize_label(outcome_label)]] %||% NULL
          if (!is.null(kp)) keep_prob <- as.numeric(kp)
        }
        upper_lim <- if (is.finite(deadline)) deadline else Inf
        res <- tryCatch(stats::integrate(integrand, lower = 0, upper = upper_lim,
                                         rel.tol = 1e-5, abs.tol = 1e-6, stop.on.error = FALSE)$value,
                        error = function(e) 0.0)
        return(as.numeric(keep_prob * res))
      }
      # General K>=2 fallback: exact guarded family integrand
      integ <- .family_integrand_for_label(prep, component, fam, outcome_label, competitor_exprs = comp_exprs_ex)
      # Cache probability integral too for repeated calls
      cache_env <- prep$.__cache_family_prob__ %||% new.env(parent = emptyenv())
      prep$.__cache_family_prob__ <- cache_env
      cache_key <- paste(component, outcome_label, fam$key, "prob", sep = "|")
      cached <- cache_env[[cache_key]]
      if (is.null(cached)) {
        upper <- if (is.finite(deadline)) deadline else Inf
        res <- tryCatch(stats::integrate(integ, lower = 0, upper = upper,
                                         rel.tol = 1e-5, abs.tol = 1e-6, stop.on.error = FALSE)$value,
                        error = function(e) 0.0)
        cache_env[[cache_key]] <- as.numeric(res)
        cached <- cache_env[[cache_key]]
      }
      keep_prob <- 1.0
      gp_comp <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp) && !is.null(gp_comp$weights)) {
        kp <- gp_comp$weights[[outcome_label]] %||% gp_comp$weights[[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob <- as.numeric(kp)
      }
      return(as.numeric(keep_prob * cached))
    }
    # If no family detected, try OR-of-guarded-ANDs with shared gate member inside this label
    or_guard_info <- .detect_or_guard_shared_k1(expr, prep)
    if (!is.null(or_guard_info)) {
      comp_labels_ex <- setdiff(names(outcome_defs), outcome_label)
      if (length(comp_labels_ex) > 0) {
        # Exclude alias-of, special, and siblings mapping to same observable target
        this_target <- options$map_outcome_to %||% outcome_label
        comp_labels_ex <- Filter(function(lbl) {
          if (!is.null(outcome_defs[[lbl]]$options$alias_of)) return(FALSE)
          expr_lbl <- outcome_defs[[lbl]]$expr
          if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
            src <- expr_lbl$source
            if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
          }
          comp_target <- outcome_defs[[lbl]]$options$map_outcome_to %||% lbl
          if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
          TRUE
        }, comp_labels_ex)
      }
      comp_exprs_ex <- if (length(comp_labels_ex) > 0) lapply(comp_labels_ex, function(lbl) outcome_defs[[lbl]]$expr) else list()
      integ <- .or_guard_shared_k1_integrand(prep, component, or_guard_info, competitor_exprs = comp_exprs_ex)
      keep_prob <- 1.0
      gp_comp <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp) && !is.null(gp_comp$weights)) {
        kp <- gp_comp$weights[[outcome_label]] %||% gp_comp$weights[[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob <- as.numeric(kp)
      }
      upper <- if (is.finite(deadline)) deadline else Inf
      res <- tryCatch(stats::integrate(integ, lower = 0, upper = upper,
                                       rel.tol = 1e-6, abs.tol = 1e-8, stop.on.error = FALSE)$value,
                      error = function(e) 0.0)
      return(as.numeric(keep_prob * res))
    }
    # If no family or intra-label shared-gate detected, proceed with generic expression integration path only
    {
      # Generic race integral with donor contributions
      integrand <- function(t) {
        vapply(t, function(tt) {
          if (!is.finite(tt) || tt < 0) return(0.0)
          base <- .eval_expr_likelihood(expr, tt, prep, component)
          if (base > 0 && length(competitor_exprs) > 0) {
            base <- base * prod(vapply(competitor_exprs, function(ce) .eval_expr_survival(ce, tt, prep, component), numeric(1)))
          }
          add <- 0.0
          if (length(donors) > 0) {
            for (d in donors) {
              dens_d <- .eval_expr_likelihood(d$expr, tt, prep, component)
              if (dens_d == 0) next
              donor_src <- NULL
              for (lblx in names(outcome_defs)) { if (identical(outcome_defs[[lblx]]$expr, d$expr)) { donor_src <- lblx; break } }
              d_comp_labels <- setdiff(names(outcome_defs), donor_src)
              d_comp_exprs <- if (length(d_comp_labels) > 0) lapply(d_comp_labels, function(lbl) outcome_defs[[lbl]]$expr) else list()
              surv_prod_d <- if (length(d_comp_exprs) > 0) prod(vapply(d_comp_exprs, function(ce) .eval_expr_survival(ce, tt, prep, component), numeric(1))) else 1.0
              add <- add + d$weight * dens_d * surv_prod_d
            }
          }
          base + add
        }, numeric(1))
      }
      res <- tryCatch(stats::integrate(integrand, lower = 0, upper = if (is.finite(deadline)) deadline else Inf,
                                        rel.tol = 1e-5, abs.tol = 1e-6, stop.on.error = FALSE)$value,
                      error = function(e) 0.0)
      return(as.numeric(res))
    }
  }
  
  # Race density at observed RT
  {
    # Prefer general shared-family fast path when available
    prep <- .ensure_shared_families(prep)
    fam_ids <- prep$shared_family_index[[outcome_label]] %||% integer(0)
    if (length(fam_ids) == 1) {
      fam <- prep$shared_families[[fam_ids[[1]]]]
      comp_labels_ex <- setdiff(names(outcome_defs), union(outcome_label, fam$labels))
      if (length(comp_labels_ex) > 0) {
        this_target <- options$map_outcome_to %||% outcome_label
        comp_labels_ex <- Filter(function(lbl) {
          if (!is.null(outcome_defs[[lbl]]$options$alias_of)) return(FALSE)
          expr_lbl <- outcome_defs[[lbl]]$expr
          if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
            src <- expr_lbl$source
            if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
          }
          comp_target <- outcome_defs[[lbl]]$options$map_outcome_to %||% lbl
          if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
          TRUE
        }, comp_labels_ex)
      }
      comp_exprs_ex <- if (length(comp_labels_ex) > 0) lapply(comp_labels_ex, function(lbl) outcome_defs[[lbl]]$expr) else list()
      integ <- .family_integrand_for_label(prep, component, fam, outcome_label, competitor_exprs = comp_exprs_ex)
      # Build a small cache over a time grid to accelerate repeated evaluations
      # when compute_loglik calls this many times.
      # Cache key per (component,label,family)
      cache_env <- prep$.__cache_family__ %||% new.env(parent = emptyenv())
      prep$.__cache_family__ <- cache_env
      cache_key <- paste(component, outcome_label, fam$key, sep = "|")
      entry <- cache_env[[cache_key]]
      if (is.null(entry)) {
        # determine an upper bound for practical support: use deadline or a heuristic
        dl <- .get_component_attr(prep, component, "deadline") %||% prep$default_deadline
        tmax <- if (is.finite(dl)) dl else 2.0
        # rough grid resolution
        ngrid <- 1024L
        xs <- seq(0, tmax, length.out = ngrid)
        ys <- integ(xs)
        # build spline-like linear interpolator
        interp <- .interp1_builder(xs, ys)
        entry <- list(xs = xs, ys = ys, f = interp)
        cache_env[[cache_key]] <- entry
      }
      base_val <- as.numeric(entry$f(rt))
      keep_prob_rt <- 1.0
      gp_comp_rt <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp_rt) && !is.null(gp_comp_rt$weights)) {
        kp <- gp_comp_rt$weights[[outcome_label]] %||% gp_comp_rt$weights[[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob_rt <- as.numeric(kp)
      }
      return(base_val * keep_prob_rt)
    }
    # Try intra-label OR-of-guarded-ANDs with shared gate member at RT
    or_guard_info <- .detect_or_guard_shared_k1(expr, prep)
    if (!is.null(or_guard_info)) {
      comp_labels_ex <- setdiff(names(outcome_defs), outcome_label)
      if (length(comp_labels_ex) > 0) {
        this_target <- options$map_outcome_to %||% outcome_label
        comp_labels_ex <- Filter(function(lbl) {
          if (!is.null(outcome_defs[[lbl]]$options$alias_of)) return(FALSE)
          expr_lbl <- outcome_defs[[lbl]]$expr
          if (!is.null(expr_lbl$kind) && identical(expr_lbl$kind, "event")) {
            src <- expr_lbl$source
            if (identical(src, "__GUESS__") || identical(src, "__DEADLINE__")) return(FALSE)
          }
          comp_target <- outcome_defs[[lbl]]$options$map_outcome_to %||% lbl
          if (identical(as.character(comp_target), as.character(this_target))) return(FALSE)
          TRUE
        }, comp_labels_ex)
      }
      comp_exprs_ex <- if (length(comp_labels_ex) > 0) lapply(comp_labels_ex, function(lbl) outcome_defs[[lbl]]$expr) else list()
      integ <- .or_guard_shared_k1_integrand(prep, component, or_guard_info, competitor_exprs = comp_exprs_ex)
      # Cache this RT density too
      cache_env <- prep$.__cache_or_guard_k1__ %||% new.env(parent = emptyenv())
      prep$.__cache_or_guard_k1__ <- cache_env
      cache_key <- paste(component, outcome_label, sep = "|")
      entry <- cache_env[[cache_key]]
      if (is.null(entry)) {
        dl <- .get_component_attr(prep, component, "deadline") %||% prep$default_deadline
        tmax <- if (is.finite(dl)) dl else 2.0
        ngrid <- 1024L
        xs <- seq(0, tmax, length.out = ngrid)
        ys <- integ(xs)
        interp <- .interp1_builder(xs, ys)
        entry <- list(xs = xs, ys = ys, f = interp)
        cache_env[[cache_key]] <- entry
      }
      base_val <- as.numeric(entry$f(rt))
      keep_prob_rt <- 1.0
      gp_comp_rt <- .get_component_attr(prep, component, "guess")
      if (!is.null(gp_comp_rt) && !is.null(gp_comp_rt$weights)) {
        kp <- gp_comp_rt$weights[[outcome_label]] %||% gp_comp_rt$weights[[normalize_label(outcome_label)]] %||% NULL
        if (!is.null(kp)) keep_prob_rt <- as.numeric(kp)
      }
      return(base_val * keep_prob_rt)
    }
    dens_r <- .eval_expr_likelihood(expr, rt, prep, component)
    if (dens_r == 0) return(0.0)
    surv_comp <- if (length(competitor_exprs) > 0) prod(vapply(competitor_exprs, function(ce) .eval_expr_survival(ce, rt, prep, component), numeric(1))) else 1.0
    base_val <- dens_r * surv_comp
  }
  # Apply component-level keep probability to base density for observed labels
  keep_prob_rt <- 1.0
  gp_comp_rt <- .get_component_attr(prep, component, "guess")
  if (!is.null(gp_comp_rt) && !is.null(gp_comp_rt$weights)) {
    kp <- gp_comp_rt$weights[[outcome_label]] %||% gp_comp_rt$weights[[normalize_label(outcome_label)]] %||% NULL
    if (!is.null(kp)) keep_prob_rt <- as.numeric(kp)
  }
  base_val <- base_val * keep_prob_rt
  # Drop donor contributions here; component-level guess is handled via GUESS outcome
  base_val
}

# ==============================================================================
# PART 6: Top-level Likelihood Function
# ==============================================================================

#' Compute log-likelihood for data given a model
#'
#' @param model Model specification from new_API.R (race_spec or race_model_spec)
#' @param data Data frame with columns: outcome, rt, (optionally: component)
#' @param cache_env Optional environment for memoising per-trial likelihoods
#' @return Log-likelihood value with per-trial attributes
compute_loglik <- function(model, data, cache_env = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }

  if (!"outcome" %in% names(data) || !"rt" %in% names(data)) {
    stop("data must contain 'outcome' and 'rt' columns")
  }

  # Shared cache for memoisation (allows reuse across repeated evaluations)
  if (is.null(cache_env)) {
    cache_env <- new.env(parent = emptyenv())
  } else if (!is.environment(cache_env)) {
    stop("cache_env must be an environment or NULL")
  }

  # Prepare model
  prep <- .prepare_model_for_likelihood(model)

  # Check if mixture model and whether weights are parameters (auto-estimate)
  comp_info <- prep$components
  is_mixture <- length(comp_info$ids) > 1
  has_weight_param <- FALSE
  if (!is.null(comp_info$has_weight_param)) {
    has_weight_param <- any(comp_info$has_weight_param)
  }

  # Helper to compute likelihood for a single data row (with memoisation key)
  eval_row_likelihood <- function(idx) {
    outcome <- as.character(data$outcome[[idx]])
    rt_val <- as.numeric(data$rt[[idx]])
    rt_key <- if (is.na(rt_val)) "NA" else format(rt_val, digits = 22, scientific = FALSE, trim = TRUE)
    component_key <- ""

    if (is_mixture) {
      if (has_weight_param || !("component" %in% names(data))) {
        # Mixture-averaged likelihood (weights may be parameters)
        comp_ids <- comp_info$ids
        weights <- comp_info$weights
        total_lik <- 0.0
        for (j in seq_along(comp_ids)) {
          comp_lik <- .outcome_likelihood(outcome, rt_val, prep, comp_ids[[j]])
          total_lik <- total_lik + weights[[j]] * comp_lik
        }
        component_key <- "__mixture__"
        lik <- total_lik
      } else {
        component <- as.character(data$component[[idx]])
        component_key <- component
        lik <- .outcome_likelihood(outcome, rt_val, prep, component)
      }
    } else {
      component <- if (length(comp_info$ids) > 0) comp_info$ids[[1]] else "__default__"
      component_key <- component
      lik <- .outcome_likelihood(outcome, rt_val, prep, component)
    }

    key <- paste(component_key, outcome, rt_key, sep = "|")

    if (!is.null(cache_env[[key]])) {
      return(list(key = key, value = cache_env[[key]]))
    }

    cache_env[[key]] <- lik
    list(key = key, value = lik)
  }

  # Pre-evaluate unique rows to make reuse/caching explicit
  n_rows <- nrow(data)
  if (n_rows == 0) {
    total_ll <- 0.0
    attr(total_ll, "contributions") <- numeric(0)
    attr(total_ll, "log_contributions") <- numeric(0)
    attr(total_ll, "cache") <- cache_env
    attr(total_ll, "keys") <- character(0)
    return(total_ll)
  }

  raw_keys <- character(n_rows)
  for (i in seq_len(n_rows)) {
    outcome <- as.character(data$outcome[[i]])
    rt_val <- as.numeric(data$rt[[i]])
    rt_key <- if (is.na(rt_val)) "NA" else format(rt_val, digits = 22, scientific = FALSE, trim = TRUE)
    component_key <- ""
    if (is_mixture) {
      if (has_weight_param || !("component" %in% names(data))) {
        component_key <- "__mixture__"
      } else {
        component_key <- as.character(data$component[[i]])
      }
    } else {
      component_key <- if (length(comp_info$ids) > 0) comp_info$ids[[1]] else "__default__"
    }
    raw_keys[[i]] <- paste(component_key, outcome, rt_key, sep = "|")
  }

  unique_keys <- unique(raw_keys)
  key_to_indices <- split(seq_len(n_rows), raw_keys)

  for (key in unique_keys) {
    idx <- key_to_indices[[key]][[1]]
    res <- eval_row_likelihood(idx)
    # Ensure cache is populated even if likelihood computed downstream encountered Inf/NaN
    if (is.null(cache_env[[res$key]])) {
      cache_env[[res$key]] <- res$value
    }
  }

  lik_values <- vapply(raw_keys, function(key) cache_env[[key]], numeric(1))
  log_likes <- ifelse(is.finite(lik_values) & lik_values > 0, log(lik_values), -Inf)

  total_ll <- sum(log_likes)
  attr(total_ll, "contributions") <- lik_values
  attr(total_ll, "log_contributions") <- log_likes
  attr(total_ll, "keys") <- raw_keys
  attr(total_ll, "cache") <- cache_env

  total_ll
}

compute_loglik_from_tables <- function(model_tables, data, cache_env = NULL) {
  compute_loglik(model_tables, data, cache_env = cache_env)
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
response_probability <- function(model, outcome_label, component = NULL, prep = NULL) {
  if (is.null(prep)) {
    if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
      model <- tables_to_model(model)
    }
    prep <- .prepare_model_for_likelihood(model)
  }
  comp_ids <- prep$components$ids
  weights <- prep$components$weights
  
  # Helper for a single component
  .prob_for_component <- function(comp_id) {
    # Alias-of outcomes: sum of referenced labels
    out_def <- prep$outcomes[[outcome_label]]
    if (!is.null(out_def$options$alias_of)) {
      refs <- out_def$options$alias_of
      refs <- as.character(refs)
      vals <- vapply(refs, function(lbl) as.numeric(.outcome_likelihood(lbl, NA_real_, prep, comp_id)), numeric(1))
      return(sum(vals))
    }
    base <- as.numeric(.outcome_likelihood(outcome_label, NA_real_, prep, comp_id))
    # Apply component-level keep for non-special labels (complement goes to GUESS)
    if (!identical(outcome_label, "GUESS")) {
      gp <- .get_component_attr(prep, comp_id, "guess")
      if (!is.null(gp) && !is.null(gp$weights)) {
        keep <- gp$weights[[outcome_label]] %||% gp$weights[[normalize_label(outcome_label)]] %||% 1.0
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
response_probabilities <- function(model, component = NULL, prep = NULL) {
  if (is.null(prep)) {
    if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
      model <- tables_to_model(model)
    }
    prep <- .prepare_model_for_likelihood(model)
  }
  labels <- names(prep$outcomes)
  vals <- vapply(labels, function(lbl) response_probability(model, lbl, component = component, prep = prep), numeric(1))
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
  labels <- names(prep$outcomes)
  # Exclude alias-of outcomes from observed set
  labels <- Filter(function(lbl) {
    is.null(prep$outcomes[[lbl]]$options$alias_of)
  }, labels)
  # Compute base probabilities for all labels
  base <- vapply(labels, function(lbl) response_probability(model, lbl, component = component, prep = prep), numeric(1))
  names(base) <- labels
  # Partition into observable vs mapped-to-NA
  obs <- numeric(0)
  na_sum <- 0.0
  for (lbl in labels) {
    out_def <- prep$outcomes[[lbl]]
    map_to <- out_def$options$map_outcome_to %||% NULL
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
  obs
}
