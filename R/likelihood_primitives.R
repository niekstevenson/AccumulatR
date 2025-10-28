.acc_density <- function(acc_def, t) {
  if (!is.finite(t) || t < 0) return(0.0)
  if (t < acc_def[['onset']]) return(0.0)

  success_prob <- 1 - acc_def[['q']]
  if (success_prob <= 0) return(0.0)

  reg <- dist_registry(acc_def[['dist']])
  if (is.null(reg) || is.null(reg[['d']])) {
    stop(sprintf("No density function for distribution '%s'", acc_def[['dist']]))
  }

  t_adj <- t - acc_def[['onset']]
  dens <- reg[['d']](t_adj, acc_def[['params']])

  success_prob * dens
}

.acc_survival <- function(acc_def, t) {
  if (!is.finite(t)) return(0.0)
  if (t < 0) return(1.0)
  if (t < acc_def[['onset']]) return(1.0)

  reg <- dist_registry(acc_def[['dist']])
  if (is.null(reg) || is.null(reg[['p']])) {
    stop(sprintf("No CDF function for distribution '%s'", acc_def[['dist']]))
  }

  t_adj <- t - acc_def[['onset']]
  surv_underlying <- 1 - reg[['p']](t_adj, acc_def[['params']])

  acc_def[['q']] + (1 - acc_def[['q']]) * surv_underlying
}

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
    if (length(forced_survive) > 0L && any(forced_survive == id_idx)) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
    if (length(forced_complete) > 0L && any(forced_complete == id_idx)) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }
  }

  pool_def <- prep[["pools"]][[id]]
  if (!is.null(pool_def)) {
    return(.pool_density(
      prep, id, component, t,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    ))
  }

  acc_def <- prep[["accumulators"]][[id]]
  if (!is.null(acc_def)) {
    dens <- .acc_density(acc_def, t)
    attr(dens, "scenarios") <- list()
    return(dens)
  }

  val <- 0.0
  attr(val, "scenarios") <- list()
  val
}

.pool_density <- function(prep, pool_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0)) {
  pool_defs <- prep[["pools"]]
  pool_def <- pool_defs[[pool_id]]
  if (is.null(pool_def)) stop(sprintf("Unknown pool '%s'", pool_id))
  members <- .pool_active_members(prep, pool_id, component)
  n <- length(members)
  if (n == 0L) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  k <- as.integer(pool_def[["k"]] %||% 1L)
  if (k < 1L || k > n) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }

  forced_complete <- .coerce_forced_ids(prep, forced_complete)
  forced_survive <- .coerce_forced_ids(prep, forced_survive)

  member_ids <- .labels_to_ids(prep, members)
  pool_idx <- .label_to_id(prep, pool_id)
  templates_info <- .build_pool_templates(pool_id, members, member_ids, pool_idx, k)
  templates <- templates_info$templates
  finisher_map <- templates_info$finisher_map

  dens_vec <- numeric(n)
  cdf_vec <- numeric(n)
  surv_vec <- numeric(n)
  acc_defs <- prep[["accumulators"]]
  pool_defs <- prep[["pools"]]
  for (i in seq_len(n)) {
    mid <- members[[i]]
    dens_vec[[i]] <- if (!is.null(acc_defs[[mid]])) {
      .acc_density(acc_defs[[mid]], t)
    } else if (!is.null(pool_defs[[mid]])) {
      .pool_density(prep, mid, component, t,
                    forced_complete = forced_complete,
                    forced_survive = forced_survive)
    } else {
      0.0
    }
    cdf_vec[[i]] <- if (!is.null(acc_defs[[mid]])) {
      1.0 - .acc_survival(acc_defs[[mid]], t)
    } else if (!is.null(pool_defs[[mid]])) {
      .event_cdf_at(prep, mid, component, t,
                    forced_complete = forced_complete,
                    forced_survive = forced_survive)
    } else {
      0.0
    }
    surv_vec[[i]] <- 1.0 - cdf_vec[[i]]
  }

  scenarios <- list()
  add_scenario <- function(weight, fcomp, fsurv) {
    sc <- .make_scenario_record(prep, weight, fcomp, fsurv)
    if (!is.null(sc)) scenarios[[length(scenarios) + 1L]] <<- sc
  }

  idx_seq <- seq_len(n)
  for (template in templates) {
    idx <- template$finisher_idx
    others_idx_list <- finisher_map
    others_idx <- others_idx_list[[idx]]
    dens_mid <- dens_vec[[idx]]
    if (!is.finite(dens_mid) || dens_mid <= 0) next

    complete_idx <- template$complete_idx
    survivor_idx <- template$survivor_idx
    forced_complete_ids <- .forced_union(prep, forced_complete, template$forced_complete_ids)
    forced_survive_ids <- .forced_union(prep, forced_survive, template$forced_survive_ids)

    if (length(complete_idx) > 0L) dens_mid <- dens_mid * prod(cdf_vec[complete_idx])
    if (length(survivor_idx) > 0L) dens_mid <- dens_mid * prod(surv_vec[survivor_idx])
    if (!is.finite(dens_mid) || dens_mid <= 0) next

    add_scenario(dens_mid, forced_complete_ids, forced_survive_ids)
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
