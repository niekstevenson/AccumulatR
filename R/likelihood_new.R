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
  prep
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

.build_pool_templates <- function(pool_id, members, k) {
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
      templates[[template_count]] <- list(
        finisher_idx = idx,
        complete_idx = if (length(combo_idx) > 0L) combo_idx else integer(0),
        survivor_idx = survivors,
        forced_complete_ids = c(members[[idx]], if (length(combo_idx) > 0L) members[combo_idx], pool_id),
        forced_survive_ids = if (length(survivors) > 0L) members[survivors] else character(0)
      )
      idx_entries[[j]] <- template_count
    }
    finisher_map[[idx]] <- idx_entries
  }
  list(templates = templates, finisher_map = finisher_map)
}

.append_unique <- function(base_set, additions, prep) {
  if (is.null(base_set) && (is.null(additions) || length(additions) == 0L)) {
    return(character(0))
  }
  combined <- c(base_set, additions)
  .normalize_forced_labels(combined, prep)
}

.normalize_forced_labels <- function(ids, prep) {
  if (is.null(ids) || length(ids) == 0L) return(character(0))
  chr <- as.character(ids)
  if (length(chr) == 0L) return(character(0))
  idx_map <- prep[[".id_index"]]
  if (is.null(idx_map)) {
    return(sort(unique(chr)))
  }
  ord <- idx_map[chr]
  na_mask <- is.na(ord)
  if (any(na_mask)) {
    max_ord <- if (all(is.na(ord[!na_mask]))) length(chr) else max(ord[!na_mask], na.rm = TRUE)
    ord[na_mask] <- max_ord + seq_len(sum(na_mask))
  }
  chr_unique <- chr[order(ord, chr)]
  unique(chr_unique)
}

.event_survival_at <- function(prep, id, component, t,
                               forced_complete = character(0),
                               forced_survive = character(0)) {
  if (id %in% forced_survive) return(1.0)
  if (id %in% forced_complete) return(0.0)
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
                          forced_complete = character(0),
                          forced_survive = character(0)) {
  1.0 - .event_survival_at(prep, id, component, t,
                           forced_complete = forced_complete,
                           forced_survive = forced_survive)
}

.event_density_at <- function(prep, id, component, t,
                              forced_complete = character(0),
                              forced_survive = character(0)) {
  if (id %in% forced_survive || id %in% forced_complete) return(0.0)
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
                          forced_complete = character(0),
                          forced_survive = character(0)) {
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
      forced_complete = as.character(fcomp),
      forced_survive = as.character(fsurv)
    )
  }

  fast_mode <- (length(forced_complete) == 0L && length(forced_survive) == 0L)
  fast_done <- FALSE
  if (fast_mode) {
    if (k == 1L) {
      if (n == 1L) {
        weight <- dens_vec[[1]]
        if (is.finite(weight) && weight > 0) {
          add_scenario(weight, c(members[[1]], pool_id), character(0))
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
          add_scenario(weight, c(members[[idx]], pool_id), members[survivors_idx])
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
        add_scenario(weight, c(members, pool_id), character(0))
        fast_done <- TRUE
      }
    } else {
      template_info <- .build_pool_templates(pool_id, members, k)
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
    forced_complete_flags <- if (length(forced_complete) > 0L) {
      !is.na(match(members, forced_complete))
    } else {
      rep(FALSE, n)
    }
    forced_survive_flags <- if (length(forced_survive) > 0L) {
      !is.na(match(members, forced_survive))
    } else {
      rep(FALSE, n)
    }
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

      forced_complete_ids <- .append_unique(forced_complete, c(members[[idx]], members[fc_idx], pool_id), prep)
      base_survive_ids <- .append_unique(forced_survive, members[fs_idx], prep)

      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .append_unique(base_survive_ids, members[survivors_idx], prep)
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
          .append_unique(forced_complete_ids, members[combo_idx], prep),
          .append_unique(base_survive_ids, members[survivors_idx], prep)
        )
      }
    }
  }

  if (!fast_done || !fast_mode) {
    others_idx_list <- lapply(idx_seq, function(i) idx_seq[idx_seq != i])
    forced_complete_flags <- if (length(forced_complete) > 0L) {
      !is.na(match(members, forced_complete))
    } else {
      rep(FALSE, n)
    }
    forced_survive_flags <- if (length(forced_survive) > 0L) {
      !is.na(match(members, forced_survive))
    } else {
      rep(FALSE, n)
    }
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

      forced_complete_ids <- .append_unique(forced_complete, c(members[[idx]], members[fc_idx], pool_id), prep)
      base_survive_ids <- .append_unique(forced_survive, members[fs_idx], prep)

      if (remaining_need == 0L) {
        survivors_idx <- free_idx
        weight <- base_factor
        if (length(survivors_idx) > 0L) weight <- weight * prod(surv_vec[survivors_idx])
        if (is.finite(weight) && weight > 0) {
          add_scenario(
            weight,
            forced_complete_ids,
            .append_unique(base_survive_ids, members[survivors_idx], prep)
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
          .append_unique(forced_complete_ids, members[combo_idx], prep),
          .append_unique(base_survive_ids, members[survivors_idx], prep)
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
                           forced_complete = character(0),
                           forced_survive = character(0)) {
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
#        forced_complete = character(),
#        forced_survive  = character())
.expr_scenarios_at <- function(expr, t, prep, component,
                               forced_complete = character(0),
                               forced_survive = character(0)) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(list())
  kind <- expr[["kind"]]
  forced_complete <- .normalize_forced_labels(forced_complete, prep)
  forced_survive <- .normalize_forced_labels(forced_survive, prep)

  make_scenario <- function(weight, fcomp, fsurv) {
    if (!is.finite(weight) || weight <= 0) return(NULL)
    list(
      weight = as.numeric(weight),
      forced_complete = as.character(fcomp),
      forced_survive = as.character(fsurv)
    )
  }

  if (identical(kind, "event")) {
    source_id <- expr[["source"]]
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
        fcomp <- .append_unique(forced_complete, c(source_id, psc$forced_complete), prep)
        fsurv <- .append_unique(forced_survive, psc$forced_survive, prep)
        sc <- make_scenario(w, fcomp, fsurv)
        if (!is.null(sc)) out[[length(out) + 1L]] <- sc
      }
      return(out)
    }
    sc <- make_scenario(
      weight,
      fcomp = .append_unique(forced_complete, source_id, prep),
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
        forced_survive = forced_survive)
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
            forced_complete = fcomp,
            forced_survive = fsurv)
          if (!is.finite(Fj) || Fj <= 0) {
            ok <- FALSE
            break
          }
          weight <- weight * Fj
          if (!is.finite(weight) || weight <= 0) {
            ok <- FALSE
            break
          }
          fcomp <- .append_unique(fcomp, .expr_sources(aj, prep), prep)
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
        forced_survive = forced_survive)
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
            fsurv <- .append_unique(fsurv, witness[[1]], prep)
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

  weight <- .eval_expr_likelihood(expr, t, prep, component)
  if (!is.finite(weight) || weight <= 0) return(list())
  sc <- make_scenario(
    weight,
    fcomp = .append_unique(forced_complete, .expr_sources(expr, prep), prep),
    fsurv = forced_survive
  )
  if (is.null(sc)) list() else list(sc)
}

.eval_expr_cdf_cond <- function(expr, t, prep, component,
                                forced_complete = character(0),
                                forced_survive = character(0)) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(0.0)
  kind <- expr[["kind"]]
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
                                forced_survive = forced_survive)
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
                                     forced_survive = forced_survive)
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
        .eval_expr_likelihood(expr, ui, prep, component)
      }, numeric(1))
    }
    val <- tryCatch(stats::integrate(dens_fun, lower = 0, upper = t,
                                     rel.tol = .integrate_rel_tol(), abs.tol = .integrate_abs_tol(),
                                     stop.on.error = FALSE)$value,
                    error = function(e) 0.0)
    out <- as.numeric(val)
    if (!is.finite(out)) out <- 0.0
    max(0.0, min(1.0, out))
  }
  .eval_expr_cdf(expr, t, prep, component)
}

.eval_expr_survival_cond <- function(expr, t, prep, component,
                                     forced_complete = character(0),
                                     forced_survive = character(0)) {
  if (is.null(expr) || is.null(expr[["kind"]])) return(1.0)
  kind <- expr[["kind"]]
  if (identical(kind, "event")) {
    src <- expr[["source"]]
    return(.event_survival_at(prep, src, component, t,
                              forced_complete = forced_complete,
                              forced_survive = forced_survive))
  }
  if (identical(kind, "and")) {
    return(1.0 - .eval_expr_cdf_cond(expr, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive))
  }
  if (identical(kind, "or")) {
    args <- expr[["args"]] %||% list()
    if (length(args) == 0) return(1.0)
    prod_val <- 1.0
    for (a in args) {
      Sa <- .eval_expr_survival_cond(a, t, prep, component,
                                     forced_complete = forced_complete,
                                     forced_survive = forced_survive)
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
        .eval_expr_likelihood(expr, uu, prep, component)
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
    max(0.0, min(1.0, out))
  }
.eval_expr_survival(expr, t, prep, component)
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
                         forced_complete = character(0),
                         forced_survive = character(0)))
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
                              forced_complete = character(0),
                              forced_survive = character(0)))
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
                             forced_complete = character(0),
                             forced_survive = character(0)))
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
                                               competitor_exprs = list()) {
  if (is.na(t) || !is.finite(t) || t < 0) {
    val <- 0.0
    attr(val, "scenarios") <- list()
    return(val)
  }
  simple_event <- function(e) {
    if (is.null(e) || is.null(e[['kind']]) || !identical(e[['kind']], "event")) return(FALSE)
    src <- e[['source']]
    !is.null(prep[["accumulators"]][[src]]) && is.null(prep[["pools"]][[src]])
  }
  fast_candidates <- function(exprs) {
    if (is.null(exprs) || length(exprs) == 0) return(TRUE)
    all(vapply(exprs, .expr_fastpath_supported, logical(1)))
  }
  if (.expr_fastpath_supported(expr) && fast_candidates(competitor_exprs)) {
    dens_fast <- .expr_density_fast(expr, t, prep, component)
    if (!is.finite(dens_fast) || dens_fast <= 0) {
      dens_fast <- max(0.0, dens_fast)
    } else if (length(competitor_exprs) > 0) {
      surv_fast <- 1.0
      for (comp_expr in competitor_exprs) {
        Sv <- .expr_survival_fast(comp_expr, t, prep, component)
        if (!is.finite(Sv) || Sv <= 0) {
          surv_fast <- 0.0
          break
        }
        surv_fast <- surv_fast * Sv
      }
      dens_fast <- dens_fast * surv_fast
    }
    if (!is.finite(dens_fast) || dens_fast <= 0) return(0.0)
    return(dens_fast)
  }
  if (simple_event(expr) && length(competitor_exprs) > 0 &&
      all(vapply(competitor_exprs, simple_event, logical(1)))) {
    src <- expr[['source']]
    dens <- .event_density_at(prep, src, component, t,
                              forced_complete = character(0),
                              forced_survive = character(0))
    if (dens <= 0) return(dens)
    val <- dens * .compute_survival_product(expr, competitor_exprs, prep, component, t)
    return(val)
  }
  scenarios <- .expr_scenarios_at(
    expr, t, prep, component,
    forced_complete = character(0),
    forced_survive = character(0))
  if (length(scenarios) == 0) {
    return(0.0)
  }
  total <- 0.0
  for (sc in scenarios) {
    if (is.null(sc) || sc$weight <= 0) next
    weight <- sc$weight
    if (length(competitor_exprs) > 0) {
      surv_prod <- 1.0
      for (comp_expr in competitor_exprs) {
        surv_val <- .eval_expr_survival_cond(
          comp_expr, t, prep, component,
          forced_complete = sc$forced_complete,
          forced_survive = sc$forced_survive)
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

.compute_survival_product <- function(outcome_expr, competitor_exprs, prep, component, t) {
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
      cluster_val <- prod(vapply(exprs, function(ce) .eval_expr_survival(ce, t, prep, component), numeric(1)))
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
                           forced_complete = character(0),
                           forced_survive = character(0)))
    }

    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = instant completion
      }
      return(.event_cdf_at(prep, source_id, component, t,
                           forced_complete = character(0),
                           forced_survive = character(0)))
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
                                forced_complete = character(0),
                                forced_survive = character(0)))
    }

    acc_def <- prep[["accumulators"]][[source_id]]
    if (!is.null(acc_def)) {
      comps <- acc_def[['components']]
      if (length(comps) > 0 && !is.null(component) &&
          !identical(component, "__default__") && !(component %in% comps)) {
        return(1.0)  # Not active = always survive (never fires)
      }
      return(.event_survival_at(prep, source_id, component, t,
                                forced_complete = character(0),
                                forced_survive = character(0)))
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
.eval_expr_likelihood <- function(expr, t, prep, component) {
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
                               forced_complete = character(0),
                               forced_survive = character(0)))
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
                               forced_complete = character(0),
                               forced_survive = character(0)))
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
      forced_complete = character(0),
      forced_survive = character(0))
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
      forced_complete = character(0),
      forced_survive = character(0))
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
      forced_complete = character(0),
      forced_survive = character(0))
    if (length(ref_scen) == 0) {
      val <- 0.0
      attr(val, "scenarios") <- list()
      return(val)
    }

    S_eff <- (function() {
      if (is.null(blocker)) return(1.0)
      if (length(unless_list) == 0) {
        return(.eval_expr_survival(blocker, t, prep, component))
      }
      f_blocker <- function(tt) .eval_expr_likelihood(blocker, tt, prep, component)
      S_prot_prod <- function(tt) {
        if (length(unless_list) == 0) return(rep(1.0, length(tt)))
        vapply(tt, function(ui) {
          vals <- vapply(unless_list, function(unl) .eval_expr_survival(unl, ui, prep, component), numeric(1))
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
    })()
    if (!is.finite(S_eff)) S_eff <- 0.0
    if (S_eff < 0) S_eff <- 0.0
    if (S_eff > 1) S_eff <- 1.0

    blocker_sources <- if (!is.null(blocker)) .expr_sources(blocker, prep) else character(0)
    scenarios <- list()
    for (sc in ref_scen) {
      if (is.null(sc)) next
      weight <- sc$weight * S_eff
      if (!is.finite(weight) || weight <= 0) next
      fcomp <- sc$forced_complete
      fsurv <- .append_unique(sc$forced_survive, blocker_sources, prep)
      scenarios[[length(scenarios) + 1L]] <- list(
        weight = as.numeric(weight),
        forced_complete = as.character(fcomp),
        forced_survive = as.character(fsurv)
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
.integrand_outcome_density <- function(t, expr, prep, component, competitor_exprs) {
  vapply(t, function(tt) {
    if (!is.finite(tt) || tt < 0) return(0.0)
    .scenario_density_with_competitors(expr, tt, prep, component,
                                       competitor_exprs = competitor_exprs)
  }, numeric(1))
}

# Scenario-conditioned probability of an expression over [0, upper_limit].
# competitor_exprs retained for backward compatibility but ignored.
.integrate_outcome_probability <- function(expr, prep, component, upper_limit = Inf, competitor_exprs = NULL) {
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
  if (is.null(expr) || is.null(expr[['kind']])) return(character(0))
  kind <- expr[['kind']]
  if (identical(kind, "event")) {
    source_id <- expr[['source']]
    if (!is.null(prep[["pools"]][[source_id]])) {
      members <- prep[["pools"]][[source_id]][['members']] %||% character(0)
      unique(unlist(lapply(members, function(mid) {
        if (!is.null(prep[["accumulators"]][[mid]])) return(mid)
        if (!is.null(prep[["pools"]][[mid]])) {
          return(.expr_sources(list(kind = "event", source = mid), prep))
        }
        character(0)
      })))
    } else {
      source_id
    }
  } else if (identical(kind, "and") || identical(kind, "or")) {
    unique(unlist(lapply(expr[['args']], .expr_sources, prep = prep)))
  } else if (identical(kind, "guard")) {
    ref_src <- .expr_sources(expr[['reference']], prep)
    blk_src <- .expr_sources(expr[['blocker']], prep)
    unless_src <- unique(unlist(lapply(expr[['unless']] %||% list(), .expr_sources, prep = prep)))
    unique(c(ref_src, blk_src, unless_src))
  } else if (identical(kind, "not")) {
    .expr_sources(expr[['arg']], prep)
  } else {
    character(0)
  }
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
            forced_complete = character(0),
            forced_survive = character(0)
          ))
        }
        if (!is.null(prep[["accumulators"]][[id]])) {
          return(.event_density_at(
            prep, id, component, ti,
            forced_complete = character(0),
            forced_survive = character(0)
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
            forced_complete = character(0),
            forced_survive = character(0)
          ))
        }
        if (!is.null(prep[["accumulators"]][[id]])) {
          return(.event_survival_at(
            prep, id, component, ti,
            forced_complete = character(0),
            forced_survive = character(0)
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
                                                           competitor_exprs = comp_exprs_guess)
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
                                                         competitor_exprs = comp_exprs_guess)
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
                                                    competitor_exprs = comp_exprs_map))
            } else {
              # Race density at rt for the mapped source outcome
              comp_exprs_map <- prep[['.competitors']][[label]] %||% list()
              dens_r <- .scenario_density_with_competitors(expr, rt, prep, component,
                                                           competitor_exprs = comp_exprs_map)
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
                          forced_complete = character(0),
                          forced_survive = character(0))
      }, numeric(1))
    }
  }
  get_S <- function(id) {
    function(tt) {
      vapply(tt, function(ti) {
        if (!is.finite(ti)) return(0.0)
        if (ti < 0) return(1.0)
        .event_survival_at(prep, id, component, ti,
                            forced_complete = character(0),
                            forced_survive = character(0))
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
                                                   competitor_exprs = competitor_exprs)
        if (base == 0) return(0.0)
        add <- 0.0
        if (length(donors) > 0) {
          for (d in donors) {
            dens_d <- .scenario_density_with_competitors(d[['expr']], tt, prep, component,
                                                          competitor_exprs = d[['competitors']] %||% list())
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
                                                 competitor_exprs = competitor_exprs)
    if (dens_r == 0) return(0.0)
    base_val <- dens_r
    # Include donor mass (e.g., timeout guesses) that keeps the observed RT
    if (length(donors) > 0) {
      donor_add <- 0.0
      for (d in donors) {
        if (identical(d[['rt_policy']], "na")) next
        dens_d <- .scenario_density_with_competitors(d[['expr']], rt, prep, component,
                                                     competitor_exprs = d[['competitors']] %||% list())
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

  for (i in seq_len(n_rows)) {
    outcome <- as.character(data[['outcome']][[i]])
    rt_val <- as.numeric(data[['rt']][[i]])
    rt_key <- .time_key(rt_val)
    component_key <- default_component
    if (is_mixture) {
      if (has_weight_param || !has_component_col) {
        component_key <- "__mixture__"
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
        raw_component <- data[['component']][[i]]
        if (is.na(raw_component)) {
          component_key <- default_component
        } else {
          component_chr <- as.character(raw_component)
          component_key <- if (nzchar(component_chr)) component_chr else default_component
        }
        lik <- .outcome_likelihood(outcome, rt_val, prep, component_key)
      }
    } else {
      lik <- .outcome_likelihood(outcome, rt_val, prep, component_key)
    }
    raw_keys[[i]] <- paste(component_key, outcome, rt_key, sep = "|")
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
