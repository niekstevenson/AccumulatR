rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_prep.R")
base::source("R/pool_math.R")
base::source("R/likelihood_primitives.R")

approx_equal <- function(a, b, tol = 1e-10) {
  diff <- abs(a - b)
  diff <- diff[is.finite(diff)]
  if (length(diff) == 0L) return(TRUE)
  max(diff) <= tol
}

pool_coeffs_ref <- function(Svec, Fvec) {
  stopifnot(length(Svec) == length(Fvec))
  coeff <- 1.0
  for (i in seq_along(Svec)) {
    coeff <- c(coeff * Svec[i], 0) + c(0, coeff * Fvec[i])
  }
  as.numeric(coeff)
}

pool_density_fast_ref <- function(prep, pool_id, component, t) {
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
      pool_density_fast_ref(prep, mid, component, t)
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
        pool_survival_fast_ref(prep, oid, component, t)
      } else {
        1.0
      }
      Svec[[j]] <- surv_j
      Fvec[[j]] <- 1.0 - surv_j
    }
    coeffs <- pool_coeffs_ref(Svec, Fvec)
    if (k <= length(coeffs)) {
      total_density <- total_density + dens_i * coeffs[[k]]
    }
  }
  total_density
}

pool_survival_fast_ref <- function(prep, pool_id, component, t) {
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
      pool_survival_fast_ref(prep, mid, component, t)
    } else {
      1.0
    }
    Svec[[i]] <- surv
    Fvec[[i]] <- 1.0 - surv
  }
  coeffs <- pool_coeffs_ref(Svec, Fvec)
  upto <- min(length(coeffs), k)
  if (upto <= 0L) return(0.0)
  sum(coeffs[seq_len(upto)])
}

# ------------------------------------------------------------------------------
# Basic coefficient parity
# ------------------------------------------------------------------------------

Svec <- c(0.8, 0.65, 0.55)
Fvec <- 1.0 - Svec
coeff_cpp <- pool_coeffs(Svec, Fvec)
coeff_ref <- pool_coeffs_ref(Svec, Fvec)
if (!approx_equal(coeff_cpp, coeff_ref)) {
  stop("pool_coeffs_cpp mismatch detected")
}

# ------------------------------------------------------------------------------
# Fast-path parity on simple and nested pools
# ------------------------------------------------------------------------------

prep <- list()
prep$accumulators <- list(
  A = list(dist = "lognormal", params = list(meanlog = -0.1, sdlog = 0.4), onset = 0.0, q = 0.1),
  B = list(dist = "gamma", params = list(shape = 2.5, rate = 2.0), onset = 0.2, q = 0.05),
  C = list(dist = "exgauss", params = list(mu = 0.5, sigma = 0.3, tau = 0.8), onset = 0.15, q = 0.0)
)
prep$pools <- list(
  pool_k1 = list(members = c("A", "B", "C"), k = 1L),
  pool_k2 = list(members = c("A", "B", "C"), k = 2L),
  pool_nested = list(members = c("pool_k1", "C"), k = 1L)
)
all_ids <- c(names(prep$accumulators), names(prep$pools))
prep$.runtime <- list(
  id_index = setNames(seq_along(all_ids), all_ids)
)

times <- c(0.25, 0.5, 0.9, 1.3, 2.0)
for (tt in times) {
  val_new <- .pool_density_fast_value(prep, "pool_k1", NULL, tt)
  val_ref <- pool_density_fast_ref(prep, "pool_k1", NULL, tt)
  if (!approx_equal(val_new, val_ref, tol = 1e-10)) {
    stop(sprintf("pool_k1 density mismatch at t=%g", tt))
  }

  surv_new <- .pool_survival_fast_value(prep, "pool_k1", NULL, tt)
  surv_ref <- pool_survival_fast_ref(prep, "pool_k1", NULL, tt)
  if (!approx_equal(surv_new, surv_ref, tol = 1e-10)) {
    stop(sprintf("pool_k1 survival mismatch at t=%g", tt))
  }
}

for (tt in times) {
  val_new <- .pool_density_fast_value(prep, "pool_k2", NULL, tt)
  val_ref <- pool_density_fast_ref(prep, "pool_k2", NULL, tt)
  if (!approx_equal(val_new, val_ref, tol = 1e-10)) {
    stop(sprintf("pool_k2 density mismatch at t=%g", tt))
  }

  surv_new <- .pool_survival_fast_value(prep, "pool_k2", NULL, tt)
  surv_ref <- pool_survival_fast_ref(prep, "pool_k2", NULL, tt)
  if (!approx_equal(surv_new, surv_ref, tol = 1e-10)) {
    stop(sprintf("pool_k2 survival mismatch at t=%g", tt))
  }
}

for (tt in times) {
  val_new <- .pool_density_fast_value(prep, "pool_nested", NULL, tt)
  val_ref <- pool_density_fast_ref(prep, "pool_nested", NULL, tt)
  if (!approx_equal(val_new, val_ref, tol = 1e-10)) {
    stop(sprintf("pool_nested density mismatch at t=%g", tt))
  }

  surv_new <- .pool_survival_fast_value(prep, "pool_nested", NULL, tt)
  surv_ref <- pool_survival_fast_ref(prep, "pool_nested", NULL, tt)
  if (!approx_equal(surv_new, surv_ref, tol = 1e-10)) {
    stop(sprintf("pool_nested survival mismatch at t=%g", tt))
  }
}

# Forced scenario propagation retains weights and forced sets
forced_ids <- .labels_to_ids(prep, "A")
forced_density <- .pool_density(prep, "pool_k1", NULL, 0.7,
                                forced_complete = forced_ids)
scenarios_fc <- attr(forced_density, "scenarios")
if (!is.list(scenarios_fc)) {
  stop("Expected scenarios list for forced pool density")
}
if (length(scenarios_fc) > 0L) {
  weights <- vapply(scenarios_fc, `[[`, numeric(1), "weight")
  if (!approx_equal(sum(weights), forced_density)) {
    stop("Scenario weights do not sum to density under forced conditions")
  }
  contains_forced <- vapply(
    scenarios_fc,
    function(sc) all(forced_ids %in% as.integer(sc$forced_complete)),
    logical(1)
  )
  if (!all(contains_forced)) {
    stop("Forced complete IDs not propagated into scenarios")
  }
}

cat("pool_fast_cpp: ok\n")
