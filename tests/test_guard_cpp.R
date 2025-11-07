rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_cache.R")
base::source("R/likelihood_prep.R")
base::source("R/pool_math.R")
base::source("R/likelihood_primitives.R")
base::source("R/likelihood_integrate.R")
base::source("R/likelihood_kernels.R")

acc_defs <- list(
  A = list(dist = "lognormal", params = list(meanlog = -0.2, sdlog = 0.45), onset = 0.1, q = 0.0),
  B = list(dist = "gamma", params = list(shape = 2.2, rate = 3.1), onset = 0.0, q = 0.0),
  C = list(dist = "exgauss", params = list(mu = 0.35, sigma = 0.25, tau = 0.6), onset = 0.05, q = 0.0)
)

guard_simple <- list(
  kind = "guard",
  reference = list(kind = "event", source = "A"),
  blocker = list(kind = "event", source = "B"),
  unless = list()
)

guard_protected <- list(
  kind = "guard",
  reference = list(kind = "event", source = "A"),
  blocker = list(kind = "event", source = "B"),
  unless = list(list(kind = "event", source = "C"))
)

prep <- list(
  accumulators = acc_defs,
  pools = list()
)
all_ids <- names(acc_defs)
id_index <- setNames(seq_along(all_ids), all_ids)

prep$outcomes <- list(
  guard_simple = list(expr = guard_simple),
  guard_protected = list(expr = guard_protected)
)
prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)
prep$.runtime <- list(
  expr_compiled = prep[[".expr_compiled"]],
  label_cache = new.env(parent = emptyenv(), hash = TRUE),
  competitor_map = list(),
  id_index = id_index,
  pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
  cache_bundle = .build_likelihood_cache_bundle(prep)
)
prep[[".id_index"]] <- id_index
prep[[".label_cache"]] <- prep$.runtime$label_cache

times <- c(0.25, 0.6, 1.1, 1.8)
for (tt in times) {
  surv_new <- .guard_effective_survival(guard_simple, tt, prep, NULL, integer(0), integer(0), .eval_state_create())
  surv_block <- .acc_survival(acc_defs$B, tt)
  if (!isTRUE(all.equal(surv_new, surv_block, tolerance = 1e-10))) {
    stop(sprintf("Simple guard survival mismatch at t=%g", tt))
  }
}

manual_protected_survival <- function(t) {
  if (t <= 0) return(1.0)
  integrand <- function(u) {
    .acc_density(acc_defs$B, u) * .acc_survival(acc_defs$C, u)
  }
  res <- tryCatch(
    stats::integrate(
      integrand,
      lower = 0,
      upper = t,
      rel.tol = 1e-8,
      abs.tol = 1e-10,
      stop.on.error = FALSE
    ),
    error = function(e) list(value = 0.0)
  )
  val <- 1.0 - as.numeric(res$value)
  if (!is.finite(val)) val <- 0.0
  max(0.0, min(1.0, val))
}

times_guard <- c(0.3, 0.8, 1.4)
native_guard <- .lik_native_fn("guard_effective_survival_cpp")
for (tt in times_guard) {
  integrand <- function(u) {
    .acc_density(acc_defs$B, u) * .acc_survival(acc_defs$C, u)
  }
  surv_cpp <- native_guard(integrand, tt, .integrate_rel_tol(), .integrate_abs_tol(), 12L)
  surv_manual <- manual_protected_survival(tt)
  if (abs(surv_cpp - surv_manual) > 5e-4) {
    stop(sprintf("Protected guard survival mismatch at t=%g (cpp=%.6f manual=%.6f)", tt, surv_cpp, surv_manual))
  }
}

cat("guard_cpp: ok\n")
