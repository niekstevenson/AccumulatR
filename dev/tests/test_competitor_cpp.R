rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_cache.R")
base::source("R/pool_math.R")
base::source("R/likelihood_primitives.R")
base::source("R/likelihood_prep.R")
base::source("R/likelihood_integrate.R")
base::source("R/likelihood_kernels.R")

acc_defs <- list(
  A = list(dist = "lognormal", params = list(meanlog = -0.2, sdlog = 0.45), onset = 0.1, q = 0.0),
  B = list(dist = "gamma", params = list(shape = 2.2, rate = 3.1), onset = 0.0, q = 0.0),
  C = list(dist = "exgauss", params = list(mu = 0.35, sigma = 0.25, tau = 0.6), onset = 0.05, q = 0.0)
)

guard_comp <- list(
  kind = "guard",
  reference = list(kind = "event", source = "B"),
  blocker = list(kind = "event", source = "C"),
  unless = list()
)

prep <- list(
  accumulators = acc_defs,
  pools = list()
)

prep$outcomes <- list(
  guard_comp = list(expr = guard_comp),
  event_B = list(expr = list(kind = "event", source = "B")),
  event_C = list(expr = list(kind = "event", source = "C")),
  event_A = list(expr = list(kind = "event", source = "A"))
)

prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)

id_index <- setNames(seq_along(names(acc_defs)), names(acc_defs))
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

ctx_ptr <- native_context_build(prep)
if (is.null(ctx_ptr) || !inherits(ctx_ptr, "externalptr")) {
  stop("Failed to build native context pointer")
}

expr_event <- function(id) list(kind = "event", source = id)
competitors_simple <- list(expr_event("B"), expr_event("C"))
competitors_guard <- list(guard_comp, expr_event("A"))

node_ids_simple <- vapply(competitors_simple, function(ex) {
  node <- .expr_lookup_compiled(ex, prep)
  if (is.null(node)) NA_integer_ else node$id
}, integer(1))

node_ids_guard <- vapply(competitors_guard, function(ex) {
  node <- .expr_lookup_compiled(ex, prep)
  if (is.null(node)) NA_integer_ else node$id
}, integer(1))

times <- c(0.2, 0.6, 1.0)
for (tt in times) {
  native_val <- native_competitor_survival_cpp(
    ctx_ptr,
    as.integer(node_ids_simple),
    as.numeric(tt),
    NULL
  )
  ref_val <- .compute_survival_product(
    outcome_expr = NULL,
    competitor_exprs = competitors_simple,
    prep = prep,
    component = NULL,
    t = tt
  )
  if (!isTRUE(all.equal(native_val, ref_val, tolerance = 1e-10))) {
    stop(sprintf("Simple competitor survival mismatch at t=%g", tt))
  }
}

for (tt in times) {
  native_val <- native_competitor_survival_cpp(
    ctx_ptr,
    as.integer(node_ids_guard),
    as.numeric(tt),
    NULL
  )
  ref_val <- .compute_survival_product(
    outcome_expr = NULL,
    competitor_exprs = competitors_guard,
    prep = prep,
    component = NULL,
    t = tt
  )
  if (!isTRUE(all.equal(native_val, ref_val, tolerance = 5e-4))) {
    stop(sprintf("Guard competitor survival mismatch at t=%g", tt))
  }
}

cat("competitor_cpp: ok\n")
