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
  B = list(dist = "gamma", params = list(shape = 2.2, rate = 3.1), onset = 0.0, q = 0.0)
)

guard_expr <- list(
  kind = "guard",
  reference = list(kind = "event", source = "A"),
  blocker = list(kind = "event", source = "B"),
  unless = list()
)

prep <- list(
  accumulators = acc_defs,
  pools = list()
)

prep$outcomes <- list(
  guard = list(expr = guard_expr),
  eventA = list(expr = list(kind = "event", source = "A"))
)

prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)
prep$.runtime <- list(
  expr_compiled = prep[[".expr_compiled"]],
  label_cache = new.env(parent = emptyenv(), hash = TRUE),
  competitor_map = list(),
  id_index = setNames(seq_along(names(acc_defs)), names(acc_defs)),
  pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
  cache_bundle = .build_likelihood_cache_bundle(prep)
)
prep[[".id_index"]] <- prep$.runtime$id_index
prep[[".label_cache"]] <- prep$.runtime$label_cache

native_ctx <- .prep_native_context(prep)
native_prob <- .lik_native_fn("native_outcome_probability_cpp")

guard_id <- attr(prep$outcomes$guard$expr, ".lik_id", exact = TRUE)
event_id <- attr(prep$outcomes$eventA$expr, ".lik_id", exact = TRUE)
if (is.na(guard_id) || is.na(event_id)) stop("Missing compiled node ids")

upper_limits <- c(0.5, 1.0, 1.5)
for (upper in upper_limits) {
  native_val <- native_prob(
    native_ctx,
    as.integer(guard_id),
    as.numeric(upper),
    NULL,
    integer(0),
    integer(0),
    integer(0),
    .integrate_rel_tol(),
    .integrate_abs_tol(),
    12L
  )
  ref_val <- .integrate_outcome_probability(
    guard_expr,
    prep,
    NULL,
    upper_limit = upper,
    competitor_exprs = NULL
  )
  if (!isTRUE(all.equal(native_val, ref_val, tolerance = 5e-4))) {
    stop(sprintf("Guard outcome probability mismatch (upper=%g)", upper))
  }
}

# Event A with competitor B up to finite upper limit
competitors <- list(list(kind = "event", source = "B"))
comp_node <- .expr_lookup_compiled(competitors[[1]], prep)
if (is.null(comp_node)) stop("Competitor node missing")
for (upper in upper_limits) {
  native_val <- native_prob(
    native_ctx,
    as.integer(event_id),
    as.numeric(upper),
    NULL,
    integer(0),
    integer(0),
    as.integer(comp_node$id),
    .integrate_rel_tol(),
    .integrate_abs_tol(),
    12L
  )
  ref_val <- .integrate_outcome_probability(
    prep$outcomes$eventA$expr,
    prep,
    NULL,
    upper_limit = upper,
    competitor_exprs = competitors
  )
  if (!isTRUE(all.equal(native_val, ref_val, tolerance = 1e-6))) {
    stop(sprintf("Event outcome probability mismatch (upper=%g)", upper))
  }
}

cat("outcome_probability_cpp: ok\n")
