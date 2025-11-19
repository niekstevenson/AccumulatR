rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_cache.R")
base::source("R/pool_math.R")
base::source("R/likelihood_primitives.R")
base::source("R/likelihood_prep.R")
base::source("R/likelihood_integrate.R")
base::source("R/likelihood_kernels.R")

expr_event <- function(src) list(kind = "event", source = src)

acc_defs <- list(
  A = list(dist = "lognormal", params = list(meanlog = -0.2, sdlog = 0.4), onset = 0.1, q = 0.0),
  B = list(dist = "gamma", params = list(shape = 2.0, rate = 3.0), onset = 0.05, q = 0.0)
)

pool_defs <- list(
  P = list(members = c("A", "B"), k = 1L)
)

prep <- list(
  accumulators = acc_defs,
  pools = pool_defs
)

prep$outcomes <- list(
  eventA = list(expr = expr_event("A")),
  eventB = list(expr = expr_event("B")),
  poolP = list(expr = expr_event("P")),
  andAB = list(expr = list(kind = "and", args = list(expr_event("A"), expr_event("B"))))
)

prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)

id_labels <- c(names(acc_defs), names(pool_defs))
id_index <- setNames(seq_along(id_labels), id_labels)

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

node_ids <- vapply(prep$outcomes, function(outcome) attr(outcome$expr, ".lik_id", exact = TRUE), integer(1))
if (any(is.na(node_ids))) stop("Missing compiled node ids")

times <- c(0.25, 0.6, 1.2)

batches <- .native_node_batch_eval(
  prep,
  node_ids,
  times,
  component = NULL,
  forced_complete = integer(0),
  forced_survive = integer(0)
)

if (length(batches) != length(node_ids)) {
  stop(sprintf("Batch result length mismatch (got %d expected %d)", length(batches), length(node_ids)))
}

native_single <- native_node_eval_cpp
ctx_ptr <- .prep_native_context(prep)

for (idx in seq_along(node_ids)) {
  res_entry <- batches[[idx]]
  if (is.null(res_entry)) stop(sprintf("Missing batch entry %d", idx))
  evals <- res_entry$evaluations
  for (j in seq_along(times)) {
    single <- native_single(
      ctx_ptr,
      as.integer(node_ids[[idx]]),
      times[[j]],
      NULL,
      integer(0),
      integer(0)
    )
    batch <- evals[[j]]
    fields <- c("density", "survival", "cdf")
    for (fld in fields) {
      bval <- as.numeric(batch[[fld]])
      sval <- as.numeric(single[[fld]])
      if (!isTRUE(all.equal(bval, sval, tolerance = 1e-10))) {
        stop(sprintf("Field %s mismatch for node %d time %g", fld, node_ids[[idx]], times[[j]]))
      }
    }
  }
}

cat("node_batch_eval: ok\n")
