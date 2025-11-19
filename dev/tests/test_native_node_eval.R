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
  A = list(
    dist = "lognormal",
    params = list(meanlog = -0.2, sdlog = 0.45),
    onset = 0.1,
    q = 0.0
  ),
  B = list(
    dist = "gamma",
    params = list(shape = 2.0, rate = 3.0),
    onset = 0.0,
    q = 0.0
  )
)

pool_defs <- list(
  P = list(
    members = c("A", "B"),
    k = 1L
  )
)

prep <- list(
  accumulators = acc_defs,
  pools = pool_defs
)

prep$outcomes <- list(
  accA = list(expr = list(kind = "event", source = "A")),
  poolP = list(expr = list(kind = "event", source = "P")),
  andAB = list(expr = list(
    kind = "and",
    args = list(
      list(kind = "event", source = "A"),
      list(kind = "event", source = "B")
    )
  )),
  orAB = list(expr = list(
    kind = "or",
    args = list(
      list(kind = "event", source = "A"),
      list(kind = "event", source = "B")
    )
  )),
  notB = list(expr = list(
    kind = "not",
    arg = list(kind = "event", source = "B")
  ))
)

prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)

id_index <- setNames(
  seq_along(c(names(acc_defs), names(pool_defs))),
  c(names(acc_defs), names(pool_defs))
)

prep$.runtime <- list(
  expr_compiled = prep[[".expr_compiled"]],
  label_cache = new.env(parent = emptyenv(), hash = TRUE),
  competitor_map = list(),
  id_index = id_index,
  pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
  cache_bundle = .build_likelihood_cache_bundle(prep)
)
prep[[".id_index"]] <- prep$.runtime$id_index
prep[[".label_cache"]] <- prep$.runtime$label_cache
prep <- .refresh_compiled_prep_refs(prep)

ctx_ptr <- native_context_build(prep)
if (is.null(ctx_ptr) || !inherits(ctx_ptr, "externalptr")) {
  stop("Failed to build native context pointer")
}

node_ids <- vapply(
  prep$outcomes,
  function(outcome) {
    attr(outcome$expr, ".lik_id", exact = TRUE)
  },
  integer(1)
)
if (any(is.na(node_ids))) stop("Missing compiled node id")

compiled_nodes <- prep[[".expr_compiled"]]$nodes
ref_values <- function(node_id, t,
                       forced_complete = integer(0),
                       forced_survive = integer(0)) {
  node <- compiled_nodes[[node_id]]
  if (is.null(node)) stop(sprintf("Missing compiled node %d", node_id))
  list(
    density = .node_density(
      node, t, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    ),
    survival = .node_survival_cond(
      node, t, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    ),
    cdf = .node_cdf_cond(
      node, t, prep, NULL,
      forced_complete = forced_complete,
      forced_survive = forced_survive
    )
  )
}

compare_results <- function(native_res, ref_res, label, tt, tol = 1e-8) {
  fields <- c("density", "survival", "cdf")
  for (fld in fields) {
    native_val <- as.numeric(native_res[[fld]])
    ref_val <- as.numeric(ref_res[[fld]])
    if (!isTRUE(all.equal(native_val, ref_val, tolerance = tol))) {
      stop(sprintf(
        "%s mismatch at t=%g for %s (native=%.12f ref=%.12f)",
        fld, tt, label, native_val, ref_val
      ))
    }
  }
}

t_values <- c(0.2, 0.5, 1.0)
test_nodes <- c(
  accA = node_ids[["accA"]],
  poolP = node_ids[["poolP"]],
  andAB = node_ids[["andAB"]],
  orAB = node_ids[["orAB"]],
  notB = node_ids[["notB"]]
)

for (label in names(test_nodes)) {
  node_id <- as.integer(test_nodes[[label]])
  for (tt in t_values) {
    native_res <- native_node_eval_cpp(
      ctx_ptr,
      node_id,
      as.numeric(tt),
      NULL,
      integer(0),
      integer(0)
    )
    ref_res <- ref_values(node_id, tt)
    compare_results(native_res, ref_res, label, tt)
  }
}

forced_survive_vec <- as.integer(prep[[".id_index"]][["A"]])
if (length(forced_survive_vec) == 0L || is.na(forced_survive_vec)) {
  stop("Failed to resolve forced ID for component A")
}
forced_survive_vec <- forced_survive_vec[!is.na(forced_survive_vec)]

forced_time <- 0.65
native_forced <- native_node_eval_cpp(
  ctx_ptr,
  as.integer(node_ids[["orAB"]]),
  forced_time,
  NULL,
  integer(0),
  forced_survive_vec
)
ref_forced <- ref_values(
  node_ids[["orAB"]],
  forced_time,
  forced_survive = forced_survive_vec
)
compare_results(native_forced, ref_forced, "orAB (forced survive A)", forced_time)

cat("native_node_eval: ok\n")
