rm(list = ls())

base::source("R/dist.R")
base::source("R/likelihood_common.R")
base::source("R/likelihood_cache.R")
base::source("R/pool_math.R")
base::source("R/likelihood_primitives.R")
base::source("R/likelihood_prep.R")
base::source("R/likelihood_integrate.R")
base::source("R/likelihood_kernels.R")
base::source("R/generator_new.R")
base::source("R/likelihood_param_interface.R")

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

native_ctx <- .prep_native_context(prep)
serialize_fn <- .lik_native_fn("native_prep_serialize_cpp")
deserialize_fn <- .lik_native_fn("native_context_from_proto_cpp")

blob <- serialize_fn(prep)
if (!is.raw(blob) || length(blob) == 0L) {
  stop("Serialized prep blob is empty")
}

ctx_from_blob <- deserialize_fn(blob)

node_ids <- vapply(prep$outcomes, function(outcome) attr(outcome$expr, ".lik_id", exact = TRUE), integer(1))
node_ids <- node_ids[!is.na(node_ids)]
if (length(node_ids) == 0L) stop("No compiled nodes found")

times <- c(0.3, 0.8)
eval_fn <- .lik_native_fn("native_node_eval_cpp")

for (node_id in node_ids) {
  for (tt in times) {
    base_res <- eval_fn(native_ctx, as.integer(node_id), tt, NULL, integer(0), integer(0))
    blob_res <- eval_fn(ctx_from_blob, as.integer(node_id), tt, NULL, integer(0), integer(0))
    fields <- c("density", "survival", "cdf")
    for (fld in fields) {
      bval <- as.numeric(base_res[[fld]])
      rval <- as.numeric(blob_res[[fld]])
      if (!isTRUE(all.equal(bval, rval, tolerance = 1e-10))) {
        stop(sprintf("Mismatch in field %s for node %d time %g", fld, node_id, tt))
      }
    }
  }
}

if (!is.null(prep$.runtime$cache_bundle)) {
  prep$.runtime$cache_bundle$native_ctx$ptr <- NULL
  rebuilt_ctx <- .prep_native_context(prep)
  sample_node <- node_ids[[1]]
  sample_time <- times[[1]]
  rebuilt_res <- eval_fn(rebuilt_ctx, as.integer(sample_node), sample_time, NULL, integer(0), integer(0))
  base_res <- eval_fn(native_ctx, as.integer(sample_node), sample_time, NULL, integer(0), integer(0))
  fields <- c("density", "survival", "cdf")
  for (fld in fields) {
    if (!isTRUE(all.equal(as.numeric(rebuilt_res[[fld]]), as.numeric(base_res[[fld]]), tolerance = 1e-10))) {
      stop(sprintf("Rehydrated context mismatch for field %s", fld))
    }
  }
}

override_rows <- data.frame(
  trial = 1L,
  accumulator_id = "A",
  dist = "lognormal",
  onset = 0.25,
  q = 0.15,
  stringsAsFactors = FALSE
)
override_rows$params <- I(list(list(meanlog = -0.05, sdlog = 0.55)))

density_params_fn <- .lik_native_fn("native_density_with_competitors_params_cpp")
density_expected_fn <- .lik_native_fn("acc_density_cpp")
event_id <- attr(prep$outcomes$eventA$expr, ".lik_id", exact = TRUE)

param_density <- density_params_fn(
  native_ctx,
  as.integer(event_id),
  0.4,
  NULL,
  integer(0),
  integer(0),
  integer(0),
  override_rows
)$density

density_expected <- density_expected_fn(
  0.4,
  0.25,
  0.15,
  "lognormal",
  list(meanlog = -0.05, sdlog = 0.55)
)
if (!isTRUE(all.equal(as.numeric(density_expected), as.numeric(param_density), tolerance = 1e-10))) {
  stop("native_density_with_competitors_params_cpp did not match expected density")
}

expr_event <- prep$outcomes$eventA$expr
dens_trial_rows <- .scenario_density_with_competitors(
  expr_event,
  0.4,
  prep,
  NULL,
  state = NULL,
  trial_rows = override_rows
)
if (!isTRUE(all.equal(as.numeric(param_density), as.numeric(dens_trial_rows), tolerance = 1e-10))) {
  stop("scenario_density_with_competitors trial_rows path diverged from native params result")
}

prob_trial_rows <- .integrate_outcome_probability(
  expr_event,
  prep,
  NULL,
  upper_limit = 0.6,
  trial_rows = override_rows
)
prob_params <- .lik_native_fn("native_outcome_probability_params_cpp")(
  native_ctx,
  as.integer(event_id),
  0.6,
  NULL,
  integer(0),
  integer(0),
  integer(0),
  .integrate_rel_tol(),
  .integrate_abs_tol(),
  12L,
  override_rows
)
if (!isTRUE(all.equal(as.numeric(prob_params), as.numeric(prob_trial_rows), tolerance = 1e-10))) {
  stop("integrate_outcome_probability trial_rows path diverged from native params result")
}

cat("native_prep_serialization: ok\n")
