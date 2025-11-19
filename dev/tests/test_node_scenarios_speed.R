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

strip_compiled_expr <- function(expr) {
  if (!is.list(expr)) {
    attr(expr, ".lik_id") <- NULL
    return(expr)
  }
  expr_copy <- expr
  attr(expr_copy, ".lik_id") <- NULL
  for (nm in names(expr_copy)) {
    expr_copy[[nm]] <- strip_compiled_expr(expr_copy[[nm]])
  }
  expr_copy
}

acc_defs <- list(
  A = list(dist = "lognormal", params = list(meanlog = -0.2, sdlog = 0.4), onset = 0.1, q = 0.0),
  B = list(dist = "gamma", params = list(shape = 2.0, rate = 3.0), onset = 0.05, q = 0.0),
  C = list(dist = "exgauss", params = list(mu = 0.4, sigma = 0.2, tau = 0.6), onset = 0.0, q = 0.0)
)

pool_defs <- list(
  P = list(members = c("A", "B", "C"), k = 2L)
)

guard_expr <- list(
  kind = "guard",
  reference = expr_event("A"),
  blocker = expr_event("B"),
  unless = list(expr_event("C"))
)

prep <- list(
  accumulators = acc_defs,
  pools = pool_defs
)

prep$outcomes <- list(
  eventA = list(expr = expr_event("A")),
  eventB = list(expr = expr_event("B")),
  eventC = list(expr = expr_event("C")),
  poolP = list(expr = expr_event("P")),
  andABC = list(expr = list(
    kind = "and",
    args = list(expr_event("A"), expr_event("B"), expr_event("C"))
  )),
  orABC = list(expr = list(
    kind = "or",
    args = list(expr_event("A"), expr_event("B"), expr_event("C"))
  )),
  guardAB = list(expr = guard_expr)
)

id_labels <- c(names(acc_defs), names(pool_defs))
id_index <- setNames(seq_along(id_labels), id_labels)

prep$.runtime <- list(
  expr_compiled = NULL,
  label_cache = new.env(parent = emptyenv(), hash = TRUE),
  competitor_map = list(),
  id_index = id_index,
  pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
  cache_bundle = NULL
)
prep[[".id_index"]] <- id_index
prep[[".label_cache"]] <- prep$.runtime$label_cache

prep <- .precompile_likelihood_expressions(prep)
prep <- .refresh_compiled_prep_refs(prep)

prep$.runtime$expr_compiled <- prep[[".expr_compiled"]]
prep$.runtime$cache_bundle <- .build_likelihood_cache_bundle(prep)
prep <- .refresh_compiled_prep_refs(prep)

native_ctx <- .prep_native_context(prep)
if (is.null(native_ctx) || !inherits(native_ctx, "externalptr")) {
  stop("Failed to acquire native context pointer")
}

node_ids <- vapply(
  prep$outcomes,
  function(outcome) attr(outcome$expr, ".lik_id", exact = TRUE),
  integer(1)
)
if (any(is.na(node_ids))) {
  stop("Missing compiled node ids for speed test")
}

expr_uncompiled <- lapply(prep$outcomes, function(outcome) strip_compiled_expr(outcome$expr))

forced_id <- function(label) {
  val <- prep$.runtime$id_index[[label]]
  if (is.null(val) || is.na(val)) return(integer(0))
  as.integer(val)
}

forced_sets <- list(
  list(fc = integer(0), fs = integer(0)),
  list(fc = forced_id("B"), fs = integer(0)),
  list(fc = integer(0), fs = forced_id("C"))
)

case_grid <- expand.grid(
  name = names(node_ids),
  forced = seq_along(forced_sets),
  stringsAsFactors = FALSE
)
times <- c(0.2, 0.5, 0.9, 1.3)

run_native <- function(iterations) {
  total_cases <- nrow(case_grid) * length(times)
  idx <- 1L
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(iterations)) {
    case_idx <- ((idx - 1L) %% nrow(case_grid)) + 1L
    time_idx <- ((idx - 1L) %% length(times)) + 1L
    idx <- idx + 1L
    name <- case_grid$name[[case_idx]]
    node_id <- node_ids[[name]]
    forced <- forced_sets[[case_grid$forced[[case_idx]]]]
    native_node_scenarios_cpp(
      native_ctx,
      as.integer(node_id),
      as.numeric(times[[time_idx]]),
      NULL,
      forced$fc,
      forced$fs
    )
  }
  proc.time()[["elapsed"]] - start
}

run_legacy <- function(iterations) {
  total_cases <- nrow(case_grid) * length(times)
  idx <- 1L
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(iterations)) {
    case_idx <- ((idx - 1L) %% nrow(case_grid)) + 1L
    time_idx <- ((idx - 1L) %% length(times)) + 1L
    idx <- idx + 1L
    name <- case_grid$name[[case_idx]]
    expr_ref <- expr_uncompiled[[name]]
    forced <- forced_sets[[case_grid$forced[[case_idx]]]]
    .expr_scenarios_at_slow(
      expr_ref,
      as.numeric(times[[time_idx]]),
      prep,
      NULL,
      forced_complete = forced$fc,
      forced_survive = forced$fs,
      state = NULL
    )
  }
  proc.time()[["elapsed"]] - start
}

benchmark_iterations <- 4000L
gc()
native_time <- run_native(benchmark_iterations)
gc()
legacy_time <- run_legacy(benchmark_iterations)

speedup <- legacy_time / native_time

cat(sprintf(
  "node_scenarios speed: native=%.3fs legacy=%.3fs (%.2fx faster)\n",
  native_time,
  legacy_time,
  speedup
))
