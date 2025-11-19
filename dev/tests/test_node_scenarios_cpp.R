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
  B = list(dist = "gamma", params = list(shape = 2.0, rate = 3.0), onset = 0.05, q = 0.0),
  C = list(dist = "exgauss", params = list(mu = 0.4, sigma = 0.2, tau = 0.6), onset = 0.0, q = 0.0)
)

pool_defs <- list(
  P = list(members = c("A", "B"), k = 1L)
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
  poolP = list(expr = expr_event("P")),
  andAB = list(expr = list(
    kind = "and",
    args = list(expr_event("A"), expr_event("B"))
  )),
  orAB = list(expr = list(
    kind = "or",
    args = list(expr_event("A"), expr_event("B"))
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
  stop("Missing compiled node ids for scenarios test")
}

compiled_nodes <- prep[[".expr_compiled"]]$nodes

normalize_scenarios <- function(lst) {
  if (length(lst) == 0L) {
    return(list())
  }
  entries <- lapply(lst, function(sc) {
    fc <- as.integer(sc$forced_complete %||% integer(0))
    fc <- fc[!is.na(fc)]
    if (length(fc) > 1L) fc <- sort(unique(fc))
    fs <- as.integer(sc$forced_survive %||% integer(0))
    fs <- fs[!is.na(fs)]
    if (length(fs) > 1L) fs <- sort(unique(fs))
    list(
      weight = as.numeric(sc$weight),
      fc = fc,
      fs = fs,
      fc_key = if (length(fc) == 0L) "." else paste(fc, collapse = ","),
      fs_key = if (length(fs) == 0L) "." else paste(fs, collapse = ",")
    )
  })
  ord <- order(
    vapply(entries, function(x) x$fc_key, character(1)),
    vapply(entries, function(x) x$fs_key, character(1)),
    vapply(entries, function(x) sprintf("%.12f", x$weight), character(1))
  )
  entries[ord]
}

compare_scenarios <- function(native, legacy, label, tt, tol = 1e-9) {
  norm_native <- normalize_scenarios(native)
  norm_legacy <- normalize_scenarios(legacy)
  if (length(norm_native) != length(norm_legacy)) {
    stop(sprintf("Scenario count mismatch for %s at t=%g", label, tt))
  }
  for (idx in seq_along(norm_native)) {
    na <- norm_native[[idx]]
    nb <- norm_legacy[[idx]]
    if (!identical(na$fc, nb$fc)) {
      stop(sprintf("Forced-complete mismatch for %s at t=%g", label, tt))
    }
    if (!identical(na$fs, nb$fs)) {
      stop(sprintf("Forced-survive mismatch for %s at t=%g", label, tt))
    }
    if (!isTRUE(all.equal(na$weight, nb$weight, tolerance = tol))) {
      stop(sprintf("Weight mismatch for %s at t=%g (native=%.12f ref=%.12f)", label, tt, na$weight, nb$weight))
    }
  }
}

legacy_scenarios <- function(expr, t, forced_complete, forced_survive) {
  original_fn <- .prep_native_context
  assign(".prep_native_context", function(...) NULL, envir = .GlobalEnv)
  on.exit(assign(".prep_native_context", original_fn, envir = .GlobalEnv), add = TRUE)
  .node_scenarios_at(
    expr,
    t,
    prep,
    NULL,
    forced_complete = forced_complete,
    forced_survive = forced_survive,
    state = .eval_state_create()
  )
}

forced_id <- function(label) {
  val <- prep$.runtime$id_index[[label]]
  if (is.null(val) || is.na(val)) return(integer(0))
  as.integer(val)
}

forced_survive_A <- forced_id("A")

test_cases <- list(
  list(name = "eventA", expr = "eventA"),
  list(name = "poolP", expr = "poolP"),
  list(name = "andAB", expr = "andAB"),
  list(name = "orAB", expr = "orAB"),
  list(name = "orAB_forced", expr = "orAB", forced_survive = forced_survive_A),
  list(name = "guardAB", expr = "guardAB")
)

times <- c(0.25, 0.6, 1.1)

for (case in test_cases) {
  node_id <- node_ids[[case$expr]]
  if (is.na(node_id)) next
  node_obj <- compiled_nodes[[as.integer(node_id)]]
  if (is.null(node_obj)) {
    stop(sprintf("Missing compiled node object for %s", case$expr))
  }
  forced_complete <- case$forced_complete %||% integer(0)
  forced_survive <- case$forced_survive %||% integer(0)
  for (tt in times) {
    native <- native_node_scenarios_cpp(
      native_ctx,
      as.integer(node_id),
      as.numeric(tt),
      NULL,
      forced_complete,
      forced_survive
    )
    legacy <- legacy_scenarios(
      node_obj,
      as.numeric(tt),
      forced_complete,
      forced_survive
    )
    compare_scenarios(native, legacy, case$name, tt)
  }
}

cat("node_scenarios_cpp: ok\n")
