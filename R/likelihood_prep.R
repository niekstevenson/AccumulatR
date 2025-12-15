.refresh_compiled_prep_refs <- function(prep) {
  comp <- prep[[".expr_compiled"]] %||% list()
  nodes <- comp$nodes %||% list()
  if (length(nodes) == 0L) return(prep)
  update_env_prep <- function(fn) {
    if (!is.function(fn)) return()
    env <- environment(fn)
    if (is.null(env)) return()
    if (exists("prep", envir = env, inherits = FALSE)) {
      env$prep <- prep
    }
  }
  lapply(nodes, function(node) {
    if (!is.list(node)) return(NULL)
    lapply(c("density_fn", "surv_fn", "cdf_fn", "scenario_fn",
             "density_fast_fn", "surv_fast_fn", "cdf_fast_fn"),
           function(name) update_env_prep(node[[name]]))
    NULL
  })
  prep
}

.prepare_model_for_likelihood <- function(model) {
  if (exists("is_model_tables", mode = "function") && is_model_tables(model)) {
    model <- tables_to_model(model)
  }
  if (!exists("prepare_model", mode = "function")) {
    stop("prepare_model function not found - source generator_new.R first")
  }
  prep <- prepare_model(model)
  acc_ids <- names(prep[["accumulators"]] %||% list())
  pool_ids <- names(prep[["pools"]] %||% list())
  all_ids <- unique(c(acc_ids, pool_ids))
  prep[[".id_index"]] <- setNames(seq_along(all_ids), all_ids)
  prep[[".label_cache"]] <- new.env(parent = emptyenv(), hash = TRUE)
  prep <- .precompile_likelihood_expressions(prep)
  prep[[".competitors"]] <- .prepare_competitor_map(prep)
  runtime <- list(
    expr_compiled = prep[[".expr_compiled"]],
    label_cache = prep[[".label_cache"]],
    competitor_map = prep[[".competitors"]],
    id_index = prep[[".id_index"]],
    pool_members_cache = new.env(parent = emptyenv(), hash = TRUE),
    cache_bundle = .build_likelihood_cache_bundle(prep)
  )
  prep[[".runtime"]] <- runtime
  prep <- .refresh_compiled_prep_refs(prep)
  prep
}

.labels_to_ids <- function(prep, labels) {
  if (is.null(labels) || length(labels) == 0L) return(integer(0))
  idx_map <- .prep_id_index(prep)
  if (is.null(idx_map)) return(integer(0))
  chr <- as.character(labels)
  if (length(chr) == 0L) return(integer(0))
  cache_env <- .prep_label_cache(prep)
  cache_key <- NULL
  if (!is.null(cache_env)) {
    cache_key <- paste(chr, collapse = "\r")
    cached <- cache_env[[cache_key]]
    if (!is.null(cached)) return(cached)
  }
  vals <- idx_map[chr]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(integer(0))
  ids <- as.integer(vals)
  if (length(ids) > 1L) ids <- sort(unique(ids))
  if (!is.null(cache_env) && nzchar(cache_key %||% "")) {
    cache_env[[cache_key]] <- ids
  }
  ids
}
