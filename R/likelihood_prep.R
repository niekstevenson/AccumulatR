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
  if (length(all_ids) > 0) {
    idx_map <- prep[[".id_index"]]
    rev_labels <- rep(NA_character_, length(idx_map))
    rev_labels[as.integer(idx_map)] <- names(idx_map)
    prep[[".id_labels"]] <- rev_labels
  } else {
    prep[[".id_labels"]] <- character(0)
  }
  prep[[".label_cache"]] <- new.env(parent = emptyenv(), hash = TRUE)
  prep <- .precompile_likelihood_expressions(prep)
  prep[[".competitors"]] <- .prepare_competitor_map(prep)
  lik_cache <- lik_cache_create()
  runtime <- list(
    cache = lik_cache,
    expr_compiled = prep[[".expr_compiled"]],
    label_cache = prep[[".label_cache"]],
    competitor_map = prep[[".competitors"]],
    id_index = prep[[".id_index"]],
    id_labels = prep[[".id_labels"]]
  )
  prep[[".runtime"]] <- runtime
  prep
}

.label_to_id <- function(prep, label) {
  if (is.null(label) || length(label) == 0L) return(NA_integer_)
  idx_map <- .prep_id_index(prep)
  if (is.null(idx_map)) return(NA_integer_)
  val <- idx_map[[as.character(label)]]
  if (is.null(val) || is.na(val)) return(NA_integer_)
  as.integer(val)
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

.ids_to_labels <- function(prep, ids) {
  if (is.null(ids) || length(ids) == 0L) return(character(0))
  labels <- .prep_id_labels(prep)
  if (is.null(labels) || length(labels) == 0L) return(character(0))
  ids <- as.integer(ids)
  out <- labels[ids]
  out[!is.na(out)]
}

.coerce_forced_ids <- function(prep, values) {
  if (is.null(values) || length(values) == 0L) return(integer(0))
  if (is.integer(values)) {
    if (length(values) > 1L) return(sort(unique(values)))
    return(values)
  }
  .labels_to_ids(prep, values)
}

.forced_union <- function(prep, base_set, additions) {
  base_ids <- .coerce_forced_ids(prep, base_set)
  add_ids <- .coerce_forced_ids(prep, additions)
  if (length(add_ids) == 0L) return(base_ids)
  if (length(base_ids) == 0L) return(add_ids)
  out <- c(base_ids, add_ids)
  out <- as.integer(out)
  if (length(out) > 1L) out <- sort(unique(out))
  out
}
