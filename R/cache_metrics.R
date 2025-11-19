.likelihood_cache_metric_fields <- c(
  "struct_hits", "struct_misses",
  "na_hits", "na_misses",
  "scratch_hits", "scratch_misses"
)

.likelihood_cache_metrics_env <- local({
  env <- new.env(parent = emptyenv())
  for (fld in .likelihood_cache_metric_fields) env[[fld]] <- 0L
  env
})

.cache_metrics_inc <- function(field, value = 1L) {
  if (!is.character(field) || length(field) != 1L) return(invisible(NULL))
  current <- .likelihood_cache_metrics_env[[field]]
  if (is.null(current)) current <- 0L
  .likelihood_cache_metrics_env[[field]] <- current + as.integer(value)
  invisible(NULL)
}

likelihood_cache_reset_stats <- function(prep = NULL) {
  for (fld in .likelihood_cache_metric_fields) {
    .likelihood_cache_metrics_env[[fld]] <- 0L
  }
  if (!is.null(prep)) {
    native_ctx <- try(.prep_native_context(prep), silent = TRUE)
    if (!inherits(native_ctx, "try-error") && inherits(native_ctx, "externalptr")) {
      try(native_cache_reset_stats_cpp(native_ctx), silent = TRUE)
    }
  }
  invisible(NULL)
}

likelihood_cache_stats <- function(prep = NULL) {
  native_stats <- NULL
  if (!is.null(prep)) {
    native_ctx <- try(.prep_native_context(prep), silent = TRUE)
    if (!inherits(native_ctx, "try-error") && inherits(native_ctx, "externalptr")) {
      native_stats <- try(native_cache_stats_cpp(native_ctx), silent = TRUE)
      if (inherits(native_stats, "try-error")) native_stats <- NULL
    }
  }
  r_stats <- sapply(.likelihood_cache_metric_fields, function(fld) {
    val <- .likelihood_cache_metrics_env[[fld]]
    if (is.null(val)) 0L else as.integer(val)
  }, simplify = FALSE, USE.NAMES = TRUE)
  list(
    r = r_stats,
    native = native_stats
  )
}
