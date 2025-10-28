# Helper utilities for memoizing likelihood computations.

lik_cache_create <- function(parent = emptyenv()) {
  env <- new.env(parent = parent, hash = TRUE)
  class(env) <- c("likelihood_cache", class(env))
  env
}

lik_cache_is <- function(x) inherits(x, "likelihood_cache")

lik_cache_resolve <- function(cache = NULL, prep = NULL) {
  if (is.environment(cache)) {
    if (!lik_cache_is(cache)) {
      class(cache) <- c("likelihood_cache", class(cache))
    }
    return(cache)
  }
  if (!is.null(prep)) {
    runtime <- prep[[".runtime"]]
    if (!is.null(runtime)) {
      existing <- runtime[["cache"]]
      if (is.environment(existing)) {
        if (!lik_cache_is(existing)) {
          class(existing) <- c("likelihood_cache", class(existing))
        }
        return(existing)
      }
    }
  }
  lik_cache_create()
}

lik_cache_time_key <- function(x) {
  if (length(x) == 0L) return(".")
  vapply(x, function(xx) {
    if (is.na(xx)) return("NA")
    if (!is.finite(xx)) {
      if (xx > 0) return("Inf")
      if (xx < 0) return("-Inf")
      return("NA")
    }
    sprintf("%.15g", xx)
  }, character(1), USE.NAMES = FALSE)[[1]]
}

lik_cache_forced_key <- function(ids) {
  if (length(ids) == 0L) return(".")
  ids <- as.integer(ids)
  if (length(ids) > 1L) ids <- sort(unique(ids))
  paste(ids, collapse = ",")
}

lik_cache_key <- function(kind, expr_id, component, t,
                          forced_complete = integer(0),
                          forced_survive = integer(0),
                          extra = NULL) {
  if (is.null(expr_id) || is.na(expr_id)) return(NULL)
  comp_key <- if (is.null(component) || length(component) == 0L) "__default__" else as.character(component)[[1]]
  t_key <- lik_cache_time_key(t)
  fc_key <- lik_cache_forced_key(forced_complete)
  fs_key <- lik_cache_forced_key(forced_survive)
  parts <- c(kind, expr_id, comp_key, t_key, fc_key, fs_key)
  if (!is.null(extra) && nzchar(extra)) {
    parts <- c(parts, extra)
  }
  paste(parts, collapse = "|")
}

lik_cache_get <- function(cache, key, default = NULL) {
  if (is.null(key) || !is.environment(cache)) return(default)
  val <- cache[[key]]
  if (is.null(val)) default else val
}

lik_cache_set <- function(cache, key, value) {
  if (is.null(key) || !is.environment(cache)) return(value)
  cache[[key]] <- value
  value
}
