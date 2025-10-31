# Evaluation state helpers for compiled likelihood execution.

.eval_state_create <- function(parent = emptyenv()) {
  env <- new.env(parent = parent, hash = TRUE)
  # Bounded LRU controls (disabled unless option set)
  env[[".lru_tick"]] <- 0L
  env[[".lru_index"]] <- new.env(parent = emptyenv(), hash = TRUE)
  env[[".key_count"]] <- 0L
  env[[".max_entries"]] <- as.integer(getOption("uuber.eval_state_max_entries", NA_integer_))
  env
}

# Internal helpers for bounded LRU environments used by eval state
.eval_state_is_control_key <- function(key) {
  is.character(key) && length(key) == 1L && startsWith(key, ".")
}

.eval_state_max_entries <- function(state) {
  mx <- state[[".max_entries"]]
  if (is.null(mx)) mx <- NA_integer_
  as.integer(mx)
}

.eval_state_touch_key <- function(state, key) {
  lru <- state[[".lru_index"]]
  if (!is.environment(lru)) return(invisible(NULL))
  tick <- as.integer(state[[".lru_tick"]] %||% 0L) + 1L
  state[[".lru_tick"]] <- tick
  lru[[key]] <- tick
  invisible(NULL)
}

.eval_state_register_insert <- function(state, key, was_present) {
  if (isTRUE(was_present)) {
    .eval_state_touch_key(state, key)
    return(invisible(NULL))
  }
  cnt <- as.integer(state[[".key_count"]] %||% 0L) + 1L
  state[[".key_count"]] <- cnt
  .eval_state_touch_key(state, key)
  .eval_state_maybe_evict(state)
  invisible(NULL)
}

.eval_state_maybe_evict <- function(state) {
  max_entries <- .eval_state_max_entries(state)
  if (is.na(max_entries) || max_entries <= 0L) return(invisible(NULL))
  cnt <- as.integer(state[[".key_count"]] %||% 0L)
  if (cnt <= max_entries) return(invisible(NULL))
  lru <- state[[".lru_index"]]
  if (!is.environment(lru)) return(invisible(NULL))
  # Evict least-recently-used entries until within bound
  while (cnt > max_entries) {
    keys <- ls(envir = lru, all.names = TRUE)
    if (length(keys) == 0L) break
    ticks <- vapply(keys, function(k) as.integer(lru[[k]] %||% .Machine$integer.max), integer(1))
    victim <- keys[[which.min(ticks)]]
    # Remove from state and LRU index
    if (!.eval_state_is_control_key(victim)) {
      rm(list = victim, envir = state)
    }
    rm(list = victim, envir = lru)
    cnt <- cnt - 1L
    state[[".key_count"]] <- cnt
  }
  invisible(NULL)
}

.eval_state_component_key <- function(component) {
  if (is.null(component) || length(component) == 0L) {
    "__default__"
  } else {
    as.character(component)[[1]]
  }
}

.eval_state_time_key <- function(x) {
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

.eval_state_ids_key <- function(ids) {
  if (length(ids) == 0L) return(".")
  ids <- as.integer(ids)
  if (length(ids) > 1L) ids <- sort(unique(ids))
  paste(ids, collapse = ",")
}

.eval_state_key <- function(node_id, component, t,
                            forced_complete = integer(0),
                            forced_survive = integer(0),
                            extra = NULL) {
  if (is.null(node_id) || is.na(node_id)) return(NULL)
  parts <- c(
    as.integer(node_id),
    .eval_state_component_key(component),
    .eval_state_time_key(t),
    .eval_state_ids_key(forced_complete),
    .eval_state_ids_key(forced_survive)
  )
  if (!is.null(extra) && nzchar(extra)) {
    parts <- c(parts, extra)
  }
  paste(parts, collapse = "|")
}

.eval_state_entry <- function(state, node_id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0),
                              extra = NULL,
                              init = TRUE) {
  if (is.null(state) || !is.environment(state)) return(NULL)
  key <- .eval_state_key(node_id, component, t, forced_complete, forced_survive, extra)
  if (is.null(key)) return(NULL)
  entry <- state[[key]]
  if (!is.null(entry) || !isTRUE(init)) {
    if (!is.null(entry)) .eval_state_touch_key(state, key)
    return(entry)
  }
  entry <- new.env(parent = emptyenv())
  entry$density <- NULL
  entry$survival <- NULL
  entry$cdf <- NULL
  entry$scenarios <- NULL
  entry$meta <- list()
  was_present <- !is.null(state[[key]])
  state[[key]] <- entry
  .eval_state_register_insert(state, key, was_present)
  entry
}

.eval_state_get_extra <- function(state, tag) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(NULL)
  key <- paste0("extra::", tag)
  val <- state[[key]]
  if (!is.null(val)) .eval_state_touch_key(state, key)
  val
}

.eval_state_set_extra <- function(state, tag, value) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(value)
  key <- paste0("extra::", tag)
  was_present <- !is.null(state[[key]])
  state[[key]] <- value
  .eval_state_register_insert(state, key, was_present)
  value
}
