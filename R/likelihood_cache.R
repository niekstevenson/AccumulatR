# Evaluation state helpers for compiled likelihood execution.

.eval_state_create <- function(parent = emptyenv()) {
  new.env(parent = parent, hash = TRUE)
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
  if (!is.null(entry) || !isTRUE(init)) return(entry)
  entry <- new.env(parent = emptyenv())
  entry$density <- NULL
  entry$survival <- NULL
  entry$cdf <- NULL
  entry$scenarios <- NULL
  entry$meta <- list()
  state[[key]] <- entry
  entry
}

.eval_state_get_extra <- function(state, tag) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(NULL)
  state[[paste0("extra::", tag)]]
}

.eval_state_set_extra <- function(state, tag, value) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(value)
  state[[paste0("extra::", tag)]] <- value
  value
}
