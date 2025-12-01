# Lightweight stubs for the old scratch cache. All functions keep their
# signatures so callers continue to work, but no memoization is performed.

.eval_state_create <- function(parent = emptyenv()) {
  new.env(parent = parent, hash = FALSE)
}

.eval_state_component_key <- function(component) {
  if (is.null(component) || length(component) == 0L) "__default__" else as.character(component)[[1]]
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
  vals <- as.integer(ids)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(".")
  if (length(vals) > 1L) vals <- sort(unique(vals))
  paste(vals, collapse = ",")
}

.eval_state_entry <- function(state, node_id, component, t,
                              forced_complete = integer(0),
                              forced_survive = integer(0),
                              extra = NULL,
                              init = TRUE) {
  NULL
}

.eval_state_get_extra <- function(state, tag) {
  NULL
}

.eval_state_set_extra <- function(state, tag, value) {
  value
}

.state_entry_is_new <- function(entry) {
  FALSE
}

.state_entry_get <- function(entry, field) {
  NULL
}

.state_entry_set <- function(entry, field, value) {
  value
}

.eval_state_time_id <- function(state, t) {
  .eval_state_time_key(t)
}
