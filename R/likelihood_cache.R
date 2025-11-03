# Evaluation state helpers for compiled likelihood execution.

.eval_state_create <- function(parent = emptyenv()) {
  state <- new.env(parent = parent, hash = FALSE)
  state$entries <- new.env(parent = emptyenv(), hash = TRUE)
  state$extras <- new.env(parent = emptyenv(), hash = TRUE)
  state$time_cache <- new.env(parent = emptyenv(), hash = TRUE)
  state$time_counter <- 0L
  state$forced_cache <- new.env(parent = emptyenv(), hash = TRUE)
  state$forced_counter <- 0L
  state
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
  if (is.null(state) || !is.environment(state)) return(NULL)
  if (is.null(node_id) || is.na(node_id)) return(NULL)
  node_store <- .state_fetch_node_store(state, node_id, init = init)
  if (is.null(node_store)) return(NULL)
  comp_id <- .state_component_id(node_store, component, init = init)
  if (is.null(comp_id)) return(NULL)
  comp_store <- .state_fetch_component_store(node_store, comp_id, init = init)
  if (is.null(comp_store)) return(NULL)
  time_id <- if (!is.null(extra) && is.list(extra) && length(extra$time_id) == 1L) {
    extra$time_id
  } else {
    .eval_state_time_id(state, t)
  }
  time_store <- .state_fetch_time_store(comp_store, time_id, init = init)
  if (is.null(time_store)) return(NULL)
  forced_info <- .state_forced_key(state, forced_complete, forced_survive, init = init)
  if (is.null(forced_info)) return(NULL)
  forced_env <- time_store$forced
  handle <- .state_entry_handle(forced_env, forced_info$key, init = init)
  if (is.null(handle)) return(NULL)
  if (.state_entry_is_new(handle)) {
    .state_entry_set(handle, "meta", list(
      node_id = as.integer(node_id),
      component_id = comp_id,
      time_id = time_id,
      forced_complete_id = forced_info$fc_id,
      forced_survive_id = forced_info$fs_id
    ))
    attr(handle, "created") <- NULL
  }
  handle
}

.eval_state_get_extra <- function(state, tag) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(NULL)
  extras <- state$extras
  extras[[tag]]
}

.eval_state_set_extra <- function(state, tag, value) {
  if (is.null(state) || !is.environment(state) || !nzchar(tag)) return(value)
  extras <- state$extras
  extras[[tag]] <- value
  value
}

.state_fetch_node_store <- function(state, node_id, init = TRUE) {
  entries <- state$entries
  key <- as.character(as.integer(node_id))
  node_store <- entries[[key]]
  if (!is.null(node_store) || !isTRUE(init)) return(node_store)
  node_store <- new.env(parent = emptyenv(), hash = FALSE)
  node_store$components <- new.env(parent = emptyenv(), hash = TRUE)
  node_store$component_map <- new.env(parent = emptyenv(), hash = TRUE)
  node_store$component_counter <- 0L
  entries[[key]] <- node_store
  node_store
}

.state_component_id <- function(node_store, component, init = TRUE) {
  comp_key <- .eval_state_component_key(component)
  comp_map <- node_store$component_map
  id <- comp_map[[comp_key]]
  if (!is.null(id)) return(id)
  if (!isTRUE(init)) return(NULL)
  id <- as.integer(node_store$component_counter) + 1L
  node_store$component_counter <- id
  comp_map[[comp_key]] <- id
  node_store$components[[as.character(id)]] <- NULL
  id
}

.state_fetch_component_store <- function(node_store, comp_id, init = TRUE) {
  comp_env <- node_store$components
  comp_key <- as.character(comp_id)
  comp_store <- comp_env[[comp_key]]
  if (!is.null(comp_store) || !isTRUE(init)) return(comp_store)
  comp_store <- new.env(parent = emptyenv(), hash = FALSE)
  comp_store$times <- new.env(parent = emptyenv(), hash = TRUE)
  comp_env[[comp_key]] <- comp_store
  comp_store
}

.state_fetch_time_store <- function(comp_store, time_id, init = TRUE) {
  times_env <- comp_store$times
  time_key <- as.character(time_id)
  time_store <- times_env[[time_key]]
  if (!is.null(time_store) || !isTRUE(init)) return(time_store)
  time_store <- new.env(parent = emptyenv(), hash = FALSE)
  time_store$forced <- new.env(parent = emptyenv(), hash = TRUE)
  times_env[[time_key]] <- time_store
  time_store
}

.state_forced_id <- function(state, values, tag, init = TRUE) {
  if (length(values) == 0L) return(0L)
  vals <- sort(unique(as.integer(values)))
  key <- paste(tag, paste(vals, collapse = ","), sep = ":")
  cache <- state$forced_cache
  id <- cache[[key]]
  if (!is.null(id)) return(id)
  if (!isTRUE(init)) return(NULL)
  state$forced_counter <- as.integer(state$forced_counter) + 1L
  id <- state$forced_counter
  cache[[key]] <- id
  cache[[paste0("values#", id)]] <- vals
  id
}

.state_forced_key <- function(state, forced_complete, forced_survive, init = TRUE) {
  fc_id <- .state_forced_id(state, forced_complete, "fc", init = init)
  if (is.null(fc_id)) return(NULL)
  fs_id <- .state_forced_id(state, forced_survive, "fs", init = init)
  if (is.null(fs_id)) return(NULL)
  list(
    key = paste(fc_id, fs_id, sep = "|"),
    fc_id = fc_id,
    fs_id = fs_id
  )
}

.state_entry_proto <- function() {
  list(
    density = NULL,
    survival = NULL,
    cdf = NULL,
    scenarios = NULL,
    meta = NULL
  )
}

.state_entry_handle <- function(env, key, init = TRUE) {
  data <- env[[key]]
  created <- FALSE
  if (is.null(data)) {
    if (!isTRUE(init)) return(NULL)
    data <- .state_entry_proto()
    env[[key]] <- data
    created <- TRUE
  }
  handle <- structure(list(env = env, key = key), class = "state_entry")
  if (created) attr(handle, "created") <- TRUE
  handle
}

.state_entry_is_new <- function(entry) {
  if (is.null(entry)) return(FALSE)
  isTRUE(attr(entry, "created", exact = TRUE))
}

.state_entry_get <- function(entry, field) {
  if (is.null(entry)) return(NULL)
  data <- entry$env[[entry$key]]
  if (is.null(data)) return(NULL)
  data[[field]]
}

.state_entry_set <- function(entry, field, value) {
  if (is.null(entry)) return(value)
  data <- entry$env[[entry$key]]
  if (is.null(data)) data <- .state_entry_proto()
  data[[field]] <- value
  entry$env[[entry$key]] <- data
  value
}

.eval_state_time_id <- function(state, t) {
  if (is.null(state) || !is.environment(state)) {
    return(.eval_state_time_key(t))
  }
  key <- .eval_state_time_key(t)
  cache <- state$time_cache
  id <- cache[[key]]
  if (!is.null(id)) return(id)
  state$time_counter <- as.integer(state$time_counter) + 1L
  id <- state$time_counter
  cache[[key]] <- id
  id
}
