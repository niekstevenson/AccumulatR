.model_view_named <- function(items, id_field, object_name) {
  items <- items %||% list()
  if (length(items) == 0L) {
    return(list())
  }
  ids <- names(items)
  if (is.null(ids) || length(ids) != length(items) ||
      any(is.na(ids) | !nzchar(ids))) {
    ids <- vapply(items, function(item) item[[id_field]] %||% "", character(1))
  }
  if (any(is.na(ids) | !nzchar(ids))) {
    stop(sprintf("Model view contains unnamed %s entries", object_name), call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop(sprintf("Model view contains duplicate %s ids", object_name), call. = FALSE)
  }
  setNames(items, ids)
}

.model_view <- function(model) {
  model_structure <- .as_model_structure(model)
  prep <- model_structure$prep %||% prepare_model(model_structure$model_spec)
  list(
    accumulators = .model_view_named(prep$accumulators, "id", "accumulator"),
    pools = .model_view_named(prep$pools, "id", "pool"),
    outcomes = .model_view_named(prep$outcomes, "label", "outcome")
  )
}

.model_view_accumulator_onset <- function(accumulator) {
  onset_spec <- accumulator$onset_spec %||% NULL
  if (is.list(onset_spec)) {
    kind <- onset_spec$kind %||% ""
    if (identical(kind, "absolute")) {
      return(as.numeric(onset_spec$value %||% 0)[1])
    }
    if (identical(kind, "after")) {
      return(list(
        kind = "after",
        source = as.character(onset_spec$source %||% ""),
        lag = as.numeric(onset_spec$lag %||% 0)[1]
      ))
    }
  }
  accumulator$onset %||% 0
}
