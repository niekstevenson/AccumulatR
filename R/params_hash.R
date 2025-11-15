.param_hash_drop_cols <- c(
  "trial", "component", "accumulator_index", "acc_idx",
  "type", "role", "outcome", "rt", "condition",
  "component_weight", "params_hash"
)

.param_hash_order_cols <- c(
  "accumulator_id", "accumulator", "accumulator_index", "acc_idx"
)

.param_hash_canonical_rows <- function(rows) {
  if (is.null(rows) || nrow(rows) == 0L) return(NULL)
  keep_cols <- setdiff(names(rows), .param_hash_drop_cols)
  if (length(keep_cols) == 0L) return(NULL)
  canon <- rows[, keep_cols, drop = FALSE]
  if (ncol(canon) == 0L) return(NULL)
  for (col in names(canon)) {
    if (is.factor(canon[[col]])) {
      canon[[col]] <- as.character(canon[[col]])
    }
  }
  if ("params" %in% names(canon)) {
    canon$params <- lapply(canon$params, function(entry) {
      if (is.null(entry) || !is.list(entry)) return(entry)
      if (is.null(names(entry))) return(entry)
      entry[order(names(entry))]
    })
  }
  order_cols <- intersect(.param_hash_order_cols, names(rows))
  if (length(order_cols) > 0L) {
    ord <- do.call(order, lapply(order_cols, function(col) rows[[col]]))
    canon <- canon[ord, , drop = FALSE]
  }
  rownames(canon) <- NULL
  canon
}

.param_rows_hash_value <- function(rows) {
  canon <- .param_hash_canonical_rows(rows)
  if (is.null(canon)) return("__base__")
  raw_bytes <- serialize(canon, connection = NULL, ascii = TRUE)
  paste(as.character(raw_bytes), collapse = "")
}

.params_hash_attach <- function(params_df) {
  if (is.null(params_df)) return(params_df)
  params_df <- as.data.frame(params_df)
  n <- nrow(params_df)
  if (n == 0L) {
    if (!"params_hash" %in% names(params_df)) {
      params_df$params_hash <- character(0)
    }
    return(params_df)
  }
  if (!"trial" %in% names(params_df)) {
    params_df$trial <- 1L
  }
  if ("component" %in% names(params_df)) {
    comp_vals <- as.character(params_df$component)
  } else {
    comp_vals <- rep(NA_character_, n)
  }
  comp_vals[is.na(comp_vals) | !nzchar(comp_vals)] <- "__default__"
  trial_vals <- params_df$trial %||% seq_len(n)
  group_keys <- paste(trial_vals, comp_vals, sep = "|")
  params_hash <- if ("params_hash" %in% names(params_df)) {
    as.character(params_df$params_hash)
  } else {
    rep(NA_character_, n)
  }
  idx_list <- split(seq_len(n), group_keys)
  for (idx in idx_list) {
    rows <- params_df[idx, , drop = FALSE]
    rows$params_hash <- NULL
    hash_val <- .param_rows_hash_value(rows)
    params_hash[idx] <- hash_val
  }
  params_df$params_hash <- params_hash
  params_df
}
