`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

.ensure_acc_param_t0 <- function(params) {
  if (is.null(params) || length(params) == 0L) {
    params <- list()
  }
  if (is.null(params$t0) || length(params$t0) == 0L) {
    params$t0 <- 0
  } else {
    params$t0 <- as.numeric(params$t0)[1]
  }
  params
}

normalize_label <- function(x) {
  idx <- is.na(x)
  out <- as.character(x)
  out[idx] <- "<NA>"
  out
}
