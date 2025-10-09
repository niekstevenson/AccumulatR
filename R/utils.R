normalize_label <- function(x) {
  idx <- is.na(x)
  out <- as.character(x)
  out[idx] <- "<NA>"
  out
}
