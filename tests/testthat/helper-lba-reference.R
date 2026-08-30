lba_denom_ref <- function(v, sv) {
  max(pnorm(v / sv), 1e-10)
}

lba_pdf_ref <- function(x, v, B, A, sv) {
  if (length(x) != 1L) {
    return(vapply(x, lba_pdf_ref, numeric(1), v = v, B = B, A = A, sv = sv))
  }
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) ||
      !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom_ref(v, sv)
  if (A <= 1e-10) {
    return(dnorm(B / x, mean = v, sd = sv) * B / (x * x * denom))
  }
  zs <- x * sv
  cmz <- B - x * v
  cz <- cmz / zs
  cz_max <- (cmz - A) / zs
  pdf <- (v * (pnorm(cz) - pnorm(cz_max)) +
    sv * (dnorm(cz_max) - dnorm(cz))) / (A * denom)
  if (is.finite(pdf) && pdf > 0) pdf else 0
}

lba_cdf_ref <- function(x, v, B, A, sv) {
  if (length(x) != 1L) {
    return(vapply(x, lba_cdf_ref, numeric(1), v = v, B = B, A = A, sv = sv))
  }
  if (!is.finite(x) || x <= 0 || !is.finite(v) || !is.finite(B) ||
      !is.finite(A) || !is.finite(sv) || sv <= 0) {
    return(0)
  }
  denom <- lba_denom_ref(v, sv)
  if (A <= 1e-10) {
    return(pnorm(B / x, mean = v, sd = sv, lower.tail = FALSE) / denom)
  }
  zs <- x * sv
  cmz <- B - x * v
  xx <- cmz - A
  cz <- cmz / zs
  cz_max <- xx / zs
  cdf <- (1 + (
    zs * (dnorm(cz_max) - dnorm(cz)) +
      xx * pnorm(cz_max) - cmz * pnorm(cz)
  ) / A) / denom
  min(max(cdf, 0), 1)
}
