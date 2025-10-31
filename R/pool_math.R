# Pool polynomial coefficients and prefix sums

# coeffs[m+1] = Prob(exactly m finished by time t) for a pool A
# given vectors Svec = tilde_S_i(t), Fvec = tilde_F_i(t) over i in A
pool_coeffs <- function(Svec, Fvec) {
  stopifnot(length(Svec) == length(Fvec))
  coeff <- 1.0
  for (i in seq_along(Svec)) {
    coeff <- c(coeff * Svec[i], 0) + c(0, coeff * Fvec[i])
  }
  as.numeric(coeff)
}
