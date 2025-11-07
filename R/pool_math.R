# Pool polynomial coefficients and prefix sums

# coeffs[m+1] = Prob(exactly m finished by time t) for a pool A
# given vectors Svec = tilde_S_i(t), Fvec = tilde_F_i(t) over i in A
pool_coeffs <- function(Svec, Fvec) {
  stopifnot(length(Svec) == length(Fvec))
  native <- .lik_native_fn("pool_coeffs_cpp")
  as.numeric(native(as.numeric(Svec), as.numeric(Fvec)))
}
