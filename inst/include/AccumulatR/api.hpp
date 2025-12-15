#pragma once

#include <Rcpp.h>

namespace accumulatr {

inline Rcpp::List cpp_loglik(
    SEXP ctx_ptr,
    const Rcpp::List& structure,
    const Rcpp::List& trial_entries,
    const Rcpp::Nullable<Rcpp::DataFrame>& component_weights,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["cpp_loglik"];
  return fn(ctx_ptr,
            structure,
            trial_entries,
            component_weights,
            rel_tol,
            abs_tol,
            max_depth);
}

inline Rcpp::NumericVector cpp_loglik_multiple(
    SEXP ctx_ptr,
    const Rcpp::List& structure,
    const Rcpp::List& entries,
    const Rcpp::Nullable<Rcpp::DataFrame>& component_weights,
    const Rcpp::List& params_list,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["cpp_loglik_multiple"];
  return fn(ctx_ptr,
            structure,
            entries,
            component_weights,
            params_list,
            rel_tol,
            abs_tol,
            max_depth);
}

inline Rcpp::List cpp_update_entries(
    SEXP ctx_ptr,
    const Rcpp::List& entries,
    SEXP params_obj) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["cpp_update_entries"];
  return fn(ctx_ptr, entries, params_obj);
}

} // namespace accumulatr
