#pragma once

#include <Rcpp.h>

namespace accumulatr {

inline Rcpp::NumericVector native_loglik_from_buffer(
    SEXP ctx_ptr,
    const Rcpp::List& structure,
    const Rcpp::List& trial_entries,
    const Rcpp::List& params_list,
    const Rcpp::Nullable<Rcpp::DataFrame>& component_weights,
    double default_deadline,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["native_loglik_from_buffer_cpp"];
  return fn(ctx_ptr,
            structure,
            trial_entries,
            params_list,
            component_weights,
            default_deadline,
            rel_tol,
            abs_tol,
            max_depth);
}

} // namespace accumulatr
EOF
