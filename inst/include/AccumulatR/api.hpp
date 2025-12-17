#pragma once

#include <Rcpp.h>

namespace accumulatr {

inline Rcpp::List cpp_loglik(
    SEXP ctx_ptr,
    const Rcpp::List& structure,
    const Rcpp::NumericMatrix& param_matrix,
    const Rcpp::DataFrame& data_df,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["cpp_loglik"];
  return fn(ctx_ptr,
            structure,
            param_matrix,
            data_df,
            rel_tol,
            abs_tol,
            max_depth);
}

inline Rcpp::NumericVector cpp_loglik_multiple(
    SEXP ctx_ptr,
    const Rcpp::List& structure,
    const Rcpp::List& params_list,
    const Rcpp::DataFrame& data_df,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  Rcpp::Function fn = Rcpp::Environment::namespace_env("AccumulatR")["cpp_loglik_multiple"];
  return fn(ctx_ptr,
            structure,
            params_list,
            data_df,
            rel_tol,
            abs_tol,
            max_depth);
}

} // namespace accumulatr
