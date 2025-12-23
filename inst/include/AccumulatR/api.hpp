#pragma once

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

namespace accumulatr {

using cpp_loglik_ccallable_t = double (*)(SEXP, SEXP, SEXP, SEXP, SEXP, double, double, double, int);
using cpp_loglik_multiple_ccallable_t = SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP, double, double, double, int);

inline cpp_loglik_ccallable_t cpp_loglik_ccallable() {
  static cpp_loglik_ccallable_t fn = nullptr;
  if (!fn) {
    fn = reinterpret_cast<cpp_loglik_ccallable_t>(R_GetCCallable("AccumulatR", "cpp_loglik"));
    if (!fn) {
      Rcpp::stop("AccumulatR C-callable 'cpp_loglik' not found (is AccumulatR loaded?)");
    }
  }
  return fn;
}

inline cpp_loglik_multiple_ccallable_t cpp_loglik_multiple_ccallable() {
  static cpp_loglik_multiple_ccallable_t fn = nullptr;
  if (!fn) {
    fn = reinterpret_cast<cpp_loglik_multiple_ccallable_t>(
      R_GetCCallable("AccumulatR", "cpp_loglik_multiple"));
    if (!fn) {
      Rcpp::stop("AccumulatR C-callable 'cpp_loglik_multiple' not found (is AccumulatR loaded?)");
    }
  }
  return fn;
}

inline double cpp_loglik(
    SEXP ctx_ptr,
    const Rcpp::NumericMatrix& param_matrix,
    const Rcpp::DataFrame& data_df,
    const Rcpp::LogicalVector& ok,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  return cpp_loglik_ccallable()(ctx_ptr,
                                param_matrix,
                                data_df,
                                ok,
                                expand,
                                min_ll,
                                rel_tol,
                                abs_tol,
                                max_depth);
}

inline Rcpp::NumericVector cpp_loglik_multiple(
    SEXP ctx_ptr,
    const Rcpp::List& params_list,
    const Rcpp::DataFrame& data_df,
    const Rcpp::LogicalVector& ok,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    double rel_tol,
    double abs_tol,
    int max_depth) {
  SEXP out = cpp_loglik_multiple_ccallable()(ctx_ptr,
                                            params_list,
                                            data_df,
                                            ok,
                                            expand,
                                            min_ll,
                                            rel_tol,
                                            abs_tol,
                                            max_depth);
  return Rcpp::as<Rcpp::NumericVector>(out);
}

} // namespace accumulatr
