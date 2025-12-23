#include <Rcpp.h>
#include <R_ext/Rdynload.h>

// Internal C++ entrypoints (implemented in src/distributions.cpp)
double cpp_loglik(SEXP ctxSEXP,
                  Rcpp::NumericMatrix params_mat,
                  Rcpp::DataFrame data_df,
                  Rcpp::LogicalVector ok,
                  Rcpp::IntegerVector expand,
                  double min_ll,
                  double rel_tol,
                  double abs_tol,
                  int max_depth);

Rcpp::NumericVector cpp_loglik_multiple(SEXP ctxSEXP,
                                        Rcpp::List params_list,
                                        Rcpp::DataFrame data_df,
                                        Rcpp::LogicalVector ok,
                                        Rcpp::IntegerVector expand,
                                        double min_ll,
                                        double rel_tol,
                                        double abs_tol,
                                        int max_depth);

extern "C" {

typedef double (*accumulatr_cpp_loglik_fn_t)(SEXP, SEXP, SEXP, SEXP, SEXP, double, double, double, int);
typedef SEXP (*accumulatr_cpp_loglik_multiple_fn_t)(SEXP, SEXP, SEXP, SEXP, SEXP, double, double, double, int);

double accumulatr_cpp_loglik_ccallable(SEXP ctxSEXP,
                                       SEXP params_matSEXP,
                                       SEXP data_dfSEXP,
                                       SEXP okSEXP,
                                       SEXP expandSEXP,
                                       double min_ll,
                                       double rel_tol,
                                       double abs_tol,
                                       int max_depth) {
  try {
    Rcpp::NumericMatrix params_mat(params_matSEXP);
    Rcpp::DataFrame data_df(data_dfSEXP);
    Rcpp::LogicalVector ok(okSEXP);
    Rcpp::IntegerVector expand(expandSEXP);
    return cpp_loglik(ctxSEXP, params_mat, data_df, ok, expand, min_ll, rel_tol, abs_tol, max_depth);
  } catch (const std::exception& e) {
    ::Rf_error("%s", e.what());
  } catch (...) {
    ::Rf_error("Unknown C++ exception in AccumulatR::cpp_loglik");
  }
  return NA_REAL;
}

SEXP accumulatr_cpp_loglik_multiple_ccallable(SEXP ctxSEXP,
                                              SEXP params_listSEXP,
                                              SEXP data_dfSEXP,
                                              SEXP okSEXP,
                                              SEXP expandSEXP,
                                              double min_ll,
                                              double rel_tol,
                                              double abs_tol,
                                              int max_depth) {
  BEGIN_RCPP
  Rcpp::List params_list(params_listSEXP);
  Rcpp::DataFrame data_df(data_dfSEXP);
  Rcpp::LogicalVector ok(okSEXP);
  Rcpp::IntegerVector expand(expandSEXP);
  Rcpp::NumericVector res =
    cpp_loglik_multiple(ctxSEXP, params_list, data_df, ok, expand, min_ll, rel_tol, abs_tol, max_depth);
  return res;
  END_RCPP
}

} // extern "C"

// [[Rcpp::export]]
SEXP register_ccallables_cpp() {
  static bool registered = false;
  if (registered) return R_NilValue;

  R_RegisterCCallable("AccumulatR",
                      "cpp_loglik",
                      reinterpret_cast<DL_FUNC>(accumulatr_cpp_loglik_ccallable));
  R_RegisterCCallable("AccumulatR",
                      "cpp_loglik_multiple",
                      reinterpret_cast<DL_FUNC>(accumulatr_cpp_loglik_multiple_ccallable));

  registered = true;
  return R_NilValue;
}
