#pragma once

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

namespace accumulatr {

using loglik_total_ccallable_t =
    double (*)(SEXP, SEXP, SEXP, SEXP, SEXP, double);

inline loglik_total_ccallable_t loglik_total_ccallable() {
  static loglik_total_ccallable_t fn = nullptr;
  if (!fn) {
    fn = reinterpret_cast<loglik_total_ccallable_t>(
        R_GetCCallable("AccumulatR", "loglik_total"));
    if (!fn) {
      Rcpp::stop(
          "AccumulatR C-callable 'loglik_total' not found (is AccumulatR loaded?)");
    }
  }
  return fn;
}

inline double loglik_total(SEXP ctx_ptr,
                           SEXP param_matrix,
                           SEXP data_df,
                           SEXP ok,
                           SEXP trial_weights,
                           double min_ll) {
  return loglik_total_ccallable()(ctx_ptr,
                                  param_matrix,
                                  data_df,
                                  ok,
                                  trial_weights,
                                  min_ll);
}

} // namespace accumulatr
