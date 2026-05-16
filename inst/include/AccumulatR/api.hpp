#pragma once

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

namespace accumulatr {

using loglik_trials_ccallable_t =
    void (*)(SEXP, SEXP, SEXP, SEXP, double, double *);

inline loglik_trials_ccallable_t loglik_trials_ccallable() {
  static loglik_trials_ccallable_t fn = nullptr;
  if (!fn) {
    fn = reinterpret_cast<loglik_trials_ccallable_t>(
        R_GetCCallable("AccumulatR", "loglik_trials"));
    if (!fn) {
      Rcpp::stop(
          "AccumulatR C-callable 'loglik_trials' not found (is AccumulatR loaded?)");
    }
  }
  return fn;
}

inline void loglik_trials(SEXP ctx_ptr,
                          SEXP param_matrix,
                          SEXP data_df,
                          SEXP ok,
                          double min_ll,
                          double *out) {
  loglik_trials_ccallable()(ctx_ptr,
                            param_matrix,
                            data_df,
                            ok,
                            min_ll,
                            out);
}

} // namespace accumulatr
