#pragma once

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

namespace accumulatr {

struct PreparedLikelihood;
struct LikelihoodWorkspace;

// Prepare and destroy on the R thread. A prepared likelihood owns its context
// and data; independent workspaces may evaluate it concurrently without R calls.
struct LikelihoodAPI {
  PreparedLikelihood *(*prepare)(SEXP context, SEXP data, int parameter_rows,
                                int parameter_columns, double min_ll);
  void (*destroy)(PreparedLikelihood *);
  LikelihoodWorkspace *(*workspace)(const PreparedLikelihood *);
  void (*destroy_workspace)(LikelihoodWorkspace *);
  // Column-major parameters; a null mask includes every compact trial.
  // C++ exceptions propagate to the caller; evaluation never invokes R.
  void (*evaluate)(const PreparedLikelihood *, LikelihoodWorkspace *,
                   const double *parameters, const int *valid_trials,
                   double *trial_loglik);
};

inline const LikelihoodAPI &likelihood_api() {
  using Getter = const LikelihoodAPI *(*)();
  static const auto get = reinterpret_cast<Getter>(
      R_GetCCallable("AccumulatR", "likelihood_api"));
  return *get();
}

} // namespace accumulatr
