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

// One-trial simulation with reusable storage. The caller owns R's RNGScope;
// use independent processes, not worker threads, for parallel simulation.
enum class SimulationDistribution { Lognormal = 0, Gamma = 1, Exgauss = 2, LBA = 3, RDM = 4 };
struct SimulationAPI {
  void* (*prepare)(SEXP prep);
  void (*destroy)(void*);
  void* (*workspace)(const void*);
  void (*destroy_workspace)(void*);
  // Labels remain valid until the next draw; null denotes a missing response.
  int (*draw)(const void*, void*, const double*, int stride, int columns,
              int component, const double* onsets, const char** labels, double* times);
  double (*sample_time)(int distribution, const double*, int stride);
};
inline const SimulationAPI& simulation_api() {
  using Getter = const SimulationAPI*(*)();
  static const auto get = reinterpret_cast<Getter>(R_GetCCallable("AccumulatR", "simulation_api_1"));
  return *get();
}

} // namespace accumulatr
