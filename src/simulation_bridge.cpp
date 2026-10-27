#include <Rcpp.h>

#include "simulation.hpp"
#include "../inst/include/AccumulatR/api.hpp"

// [[Rcpp::export(rng = false)]]
SEXP make_simulation_context_cpp() {
  return R_MakeExternalPtr(nullptr, R_NilValue, R_NilValue);
}

// [[Rcpp::export]]
SEXP simulate_cpp(SEXP prep, SEXP context,
                  SEXP parameters, SEXP component, SEXP onset,
                  const bool keep_detail, const bool keep_component) {
  using namespace accumulatr::simulation;
  auto *native = static_cast<Program *>(R_ExternalPtrAddr(context));
  // Serialization drops the native address and finalizer; initialize both together.
  if (native == nullptr) {
    native = new Program(Rcpp::List(prep));
    R_SetExternalPtrAddr(context, native);
    Rcpp::XPtr<Program>(context).setDeleteFinalizer();
  }
  const auto &program = *native;
  const auto &model = program.model;
  const int n_leaves = model.leaves.size();
  const int stride = Rf_nrows(parameters);
  const int n_trials = stride / n_leaves;
  const int weight_start = Rf_ncols(parameters) - program.weight_param_count;
  const double *parameter_values = REAL(parameters);
  const int *component_codes = component == R_NilValue ? nullptr : INTEGER(component);
  const double *onsets = onset == R_NilValue ? nullptr : REAL(onset);
  Rcpp::List output;
  Rcpp::IntegerVector trials(n_trials);
  std::iota(trials.begin(), trials.end(), 1);
  output.push_back(trials, "trials");
  std::vector<Rcpp::CharacterVector> responses;
  std::vector<Rcpp::NumericVector> times;
  for (int rank = 0; rank < program.readouts; ++rank) {
    responses.emplace_back(n_trials, NA_STRING);
    times.emplace_back(n_trials, NA_REAL);
    const auto suffix = rank == 0 ? "" : std::to_string(rank + 1);
    output.push_back(responses.back(), "R" + suffix);
    output.push_back(times.back(), "rt" + suffix);
  }
  Rcpp::CharacterVector components(keep_component ? n_trials : 0);
  Rcpp::List details(keep_detail ? n_trials : 0);
  Trial trial(program);
  for (int i = 0; i < n_trials; ++i) {
    const double *params = parameter_values + i * n_leaves;
    int chosen = trial.sample(params, stride, weight_start,
      component_codes ? component_codes[i] : NA_INTEGER, onsets ? onsets + i*n_leaves : nullptr);
    const auto& plan = program.components[chosen];
    if (keep_component) components[i] = model.components[chosen].id;
    if (keep_detail) details[i] = trial.detail(chosen);
    for (int rank = 0; rank < program.readouts; ++rank) {
      if (!trial.labels[rank].empty()) responses[rank][i] = trial.labels[rank];
      times[rank][i] = trial.times[rank];
    }
    if (keep_detail && plan.readouts > 1) {
      Rcpp::IntegerVector ranks(plan.readouts);
      Rcpp::CharacterVector labels(plan.readouts, NA_STRING);
      Rcpp::NumericVector values(plan.readouts, NA_REAL);
      for (int rank = 0; rank < plan.readouts; ++rank) {
        ranks[rank] = rank + 1;
        labels[rank] = responses[rank][i];
        values[rank] = times[rank][i];
      }
      Rcpp::List detail(details[i]);
      Rcpp::List ranked = Rcpp::List::create(
          Rcpp::Named("rank") = ranks, Rcpp::Named("label") = labels,
          Rcpp::Named("time") = values);
      ranked.attr("class") = "data.frame";
      ranked.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -plan.readouts);
      detail["ranked_outcomes"] = ranked;
      details[i] = detail;
    }
  }
  if (keep_component) output.push_back(components, "component");
  output.attr("class") = "data.frame";
  output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_trials);
  if (keep_detail) output.attr("details") = details;
  return output;
}

namespace {
const accumulatr::SimulationAPI* native_simulation_api() {
  using namespace accumulatr::simulation;
  static const accumulatr::SimulationAPI api{
    [](SEXP prep) -> void* { return new Program(Rcpp::List(prep)); },
    [](void* p) { delete static_cast<Program*>(p); },
    [](const void* p) -> void* { return new Trial(*static_cast<const Program*>(p)); },
    [](void* w) { delete static_cast<Trial*>(w); },
    [](const void* p, void* w, const double* parameters, int stride, int columns,
       int component, const double* onsets, const char** labels, double* times) {
      const auto& program = *static_cast<const Program*>(p);
      auto& trial = *static_cast<Trial*>(w);
      int chosen = trial.sample(parameters,stride,columns-program.weight_param_count,component,onsets);
      for (int rank = 0; rank < program.readouts; ++rank) {
        labels[rank] = trial.labels[rank].empty() ? nullptr : trial.labels[rank].c_str();
        times[rank] = trial.times[rank];
      }
      return chosen;
    },
    [](int distribution, const double* row, int stride) {
      return accumulatr::leaf::sample_time(static_cast<accumulatr::leaf::DistKind>(distribution),row,stride);
    }
  };
  return &api;
}
}
// [[Rcpp::init]]
void register_simulation_api(DllInfo* dll) {
  R_RegisterCCallable("AccumulatR", "simulation_api_1", reinterpret_cast<DL_FUNC>(&native_simulation_api));
}
