#include <Rcpp.h>

#include "simulation.hpp"

// [[Rcpp::export]]
SEXP simulate_cpp(SEXP prep, SEXP cache,
                  SEXP parameters, SEXP component, SEXP onset,
                  const bool keep_detail, const bool keep_component) {
  using namespace accumulatr::simulation;
  const SEXP native_symbol = Rf_install("native");
  SEXP native = Rf_findVarInFrame(cache, native_symbol);
  // Serialization clears external pointers; rebuild once in each receiving worker.
  if (native == R_UnboundValue || R_ExternalPtrAddr(native) == nullptr) {
    Rcpp::XPtr<Program> compiled(new Program(Rcpp::List(prep)), true);
    Rf_defineVar(native_symbol, compiled, cache);
    native = compiled;
  }
  const auto &program = *static_cast<const Program *>(R_ExternalPtrAddr(native));
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
    int chosen = 0;
    if (component_codes == nullptr || component_codes[i] == NA_INTEGER) {
      if (model.components.size() > 1) {
        double sampled_total = 0.0;
        for (std::size_t c = 0; c < model.components.size(); ++c) {
          const int weight_index = program.components[c].weight_param_index;
          trial.weights[c] = weight_index < 0 ? model.components[c].weight
              : params[(weight_start + weight_index) * stride];
          if (weight_index >= 0) sampled_total += trial.weights[c];
        }
        for (std::size_t c = 0; c < model.components.size(); ++c) {
          if (model.component_mode == "sample" && program.components[c].weight_param_index < 0)
            trial.weights[c] = 1.0 - sampled_total;
        }
        chosen = choose(trial.weights);
      }
    } else {
      chosen = component_codes[i] - 1;
    }
    const auto &plan = program.components[chosen];
    trial.run(plan, params, stride, onsets == nullptr ? nullptr : onsets + i * n_leaves);
    if (keep_component) components[i] = model.components[chosen].id;
    if (keep_detail) details[i] = trial.detail(chosen);
    if (trial.candidates.empty()) continue;
    const auto candidate_less = [&](const int a, const int b) {
      return plan.readouts == 1 ? earlier(trial.result(a), trial.result(b))
          : trial.result(a).time < trial.result(b).time ||
              (trial.result(a).time == trial.result(b).time && a < b);
    };
    const int count = std::min<int>(plan.readouts, trial.candidates.size());
    if (count == 1) {
      const auto winner = std::min_element(trial.candidates.begin(), trial.candidates.end(), candidate_less);
      std::iter_swap(trial.candidates.begin(), winner);
    } else {
      std::partial_sort(trial.candidates.begin(), trial.candidates.begin() + count,
                        trial.candidates.end(), candidate_less);
    }
    for (int rank = 0; rank < count; ++rank) {
      const auto winner = trial.candidates[rank];
      const auto &outcome = model.outcomes[winner];
      const auto &policy = program.readout[winner];
      std::string label = outcome.label;
      double time = trial.result(winner).time;
      if (!policy.weights.empty()) label = policy.labels[choose(policy.weights)];
      if (policy.missing_rt) time = NA_REAL;
      if (outcome.mapping.maps_to_missing) continue;
      if (!outcome.mapping.observed_label.empty()) label = outcome.mapping.observed_label;
      responses[rank][i] = label;
      times[rank][i] = time;
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
