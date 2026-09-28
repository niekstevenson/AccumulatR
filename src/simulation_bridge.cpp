#include <Rcpp.h>

#include "simulation.hpp"

// [[Rcpp::export]]
Rcpp::DataFrame simulate_cpp(const Rcpp::List &prep, Rcpp::Environment cache,
                            const Rcpp::NumericMatrix &parameters,
                            const Rcpp::IntegerVector &component,
                            Rcpp::Nullable<Rcpp::NumericVector> onset,
                            const bool keep_detail, const bool keep_component) {
  using namespace accumulatr::simulation;
  SEXP native = cache.exists("native") ? static_cast<SEXP>(cache["native"]) : R_NilValue;
  // Serialization clears external pointers; rebuild once in each receiving worker.
  if (native == R_NilValue || R_ExternalPtrAddr(native) == nullptr) {
    Rcpp::XPtr<Program> compiled(new Program(prep), true);
    cache["native"] = compiled;
    native = compiled;
  }
  const auto &program = *Rcpp::XPtr<Program>(native);
  const auto &model = program.model;
  const int n_leaves = model.leaves.size();
  const int n_trials = parameters.nrow() / n_leaves;
  const int stride = parameters.nrow();
  const auto column_names = Rcpp::as<std::vector<std::string>>(
      Rcpp::colnames(parameters));
  std::vector<int> weight_columns;
  for (const auto &definition : model.components) {
    const auto found = std::find(column_names.begin(), column_names.end(), definition.weight_name);
    if (!definition.weight_name.empty() && found == column_names.end())
      Rcpp::stop("Missing mixture parameter '%s'", definition.weight_name);
    weight_columns.push_back(definition.weight_name.empty() ? -1 : found - column_names.begin());
  }
  const double *onsets = onset.isNotNull() ? REAL(onset.get()) : nullptr;
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
    if (i % 4096 == 0) Rcpp::checkUserInterrupt();
    const double *params = parameters.begin() + i * n_leaves;
    int chosen = 0;
    if (component[i] == NA_INTEGER) {
      if (model.components.size() > 1) {
        double sampled_total = 0.0;
        for (std::size_t c = 0; c < model.components.size(); ++c) {
          trial.weights[c] = weight_columns[c] < 0 ? model.components[c].weight
              : params[weight_columns[c] * stride];
          if (weight_columns[c] >= 0) sampled_total += trial.weights[c];
        }
        for (std::size_t c = 0; c < model.components.size(); ++c) {
          if (model.component_mode == "sample" && weight_columns[c] < 0)
            trial.weights[c] = 1.0 - sampled_total;
          if (!std::isfinite(trial.weights[c]) || trial.weights[c] < 0.0)
            Rcpp::stop("Mixture weights must be finite, non-negative, and sum to at most one");
        }
        chosen = choose(trial.weights);
      }
    } else {
      chosen = component[i] - 1;
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
      detail["ranked_outcomes"] = Rcpp::DataFrame::create(
          Rcpp::Named("rank") = ranks, Rcpp::Named("label") = labels,
          Rcpp::Named("time") = values);
      details[i] = detail;
    }
  }
  if (keep_component) output.push_back(components, "component");
  Rcpp::DataFrame result(output);
  if (keep_detail) result.attr("details") = details;
  return result;
}
