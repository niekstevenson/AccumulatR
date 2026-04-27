#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "eval_query.hpp"
#include "exact_sequence.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

inline Rcpp::NumericVector evaluate_loglik_records(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<ObservedRecord> &records,
    const double min_ll,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericVector out(records.size(), R_NegInf);
  if (records.empty()) {
    return out;
  }

  std::vector<EvalLoglikQuery> queries;
  queries.reserve(records.size());
  for (const auto &record : records) {
    if (record.backend != compile::BackendKind::Exact) {
      throw std::runtime_error(
          "observed likelihood received a non-exact record after exact planning");
    }
    queries.push_back(EvalLoglikQuery{
        record.trial_index,
        record.variant_index,
        record.semantic_code,
        record.observed_rt,
        record.row_map_index});
  }
  const auto loglik = evaluate_exact_loglik_queries_cached(
      exact_plans,
      layout,
      paramsSEXP,
      queries,
      min_ll,
      row_maps);
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    out[i] = loglik[i];
  }

  return out;
}

inline Rcpp::NumericVector evaluate_probability_records(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<ObservedRecord> &records,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericVector out(records.size(), 0.0);
  if (records.empty()) {
    return out;
  }

  std::vector<EvalProbabilityQuery> queries;
  queries.reserve(records.size());
  for (const auto &record : records) {
    if (record.backend != compile::BackendKind::Exact) {
      throw std::runtime_error(
          "observed probability received a non-exact record after exact planning");
    }
    queries.push_back(EvalProbabilityQuery{
        record.trial_index,
        record.variant_index,
        record.semantic_code,
        record.row_map_index});
  }
  const auto prob = evaluate_exact_probability_queries_cached(
      exact_plans,
      layout,
      paramsSEXP,
      queries,
      row_maps);
  for (R_xlen_t i = 0; i < prob.size(); ++i) {
    out[i] = prob[i];
  }

  return out;
}

inline double logsumexp_records(const std::vector<double> &values) {
  if (values.empty()) {
    return R_NegInf;
  }
  double anchor = R_NegInf;
  for (const auto value : values) {
    if (std::isfinite(value) && value > anchor) {
      anchor = value;
    }
  }
  if (!std::isfinite(anchor)) {
    return R_NegInf;
  }
  double sum = 0.0;
  for (const auto value : values) {
    if (std::isfinite(value)) {
      sum += std::exp(value - anchor);
    }
  }
  return sum > 0.0 ? anchor + std::log(sum) : R_NegInf;
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code);

struct NamedParamMatrix {
  const double *base{nullptr};
  int nrow{0};
  std::unordered_map<std::string, int> column_by_name;

  explicit NamedParamMatrix(SEXP paramsSEXP)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)) {
    const SEXP dimnames = Rf_getAttrib(paramsSEXP, R_DimNamesSymbol);
    if (TYPEOF(dimnames) != VECSXP || Rf_length(dimnames) < 2) {
      return;
    }
    const SEXP colnames = VECTOR_ELT(dimnames, 1);
    if (TYPEOF(colnames) != STRSXP) {
      return;
    }
    const auto n_cols = Rf_length(colnames);
    column_by_name.reserve(static_cast<std::size_t>(n_cols));
    for (int i = 0; i < n_cols; ++i) {
      if (STRING_ELT(colnames, i) == NA_STRING) {
        continue;
      }
      column_by_name.emplace(CHAR(STRING_ELT(colnames, i)), i);
    }
  }

  inline bool has(const std::string &name) const {
    return column_by_name.find(name) != column_by_name.end();
  }

  inline double value(const int row, const std::string &name) const {
    const auto it = column_by_name.find(name);
    if (it == column_by_name.end()) {
      return NA_REAL;
    }
    return base[static_cast<R_xlen_t>(it->second) * nrow + row];
  }
};

inline bool trials_have_observed_components(
    const PreparedTrialLayout &layout,
    const Rcpp::IntegerVector &component,
    const std::vector<unsigned char> *selected = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(selected, trial_index)) {
      continue;
    }
    if (integer_cell_is_na(
            component,
            static_cast<R_xlen_t>(layout.spans[trial_index].start_row))) {
      return false;
    }
  }
  return true;
}

inline std::vector<double> resolve_component_weights(
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &component_codes,
    const NamedParamMatrix &params,
    const int row) {
  std::vector<double> weights(component_codes.size(), 0.0);
  if (component_codes.empty()) {
    return weights;
  }

  if (model.component_mode != "sample") {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto weight = model.components[component_index].weight;
      if (!(std::isfinite(weight) && weight >= 0.0)) {
        throw std::runtime_error(
            "Component weights must be non-negative finite numbers");
      }
      weights[i] = weight;
    }
  } else {
    const auto reference_id =
        !model.component_reference.empty() && std::any_of(
                                                 component_codes.begin(),
                                                 component_codes.end(),
                                                 [&](const auto code) {
                                                   return model.components
                                                              [static_cast<std::size_t>(code - 1)]
                                                                  .id ==
                                                          model.component_reference;
                                                 })
            ? model.component_reference
            : model.components[static_cast<std::size_t>(component_codes.front() - 1)]
                  .id;
    double sum_nonref = 0.0;
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id == reference_id) {
        continue;
      }
      double weight = component.weight;
      if (!component.weight_name.empty() &&
          params.has(component.weight_name)) {
        weight = params.value(row, component.weight_name);
      }
      if (!(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0)) {
        throw std::runtime_error(
            "Mixture weight '" + component.weight_name +
            "' must be a probability in [0,1]");
      }
      weights[i] = weight;
      sum_nonref += weight;
    }
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      const auto &component = model.components[component_index];
      if (component.id != reference_id) {
        continue;
      }
      double ref_weight = 1.0 - sum_nonref;
      if (ref_weight < -1e-8) {
        throw std::runtime_error(
            "Mixture weights sum to >1; check non-reference weight params");
      }
      if (ref_weight < 0.0) {
        ref_weight = 0.0;
      }
      weights[i] = ref_weight;
      break;
    }
  }

  double total = 0.0;
  for (const auto weight : weights) {
    total += weight;
  }
  if (!(std::isfinite(total) && total > 0.0)) {
    const auto uniform = 1.0 / static_cast<double>(weights.size());
    std::fill(weights.begin(), weights.end(), uniform);
    return weights;
  }
  for (auto &weight : weights) {
    weight /= total;
  }
  return weights;
}

struct TrialComponentChoice {
  semantic::Index component_code{semantic::kInvalidIndex};
  double weight{0.0};
};

struct MissingAllAccumulator {
  std::vector<ObservedRecord> exact_records;
  std::vector<double> exact_component_weight_sum;

  explicit MissingAllAccumulator(const std::size_t n_trials)
      : exact_component_weight_sum(n_trials, 0.0) {}
};

inline std::unordered_map<std::string, int> build_trial_row_map(
    const Rcpp::CharacterVector &accumulator,
    const PreparedTrialSpan &span) {
  std::unordered_map<std::string, int> row_by_accumulator;
  row_by_accumulator.reserve(
      static_cast<std::size_t>(span.end_row - span.start_row + 1));
  for (auto row = static_cast<int>(span.start_row);
       row <= static_cast<int>(span.end_row);
       ++row) {
    row_by_accumulator.emplace(Rcpp::as<std::string>(accumulator[row]), row);
  }
  return row_by_accumulator;
}

inline std::vector<int> build_component_row_map(
    const std::unordered_map<std::string, int> &row_by_accumulator,
    const std::vector<std::string> &leaf_ids) {
  std::vector<int> row_map;
  row_map.reserve(leaf_ids.size());
  for (const auto &leaf_id : leaf_ids) {
    const auto it = row_by_accumulator.find(leaf_id);
    if (it == row_by_accumulator.end()) {
      throw std::runtime_error(
          "prepared data are missing accumulator row for '" + leaf_id + "'");
    }
    row_map.push_back(it->second);
  }
  return row_map;
}

inline std::vector<TrialComponentChoice> resolve_trial_components(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const semantic::SemanticModel &model,
    const PreparedTrialLayout &layout,
    const Rcpp::IntegerVector &component,
    const NamedParamMatrix &params,
    const std::size_t trial_index) {
  const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
  std::vector<TrialComponentChoice> resolved;

  if (!integer_cell_is_na(component, row)) {
    resolved.push_back(TrialComponentChoice{
        static_cast<semantic::Index>(component[row]),
        1.0});
    return resolved;
  }

  std::vector<semantic::Index> component_codes;
  const auto max_component_code = std::min(
      component_plans_by_code.size(),
      model.components.size() + 1U);
  for (std::size_t component_code = 1; component_code < max_component_code;
       ++component_code) {
    if (component_plans_by_code[component_code].present) {
      component_codes.push_back(static_cast<semantic::Index>(component_code));
    }
  }
  if (component_codes.empty()) {
    return resolved;
  }

  const auto weights = resolve_component_weights(
      model,
      component_codes,
      params,
      static_cast<int>(row));
  resolved.reserve(component_codes.size());
  for (std::size_t i = 0; i < component_codes.size(); ++i) {
    if (!(std::isfinite(weights[i]) && weights[i] > 0.0)) {
      continue;
    }
    resolved.push_back(TrialComponentChoice{component_codes[i], weights[i]});
  }
  return resolved;
}

inline void append_weighted_observed_branches(
    std::vector<ObservedRecord> *out,
    const semantic::Index trial_index,
    const semantic::Index variant_index,
    const compile::BackendKind backend,
    const double observed_rt,
    const double scale,
    const std::vector<ObservedBranch> &branches,
    const semantic::Index row_map_index) {
  if (!(std::isfinite(scale) && scale > 0.0)) {
    return;
  }
  for (const auto &branch : branches) {
    const double weight = scale * branch.weight;
    if (!(std::isfinite(weight) && weight > 0.0)) {
      continue;
    }
    out->push_back(ObservedRecord{
        trial_index,
        variant_index,
        backend,
        branch.semantic_code,
        observed_rt,
        weight,
        row_map_index});
  }
}

inline void append_trial_component_records(
    const semantic::Index trial_index,
    const ObservedTrialKind kind,
    const semantic::Index observed_label_code,
    const double observed_rt,
    const semantic::Index variant_index,
    const double component_weight,
    const semantic::Index row_map_index,
    const ComponentObservationPlan &component_plan,
    std::vector<ObservedRecord> *finite_records,
    std::vector<ObservedRecord> *missing_rt_records,
    MissingAllAccumulator *missing_all) {
  if (variant_index == semantic::kInvalidIndex ||
      !(std::isfinite(component_weight) && component_weight > 0.0)) {
    return;
  }

  switch (kind) {
  case ObservedTrialKind::Finite: {
    if (observed_label_code == semantic::kInvalidIndex) {
      return;
    }
    const auto observed_code = static_cast<std::size_t>(observed_label_code);
    if (observed_code >= component_plan.keep_by_code.size()) {
      return;
    }
    append_weighted_observed_branches(
        finite_records,
        trial_index,
        variant_index,
        component_plan.semantic_backend,
        observed_rt,
        component_weight,
        component_plan.keep_by_code[observed_code],
        row_map_index);
    return;
  }
  case ObservedTrialKind::MissingRt: {
    if (observed_label_code == semantic::kInvalidIndex) {
      return;
    }
    const auto observed_code = static_cast<std::size_t>(observed_label_code);
    if (observed_code >= component_plan.missing_rt_by_code.size()) {
      return;
    }
    append_weighted_observed_branches(
        missing_rt_records,
        trial_index,
        variant_index,
        component_plan.semantic_backend,
        observed_rt,
        component_weight,
        component_plan.missing_rt_by_code[observed_code],
        row_map_index);
    return;
  }
  case ObservedTrialKind::MissingAll:
    missing_all->exact_component_weight_sum[static_cast<std::size_t>(trial_index)] +=
        component_weight;
    append_weighted_observed_branches(
        &missing_all->exact_records,
        trial_index,
        variant_index,
        component_plan.semantic_backend,
        observed_rt,
        component_weight,
        component_plan.finite_observed_branches,
        row_map_index);
    return;
  }
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  const auto idx = static_cast<std::size_t>(component_code);
  return exact_variant_index_by_component_code[idx];
}

inline bool identity_trials_are_all_finite(
    const PreparedTrialLayout &layout,
    const Rcpp::IntegerVector &label,
    const Rcpp::NumericVector &rt,
    const std::vector<unsigned char> *selected = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(selected, trial_index)) {
      continue;
    }
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    if (integer_cell_is_na(label, row) || Rcpp::NumericVector::is_na(rt[row])) {
      return false;
    }
  }
  return true;
}

inline Rcpp::List evaluate_outcome_queries_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  Rcpp::DataFrame data(dataSEXP);
  const auto table = read_prepared_data_view(data);
  const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector probability(n_trials, 0.0);
  std::vector<ObservedRecord> records;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    if (integer_cell_is_na(table.component, row) || integer_cell_is_na(label, row)) {
      continue;
    }

    const auto component_code =
        static_cast<semantic::Index>(table.component[row]);
    if (component_code <= 0 ||
        component_code >=
            static_cast<semantic::Index>(component_plans_by_code.size())) {
      continue;
    }

    const auto &component_plan =
        component_plans_by_code[static_cast<std::size_t>(component_code)];
    if (!component_plan.present) {
      continue;
    }

    const auto variant_index = resolve_variant_index_by_component_code(
        component_code,
        exact_variant_index_by_component_code);
    if (variant_index == semantic::kInvalidIndex) {
      continue;
    }

    const auto observed_code =
        static_cast<std::size_t>(static_cast<semantic::Index>(label[row]));
    const auto append_records =
        [&](const std::vector<ObservedBranch> &branches) {
          for (const auto &branch : branches) {
            records.push_back(ObservedRecord{
                static_cast<semantic::Index>(trial_index),
                variant_index,
                component_plan.semantic_backend,
                branch.semantic_code,
                NA_REAL,
                branch.weight});
          }
        };

    if (observed_code < component_plan.keep_by_code.size()) {
      append_records(component_plan.keep_by_code[observed_code]);
    }
    if (observed_code < component_plan.missing_rt_by_code.size()) {
      append_records(component_plan.missing_rt_by_code[observed_code]);
    }
  }

  if (!records.empty()) {
    const auto donor_prob = evaluate_probability_records(
        exact_plans,
        layout,
        paramsSEXP,
        records);
    std::size_t i = 0;
    while (i < records.size()) {
      const auto trial_index = records[i].trial_index;
      double total = 0.0;
      while (i < records.size() && records[i].trial_index == trial_index) {
        total += records[i].weight * donor_prob[static_cast<R_xlen_t>(i)];
        ++i;
      }
      probability[static_cast<R_xlen_t>(trial_index)] =
          std::isfinite(total) && total > 0.0 ? total : 0.0;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = probability,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

inline SEXP evaluate_identity_trials_with_missing_all(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const std::vector<unsigned char> *selected = nullptr) {
  Rcpp::DataFrame data(dataSEXP);
  const auto table = read_prepared_data_view(data);
  const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector loglik(n_trials, min_ll);
  std::vector<ObservedRecord> finite_records;
  std::vector<ObservedRecord> missing_all_records;
  std::vector<bool> missing_all_trial(n_trials, false);

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    if (!trial_is_selected(selected, trial_index)) {
      continue;
    }
    const auto component_code = static_cast<semantic::Index>(table.component[row]);
    const auto &component_plan =
        component_plans_by_code[static_cast<std::size_t>(component_code)];
    const auto variant_index = resolve_variant_index_by_component_code(
        component_code,
        exact_variant_index_by_component_code);

    if (integer_cell_is_na(label, row)) {
      missing_all_trial[trial_index] = true;
      for (const auto &branch : component_plan.finite_observed_branches) {
        missing_all_records.push_back(ObservedRecord{
            static_cast<semantic::Index>(trial_index),
            variant_index,
            component_plan.semantic_backend,
            branch.semantic_code,
            NA_REAL,
            branch.weight});
      }
      continue;
    }

    finite_records.push_back(ObservedRecord{
        static_cast<semantic::Index>(trial_index),
        variant_index,
        component_plan.semantic_backend,
        static_cast<semantic::Index>(label[row]),
        rt[row],
        1.0});
  }

  if (!finite_records.empty()) {
    const auto donor_loglik = evaluate_loglik_records(
        exact_plans,
        layout,
        paramsSEXP,
        finite_records,
        min_ll);
    for (std::size_t i = 0; i < finite_records.size(); ++i) {
      const auto trial_index =
          static_cast<R_xlen_t>(finite_records[i].trial_index);
      loglik[trial_index] = donor_loglik[static_cast<R_xlen_t>(i)];
    }
  }

  if (!missing_all_records.empty()) {
    const auto finite_prob = evaluate_probability_records(
        exact_plans,
        layout,
        paramsSEXP,
        missing_all_records);
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (!missing_all_trial[trial_index]) {
        continue;
      }
      double total_finite = 0.0;
      for (std::size_t i = 0; i < missing_all_records.size(); ++i) {
        if (missing_all_records[i].trial_index !=
            static_cast<semantic::Index>(trial_index)) {
          continue;
        }
        total_finite += missing_all_records[i].weight *
                        finite_prob[static_cast<R_xlen_t>(i)];
      }
      const double missing_prob =
          std::max(0.0, 1.0 - (std::isfinite(total_finite) ? total_finite : 1.0));
      loglik[static_cast<R_xlen_t>(trial_index)] =
          missing_prob > 0.0 ? std::log(missing_prob) : min_ll;
    }
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    if (!std::isfinite(loglik[i])) {
      loglik[i] = min_ll;
    }
    total_loglik += static_cast<double>(loglik[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

inline SEXP evaluate_observed_trials_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observed_identity,
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const std::vector<unsigned char> *selected = nullptr) {
  Rcpp::DataFrame data(dataSEXP);
  const auto table = read_prepared_data_view(data);

  if (observed_identity) {
    const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
    const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
    if (trials_have_observed_components(layout, table.component, selected)) {
      if (identity_trials_are_all_finite(layout, label, rt, selected)) {
        (void)model;
        return evaluate_exact_trials_cached(
            exact_variant_index_by_component_code,
            exact_plans,
            layout,
            paramsSEXP,
            dataSEXP,
            min_ll,
            selected);
      }
      return evaluate_identity_trials_with_missing_all(
          component_plans_by_code,
          exact_variant_index_by_component_code,
          exact_plans,
          layout,
          paramsSEXP,
          dataSEXP,
          min_ll,
          selected);
    }
  }

  const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  const auto accumulator = Rcpp::as<Rcpp::CharacterVector>(data["accumulator"]);
  const NamedParamMatrix named_params(paramsSEXP);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector loglik(n_trials, min_ll);
  std::vector<ObservedTrialKind> trial_kinds(n_trials, ObservedTrialKind::Finite);
  std::vector<ObservedRecord> finite_records;
  std::vector<ObservedRecord> missing_rt_records;
  MissingAllAccumulator missing_all(n_trials);
  std::vector<std::vector<int>> row_maps;
  bool has_missing_all = false;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    if (!trial_is_selected(selected, trial_index)) {
      continue;
    }
    const auto observed_label_code =
        integer_cell_is_na(label, row)
            ? semantic::kInvalidIndex
            : static_cast<semantic::Index>(label[row]);
    const double observed_rt = rt[row];

    ObservedTrialKind kind = ObservedTrialKind::Finite;
    if (observed_label_code == semantic::kInvalidIndex &&
        Rcpp::NumericVector::is_na(observed_rt)) {
      kind = ObservedTrialKind::MissingAll;
    } else if (observed_label_code != semantic::kInvalidIndex &&
               Rcpp::NumericVector::is_na(observed_rt)) {
      kind = ObservedTrialKind::MissingRt;
    }
    trial_kinds[trial_index] = kind;

    if (kind == ObservedTrialKind::MissingAll) {
      has_missing_all = true;
    }

    const auto components = resolve_trial_components(
        component_plans_by_code,
        model,
        layout,
        table.component,
        named_params,
        trial_index);
    const auto latent_trial = integer_cell_is_na(table.component, row);
    std::unordered_map<std::string, int> row_by_accumulator;
    if (latent_trial) {
      row_by_accumulator = build_trial_row_map(accumulator, span);
    }
    for (const auto &choice : components) {
      if (choice.component_code <= 0 ||
          choice.component_code >=
              static_cast<semantic::Index>(component_plans_by_code.size())) {
        continue;
      }
      const auto &component_plan =
          component_plans_by_code[static_cast<std::size_t>(choice.component_code)];
      if (!component_plan.present) {
        continue;
      }
      const auto variant_index = resolve_variant_index_by_component_code(
          choice.component_code,
          exact_variant_index_by_component_code);
      semantic::Index row_map_index = semantic::kInvalidIndex;
      if (latent_trial && variant_index != semantic::kInvalidIndex) {
        const auto &leaf_ids =
            exact_plans.at(static_cast<std::size_t>(variant_index))
                .lowered.leaf_ids;
        row_maps.push_back(build_component_row_map(row_by_accumulator, leaf_ids));
        row_map_index =
            static_cast<semantic::Index>(row_maps.size() - 1U);
      }
      append_trial_component_records(
          static_cast<semantic::Index>(trial_index),
          kind,
          observed_label_code,
          observed_rt,
          variant_index,
          choice.weight,
          row_map_index,
          component_plan,
          &finite_records,
          &missing_rt_records,
          &missing_all);
    }
  }

  if (!finite_records.empty()) {
    const auto donor_loglik = evaluate_loglik_records(
        exact_plans,
        layout,
        paramsSEXP,
        finite_records,
        min_ll,
        &row_maps);
    std::size_t i = 0;
    while (i < finite_records.size()) {
      const auto trial_index = finite_records[i].trial_index;
      std::vector<double> values;
      while (i < finite_records.size() &&
             finite_records[i].trial_index == trial_index) {
        if (finite_records[i].weight > 0.0) {
          values.push_back(std::log(finite_records[i].weight) +
                           donor_loglik[static_cast<R_xlen_t>(i)]);
        }
        ++i;
      }
      const auto value = logsumexp_records(values);
      loglik[static_cast<R_xlen_t>(trial_index)] =
          std::isfinite(value) ? value : min_ll;
    }
  }

  if (!missing_rt_records.empty()) {
    const auto donor_prob = evaluate_probability_records(
        exact_plans,
        layout,
        paramsSEXP,
        missing_rt_records,
        &row_maps);
    std::size_t i = 0;
    while (i < missing_rt_records.size()) {
      const auto trial_index = missing_rt_records[i].trial_index;
      double total = 0.0;
      while (i < missing_rt_records.size() &&
             missing_rt_records[i].trial_index == trial_index) {
        total += missing_rt_records[i].weight *
                 donor_prob[static_cast<R_xlen_t>(i)];
        ++i;
      }
      loglik[static_cast<R_xlen_t>(trial_index)] =
          std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
    }
  }

  if (has_missing_all) {
    Rcpp::NumericVector exact_finite_prob;
    if (!missing_all.exact_records.empty()) {
      exact_finite_prob = evaluate_probability_records(
          exact_plans,
          layout,
          paramsSEXP,
          missing_all.exact_records,
          &row_maps);
    }
    std::size_t exact_i = 0;
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (trial_kinds[trial_index] != ObservedTrialKind::MissingAll) {
        continue;
      }
      double missing_prob = 0.0;
      double exact_total_finite = 0.0;
      while (exact_i < missing_all.exact_records.size() &&
             missing_all.exact_records[exact_i].trial_index ==
                 static_cast<semantic::Index>(trial_index)) {
        exact_total_finite += missing_all.exact_records[exact_i].weight *
                              exact_finite_prob[static_cast<R_xlen_t>(exact_i)];
        ++exact_i;
      }
      const double exact_missing = std::max(
          0.0,
          missing_all.exact_component_weight_sum[trial_index] -
              (std::isfinite(exact_total_finite) ? exact_total_finite
                                                 : missing_all.exact_component_weight_sum[trial_index]));
      missing_prob += exact_missing;
      loglik[static_cast<R_xlen_t>(trial_index)] =
          missing_prob > 0.0 ? std::log(missing_prob) : min_ll;
    }
  }

  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    if (!std::isfinite(loglik[i])) {
      loglik[i] = min_ll;
    }
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    total_loglik += static_cast<double>(loglik[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

} // namespace detail
} // namespace accumulatr::eval
