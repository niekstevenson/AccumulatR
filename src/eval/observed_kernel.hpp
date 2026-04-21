#pragma once

#include <Rcpp.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "direct_kernel.hpp"
#include "exact_kernel.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

inline Rcpp::NumericVector evaluate_loglik_records(
    const std::vector<VariantPlan> &direct_plans,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<ObservedRecord> &records,
    const double min_ll) {
  Rcpp::NumericVector out(records.size(), R_NegInf);
  if (records.empty()) {
    return out;
  }

  for (const auto backend :
       {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
    std::vector<semantic::Index> idxs;
    idxs.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
      if (records[i].backend == backend) {
        idxs.push_back(static_cast<semantic::Index>(i));
      }
    }
    if (idxs.empty()) {
      continue;
    }

    std::vector<DirectLoglikQuery> queries;
    queries.reserve(idxs.size());
    for (const auto idx : idxs) {
      const auto &record = records[static_cast<std::size_t>(idx)];
      queries.push_back(DirectLoglikQuery{
          record.trial_index,
          record.variant_index,
          record.semantic_code,
          record.observed_rt});
    }

    Rcpp::NumericVector loglik;
    if (backend == compile::BackendKind::Direct) {
      loglik = evaluate_direct_loglik_queries_cached(
          direct_plans,
          layout,
          paramsSEXP,
          queries,
          min_ll);
    } else {
      loglik = evaluate_exact_loglik_queries_cached(
          exact_plans,
          layout,
          paramsSEXP,
          queries,
          min_ll);
    }
    for (std::size_t i = 0; i < idxs.size(); ++i) {
      out[static_cast<R_xlen_t>(idxs[i])] = loglik[static_cast<R_xlen_t>(i)];
    }
  }

  return out;
}

inline Rcpp::NumericVector evaluate_probability_records(
    const std::vector<VariantPlan> &direct_plans,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<ObservedRecord> &records) {
  Rcpp::NumericVector out(records.size(), 0.0);
  if (records.empty()) {
    return out;
  }

  for (const auto backend :
       {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
    std::vector<semantic::Index> idxs;
    idxs.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
      if (records[i].backend == backend) {
        idxs.push_back(static_cast<semantic::Index>(i));
      }
    }
    if (idxs.empty()) {
      continue;
    }

    std::vector<DirectProbabilityQuery> queries;
    queries.reserve(idxs.size());
    for (const auto idx : idxs) {
      const auto &record = records[static_cast<std::size_t>(idx)];
      queries.push_back(DirectProbabilityQuery{
          record.trial_index,
          record.variant_index,
          record.semantic_code});
    }

    Rcpp::NumericVector prob;
    if (backend == compile::BackendKind::Direct) {
      prob = evaluate_direct_probability_queries_cached(
          direct_plans,
          layout,
          paramsSEXP,
          queries);
    } else {
      prob = evaluate_exact_probability_queries_cached(
          exact_plans,
          layout,
          paramsSEXP,
          queries);
    }
    for (std::size_t i = 0; i < idxs.size(); ++i) {
      out[static_cast<R_xlen_t>(idxs[i])] = prob[static_cast<R_xlen_t>(i)];
    }
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
    const compile::BackendKind backend,
    const std::vector<semantic::Index> &direct_variant_index_by_component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  const auto idx = static_cast<std::size_t>(component_code);
  return backend == compile::BackendKind::Direct
             ? direct_variant_index_by_component_code[idx]
             : exact_variant_index_by_component_code[idx];
}

inline bool identity_trials_are_all_finite(
    const PreparedTrialLayout &layout,
    const Rcpp::IntegerVector &label,
    const Rcpp::NumericVector &rt) {
  for (const auto &span : layout.spans) {
    const auto row = static_cast<R_xlen_t>(span.start_row);
    if (integer_cell_is_na(label, row) || Rcpp::NumericVector::is_na(rt[row])) {
      return false;
    }
  }
  return true;
}

inline SEXP evaluate_identity_trials_with_missing_all(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &direct_variant_index_by_component_code,
    const std::vector<VariantPlan> &direct_plans,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll) {
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
    const auto component_code = static_cast<semantic::Index>(table.component[row]);
    const auto &component_plan =
        component_plans_by_code[static_cast<std::size_t>(component_code)];
    const auto variant_index = resolve_variant_index_by_component_code(
        component_code,
        component_plan.semantic_backend,
        direct_variant_index_by_component_code,
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
        direct_plans,
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
        direct_plans,
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
    const compile::BackendKind identity_backend,
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &direct_variant_index_by_component_code,
    const std::vector<VariantPlan> &direct_plans,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll) {
  Rcpp::DataFrame data(dataSEXP);

  if (observed_identity) {
    const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
    const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
    if (identity_trials_are_all_finite(layout, label, rt)) {
      if (identity_backend == compile::BackendKind::Direct) {
        return evaluate_direct_trials_cached(
            model,
            direct_variant_index_by_component_code,
            direct_plans,
            layout,
            paramsSEXP,
            dataSEXP,
            min_ll);
      }
      return evaluate_exact_trials_cached(
          exact_variant_index_by_component_code,
          exact_plans,
          layout,
          paramsSEXP,
          dataSEXP,
          min_ll);
    }
    return evaluate_identity_trials_with_missing_all(
        component_plans_by_code,
        direct_variant_index_by_component_code,
        direct_plans,
        exact_variant_index_by_component_code,
        exact_plans,
        layout,
        paramsSEXP,
        dataSEXP,
        min_ll);
  }

  const auto table = read_prepared_data_view(data);
  const auto label = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector loglik(n_trials, min_ll);
  std::vector<ObservedTrialKind> trial_kinds(n_trials, ObservedTrialKind::Finite);
  std::vector<ObservedRecord> finite_records;
  std::vector<ObservedRecord> missing_rt_records;
  std::vector<ObservedRecord> missing_all_records;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    const int component_code_int = table.component[row];
    const auto component_code = static_cast<semantic::Index>(component_code_int);
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

    const auto &component_plan =
        component_plans_by_code[static_cast<std::size_t>(component_code)];
    const auto variant_index = resolve_variant_index_by_component_code(
        component_code,
        component_plan.semantic_backend,
        direct_variant_index_by_component_code,
        exact_variant_index_by_component_code);

    switch (kind) {
    case ObservedTrialKind::Finite: {
      const auto observed_code = static_cast<std::size_t>(observed_label_code);
      if (observed_code < component_plan.keep_by_code.size()) {
        for (const auto &branch : component_plan.keep_by_code[observed_code]) {
          finite_records.push_back(ObservedRecord{
              static_cast<semantic::Index>(trial_index),
              variant_index,
              component_plan.semantic_backend,
              branch.semantic_code,
              observed_rt,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingRt: {
      const auto observed_code = static_cast<std::size_t>(observed_label_code);
      if (observed_code < component_plan.missing_rt_by_code.size()) {
        for (const auto &branch :
             component_plan.missing_rt_by_code[observed_code]) {
          missing_rt_records.push_back(ObservedRecord{
              static_cast<semantic::Index>(trial_index),
              variant_index,
              component_plan.semantic_backend,
              branch.semantic_code,
              observed_rt,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingAll:
      for (const auto &branch : component_plan.finite_observed_branches) {
        missing_all_records.push_back(ObservedRecord{
            static_cast<semantic::Index>(trial_index),
            variant_index,
            component_plan.semantic_backend,
            branch.semantic_code,
            observed_rt,
            branch.weight});
      }
      break;
    }
  }

  if (!finite_records.empty()) {
    const auto donor_loglik = evaluate_loglik_records(
        direct_plans,
        exact_plans,
        layout,
        paramsSEXP,
        finite_records,
        min_ll);
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (trial_kinds[trial_index] != ObservedTrialKind::Finite) {
        continue;
      }
      std::vector<double> values;
      for (std::size_t i = 0; i < finite_records.size(); ++i) {
        if (finite_records[i].trial_index !=
                static_cast<semantic::Index>(trial_index) ||
            !(finite_records[i].weight > 0.0)) {
          continue;
        }
        values.push_back(std::log(finite_records[i].weight) +
                         donor_loglik[static_cast<R_xlen_t>(i)]);
      }
      const auto value = logsumexp_records(values);
      loglik[static_cast<R_xlen_t>(trial_index)] =
          std::isfinite(value) ? value : min_ll;
    }
  }

  if (!missing_rt_records.empty()) {
    const auto donor_prob = evaluate_probability_records(
        direct_plans,
        exact_plans,
        layout,
        paramsSEXP,
        missing_rt_records);
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (trial_kinds[trial_index] != ObservedTrialKind::MissingRt) {
        continue;
      }
      double total = 0.0;
      for (std::size_t i = 0; i < missing_rt_records.size(); ++i) {
        if (missing_rt_records[i].trial_index !=
            static_cast<semantic::Index>(trial_index)) {
          continue;
        }
        total += missing_rt_records[i].weight *
                 donor_prob[static_cast<R_xlen_t>(i)];
      }
      loglik[static_cast<R_xlen_t>(trial_index)] =
          std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
    }
  }

  if (!missing_all_records.empty()) {
    const auto finite_prob = evaluate_probability_records(
        direct_plans,
        exact_plans,
        layout,
        paramsSEXP,
        missing_all_records);
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (trial_kinds[trial_index] != ObservedTrialKind::MissingAll) {
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
  } else {
    for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
      if (trial_kinds[trial_index] == ObservedTrialKind::MissingAll) {
        loglik[static_cast<R_xlen_t>(trial_index)] = 0.0;
      }
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
