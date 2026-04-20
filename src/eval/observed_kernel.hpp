#pragma once

#include <Rcpp.h>

#include "direct_kernel.hpp"
#include "exact_kernel.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

inline std::string read_string_cell(const Rcpp::CharacterVector &column,
                                    const R_xlen_t row) {
  if (string_cell_is_na(column, row)) {
    return {};
  }
  return Rcpp::as<std::string>(column[row]);
}

inline std::vector<ObservedTrialPlan> build_observed_trial_plan(
    const Rcpp::DataFrame &data,
    const std::unordered_map<std::string, ComponentObservationPlan> &component_plans) {
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return {};
  }
  if (!data.containsElementNamed("R") || !data.containsElementNamed("rt")) {
    throw std::runtime_error("observed likelihood evaluator requires columns R and rt");
  }

  const int max_rank = detect_rank_count(data, "");

  const auto table = read_trial_table_view(data);
  const auto label = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);

  std::vector<Rcpp::CharacterVector> rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  std::vector<ObservedTrialPlan> out;
  int last_trial = NA_INTEGER;
  for (R_xlen_t i = 0; i < n_rows; ++i) {
    if (i > 0 && table.trial[i] == last_trial) {
      continue;
    }
    const auto start = static_cast<semantic::Index>(i);
    semantic::Index end = start;
    while (static_cast<R_xlen_t>(end + 1) < n_rows &&
           table.trial[end + 1] == table.trial[i]) {
      ++end;
    }
    last_trial = table.trial[i];

    std::string component_id =
        component_id_at_row(table, i, component_plans.size() == 1U
                                          ? component_plans.begin()->first
                                          : "__default__");
    const auto component_it = component_plans.find(component_id);
    if (component_it == component_plans.end()) {
      throw std::runtime_error(
          "observed likelihood evaluator found no variant for component '" +
          component_id + "'");
    }

    const std::string observed_label = read_string_cell(label, i);
    const double observed_rt = rt[i];
    bool has_ranked = false;
    for (int rank = 2; rank <= max_rank; ++rank) {
      const bool label_missing =
          STRING_ELT(rank_labels[static_cast<std::size_t>(rank)], i) == NA_STRING;
      const bool rt_missing =
          Rcpp::NumericVector::is_na(rank_times[static_cast<std::size_t>(rank)][i]);
      if (!label_missing || !rt_missing) {
        has_ranked = true;
        break;
      }
    }

    ObservedTrialPlan plan;
    plan.trial_index = static_cast<semantic::Index>(out.size());
    plan.start_row = start;
    plan.end_row = end;
    plan.component_id = component_id;
    plan.observed_label = observed_label;
    plan.observed_rt = observed_rt;

    if (has_ranked) {
      if (component_it->second.observed_backend != compile::BackendKind::Exact) {
        throw std::runtime_error("ranked observations require the exact backend");
      }
      plan.kind = ObservedTrialKind::Ranked;
    } else if (observed_label.empty() && Rcpp::NumericVector::is_na(observed_rt)) {
      plan.kind = ObservedTrialKind::MissingAll;
    } else if (!observed_label.empty() && Rcpp::NumericVector::is_na(observed_rt)) {
      plan.kind = ObservedTrialKind::MissingRt;
    } else if (observed_label.empty()) {
      plan.kind = ObservedTrialKind::UnsupportedMissingLabel;
    } else {
      plan.kind = ObservedTrialKind::Finite;
    }

    out.push_back(std::move(plan));
  }
  return out;
}

inline Rcpp::NumericMatrix subset_param_rows(const Rcpp::NumericMatrix &params,
                                             const std::vector<int> &rows) {
  Rcpp::NumericMatrix out(static_cast<R_xlen_t>(rows.size()), params.ncol());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    for (int j = 0; j < params.ncol(); ++j) {
      out(static_cast<R_xlen_t>(i), j) = params(rows[i], j);
    }
  }
  out.attr("dimnames") = params.attr("dimnames");
  return out;
}

inline Rcpp::DataFrame copy_trial_blocks(
    const Rcpp::DataFrame &data,
    const std::vector<ObservedTrialPlan> &trials,
    const std::vector<semantic::Index> &trial_indices,
    const std::vector<std::string> *replacement_labels,
    Rcpp::NumericMatrix *out_params,
    const Rcpp::NumericMatrix &params) {
  int max_rank = 1;
  for (int rank = 2;; ++rank) {
    const std::string r_col = "R" + std::to_string(rank);
    const std::string rt_col = "rt" + std::to_string(rank);
    const bool has_r = data.containsElementNamed(r_col.c_str());
    const bool has_rt = data.containsElementNamed(rt_col.c_str());
    if (!has_r || !has_rt) {
      break;
    }
    max_rank = rank;
  }

  const auto source_trial =
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(data.nrows());
  const auto source_r = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  const auto source_rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  const bool has_component = data.containsElementNamed("component");
  const auto source_component =
      has_component ? Rcpp::as<Rcpp::CharacterVector>(data["component"])
                    : Rcpp::CharacterVector(data.nrows(), NA_STRING);
  std::vector<Rcpp::CharacterVector> rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  std::vector<int> row_indices;
  for (const auto trial_index : trial_indices) {
    const auto &trial = trials[static_cast<std::size_t>(trial_index)];
    for (semantic::Index row = trial.start_row; row <= trial.end_row; ++row) {
      row_indices.push_back(static_cast<int>(row));
    }
  }
  *out_params = subset_param_rows(params, row_indices);

  Rcpp::IntegerVector new_trial(row_indices.size());
  Rcpp::CharacterVector new_r(row_indices.size(), NA_STRING);
  Rcpp::NumericVector new_rt(row_indices.size(), NA_REAL);
  Rcpp::CharacterVector new_component(row_indices.size(), NA_STRING);
  std::vector<Rcpp::CharacterVector> new_rank_labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> new_rank_times(static_cast<std::size_t>(max_rank + 1));
  for (int rank = 2; rank <= max_rank; ++rank) {
    new_rank_labels[static_cast<std::size_t>(rank)] =
        Rcpp::CharacterVector(row_indices.size(), NA_STRING);
    new_rank_times[static_cast<std::size_t>(rank)] =
        Rcpp::NumericVector(row_indices.size(), NA_REAL);
  }

  std::size_t out_row = 0;
  for (std::size_t i = 0; i < trial_indices.size(); ++i) {
    const auto &trial = trials[static_cast<std::size_t>(trial_indices[i])];
    for (semantic::Index row = trial.start_row; row <= trial.end_row; ++row, ++out_row) {
      new_trial[static_cast<R_xlen_t>(out_row)] = static_cast<int>(i + 1U);
      new_r[static_cast<R_xlen_t>(out_row)] = source_r[static_cast<R_xlen_t>(row)];
      new_rt[static_cast<R_xlen_t>(out_row)] = source_rt[static_cast<R_xlen_t>(row)];
      if (has_component) {
        new_component[static_cast<R_xlen_t>(out_row)] =
            source_component[static_cast<R_xlen_t>(row)];
      }
      for (int rank = 2; rank <= max_rank; ++rank) {
        new_rank_labels[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(out_row)] =
            rank_labels[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(row)];
        new_rank_times[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(out_row)] =
            rank_times[static_cast<std::size_t>(rank)][static_cast<R_xlen_t>(row)];
      }
    }
    if (replacement_labels != nullptr) {
      const auto replacement = (*replacement_labels)[i];
      out_row -= static_cast<std::size_t>(trial.end_row - trial.start_row + 1);
      new_r[static_cast<R_xlen_t>(out_row)] = replacement;
      out_row += static_cast<std::size_t>(trial.end_row - trial.start_row + 1);
    }
  }

  Rcpp::List out;
  Rcpp::CharacterVector names;
  out.push_back(new_trial, "trial");
  out.push_back(new_r, "R");
  out.push_back(new_rt, "rt");
  if (has_component) {
    out.push_back(new_component, "component");
  }
  for (int rank = 2; rank <= max_rank; ++rank) {
    out.push_back(new_rank_labels[static_cast<std::size_t>(rank)],
                  ("R" + std::to_string(rank)).c_str());
    out.push_back(new_rank_times[static_cast<std::size_t>(rank)],
                  ("rt" + std::to_string(rank)).c_str());
  }
  out.attr("class") = "data.frame";
  out.attr("row.names") =
      Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(row_indices.size()));
  return Rcpp::DataFrame(out);
}

inline Rcpp::NumericVector evaluate_loglik_records(
    const compile::CompiledModel &compiled,
    const semantic::SemanticModel &model,
    const std::vector<ObservedTrialPlan> &trials,
    const Rcpp::DataFrame &data,
    const Rcpp::NumericMatrix &params,
    const std::vector<ObservedRecord> &records,
    const double min_ll) {
  Rcpp::NumericVector out(records.size(), R_NegInf);
  if (records.empty()) {
    return out;
  }

  for (const auto backend : {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
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

    std::vector<semantic::Index> trial_indices;
    std::vector<std::string> labels;
    trial_indices.reserve(idxs.size());
    labels.reserve(idxs.size());
    for (const auto idx : idxs) {
      trial_indices.push_back(records[static_cast<std::size_t>(idx)].trial_index);
      labels.push_back(records[static_cast<std::size_t>(idx)].semantic_label);
    }

    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, trial_indices, &labels, &subset_params, params);
    SEXP eval_sexp = R_NilValue;
    if (backend == compile::BackendKind::Direct) {
      eval_sexp = evaluate_direct_trials(
          compiled, model, subset_params, subset_data, min_ll);
    } else {
      eval_sexp = evaluate_exact_trials(
          compiled, subset_params, subset_data, min_ll);
    }
    Rcpp::List eval(eval_sexp);
    const auto loglik = Rcpp::as<Rcpp::NumericVector>(eval["loglik"]);
    for (std::size_t i = 0; i < idxs.size(); ++i) {
      out[static_cast<R_xlen_t>(idxs[i])] = loglik[static_cast<R_xlen_t>(i)];
    }
  }

  return out;
}

inline Rcpp::NumericVector evaluate_probability_records(
    const compile::CompiledModel &compiled,
    const semantic::SemanticModel &model,
    const std::vector<ObservedTrialPlan> &trials,
    const Rcpp::DataFrame &data,
    const Rcpp::NumericMatrix &params,
    const std::vector<ObservedRecord> &records) {
  Rcpp::NumericVector out(records.size(), 0.0);
  if (records.empty()) {
    return out;
  }

  for (const auto backend : {compile::BackendKind::Direct, compile::BackendKind::Exact}) {
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

    std::vector<semantic::Index> trial_indices;
    std::vector<std::string> labels;
    trial_indices.reserve(idxs.size());
    labels.reserve(idxs.size());
    for (const auto idx : idxs) {
      trial_indices.push_back(records[static_cast<std::size_t>(idx)].trial_index);
      labels.push_back(records[static_cast<std::size_t>(idx)].semantic_label);
    }

    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, trial_indices, &labels, &subset_params, params);
    SEXP eval_sexp = R_NilValue;
    if (backend == compile::BackendKind::Direct) {
      eval_sexp = evaluate_direct_outcome_probabilities(
          compiled, model, subset_params, subset_data);
    } else {
      eval_sexp = evaluate_exact_outcome_probabilities(
          compiled, subset_params, subset_data);
    }
    Rcpp::List eval(eval_sexp);
    const auto prob = Rcpp::as<Rcpp::NumericVector>(eval["probability"]);
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

inline SEXP evaluate_observed_trials(SEXP observedPlanSEXP,
                                     const semantic::SemanticModel &model,
                                     const compile::CompiledModel &compiled,
                                     SEXP paramsSEXP,
                                     SEXP dataSEXP,
                                     const double min_ll) {
  Rcpp::DataFrame data(dataSEXP);
  Rcpp::NumericMatrix params(Rcpp::as<Rcpp::NumericMatrix>(paramsSEXP));
  const auto component_plans = observed_plan_from_r(observedPlanSEXP);
  const auto trials = build_observed_trial_plan(data, component_plans);

  Rcpp::NumericVector loglik(trials.size(), min_ll);

  std::vector<semantic::Index> ranked_trial_indices;
  ranked_trial_indices.reserve(trials.size());
  std::vector<ObservedRecord> finite_records;
  std::vector<ObservedRecord> missing_rt_records;
  std::vector<ObservedRecord> missing_all_records;

  for (const auto &trial : trials) {
    const auto component_it = component_plans.find(trial.component_id);
    if (component_it == component_plans.end()) {
      throw std::runtime_error(
          "observed likelihood evaluator found no observation plan for component '" +
          trial.component_id + "'");
    }
    const auto &component_plan = component_it->second;

    switch (trial.kind) {
    case ObservedTrialKind::Ranked:
      ranked_trial_indices.push_back(trial.trial_index);
      break;
    case ObservedTrialKind::Finite: {
      const auto keep_it = component_plan.keep.find(trial.observed_label);
      if (keep_it != component_plan.keep.end()) {
        for (const auto &branch : keep_it->second) {
          finite_records.push_back(ObservedRecord{
              trial.trial_index,
              component_plan.semantic_backend,
              branch.semantic_label,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingRt: {
      const auto missing_it = component_plan.missing_rt.find(trial.observed_label);
      if (missing_it != component_plan.missing_rt.end()) {
        for (const auto &branch : missing_it->second) {
          missing_rt_records.push_back(ObservedRecord{
              trial.trial_index,
              component_plan.semantic_backend,
              branch.semantic_label,
              branch.weight});
        }
      }
      break;
    }
    case ObservedTrialKind::MissingAll:
      for (const auto &branch : component_plan.finite_observed_branches) {
        missing_all_records.push_back(ObservedRecord{
            trial.trial_index,
            component_plan.semantic_backend,
            branch.semantic_label,
            branch.weight});
      }
      break;
    case ObservedTrialKind::UnsupportedMissingLabel:
      break;
    }
  }

  if (!ranked_trial_indices.empty()) {
    Rcpp::NumericMatrix subset_params;
    const auto subset_data =
        copy_trial_blocks(data, trials, ranked_trial_indices, nullptr, &subset_params, params);
    Rcpp::List eval(
        evaluate_exact_trials(compiled, subset_params, subset_data, min_ll));
    const auto ranked_loglik = Rcpp::as<Rcpp::NumericVector>(eval["loglik"]);
    for (std::size_t i = 0; i < ranked_trial_indices.size(); ++i) {
      loglik[static_cast<R_xlen_t>(ranked_trial_indices[i])] =
          ranked_loglik[static_cast<R_xlen_t>(i)];
    }
  }

  if (!finite_records.empty()) {
    const auto donor_loglik = evaluate_loglik_records(
        compiled, model, trials, data, params, finite_records, min_ll);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::Finite) {
        continue;
      }
      std::vector<double> values;
      for (std::size_t i = 0; i < finite_records.size(); ++i) {
        if (finite_records[i].trial_index != trial.trial_index) {
          continue;
        }
        if (!(finite_records[i].weight > 0.0)) {
          continue;
        }
        values.push_back(std::log(finite_records[i].weight) + donor_loglik[static_cast<R_xlen_t>(i)]);
      }
      const auto value = logsumexp_records(values);
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          std::isfinite(value) ? value : min_ll;
    }
  }

  if (!missing_rt_records.empty()) {
    const auto donor_prob = evaluate_probability_records(
        compiled, model, trials, data, params, missing_rt_records);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::MissingRt) {
        continue;
      }
      double total = 0.0;
      for (std::size_t i = 0; i < missing_rt_records.size(); ++i) {
        if (missing_rt_records[i].trial_index != trial.trial_index) {
          continue;
        }
        total += missing_rt_records[i].weight * donor_prob[static_cast<R_xlen_t>(i)];
      }
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
    }
  }

  if (!missing_all_records.empty()) {
    const auto finite_prob = evaluate_probability_records(
        compiled, model, trials, data, params, missing_all_records);
    for (const auto &trial : trials) {
      if (trial.kind != ObservedTrialKind::MissingAll) {
        continue;
      }
      double total_finite = 0.0;
      for (std::size_t i = 0; i < missing_all_records.size(); ++i) {
        if (missing_all_records[i].trial_index != trial.trial_index) {
          continue;
        }
        total_finite += missing_all_records[i].weight * finite_prob[static_cast<R_xlen_t>(i)];
      }
      const double missing_prob =
          std::max(0.0, 1.0 - (std::isfinite(total_finite) ? total_finite : 1.0));
      loglik[static_cast<R_xlen_t>(trial.trial_index)] =
          missing_prob > 0.0 ? std::log(missing_prob) : min_ll;
    }
  } else {
    for (const auto &trial : trials) {
      if (trial.kind == ObservedTrialKind::MissingAll) {
        loglik[static_cast<R_xlen_t>(trial.trial_index)] = 0.0;
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
      Rcpp::Named("n_trials") = static_cast<int>(trials.size()));
}

} // namespace detail
} // namespace accumulatr::eval
