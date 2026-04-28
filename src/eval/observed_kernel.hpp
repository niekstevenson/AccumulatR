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

struct TrustedParamMatrix {
  const double *base{nullptr};
  int nrow{0};
  int component_weight_start{0};

  TrustedParamMatrix(SEXP paramsSEXP, const int component_weight_param_count)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        component_weight_start(
            Rf_ncols(paramsSEXP) - component_weight_param_count) {}

  inline double component_weight(const int row,
                                 const int weight_param_index) const {
    return base[
        static_cast<R_xlen_t>(component_weight_start + weight_param_index) *
            nrow +
        row];
  }

};

inline bool trials_have_observed_components(
    const PreparedTrialLayout &layout,
    const int *component,
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
    const TrustedParamMatrix &params,
    const int row) {
  std::vector<double> weights(component_codes.size(), 0.0);
  if (component_codes.empty()) {
    return weights;
  }

  if (model.component_mode != "sample") {
    for (std::size_t i = 0; i < component_codes.size(); ++i) {
      const auto component_index =
          static_cast<std::size_t>(component_codes[i] - 1);
      weights[i] = model.components[component_index].weight;
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
      if (component.weight_param_index >= 0) {
        weight = params.component_weight(row, component.weight_param_index);
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
      weights[i] = ref_weight;
      break;
    }
  }

  double total = 0.0;
  for (const auto weight : weights) {
    total += weight;
  }
  const double inv_total = 1.0 / total;
  for (auto &weight : weights) {
    weight *= inv_total;
  }
  return weights;
}

struct TrialComponentChoice {
  semantic::Index component_code{semantic::kInvalidIndex};
  double weight{0.0};
};

struct ObservationPlanRequest {
  const ObservationProbabilityPlan *plan{nullptr};
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  double observed_rt{NA_REAL};
  double component_weight{1.0};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

struct ObservationLeafTarget {
  std::size_t value_index{0};
  double weight{1.0};
};

inline std::unordered_map<std::string, int> build_trial_row_map(
    SEXP accumulator,
    const PreparedTrialSpan &span) {
  std::unordered_map<std::string, int> row_by_accumulator;
  row_by_accumulator.reserve(
      static_cast<std::size_t>(span.end_row - span.start_row + 1));
  for (auto row = static_cast<int>(span.start_row);
       row <= static_cast<int>(span.end_row);
       ++row) {
    row_by_accumulator.emplace(CHAR(STRING_ELT(accumulator, row)), row);
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
    const int *component,
    const TrustedParamMatrix &params,
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
    if (!(weights[i] > 0.0)) {
      continue;
    }
    resolved.push_back(TrialComponentChoice{component_codes[i], weights[i]});
  }
  return resolved;
}

inline void append_observation_plan_request(
    std::vector<ObservationPlanRequest> *requests,
    const semantic::Index trial_index,
    const semantic::Index variant_index,
    const double observed_rt,
    const double component_weight,
    const semantic::Index row_map_index,
    const ObservationProbabilityPlan &plan) {
  if (requests == nullptr ||
      plan.empty() ||
      variant_index == semantic::kInvalidIndex ||
      !(std::isfinite(component_weight) && component_weight > 0.0)) {
    return;
  }
  requests->push_back(ObservationPlanRequest{
      &plan,
      trial_index,
      variant_index,
      observed_rt,
      component_weight,
      row_map_index});
}

inline semantic::Index resolve_variant_index_by_component_code(
    const semantic::Index component_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code) {
  const auto idx = static_cast<std::size_t>(component_code);
  return exact_variant_index_by_component_code[idx];
}

inline bool identity_trials_are_all_finite(
    const PreparedTrialLayout &layout,
    const int *label,
    const double *rt,
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

inline std::vector<double> evaluate_observation_plan_requests(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<ObservationPlanRequest> &requests,
    const double min_ll,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  std::vector<double> roots(requests.size(), min_ll);
  if (requests.empty()) {
    return roots;
  }

  std::vector<std::size_t> request_offsets(requests.size(), 0U);
  std::size_t value_count = 0;
  for (std::size_t i = 0; i < requests.size(); ++i) {
    request_offsets[i] = value_count;
    value_count += requests[i].plan == nullptr ? 0U : requests[i].plan->ops.size();
  }

  std::vector<double> values(value_count, 0.0);
  std::vector<ObservedRecord> log_records;
  std::vector<ObservationLeafTarget> log_targets;
  std::vector<ObservedRecord> probability_records;
  std::vector<ObservationLeafTarget> probability_targets;
  std::vector<EvalNoResponseQuery> no_response_queries;
  std::vector<ObservationLeafTarget> no_response_targets;

  for (std::size_t request_index = 0; request_index < requests.size();
       ++request_index) {
    const auto &request = requests[request_index];
    const auto *plan = request.plan;
    if (plan == nullptr || plan->empty()) {
      continue;
    }
    const auto offset = request_offsets[request_index];
    for (std::size_t op_index = 0; op_index < plan->ops.size(); ++op_index) {
      const auto &op = plan->ops[op_index];
      const auto value_index = offset + op_index;
      switch (op.kind) {
      case ObservationPlanOpKind::Constant:
        values[value_index] = op.constant;
        break;
      case ObservationPlanOpKind::LogDensity:
        if (std::isfinite(op.weight) && op.weight > 0.0) {
          log_targets.push_back(ObservationLeafTarget{value_index, op.weight});
          log_records.push_back(ObservedRecord{
              request.trial_index,
              request.variant_index,
              op.backend,
              op.semantic_code,
              request.observed_rt,
              op.weight,
              request.row_map_index});
        } else {
          values[value_index] = min_ll;
        }
        break;
      case ObservationPlanOpKind::FiniteOutcomeProbability:
        if (std::isfinite(op.weight) && op.weight > 0.0) {
          probability_targets.push_back(
              ObservationLeafTarget{value_index, op.weight});
          probability_records.push_back(ObservedRecord{
              request.trial_index,
              request.variant_index,
              op.backend,
              op.semantic_code,
              NA_REAL,
              op.weight,
              request.row_map_index});
        }
        break;
      case ObservationPlanOpKind::NoResponseProbability:
        if (std::isfinite(op.weight) && op.weight > 0.0) {
          no_response_targets.push_back(
              ObservationLeafTarget{value_index, op.weight});
          no_response_queries.push_back(EvalNoResponseQuery{
              request.trial_index,
              request.variant_index,
              request.row_map_index});
        }
        break;
      case ObservationPlanOpKind::Complement:
      case ObservationPlanOpKind::WeightedSum:
      case ObservationPlanOpKind::Log:
        break;
      }
    }
  }

  if (!log_records.empty()) {
    const auto loglik = evaluate_loglik_records(
        exact_plans,
        layout,
        paramsSEXP,
        log_records,
        min_ll,
        row_maps);
    for (std::size_t i = 0; i < log_targets.size(); ++i) {
      const auto &target = log_targets[i];
      const double value = loglik[static_cast<R_xlen_t>(i)];
      values[target.value_index] =
          std::isfinite(value) && target.weight > 0.0
              ? std::log(target.weight) + value
              : min_ll;
    }
  }

  if (!probability_records.empty()) {
    const auto probabilities = evaluate_probability_records(
        exact_plans,
        layout,
        paramsSEXP,
        probability_records,
        row_maps);
    for (std::size_t i = 0; i < probability_targets.size(); ++i) {
      const auto &target = probability_targets[i];
      const double probability = probabilities[static_cast<R_xlen_t>(i)];
      values[target.value_index] =
          std::isfinite(probability) && probability > 0.0
              ? target.weight * probability
              : 0.0;
    }
  }

  if (!no_response_queries.empty()) {
    const auto probabilities = evaluate_exact_no_response_queries_cached(
        exact_plans,
        layout,
        paramsSEXP,
        no_response_queries,
        row_maps);
    for (std::size_t i = 0; i < no_response_targets.size(); ++i) {
      const double probability = probabilities[static_cast<R_xlen_t>(i)];
      if (Rcpp::NumericVector::is_na(probability) ||
          !std::isfinite(probability)) {
        throw std::runtime_error(
            "compiled no-response observation probability was unavailable");
      }
      values[no_response_targets[i].value_index] =
          no_response_targets[i].weight * probability;
    }
  }

  for (std::size_t request_index = 0; request_index < requests.size();
       ++request_index) {
    const auto &request = requests[request_index];
    const auto *plan = request.plan;
    if (plan == nullptr || plan->empty()) {
      continue;
    }
    const auto offset = request_offsets[request_index];
    for (std::size_t op_index = 0; op_index < plan->ops.size(); ++op_index) {
      const auto &op = plan->ops[op_index];
      const auto value_index = offset + op_index;
      switch (op.kind) {
      case ObservationPlanOpKind::WeightedSum: {
        if (op.value_kind == ObservationPlanValueKind::Log) {
          std::vector<double> child_values;
          child_values.reserve(static_cast<std::size_t>(op.children.size));
          for (semantic::Index i = 0; i < op.children.size; ++i) {
            const auto child = plan->child_ops[
                static_cast<std::size_t>(op.children.offset + i)];
            if (child != semantic::kInvalidIndex) {
              child_values.push_back(values[offset + static_cast<std::size_t>(child)]);
            }
          }
          const double value = logsumexp_records(child_values);
          values[value_index] = std::isfinite(value) ? value : min_ll;
        } else {
          double sum = 0.0;
          for (semantic::Index i = 0; i < op.children.size; ++i) {
            const auto child = plan->child_ops[
                static_cast<std::size_t>(op.children.offset + i)];
            if (child == semantic::kInvalidIndex) {
              continue;
            }
            const double value = values[offset + static_cast<std::size_t>(child)];
            if (std::isfinite(value) && value > 0.0) {
              sum += value;
            }
          }
          values[value_index] = sum;
        }
        break;
      }
      case ObservationPlanOpKind::Complement: {
        double sum = 0.0;
        bool unavailable = false;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = plan->child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double value = values[offset + static_cast<std::size_t>(child)];
          if (!std::isfinite(value)) {
            unavailable = true;
            break;
          }
          sum += value;
        }
        values[value_index] = unavailable ? 0.0 : std::max(0.0, 1.0 - sum);
        break;
      }
      case ObservationPlanOpKind::Log: {
        double probability = 0.0;
        if (op.children.size > 0) {
          const auto child = plan->child_ops[
              static_cast<std::size_t>(op.children.offset)];
          if (child != semantic::kInvalidIndex) {
            probability = values[offset + static_cast<std::size_t>(child)];
          }
        }
        values[value_index] =
            std::isfinite(probability) && probability > 0.0
                ? std::log(probability)
                : min_ll;
        break;
      }
      case ObservationPlanOpKind::Constant:
      case ObservationPlanOpKind::LogDensity:
      case ObservationPlanOpKind::FiniteOutcomeProbability:
      case ObservationPlanOpKind::NoResponseProbability:
        break;
      }
    }
    if (plan->root != semantic::kInvalidIndex) {
      roots[request_index] =
          values[offset + static_cast<std::size_t>(plan->root)];
    }
  }

  return roots;
}

inline double evaluate_observation_plan_direct(
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const ObservationProbabilityPlan &obs_plan,
    const semantic::Index trial_index,
    const semantic::Index variant_index,
    const double observed_rt,
    const double min_ll,
    const int *row_map,
    const int row_offset,
    ExactStepWorkspacePool *workspace_pool,
    std::vector<double> *values) {
  if (values == nullptr || workspace_pool == nullptr || obs_plan.empty()) {
    return min_ll;
  }
  values->assign(obs_plan.ops.size(), 0.0);
  const auto &exact_plan = exact_plans[static_cast<std::size_t>(variant_index)];
  ParamView params(paramsSEXP, row_map, row_offset);
  const int first_param_row =
      row_map == nullptr
          ? static_cast<int>(
                layout.spans[static_cast<std::size_t>(trial_index)].start_row)
          : 0;
  auto &workspace = workspace_pool->get(exact_plans, variant_index);

  for (std::size_t op_index = 0; op_index < obs_plan.ops.size(); ++op_index) {
    const auto &op = obs_plan.ops[op_index];
    double value = 0.0;
    switch (op.kind) {
    case ObservationPlanOpKind::Constant:
      value = op.constant;
      break;
    case ObservationPlanOpKind::LogDensity: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      const double density = exact_unranked_target_density(
          exact_plan,
          params,
          first_param_row,
          target_idx,
          observed_rt,
          &workspace);
      value =
          std::isfinite(density) && density > 0.0 && op.weight > 0.0
              ? std::log(op.weight) + std::log(density)
              : min_ll;
      break;
    }
    case ObservationPlanOpKind::FiniteOutcomeProbability: {
      const auto target_idx =
          exact_plan.outcome_index_by_code[static_cast<std::size_t>(
              op.semantic_code)];
      const double probability = integrate_to_infinity(
          [&](const double rt) {
            return exact_unranked_target_density(
                exact_plan,
                params,
                first_param_row,
                target_idx,
                rt,
                &workspace);
          });
      value =
          std::isfinite(probability) && probability > 0.0 && op.weight > 0.0
              ? op.weight * probability
              : 0.0;
      break;
    }
    case ObservationPlanOpKind::NoResponseProbability: {
      double probability = 0.0;
      if (!exact_no_response_probability(
              exact_plan,
              params,
              first_param_row,
              &probability,
              &workspace)) {
        throw std::runtime_error(
            "compiled no-response observation probability was unavailable");
      }
      value =
          std::isfinite(probability) && probability > 0.0 && op.weight > 0.0
              ? op.weight * probability
              : 0.0;
      break;
    }
    case ObservationPlanOpKind::WeightedSum:
      if (op.value_kind == ObservationPlanValueKind::Log) {
        double anchor = R_NegInf;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value) && child_value > anchor) {
            anchor = child_value;
          }
        }
        if (!std::isfinite(anchor)) {
          value = min_ll;
          break;
        }
        double sum = 0.0;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value)) {
            sum += std::exp(child_value - anchor);
          }
        }
        value = sum > 0.0 ? anchor + std::log(sum) : min_ll;
      } else {
        double sum = 0.0;
        for (semantic::Index i = 0; i < op.children.size; ++i) {
          const auto child = obs_plan.child_ops[
              static_cast<std::size_t>(op.children.offset + i)];
          if (child == semantic::kInvalidIndex) {
            continue;
          }
          const double child_value =
              (*values)[static_cast<std::size_t>(child)];
          if (std::isfinite(child_value) && child_value > 0.0) {
            sum += child_value;
          }
        }
        value = sum;
      }
      break;
    case ObservationPlanOpKind::Complement: {
      double sum = 0.0;
      for (semantic::Index i = 0; i < op.children.size; ++i) {
        const auto child = obs_plan.child_ops[
            static_cast<std::size_t>(op.children.offset + i)];
        if (child == semantic::kInvalidIndex) {
          continue;
        }
        const double child_value =
            (*values)[static_cast<std::size_t>(child)];
        if (std::isfinite(child_value)) {
          sum += child_value;
        }
      }
      value = std::max(0.0, 1.0 - sum);
      break;
    }
    case ObservationPlanOpKind::Log: {
      double probability = 0.0;
      if (op.children.size > 0) {
        const auto child =
            obs_plan.child_ops[static_cast<std::size_t>(op.children.offset)];
        if (child != semantic::kInvalidIndex) {
          probability = (*values)[static_cast<std::size_t>(child)];
        }
      }
      value =
          std::isfinite(probability) && probability > 0.0
              ? std::log(probability)
              : min_ll;
      break;
    }
    }
    (*values)[op_index] = value;
  }

  return obs_plan.root == semantic::kInvalidIndex
             ? min_ll
             : (*values)[static_cast<std::size_t>(obs_plan.root)];
}

inline Rcpp::List evaluate_outcome_queries_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector probability(n_trials, 0.0);
  std::vector<ObservationPlanRequest> requests;
  requests.reserve(n_trials);

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
        static_cast<semantic::Index>(label[row]);
    if (observed_code <= 0 ||
        static_cast<std::size_t>(observed_code) >=
            component_plan.missing_rt_probability_plans_by_code.size()) {
      continue;
    }
    append_observation_plan_request(
        &requests,
        static_cast<semantic::Index>(trial_index),
        variant_index,
        NA_REAL,
        1.0,
        semantic::kInvalidIndex,
        component_plan.missing_rt_probability_plans_by_code[
            static_cast<std::size_t>(observed_code)]);
  }

  const auto probabilities = evaluate_observation_plan_requests(
      exact_plans,
      layout,
      paramsSEXP,
      requests,
      std::log(1e-10));
  for (std::size_t i = 0; i < requests.size(); ++i) {
    const auto trial_index = static_cast<R_xlen_t>(requests[i].trial_index);
    const double value = probabilities[i] * requests[i].component_weight;
    if (std::isfinite(value) && value > 0.0) {
      probability[trial_index] += value;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = probability,
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
  const auto table = read_prepared_data_view(dataSEXP, layout);

  if (observed_identity) {
    const int *label =
        INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
    const double *rt =
        REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
    if (trials_have_observed_components(layout, table.component, selected)) {
      if (identity_trials_are_all_finite(layout, label, rt, selected)) {
        return evaluate_exact_trials_cached(
            exact_variant_index_by_component_code,
            exact_plans,
            layout,
            paramsSEXP,
            dataSEXP,
            min_ll,
            selected);
      }
    }
  }

  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
  const double *rt =
      REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
  const TrustedParamMatrix trusted_params(
      paramsSEXP,
      model.component_weight_param_count);

  const auto n_trials = layout.spans.size();
  Rcpp::NumericVector loglik(n_trials, min_ll);
  ExactStepWorkspacePool workspace_pool(exact_plans.size());
  std::vector<double> trial_values;
  std::vector<double> plan_values;

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    if (!trial_is_selected(selected, trial_index)) {
      continue;
    }
    trial_values.clear();
    const auto observed_label_code =
        integer_cell_is_na(label, row)
            ? semantic::kInvalidIndex
            : static_cast<semantic::Index>(label[row]);
    const double observed_rt = rt[row];

    const auto components = resolve_trial_components(
        component_plans_by_code,
        model,
        layout,
        table.component,
        trusted_params,
        trial_index);
    const auto latent_trial = integer_cell_is_na(table.component, row);
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
      if (variant_index == semantic::kInvalidIndex) {
        continue;
      }

      const ObservationProbabilityPlan *plan = nullptr;
      double plan_rt = NA_REAL;
      if (observed_label_code == semantic::kInvalidIndex &&
          Rcpp::NumericVector::is_na(observed_rt)) {
        plan = &component_plan.missing_all_log_plan;
      } else if (observed_label_code != semantic::kInvalidIndex &&
                 Rcpp::NumericVector::is_na(observed_rt)) {
        const auto observed_code = static_cast<std::size_t>(observed_label_code);
        if (observed_code < component_plan.missing_rt_log_plans_by_code.size()) {
          plan = &component_plan.missing_rt_log_plans_by_code[observed_code];
        }
      } else if (observed_label_code != semantic::kInvalidIndex) {
        const auto observed_code = static_cast<std::size_t>(observed_label_code);
        if (observed_code < component_plan.finite_log_plans_by_code.size()) {
          plan = &component_plan.finite_log_plans_by_code[observed_code];
          plan_rt = observed_rt;
        }
      }
      if (plan == nullptr || plan->empty()) {
        continue;
      }

      const int *row_map_ptr = nullptr;
      int row_offset = 0;
      if (latent_trial) {
        const auto &leaf_offsets =
            exact_plans[static_cast<std::size_t>(variant_index)]
                .leaf_row_offsets;
        row_map_ptr = leaf_offsets.data();
        row_offset = static_cast<int>(span.start_row);
      }
      const double value = evaluate_observation_plan_direct(
          exact_plans,
          layout,
          paramsSEXP,
          *plan,
          static_cast<semantic::Index>(trial_index),
          variant_index,
          plan_rt,
          min_ll,
          row_map_ptr,
          row_offset,
          &workspace_pool,
          &plan_values);
      if (std::isfinite(value) && choice.weight > 0.0) {
        trial_values.push_back(std::log(choice.weight) + value);
      }
    }
    const auto value = logsumexp_records(trial_values);
    loglik[static_cast<R_xlen_t>(trial_index)] =
        std::isfinite(value) ? value : min_ll;
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
