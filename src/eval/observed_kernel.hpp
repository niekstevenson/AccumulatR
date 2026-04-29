#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "eval_query.hpp"
#include "exact_sequence.hpp"
#include "observed_plan.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

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
    const int *ok = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
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
    const int *ok = nullptr) {
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    const auto row = static_cast<R_xlen_t>(layout.spans[trial_index].start_row);
    if (integer_cell_is_na(label, row) || Rcpp::NumericVector::is_na(rt[row])) {
      return false;
    }
  }
  return true;
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
  if (values == nullptr || workspace_pool == nullptr) {
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
      const double probability = exact_finite_probability_program_value(
          exact_plan,
          params,
          first_param_row,
          target_idx,
          &workspace);
      value =
          std::isfinite(probability) && probability > 0.0 && op.weight > 0.0
              ? op.weight * probability
              : 0.0;
      break;
    }
    case ObservationPlanOpKind::NoResponseProbability: {
      const double probability = exact_no_response_program_value(
          exact_plan,
          params,
          first_param_row,
          &workspace);
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
  ExactStepWorkspacePool workspace_pool(exact_plans.size());
  std::vector<double> plan_values;

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
    if (observed_code <= 0) {
      continue;
    }
    const auto state_code =
        missing_rt_observation_state_code(component_plan, observed_code);
    const auto &plan =
        observation_probability_plan_for_state(component_plan, state_code);
    const double value = evaluate_observation_plan_direct(
        exact_plans,
        layout,
        paramsSEXP,
        plan,
        static_cast<semantic::Index>(trial_index),
        variant_index,
        NA_REAL,
        std::log(1e-10),
        nullptr,
        0,
        &workspace_pool,
        &plan_values);
    if (std::isfinite(value) && value > 0.0) {
      probability[static_cast<R_xlen_t>(trial_index)] = value;
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
    const int *ok = nullptr) {
  const auto table = read_prepared_data_view(dataSEXP, layout);

  if (observed_identity) {
    const int *label =
        INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
    const double *rt =
        REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
    if (trials_have_observed_components(layout, table.component, ok)) {
      if (identity_trials_are_all_finite(layout, label, rt, ok)) {
        return evaluate_exact_trials_cached(
            exact_variant_index_by_component_code,
            exact_plans,
            layout,
            paramsSEXP,
            dataSEXP,
            min_ll,
            ok);
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
    if (!trial_is_selected(ok, trial_index)) {
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

      const auto state_code = observation_state_code(
          component_plan,
          observed_label_code,
          observed_rt);
      if (state_code == semantic::kInvalidIndex) {
        continue;
      }
      const auto &plan =
          observation_log_plan_for_state(component_plan, state_code);
      const double plan_rt =
          observation_state_uses_rt(component_plan, state_code)
              ? observed_rt
              : NA_REAL;

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
          plan,
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

  const double total_loglik = aggregate_trial_loglik(loglik, layout);

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

} // namespace detail
} // namespace accumulatr::eval
