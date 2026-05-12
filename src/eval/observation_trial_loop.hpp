#pragma once

#include <Rcpp.h>

#include <cmath>
#include <unordered_map>
#include <vector>

#include "exact_sequence.hpp"
#include "observation_component_mixture.hpp"
#include "observation_plan_eval.hpp"
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

inline Rcpp::List evaluate_observation_probability_queries_cached(
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

inline SEXP evaluate_observation_likelihood_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observation_is_identity,
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    SEXP expandSEXP,
    const int *ok = nullptr) {
  const auto table = read_prepared_data_view(dataSEXP, layout);

  if (observation_is_identity) {
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
            expandSEXP,
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
  std::unordered_map<
      RtFreeObservationPlanCacheKey,
      std::vector<RtFreeObservationPlanCacheEntry>,
      RtFreeObservationPlanCacheKeyHash> rt_free_plan_cache;

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
      const auto &exact_plan =
          exact_plans[static_cast<std::size_t>(variant_index)];
      const auto leaf_count =
          static_cast<std::size_t>(exact_plan.lowered.program.layout.n_leaves);
      if (latent_trial) {
        const auto &leaf_offsets =
            exact_leaf_row_offsets_by_variant[
                static_cast<std::size_t>(variant_index)];
        row_map_ptr = leaf_offsets.data();
        row_offset = static_cast<int>(span.start_row);
      }
      const bool cacheable_rt_free_plan =
          row_map_ptr == nullptr && !std::isfinite(plan_rt);
      const int first_param_row = static_cast<int>(span.start_row);
      const auto cache_key = RtFreeObservationPlanCacheKey{
          choice.component_code,
          variant_index,
          state_code,
          leaf_count,
          cacheable_rt_free_plan
              ? param_leaf_block_hash(paramsSEXP, first_param_row, leaf_count)
              : 0U};
      double value = 0.0;
      bool cache_hit = false;
      if (cacheable_rt_free_plan) {
        cache_hit = rt_free_observation_cache_lookup(
            rt_free_plan_cache,
            paramsSEXP,
            cache_key,
            first_param_row,
            &value);
      }
      if (!cache_hit) {
        value = evaluate_observation_plan_direct(
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
        if (cacheable_rt_free_plan && std::isfinite(value)) {
          rt_free_plan_cache[cache_key].push_back(
              RtFreeObservationPlanCacheEntry{first_param_row, value});
        }
      }
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

  const double total_loglik = aggregate_trial_loglik(loglik, expandSEXP);

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("n_trials") = static_cast<int>(n_trials));
}

} // namespace detail
} // namespace accumulatr::eval
