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
  for (std::size_t trial_index = 0; trial_index < layout.trials.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    if (integer_cell_is_na(
            component,
            static_cast<R_xlen_t>(layout.trials[trial_index].start_row))) {
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
  for (std::size_t trial_index = 0; trial_index < layout.trials.size(); ++trial_index) {
    if (!trial_is_selected(ok, trial_index)) {
      continue;
    }
    const auto row = static_cast<R_xlen_t>(layout.trials[trial_index].start_row);
    if (integer_cell_is_na(label, row) || Rcpp::NumericVector::is_na(rt[row])) {
      return false;
    }
  }
  return true;
}

inline Rcpp::NumericVector evaluate_response_probabilities_cached(
    const ComponentMixturePlan &component_mixture,
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    const std::size_t n_outcomes,
    SEXP paramsSEXP,
    SEXP layoutSEXP) {
  const Rcpp::List layout(layoutSEXP);
  const Rcpp::IntegerVector full_start_rows(layout["full_start_rows"]);
  const Rcpp::IntegerMatrix component_start_rows(layout["component_start_rows"]);
  const auto n_trials = static_cast<std::size_t>(full_start_rows.size());

  Rcpp::NumericVector probability(static_cast<R_xlen_t>(n_outcomes));
  if (n_trials == 0U || n_outcomes == 0U) {
    return probability;
  }

  const TrustedParamMatrix trusted_params(
      paramsSEXP,
      component_mixture.weight_param_count);
  ExactStepWorkspacePool workspace_pool(exact_plans.size());
  std::vector<double> plan_values;
  std::vector<semantic::Index> available_component_codes;

  const auto component_start =
      [&](const std::size_t trial_index,
          const semantic::Index component_code) -> int {
    const auto col = static_cast<int>(component_code - 1);
    if (col < 0 || col >= component_start_rows.ncol()) {
      return NA_INTEGER;
    }
    return component_start_rows(
        static_cast<int>(trial_index),
        col);
  };

  for (std::size_t trial_index = 0; trial_index < n_trials; ++trial_index) {
    available_component_codes.clear();
    const int full_start_row = full_start_rows[static_cast<R_xlen_t>(trial_index)];
    const bool has_full_trial = full_start_row != NA_INTEGER;
    int weight_row = NA_INTEGER;

    if (has_full_trial) {
      available_component_codes = component_mixture.present_component_codes;
      weight_row = full_start_row - 1;
    } else {
      for (std::size_t component_code = 1;
           component_code < component_plans_by_code.size() &&
           component_code < exact_variant_index_by_component_code.size();
           ++component_code) {
        const int start_row =
            component_start(
                trial_index,
                static_cast<semantic::Index>(component_code));
        if (start_row == NA_INTEGER ||
            !component_plans_by_code[component_code].present ||
            exact_variant_index_by_component_code[component_code] ==
                semantic::kInvalidIndex) {
          continue;
        }
        if (weight_row == NA_INTEGER) {
          weight_row = start_row - 1;
        }
        available_component_codes.push_back(
            static_cast<semantic::Index>(component_code));
      }
    }

    if (available_component_codes.empty() || weight_row == NA_INTEGER) {
      continue;
    }

    const auto weights = resolve_component_weights(
        component_mixture,
        available_component_codes,
        trusted_params,
        weight_row);
    for (std::size_t choice_index = 0;
         choice_index < available_component_codes.size();
         ++choice_index) {
      const auto component_code = available_component_codes[choice_index];
      if (!(weights[choice_index] > 0.0)) {
        continue;
      }
      const auto variant_index = resolve_variant_index_by_component_code(
          component_code,
          exact_variant_index_by_component_code);
      if (variant_index == semantic::kInvalidIndex ||
          static_cast<std::size_t>(variant_index) >= exact_plans.size()) {
        continue;
      }
      const auto &component_plan =
          component_plans_by_code[static_cast<std::size_t>(component_code)];
      if (!component_plan.present) {
        continue;
      }

      const int *row_map = nullptr;
      int row_offset = 0;
      int first_param_row = 0;
      if (has_full_trial) {
        const auto &leaf_offsets =
            exact_leaf_row_offsets_by_variant[
                static_cast<std::size_t>(variant_index)];
        row_map = leaf_offsets.data();
        row_offset = full_start_row - 1;
      } else {
        const int start_row = component_start(trial_index, component_code);
        if (start_row == NA_INTEGER) {
          continue;
        }
        first_param_row = start_row - 1;
      }

      for (std::size_t observed_pos = 1U;
           observed_pos <= n_outcomes;
           ++observed_pos) {
        const auto observed_code =
            static_cast<semantic::Index>(observed_pos);
        const auto state_code =
            missing_rt_observation_state_code(component_plan, observed_code);
        if (state_code == semantic::kInvalidIndex ||
            static_cast<std::size_t>(state_code) >=
                component_plan.probability_plans_by_state_code.size()) {
          continue;
        }
        const auto &plan =
            observation_probability_plan_for_state(component_plan, state_code);
        if (plan.empty()) {
          continue;
        }
        const double value = evaluate_observation_plan_at_row(
            exact_plans,
            paramsSEXP,
            nullptr,
            plan,
            variant_index,
            NA_REAL,
            std::log(1e-10),
            row_map,
            row_offset,
            first_param_row,
            &workspace_pool,
            &plan_values);
        if (std::isfinite(value) && value > 0.0) {
          probability[static_cast<R_xlen_t>(observed_pos - 1U)] +=
              weights[choice_index] * value;
        }
      }
    }
  }

  const double inv_trials = 1.0 / static_cast<double>(n_trials);
  for (R_xlen_t i = 0; i < probability.size(); ++i) {
    probability[i] *= inv_trials;
  }
  return probability;
}

template <typename TrialSink>
inline void evaluate_observation_likelihood_trials_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const ComponentMixturePlan &component_mixture,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const int *ok,
    TrialSink &&sink) {
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
  const double *rt =
      REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
  const double *onset =
      layout.onset_col >= 0
          ? REAL(trusted_data_column(dataSEXP, layout.onset_col))
          : nullptr;
  const TrustedParamMatrix trusted_params(
      paramsSEXP,
      component_mixture.weight_param_count);

  ExactStepWorkspacePool workspace_pool(exact_plans.size());
  std::vector<double> trial_values;
  std::vector<double> plan_values;
  std::unordered_map<
      RtFreeObservationPlanCacheKey,
      std::vector<RtFreeObservationPlanCacheEntry>,
      RtFreeObservationPlanCacheKeyHash> rt_free_plan_cache;

  for (std::size_t trial_index = 0; trial_index < layout.trials.size(); ++trial_index) {
    const auto &trial_row = layout.trials[trial_index];
    const auto row = static_cast<R_xlen_t>(trial_row.start_row);
    double trial_loglik = min_ll;
    if (trial_is_selected(ok, trial_index)) {
      trial_values.clear();
      const auto observed_label_code =
          integer_cell_is_na(label, row)
              ? semantic::kInvalidIndex
              : static_cast<semantic::Index>(label[row]);
      const double observed_rt = rt[row];

      const auto components = resolve_trial_components(
          component_mixture,
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
            static_cast<std::size_t>(exact_plan.program.layout.n_leaves);
        if (latent_trial) {
          const auto &leaf_offsets =
              exact_leaf_row_offsets_by_variant[
                  static_cast<std::size_t>(variant_index)];
          row_map_ptr = leaf_offsets.data();
          row_offset = static_cast<int>(trial_row.start_row);
        }
        const bool cacheable_rt_free_plan =
            row_map_ptr == nullptr && !std::isfinite(plan_rt);
        const int first_param_row = static_cast<int>(trial_row.start_row);
        const auto cache_key = RtFreeObservationPlanCacheKey{
            choice.component_code,
            variant_index,
            state_code,
            leaf_count,
            cacheable_rt_free_plan
                ? param_leaf_block_hash(
                      paramsSEXP,
                      onset,
                      first_param_row,
                      leaf_count)
                : 0U};
        double value = 0.0;
        bool cache_hit = false;
        if (cacheable_rt_free_plan) {
          cache_hit = rt_free_observation_cache_lookup(
              rt_free_plan_cache,
              paramsSEXP,
              onset,
              cache_key,
              first_param_row,
              &value);
        }
        if (!cache_hit) {
          value = evaluate_observation_plan_direct(
              exact_plans,
              layout,
              onset,
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
      if (std::isfinite(value)) {
        trial_loglik = value;
      }
    }

    sink(trial_index, trial_loglik);
  }
}

inline void evaluate_observation_likelihood_trial_values_cached(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observation_is_identity,
    const ComponentMixturePlan &component_mixture,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const int *ok,
    double *trial_loglik) {
  const auto table = read_prepared_data_view(dataSEXP, layout);

  if (observation_is_identity) {
    const int *label =
        INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
    const double *rt =
        REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));
    if (trials_have_observed_components(layout, table.component, ok) &&
        identity_trials_are_all_finite(layout, label, rt, ok)) {
      evaluate_exact_trials_cached(
          exact_variant_index_by_component_code,
          exact_plans,
          layout,
          paramsSEXP,
          dataSEXP,
          min_ll,
          ok,
          [&](const std::size_t trial_index, const double value) {
            trial_loglik[trial_index] = value;
          });
      return;
    }
  }

  evaluate_observation_likelihood_trials_cached(
      component_plans_by_code,
      component_mixture,
      exact_variant_index_by_component_code,
      exact_plans,
      exact_leaf_row_offsets_by_variant,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll,
      ok,
      [&](const std::size_t trial_index, const double value) {
        trial_loglik[trial_index] = value;
      });
}

} // namespace detail
} // namespace accumulatr::eval
