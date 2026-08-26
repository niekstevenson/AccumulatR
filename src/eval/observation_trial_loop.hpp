#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "exact_sequence.hpp"
#include "observation_component_mixture.hpp"
#include "observation_plan_eval.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ObservationScheduleEntry {
  std::size_t trial_index{0U};
  R_xlen_t row{0};
  const int *row_map{nullptr};
  int row_offset{0};
  double observed_time{NA_REAL};
  semantic::Index component_weight_index{semantic::kInvalidIndex};
  int rank_count{0};
};

struct ObservationScheduleGroup {
  semantic::Index variant_index{semantic::kInvalidIndex};
  const ObservationProbabilityPlan *plan{nullptr};
  std::vector<ObservationScheduleEntry> entries;
  ObservationLaneBatch direct_lanes;
};

struct LatentTrialScheduleEntry {
  std::size_t trial_index{0U};
  int parameter_row{0};
};

struct ObservationLikelihoodSchedule {
  bool matches(SEXP dataSEXP) const {
    return data != R_NilValue && data == dataSEXP;
  }

  void identify(SEXP dataSEXP) {
    data_owner = dataSEXP;
    data = dataSEXP;
  }

  Rcpp::RObject data_owner;
  SEXP data{R_NilValue};
  std::size_t trial_count{0U};
  std::size_t component_count{0U};
  bool direct_trial_values{false};
  const double *onset{nullptr};
  ExactTrialColumns ranked_columns;
  std::vector<ObservationScheduleGroup> groups;
  std::vector<std::vector<ObservationScheduleEntry>> ranked_by_variant;
  std::vector<LatentTrialScheduleEntry> latent_trials;
};

struct ObservationLikelihoodLaneWorkspace {
  explicit ObservationLikelihoodLaneWorkspace(const std::size_t plan_count)
      : observation(plan_count), ranked(plan_count) {
    observation_lanes.reserve(kExactLaneTileSize);
    ranked_lanes.reserve(kExactRankedTrialTileSize);
    active_entries.reserve(kExactLaneTileSize);
  }

  ObservationLaneWorkspace observation;
  ExactRankedLaneWorkspacePool ranked;
  ObservationLikelihoodSchedule schedule;
  ObservationLaneBatch observation_lanes;
  std::vector<ExactRankedLane> ranked_lanes;
  std::vector<const ObservationScheduleEntry *> active_entries;
  std::vector<double> scaled_sums;
  std::vector<double> latent_weights;
  std::vector<double> component_weights;
};

inline ObservationScheduleGroup &observation_schedule_group(
    std::vector<ObservationScheduleGroup> *groups,
    const semantic::Index variant_index,
    const ObservationProbabilityPlan &plan) {
  for (auto &group : *groups) {
    if (group.variant_index == variant_index && group.plan == &plan) {
      return group;
    }
  }
  groups->push_back(ObservationScheduleGroup{variant_index, &plan, {}});
  return groups->back();
}

inline Rcpp::NumericVector evaluate_response_probabilities_cached(
    const ComponentMixturePlan &component_mixture,
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    const std::size_t n_outcomes,
    SEXP paramsSEXP,
    SEXP layoutSEXP,
    ObservationLaneWorkspace *workspace) {
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
  const ParamMatrixView parameter_matrix(paramsSEXP);
  workspace->begin();
  auto &groups = workspace->groups;
  auto &values = workspace->group_values;
  auto &available_component_codes = workspace->component_codes;
  auto &weights = workspace->component_weights;

  const auto flush_group = [&](ObservationLaneGroup *group) {
    if (group->lanes.empty()) {
      return;
    }
    evaluate_observation_lane_group(
        exact_plans,
        std::log(1e-10),
        *group,
        parameter_matrix,
        workspace,
        &values);
    for (std::size_t lane_index = 0U;
         lane_index < group->lanes.size();
         ++lane_index) {
      const double value = values[lane_index];
      if (std::isfinite(value) && value > 0.0) {
        probability[static_cast<R_xlen_t>(
            group->destination_indices[lane_index])] +=
            group->component_weights[lane_index] * value;
      }
    }
    group->clear();
  };

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

    resolve_component_weights(
        component_mixture,
        available_component_codes,
        trusted_params,
        weight_row,
        &weights);
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
        row_offset = start_row - 1;
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
        auto &group = find_observation_lane_group(
            &groups, variant_index, plan);
        group.lanes.emplace_back(row_map, row_offset, NA_REAL);
        group.destination_indices.push_back(observed_pos - 1U);
        group.component_weights.push_back(weights[choice_index]);
        if (group.lanes.size() == kExactLaneTileSize) {
          flush_group(&group);
        }
      }
    }
  }

  for (auto &group : groups) {
    flush_group(&group);
  }

  const double inv_trials = 1.0 / static_cast<double>(n_trials);
  for (R_xlen_t i = 0; i < probability.size(); ++i) {
    probability[i] *= inv_trials;
  }
  return probability;
}

inline void build_observation_likelihood_schedule(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observation_is_identity,
    const ComponentMixturePlan &component_mixture,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    SEXP dataSEXP,
    ObservationLikelihoodSchedule *schedule) {
  const auto layout = read_prepared_trial_layout(dataSEXP);
  const auto table = read_prepared_data_view(dataSEXP, layout);
  const int *label =
      INTEGER(trusted_data_column(dataSEXP, layout.label_cols[1]));
  const double *rt =
      REAL(trusted_data_column(dataSEXP, layout.time_cols[1]));

  schedule->trial_count = layout.trials.size();
  schedule->component_count = component_mixture.present_component_codes.size();
  schedule->onset = layout.onset_col >= 0
                        ? REAL(trusted_data_column(
                              dataSEXP, layout.onset_col))
                        : nullptr;
  schedule->groups.clear();
  schedule->ranked_by_variant.assign(exact_plans.size(), {});
  schedule->latent_trials.clear();
  schedule->ranked_columns = {};
  std::vector<std::size_t> contribution_count(
      schedule->trial_count, 0U);
  bool unit_contributions = true;

  const int *second_rank_labels = nullptr;
  const double *second_rank_times = nullptr;
  if (observation_is_identity && layout.max_rank > 1) {
    schedule->ranked_columns = make_exact_trial_columns(dataSEXP, layout);
    second_rank_labels = schedule->ranked_columns.labels[2U];
    second_rank_times = schedule->ranked_columns.times[2U];
  }

  for (std::size_t trial_index = 0U;
       trial_index < schedule->trial_count;
       ++trial_index) {
    const auto trial_row = layout.trials[trial_index];
    const auto row = static_cast<R_xlen_t>(trial_row.start_row);
    const auto observed_label =
        integer_cell_is_na(label, row)
            ? semantic::kInvalidIndex
            : static_cast<semantic::Index>(label[row]);
    const double observed_rt = rt[row];
    const bool latent_trial = integer_cell_is_na(table.component, row);

    int rank_count = 0;
    if (second_rank_labels != nullptr &&
        (!integer_cell_is_na(second_rank_labels, row) ||
         !Rcpp::NumericVector::is_na(second_rank_times[row]))) {
      for (int rank = 1; rank <= layout.max_rank; ++rank) {
        const auto *rank_labels =
            schedule->ranked_columns.labels[static_cast<std::size_t>(rank)];
        const auto *rank_times =
            schedule->ranked_columns.times[static_cast<std::size_t>(rank)];
        if (integer_cell_is_na(rank_labels, row) &&
            Rcpp::NumericVector::is_na(rank_times[row])) {
          break;
        }
        ++rank_count;
      }
    }

    if (latent_trial && schedule->component_count > 1U) {
      schedule->latent_trials.push_back(
          LatentTrialScheduleEntry{trial_index, static_cast<int>(row)});
    }
    const std::size_t choice_count =
        latent_trial ? component_mixture.present_component_codes.size() : 1U;
    for (std::size_t choice_index = 0U;
         choice_index < choice_count;
         ++choice_index) {
      const auto component_code =
          latent_trial
              ? component_mixture.present_component_codes[choice_index]
              : static_cast<semantic::Index>(table.component[row]);
      if (component_code <= 0 ||
          component_code >= static_cast<semantic::Index>(
                                component_plans_by_code.size())) {
        continue;
      }
      const auto &component_plan =
          component_plans_by_code[static_cast<std::size_t>(component_code)];
      if (!component_plan.present) {
        continue;
      }
      const auto variant_index = resolve_variant_index_by_component_code(
          component_code, exact_variant_index_by_component_code);
      if (variant_index == semantic::kInvalidIndex ||
          static_cast<std::size_t>(variant_index) >= exact_plans.size()) {
        continue;
      }

      const int *row_map = nullptr;
      if (latent_trial) {
        row_map = exact_leaf_row_offsets_by_variant[
                      static_cast<std::size_t>(variant_index)]
                      .data();
      }
      ObservationScheduleEntry entry{
          trial_index,
          row,
          row_map,
          static_cast<int>(trial_row.start_row),
          observed_rt,
          latent_trial && schedule->component_count > 1U
              ? static_cast<semantic::Index>(choice_index)
              : semantic::kInvalidIndex,
          rank_count};

      if (rank_count > 1) {
        schedule->ranked_by_variant[
            static_cast<std::size_t>(variant_index)]
            .push_back(entry);
        ++contribution_count[trial_index];
        unit_contributions =
            unit_contributions &&
            entry.component_weight_index == semantic::kInvalidIndex;
        continue;
      }

      const auto state_code = observation_state_code(
          component_plan, observed_label, observed_rt);
      if (state_code == semantic::kInvalidIndex) {
        continue;
      }
      const auto &plan =
          observation_log_plan_for_state(component_plan, state_code);
      entry.observed_time =
          observation_state_uses_rt(component_plan, state_code)
              ? observed_rt
              : NA_REAL;
      observation_schedule_group(
          &schedule->groups, variant_index, plan)
          .entries.push_back(entry);
      ++contribution_count[trial_index];
      unit_contributions =
          unit_contributions &&
          entry.component_weight_index == semantic::kInvalidIndex;
    }
  }
  schedule->direct_trial_values =
      unit_contributions &&
      std::all_of(
          contribution_count.begin(),
          contribution_count.end(),
          [](const std::size_t count) { return count == 1U; });
  schedule->identify(dataSEXP);
}

inline void evaluate_observation_likelihood_trial_values_lanes(
    const std::vector<ComponentObservationPlan> &component_plans_by_code,
    const bool observation_is_identity,
    const ComponentMixturePlan &component_mixture,
    const std::vector<semantic::Index> &exact_variant_index_by_component_code,
    const std::vector<ExactVariantPlan> &exact_plans,
    const std::vector<std::vector<int>> &exact_leaf_row_offsets_by_variant,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const int *ok,
    ObservationLikelihoodLaneWorkspace *lane_workspace,
    double *trial_loglik) {
  auto &schedule = lane_workspace->schedule;
  bool rebuilt_schedule = false;
  if (!schedule.matches(dataSEXP)) {
    ObservationLikelihoodSchedule replacement;
    build_observation_likelihood_schedule(
        component_plans_by_code,
        observation_is_identity,
        component_mixture,
        exact_variant_index_by_component_code,
        exact_plans,
        exact_leaf_row_offsets_by_variant,
        dataSEXP,
        &replacement);
    schedule = std::move(replacement);
    rebuilt_schedule = true;
  }
  const ParamMatrixView parameter_matrix(paramsSEXP, schedule.onset);
  if (rebuilt_schedule && schedule.direct_trial_values) {
    for (auto &group : schedule.groups) {
      group.direct_lanes.clear();
      group.direct_lanes.reserve(group.entries.size());
      for (const auto &entry : group.entries) {
        group.direct_lanes.emplace_back(
            entry.row_map, entry.row_offset, entry.observed_time);
      }
      group.direct_lanes.materialize_physical_rows(
          exact_plans[static_cast<std::size_t>(group.variant_index)]
              .leaf_descriptors.size());
    }
  }
  const auto trial_count = schedule.trial_count;
  auto &scaled_sums = lane_workspace->scaled_sums;
  auto &workspace = lane_workspace->observation;
  auto &ranked_workspaces = lane_workspace->ranked;
  auto &observation_lanes = lane_workspace->observation_lanes;
  auto &ranked_lanes = lane_workspace->ranked_lanes;
  auto &active_entries = lane_workspace->active_entries;
  auto &values = workspace.group_values;
  const bool direct_trial_values = schedule.direct_trial_values;
  const bool direct_unmasked = direct_trial_values && ok == nullptr;

  if (!direct_trial_values) {
    scaled_sums.assign(trial_count, 0.0);
    for (std::size_t trial = 0U; trial < trial_count; ++trial) {
      trial_loglik[trial] = !trial_is_selected(ok, trial)
                                 ? min_ll
                                 : R_NegInf;
    }
  }

  auto &latent_weights = lane_workspace->latent_weights;
  auto &component_weights = lane_workspace->component_weights;
  if (!direct_trial_values) {
    const TrustedParamMatrix trusted_params(
        paramsSEXP, component_mixture.weight_param_count);
    latent_weights.assign(
        schedule.latent_trials.empty()
            ? 0U
            : trial_count * schedule.component_count,
        0.0);
    for (const auto &latent : schedule.latent_trials) {
      if (!trial_is_selected(ok, latent.trial_index)) {
        continue;
      }
      resolve_component_weights(
          component_mixture,
          component_mixture.present_component_codes,
          trusted_params,
          latent.parameter_row,
          &component_weights);
      std::copy_n(
          component_weights.begin(),
          schedule.component_count,
          latent_weights.begin() +
              static_cast<std::ptrdiff_t>(
                  latent.trial_index * schedule.component_count));
    }
  }

  const auto entry_weight = [&](const ObservationScheduleEntry &entry) {
    if (entry.component_weight_index == semantic::kInvalidIndex) {
      return 1.0;
    }
    return latent_weights[
        entry.trial_index * schedule.component_count +
        static_cast<std::size_t>(entry.component_weight_index)];
  };

  const auto accumulate = [&](const std::size_t trial,
                              const double value) {
    if (!std::isfinite(value)) {
      return;
    }
    double &anchor = trial_loglik[trial];
    double &sum = scaled_sums[trial];
    if (!std::isfinite(anchor)) {
      anchor = value;
      sum = 1.0;
    } else if (value > anchor) {
      sum = sum * std::exp(anchor - value) + 1.0;
      anchor = value;
    } else {
      sum += std::exp(value - anchor);
    }
  };
  const auto consume_values = [&](const std::size_t lane_count,
                                  const ObservationScheduleEntry *entries) {
    if (entries != nullptr) {
      for (std::size_t lane = 0U; lane < lane_count; ++lane) {
        const double value = values[lane];
        trial_loglik[entries[lane].trial_index] =
            std::isfinite(value) ? value : min_ll;
      }
      return;
    }
    if (direct_trial_values) {
      for (std::size_t lane = 0U; lane < lane_count; ++lane) {
        const double value = values[lane];
        trial_loglik[active_entries[lane]->trial_index] =
            std::isfinite(value) ? value : min_ll;
      }
      return;
    }
    for (std::size_t lane = 0U; lane < lane_count; ++lane) {
      const auto &entry = *active_entries[lane];
      const double value = values[lane];
      const double component_weight = entry_weight(entry);
      accumulate(
          entry.trial_index,
          component_weight == 1.0
              ? value
              : std::log(component_weight) + value);
    }
  };
  const auto flush_observation_group = [&]
      (const ObservationScheduleGroup &group,
       const ObservationLaneBatchView lanes,
       const ObservationScheduleEntry *entries) {
    const auto lane_count = lanes.size;
    if (lane_count == 0U) {
      return;
    }
    evaluate_observation_lanes(
        exact_plans,
        min_ll,
        group.variant_index,
        *group.plan,
        lanes,
        &workspace,
        &values);
    consume_values(lane_count, entries);
    active_entries.clear();
  };
  const auto flush_ranked_group = [&]
      (const std::size_t variant,
       const ObservationScheduleEntry *entries) {
    const auto lane_count = ranked_lanes.size();
    if (lane_count == 0U) {
      return;
    }
    exact_ranked_loglik_lanes(
        exact_plans[variant],
        ranked_lanes.data(),
        lane_count,
        min_ll,
        &workspace.exact_workspaces.get(
            exact_plans, static_cast<semantic::Index>(variant)),
        &ranked_workspaces.get(
            exact_plans, static_cast<semantic::Index>(variant)),
        &values);
    consume_values(lane_count, entries);
    ranked_lanes.clear();
    active_entries.clear();
  };

  for (const auto &group : schedule.groups) {
    if (direct_unmasked) {
      for (std::size_t begin = 0U;
           begin < group.entries.size();
           begin += kExactLaneTileSize) {
        const auto lane_count = std::min(
            kExactLaneTileSize, group.entries.size() - begin);
        flush_observation_group(
            group,
            group.direct_lanes.view(parameter_matrix, begin, lane_count),
            group.entries.data() + begin);
      }
      continue;
    }
    active_entries.clear();
    observation_lanes.clear();
    std::size_t lane_count = 0U;
    for (const auto &entry : group.entries) {
      if (!trial_is_selected(ok, entry.trial_index)) {
        if (direct_trial_values) {
          trial_loglik[entry.trial_index] = min_ll;
        }
        continue;
      }
      if (!direct_trial_values && !(entry_weight(entry) > 0.0)) {
        continue;
      }
      observation_lanes.emplace_back(
          entry.row_map, entry.row_offset, entry.observed_time);
      active_entries.push_back(&entry);
      ++lane_count;
      if (lane_count == kExactLaneTileSize) {
        flush_observation_group(
            group, observation_lanes.view(parameter_matrix), nullptr);
        observation_lanes.clear();
        lane_count = 0U;
      }
    }
    flush_observation_group(
        group, observation_lanes.view(parameter_matrix), nullptr);
  }

  for (std::size_t variant = 0U;
       variant < schedule.ranked_by_variant.size();
       ++variant) {
    active_entries.clear();
    ranked_lanes.clear();
    const ObservationScheduleEntry *entries = nullptr;
    for (const auto &entry : schedule.ranked_by_variant[variant]) {
      if (!direct_unmasked) {
        if (!trial_is_selected(ok, entry.trial_index)) {
          if (direct_trial_values) {
            trial_loglik[entry.trial_index] = min_ll;
          }
          continue;
        }
        if (!direct_trial_values && !(entry_weight(entry) > 0.0)) {
          continue;
        }
      } else if (ranked_lanes.empty()) {
        entries = &entry;
      }
      ranked_lanes.emplace_back(
          parameter_matrix,
          entry.row_map,
          entry.row_offset,
          ExactTrialView{
              entry.row, entry.rank_count, &schedule.ranked_columns});
      if (!direct_unmasked) {
        active_entries.push_back(&entry);
      }
      if (ranked_lanes.size() == kExactRankedTrialTileSize) {
        flush_ranked_group(variant, entries);
        entries = nullptr;
      }
    }
    flush_ranked_group(variant, entries);
  }

  if (!direct_trial_values) {
    for (std::size_t trial = 0; trial < trial_count; ++trial) {
      if (!trial_is_selected(ok, trial)) {
        continue;
      }
      trial_loglik[trial] =
          std::isfinite(trial_loglik[trial]) && scaled_sums[trial] > 0.0
              ? (scaled_sums[trial] == 1.0
                     ? trial_loglik[trial]
                     : trial_loglik[trial] + std::log(scaled_sums[trial]))
              : min_ll;
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
