#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "exact_interval.hpp"
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
  ObservationLaneBatch lanes;
};

struct ExactResponseScheduleEntry {
  ObservationScheduleEntry observation;
  double lower{0.0};
  double upper{R_PosInf};
  bool denominator{false};
};

struct ExactResponseScheduleGroup {
  semantic::Index variant_index{semantic::kInvalidIndex};
  ExactResponseMeasure measure{ExactResponseMeasure::ObservableResponses};
  const ObservationProbabilityPlan *plan{nullptr};
  std::vector<ExactOutcomeTerm> terms;
  bool complete_outcome_partition{false};
  std::vector<ExactResponseScheduleEntry> entries;
  ObservationLaneBatch lanes;
  std::vector<double> lower;
  std::vector<double> upper;
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
  std::vector<ExactResponseScheduleGroup> response_groups;
  std::vector<std::uint8_t> truncated_trials;
  std::vector<std::vector<ObservationScheduleEntry>> ranked_by_variant;
  std::vector<LatentTrialScheduleEntry> latent_trials;
};

struct ObservationLikelihoodLaneWorkspace {
  explicit ObservationLikelihoodLaneWorkspace(const std::size_t plan_count)
      : observation(plan_count), ranked(plan_count) {
    observation_lanes.reserve(kExactLaneTileSize);
    ranked_lanes.reserve(kExactRankedTrialTileSize);
    active_entries.reserve(kExactLaneTileSize);
    active_response_entries.reserve(kExactLaneTileSize);
    interval_lower.reserve(kExactLaneTileSize);
    interval_upper.reserve(kExactLaneTileSize);
  }

  ObservationLaneWorkspace observation;
  ExactIntervalLaneWorkspace response_interval;
  ExactRankedLaneWorkspacePool ranked;
  ObservationLikelihoodSchedule schedule;
  ObservationLaneBatch observation_lanes;
  std::vector<ExactRankedLane> ranked_lanes;
  std::vector<const ObservationScheduleEntry *> active_entries;
  std::vector<const ExactResponseScheduleEntry *> active_response_entries;
  std::vector<double> scaled_sums;
  std::vector<double> denominator_loglik;
  std::vector<double> denominator_scaled_sums;
  std::vector<double> latent_weights;
  std::vector<double> component_weights;
  std::vector<double> interval_lower;
  std::vector<double> interval_upper;
};

inline void append_exact_outcome_term(
    const ExactVariantPlan &exact_plan,
    const ObservationProbabilityPlan &observation,
    const semantic::Index op_index,
    std::vector<ExactOutcomeTerm> *terms) {
  const auto &op = observation.ops[static_cast<std::size_t>(op_index)];
  if (op.kind == ObservationPlanOpKind::FiniteOutcomeProbability) {
    const auto target = exact_plan.outcome_index_by_code[
        static_cast<std::size_t>(op.semantic_code)];
    for (auto &term : *terms) {
      if (term.target == target) {
        term.weight += op.weight;
        return;
      }
    }
    terms->push_back(ExactOutcomeTerm{target, op.weight});
    return;
  }
  if (op.kind != ObservationPlanOpKind::WeightedSum) {
    return;
  }
  for (semantic::Index child = 0; child < op.children.size; ++child) {
    append_exact_outcome_term(
        exact_plan,
        observation,
        observation.child_ops[
            static_cast<std::size_t>(op.children.offset + child)],
        terms);
  }
}

inline std::vector<ExactOutcomeTerm> exact_terms_from_observation_plan(
    const ExactVariantPlan &exact_plan,
    const ObservationProbabilityPlan &observation) {
  std::vector<ExactOutcomeTerm> terms;
  if (!observation.empty()) {
    append_exact_outcome_term(
        exact_plan, observation, observation.root, &terms);
  }
  return terms;
}

inline bool exact_terms_form_complete_outcome_partition(
    const ExactVariantPlan &exact_plan,
    const std::vector<ExactOutcomeTerm> &terms) {
  if (terms.size() != exact_plan.compiled_outcomes.size()) {
    return false;
  }
  std::vector<std::uint8_t> seen(terms.size(), 0U);
  for (const auto &term : terms) {
    if (term.target < 0 ||
        static_cast<std::size_t>(term.target) >= seen.size() ||
        term.weight != 1.0 ||
        seen[static_cast<std::size_t>(term.target)] != 0U) {
      return false;
    }
    seen[static_cast<std::size_t>(term.target)] = 1U;
  }
  return true;
}

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

inline ExactResponseScheduleGroup &exact_response_schedule_group(
    std::vector<ExactResponseScheduleGroup> *groups,
    const semantic::Index variant_index,
    const ExactResponseMeasure measure,
    const ObservationProbabilityPlan &plan,
    const ExactVariantPlan &exact_plan) {
  for (auto &group : *groups) {
    if (group.variant_index == variant_index &&
        group.measure == measure && group.plan == &plan) {
      return group;
    }
  }
  auto terms = exact_terms_from_observation_plan(exact_plan, plan);
  const bool complete =
      exact_terms_form_complete_outcome_partition(exact_plan, terms);
  groups->push_back(ExactResponseScheduleGroup{
      variant_index, measure, &plan, std::move(terms), complete});
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
  const bool has_observation_bounds = layout.observation.present();
  const auto observation_data =
      has_observation_bounds
          ? read_prepared_observation_data_view(dataSEXP, layout)
          : PreparedObservationDataView{};
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
  schedule->response_groups.clear();
  schedule->truncated_trials.clear();
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
    const auto bounds = has_observation_bounds
                            ? observation_bounds_for_row(observation_data, row)
                            : ObservationBounds{};
    if (bounds.truncates()) {
      if (schedule->truncated_trials.empty()) {
        schedule->truncated_trials.assign(schedule->trial_count, 0U);
      }
      schedule->truncated_trials[trial_index] = 1U;
    }
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

      if (bounds.truncates()) {
        exact_response_schedule_group(
            &schedule->response_groups,
            variant_index,
            ExactResponseMeasure::ObservableResponses,
            component_plan.finite_response_plan,
            exact_plans[static_cast<std::size_t>(variant_index)])
            .entries.push_back(ExactResponseScheduleEntry{
                entry, bounds.trunc_lower, bounds.trunc_upper, true});
      }
      const bool interval_numerator =
          Rcpp::NumericVector::is_na(observed_rt) && bounds.censored();
      if (interval_numerator) {
        const bool known_response =
            observed_label != semantic::kInvalidIndex;
        const ObservationProbabilityPlan *selected_plan = nullptr;
        if (known_response) {
          const auto state_code = missing_rt_observation_state_code(
              component_plan, observed_label);
          if (state_code == semantic::kInvalidIndex ||
              static_cast<std::size_t>(state_code) >=
                  component_plan.probability_plans_by_state_code.size()) {
            continue;
          }
          selected_plan = &observation_probability_plan_for_state(
              component_plan, state_code);
        }
        const auto measure =
            known_response ? ExactResponseMeasure::SelectedResponse
                           : ExactResponseMeasure::ObservableResponses;
        const auto &response_plan =
            known_response ? *selected_plan
                           : component_plan.finite_response_plan;
        auto &group = exact_response_schedule_group(
            &schedule->response_groups,
            variant_index,
            measure,
            response_plan,
            exact_plans[static_cast<std::size_t>(variant_index)]);
        const auto append_interval = [&](const double lower,
                                         const double upper) {
          group.entries.push_back(ExactResponseScheduleEntry{
              entry, lower, upper});
        };
        if (bounds.missingness == 1 || bounds.missingness == 3) {
          append_interval(bounds.trunc_lower, bounds.censor_lower);
        }
        if (bounds.missingness == 2 || bounds.missingness == 3) {
          append_interval(bounds.censor_upper, bounds.trunc_upper);
        }
      } else {
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
      }
      unit_contributions =
          unit_contributions &&
          entry.component_weight_index == semantic::kInvalidIndex;
    }
  }
  schedule->direct_trial_values =
      schedule->response_groups.empty() && unit_contributions &&
      std::all_of(
          contribution_count.begin(),
          contribution_count.end(),
          [](const std::size_t count) { return count == 1U; });

  for (auto &group : schedule->groups) {
    const bool unit_weights = std::all_of(
        group.entries.begin(),
        group.entries.end(),
        [](const ObservationScheduleEntry &entry) {
          return entry.component_weight_index == semantic::kInvalidIndex;
        });
    if (!unit_weights) {
      continue;
    }
    group.lanes.reserve(group.entries.size());
    for (const auto &entry : group.entries) {
      group.lanes.emplace_back(
          entry.row_map, entry.row_offset, entry.observed_time);
    }
    group.lanes.materialize_physical_rows(
        exact_plans[static_cast<std::size_t>(group.variant_index)]
            .leaf_descriptors.size());
  }
  for (auto &group : schedule->response_groups) {
    const bool unit_weights = std::all_of(
        group.entries.begin(),
        group.entries.end(),
        [](const ExactResponseScheduleEntry &entry) {
          return entry.observation.component_weight_index ==
                 semantic::kInvalidIndex;
        });
    if (!unit_weights) {
      continue;
    }
    group.lanes.reserve(group.entries.size());
    group.lower.reserve(group.entries.size());
    group.upper.reserve(group.entries.size());
    for (const auto &entry : group.entries) {
      group.lanes.emplace_back(
          entry.observation.row_map,
          entry.observation.row_offset,
          NA_REAL);
      group.lower.push_back(entry.lower);
      group.upper.push_back(entry.upper);
    }
    group.lanes.materialize_physical_rows(
        exact_plans[static_cast<std::size_t>(group.variant_index)]
            .leaf_descriptors.size());
  }
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
  }
  const ParamMatrixView parameter_matrix(paramsSEXP, schedule.onset);
  const auto trial_count = schedule.trial_count;
  auto &scaled_sums = lane_workspace->scaled_sums;
  auto &denominator_loglik = lane_workspace->denominator_loglik;
  auto &denominator_scaled_sums =
      lane_workspace->denominator_scaled_sums;
  auto &workspace = lane_workspace->observation;
  auto &response_interval = lane_workspace->response_interval;
  auto &ranked_workspaces = lane_workspace->ranked;
  auto &observation_lanes = lane_workspace->observation_lanes;
  auto &ranked_lanes = lane_workspace->ranked_lanes;
  auto &active_entries = lane_workspace->active_entries;
  auto &active_response_entries =
      lane_workspace->active_response_entries;
  auto &interval_lower = lane_workspace->interval_lower;
  auto &interval_upper = lane_workspace->interval_upper;
  auto &values = workspace.group_values;
  const bool direct_trial_values = schedule.direct_trial_values;
  const bool overwrite_unmasked = direct_trial_values && ok == nullptr;

  if (!direct_trial_values) {
    scaled_sums.assign(trial_count, 0.0);
    for (std::size_t trial = 0U; trial < trial_count; ++trial) {
      trial_loglik[trial] = !trial_is_selected(ok, trial)
                                 ? min_ll
                                 : R_NegInf;
    }
    if (!schedule.truncated_trials.empty()) {
      denominator_loglik.assign(trial_count, R_NegInf);
      denominator_scaled_sums.assign(trial_count, 0.0);
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

  const auto accumulate = [&](double *anchors,
                              double *sums,
                              const std::size_t trial,
                              const double value) {
    if (!std::isfinite(value)) {
      return;
    }
    double &anchor = anchors[trial];
    double &sum = sums[trial];
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
    if (entries != nullptr && direct_trial_values) {
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
      const auto &entry = entries != nullptr
                              ? entries[lane]
                              : *active_entries[lane];
      const double value = values[lane];
      const double component_weight =
          entries != nullptr ? 1.0 : entry_weight(entry);
      accumulate(
          trial_loglik,
          scaled_sums.data(),
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
  const auto flush_response_group = [&]
      (const ExactResponseScheduleGroup &group,
       const ObservationLaneBatchView lanes,
       const double *lower,
       const double *upper,
       const ExactResponseScheduleEntry *entries) {
    const auto lane_count = lanes.size;
    if (lane_count == 0U) {
      return;
    }
    const auto &exact_plan =
        exact_plans[static_cast<std::size_t>(group.variant_index)];
    exact_response_probability_between_lanes(
        exact_plan,
        group.measure,
        group.complete_outcome_partition,
        group.terms,
        lanes,
        lower,
        upper,
        &workspace.exact_workspaces.get(
            exact_plans, group.variant_index),
        &response_interval,
        &values);
    for (std::size_t lane = 0U; lane < lane_count; ++lane) {
      const double probability = values[lane];
      if (!(std::isfinite(probability) && probability > 0.0)) {
        continue;
      }
      const auto &entry = entries != nullptr
                              ? entries[lane]
                              : *active_response_entries[lane];
      const double component_weight =
          entries != nullptr ? 1.0 : entry_weight(entry.observation);
      const double value =
          std::log(probability) +
          (component_weight == 1.0 ? 0.0 : std::log(component_weight));
      if (entry.denominator) {
        accumulate(
            denominator_loglik.data(),
            denominator_scaled_sums.data(),
            entry.observation.trial_index,
            value);
      } else {
        accumulate(
            trial_loglik,
            scaled_sums.data(),
            entry.observation.trial_index,
            value);
      }
    }
    active_response_entries.clear();
    interval_lower.clear();
    interval_upper.clear();
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
    if (ok == nullptr && group.lanes.size() == group.entries.size()) {
      for (std::size_t begin = 0U;
           begin < group.entries.size();
           begin += kExactLaneTileSize) {
        const auto lane_count = std::min(
            kExactLaneTileSize, group.entries.size() - begin);
        flush_observation_group(
            group,
            group.lanes.view(parameter_matrix, begin, lane_count),
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

  for (const auto &group : schedule.response_groups) {
    if (ok == nullptr && group.lanes.size() == group.entries.size()) {
      for (std::size_t begin = 0U;
           begin < group.entries.size();
           begin += kExactLaneTileSize) {
        const auto lane_count = std::min(
            kExactLaneTileSize, group.entries.size() - begin);
        flush_response_group(
            group,
            group.lanes.view(parameter_matrix, begin, lane_count),
            group.lower.data() + begin,
            group.upper.data() + begin,
            group.entries.data() + begin);
      }
      continue;
    }
    const auto leaf_count =
        exact_plans[static_cast<std::size_t>(group.variant_index)]
            .leaf_descriptors.size();
    active_response_entries.clear();
    observation_lanes.clear();
    interval_lower.clear();
    interval_upper.clear();
    for (const auto &entry : group.entries) {
      const auto &observation = entry.observation;
      if (!trial_is_selected(ok, observation.trial_index) ||
          !(entry_weight(observation) > 0.0)) {
        continue;
      }
      observation_lanes.emplace_back(
          observation.row_map,
          observation.row_offset,
          NA_REAL);
      active_response_entries.push_back(&entry);
      interval_lower.push_back(entry.lower);
      interval_upper.push_back(entry.upper);
      if (observation_lanes.size() == kExactLaneTileSize) {
        observation_lanes.materialize_physical_rows(leaf_count);
        flush_response_group(
            group,
            observation_lanes.view(parameter_matrix),
            interval_lower.data(),
            interval_upper.data(),
            nullptr);
        observation_lanes.clear();
      }
    }
    observation_lanes.materialize_physical_rows(leaf_count);
    flush_response_group(
        group,
        observation_lanes.view(parameter_matrix),
        interval_lower.data(),
        interval_upper.data(),
        nullptr);
  }

  for (std::size_t variant = 0U;
       variant < schedule.ranked_by_variant.size();
       ++variant) {
    active_entries.clear();
    ranked_lanes.clear();
    const ObservationScheduleEntry *entries = nullptr;
    for (const auto &entry : schedule.ranked_by_variant[variant]) {
      if (!overwrite_unmasked) {
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
      if (!overwrite_unmasked) {
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
      const bool has_numerator =
          std::isfinite(trial_loglik[trial]) && scaled_sums[trial] > 0.0;
      double value =
          has_numerator
              ? (scaled_sums[trial] == 1.0
                     ? trial_loglik[trial]
                     : trial_loglik[trial] + std::log(scaled_sums[trial]))
              : min_ll;
      if (!schedule.truncated_trials.empty() &&
          schedule.truncated_trials[trial] != 0U) {
        const double denominator =
            trial < denominator_loglik.size() &&
                    std::isfinite(denominator_loglik[trial]) &&
                    denominator_scaled_sums[trial] > 0.0
                ? (denominator_scaled_sums[trial] == 1.0
                       ? denominator_loglik[trial]
                       : denominator_loglik[trial] +
                             std::log(denominator_scaled_sums[trial]))
                : R_NegInf;
        value = has_numerator && std::isfinite(denominator)
                    ? value - denominator
                    : min_ll;
      }
      trial_loglik[trial] =
          std::isfinite(value) ? std::max(value, min_ll) : min_ll;
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
