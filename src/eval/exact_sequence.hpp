#pragma once

#include <limits>
#include <utility>

#include "eval_query.hpp"
#include "exact_step_distribution.hpp"
#include "trial_data.hpp"
#include "../compile/exact_evaluation_program_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactTrialColumns {
  std::vector<const int *> labels;
  std::vector<const double *> times;
};

struct ExactTrialView {
  R_xlen_t row{0};
  int rank_count{0};
  const ExactTrialColumns *columns{nullptr};
};

inline ExactTrialColumns make_exact_trial_columns(
    SEXP dataSEXP,
    const PreparedTrialLayout &layout) {
  ExactTrialColumns columns;
  const int max_rank = layout.max_rank;
  columns.labels.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  columns.times.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  for (int rank = 1; rank <= max_rank; ++rank) {
    columns.labels[static_cast<std::size_t>(rank)] =
        INTEGER(VECTOR_ELT(
            dataSEXP,
            layout.label_cols[static_cast<std::size_t>(rank)]));
    columns.times[static_cast<std::size_t>(rank)] =
        REAL(VECTOR_ELT(
            dataSEXP,
            layout.time_cols[static_cast<std::size_t>(rank)]));
  }
  return columns;
}

inline semantic::Index exact_trial_view_outcome_code(const ExactTrialView &view,
                                                     const std::size_t rank_idx) {
  return view.columns->labels[rank_idx + 1U][view.row];
}

inline double exact_trial_view_rt(const ExactTrialView &view,
                                  const std::size_t rank_idx) {
  return view.columns->times[rank_idx + 1U][view.row];
}

inline void build_exact_plan_cache(
    const compile::CompiledModel &compiled,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_component_codes,
    const std::size_t n_outcome_codes,
    std::vector<semantic::Index> *variant_index_by_component_code,
    std::vector<ExactVariantPlan> *plans,
    std::vector<ExactComplexityMetrics> *complexity_metrics_by_variant = nullptr) {
  variant_index_by_component_code->assign(
      n_component_codes + 1U,
      semantic::kInvalidIndex);
  plans->clear();
  plans->reserve(compiled.variants.size());
  if (complexity_metrics_by_variant != nullptr) {
    complexity_metrics_by_variant->clear();
    complexity_metrics_by_variant->reserve(compiled.variants.size());
  }
  for (const auto &variant : compiled.variants) {
    const auto plan_index = static_cast<semantic::Index>(plans->size());
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared component code for '" +
          variant.component_id + "'");
    }
    (*variant_index_by_component_code)[static_cast<std::size_t>(component_it->second)] =
        plan_index;
    auto evaluation_program =
        accumulatr::compile::lower_exact_evaluation_program(
            variant,
            outcome_code_by_label);
    ExactComplexityMetrics metrics;
    plans->push_back(
        make_exact_variant_plan(
            std::move(evaluation_program),
            n_outcome_codes,
            complexity_metrics_by_variant == nullptr ? nullptr : &metrics));
    if (complexity_metrics_by_variant != nullptr) {
      complexity_metrics_by_variant->push_back(metrics);
    }
  }
}

inline bool exact_has_single_unweighted_trigger_state(
    const ExactVariantPlan &plan) noexcept {
  const auto &states = plan.trigger_state_table.states;
  return states.size() == 1U && states.front().weight_terms.empty();
}

inline void exact_unranked_target_density_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const semantic::Index target_idx,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *out) {
  const auto lane_count = lanes.size;
  const auto &trigger_states = plan.trigger_state_table.states;
  out->resize(lane_count);
  if (lane_count == 0U) {
    return;
  }
  auto &values = workspace->input_values;
  if (exact_has_single_unweighted_trigger_state(plan)) {
    const auto &compiled_state = trigger_states.front();
    const auto *shared_started =
        exact_compiled_trigger_shared_started(plan, compiled_state);
    auto &frame = workspace->prepare_initial(
        lanes,
        shared_started);
    evaluate_exact_outcome_density_lanes(
        plan,
        &frame,
        lane_count,
        target_idx,
        workspace,
        out);
    return;
  }

  for (std::size_t state_index = 0U;
       state_index < trigger_states.size();
       ++state_index) {
    const auto &compiled_state = trigger_states[state_index];
    const bool fixed_weight = compiled_state.weight_terms.empty();
    if (!fixed_weight) {
      exact_compiled_trigger_state_weights_lanes(
          plan,
          lanes,
          lane_count,
          compiled_state,
          &workspace->trigger_weights);
    }
    const auto *shared_started =
        exact_compiled_trigger_shared_started(plan, compiled_state);
    auto &frame = workspace->prepare_initial(
        lanes,
        shared_started);
    evaluate_exact_outcome_density_lanes(
        plan,
        &frame,
        lane_count,
        target_idx,
        workspace,
        &values);
    for (std::size_t lane = 0U; lane < lane_count; ++lane) {
      const double value = values[lane];
      const double weight = fixed_weight
                                ? 1.0
                                : workspace->trigger_weights[lane];
      double &total = (*out)[lane];
      if (state_index == 0U) {
        total = weight > 0.0 ? weight * value : 0.0;
      } else if (weight > 0.0) {
        total += weight * value;
      }
      if (state_index + 1U == trigger_states.size() &&
          (!std::isfinite(total) || !(total > 0.0))) {
        total = 0.0;
      }
    }
  }
}

inline void exact_finite_outcome_probability_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
    const semantic::Index target_idx,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *out) {
  const auto lane_count = lanes.size;
  out->assign(lane_count, 0.0);
  if (lane_count == 0U) {
    return;
  }
  const auto &tail = quadrature::canonical_tail_batch().nodes;
  auto &input_lanes = workspace->input_lanes;
  auto &positions = workspace->input_positions;
  auto &weights = workspace->input_weights;
  auto &values = workspace->expanded_values;
  auto &density_totals = workspace->finite_density_totals;

  for (std::size_t tile_start = 0U;
       tile_start < lane_count;
       tile_start += kExactLaneTileSize) {
    const auto tile_count =
        std::min(kExactLaneTileSize, lane_count - tile_start);
    density_totals.assign(
        quadrature::kDefaultTailOrder * tile_count, 0.0);
    for (const auto &compiled_state : plan.trigger_state_table.states) {
      const bool fixed_weight = compiled_state.weight_terms.empty();
      if (!fixed_weight) {
        exact_compiled_trigger_state_weights_lanes(
            plan,
            lanes + tile_start,
            tile_count,
            compiled_state,
            &workspace->trigger_weights);
      }
      const auto *shared_started =
          exact_compiled_trigger_shared_started(plan, compiled_state);
      input_lanes.clear();
      positions.clear();
      weights.clear();
      for (std::size_t tile_lane = 0U;
           tile_lane < tile_count;
           ++tile_lane) {
        const auto lane = tile_start + tile_lane;
        const double weight = fixed_weight ? 1.0 : workspace->trigger_weights[tile_lane];
        if (!(weight > 0.0)) {
          continue;
        }
        input_lanes.emplace_back(lanes.row_maps[lane], lanes.row_offsets[lane], NA_REAL);
        positions.push_back(tile_lane);
        weights.push_back(weight);
      }
      const auto active_count = input_lanes.size();
      if (active_count == 0U) {
        continue;
      }
      workspace->bind_initial_sources(input_lanes.view(*lanes.matrix), shared_started);
      const auto times_per_block = std::max<std::size_t>(
          1U, kExactExpandedLaneTileSize / active_count);
      for (std::size_t q_start = 0U;
           q_start < quadrature::kDefaultTailOrder;
           q_start += times_per_block) {
        const auto q_count = std::min(
            times_per_block,
            quadrature::kDefaultTailOrder - q_start);
        auto &frame = workspace->prepare_repeated_times(
            active_count,
            tail.nodes.data() + q_start,
            q_count);
        const auto expanded_count = active_count * q_count;
        evaluate_exact_outcome_density_lanes(
            plan,
            &frame,
            expanded_count,
            target_idx,
            workspace,
            &values);
        for (std::size_t q = 0; q < q_count; ++q) {
          const auto total_offset = (q_start + q) * tile_count;
          const auto value_offset = q * active_count;
          for (std::size_t lane = 0; lane < active_count; ++lane) {
            const double value = values[value_offset + lane];
            if (std::isfinite(value) && value > 0.0) {
              density_totals[total_offset + positions[lane]] +=
                  weights[lane] * value;
            }
          }
        }
      }
    }
    for (std::size_t q = 0; q < quadrature::kDefaultTailOrder; ++q) {
      const auto total_offset = q * tile_count;
      for (std::size_t lane = 0; lane < tile_count; ++lane) {
        const double value = density_totals[total_offset + lane];
        if (std::isfinite(value) && value > 0.0) {
          (*out)[tile_start + lane] += tail.weights[q] * value;
        }
      }
    }
  }
  for (auto &value : *out) {
    value = std::isfinite(value) ? clamp_probability(value) : 0.0;
  }
}

inline void exact_terminal_no_response_probability_lanes(
    const ExactVariantPlan &plan,
    const ObservationLaneBatchView lanes,
  ExactStepLaneWorkspace *workspace,
  std::vector<double> *out) {
  const auto lane_count = lanes.size;
  out->assign(lane_count, 0.0);
  auto &products = workspace->input_weights;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    exact_compiled_trigger_state_weights_lanes(
        plan, lanes, lane_count, compiled_state, &products);
    const auto *shared_started =
        exact_compiled_trigger_shared_started(plan, compiled_state);
    for (const auto leaf_index : plan.no_response.leaf_indices) {
      for (std::size_t lane_index = 0; lane_index < lane_count; ++lane_index) {
        if (!(products[lane_index] > 0.0)) {
          continue;
        }
        products[lane_index] *= clamp_probability(exact_leaf_q_for_trigger_state(
            plan.leaf_trigger_index,
            shared_started,
            leaf_index,
            lanes.q(leaf_index, lane_index)));
      }
    }
    for (std::size_t lane_index = 0; lane_index < lane_count; ++lane_index) {
      (*out)[lane_index] += products[lane_index];
    }
  }
  for (auto &value : *out) {
    value = std::isfinite(value) ? clamp_probability(value) : 0.0;
  }
}

inline void advance_exact_sequence_state(
    ExactSequenceState &state,
    const ExactCompiledTransitionPlan &transition,
    const double observed_time,
    const std::vector<double> &ready_expr_normalizers) {
  state.has_history = true;
  state.lower_bound = observed_time;
  if (transition.release_source_id != semantic::kInvalidIndex) {
    state.exact_times[static_cast<std::size_t>(transition.release_source_id)] =
        observed_time;
  }
  for (const auto source_id : transition.readiness_source_ids) {
    auto &upper = state.upper_bounds[static_cast<std::size_t>(source_id)];
    upper = std::isfinite(upper) ? std::min(upper, observed_time)
                                 : observed_time;
  }
  for (std::size_t i = 0; i < transition.readiness_expr_ids.size(); ++i) {
    const auto expr_id = transition.readiness_expr_ids[i];
    const double normalizer = ready_expr_normalizers[i];
    if (!(normalizer > 0.0) || !std::isfinite(normalizer)) {
      continue;
    }
    auto &upper =
        state.expr_upper_bounds[static_cast<std::size_t>(expr_id)];
    if (std::isfinite(upper) && upper <= observed_time) {
      continue;
    }
    upper = observed_time;
    state.expr_upper_normalizers[static_cast<std::size_t>(expr_id)] =
        normalizer;
  }
}

inline bool exact_sequence_states_equal(const ExactSequenceState &lhs,
                                        const ExactSequenceState &rhs) {
  if (lhs.has_history != rhs.has_history ||
      lhs.lower_bound != rhs.lower_bound) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.exact_times.size(); ++i) {
    const bool lhs_na = std::isnan(lhs.exact_times[i]);
    const bool rhs_na = std::isnan(rhs.exact_times[i]);
    if (lhs_na || rhs_na) {
      if (lhs_na != rhs_na) {
        return false;
      }
      continue;
    }
    if (lhs.exact_times[i] != rhs.exact_times[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.upper_bounds.size(); ++i) {
    if (lhs.upper_bounds[i] != rhs.upper_bounds[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.expr_upper_bounds.size(); ++i) {
    if (lhs.expr_upper_bounds[i] != rhs.expr_upper_bounds[i]) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.expr_upper_normalizers.size(); ++i) {
    if (lhs.expr_upper_normalizers[i] != rhs.expr_upper_normalizers[i]) {
      return false;
    }
  }
  return true;
}

inline ExactSequenceState &ranked_sequence_state_slot(
    std::vector<ExactSequenceState> *states,
    const ExactVariantPlan &plan,
    const std::size_t index) {
  while (states->size() <= index) {
    states->push_back(make_exact_sequence_state(plan));
  }
  return (*states)[index];
}

struct ExactRankedLane {
  ExactRankedLane(const ParamMatrixView &parameter_matrix,
                  const int *row_map,
                  const int row_offset,
                  const ExactTrialView &observation_)
      : params(parameter_matrix, row_map, row_offset),
        observation(observation_) {}

  ParamView params;
  ExactTrialView observation;
};

struct ExactRankedLaneState {
  ExactTriggerState trigger;
  std::vector<std::uint8_t> used_outcomes;
  semantic::Index pending_target{semantic::kInvalidIndex};
  double pending_time{0.0};
  bool trigger_active{false};
  bool rank_pending{false};
};

struct ExactRankedLaneFrontierEntry {
  std::size_t trial_index{0U};
  double probability{0.0};
  semantic::Index state_index{semantic::kInvalidIndex};
};

struct ExactRankedWorkItem {
  std::size_t trial_index{0U};
  semantic::Index state_index{semantic::kInvalidIndex};
  double frontier_probability{0.0};
};

constexpr std::size_t kExactRankedTrialTileSize = 256U;

struct ExactRankedLaneWorkspace {
  explicit ExactRankedLaneWorkspace(const ExactVariantPlan &plan)
      : candidate_state(make_exact_sequence_state(plan)),
        trial_states(kExactRankedTrialTileSize),
        work_by_outcome(plan.compiled_outcomes.size()) {
    totals.reserve(kExactRankedTrialTileSize);
    next_frontier_head_by_trial.resize(
        kExactRankedTrialTileSize, std::numeric_limits<std::size_t>::max());
    conditional_totals.resize(kExactRankedTrialTileSize, 0.0);
  }

  void ensure_trials(const std::size_t lane_count) {
    totals.assign(lane_count, 0.0);
  }

  ExactSequenceState candidate_state;
  std::vector<ExactRankedLaneState> trial_states;
  std::vector<ExactRankedLaneFrontierEntry> frontier;
  std::vector<ExactRankedLaneFrontierEntry> next_frontier;
  std::vector<ExactSequenceState> states;
  std::vector<ExactSequenceState> next_states;
  std::size_t next_state_count{0U};
  std::vector<std::size_t> next_frontier_head_by_trial;
  std::vector<std::size_t> next_frontier_links;
  std::vector<double> conditional_totals;
  std::vector<ExactRankedWorkItem> work;
  std::vector<std::vector<std::size_t>> work_by_outcome;
  std::vector<double> totals;
  std::vector<double> transition_values;
  std::vector<double> readiness_values;
  std::vector<double> transition_normalizers;
  std::vector<double> trigger_weights;
};

struct ExactRankedLaneWorkspacePool {
  explicit ExactRankedLaneWorkspacePool(const std::size_t plan_count)
      : workspaces(plan_count) {}

  ExactRankedLaneWorkspace &get(
      const std::vector<ExactVariantPlan> &plans,
      const semantic::Index variant_index) {
    const auto position = static_cast<std::size_t>(variant_index);
    if (!workspaces[position]) {
      workspaces[position] =
          std::make_unique<ExactRankedLaneWorkspace>(plans[position]);
    }
    return *workspaces[position];
  }

  std::vector<std::unique_ptr<ExactRankedLaneWorkspace>> workspaces;
};

inline void reset_exact_ranked_tile(
    const ExactVariantPlan &plan,
    const ExactRankedLane *lanes,
    const std::size_t lane_count,
    const ExactCompiledTriggerState &compiled_trigger,
    const ExactSequenceState &initial_state,
    ExactRankedLaneWorkspace *workspace) {
  workspace->frontier.clear();
  workspace->next_frontier.clear();
  workspace->next_state_count = 0U;
  const bool fixed_weight = compiled_trigger.weight_terms.empty();
  if (!fixed_weight) {
    exact_compiled_trigger_state_weights_lanes(
        plan,
        lanes,
        lane_count,
        compiled_trigger,
        &workspace->trigger_weights);
  }
  const auto *shared_started =
      exact_compiled_trigger_shared_started(plan, compiled_trigger);
  for (std::size_t lane_index = 0U;
       lane_index < lane_count;
       ++lane_index) {
    auto &state = workspace->trial_states[lane_index];
    state.trigger = ExactTriggerState{
        fixed_weight ? 1.0
                     : workspace->trigger_weights[lane_index],
        shared_started};
    state.trigger_active = state.trigger.weight > 0.0;
    state.rank_pending = false;
    state.pending_target = semantic::kInvalidIndex;
    state.used_outcomes.assign(plan.compiled_outcomes.size(), 0U);
    if (!state.trigger_active) {
      continue;
    }
    const auto state_index = workspace->frontier.size();
    ranked_sequence_state_slot(
        &workspace->states, plan, state_index) = initial_state;
    workspace->frontier.push_back(ExactRankedLaneFrontierEntry{
        lane_index,
        1.0,
        static_cast<semantic::Index>(state_index)});
  }
}

inline void build_exact_ranked_work(
    const ExactVariantPlan &plan,
    const ExactRankedLane *lanes,
    const std::size_t lane_count,
    const std::size_t rank_index,
    ExactRankedLaneWorkspace *workspace) {
  workspace->work.clear();
  workspace->next_frontier.clear();
  workspace->next_frontier_links.clear();
  workspace->next_state_count = 0U;
  std::fill_n(
      workspace->next_frontier_head_by_trial.begin(),
      lane_count,
      std::numeric_limits<std::size_t>::max());
  std::fill_n(
      workspace->conditional_totals.begin(), lane_count, 0.0);
  for (auto &group : workspace->work_by_outcome) {
    group.clear();
  }
  for (const auto &entry : workspace->frontier) {
    if (entry.probability > 0.0) {
      workspace->conditional_totals[entry.trial_index] += entry.probability;
    }
  }
  for (std::size_t lane_index = 0; lane_index < lane_count; ++lane_index) {
    auto &state = workspace->trial_states[lane_index];
    state.rank_pending = false;
    if (!state.trigger_active) {
      continue;
    }
    if (rank_index >=
        static_cast<std::size_t>(
            lanes[lane_index].observation.rank_count)) {
      workspace->totals[lane_index] +=
          state.trigger.weight * workspace->conditional_totals[lane_index];
      state.trigger_active = false;
      continue;
    }
    const auto outcome_code = exact_trial_view_outcome_code(
        lanes[lane_index].observation, rank_index);
    const auto target = plan.outcome_index_by_code[
        static_cast<std::size_t>(outcome_code)];
    if (target == semantic::kInvalidIndex) {
      state.trigger_active = false;
      continue;
    }
    state.rank_pending = true;
    state.pending_target = target;
    state.pending_time = exact_trial_view_rt(
        lanes[lane_index].observation, rank_index);
  }
  for (const auto &entry : workspace->frontier) {
    const auto lane_index = entry.trial_index;
    if (!(entry.probability > 0.0)) {
      continue;
    }
    const auto &state = workspace->trial_states[lane_index];
    if (!state.trigger_active || !state.rank_pending) {
      continue;
    }
    const auto work_index = workspace->work.size();
    workspace->work.push_back(ExactRankedWorkItem{
        lane_index, entry.state_index, entry.probability});
    workspace->work_by_outcome[
        static_cast<std::size_t>(state.pending_target)].push_back(work_index);
  }
}

inline void evaluate_exact_ranked_work_group(
    const ExactVariantPlan &plan,
    const ExactRankedLane *lanes,
    const semantic::Index target,
    const std::vector<std::size_t> &group,
    ExactStepLaneWorkspace *step,
    ExactRankedLaneWorkspace *workspace) {
  const auto target_position = static_cast<std::size_t>(target);
  const auto &outcome = plan.compiled_outcomes[target_position];
  auto &step_inputs = step->step_inputs;
  for (std::size_t tile_start = 0U;
       tile_start < group.size();
       tile_start += kExactLaneTileSize) {
    const auto tile_count =
        std::min(kExactLaneTileSize, group.size() - tile_start);
    step_inputs.clear();
    for (std::size_t tile_lane = 0U;
         tile_lane < tile_count;
         ++tile_lane) {
      const auto &work = workspace->work[group[tile_start + tile_lane]];
      auto &trial_state = workspace->trial_states[work.trial_index];
      step_inputs.push_back(ExactStepLaneInput{
          &lanes[work.trial_index].params,
          trial_state.trigger.shared_started,
          &workspace->states[static_cast<std::size_t>(work.state_index)],
          trial_state.pending_time,
          &trial_state.used_outcomes});
    }
    evaluate_exact_ranked_step_lanes(
        plan,
        step_inputs.data(),
        tile_count,
        target,
        step,
        &workspace->transition_values,
        &workspace->readiness_values);

    for (std::size_t transition_index = 0U;
         transition_index < outcome.transitions.size();
         ++transition_index) {
      const auto &transition = outcome.transitions[transition_index];
      const auto readiness_span =
          outcome.transition_readiness_slots[transition_index];
      for (std::size_t tile_lane = 0U;
           tile_lane < tile_count;
           ++tile_lane) {
        const auto &work = workspace->work[group[tile_start + tile_lane]];
        auto &trial_state = workspace->trial_states[work.trial_index];
        const double transition_probability =
            workspace->transition_values[
                transition_index * tile_count + tile_lane];
        if (!(transition_probability > 0.0)) {
          continue;
        }
        const double probability =
            work.frontier_probability * transition_probability;
        if (!(probability > 0.0)) {
          continue;
        }
        const auto &current_state = workspace->states[
            static_cast<std::size_t>(work.state_index)];
        workspace->candidate_state = current_state;

        workspace->transition_normalizers.clear();
        workspace->transition_normalizers.reserve(
            static_cast<std::size_t>(readiness_span.size));
        for (semantic::Index i = 0; i < readiness_span.size; ++i) {
          const auto root_slot = outcome.readiness_root_slot_by_item[
              static_cast<std::size_t>(readiness_span.offset + i)];
          workspace->transition_normalizers.push_back(
              root_slot == semantic::kInvalidIndex
                  ? 0.0
                  : workspace->readiness_values[
                        static_cast<std::size_t>(root_slot) * tile_count +
                        tile_lane]);
        }
        advance_exact_sequence_state(
            workspace->candidate_state,
            transition,
            trial_state.pending_time,
            workspace->transition_normalizers);
        bool merged = false;
        auto existing_index =
            workspace->next_frontier_head_by_trial[work.trial_index];
        while (existing_index != std::numeric_limits<std::size_t>::max()) {
          auto &existing = workspace->next_frontier[existing_index];
          if (exact_sequence_states_equal(
                  workspace->next_states[
                      static_cast<std::size_t>(existing.state_index)],
                  workspace->candidate_state)) {
            existing.probability += probability;
            merged = true;
            break;
          }
          existing_index = workspace->next_frontier_links[existing_index];
        }
        if (!merged) {
          const auto candidate_index = workspace->next_state_count;
          auto &candidate_state = ranked_sequence_state_slot(
              &workspace->next_states, plan, candidate_index);
          std::swap(candidate_state, workspace->candidate_state);
          ++workspace->next_state_count;
          const auto new_index = workspace->next_frontier.size();
          workspace->next_frontier.push_back(
              ExactRankedLaneFrontierEntry{
                  work.trial_index,
                  probability,
                  static_cast<semantic::Index>(candidate_index)});
          workspace->next_frontier_links.push_back(
              workspace->next_frontier_head_by_trial[work.trial_index]);
          workspace->next_frontier_head_by_trial[work.trial_index] = new_index;
        }
      }
    }
  }
}

inline void finish_exact_ranked_rank(
    const std::size_t lane_count,
    ExactRankedLaneWorkspace *workspace) {
  for (std::size_t lane_index = 0U;
       lane_index < lane_count;
       ++lane_index) {
    auto &state = workspace->trial_states[lane_index];
    if (!state.trigger_active || !state.rank_pending) {
      continue;
    }
    if (workspace->next_frontier_head_by_trial[lane_index] ==
        std::numeric_limits<std::size_t>::max()) {
      state.trigger_active = false;
      continue;
    }
    state.used_outcomes[
        static_cast<std::size_t>(state.pending_target)] = 1U;
  }
  workspace->frontier.swap(workspace->next_frontier);
  workspace->states.swap(workspace->next_states);
  workspace->next_state_count = 0U;
}

inline void exact_ranked_loglik_lanes(
    const ExactVariantPlan &plan,
    const ExactRankedLane *lanes,
    const std::size_t lane_count,
    const double min_ll,
    ExactStepLaneWorkspace *step,
    ExactRankedLaneWorkspace *workspace,
    std::vector<double> *out) {
  out->assign(lane_count, min_ll);
  if (lane_count == 0U) {
    return;
  }
  for (std::size_t tile_start = 0U;
       tile_start < lane_count;
       tile_start += kExactRankedTrialTileSize) {
    const auto tile_count =
        std::min(kExactRankedTrialTileSize, lane_count - tile_start);
    const auto *tile_lanes = lanes + tile_start;
    workspace->ensure_trials(tile_count);
    std::size_t max_rank = 0U;
    for (std::size_t lane_index = 0U;
         lane_index < tile_count;
         ++lane_index) {
      max_rank = std::max(
          max_rank,
          static_cast<std::size_t>(
              tile_lanes[lane_index].observation.rank_count));
    }

    for (const auto &compiled_trigger : plan.trigger_state_table.states) {
      reset_exact_ranked_tile(
          plan,
          tile_lanes,
          tile_count,
          compiled_trigger,
          step->initial_state,
          workspace);

      for (std::size_t rank_index = 0U;
           rank_index < max_rank;
           ++rank_index) {
        build_exact_ranked_work(
            plan, tile_lanes, tile_count, rank_index, workspace);
        for (std::size_t target_position = 0U;
             target_position < workspace->work_by_outcome.size();
             ++target_position) {
          const auto &group = workspace->work_by_outcome[target_position];
          if (group.empty()) {
            continue;
          }
          evaluate_exact_ranked_work_group(
              plan,
              tile_lanes,
              static_cast<semantic::Index>(target_position),
              group,
              step,
              workspace);
        }
        finish_exact_ranked_rank(tile_count, workspace);
      }

      std::fill_n(
          workspace->conditional_totals.begin(), tile_count, 0.0);
      for (const auto &entry : workspace->frontier) {
        if (entry.probability > 0.0) {
          workspace->conditional_totals[entry.trial_index] +=
              entry.probability;
        }
      }
      for (std::size_t lane_index = 0U;
           lane_index < tile_count;
           ++lane_index) {
        const auto &state = workspace->trial_states[lane_index];
        if (state.trigger_active) {
          workspace->totals[lane_index] +=
              state.trigger.weight *
              workspace->conditional_totals[lane_index];
        }
      }
    }

    for (std::size_t lane_index = 0U;
         lane_index < tile_count;
         ++lane_index) {
      const double total = workspace->totals[lane_index];
      if (std::isfinite(total) && total > 0.0) {
        (*out)[tile_start + lane_index] = std::log(total);
      }
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
