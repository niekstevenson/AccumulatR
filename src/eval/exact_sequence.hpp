#pragma once

#include "eval_query.hpp"
#include "exact_competitor_union.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactTrialColumns {
  std::vector<const int *> labels;
  std::vector<const double *> times;
};

inline ExactTrialColumns read_exact_trial_columns(
    const PreparedTrialLayout &layout,
    SEXP dataSEXP) {
  ExactTrialColumns columns;
  const int max_rank = layout.max_rank;
  columns.labels.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  columns.times.assign(static_cast<std::size_t>(max_rank + 1), nullptr);
  for (int rank = 1; rank <= max_rank; ++rank) {
    columns.labels[static_cast<std::size_t>(rank)] =
        INTEGER(trusted_data_column(
            dataSEXP,
            layout.label_cols[static_cast<std::size_t>(rank)]));
    columns.times[static_cast<std::size_t>(rank)] =
        REAL(trusted_data_column(
            dataSEXP,
            layout.time_cols[static_cast<std::size_t>(rank)]));
  }
  return columns;
}

struct ExactTrialView {
  semantic::Index variant_index{semantic::kInvalidIndex};
  R_xlen_t row{0};
  int rank_count{0};
  const ExactTrialColumns *columns{nullptr};
};

inline semantic::Index exact_trial_view_outcome_code(const ExactTrialView &view,
                                                     const std::size_t rank_idx) {
  return view.columns->labels[rank_idx + 1U][view.row];
}

inline double exact_trial_view_rt(const ExactTrialView &view,
                                  const std::size_t rank_idx) {
  return view.columns->times[rank_idx + 1U][view.row];
}

inline ExactTrialView read_exact_observation_view(
    const PreparedDataView &table,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const PreparedTrialLayout &layout,
    const std::size_t trial_index,
    const ExactTrialColumns &columns) {
  const auto &span = layout.spans[trial_index];
  const auto row = static_cast<R_xlen_t>(span.start_row);
  const int component_code = table.component[row];
  ExactTrialView obs;
  obs.variant_index =
      variant_index_by_component_code[static_cast<std::size_t>(component_code)];
  obs.row = row;
  obs.columns = &columns;
  for (int rank = 1; rank <= layout.max_rank; ++rank) {
    const auto *rank_labels = columns.labels[static_cast<std::size_t>(rank)];
    const auto *rank_times = columns.times[static_cast<std::size_t>(rank)];
    if (integer_cell_is_na(rank_labels, row) &&
        Rcpp::NumericVector::is_na(rank_times[row])) {
      break;
    }
    ++obs.rank_count;
  }
  return obs;
}

inline void build_exact_plan_cache(
    const compile::CompiledModel &compiled,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_component_codes,
    const std::size_t n_outcome_codes,
    std::vector<semantic::Index> *variant_index_by_component_code,
    std::vector<ExactVariantPlan> *plans) {
  variant_index_by_component_code->assign(
      n_component_codes + 1U,
      semantic::kInvalidIndex);
  plans->clear();
  plans->reserve(compiled.variants.size());
  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    const auto plan_index = static_cast<semantic::Index>(plans->size());
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared component code for '" +
          variant.component_id + "'");
    }
    (*variant_index_by_component_code)[static_cast<std::size_t>(component_it->second)] =
        plan_index;
    plans->push_back(
      make_exact_variant_plan(
          runtime::lower_exact_variant(variant),
          outcome_code_by_label,
          n_outcome_codes));
  }
}

struct ExactProbabilityBatchLane {
  ParamView params;
  int first_param_row{0};
  double observed_time{NA_REAL};
  ExactStepWorkspace *workspace{nullptr};

  ExactProbabilityBatchLane(SEXP paramsSEXP,
                            const int *row_map,
                            const int row_offset,
                            const int first_param_row_,
                            const double observed_time_,
                            ExactStepWorkspace *workspace_)
      : params(paramsSEXP, row_map, row_offset),
        first_param_row(first_param_row_),
        observed_time(observed_time_),
        workspace(workspace_) {}
};

struct ExactProbabilityBatchWorkspace {
  std::vector<CompiledMathBatchLane> compiled_lanes;
  std::vector<ExactTriggerState> trigger_states;
  std::vector<ExactProbabilityBatchLane> sub_lanes;
  std::vector<std::size_t> sub_lane_indices;
  std::vector<ExactSourceChannels *> source_channels_by_lane;
  std::vector<semantic::Index> active_lanes;
  std::vector<double> observed_times;
  std::vector<double> values_a;
  std::vector<double> values_b;
  CompiledMathBatchSourceProductScratch source_scratch;
  CompiledMathBatchWorkspace compiled_workspace;
};

struct ExactRankedBatchLane {
  ParamView params;
  int first_param_row{0};
  ExactTrialView obs{};
  ExactStepWorkspace *workspace{nullptr};

  ExactRankedBatchLane(SEXP paramsSEXP,
                       const int *row_map,
                       const int row_offset,
                       const int first_param_row_,
                       const ExactTrialView &obs_,
                       ExactStepWorkspace *workspace_)
      : params(paramsSEXP, row_map, row_offset),
        first_param_row(first_param_row_),
        obs(obs_),
        workspace(workspace_) {}
};

struct ExactRankedStepBatchLane {
  ParamView params;
  int first_param_row{0};
  double observed_time{NA_REAL};
  ExactStepWorkspace *workspace{nullptr};
  const ExactSequenceState *sequence_state{nullptr};
  const ExactTriggerState *trigger_state{nullptr};
  const std::vector<std::uint8_t> *used_outcomes{nullptr};

  ExactRankedStepBatchLane(const ParamView &params_,
                           const int first_param_row_,
                           const double observed_time_,
                           ExactStepWorkspace *workspace_,
                           const ExactSequenceState *sequence_state_,
                           const ExactTriggerState *trigger_state_,
                           const std::vector<std::uint8_t> *used_outcomes_)
      : params(params_),
        first_param_row(first_param_row_),
        observed_time(observed_time_),
        workspace(workspace_),
        sequence_state(sequence_state_),
        trigger_state(trigger_state_),
        used_outcomes(used_outcomes_) {}
};

struct ExactRankedStepItem {
  std::size_t lane{0};
  std::size_t frontier_entry{0};
};

struct ExactRankedBatchWorkspace {
  std::vector<ExactRankedBatchLane> sub_lanes;
  std::vector<std::size_t> sub_lane_indices;
  std::vector<ExactTriggerState> trigger_states;
  std::vector<double> trigger_weights;
  std::vector<double> trigger_values;
  std::vector<double> lane_totals;
  std::vector<double> single_path_totals;

  std::vector<std::uint8_t> lane_alive;
  std::vector<semantic::Index> target_by_lane;
  std::vector<semantic::Index> unique_targets;
  std::vector<std::size_t> next_state_counts;

  std::vector<ExactRankedStepBatchLane> step_lanes;
  std::vector<ExactRankedStepItem> step_items;
  std::vector<CompiledMathBatchLane> compiled_lanes;
  std::vector<double> root_values;
  std::vector<double> transition_values;
  std::vector<double> ready_values;
  std::vector<double> candidate_ready_expr_normalizers;
  CompiledMathBatchWorkspace compiled_workspace;
};

inline bool prepare_exact_probability_leaf_batch_state(
    const ExactVariantPlan &plan,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    CompiledMathLaneBatchState *state) {
  (void)plan;
  if (batch_workspace == nullptr || state == nullptr ||
      trigger_states.size() < lanes.size()) {
    return false;
  }
  const auto lane_count = lanes.size();
  batch_workspace->source_channels_by_lane.assign(lane_count, nullptr);
  batch_workspace->active_lanes.resize(lane_count);
  batch_workspace->observed_times.resize(lane_count);
  batch_workspace->values_a.assign(lane_count, 0.0);
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr) {
      return false;
    }
    workspace->source_channels.reset(
        lanes[lane].params,
        lanes[lane].first_param_row,
        trigger_states[lane],
        workspace->initial_state,
        lanes[lane].observed_time);
    batch_workspace->source_channels_by_lane[lane] =
        &workspace->source_channels;
    batch_workspace->active_lanes[lane] =
        static_cast<semantic::Index>(lane);
    batch_workspace->observed_times[lane] = lanes[lane].observed_time;
  }
  *state = CompiledMathLaneBatchState{
      lane_count,
      0U,
      BatchTimeSlotView{},
      nullptr,
      &batch_workspace->source_channels_by_lane,
      nullptr,
      nullptr,
      nullptr};
  return true;
}

inline bool evaluate_exact_probability_top1_leaf_race_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  if (batch_workspace == nullptr || out_by_lane == nullptr) {
    return false;
  }
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return true;
  }
  CompiledMathLaneBatchState state;
  if (!prepare_exact_probability_leaf_batch_state(
          plan,
          lanes,
          trigger_states,
          batch_workspace,
          &state)) {
    return false;
  }
  std::fill(out_by_lane->begin(), out_by_lane->end(), 1.0);
  for (semantic::Index i = 0; i < op.source_span.size; ++i) {
    const auto source_id =
        programs.source_ids[
            static_cast<std::size_t>(op.source_span.offset + i)];
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= plan.source_kernels.size()) {
      return false;
    }
    const auto &kernel =
        plan.source_kernels[static_cast<std::size_t>(source_id)];
    if (kernel.leaf_index == semantic::kInvalidIndex ||
        static_cast<std::size_t>(kernel.leaf_index) >=
            plan.lowered.program.leaf_descriptors.size()) {
      return false;
    }
    const auto &leaf =
        plan.lowered.program.leaf_descriptors[
            static_cast<std::size_t>(kernel.leaf_index)];
    const auto channel_mask =
        source_id == op.target_source_id
            ? kLeafChannelPdf
            : kLeafChannelSurvival;
    if (!compiled_math_batch_leaf_values_from_times(
            leaf.dist_kind,
            kernel.leaf_index,
            leaf.onset_abs_value,
            state,
            batch_workspace->active_lanes.data(),
            lanes.size(),
            batch_workspace->observed_times.data(),
            channel_mask,
            batch_workspace->values_a.data(),
            &batch_workspace->source_scratch)) {
      return false;
    }
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const double product =
          (*out_by_lane)[lane] * batch_workspace->values_a[lane];
      (*out_by_lane)[lane] =
          std::isfinite(product) && product > 0.0 ? product : 0.0;
    }
  }
  return true;
}

inline bool evaluate_exact_probability_terminal_no_response_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  if (batch_workspace == nullptr || out_by_lane == nullptr) {
    return false;
  }
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return true;
  }
  CompiledMathLaneBatchState state;
  if (!prepare_exact_probability_leaf_batch_state(
          plan,
          lanes,
          trigger_states,
          batch_workspace,
          &state)) {
    return false;
  }
  std::fill(out_by_lane->begin(), out_by_lane->end(), 1.0);
  for (semantic::Index i = 0; i < op.source_span.size; ++i) {
    const auto source_id =
        programs.source_ids[
            static_cast<std::size_t>(op.source_span.offset + i)];
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= plan.source_kernels.size()) {
      return false;
    }
    const auto &kernel =
        plan.source_kernels[static_cast<std::size_t>(source_id)];
    if (kernel.leaf_index == semantic::kInvalidIndex) {
      return false;
    }
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const auto *source_channels = state.source_channels(
          static_cast<semantic::Index>(lane));
      if (source_channels == nullptr) {
        return false;
      }
      const double q =
          clamp_probability(
              source_channels->source_product_leaf_input(
                  kernel.leaf_index).q);
      const double product = (*out_by_lane)[lane] * q;
      (*out_by_lane)[lane] =
          std::isfinite(product) && product > 0.0 ? product : 0.0;
    }
  }
  for (auto &value : *out_by_lane) {
    value = clamp_probability(value);
  }
  return true;
}

inline void evaluate_exact_probability_trigger_op_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane);

inline void evaluate_exact_probability_program_op_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane);

inline void exact_probability_batch_lanes_at_time(
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const double observed_time,
    std::vector<ExactProbabilityBatchLane> *out) {
  out->clear();
  out->reserve(lanes.size());
  for (const auto &lane : lanes) {
    out->push_back(lane);
    out->back().observed_time = observed_time;
  }
}

inline void exact_probability_accumulate_positive_values(
    const std::vector<double> &values,
    const double weight,
    std::vector<double> *out_by_lane) {
  for (std::size_t lane = 0; lane < values.size(); ++lane) {
    const double value = values[lane];
    if (std::isfinite(value) && value > 0.0) {
      (*out_by_lane)[lane] += weight * value;
    }
  }
}

inline void prepare_exact_ranked_step_compiled_lanes(
    const std::vector<ExactRankedStepBatchLane> &lanes,
    ExactRankedBatchWorkspace *batch_workspace) {
  auto &compiled_lanes = batch_workspace->compiled_lanes;
  compiled_lanes.clear();
  compiled_lanes.reserve(lanes.size());
  for (const auto &lane : lanes) {
    if (lane.workspace == nullptr ||
        lane.sequence_state == nullptr ||
        lane.trigger_state == nullptr) {
      compiled_lanes.push_back(CompiledMathBatchLane{});
      continue;
    }
    lane.workspace->reset(
        lane.params,
        lane.first_param_row,
        *lane.trigger_state,
        *lane.sequence_state,
        lane.observed_time);
    lane.workspace->target_workspace.compiled_math.used_outcomes =
        lane.used_outcomes;
    compiled_lanes.push_back(CompiledMathBatchLane{
        &lane.workspace->target_workspace.compiled_math,
        &lane.workspace->target_evaluator,
        nullptr,
        &lane.workspace->target_workspace});
  }
}

inline void evaluate_exact_ranked_step_transition_probabilities_batch(
    const ExactVariantPlan &plan,
    const semantic::Index outcome_index,
    const std::vector<ExactRankedStepBatchLane> &lanes,
    ExactRankedBatchWorkspace *batch_workspace,
    std::vector<double> *transition_values) {
  transition_values->clear();
  if (lanes.empty() || outcome_index == semantic::kInvalidIndex) {
    return;
  }
  const auto outcome_pos = static_cast<std::size_t>(outcome_index);
  const auto &successors =
      plan.runtime.outcomes[outcome_pos].successor_distribution;
  const auto transition_count = successors.transitions.size();
  transition_values->assign(transition_count * lanes.size(), 0.0);
  if (transition_count == 0U) {
    return;
  }

  prepare_exact_ranked_step_compiled_lanes(lanes, batch_workspace);
  for (std::size_t transition_idx = 0; transition_idx < transition_count;
       ++transition_idx) {
    const auto root_id =
        successors.transitions[transition_idx].probability_root_id;
    if (root_id == semantic::kInvalidIndex) {
      continue;
    }
    evaluate_compiled_math_root_batch(
        plan.compiled_math,
        root_id,
        batch_workspace->compiled_lanes,
        &batch_workspace->compiled_workspace,
        &batch_workspace->root_values);
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const double value = batch_workspace->root_values[lane];
      (*transition_values)[transition_idx * lanes.size() + lane] =
          std::isfinite(value) && value > 0.0 ? value : 0.0;
    }
  }
}

inline void exact_sequence_ready_expr_normalizers_for_prepared_ranked_steps(
    const ExactVariantPlan &plan,
    const ExactRuntimeScenarioTransitionPlan &transition,
    const std::size_t lane_count,
    ExactRankedBatchWorkspace *batch_workspace,
    std::vector<double> *normalizers) {
  const auto ready_count =
      static_cast<std::size_t>(transition.ready_expr_span.size);
  normalizers->assign(ready_count * lane_count, 0.0);
  if (ready_count == 0U || lane_count == 0U) {
    return;
  }

  for (semantic::Index i = 0; i < transition.ready_expr_span.size; ++i) {
    const auto expr_id =
        plan.scenario_expr_ids[
            static_cast<std::size_t>(
                transition.ready_expr_span.offset + i)];
    if (expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(expr_id) >=
            plan.sequence_expr_cdf_roots.size()) {
      continue;
    }
    const auto root_id =
        plan.sequence_expr_cdf_roots[static_cast<std::size_t>(expr_id)];
    if (root_id == semantic::kInvalidIndex) {
      continue;
    }
    evaluate_compiled_math_root_batch(
        plan.compiled_math,
        root_id,
        batch_workspace->compiled_lanes,
        &batch_workspace->compiled_workspace,
        &batch_workspace->root_values);
    const auto ready_pos = static_cast<std::size_t>(i);
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const double value = batch_workspace->root_values[lane];
      (*normalizers)[ready_pos * lane_count + lane] =
          std::isfinite(value) ? clamp_probability(value) : 0.0;
    }
  }
}

inline void evaluate_exact_step_distribution_total_probability_batch(
    const ExactVariantPlan &plan,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const std::vector<ExactTriggerState> &trigger_states,
    const semantic::Index outcome_index,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty() || outcome_index == semantic::kInvalidIndex) {
    return;
  }
  const auto target_pos = static_cast<std::size_t>(outcome_index);
  const auto &runtime_outcome = plan.runtime.outcomes[target_pos];
  const auto &successors = runtime_outcome.successor_distribution;
  const auto root_id = successors.total_probability_root_id;
  if (root_id == semantic::kInvalidIndex) {
    return;
  }

  auto &compiled_lanes = batch_workspace->compiled_lanes;
  compiled_lanes.clear();
  compiled_lanes.reserve(lanes.size());
  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr) {
      compiled_lanes.push_back(CompiledMathBatchLane{});
      continue;
    }
    workspace->reset(
        lanes[lane].params,
        lanes[lane].first_param_row,
        trigger_states[lane],
        workspace->initial_state,
        lanes[lane].observed_time);
    workspace->target_workspace.compiled_math.used_outcomes = nullptr;
    compiled_lanes.push_back(CompiledMathBatchLane{
        &workspace->target_workspace.compiled_math,
        &workspace->target_evaluator,
        nullptr,
        &workspace->target_workspace});
  }

  evaluate_compiled_math_root_batch(
      plan.compiled_math,
      root_id,
      compiled_lanes,
      &batch_workspace->compiled_workspace,
      out_by_lane);
  for (auto &value : *out_by_lane) {
    if (!std::isfinite(value) || !(value > 0.0)) {
      value = 0.0;
    }
  }
}

inline void evaluate_exact_probability_trigger_children_sum_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &root,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  for (semantic::Index i = 0; i < root.children.size; ++i) {
    const auto child =
        program.child_ops[
            static_cast<std::size_t>(root.children.offset + i)];
    if (child == semantic::kInvalidIndex) {
      continue;
    }
    std::vector<double> child_values;
    evaluate_exact_probability_trigger_op_batch(
        plan,
        program,
        programs,
        lanes,
        program.ops[static_cast<std::size_t>(child)],
        trigger_states,
        batch_workspace,
        &child_values);
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      (*out_by_lane)[lane] += child_values[lane];
    }
  }
}

inline void evaluate_exact_probability_weighted_trigger_sum_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &root,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return;
  }

  auto evaluate_for_states =
      [&](const std::vector<ExactProbabilityBatchLane> &state_lanes,
          const std::vector<ExactTriggerState> &states,
          std::vector<double> *state_values) {
        if (program.root_child != semantic::kInvalidIndex) {
          evaluate_exact_probability_trigger_op_batch(
              plan,
              program,
              programs,
              state_lanes,
              program.ops[static_cast<std::size_t>(program.root_child)],
              states,
              batch_workspace,
              state_values);
          return;
        }
        evaluate_exact_probability_trigger_children_sum_batch(
            plan,
            program,
            programs,
            state_lanes,
            root,
            states,
            batch_workspace,
            state_values);
      };

  if (!program.requires_trigger_enumeration) {
    auto &states = batch_workspace->trigger_states;
    states.clear();
    states.reserve(lanes.size());
    for (const auto &lane : lanes) {
      states.push_back(lane.workspace->default_trigger_state);
    }
    evaluate_for_states(lanes, states, out_by_lane);
    for (auto &value : *out_by_lane) {
      value = root.value_kind == ExactProbabilityValueKind::Probability
                  ? clamp_probability(value)
                  : (std::isfinite(value) && value > 0.0 ? value : 0.0);
    }
    return;
  }

  auto &sub_lanes = batch_workspace->sub_lanes;
  auto &sub_indices = batch_workspace->sub_lane_indices;
  auto &states = batch_workspace->trigger_states;
  std::vector<double> state_weights;
  for (const auto &compiled_state : plan.trigger_state_table.states) {
    sub_lanes.clear();
    sub_indices.clear();
    states.clear();
    state_weights.clear();
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const auto state =
          exact_compiled_trigger_state_view(
              plan,
              lanes[lane].params,
              lanes[lane].first_param_row,
              compiled_state);
      if (!(state.weight > 0.0)) {
        continue;
      }
      sub_indices.push_back(lane);
      sub_lanes.push_back(lanes[lane]);
      states.push_back(state);
      state_weights.push_back(state.weight);
    }
    if (sub_lanes.empty()) {
      continue;
    }
    std::vector<double> state_values;
    evaluate_for_states(sub_lanes, states, &state_values);
    for (std::size_t i = 0; i < sub_lanes.size(); ++i) {
      (*out_by_lane)[sub_indices[i]] += state_weights[i] * state_values[i];
    }
  }

  for (auto &value : *out_by_lane) {
    value = root.value_kind == ExactProbabilityValueKind::Probability
                ? clamp_probability(value)
                : (std::isfinite(value) && value > 0.0 ? value : 0.0);
  }
}

inline void evaluate_exact_probability_generic_transition_probability_batch(
    const ExactVariantPlan &plan,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const std::vector<ExactTriggerState> &trigger_states,
    const semantic::Index outcome_index,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return;
  }
  const auto &rule = quadrature::canonical_tail_batch().nodes;
  std::vector<ExactProbabilityBatchLane> time_lanes;
  std::vector<double> sample_values;
  time_lanes.reserve(lanes.size());
  for (std::size_t q = 0; q < quadrature::kDefaultTailOrder; ++q) {
    exact_probability_batch_lanes_at_time(lanes, rule.nodes[q], &time_lanes);
    evaluate_exact_step_distribution_total_probability_batch(
        plan,
        time_lanes,
        trigger_states,
        outcome_index,
        batch_workspace,
        &sample_values);
    exact_probability_accumulate_positive_values(
        sample_values,
        rule.weights[q],
        out_by_lane);
  }
  for (auto &value : *out_by_lane) {
    value = clamp_probability(value);
  }
}

inline void evaluate_exact_probability_trigger_integral_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty() || op.children.size == 0) {
    return;
  }
  const auto child =
      program.child_ops[static_cast<std::size_t>(op.children.offset)];
  const auto &child_op = program.ops[static_cast<std::size_t>(child)];
  const auto &rule = quadrature::canonical_tail_batch().nodes;
  std::vector<ExactProbabilityBatchLane> time_lanes;
  std::vector<double> sample_values;
  time_lanes.reserve(lanes.size());
  for (std::size_t q = 0; q < quadrature::kDefaultTailOrder; ++q) {
    exact_probability_batch_lanes_at_time(lanes, rule.nodes[q], &time_lanes);
    evaluate_exact_probability_trigger_op_batch(
        plan,
        program,
        programs,
        time_lanes,
        child_op,
        trigger_states,
        batch_workspace,
        &sample_values);
    exact_probability_accumulate_positive_values(
        sample_values,
        rule.weights[q],
        out_by_lane);
  }
  if (op.value_kind == ExactProbabilityValueKind::Probability) {
    for (auto &value : *out_by_lane) {
      value = clamp_probability(value);
    }
  }
}

inline void evaluate_exact_probability_trigger_op_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const ExactProbabilityProgramSet &programs,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return;
  }
  switch (op.kind) {
  case ExactProbabilityOpKind::Constant:
    std::fill(out_by_lane->begin(), out_by_lane->end(), op.constant);
    return;

  case ExactProbabilityOpKind::Top1LeafRaceDensity:
    if (!evaluate_exact_probability_top1_leaf_race_batch(
            plan,
            programs,
            lanes,
            op,
            trigger_states,
            batch_workspace,
            out_by_lane)) {
      throw std::runtime_error(
          "top-1 leaf race probability op has no batch implementation");
    }
    return;

  case ExactProbabilityOpKind::TerminalNoResponseProbability:
    if (!evaluate_exact_probability_terminal_no_response_batch(
            plan,
            programs,
            lanes,
            op,
            trigger_states,
            batch_workspace,
            out_by_lane)) {
      throw std::runtime_error(
          "terminal no-response probability op has no batch implementation");
    }
    return;

  case ExactProbabilityOpKind::GenericTransitionDensity:
    evaluate_exact_step_distribution_total_probability_batch(
        plan,
        lanes,
        trigger_states,
        op.outcome_index,
        batch_workspace,
        out_by_lane);
    return;

  case ExactProbabilityOpKind::GenericTransitionProbability:
    evaluate_exact_probability_generic_transition_probability_batch(
        plan,
        lanes,
        trigger_states,
        op.outcome_index,
        batch_workspace,
        out_by_lane);
    return;

  case ExactProbabilityOpKind::Integral:
    evaluate_exact_probability_trigger_integral_batch(
        plan,
        program,
        programs,
        lanes,
        op,
        trigger_states,
        batch_workspace,
        out_by_lane);
    return;

  case ExactProbabilityOpKind::Log: {
    const auto child =
        program.child_ops[static_cast<std::size_t>(op.children.offset)];
    std::vector<double> child_values;
    evaluate_exact_probability_trigger_op_batch(
        plan,
        program,
        programs,
        lanes,
        program.ops[static_cast<std::size_t>(child)],
        trigger_states,
        batch_workspace,
        &child_values);
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const double probability = child_values[lane];
      (*out_by_lane)[lane] =
          std::isfinite(probability) && probability > 0.0
              ? std::log(probability)
              : R_NegInf;
    }
    return;
  }

  case ExactProbabilityOpKind::WeightedTriggerSum:
  case ExactProbabilityOpKind::RankedTransitionSequence:
    return;
  }
}

inline void evaluate_exact_probability_program_integral_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty() || op.children.size == 0) {
    return;
  }
  const auto child =
      program.child_ops[static_cast<std::size_t>(op.children.offset)];
  const auto &child_op = program.ops[static_cast<std::size_t>(child)];
  const auto &rule = quadrature::canonical_tail_batch().nodes;
  std::vector<ExactProbabilityBatchLane> time_lanes;
  std::vector<double> sample_values;
  time_lanes.reserve(lanes.size());
  for (std::size_t q = 0; q < quadrature::kDefaultTailOrder; ++q) {
    exact_probability_batch_lanes_at_time(lanes, rule.nodes[q], &time_lanes);
    evaluate_exact_probability_program_op_batch(
        plan,
        program,
        time_lanes,
        child_op,
        batch_workspace,
        &sample_values);
    exact_probability_accumulate_positive_values(
        sample_values,
        rule.weights[q],
        out_by_lane);
  }
  if (op.value_kind == ExactProbabilityValueKind::Probability) {
    for (auto &value : *out_by_lane) {
      value = clamp_probability(value);
    }
  }
}

inline void evaluate_exact_probability_program_op_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    const ExactProbabilityOp &op,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (lanes.empty()) {
    return;
  }
  switch (op.kind) {
  case ExactProbabilityOpKind::WeightedTriggerSum:
    evaluate_exact_probability_weighted_trigger_sum_batch(
        plan,
        program,
        plan.probability_programs,
        lanes,
        op,
        batch_workspace,
        out_by_lane);
    return;

  case ExactProbabilityOpKind::Integral:
    evaluate_exact_probability_program_integral_batch(
        plan,
        program,
        lanes,
        op,
        batch_workspace,
        out_by_lane);
    return;

  case ExactProbabilityOpKind::Log: {
    const auto child =
        program.child_ops[static_cast<std::size_t>(op.children.offset)];
    std::vector<double> child_values;
    evaluate_exact_probability_program_op_batch(
        plan,
        program,
        lanes,
        program.ops[static_cast<std::size_t>(child)],
        batch_workspace,
        &child_values);
    for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
      const double probability = child_values[lane];
      (*out_by_lane)[lane] =
          std::isfinite(probability) && probability > 0.0
              ? std::log(probability)
              : R_NegInf;
    }
    return;
  }

  case ExactProbabilityOpKind::Top1LeafRaceDensity:
  case ExactProbabilityOpKind::TerminalNoResponseProbability:
  case ExactProbabilityOpKind::Constant:
  case ExactProbabilityOpKind::GenericTransitionDensity:
  case ExactProbabilityOpKind::GenericTransitionProbability: {
    auto &states = batch_workspace->trigger_states;
    states.clear();
    states.reserve(lanes.size());
    for (const auto &lane : lanes) {
      states.push_back(lane.workspace->default_trigger_state);
    }
    evaluate_exact_probability_trigger_op_batch(
        plan,
        program,
        plan.probability_programs,
        lanes,
        op,
        states,
        batch_workspace,
        out_by_lane);
    return;
  }

  case ExactProbabilityOpKind::RankedTransitionSequence:
    return;
  }
}

inline void evaluate_exact_probability_program_batch(
    const ExactVariantPlan &plan,
    const ExactProbabilityProgram &program,
    const std::vector<ExactProbabilityBatchLane> &lanes,
    ExactProbabilityBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  out_by_lane->assign(lanes.size(), 0.0);
  if (program.empty() || lanes.empty()) {
    return;
  }
  const auto &root = program.ops[static_cast<std::size_t>(program.root)];
  evaluate_exact_probability_program_op_batch(
      plan,
      program,
      lanes,
      root,
      batch_workspace,
      out_by_lane);
}

inline void advance_exact_sequence_state(
    ExactSequenceState *state,
    const ExactRuntimeScenarioTransitionPlan &transition,
    const ExactVariantPlan &plan,
    const double observed_time,
    const std::vector<double> *ready_expr_normalizers = nullptr) {
  if (state == nullptr) {
    return;
  }
  state->lower_bound = observed_time;
  if (transition.active_source_id != semantic::kInvalidIndex &&
      static_cast<std::size_t>(transition.active_source_id) <
          state->exact_times.size()) {
    state->exact_times[static_cast<std::size_t>(transition.active_source_id)] =
        observed_time;
  }
  for (semantic::Index i = 0; i < transition.before_source_span.size; ++i) {
    const auto source_id =
        plan.scenario_source_ids[
            static_cast<std::size_t>(
                transition.before_source_span.offset + i)];
    if (source_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(source_id) >= state->upper_bounds.size()) {
      continue;
    }
    auto &upper = state->upper_bounds[static_cast<std::size_t>(source_id)];
    upper = std::isfinite(upper) ? std::min(upper, observed_time)
                                 : observed_time;
  }
  for (semantic::Index i = 0; i < transition.ready_expr_span.size; ++i) {
    const auto expr_id =
        plan.scenario_expr_ids[
            static_cast<std::size_t>(
                transition.ready_expr_span.offset + i)];
    if (expr_id == semantic::kInvalidIndex ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_bounds.size() ||
        static_cast<std::size_t>(expr_id) >=
            state->expr_upper_normalizers.size()) {
      continue;
    }
    const std::size_t ready_pos = static_cast<std::size_t>(i);
    if (ready_expr_normalizers == nullptr ||
        ready_pos >= ready_expr_normalizers->size()) {
      continue;
    }
    const double normalizer = (*ready_expr_normalizers)[ready_pos];
    if (!(normalizer > 0.0) || !std::isfinite(normalizer)) {
      continue;
    }
    auto &upper =
        state->expr_upper_bounds[static_cast<std::size_t>(expr_id)];
    if (std::isfinite(upper) && upper <= observed_time) {
      continue;
    }
    upper = observed_time;
    state->expr_upper_normalizers[static_cast<std::size_t>(expr_id)] =
        normalizer;
  }
}

inline bool exact_sequence_states_equal(const ExactSequenceState &lhs,
                                        const ExactSequenceState &rhs) {
  if (lhs.lower_bound != rhs.lower_bound ||
      lhs.exact_times.size() != rhs.exact_times.size() ||
      lhs.upper_bounds.size() != rhs.upper_bounds.size()) {
    return false;
  }
  if (lhs.expr_upper_bounds.size() != rhs.expr_upper_bounds.size() ||
      lhs.expr_upper_normalizers.size() !=
          rhs.expr_upper_normalizers.size()) {
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

inline void append_ranked_frontier_entry(
    std::vector<ExactRankedFrontierEntry> *frontier,
    const std::vector<ExactSequenceState> &states,
    const semantic::Index state_index,
    const double probability) {
  if (!(probability > 0.0)) {
    return;
  }
  const auto state_pos = static_cast<std::size_t>(state_index);
  for (auto &existing : *frontier) {
    if (exact_sequence_states_equal(
            states[static_cast<std::size_t>(existing.state_index)],
            states[state_pos])) {
      existing.probability += probability;
      return;
    }
  }
  frontier->push_back(ExactRankedFrontierEntry{probability, state_index});
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

inline bool exact_ranked_add_unique_target(
    std::vector<semantic::Index> *targets,
    const semantic::Index target) {
  for (const auto existing : *targets) {
    if (existing == target) {
      return false;
    }
  }
  targets->push_back(target);
  return true;
}

inline bool exact_ranked_trigger_probability_single_path_batch(
    const ExactVariantPlan &plan,
    const std::vector<ExactRankedBatchLane> &lanes,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactRankedBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  const auto lane_count = lanes.size();
  if (batch_workspace == nullptr || out_by_lane == nullptr) {
    return false;
  }
  out_by_lane->assign(lane_count, 0.0);
  if (lane_count == 0U) {
    return true;
  }

  auto &alive = batch_workspace->lane_alive;
  auto &target_by_lane = batch_workspace->target_by_lane;
  auto &unique_targets = batch_workspace->unique_targets;
  auto &totals = batch_workspace->single_path_totals;
  alive.assign(lane_count, 1U);
  target_by_lane.assign(lane_count, semantic::kInvalidIndex);
  totals.assign(lane_count, 1.0);

  int max_rank_count = 0;
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr || lane >= trigger_states.size()) {
      alive[lane] = 0U;
      totals[lane] = 0.0;
      continue;
    }
    workspace->ranked_used_outcomes.assign(plan.outcomes.size(), 0U);
    workspace->ranked_states.clear();
    ranked_sequence_state_slot(&workspace->ranked_states, plan, 0) =
        workspace->initial_state;
    max_rank_count = std::max(max_rank_count, lanes[lane].obs.rank_count);
  }

  for (std::size_t rank_idx = 0;
       rank_idx < static_cast<std::size_t>(max_rank_count);
       ++rank_idx) {
    bool any_active = false;
    unique_targets.clear();
    std::fill(
        target_by_lane.begin(),
        target_by_lane.end(),
        semantic::kInvalidIndex);

    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      if (alive[lane] == 0U ||
          rank_idx >= static_cast<std::size_t>(lanes[lane].obs.rank_count)) {
        continue;
      }
      const auto outcome_code =
          exact_trial_view_outcome_code(lanes[lane].obs, rank_idx);
      const auto target_outcome_index =
          plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
      if (target_outcome_index == semantic::kInvalidIndex) {
        alive[lane] = 0U;
        totals[lane] = 0.0;
        continue;
      }
      auto *workspace = lanes[lane].workspace;
      const auto target_pos = static_cast<std::size_t>(target_outcome_index);
      if (target_pos >= workspace->ranked_used_outcomes.size() ||
          workspace->ranked_used_outcomes[target_pos] != 0U) {
        alive[lane] = 0U;
        totals[lane] = 0.0;
        continue;
      }
      const auto &successors =
          plan.runtime.outcomes[target_pos].successor_distribution;
      if (successors.transitions.size() != 1U) {
        return false;
      }
      target_by_lane[lane] = target_outcome_index;
      exact_ranked_add_unique_target(&unique_targets, target_outcome_index);
      any_active = true;
    }
    if (!any_active) {
      continue;
    }

    for (const auto target_outcome_index : unique_targets) {
      auto &step_lanes = batch_workspace->step_lanes;
      auto &step_items = batch_workspace->step_items;
      step_lanes.clear();
      step_items.clear();
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        if (target_by_lane[lane] != target_outcome_index) {
          continue;
        }
        auto *workspace = lanes[lane].workspace;
        if (workspace->ranked_states.empty()) {
          return false;
        }
        const double observed_time =
            exact_trial_view_rt(lanes[lane].obs, rank_idx);
        step_items.push_back(ExactRankedStepItem{lane, 0U});
        step_lanes.emplace_back(
            lanes[lane].params,
            lanes[lane].first_param_row,
            observed_time,
            workspace,
            &workspace->ranked_states[0],
            &trigger_states[lane],
            &workspace->ranked_used_outcomes);
      }
      if (step_lanes.empty()) {
        continue;
      }

      evaluate_exact_ranked_step_transition_probabilities_batch(
          plan,
          target_outcome_index,
          step_lanes,
          batch_workspace,
          &batch_workspace->transition_values);
      const auto &successors =
          plan.runtime.outcomes[static_cast<std::size_t>(target_outcome_index)]
              .successor_distribution;
      const auto &transition = successors.transitions[0];
      exact_sequence_ready_expr_normalizers_for_prepared_ranked_steps(
          plan,
          transition,
          step_lanes.size(),
          batch_workspace,
          &batch_workspace->ready_values);
      const auto ready_count =
          static_cast<std::size_t>(transition.ready_expr_span.size);
      for (std::size_t step_lane = 0; step_lane < step_lanes.size();
           ++step_lane) {
        const auto lane = step_items[step_lane].lane;
        const double transition_probability =
            batch_workspace->transition_values[step_lane];
        if (!(transition_probability > 0.0)) {
          alive[lane] = 0U;
          totals[lane] = 0.0;
          continue;
        }
        totals[lane] *= transition_probability;
        auto &candidate_normalizers =
            batch_workspace->candidate_ready_expr_normalizers;
        if (ready_count != 0U) {
          candidate_normalizers.resize(ready_count);
          for (std::size_t ready_pos = 0; ready_pos < ready_count;
               ++ready_pos) {
            candidate_normalizers[ready_pos] =
                batch_workspace->ready_values[
                    ready_pos * step_lanes.size() + step_lane];
          }
        }
        advance_exact_sequence_state(
            lanes[lane].workspace->ranked_states.data(),
            transition,
            plan,
            step_lanes[step_lane].observed_time,
            ready_count == 0U ? nullptr : &candidate_normalizers);
      }
    }

    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto target_outcome_index = target_by_lane[lane];
      if (target_outcome_index == semantic::kInvalidIndex) {
        continue;
      }
      if (alive[lane] == 0U) {
        continue;
      }
      lanes[lane].workspace->ranked_used_outcomes[
          static_cast<std::size_t>(target_outcome_index)] = 1U;
    }
  }

  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    (*out_by_lane)[lane] =
        alive[lane] != 0U && std::isfinite(totals[lane]) && totals[lane] > 0.0
            ? totals[lane]
            : 0.0;
  }
  return true;
}

inline void exact_ranked_trigger_probability_batch(
    const ExactVariantPlan &plan,
    const std::vector<ExactRankedBatchLane> &lanes,
    const std::vector<ExactTriggerState> &trigger_states,
    ExactRankedBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  const auto lane_count = lanes.size();
  out_by_lane->assign(lane_count, 0.0);
  if (lane_count == 0U) {
    return;
  }
  if (exact_ranked_trigger_probability_single_path_batch(
          plan,
          lanes,
          trigger_states,
          batch_workspace,
          out_by_lane)) {
    return;
  }

  auto &alive = batch_workspace->lane_alive;
  auto &target_by_lane = batch_workspace->target_by_lane;
  auto &unique_targets = batch_workspace->unique_targets;
  auto &next_state_counts = batch_workspace->next_state_counts;
  alive.assign(lane_count, 1U);
  target_by_lane.assign(lane_count, semantic::kInvalidIndex);
  next_state_counts.assign(lane_count, 0U);

  int max_rank_count = 0;
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    auto *workspace = lanes[lane].workspace;
    if (workspace == nullptr || lane >= trigger_states.size()) {
      alive[lane] = 0U;
      continue;
    }
    workspace->ranked_used_outcomes.assign(plan.outcomes.size(), 0U);
    workspace->ranked_frontier.clear();
    workspace->ranked_next_frontier.clear();
    workspace->ranked_states.clear();
    workspace->ranked_next_states.clear();
    ranked_sequence_state_slot(&workspace->ranked_states, plan, 0) =
        workspace->initial_state;
    workspace->ranked_frontier.push_back(ExactRankedFrontierEntry{1.0, 0});
    max_rank_count = std::max(max_rank_count, lanes[lane].obs.rank_count);
  }

  for (std::size_t rank_idx = 0;
       rank_idx < static_cast<std::size_t>(max_rank_count);
       ++rank_idx) {
    bool any_active = false;
    unique_targets.clear();
    std::fill(
        target_by_lane.begin(),
        target_by_lane.end(),
        semantic::kInvalidIndex);

    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      if (alive[lane] == 0U ||
          rank_idx >= static_cast<std::size_t>(lanes[lane].obs.rank_count)) {
        continue;
      }
      const auto outcome_code =
          exact_trial_view_outcome_code(lanes[lane].obs, rank_idx);
      const auto target_outcome_index =
          plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
      if (target_outcome_index == semantic::kInvalidIndex) {
        alive[lane] = 0U;
        continue;
      }
      auto *workspace = lanes[lane].workspace;
      const auto target_pos = static_cast<std::size_t>(target_outcome_index);
      if (target_pos >= workspace->ranked_used_outcomes.size() ||
          workspace->ranked_used_outcomes[target_pos] != 0U) {
        alive[lane] = 0U;
        workspace->ranked_frontier.clear();
        continue;
      }
      target_by_lane[lane] = target_outcome_index;
      exact_ranked_add_unique_target(&unique_targets, target_outcome_index);
      workspace->ranked_next_frontier.clear();
      next_state_counts[lane] = 0U;
      any_active = true;
    }

    if (!any_active) {
      continue;
    }

    for (const auto target_outcome_index : unique_targets) {
      auto &step_lanes = batch_workspace->step_lanes;
      auto &step_items = batch_workspace->step_items;
      step_lanes.clear();
      step_items.clear();
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        if (target_by_lane[lane] != target_outcome_index) {
          continue;
        }
        auto *workspace = lanes[lane].workspace;
        const double observed_time =
            exact_trial_view_rt(lanes[lane].obs, rank_idx);
        for (std::size_t entry_idx = 0;
             entry_idx < workspace->ranked_frontier.size();
             ++entry_idx) {
          const auto &entry = workspace->ranked_frontier[entry_idx];
          if (!(entry.probability > 0.0)) {
            continue;
          }
          const auto state_pos = static_cast<std::size_t>(entry.state_index);
          if (state_pos >= workspace->ranked_states.size()) {
            continue;
          }
          step_items.push_back(ExactRankedStepItem{lane, entry_idx});
          step_lanes.emplace_back(
              lanes[lane].params,
              lanes[lane].first_param_row,
              observed_time,
              workspace,
              &workspace->ranked_states[state_pos],
              &trigger_states[lane],
              &workspace->ranked_used_outcomes);
        }
      }
      if (step_lanes.empty()) {
        continue;
      }

      evaluate_exact_ranked_step_transition_probabilities_batch(
          plan,
          target_outcome_index,
          step_lanes,
          batch_workspace,
          &batch_workspace->transition_values);
      const auto &successors =
          plan.runtime.outcomes[static_cast<std::size_t>(target_outcome_index)]
              .successor_distribution;
      const auto transition_count = successors.transitions.size();
      for (std::size_t transition_idx = 0;
           transition_idx < transition_count;
           ++transition_idx) {
        const auto &transition = successors.transitions[transition_idx];
        exact_sequence_ready_expr_normalizers_for_prepared_ranked_steps(
            plan,
            transition,
            step_lanes.size(),
            batch_workspace,
            &batch_workspace->ready_values);
        const auto ready_count =
            static_cast<std::size_t>(transition.ready_expr_span.size);
        for (std::size_t step_lane = 0; step_lane < step_lanes.size();
             ++step_lane) {
          const double transition_probability =
              batch_workspace->transition_values[
                  transition_idx * step_lanes.size() + step_lane];
          if (!(transition_probability > 0.0)) {
            continue;
          }
          const auto &item = step_items[step_lane];
          auto *workspace = lanes[item.lane].workspace;
          const auto &entry =
              workspace->ranked_frontier[item.frontier_entry];
          const auto candidate_index = next_state_counts[item.lane]++;
          auto &candidate_state =
              ranked_sequence_state_slot(
                  &workspace->ranked_next_states,
                  plan,
                  candidate_index);
          candidate_state =
              workspace->ranked_states[
                  static_cast<std::size_t>(entry.state_index)];
          auto &candidate_normalizers =
              batch_workspace->candidate_ready_expr_normalizers;
          candidate_normalizers.assign(ready_count, 0.0);
          for (std::size_t ready_pos = 0; ready_pos < ready_count;
               ++ready_pos) {
            candidate_normalizers[ready_pos] =
                batch_workspace->ready_values[
                    ready_pos * step_lanes.size() + step_lane];
          }
          advance_exact_sequence_state(
              &candidate_state,
              transition,
              plan,
              step_lanes[step_lane].observed_time,
              ready_count == 0U ? nullptr : &candidate_normalizers);
          append_ranked_frontier_entry(
              &workspace->ranked_next_frontier,
              workspace->ranked_next_states,
              static_cast<semantic::Index>(candidate_index),
              entry.probability * transition_probability);
        }
      }
    }

    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto target_outcome_index = target_by_lane[lane];
      if (target_outcome_index == semantic::kInvalidIndex) {
        continue;
      }
      auto *workspace = lanes[lane].workspace;
      if (workspace->ranked_next_frontier.empty()) {
        alive[lane] = 0U;
        workspace->ranked_frontier.clear();
        continue;
      }
      workspace->ranked_used_outcomes[
          static_cast<std::size_t>(target_outcome_index)] = 1U;
      workspace->ranked_frontier.swap(workspace->ranked_next_frontier);
      workspace->ranked_states.swap(workspace->ranked_next_states);
    }
  }

  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    if (alive[lane] == 0U || lanes[lane].workspace == nullptr) {
      continue;
    }
    double total = 0.0;
    for (const auto &entry : lanes[lane].workspace->ranked_frontier) {
      total += entry.probability;
    }
    (*out_by_lane)[lane] =
        std::isfinite(total) && total > 0.0 ? total : 0.0;
  }
}

inline void exact_ranked_loglik_batch(
    const ExactVariantPlan &plan,
    const std::vector<ExactRankedBatchLane> &lanes,
    const double min_ll,
    ExactRankedBatchWorkspace *batch_workspace,
    std::vector<double> *out_by_lane) {
  const auto lane_count = lanes.size();
  out_by_lane->assign(lane_count, min_ll);
  if (lane_count == 0U) {
    return;
  }
  auto &lane_totals = batch_workspace->lane_totals;
  lane_totals.assign(lane_count, 0.0);

  for (const auto &compiled_state : plan.trigger_state_table.states) {
    auto &sub_lanes = batch_workspace->sub_lanes;
    auto &sub_indices = batch_workspace->sub_lane_indices;
    auto &states = batch_workspace->trigger_states;
    auto &weights = batch_workspace->trigger_weights;
    sub_lanes.clear();
    sub_indices.clear();
    states.clear();
    weights.clear();
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto state =
          exact_compiled_trigger_state_view(
              plan,
              lanes[lane].params,
              lanes[lane].first_param_row,
              compiled_state);
      if (!(state.weight > 0.0)) {
        continue;
      }
      sub_indices.push_back(lane);
      sub_lanes.push_back(lanes[lane]);
      states.push_back(state);
      weights.push_back(state.weight);
    }
    if (sub_lanes.empty()) {
      continue;
    }
    exact_ranked_trigger_probability_batch(
        plan,
        sub_lanes,
        states,
        batch_workspace,
        &batch_workspace->trigger_values);
    for (std::size_t i = 0; i < sub_lanes.size(); ++i) {
      lane_totals[sub_indices[i]] +=
          weights[i] * batch_workspace->trigger_values[i];
    }
  }

  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    const double total = lane_totals[lane];
    (*out_by_lane)[lane] =
        std::isfinite(total) && total > 0.0 ? std::log(total) : min_ll;
  }
}

} // namespace detail
} // namespace accumulatr::eval
