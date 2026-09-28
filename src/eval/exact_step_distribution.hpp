#pragma once

#include <memory>

#include "exact_compiled_lane_eval.hpp"

namespace accumulatr::eval {
namespace detail {


constexpr std::size_t kExactLaneTileSize = 512U;

struct ExactStepLaneInput {
  const ParamView *params{nullptr};
  const std::uint8_t *shared_started{nullptr};
  const ExactSequenceState *sequence_state{nullptr};
  double observed_time{0.0};
  const std::vector<std::uint8_t> *used_outcomes{nullptr};
};

struct ExactStepLaneWorkspace {
  explicit ExactStepLaneWorkspace(const ExactVariantPlan &plan)
      : executor(plan.compiled_math),
        initial_state(make_exact_sequence_state(plan)),
        source_state(plan, kExactLaneTileSize) {
    step_inputs.reserve(kExactLaneTileSize);
    input_lanes.reserve(kExactLaneTileSize);
    input_positions.reserve(kExactLaneTileSize);
    input_weights.reserve(kExactLaneTileSize);
    executor.lanes.set_source_state(&source_state);
  }

  ExactStepLaneWorkspace(const ExactStepLaneWorkspace &) = delete;
  ExactStepLaneWorkspace &operator=(const ExactStepLaneWorkspace &) = delete;
  ExactStepLaneWorkspace(ExactStepLaneWorkspace &&) = delete;
  ExactStepLaneWorkspace &operator=(ExactStepLaneWorkspace &&) = delete;

  CompiledLaneFrame &prepare(
      const ExactStepLaneInput *lanes,
      const std::size_t lane_count) {
    source_state.bind_matrix(*lanes[0].params->matrix, lane_count);
    auto &frame = executor.lanes.top(lane_count);
    const auto observed = static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed) * frame.stride;
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto &input = lanes[lane];
      frame.has_sequence_history |= input.sequence_state->has_history;
      source_state.bind_lane(lane, *input.sequence_state, *input.params, input.shared_started);
      frame.time_values[observed + lane] = input.observed_time;
      frame.used_outcomes[lane] = input.used_outcomes;
    }
    frame.time_valid[static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed)] = 1U;
    return frame;
  }

  CompiledLaneFrame &prepare_initial(
    const ObservationLaneBatchView lanes,
    const std::uint8_t *shared_started) {
    const auto lane_count = lanes.size;
    source_state.bind_initial_batch(lanes, shared_started);
    auto &frame = executor.lanes.top(lane_count);
    const auto observed = static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed) * frame.stride;
    std::copy_n(
        lanes.observed_times,
        lane_count,
        frame.time_values.data() + observed);
    frame.time_valid[static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed)] = 1U;
    return frame;
  }

  void bind_initial_sources(
      const ObservationLaneBatchView lanes,
      const std::uint8_t *shared_started) {
    source_state.bind_initial_batch(lanes, shared_started);
  }

  template <typename LaneIndex>
  CompiledLaneFrame &prepare_mapped_times(
      const double *times,
      const LaneIndex *source_lanes,
      const std::size_t lane_count) {
    auto &frame = executor.lanes.top(lane_count);
    const auto observed = static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed) * frame.stride;
    std::copy_n(times, lane_count, frame.time_values.data() + observed);
    for (std::size_t lane = 0U; lane < lane_count; ++lane) {
      const auto source =
          static_cast<semantic::Index>(source_lanes[lane]);
      frame.source_lanes[lane] = source;
      frame.source_lanes_identity =
          frame.source_lanes_identity &&
          source == static_cast<semantic::Index>(lane);
    }
    frame.time_valid[static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed)] = 1U;
    return frame;
  }

  CompiledLaneFrame &prepare_repeated_times(
      const std::size_t lane_count,
      const double *times,
      const std::size_t time_count) {
    const auto expanded_count = lane_count * time_count;
    auto &frame = executor.lanes.top(expanded_count);
    frame.source_lanes_identity = time_count == 1U;
    const auto observed = static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed) * frame.stride;
    for (std::size_t time = 0; time < time_count; ++time) {
      std::fill_n(
          frame.time_values.data() + observed + time * lane_count,
          lane_count,
          times[time]);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const auto expanded_lane = time * lane_count + lane;
        if (!frame.source_lanes_identity) {
          frame.source_lanes[expanded_lane] =
              static_cast<semantic::Index>(lane);
        }
      }
    }
    frame.time_valid[static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed)] = 1U;
    return frame;
  }

  CompiledLaneExecutor executor;
  ExactSequenceState initial_state;
  ExactLaneSourceState source_state;
  std::vector<double> root_values;
  std::vector<ExactStepLaneInput> step_inputs;
  ObservationLaneBatch input_lanes;
  std::vector<std::size_t> input_positions;
  std::vector<double> input_weights;
  std::vector<double> trigger_weights;
  std::vector<double> input_values;
  std::vector<double> expanded_values;
  std::vector<double> finite_density_totals;
};

struct ExactStepLaneWorkspacePool {
  explicit ExactStepLaneWorkspacePool(const std::size_t plan_count)
      : workspaces(plan_count) {}

  ExactStepLaneWorkspace &get(
      const std::vector<ExactVariantPlan> &plans,
    const semantic::Index variant_index) {
    const auto position = static_cast<std::size_t>(variant_index);
    if (!workspaces[position]) {
      workspaces[position] =
          std::make_unique<ExactStepLaneWorkspace>(plans[position]);
    }
    return *workspaces[position];
  }

  std::vector<std::unique_ptr<ExactStepLaneWorkspace>> workspaces;
};

inline void evaluate_exact_outcome_density_lanes(
    const ExactVariantPlan &plan,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    const semantic::Index target_idx,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *out) {
  const auto root_id = plan.compiled_outcomes[
      static_cast<std::size_t>(target_idx)].total_probability_root_id;
  const auto &root = plan.compiled_math.roots[static_cast<std::size_t>(root_id)];
  const auto &execution = frame->has_sequence_history
                              ? root.execution : root.initial_execution;
  evaluate_compiled_lane_root(
      plan, root_id, &workspace->executor, frame, lane_count, out);
  if (execution.kind != CompiledMathExecutionKind::SourceProduct) {
    for (auto &value : *out) value = safe_density(value);
  }
}

inline void evaluate_exact_ranked_step_lanes(
    const ExactVariantPlan &plan,
    const ExactStepLaneInput *lanes,
    const std::size_t lane_count,
    const semantic::Index target_idx,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *transition_values,
    std::vector<double> *readiness_values) {
  auto &frame = workspace->prepare(lanes, lane_count);
  const auto &outcome = plan.compiled_outcomes[static_cast<std::size_t>(target_idx)];
  transition_values->resize(outcome.transitions.size() * lane_count);
  for (std::size_t i = 0; i < outcome.transitions.size(); ++i) {
    evaluate_compiled_lane_root(
        plan, outcome.transitions[i].probability_root_id,
        &workspace->executor, &frame, lane_count, &workspace->root_values);
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      (*transition_values)[i * lane_count + lane] =
          safe_density(workspace->root_values[lane]);
    }
  }
  readiness_values->resize(outcome.readiness_root_ids.size() * lane_count);
  for (std::size_t i = 0; i < outcome.readiness_root_ids.size(); ++i) {
    evaluate_compiled_lane_root(
        plan, outcome.readiness_root_ids[i],
        &workspace->executor, &frame, lane_count, &workspace->root_values);
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      (*readiness_values)[i * lane_count + lane] =
          clamp_probability(workspace->root_values[lane]);
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
