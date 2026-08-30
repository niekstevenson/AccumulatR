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
    input_params.reserve(kExactLaneTileSize);
    input_positions.reserve(kExactLaneTileSize);
    input_weights.reserve(kExactLaneTileSize);
    executor.lanes.set_source_state(&source_state);
  }

  ExactStepLaneWorkspace(const ExactStepLaneWorkspace &) = delete;
  ExactStepLaneWorkspace &operator=(const ExactStepLaneWorkspace &) = delete;
  ExactStepLaneWorkspace(ExactStepLaneWorkspace &&) = delete;
  ExactStepLaneWorkspace &operator=(ExactStepLaneWorkspace &&) = delete;

  bool bind_sources(
      const ExactStepLaneInput *lanes,
      const std::size_t lane_count) {
    bool has_sequence_history = false;
    if (lane_count > 0U) {
      source_state.bind_matrix(*lanes[0].params->matrix, lane_count);
    }
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const auto &input = lanes[lane];
      const auto &sequence =
          input.sequence_state == nullptr ? initial_state
                                          : *input.sequence_state;
      has_sequence_history = has_sequence_history || sequence.has_history;
      source_state.bind_lane(
          lane,
          sequence,
          *input.params,
          input.shared_started);
    }
    return has_sequence_history;
  }

  CompiledLaneFrame &prepare(
      const ExactStepLaneInput *lanes,
      const std::size_t lane_count) {
    const bool has_sequence_history = bind_sources(lanes, lane_count);
    auto &frame = executor.lanes.top(lane_count);
    frame.has_sequence_history = has_sequence_history;
    const auto observed = static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed) * frame.stride;
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      frame.time_values[observed + lane] = lanes[lane].observed_time;
      frame.used_outcomes[lane] = lanes[lane].used_outcomes;
    }
    frame.time_valid[static_cast<std::size_t>(
        CompiledMathTimeSlot::Observed)] = 1U;
    return frame;
  }

  CompiledLaneFrame &prepare_initial(
    const ObservationLaneBatchView lanes,
    const std::uint8_t *shared_started) {
    const auto lane_count = lanes.size;
    if (lane_count > 0U) {
      source_state.bind_initial_batch(
          lanes, shared_started);
    }
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
    if (lanes.size > 0U) {
      source_state.bind_initial_batch(lanes, shared_started);
    }
  }

  template <typename LaneIndex>
  CompiledLaneFrame &prepare_mapped_times(
      const double *times,
      const LaneIndex *source_lanes,
      const std::size_t lane_count) {
    auto &frame = executor.lanes.top(lane_count);
    frame.has_sequence_history = false;
    frame.source_lanes_identity = true;
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
    frame.has_sequence_history = false;
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
  std::vector<ParamView> input_params;
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

inline void evaluate_exact_step_distribution_prepared_lanes(
    const ExactVariantPlan &plan,
    CompiledLaneFrame *frame,
    const std::size_t lane_count,
    const semantic::Index target_idx,
    const bool collect_successors,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *total_probability,
    std::vector<double> *transition_probability = nullptr,
    const std::vector<semantic::Index> *additional_root_ids = nullptr,
    std::vector<double> *additional_root_values = nullptr) {
  if (lane_count == 0U || target_idx == semantic::kInvalidIndex) {
    total_probability->assign(lane_count, 0.0);
    if (transition_probability != nullptr) {
      transition_probability->clear();
    }
    if (additional_root_values != nullptr) {
      additional_root_values->clear();
    }
    return;
  }
  const auto &outcome =
      plan.compiled_outcomes[static_cast<std::size_t>(target_idx)];
  const bool collect_transition_values =
      collect_successors && transition_probability != nullptr;
  bool total_probability_is_clean = false;
  if (collect_transition_values) {
    total_probability->assign(lane_count, 0.0);
    transition_probability->assign(
        outcome.transitions.size() * lane_count, 0.0);
    for (std::size_t transition_index = 0;
         transition_index < outcome.transitions.size();
         ++transition_index) {
      const auto root_id =
          outcome.transitions[transition_index].probability_root_id;
      if (root_id == semantic::kInvalidIndex) {
        continue;
      }
      evaluate_compiled_lane_root(
          plan,
          root_id,
          &workspace->executor,
          frame,
          lane_count,
          &workspace->root_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double value = workspace->root_values[lane];
        (*total_probability)[lane] += value;
        (*transition_probability)[transition_index * lane_count + lane] =
            std::isfinite(value) && value > 0.0 ? value : 0.0;
      }
    }
    for (auto &value : *total_probability) {
      value = clean_signed_value(value);
      if (!(value > 0.0)) {
        value = 0.0;
      }
    }
    total_probability_is_clean = true;
  } else if (outcome.total_probability_root_id != semantic::kInvalidIndex) {
    const auto &root = plan.compiled_math.roots[static_cast<std::size_t>(
        outcome.total_probability_root_id)];
    const auto &execution = frame->has_sequence_history
                                ? root.execution
                                : root.initial_execution;
    evaluate_compiled_lane_root(
        plan,
        outcome.total_probability_root_id,
        &workspace->executor,
        frame,
        lane_count,
        total_probability);
    total_probability_is_clean =
        execution.kind == CompiledMathExecutionKind::SourceProduct;
  } else {
    total_probability->assign(lane_count, 0.0);
    for (const auto &transition : outcome.transitions) {
      if (transition.probability_root_id == semantic::kInvalidIndex) {
        continue;
      }
      evaluate_compiled_lane_root(
          plan,
          transition.probability_root_id,
          &workspace->executor,
          frame,
          lane_count,
          &workspace->root_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        (*total_probability)[lane] += workspace->root_values[lane];
      }
    }
  }
  if (!total_probability_is_clean) {
    for (auto &value : *total_probability) {
      if (!std::isfinite(value) || !(value > 0.0)) {
        value = 0.0;
      }
    }
  }
  if (additional_root_ids != nullptr &&
      additional_root_values != nullptr) {
    additional_root_values->assign(
        additional_root_ids->size() * lane_count, 0.0);
    for (std::size_t root_index = 0;
         root_index < additional_root_ids->size();
         ++root_index) {
      const auto root_id = (*additional_root_ids)[root_index];
      if (root_id == semantic::kInvalidIndex) {
        continue;
      }
      evaluate_compiled_lane_root(
          plan,
          root_id,
          &workspace->executor,
          frame,
          lane_count,
          &workspace->root_values);
      for (std::size_t lane = 0; lane < lane_count; ++lane) {
        const double value = workspace->root_values[lane];
        (*additional_root_values)[root_index * lane_count + lane] =
            std::isfinite(value) ? clamp_probability(value) : 0.0;
      }
    }
  }
}

inline void evaluate_exact_step_distribution_lanes(
    const ExactVariantPlan &plan,
    const ExactStepLaneInput *lanes,
    const std::size_t lane_count,
    const semantic::Index target_idx,
    const bool collect_successors,
    ExactStepLaneWorkspace *workspace,
    std::vector<double> *total_probability,
    std::vector<double> *transition_probability = nullptr,
    const std::vector<semantic::Index> *additional_root_ids = nullptr,
    std::vector<double> *additional_root_values = nullptr) {
  if (lane_count == 0U || target_idx == semantic::kInvalidIndex) {
    total_probability->assign(lane_count, 0.0);
    if (transition_probability != nullptr) {
      transition_probability->clear();
    }
    if (additional_root_values != nullptr) {
      additional_root_values->clear();
    }
    return;
  }
  auto &frame = workspace->prepare(lanes, lane_count);
  evaluate_exact_step_distribution_prepared_lanes(
      plan,
      &frame,
      lane_count,
      target_idx,
      collect_successors,
      workspace,
      total_probability,
      transition_probability,
      additional_root_ids,
      additional_root_values);
}

} // namespace detail
} // namespace accumulatr::eval
