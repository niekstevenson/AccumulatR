#pragma once

#include <optional>

#include "exact_truth.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactStepWorkspace {
  explicit ExactStepWorkspace(const ExactVariantPlan &plan)
      : source_channels(plan),
        target_evaluator(plan),
        target_workspace(plan) {}

  void reset(const ParamView &params,
             const int first_param_row,
             const ExactTriggerState &trigger_state,
             const ExactSequenceState &sequence_state,
             const double observed_time) {
    source_channels.reset(
        params, first_param_row, trigger_state, sequence_state, observed_time);
    target_evaluator.reset(&source_channels, RelationView{});
    target_workspace.reset(&source_channels, RelationView{});
    target_workspace.compiled_math.set_time(
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
        observed_time);
    target_workspace.compiled_math.set_time(
        static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
        observed_time);
    target_workspace.compiled_math.set_time(
        static_cast<semantic::Index>(CompiledMathTimeSlot::Zero),
        0.0);
    target_workspace.reset_planned_caches();
  }

  ExactSourceChannels source_channels;
  CompiledSourceView target_evaluator;
  CompiledEvalWorkspace target_workspace;
  std::vector<double> transition_probabilities;
};

inline ExactStepDistributionView evaluate_exact_step_distribution(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const ExactSequenceState &sequence_state,
    const semantic::Index target_idx,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes = nullptr,
    const bool collect_successors = false,
    ExactStepWorkspace *workspace = nullptr) {
  ExactStepDistributionView result;
  const auto target_pos = static_cast<std::size_t>(target_idx);
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  auto &step_workspace = *workspace;
  step_workspace.reset(
      params, first_param_row, trigger_state, sequence_state, observed_time);
  auto &target_evaluator = step_workspace.target_evaluator;
  auto &target_workspace = step_workspace.target_workspace;
  target_workspace.compiled_math.used_outcomes = used_outcomes;
  const auto &runtime_outcome = plan.runtime.outcomes[target_pos];
  const auto &successors = runtime_outcome.successor_distribution;
  step_workspace.transition_probabilities.clear();

  if (!collect_successors &&
      successors.total_probability_root_id != semantic::kInvalidIndex) {
    result.total_probability =
        evaluate_compiled_math_root(
            plan.compiled_math,
            successors.total_probability_root_id,
            &target_workspace.compiled_math,
            &target_evaluator,
            nullptr,
            &target_workspace);
    if (!std::isfinite(result.total_probability) ||
        !(result.total_probability > 0.0)) {
      result.total_probability = 0.0;
    }
    return result;
  }

  step_workspace.transition_probabilities.assign(
      successors.transitions.size(), 0.0);
  for (std::size_t transition_idx = 0;
       transition_idx < successors.transitions.size();
       ++transition_idx) {
    const auto &transition = successors.transitions[transition_idx];
    if (collect_successors && !transition.ranked_supported) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }
    if (transition.probability_root_id == semantic::kInvalidIndex) {
      continue;
    }

    const double scenario_prob =
        evaluate_compiled_math_root(
            plan.compiled_math,
            transition.probability_root_id,
            &target_workspace.compiled_math,
            &target_evaluator,
            nullptr,
            &target_workspace);
    if (!(scenario_prob > 0.0)) {
      continue;
    }
    step_workspace.transition_probabilities[transition_idx] = scenario_prob;
    result.total_probability += scenario_prob;
  }

  if (collect_successors) {
    result.transition_probabilities =
        &step_workspace.transition_probabilities;
  }
  return result;
}

} // namespace detail
} // namespace accumulatr::eval
