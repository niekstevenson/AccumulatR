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
};

inline ExactSequenceState advance_exact_sequence_state(
    const ExactVariantPlan &plan,
    const ExactSequenceState &sequence_state,
    const ExactTransitionScenario &scenario,
    const double observed_time) {
  ExactSequenceState next_state = sequence_state;
  next_state.lower_bound = observed_time;
  if (next_state.exact_times.empty()) {
    next_state.exact_times.assign(
        static_cast<std::size_t>(plan.source_count),
        std::numeric_limits<double>::quiet_NaN());
  }
  next_state.exact_times[static_cast<std::size_t>(scenario.active_source_id)] =
      observed_time;
  return next_state;
}

inline double evaluate_compiled_scenario_probability(
    const ExactRuntimeScenarioFormula &scenario_formula,
    CompiledSourceView *target_evaluator,
    CompiledEvalWorkspace *target_workspace) {
  if (scenario_formula.probability_root_id == semantic::kInvalidIndex) {
    throw std::runtime_error(
        "runtime scenario reached evaluation without a compiled probability root");
  }
  const double value =
      evaluate_compiled_math_root(
          target_evaluator->plan().compiled_math,
          scenario_formula.probability_root_id,
          &target_workspace->compiled_math,
          target_evaluator,
          nullptr,
          target_workspace);
  return std::isfinite(value) ? value : 0.0;
}

inline ExactStepResult evaluate_exact_step(
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
  ExactStepResult result;
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

  for (std::size_t scenario_idx = 0;
       scenario_idx < plan.outcomes[target_pos].scenarios.size();
       ++scenario_idx) {
    const auto &scenario =
        plan.outcomes[target_pos].scenarios[scenario_idx];
    const auto &runtime_scenario =
        runtime_outcome.scenarios[scenario_idx];
    if (collect_successors && !scenario_supports_ranked_sequence(scenario)) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }

    const double scenario_prob = evaluate_compiled_scenario_probability(
        runtime_scenario,
        &target_evaluator,
        &target_workspace);
    if (!(scenario_prob > 0.0)) {
      continue;
    }
    result.total_probability += scenario_prob;
    if (collect_successors) {
      result.branches.push_back(ExactStepBranch{
          scenario_prob,
          advance_exact_sequence_state(
              plan, sequence_state, scenario, observed_time)});
    }
  }

  return result;
}

} // namespace detail
} // namespace accumulatr::eval
