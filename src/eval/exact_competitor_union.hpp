#pragma once

#include <optional>

#include "exact_truth.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactStepWorkspace {
  explicit ExactStepWorkspace(const ExactVariantPlan &plan)
      : oracle(plan),
        target_evaluator(plan),
        competitor_evaluator(plan),
        target_workspace(plan),
        competitor_workspace(plan) {}

  void reset(const ParamView &params,
             const int first_param_row,
             const ExactTriggerState &trigger_state,
             const ExactSequenceState &sequence_state,
             const double observed_time) {
    oracle.reset(
        params, first_param_row, trigger_state, sequence_state, observed_time);
    target_evaluator.reset(&oracle, RelationView{});
    competitor_evaluator.reset(&oracle, RelationView{});
    target_workspace.reset(&oracle, RelationView{});
    competitor_workspace.reset(&oracle, RelationView{});
  }

  ExactSourceOracle oracle;
  ForcedExprEvaluator target_evaluator;
  ForcedExprEvaluator competitor_evaluator;
  ForcedExprWorkspace target_workspace;
  ForcedExprWorkspace competitor_workspace;
};

inline double competitor_subset_win_mass(
    const ExactCompetitorSubsetPlan &subset_plan,
    const semantic::Index target_active_source_id,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  double win_prob = 0.0;
  if (subset_plan.outcome_indices.size() == 1U &&
      subset_plan.singleton_expr_root != semantic::kInvalidIndex) {
    win_prob = evaluator->expr_cdf(subset_plan.singleton_expr_root);
  } else {
    for (const auto &scenario : subset_plan.scenarios) {
      const auto scenario_view =
          make_exact_scenario_runtime_view(evaluator->plan(), scenario);
      win_prob += scenario_truth_cdf(scenario_view, evaluator, t, workspace);
      if (!std::isfinite(win_prob)) {
        return 1.0;
      }
    }
    win_prob = clamp_probability(win_prob);
  }
  double same_active_t = 0.0;
  double same_active_ready = 0.0;
  for (const auto &scenario : subset_plan.scenarios) {
    if (scenario.active_source_id != target_active_source_id) {
      continue;
    }
    const auto scenario_view =
        make_exact_scenario_runtime_view(evaluator->plan(), scenario);
    same_active_t +=
        same_active_win_mass(scenario_view, evaluator, t, t, workspace);
    same_active_ready +=
        same_active_win_mass(
            scenario_view, evaluator, t, readiness_upper, workspace);
  }
  return clamp_probability(win_prob - same_active_t + same_active_ready);
}

inline double competitor_non_win_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const semantic::Index target_active_source_id,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const std::vector<std::uint8_t> *used_outcomes = nullptr) {
  double non_win = 1.0;
  for (const auto &block : competitor_plan.blocks) {
    double union_prob = 0.0;
    for (const auto &subset : block.subsets) {
      if (used_outcomes != nullptr) {
        bool skip_subset = false;
        for (const auto outcome_idx : subset.outcome_indices) {
          if ((*used_outcomes)[static_cast<std::size_t>(outcome_idx)] != 0U) {
            skip_subset = true;
            break;
          }
        }
        if (skip_subset) {
          continue;
        }
      }
      union_prob += static_cast<double>(subset.inclusion_sign) *
                    competitor_subset_win_mass(
                        subset,
                        target_active_source_id,
                        readiness_upper,
                        evaluator,
                        workspace,
                        t);
    }
    union_prob = clamp_probability(union_prob);
    non_win *= clamp_probability(1.0 - union_prob);
    if (!(non_win > 0.0)) {
      return 0.0;
    }
  }
  return clamp_probability(non_win);
}

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

inline double evaluate_scenario_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const ExactScenarioRuntimeView &scenario_view,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes,
    const quadrature::FiniteBatch *step_batch,
    ForcedExprEvaluator *target_evaluator,
    ForcedExprWorkspace *target_workspace,
    ForcedExprEvaluator *competitor_evaluator,
    ForcedExprWorkspace *competitor_workspace) {
  const double active_pdf = source_pdf_at(
      target_evaluator->oracle(),
      scenario_view.scenario.active_source_id,
      observed_time);
  const double tail = after_survival(
      scenario_view, target_evaluator, observed_time, target_workspace);
  if (!(active_pdf > 0.0) || !(tail > 0.0)) {
    return 0.0;
  }

  const auto competitor_non_win = [&](const double readiness_upper) {
    return competitor_non_win_probability(
        competitor_plan,
        scenario_view.scenario.active_source_id,
        readiness_upper,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        used_outcomes);
  };

  if (!scenario_view.has_readiness()) {
    return active_pdf * tail * competitor_non_win(0.0);
  }

  double total = 0.0;
  const double initial_ready =
      readiness_cdf(scenario_view, target_evaluator, 0.0, target_workspace);
  if (initial_ready > 0.0) {
    total += initial_ready * active_pdf * tail * competitor_non_win(0.0);
  }

  const auto owned_batch =
      step_batch == nullptr
          ? quadrature::build_finite_batch(0.0, observed_time)
          : quadrature::FiniteBatch{};
  const auto &batch = step_batch != nullptr ? *step_batch : owned_batch;
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double readiness_time = batch.nodes.nodes[i];
    double value = readiness_density(
        scenario_view, target_evaluator, readiness_time, target_workspace);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    value *= active_pdf * tail;
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    value *= competitor_non_win(readiness_time);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    total += batch.nodes.weights[i] * value;
  }
  return total;
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
    const quadrature::FiniteBatch *step_batch = nullptr,
    ExactStepWorkspace *workspace = nullptr) {
  ExactStepResult result;
  const auto target_pos = static_cast<std::size_t>(target_idx);
  const auto &competitor_plan = plan.competitor_plans[target_pos];
  std::optional<ExactStepWorkspace> local_workspace;
  if (workspace == nullptr) {
    local_workspace.emplace(plan);
    workspace = &*local_workspace;
  }
  auto &step_workspace = *workspace;
  step_workspace.reset(
      params, first_param_row, trigger_state, sequence_state, observed_time);
  auto &oracle = step_workspace.oracle;
  auto &target_evaluator = step_workspace.target_evaluator;
  auto &competitor_evaluator = step_workspace.competitor_evaluator;
  auto &target_workspace = step_workspace.target_workspace;
  auto &competitor_workspace = step_workspace.competitor_workspace;

  for (const auto &scenario : plan.outcomes[target_pos].scenarios) {
    if (collect_successors && !scenario_supports_ranked_sequence(scenario)) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }
    const auto scenario_view = make_exact_scenario_runtime_view(plan, scenario);

    RelationView competitor_view;
    if (!relation_view_with_overlay(
            RelationView{}, scenario.relation_template, &competitor_view)) {
      continue;
    }
    competitor_evaluator.reset(&oracle, competitor_view);

    const double scenario_prob = evaluate_scenario_probability(
        competitor_plan,
        scenario_view,
        observed_time,
        used_outcomes,
        step_batch,
        &target_evaluator,
        &target_workspace,
        &competitor_evaluator,
        &competitor_workspace);
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
