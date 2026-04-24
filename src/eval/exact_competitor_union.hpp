#pragma once

#include <optional>

#include "exact_truth.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactRuntimePreparedCompetitorSubset {
  semantic::Index block_index{0};
  semantic::Index subset_index{0};
  double fixed_mass{0.0};
};

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
  std::vector<semantic::Index> prepared_competitor_block_offsets;
  std::vector<ExactRuntimePreparedCompetitorSubset> prepared_competitor_subsets;
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

inline bool runtime_subset_is_used(
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const std::vector<std::uint8_t> *used_outcomes) {
  if (used_outcomes == nullptr) {
    return false;
  }
  for (const auto outcome_idx : subset.outcome_indices) {
    if ((*used_outcomes)[static_cast<std::size_t>(outcome_idx)] != 0U) {
      return true;
    }
  }
  return false;
}

inline double runtime_competitor_subset_win_at_t(
    const ExactRuntimeCompetitorSubsetPlan &subset,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t) {
  const auto guard = evaluator->oracle()->conditional_time_guard(t);
  if (subset.outcome_indices.size() == 1U &&
      subset.singleton_expr_root != semantic::kInvalidIndex) {
    return evaluator->expr_cdf(subset.singleton_expr_root);
  }

  double win_prob = 0.0;
  for (const auto &scenario : subset.scenarios) {
    win_prob += runtime_scenario_truth_cdf(
        scenario,
        evaluator,
        t,
        workspace);
    if (!std::isfinite(win_prob)) {
      return 1.0;
    }
  }
  return clamp_probability(win_prob);
}

inline double runtime_same_active_mass(
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactRuntimeScenarioSubsetView &subset_view,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const double ready_upper) {
  double total = 0.0;
  for (const auto scenario_idx : subset_view.same_active_scenario_indices) {
    total += runtime_same_active_win_mass(
        subset.scenarios[static_cast<std::size_t>(scenario_idx)],
        evaluator,
        t,
        ready_upper,
        workspace);
  }
  return total;
}

inline void prepare_runtime_competitor_non_win(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    const std::vector<std::uint8_t> *used_outcomes,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    ExactStepWorkspace *step_workspace) {
  auto &offsets = step_workspace->prepared_competitor_block_offsets;
  auto &prepared = step_workspace->prepared_competitor_subsets;
  offsets.clear();
  prepared.clear();
  offsets.reserve(runtime_outcome.competitor_blocks.size() + 1U);
  offsets.push_back(0);

  for (std::size_t block_idx = 0;
       block_idx < runtime_outcome.competitor_blocks.size();
       ++block_idx) {
    const auto &block = runtime_outcome.competitor_blocks[block_idx];
    const auto &block_view = scenario_competitor.blocks[block_idx];
    for (std::size_t subset_idx = 0; subset_idx < block.subsets.size();
         ++subset_idx) {
      const auto &subset = block.subsets[subset_idx];
      if (runtime_subset_is_used(subset, used_outcomes)) {
        continue;
      }
      const auto &subset_view = block_view.subsets[subset_idx];
      const double win_prob =
          runtime_competitor_subset_win_at_t(subset, evaluator, workspace, t);
      const double same_active_t =
          runtime_same_active_mass(subset, subset_view, evaluator, workspace, t, t);
      prepared.push_back(ExactRuntimePreparedCompetitorSubset{
          static_cast<semantic::Index>(block_idx),
          static_cast<semantic::Index>(subset_idx),
          win_prob - same_active_t});
    }
    offsets.push_back(static_cast<semantic::Index>(prepared.size()));
  }
}

inline double runtime_competitor_non_win_probability(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t,
    const double readiness_upper,
    const ExactStepWorkspace &step_workspace) {
  double non_win = 1.0;
  const auto &offsets = step_workspace.prepared_competitor_block_offsets;
  const auto &prepared = step_workspace.prepared_competitor_subsets;
  for (std::size_t block_idx = 0;
       block_idx < runtime_outcome.competitor_blocks.size();
       ++block_idx) {
    double union_prob = 0.0;
    const auto begin = static_cast<std::size_t>(offsets[block_idx]);
    const auto end = static_cast<std::size_t>(offsets[block_idx + 1U]);
    for (std::size_t i = begin; i < end; ++i) {
      const auto &prepared_subset = prepared[i];
      const auto subset_idx = static_cast<std::size_t>(
          prepared_subset.subset_index);
      const auto &subset =
          runtime_outcome.competitor_blocks[block_idx].subsets[subset_idx];
      const auto &subset_view =
          scenario_competitor.blocks[block_idx].subsets[subset_idx];
      const double same_active_ready = runtime_same_active_mass(
          subset,
          subset_view,
          evaluator,
          workspace,
          t,
          readiness_upper);
      const double subset_mass =
          clamp_probability(prepared_subset.fixed_mass + same_active_ready);
      union_prob += static_cast<double>(subset.inclusion_sign) * subset_mass;
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

  const auto batch = quadrature::build_finite_batch(0.0, observed_time);
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

inline double evaluate_runtime_scenario_probability(
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeScenarioFormula &scenario_formula,
    const ExactRuntimeScenarioCompetitorView &scenario_competitor,
    const double observed_time,
    const std::vector<std::uint8_t> *used_outcomes,
    ForcedExprEvaluator *target_evaluator,
    ForcedExprWorkspace *target_workspace,
    ForcedExprEvaluator *competitor_evaluator,
    ForcedExprWorkspace *competitor_workspace,
    ExactStepWorkspace *step_workspace) {
  const double active_pdf = source_pdf_at(
      target_evaluator->oracle(),
      scenario_formula.active_source_id,
      observed_time);
  const double tail = runtime_after_survival(
      scenario_formula, target_evaluator, observed_time, target_workspace);
  if (!(active_pdf > 0.0) || !(tail > 0.0)) {
    return 0.0;
  }

  prepare_runtime_competitor_non_win(
      runtime_outcome,
      scenario_competitor,
      used_outcomes,
      competitor_evaluator,
      competitor_workspace,
      observed_time,
      step_workspace);

  const auto competitor_non_win = [&](const double readiness_upper) {
    return runtime_competitor_non_win_probability(
        runtime_outcome,
        scenario_competitor,
        competitor_evaluator,
        competitor_workspace,
        observed_time,
        readiness_upper,
        *step_workspace);
  };

  if (!scenario_formula.has_readiness) {
    return active_pdf * tail * competitor_non_win(0.0);
  }

  double total = 0.0;
  const double initial_ready = runtime_readiness_cdf(
      scenario_formula, target_evaluator, 0.0, target_workspace);
  if (initial_ready > 0.0) {
    total += initial_ready * active_pdf * tail * competitor_non_win(0.0);
  }

  const auto batch = quadrature::build_finite_batch(0.0, observed_time);
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double readiness_time = batch.nodes.nodes[i];
    double value = runtime_readiness_density(
        scenario_formula, target_evaluator, readiness_time, target_workspace);
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
  auto &oracle = step_workspace.oracle;
  auto &target_evaluator = step_workspace.target_evaluator;
  auto &competitor_evaluator = step_workspace.competitor_evaluator;
  auto &target_workspace = step_workspace.target_workspace;
  auto &competitor_workspace = step_workspace.competitor_workspace;
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

    RelationView competitor_view;
    if (!relation_view_with_overlay(
            RelationView{}, runtime_scenario.relation_template, &competitor_view)) {
      continue;
    }
    competitor_evaluator.reset(&oracle, competitor_view);

    const double scenario_prob = evaluate_runtime_scenario_probability(
        runtime_outcome,
        runtime_scenario,
        runtime_outcome.competitor_by_scenario[scenario_idx],
        observed_time,
        used_outcomes,
        &target_evaluator,
        &target_workspace,
        &competitor_evaluator,
        &competitor_workspace,
        &step_workspace);
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
