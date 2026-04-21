#pragma once

#include "exact_truth.hpp"

namespace accumulatr::eval {
namespace detail {

inline double competitor_subset_win_mass(
    const ExactCompetitorSubsetPlan &subset_plan,
    const ExactSourceKey target_active_key,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    ForcedExprWorkspace *workspace,
    const double t) {
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  double win_prob = 0.0;
  if (subset_plan.outcome_indices.size() == 1U &&
      subset_plan.singleton_expr_root != semantic::kInvalidIndex) {
    win_prob = evaluator->expr_cdf(subset_plan.singleton_expr_root);
  } else {
    for (const auto &scenario : subset_plan.scenarios) {
      win_prob += scenario_truth_cdf(scenario, evaluator, t, workspace);
      if (!std::isfinite(win_prob)) {
        evaluator->oracle_time_ = current_time;
        return 1.0;
      }
    }
    win_prob = clamp_probability(win_prob);
  }
  double same_active_t = 0.0;
  double same_active_ready = 0.0;
  for (const auto &scenario : subset_plan.scenarios) {
    if (!(scenario.active_key == target_active_key)) {
      continue;
    }
    same_active_t += same_active_win_mass(scenario, evaluator, t, t, workspace);
    same_active_ready +=
        same_active_win_mass(scenario, evaluator, t, readiness_upper, workspace);
  }
  evaluator->oracle_time_ = current_time;
  return clamp_probability(win_prob - same_active_t + same_active_ready);
}

inline double competitor_non_win_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const ExactSourceKey target_active_key,
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
                        target_active_key,
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
    const quadrature::FiniteBatch *step_batch = nullptr) {
  ExactStepResult result;
  const auto target_pos = static_cast<std::size_t>(target_idx);
  const auto &competitor_plan = plan.competitor_plans[target_pos];
  ExactSourceOracle oracle(
      plan, params, first_param_row, trigger_state, sequence_state, observed_time);
  ForcedExprEvaluator target_evaluator(plan, &oracle, nullptr);
  std::vector<ExactRelation> competitor_relations(
      static_cast<std::size_t>(plan.source_count),
      ExactRelation::Unknown);
  ForcedExprEvaluator competitor_evaluator(plan, &oracle, nullptr);
  ForcedExprWorkspace target_workspace(plan, &oracle);
  ForcedExprWorkspace competitor_workspace(plan, &oracle);

  for (const auto &scenario : plan.outcomes[target_pos].scenarios) {
    if (collect_successors && !scenario_supports_ranked_sequence(scenario)) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }

    if (!merge_forced_relations(nullptr, scenario, &competitor_relations)) {
      continue;
    }
    target_evaluator.reset_forced_relations(nullptr);
    target_evaluator.oracle_time_ = observed_time;
    competitor_evaluator.reset_forced_relations(&competitor_relations);
    competitor_evaluator.oracle_time_ = observed_time;

    const double active_pdf =
        source_pdf_at(&oracle, scenario.active_key, observed_time);
    const double tail =
        after_survival(scenario, &target_evaluator, observed_time, &target_workspace);
    if (!(active_pdf > 0.0) || !(tail > 0.0)) {
      continue;
    }

    if (scenario.before_keys.empty() && scenario.ready_exprs.empty()) {
      const double scenario_prob =
          active_pdf * tail *
          competitor_non_win_probability(
              competitor_plan,
              scenario.active_key,
              0.0,
              &competitor_evaluator,
              &competitor_workspace,
              observed_time,
              used_outcomes);
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
      continue;
    }

    const double initial_ready =
        readiness_cdf(scenario, &target_evaluator, 0.0, &target_workspace);
    if (initial_ready > 0.0) {
      result.total_probability += initial_ready * active_pdf * tail *
                                  competitor_non_win_probability(
                                      competitor_plan,
                                      scenario.active_key,
                                      0.0,
                                      &competitor_evaluator,
                                      &competitor_workspace,
                                      observed_time,
                                      used_outcomes);
    }

    const auto owned_batch =
        step_batch == nullptr
            ? quadrature::build_finite_batch(0.0, observed_time)
            : quadrature::FiniteBatch{};
    const auto &batch = step_batch != nullptr ? *step_batch : owned_batch;
    double integrated = 0.0;
    for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
      const double readiness_time = batch.nodes.nodes[i];
      double value =
          readiness_density(
              scenario, &target_evaluator, readiness_time, &target_workspace);
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      value *= active_pdf * tail;
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      value *= competitor_non_win_probability(
          competitor_plan,
          scenario.active_key,
          readiness_time,
          &competitor_evaluator,
          &competitor_workspace,
          observed_time,
          used_outcomes);
      if (!std::isfinite(value) || value == 0.0) {
        continue;
      }
      integrated += batch.nodes.weights[i] * value;
    }
    result.total_probability += integrated;
  }

  return result;
}

} // namespace detail
} // namespace accumulatr::eval
