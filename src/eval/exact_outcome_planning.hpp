#pragma once

#include <algorithm>

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_expr_distribution.hpp"
#include "exact_transition_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index exact_terminal_leaf_release(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionScenario &scenario) {
  const auto source_id =
      scenario.transition.release_source_id;
  if (source_id == semantic::kInvalidIndex ||
      source_id >= plan.program.layout.n_leaves ||
      plan.source_count != plan.program.layout.n_leaves) {
    return semantic::kInvalidIndex;
  }
  if (!scenario.transition.readiness.empty() ||
      !scenario.transition.guards.empty() ||
      !scenario.transition.source_order_facts.empty()) {
    return semantic::kInvalidIndex;
  }
  const auto &relations = scenario.transition.relation_template;
  if (!relations.empty() &&
      !(relations.source_ids.size() == 1U &&
        relations.relations.size() == 1U &&
        relations.source_ids.front() == source_id &&
        relations.relations.front() == ExactRelation::At)) {
    return semantic::kInvalidIndex;
  }
  return source_id;
}

inline ExactTerminalNoResponsePlan compile_terminal_no_response_plan(
    const ExactVariantBuildState &plan) {
  ExactTerminalNoResponsePlan no_response;
  const auto leaf_count = plan.program.layout.n_leaves;
  if (leaf_count <= 0 ||
      plan.source_count != leaf_count ||
      plan.outcomes.empty()) {
    return no_response;
  }

  std::vector<std::uint8_t> covered(static_cast<std::size_t>(leaf_count), 0U);
  for (const auto &outcome : plan.outcomes) {
    if (outcome.scenarios.size() != 1U) {
      return ExactTerminalNoResponsePlan{};
    }
    const auto leaf_index =
        exact_terminal_leaf_release(plan, outcome.scenarios.front());
    if (leaf_index == semantic::kInvalidIndex) {
      return ExactTerminalNoResponsePlan{};
    }
    const auto pos = static_cast<std::size_t>(leaf_index);
    if (covered[pos] != 0U) {
      return ExactTerminalNoResponsePlan{};
    }
    covered[pos] = 1U;
  }

  no_response.leaf_indices.reserve(covered.size());
  for (std::size_t i = 0; i < covered.size(); ++i) {
    if (covered[i] == 0U) {
      return ExactTerminalNoResponsePlan{};
    }
    no_response.leaf_indices.push_back(static_cast<semantic::Index>(i));
  }
  no_response.direct_leaf_failure_product = true;
  return no_response;
}

inline semantic::Index compile_outcome_probability_root(
    ExactVariantBuildState *plan,
    const ExactOutcomeRegionCompileContext &outcome_context) {
  if (outcome_context.scenarios.empty()) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
  }
  std::vector<semantic::Index> scenario_nodes;
  scenario_nodes.reserve(outcome_context.scenarios.size());
  for (const auto &scenario : outcome_context.scenarios) {
    scenario_nodes.push_back(
        compiled_math_root_node_id(
            plan->compiled_math,
            scenario.probability_root_id));
  }
  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::move(scenario_nodes),
          CompiledMathValueKind::Scalar));
}

inline void mark_sequence_expr_upper_bounds_for_scenario(
    ExactVariantBuildState *plan,
    const ExactSymbolicTransitionScenario &scenario) {
  for (const auto &guard :
       scenario.transition.readiness.guards) {
    if (guard.kind != ExactTransitionGuardKind::ExprBefore) {
      continue;
    }
    const auto expr_id = guard.subject_id;
    if (expr_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(expr_id) <
            plan->sequence.expr_upper_bound_used.size()) {
      plan->sequence.expr_upper_bound_used[
          static_cast<std::size_t>(expr_id)] = 1U;
    }
  }
}

inline void compile_sequence_plan(
    ExactVariantBuildState *plan,
    const std::vector<ExactTargetCompetitorPlan> &competitor_plans) {
  const auto expr_count = plan->program.expr_kind.size();
  plan->sequence.expr_upper_bound_used.assign(expr_count, 0U);
  plan->sequence.expr_cdf_roots.assign(
      expr_count,
      semantic::kInvalidIndex);
  for (const auto &outcome : plan->outcomes) {
    for (const auto &scenario : outcome.scenarios) {
      mark_sequence_expr_upper_bounds_for_scenario(plan, scenario);
    }
  }
  for (const auto &target_plan : competitor_plans) {
    for (const auto &competitor : target_plan.competitors) {
      for (const auto &scenario : competitor.scenarios) {
        mark_sequence_expr_upper_bounds_for_scenario(plan, scenario);
      }
    }
  }
  for (semantic::Index expr_id = 0;
       expr_id < static_cast<semantic::Index>(expr_count);
       ++expr_id) {
    if (plan->sequence.expr_upper_bound_used[
            static_cast<std::size_t>(expr_id)] == 0U) {
      continue;
    }
    const auto node_id =
        compile_expr_value_node_raw(
            plan,
            expr_id,
            CompiledMathValueKind::Cdf,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
            0);
    plan->sequence.expr_cdf_roots[static_cast<std::size_t>(expr_id)] =
        compiled_math_make_root(&plan->compiled_math, node_id);
  }
}

inline void compile_finite_response_distribution_roots(
    ExactVariantBuildState *plan) {
  if (plan->compiled_outcomes.empty()) {
    return;
  }
  if (plan->no_response.direct_leaf_failure_product) {
    std::vector<semantic::Index> survival_nodes;
    survival_nodes.reserve(plan->no_response.leaf_indices.size());
    for (const auto leaf : plan->no_response.leaf_indices) {
      survival_nodes.push_back(compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceSurvival,
          leaf));
    }
    plan->finite_response_survival_root_id = compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(survival_nodes),
            CompiledMathValueKind::Survival));
    return;
  }
  std::vector<semantic::Index> outcome_density_nodes;
  outcome_density_nodes.reserve(plan->compiled_outcomes.size());
  for (const auto &outcome : plan->compiled_outcomes) {
    outcome_density_nodes.push_back(compiled_math_root_node_id(
        plan->compiled_math, outcome.total_probability_root_id));
  }
  plan->finite_response_density_root_id = compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::move(outcome_density_nodes),
          CompiledMathValueKind::Density));
}

inline std::vector<ExactCompiledOutcomePlan> compile_exact_outcome_plans(
    ExactVariantBuildState *plan,
    const std::vector<ExactTargetCompetitorPlan> &competitor_plans) {
  const ExactVariantBuildState &plan_ref = *plan;
  std::vector<ExactCompiledOutcomePlan> compiled_outcomes;
  compiled_outcomes.reserve(plan_ref.outcomes.size());

  for (semantic::Index target_idx = 0;
       target_idx < static_cast<semantic::Index>(plan_ref.outcomes.size());
       ++target_idx) {
    const auto target_pos = static_cast<std::size_t>(target_idx);
    const auto &outcome = plan_ref.outcomes[target_pos];
    const auto &competitor_plan = competitor_plans[target_pos];

    ExactOutcomeRegionCompileContext compile_context;
    compile_context.scenarios = outcome.scenarios;
    compile_context.competitors = competitor_plan.competitors;

    for (auto &scenario : compile_context.scenarios) {
      scenario.probability_root_id =
          exact_order_region_probability_root(plan, compile_context, scenario);
    }
    ExactCompiledOutcomePlan compiled_outcome;
    compiled_outcome.total_probability_root_id =
        compile_outcome_probability_root(plan, compile_context);
    compiled_outcome.transitions.reserve(compile_context.scenarios.size());
    compiled_outcome.transition_readiness_slots.reserve(
        compile_context.scenarios.size());
    for (std::size_t scenario_idx = 0;
         scenario_idx < compile_context.scenarios.size();
         ++scenario_idx) {
      ExactCompiledTransitionPlan transition;
      transition.probability_root_id =
          compile_context.scenarios[scenario_idx].probability_root_id;
      transition.release_source_id =
          compile_context.scenarios[scenario_idx].transition.release_source_id;
      const auto readiness_offset = static_cast<semantic::Index>(
          compiled_outcome.readiness_root_slot_by_item.size());
      for (const auto &guard :
           compile_context.scenarios[scenario_idx]
               .transition.readiness.guards) {
        if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
          transition.readiness_source_ids.push_back(guard.subject_id);
        } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
          transition.readiness_expr_ids.push_back(guard.subject_id);
          const auto root_id =
              guard.subject_id == semantic::kInvalidIndex ||
                      static_cast<std::size_t>(guard.subject_id) >=
                          plan->sequence.expr_cdf_roots.size()
                  ? semantic::kInvalidIndex
                  : plan->sequence.expr_cdf_roots[
                        static_cast<std::size_t>(guard.subject_id)];
          semantic::Index root_slot = semantic::kInvalidIndex;
          if (root_id != semantic::kInvalidIndex) {
            const auto existing = std::find(
                compiled_outcome.readiness_root_ids.begin(),
                compiled_outcome.readiness_root_ids.end(),
                root_id);
            if (existing == compiled_outcome.readiness_root_ids.end()) {
              root_slot = static_cast<semantic::Index>(
                  compiled_outcome.readiness_root_ids.size());
              compiled_outcome.readiness_root_ids.push_back(root_id);
            } else {
              root_slot = static_cast<semantic::Index>(std::distance(
                  compiled_outcome.readiness_root_ids.begin(), existing));
            }
          }
          compiled_outcome.readiness_root_slot_by_item.push_back(root_slot);
        }
      }
      compiled_outcome.transition_readiness_slots.push_back(ExactIndexSpan{
          readiness_offset,
          static_cast<semantic::Index>(transition.readiness_expr_ids.size())});
      compiled_outcome.transitions.push_back(std::move(transition));
    }

    compiled_outcomes.push_back(std::move(compiled_outcome));
  }

  return compiled_outcomes;
}


} // namespace detail
} // namespace accumulatr::eval
