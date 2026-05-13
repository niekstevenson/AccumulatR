#pragma once

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_expr_distribution.hpp"
#include "exact_transition_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

inline bool exact_scenario_is_terminal_leaf_release(
    const ExactVariantBuildState &plan,
    const ExactSymbolicTransitionScenario &scenario,
    semantic::Index *leaf_index) {
  const auto source_id =
      exact_symbolic_transition_release_source_id(scenario.transition);
  if (source_id == semantic::kInvalidIndex ||
      source_id >= plan.lowered.program.layout.n_leaves ||
      plan.source_count != plan.lowered.program.layout.n_leaves) {
    return false;
  }
  if (!scenario.transition.readiness_time_expr.requirements.empty() ||
      !scenario.transition.guards.empty() ||
      !scenario.transition.order_region.source_order_facts.empty()) {
    return false;
  }
  const auto &relations = scenario.transition.relation_template;
  if (!relations.empty() &&
      !(relations.source_ids.size() == 1U &&
        relations.relations.size() == 1U &&
        relations.source_ids.front() == source_id &&
        relations.relations.front() == ExactRelation::At)) {
    return false;
  }
  if (!scenario.transition.active_sources.empty() &&
      (scenario.transition.active_sources.size() != 1U ||
       scenario.transition.active_sources.front() != source_id)) {
    return false;
  }
  if (leaf_index != nullptr) {
    *leaf_index = source_id;
  }
  return true;
}

inline ExactTerminalNoResponsePlan compile_terminal_no_response_plan(
    const ExactVariantBuildState &plan) {
  ExactTerminalNoResponsePlan no_response;
  const auto leaf_count = plan.lowered.program.layout.n_leaves;
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
    semantic::Index leaf_index{semantic::kInvalidIndex};
    if (!exact_scenario_is_terminal_leaf_release(
            plan, outcome.scenarios.front(), &leaf_index)) {
      return ExactTerminalNoResponsePlan{};
    }
    const auto pos = static_cast<std::size_t>(leaf_index);
    if (pos >= covered.size() || covered[pos] != 0U) {
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

inline semantic::Index compile_scenario_probability_root(
    ExactVariantBuildState *plan,
    const ExactOutcomeRegionCompileContext &outcome_context,
    ExactSymbolicTransitionScenario *formula) {
  semantic::Index region_root{semantic::kInvalidIndex};
  if (exact_order_region_probability_root(
          plan, outcome_context, *formula, &region_root)) {
    return region_root;
  }
  throw std::runtime_error(
      "exact order-region compiler could not lower scenario probability");
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
    if (scenario.probability_root_id == semantic::kInvalidIndex) {
      continue;
    }
    const auto node_id =
        compiled_math_root_node_id(
            plan->compiled_math,
            scenario.probability_root_id);
    if (node_id != semantic::kInvalidIndex) {
      scenario_nodes.push_back(node_id);
    }
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
       scenario.transition.readiness_time_expr.requirements.guards) {
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
  const auto expr_count = plan->lowered.program.expr_kind.size();
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
            CompiledMathNodeKind::ExprCdf,
            0,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
            0);
    plan->sequence.expr_cdf_roots[static_cast<std::size_t>(expr_id)] =
        compiled_math_make_root(&plan->compiled_math, node_id);
  }
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
    compile_context.scenarios.reserve(outcome.scenarios.size());
    for (const auto &scenario : outcome.scenarios) {
      compile_context.scenarios.push_back(scenario);
    }

    compile_context.competitors = competitor_plan.competitors;

    for (std::size_t scenario_idx = 0;
         scenario_idx < compile_context.scenarios.size();
         ++scenario_idx) {
      compile_context.scenarios[scenario_idx].probability_root_id =
          compile_scenario_probability_root(
              plan,
              compile_context,
              &compile_context.scenarios[scenario_idx]);
    }
    ExactCompiledOutcomePlan compiled_outcome;
    compiled_outcome.total_probability_root_id =
        compile_outcome_probability_root(plan, compile_context);
    compiled_outcome.transitions.reserve(compile_context.scenarios.size());
    for (std::size_t scenario_idx = 0;
         scenario_idx < compile_context.scenarios.size();
         ++scenario_idx) {
      ExactCompiledTransitionPlan transition;
      transition.probability_root_id =
          compile_context.scenarios[scenario_idx].probability_root_id;
      transition.release_source_id =
          exact_symbolic_transition_release_source_id(
              compile_context.scenarios[scenario_idx].transition);
      for (const auto &guard :
           compile_context.scenarios[scenario_idx]
               .transition.readiness_time_expr.requirements.guards) {
        if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
          transition.readiness_source_ids.push_back(guard.subject_id);
        } else if (guard.kind == ExactTransitionGuardKind::ExprBefore) {
          transition.readiness_expr_ids.push_back(guard.subject_id);
        }
      }
      compiled_outcome.transitions.push_back(std::move(transition));
    }

    compiled_outcomes.push_back(std::move(compiled_outcome));
  }

  return compiled_outcomes;
}


} // namespace detail
} // namespace accumulatr::eval
