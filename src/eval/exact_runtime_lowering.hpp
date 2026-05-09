#pragma once

#include "exact_types.hpp"
#include "exact_compiled_math_lowering.hpp"
#include "exact_region_compiler.hpp"
#include "exact_transition_lowering.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index compile_runtime_scenario_probability_root(
    ExactVariantPlan *plan,
    const ExactRuntimeOutcomeCompileContext &runtime_outcome,
    ExactSymbolicTransitionScenario *formula) {
  semantic::Index region_root{semantic::kInvalidIndex};
  if (exact_order_region_probability_root(
          plan, runtime_outcome, *formula, &region_root)) {
    return region_root;
  }
  throw std::runtime_error(
      "exact order-region compiler could not lower scenario probability");
}

inline semantic::Index compile_runtime_outcome_probability_root(
    ExactVariantPlan *plan,
    const ExactRuntimeOutcomeCompileContext &runtime_outcome) {
  if (runtime_outcome.scenarios.empty()) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
  }
  std::vector<semantic::Index> scenario_nodes;
  scenario_nodes.reserve(runtime_outcome.scenarios.size());
  for (const auto &scenario : runtime_outcome.scenarios) {
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
    ExactVariantPlan *plan,
    const ExactSymbolicTransitionScenario &scenario) {
  for (const auto &guard :
       scenario.transition.readiness_time_expr.requirements.guards) {
    if (guard.kind != ExactTransitionGuardKind::ExprBefore) {
      continue;
    }
    const auto expr_id = guard.subject_id;
    if (expr_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(expr_id) <
            plan->sequence_expr_upper_bound_used.size()) {
      plan->sequence_expr_upper_bound_used[
          static_cast<std::size_t>(expr_id)] = 1U;
    }
  }
}

inline void compile_sequence_expr_upper_bound_roots(
    ExactVariantPlan *plan,
    const std::vector<ExactTargetCompetitorPlan> &competitor_plans) {
  const auto expr_count = plan->lowered.program.expr_kind.size();
  plan->sequence_expr_upper_bound_used.assign(expr_count, 0U);
  plan->sequence_expr_cdf_roots.assign(
      expr_count,
      semantic::kInvalidIndex);
  for (const auto &outcome : plan->outcomes) {
    for (const auto &scenario : outcome.scenarios) {
      mark_sequence_expr_upper_bounds_for_scenario(plan, scenario);
    }
  }
  for (const auto &target_plan : competitor_plans) {
    for (const auto &block : target_plan.blocks) {
      for (const auto &subset : block.subsets) {
        for (const auto &scenario : subset.scenarios) {
          mark_sequence_expr_upper_bounds_for_scenario(plan, scenario);
        }
      }
    }
  }
  for (semantic::Index expr_id = 0;
       expr_id < static_cast<semantic::Index>(expr_count);
       ++expr_id) {
    if (plan->sequence_expr_upper_bound_used[
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
    plan->sequence_expr_cdf_roots[static_cast<std::size_t>(expr_id)] =
        compiled_math_make_root(&plan->compiled_math, node_id);
  }
}

inline ExactRuntimeVariantPlan compile_exact_runtime_plan(
    ExactVariantPlan *plan,
    const std::vector<ExactTargetCompetitorPlan> &competitor_plans) {
  const ExactVariantPlan &plan_ref = *plan;
  ExactRuntimeVariantPlan runtime;
  runtime.outcomes.reserve(plan_ref.outcomes.size());

  for (semantic::Index target_idx = 0;
       target_idx < static_cast<semantic::Index>(plan_ref.outcomes.size());
       ++target_idx) {
    const auto target_pos = static_cast<std::size_t>(target_idx);
    const auto &outcome = plan_ref.outcomes[target_pos];
    const auto &competitor_plan = competitor_plans[target_pos];

    ExactRuntimeOutcomeCompileContext compile_context;
    compile_context.scenarios.reserve(outcome.scenarios.size());
    for (const auto &scenario : outcome.scenarios) {
      compile_context.scenarios.push_back(scenario);
    }

    compile_context.competitor_blocks.reserve(competitor_plan.blocks.size());
    for (const auto &block : competitor_plan.blocks) {
      ExactRuntimeCompetitorBlockPlan runtime_block;
      runtime_block.subsets.reserve(block.subsets.size());
      for (const auto &subset : block.subsets) {
        ExactRuntimeCompetitorSubsetPlan runtime_subset;
        runtime_subset.outcome_indices = subset.outcome_indices;
        runtime_subset.expr_roots = subset.expr_roots;
        runtime_subset.inclusion_sign = subset.inclusion_sign;
        runtime_subset.scenarios.reserve(subset.scenarios.size());
        for (const auto &scenario : subset.scenarios) {
          runtime_subset.scenarios.push_back(scenario);
        }
        runtime_block.subsets.push_back(std::move(runtime_subset));
      }
      compile_context.competitor_blocks.push_back(std::move(runtime_block));
    }

    for (std::size_t scenario_idx = 0;
         scenario_idx < compile_context.scenarios.size();
         ++scenario_idx) {
      compile_context.scenarios[scenario_idx].probability_root_id =
          compile_runtime_scenario_probability_root(
              plan,
              compile_context,
              &compile_context.scenarios[scenario_idx]);
    }
    ExactRuntimeOutcomePlan runtime_outcome;
    runtime_outcome.successor_distribution.total_probability_root_id =
        compile_runtime_outcome_probability_root(plan, compile_context);
    runtime_outcome.successor_distribution.transitions.reserve(
        compile_context.scenarios.size());
    for (std::size_t scenario_idx = 0;
         scenario_idx < compile_context.scenarios.size();
         ++scenario_idx) {
      ExactRuntimeScenarioTransitionPlan transition;
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
      runtime_outcome.successor_distribution.transitions.push_back(
          std::move(transition));
    }

    runtime.outcomes.push_back(std::move(runtime_outcome));
  }

  return runtime;
}


} // namespace detail
} // namespace accumulatr::eval
