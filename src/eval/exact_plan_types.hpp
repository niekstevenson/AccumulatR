#pragma once

#include <limits>
#include <utility>
#include <vector>

#include "compiled_math_types.hpp"
#include "exact_common_types.hpp"
#include "exact_metrics.hpp"
#include "exact_region_types.hpp"
#include "exact_support_types.hpp"
#include "exact_transition_types.hpp"
#include "../runtime/exact_program.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactOutcomePlan {
  std::vector<ExactSymbolicTransitionScenario> scenarios;
};

struct ExactCompetitorRegionPlan {
  semantic::Index expr_root{semantic::kInvalidIndex};
  std::vector<semantic::Index> outcome_indices;
  std::vector<ExactSymbolicTransitionScenario> scenarios;
};

struct ExactCompiledTransitionPlan {
  semantic::Index probability_root_id{semantic::kInvalidIndex};
  semantic::Index release_source_id{semantic::kInvalidIndex};
  std::vector<semantic::Index> readiness_source_ids;
  std::vector<semantic::Index> readiness_expr_ids;
};

struct ExactOutcomeRegionCompileContext {
  std::vector<ExactSymbolicTransitionScenario> scenarios;
  std::vector<ExactCompetitorRegionPlan> competitors;
};

struct ExactCompiledOutcomePlan {
  std::vector<ExactCompiledTransitionPlan> transitions;
  semantic::Index total_probability_root_id{semantic::kInvalidIndex};
};

struct ExactSequencePlan {
  std::vector<std::uint8_t> expr_upper_bound_used;
  std::vector<semantic::Index> expr_cdf_roots;
};

struct ExactExprDistributionKey {
  semantic::Index expr_id{semantic::kInvalidIndex};
  CompiledMathNodeKind value_kind{CompiledMathNodeKind::ExprCdf};
  semantic::Index condition_id{0};
  semantic::Index time_id{
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed)};
  semantic::Index source_view_id{0};
};

struct ExactExprDistributionPlan {
  ExactExprDistributionKey key;
  semantic::Index root_id{semantic::kInvalidIndex};
  bool compiling{false};
};

struct ExactVariantBuildState {
  runtime::LoweredExactVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactOutcomePlan> outcomes;
  std::vector<ExactCompiledOutcomePlan> compiled_outcomes;
  ExactSequencePlan sequence;
  std::vector<ExactExprDistributionPlan> expr_distributions;
  ExactComplexityMetrics complexity;
  CompiledMathProgram compiled_math;
  std::vector<ExactRelationTemplate> compiled_source_views;
  std::vector<ExactExprKernel> expr_kernels;
  std::vector<ExactSourceKernel> source_kernels;
  std::vector<std::vector<semantic::Index>> leaf_supports;
  std::vector<std::vector<semantic::Index>> pool_supports;
  std::vector<std::vector<semantic::Index>> expr_supports;
  std::vector<semantic::Index> compiled_outcome_gate_indices;
  semantic::Index source_count{0};
  std::vector<semantic::Index> leaf_source_ids;
  std::vector<semantic::Index> pool_source_ids;
  std::vector<semantic::Index> shared_trigger_indices;
  ExactCompiledTriggerStateTable trigger_state_table;
  std::vector<std::uint8_t> compiled_source_view_relations;
  semantic::Index compiled_source_view_source_count{0};
  semantic::Index next_compiled_time_id{
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) + 1U};
};

struct ExactVariantPlan {
  runtime::LoweredExactVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactCompiledOutcomePlan> compiled_outcomes;
  ExactSequencePlan sequence;
  ExactComplexityMetrics complexity;
  CompiledMathProgram compiled_math;
  std::vector<ExactRelationTemplate> compiled_source_views;
  std::vector<semantic::Index> compiled_outcome_gate_indices;
  semantic::Index source_count{0};
  ExactCompiledTriggerStateTable trigger_state_table;
  std::vector<std::uint8_t> compiled_source_view_relations;
  semantic::Index compiled_source_view_source_count{0};
};

inline ExactVariantPlan finalize_exact_variant_plan(
    ExactVariantBuildState &&build) {
  ExactVariantPlan plan;
  plan.lowered = std::move(build.lowered);
  plan.outcome_index_by_code = std::move(build.outcome_index_by_code);
  plan.compiled_outcomes = std::move(build.compiled_outcomes);
  plan.sequence = std::move(build.sequence);
  plan.complexity = build.complexity;
  plan.compiled_math = std::move(build.compiled_math);
  plan.compiled_source_views = std::move(build.compiled_source_views);
  plan.compiled_outcome_gate_indices =
      std::move(build.compiled_outcome_gate_indices);
  plan.source_count = build.source_count;
  plan.trigger_state_table = std::move(build.trigger_state_table);
  plan.compiled_source_view_relations =
      std::move(build.compiled_source_view_relations);
  plan.compiled_source_view_source_count =
      build.compiled_source_view_source_count;
  return plan;
}

inline ExactSequenceState make_exact_sequence_state(const ExactVariantPlan &plan) {
  ExactSequenceState state;
  state.exact_times.assign(
      static_cast<std::size_t>(plan.source_count),
      std::numeric_limits<double>::quiet_NaN());
  state.upper_bounds.assign(
      static_cast<std::size_t>(plan.source_count),
      std::numeric_limits<double>::infinity());
  state.expr_upper_bounds.assign(
      plan.lowered.program.expr_kind.size(),
      std::numeric_limits<double>::infinity());
  state.expr_upper_normalizers.assign(
      plan.lowered.program.expr_kind.size(),
      0.0);
  return state;
}

template <typename PlanLike>
inline ExactRelation exact_compiled_source_view_relation(
    const PlanLike &plan,
    const semantic::Index source_view_id,
    const semantic::Index source_id) noexcept {
  if (source_view_id == 0 ||
      source_view_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex ||
      plan.compiled_source_view_source_count <= 0) {
    return ExactRelation::Unknown;
  }
  const auto source_count =
      static_cast<std::size_t>(plan.compiled_source_view_source_count);
  const auto view_pos = static_cast<std::size_t>(source_view_id - 1U);
  const auto source_pos = static_cast<std::size_t>(source_id);
  const auto offset = view_pos * source_count + source_pos;
  if (source_pos >= source_count ||
      offset >= plan.compiled_source_view_relations.size()) {
    return ExactRelation::Unknown;
  }
  return static_cast<ExactRelation>(
      plan.compiled_source_view_relations[offset]);
}

inline bool expr_support_contains_source(const ExactVariantBuildState &plan,
                                         const semantic::Index expr_idx,
                                         const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  return support_contains_source(
      plan.expr_supports[static_cast<std::size_t>(expr_idx)], source_id);
}

inline bool expr_supports_overlap(const ExactVariantBuildState &plan,
                                  const semantic::Index lhs_expr_idx,
                                  const semantic::Index rhs_expr_idx) {
  if (lhs_expr_idx == semantic::kInvalidIndex ||
      rhs_expr_idx == semantic::kInvalidIndex) {
    return false;
  }
  return supports_overlap(
      plan.expr_supports[static_cast<std::size_t>(lhs_expr_idx)],
      plan.expr_supports[static_cast<std::size_t>(rhs_expr_idx)]);
}

inline semantic::Index source_ordinal(const ExactVariantBuildState &plan,
                                      const semantic::SourceKind kind,
                                      const semantic::Index index) {
  if (kind == semantic::SourceKind::Leaf) {
    return plan.leaf_source_ids[static_cast<std::size_t>(index)];
  }
  if (kind == semantic::SourceKind::Pool) {
    return plan.pool_source_ids[static_cast<std::size_t>(index)];
  }
  return semantic::kInvalidIndex;
}

} // namespace detail
} // namespace accumulatr::eval
