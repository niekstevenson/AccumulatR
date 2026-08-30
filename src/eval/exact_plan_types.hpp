#pragma once

#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "compiled_math_types.hpp"
#include "exact_common_types.hpp"
#include "exact_metrics.hpp"
#include "exact_region_types.hpp"
#include "exact_support_types.hpp"
#include "exact_transition_types.hpp"
#include "../runtime/exact_evaluation_program.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactProjectionPlannerState;

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
  std::vector<semantic::Index> readiness_root_ids;
  std::vector<ExactIndexSpan> transition_readiness_slots;
  std::vector<semantic::Index> readiness_root_slot_by_item;
};

struct ExactSequencePlan {
  std::vector<std::uint8_t> expr_upper_bound_used;
  std::vector<semantic::Index> expr_cdf_roots;
};

struct ExactTerminalNoResponsePlan {
  bool direct_leaf_failure_product{false};
  std::vector<semantic::Index> leaf_indices;
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

struct ExactSourceProgramCompileKey {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index source_view_id{0};

  bool operator==(const ExactSourceProgramCompileKey &other) const noexcept {
    return source_id == other.source_id &&
           condition_id == other.condition_id &&
           source_view_id == other.source_view_id;
  }
};

struct ExactSourceProgramCompileKeyHash {
  std::size_t operator()(
      const ExactSourceProgramCompileKey &key) const noexcept {
    std::size_t seed = static_cast<std::size_t>(key.source_id);
    seed ^= static_cast<std::size_t>(key.condition_id) +
            0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    seed ^= static_cast<std::size_t>(key.source_view_id) +
            0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
  }
};

struct ExactConditionedSourceProgramCompileKey {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index condition_id{0};
  semantic::Index source_view_id{0};
  semantic::Index time_id{0};
  semantic::Index time_cap_id{semantic::kInvalidIndex};

  bool operator==(
      const ExactConditionedSourceProgramCompileKey &other) const noexcept {
    return source_id == other.source_id &&
           condition_id == other.condition_id &&
           source_view_id == other.source_view_id &&
           time_id == other.time_id &&
           time_cap_id == other.time_cap_id;
  }
};

struct ExactConditionedSourceProgramCompileKeyHash {
  std::size_t operator()(
      const ExactConditionedSourceProgramCompileKey &key) const noexcept {
    std::size_t seed = static_cast<std::size_t>(key.source_id);
    hash_combine(&seed, static_cast<std::size_t>(key.condition_id));
    hash_combine(&seed, static_cast<std::size_t>(key.source_view_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_id));
    hash_combine(&seed, static_cast<std::size_t>(key.time_cap_id));
    return seed;
  }

private:
  static void hash_combine(std::size_t *seed, const std::size_t value) noexcept {
    *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
  }
};

struct ExactVariantBuildState {
  runtime::ExactEvaluationProgram program;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactOutcomePlan> outcomes;
  std::vector<ExactCompiledOutcomePlan> compiled_outcomes;
  ExactSequencePlan sequence;
  ExactTerminalNoResponsePlan no_response;
  semantic::Index finite_response_density_root_id{
      semantic::kInvalidIndex};
  semantic::Index finite_response_survival_root_id{
      semantic::kInvalidIndex};
  std::vector<ExactExprDistributionPlan> expr_distributions;
  mutable std::shared_ptr<ExactProjectionPlannerState> projection_planner;
  ExactComplexityMetrics *complexity_metrics{nullptr};
  CompiledMathProgram compiled_math;
  std::vector<ExactRelationTemplate> compiled_source_views;
  std::vector<ExactExprKernel> expr_kernels;
  std::vector<ExactSourceKernel> source_kernels;
  std::vector<std::vector<semantic::Index>> leaf_supports;
  std::vector<std::vector<semantic::Index>> pool_supports;
  std::vector<std::vector<semantic::Index>> expr_supports;
  std::vector<std::uint8_t> pool_transition_can_stay_aggregate;
  std::vector<semantic::Index> compiled_outcome_gate_indices;
  semantic::Index source_count{0};
  std::vector<semantic::Index> leaf_source_ids;
  std::vector<semantic::Index> pool_source_ids;
  std::vector<semantic::Index> shared_trigger_indices;
  ExactCompiledTriggerStateTable trigger_state_table;
  std::vector<std::uint8_t> compiled_source_view_relations;
  semantic::Index compiled_source_view_source_count{0};
  std::unordered_map<
      ExactSourceProgramCompileKey,
      semantic::Index,
      ExactSourceProgramCompileKeyHash>
      source_kernel_program_index;
  std::unordered_map<
      ExactSourceProgramCompileKey,
      semantic::Index,
      ExactSourceProgramCompileKeyHash>
      source_base_program_index;
  std::unordered_map<
      ExactConditionedSourceProgramCompileKey,
      semantic::Index,
      ExactConditionedSourceProgramCompileKeyHash>
      source_conditioned_program_index;
};

struct ExactVariantPlan {
  std::vector<runtime::LeafRuntimeDescriptor> leaf_descriptors;
  std::vector<semantic::Index> leaf_trigger_index;
  semantic::Index expr_count{0};
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactCompiledOutcomePlan> compiled_outcomes;
  ExactTerminalNoResponsePlan no_response;
  semantic::Index finite_response_density_root_id{
      semantic::kInvalidIndex};
  semantic::Index finite_response_survival_root_id{
      semantic::kInvalidIndex};
  CompiledMathProgram compiled_math;
  std::vector<semantic::Index> compiled_outcome_gate_indices;
  semantic::Index source_count{0};
  ExactCompiledTriggerStateTable trigger_state_table;
};

inline ExactVariantPlan finalize_exact_variant_plan(
    ExactVariantBuildState &&build) {
  ExactVariantPlan plan;
  plan.leaf_descriptors = std::move(build.program.leaf_descriptors);
  plan.leaf_trigger_index = std::move(build.program.leaf_trigger_index);
  plan.expr_count =
      static_cast<semantic::Index>(build.program.expr_kind.size());
  plan.outcome_index_by_code = std::move(build.outcome_index_by_code);
  plan.compiled_outcomes = std::move(build.compiled_outcomes);
  plan.no_response = std::move(build.no_response);
  plan.finite_response_density_root_id =
      build.finite_response_density_root_id;
  plan.finite_response_survival_root_id =
      build.finite_response_survival_root_id;
  plan.compiled_math = std::move(build.compiled_math);
  plan.compiled_outcome_gate_indices =
      std::move(build.compiled_outcome_gate_indices);
  plan.source_count = build.source_count;
  plan.trigger_state_table = std::move(build.trigger_state_table);
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
      static_cast<std::size_t>(plan.expr_count),
      std::numeric_limits<double>::infinity());
  state.expr_upper_normalizers.assign(
      static_cast<std::size_t>(plan.expr_count),
      0.0);
  return state;
}

inline ExactRelation exact_compiled_source_view_relation(
    const ExactVariantBuildState &plan,
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
  if (kind == semantic::SourceKind::Leaf && index >= 0 &&
      static_cast<std::size_t>(index) < plan.leaf_source_ids.size()) {
    return plan.leaf_source_ids[static_cast<std::size_t>(index)];
  }
  if (kind == semantic::SourceKind::Pool && index >= 0 &&
      static_cast<std::size_t>(index) < plan.pool_source_ids.size()) {
    return plan.pool_source_ids[static_cast<std::size_t>(index)];
  }
  return semantic::kInvalidIndex;
}

} // namespace detail
} // namespace accumulatr::eval
