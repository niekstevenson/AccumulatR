#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "compiled_math.hpp"
#include "../runtime/exact_program.hpp"

namespace accumulatr::eval {
namespace detail {

enum class ExactRelation : std::uint8_t {
  Unknown = 0,
  Before = 1,
  At = 2,
  After = 3
};

struct ExactIndexSpan {
  semantic::Index offset{0};
  semantic::Index size{0};

  [[nodiscard]] bool empty() const noexcept {
    return size == 0;
  }
};

struct ExactSourceOrderFact {
  semantic::Index before_source_id{semantic::kInvalidIndex};
  semantic::Index after_source_id{semantic::kInvalidIndex};
};

struct ExactTimedExprUpperBound {
  semantic::Index expr_id{semantic::kInvalidIndex};
  double time{std::numeric_limits<double>::quiet_NaN()};
  double normalizer{0.0};
};

struct ExactRelationTemplate {
  std::vector<semantic::Index> source_ids;
  std::vector<ExactRelation> relations;

  bool empty() const noexcept {
    return source_ids.empty();
  }
};

struct ExactTriggerState {
  double weight{1.0};
  const std::uint8_t *shared_started{nullptr};
};

struct ExactCompiledTriggerWeightTerm {
  semantic::Index leaf_index{semantic::kInvalidIndex};
  std::uint8_t shared_started{2U};
};

struct ExactCompiledTriggerState {
  double fixed_weight{1.0};
  ExactIndexSpan weight_terms{};
  semantic::Index shared_started_offset{0};
};

struct ExactCompiledTriggerStateTable {
  std::vector<ExactCompiledTriggerState> states;
  std::vector<ExactCompiledTriggerWeightTerm> weight_terms;
  std::vector<std::uint8_t> shared_started_values;
  semantic::Index trigger_count{0};
};

struct ExactSequenceState {
  double lower_bound{0.0};
  std::vector<double> exact_times;
  std::vector<double> upper_bounds;
  std::vector<double> expr_upper_bounds;
  std::vector<double> expr_upper_normalizers;
};

struct ExactRankedFrontierEntry {
  double probability{0.0};
  semantic::Index state_index{semantic::kInvalidIndex};
};

struct ExactStepDistributionView {
  double total_probability{0.0};
  const std::vector<double> *transition_probabilities{nullptr};
};

enum class ExactSymbolicTransitionTimeKind : std::uint8_t {
  SourceRelease = 0
};

struct ExactSymbolicTransitionTimeExpr {
  ExactSymbolicTransitionTimeKind kind{
      ExactSymbolicTransitionTimeKind::SourceRelease};
  semantic::Index source_id{semantic::kInvalidIndex};
};

enum class ExactTransitionGuardKind : std::uint8_t {
  SourceBefore = 0,
  SourceAfter = 1,
  ExprBefore = 2,
  ExprAfter = 3
};

struct ExactTransitionGuard {
  ExactTransitionGuardKind kind{ExactTransitionGuardKind::SourceBefore};
  semantic::Index subject_id{semantic::kInvalidIndex};
};

struct ExactTransitionGuardSet {
  std::vector<ExactTransitionGuard> guards;
  double empty_value{1.0};

  [[nodiscard]] bool empty() const noexcept {
    return guards.empty();
  }
};

struct ExactSymbolicReadinessTimeExpr {
  ExactTransitionGuardSet requirements;

  [[nodiscard]] bool present() const noexcept {
    return !requirements.empty();
  }
};

struct ExactSymbolicTransitionOrderRegion {
  std::vector<ExactSourceOrderFact> source_order_facts;
};

struct ExactSymbolicTransitionTime {
  ExactSymbolicTransitionTimeExpr transition_time_expr;
  ExactSymbolicReadinessTimeExpr readiness_time_expr;
  std::vector<semantic::Index> active_sources;
  ExactTransitionGuardSet guards;
  ExactSymbolicTransitionOrderRegion order_region;
  semantic::Index source_view_id{0};
  ExactRelationTemplate relation_template;
};

inline semantic::Index exact_symbolic_transition_release_source_id(
    const ExactSymbolicTransitionTime &transition) noexcept {
  return transition.transition_time_expr.kind ==
                 ExactSymbolicTransitionTimeKind::SourceRelease
             ? transition.transition_time_expr.source_id
             : semantic::kInvalidIndex;
}

inline bool exact_symbolic_transition_has_readiness(
    const ExactSymbolicTransitionTime &transition) noexcept {
  return transition.readiness_time_expr.present();
}

struct ExactSymbolicTransitionRelation {
  bool competitor_can_strictly_precede{false};
  bool competitor_can_positively_coincide{false};
};

inline ExactSymbolicTransitionRelation exact_symbolic_transition_relation(
    const ExactSymbolicTransitionTime &target,
    const ExactSymbolicTransitionTime &competitor) {
  const auto target_source_id =
      exact_symbolic_transition_release_source_id(target);
  const auto competitor_source_id =
      exact_symbolic_transition_release_source_id(competitor);
  if (target_source_id == semantic::kInvalidIndex ||
      competitor_source_id == semantic::kInvalidIndex) {
    return {};
  }
  if (target_source_id == competitor_source_id) {
    return ExactSymbolicTransitionRelation{false, true};
  }
  return ExactSymbolicTransitionRelation{true, false};
}

struct ExactSymbolicTransitionScenario {
  semantic::Index visible_outcome{semantic::kInvalidIndex};
  ExactSymbolicTransitionTime transition;
  semantic::Index probability_root_id{semantic::kInvalidIndex};
};

struct ExactOutcomePlan {
  semantic::Index expr_root{semantic::kInvalidIndex};
  std::vector<ExactSymbolicTransitionScenario> scenarios;
};

struct ExactRuntimeCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
  std::vector<semantic::Index> expr_roots;
  int inclusion_sign{1};
  std::vector<ExactSymbolicTransitionScenario> scenarios;
};

struct ExactRuntimeCompetitorBlockPlan {
  std::vector<ExactRuntimeCompetitorSubsetPlan> subsets;
};

struct ExactRuntimeScenarioTransitionPlan {
  semantic::Index probability_root_id{semantic::kInvalidIndex};
  semantic::Index release_source_id{semantic::kInvalidIndex};
  std::vector<semantic::Index> readiness_source_ids;
  std::vector<semantic::Index> readiness_expr_ids;
};

struct ExactRuntimeSuccessorDistributionPlan {
  std::vector<ExactRuntimeScenarioTransitionPlan> transitions;
  semantic::Index total_probability_root_id{semantic::kInvalidIndex};
};

struct ExactRuntimeOutcomeCompileContext {
  std::vector<ExactSymbolicTransitionScenario> scenarios;
  std::vector<ExactRuntimeCompetitorBlockPlan> competitor_blocks;
};

struct ExactRuntimeOutcomePlan {
  ExactRuntimeSuccessorDistributionPlan successor_distribution;
};

struct ExactRuntimeVariantPlan {
  std::vector<ExactRuntimeOutcomePlan> outcomes;
};

struct ExactComplexityMetrics {
  semantic::Index symbolic_region_count{0};
  semantic::Index symbolic_cell_count{0};
  semantic::Index max_symbolic_cells_per_region{0};
  semantic::Index negative_symbolic_cell_count{0};
  semantic::Index overlapping_symbolic_cell_pair_count{0};
  semantic::Index expr_relation_atom_count{0};
  semantic::Index compiled_root_count{0};
  semantic::Index compiled_node_count{0};
  semantic::Index integral_node_count{0};
  semantic::Index integral_kernel_count{0};
  semantic::Index source_product_integral_kernel_count{0};
  semantic::Index generic_integral_kernel_count{0};
  semantic::Index max_integral_depth{0};
};

struct ExactOrderRegionTimeOrder {
  semantic::Index before_time_id{semantic::kInvalidIndex};
  semantic::Index after_time_id{semantic::kInvalidIndex};
  bool strict{true};
};

struct ExactOrderRegionSourceTime {
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index time_id{semantic::kInvalidIndex};
  bool inclusive{false};
};

struct ExactOrderRegionExprValueFactor {
  semantic::Index expr_id{semantic::kInvalidIndex};
  semantic::Index time_id{semantic::kInvalidIndex};
  bool before{true};
  bool inclusive{false};
  bool density{false};
};

enum class ExactRegionVarKind : std::uint8_t {
  SourceTime = 0,
  ExprTime = 1,
  Time = 2
};

struct ExactRegionVar {
  ExactRegionVarKind kind{ExactRegionVarKind::Time};
  semantic::Index id{semantic::kInvalidIndex};
};

enum class ExactRegionAtomKind : std::uint8_t {
  SourceExact = 0,
  SourceLower = 1,
  SourceUpper = 2,
  ExprBefore = 3,
  ExprNotBefore = 4,
  ExprDensity = 5,
  TimeOrder = 6,
  OutcomeUnused = 7,
  OutcomeUsed = 8
};

enum class ExactRegionEqualityMass : std::uint8_t {
  MeasureZero = 0,
  PositiveMass = 1
};

enum class ExactRegionEqualityOrigin : std::uint8_t {
  ContinuousBoundary = 0,
  ModelTie = 1,
  SharedLatentIdentity = 2
};

struct ExactRegionAtom {
  ExactRegionAtomKind kind{ExactRegionAtomKind::SourceLower};
  ExactRegionVar lhs{};
  ExactRegionVar rhs{};
  std::vector<semantic::Index> outcome_indices;
  bool inclusive{false};
  bool strict{true};
};

struct ExactRegionEquality {
  semantic::Index lhs_time_id{semantic::kInvalidIndex};
  semantic::Index rhs_time_id{semantic::kInvalidIndex};
  ExactRegionEqualityMass mass{ExactRegionEqualityMass::MeasureZero};
  ExactRegionEqualityOrigin origin{
      ExactRegionEqualityOrigin::ContinuousBoundary};
};

struct ExactRegionCell {
  double sign{1.0};
  std::vector<ExactRegionAtom> atoms;
  std::vector<ExactRegionEquality> equalities;
  bool impossible{false};
};

struct ExactOrderRegionExpr {
  std::vector<ExactRegionCell> terms;
};

struct ExactOrderRegionBuilder {
  semantic::Index next_time_id{
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) + 1U};
};

struct ExactOrderRegionTimeClosure {
  std::vector<semantic::Index> time_ids;
  std::vector<std::uint8_t> relation;
  std::vector<semantic::Index> representative;
  bool impossible{false};
};

struct ExactExprUnionSubset {
  ExactIndexSpan children;
  ExactIndexSpan child_positions;
  int sign{1};
};

struct ExactExprKernel {
  semantic::ExprKind kind{semantic::ExprKind::Impossible};
  ExactIndexSpan children;
  ExactIndexSpan union_subset_span;
  ExactIndexSpan union_absorption_candidate_span;
  ExactIndexSpan guard_ref_child_span;
  semantic::Index event_source_id{semantic::kInvalidIndex};
  semantic::Index guard_ref_expr_id{semantic::kInvalidIndex};
  semantic::Index guard_blocker_expr_id{semantic::kInvalidIndex};
  semantic::Index guard_ref_source_id{semantic::kInvalidIndex};
  semantic::Index guard_blocker_source_id{semantic::kInvalidIndex};
  semantic::ExprKind guard_ref_kind{semantic::ExprKind::Impossible};
  semantic::ExprKind guard_blocker_kind{semantic::ExprKind::Impossible};
  semantic::Index union_kernel_cache_slot{semantic::kInvalidIndex};
  bool children_overlap{false};
  bool simple_event_guard{false};
  bool has_unless{false};
};

struct ExactSourceKernel {
  CompiledSourceChannelKernelKind kind{
      CompiledSourceChannelKernelKind::Invalid};
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index leaf_index{semantic::kInvalidIndex};
  semantic::Index pool_index{semantic::kInvalidIndex};
  semantic::Index onset_source_id{semantic::kInvalidIndex};
  semantic::Index pool_member_offset{0};
  semantic::Index pool_member_count{0};
  semantic::Index pool_k{0};
};

struct ExactVariantPlan {
  runtime::LoweredExactVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactOutcomePlan> outcomes;
  ExactRuntimeVariantPlan runtime;
  ExactComplexityMetrics complexity;
  std::vector<int> leaf_row_offsets;
  CompiledMathProgram compiled_math;
  std::vector<ExactRelationTemplate> compiled_source_views;
  std::vector<ExactExprKernel> expr_kernels;
  std::vector<ExactSourceKernel> source_kernels;
  std::vector<ExactExprUnionSubset> expr_union_subsets;
  std::vector<semantic::Index> expr_union_subset_children;
  std::vector<semantic::Index> expr_union_subset_child_positions;
  std::vector<semantic::Index> expr_union_absorption_sources;
  std::vector<std::vector<semantic::Index>> leaf_supports;
  std::vector<std::vector<semantic::Index>> pool_supports;
  std::vector<std::vector<semantic::Index>> expr_supports;
  std::vector<std::uint8_t> sequence_expr_upper_bound_used;
  std::vector<semantic::Index> sequence_expr_cdf_roots;
  std::vector<semantic::Index> compiled_outcome_gate_indices;
  semantic::Index source_count{0};
  std::vector<semantic::Index> leaf_source_ids;
  std::vector<semantic::Index> pool_source_ids;
  std::vector<semantic::Index> shared_trigger_indices;
  ExactCompiledTriggerStateTable trigger_state_table;
  std::vector<std::uint8_t> compiled_source_view_relations;
  semantic::Index compiled_source_view_source_count{0};
  semantic::Index expr_union_kernel_cache_slot_count{0};
};

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

inline ExactRelation exact_compiled_source_view_relation(
    const ExactVariantPlan &plan,
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

inline std::vector<semantic::Index> merge_sorted_support(
    std::vector<semantic::Index> merged,
    const std::vector<semantic::Index> &rhs) {
  for (const auto idx : rhs) {
    if (std::find(merged.begin(), merged.end(), idx) == merged.end()) {
      merged.push_back(idx);
    }
  }
  std::sort(merged.begin(), merged.end());
  return merged;
}

inline bool supports_overlap(const std::vector<semantic::Index> &lhs,
                             const std::vector<semantic::Index> &rhs) {
  std::vector<semantic::Index> overlap;
  std::set_intersection(lhs.begin(),
                        lhs.end(),
                        rhs.begin(),
                        rhs.end(),
                        std::back_inserter(overlap));
  return !overlap.empty();
}

inline bool support_contains_source(const std::vector<semantic::Index> &support,
                                    const semantic::Index source_id) {
  return source_id != semantic::kInvalidIndex &&
         std::binary_search(support.begin(), support.end(), source_id);
}

inline bool expr_support_contains_source(const ExactVariantPlan &plan,
                                         const semantic::Index expr_idx,
                                         const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  return support_contains_source(
      plan.expr_supports[static_cast<std::size_t>(expr_idx)], source_id);
}

inline bool expr_supports_overlap(const ExactVariantPlan &plan,
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

inline bool simple_event_guard_sources(const ExactVariantPlan &plan,
                                       const semantic::Index expr_id,
                                       semantic::Index *ref_source_id,
                                       semantic::Index *blocker_source_id) {
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan.expr_kernels.size()) {
    return false;
  }
  const auto &kernel = plan.expr_kernels[static_cast<std::size_t>(expr_id)];
  if (!kernel.simple_event_guard || kernel.has_unless) {
    return false;
  }
  if (ref_source_id != nullptr) {
    *ref_source_id = kernel.guard_ref_source_id;
  }
  if (blocker_source_id != nullptr) {
    *blocker_source_id = kernel.guard_blocker_source_id;
  }
  return true;
}

inline bool expr_contains_simple_guard_pair(const ExactVariantPlan &plan,
                                            const semantic::Index expr_idx,
                                            const semantic::Index ref_source_id,
                                            const semantic::Index blocker_source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Guard) {
    semantic::Index guard_ref_source_id{semantic::kInvalidIndex};
    semantic::Index guard_blocker_source_id{semantic::kInvalidIndex};
    if (simple_event_guard_sources(
            plan, expr_idx, &guard_ref_source_id, &guard_blocker_source_id) &&
        guard_ref_source_id == ref_source_id &&
        guard_blocker_source_id == blocker_source_id) {
      return true;
    }
    if (expr_contains_simple_guard_pair(
            plan, program.expr_ref_child[pos], ref_source_id, blocker_source_id) ||
        expr_contains_simple_guard_pair(
            plan,
            program.expr_blocker_child[pos],
            ref_source_id,
            blocker_source_id)) {
      return true;
    }
  }
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    if (expr_contains_simple_guard_pair(
            plan,
            program.expr_args[static_cast<std::size_t>(i)],
            ref_source_id,
            blocker_source_id)) {
      return true;
    }
  }
  return false;
}

inline bool has_reason(const std::vector<std::string> &reasons,
                       const std::string &needle) {
  return std::find(reasons.begin(), reasons.end(), needle) != reasons.end();
}

class ExactSupportBuilder {
public:
  explicit ExactSupportBuilder(const runtime::LoweredExactVariant &lowered)
      : program_(lowered.program),
        leaf_ready_(static_cast<std::size_t>(program_.layout.n_leaves), 0U),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U),
        expr_ready_(program_.expr_kind.size(), 0U),
        leaf_supports_(static_cast<std::size_t>(program_.layout.n_leaves)),
        pool_supports_(static_cast<std::size_t>(program_.layout.n_pools)),
        expr_supports_(program_.expr_kind.size()) {}

  std::vector<std::vector<semantic::Index>> build_leaf_supports() {
    for (semantic::Index i = 0; i < program_.layout.n_leaves; ++i) {
      leaf_support(i);
    }
    return leaf_supports_;
  }

  std::vector<std::vector<semantic::Index>> build_pool_supports() {
    for (semantic::Index i = 0; i < program_.layout.n_pools; ++i) {
      pool_support(i);
    }
    return pool_supports_;
  }

  std::vector<std::vector<semantic::Index>> build_expr_supports() {
    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(program_.expr_kind.size());
         ++i) {
      expr_support(i);
    }
    return expr_supports_;
  }

private:
  const runtime::ExactProgram &program_;
  std::vector<std::uint8_t> leaf_ready_;
  std::vector<std::uint8_t> pool_ready_;
  std::vector<std::uint8_t> expr_ready_;
  std::vector<std::vector<semantic::Index>> leaf_supports_;
  std::vector<std::vector<semantic::Index>> pool_supports_;
  std::vector<std::vector<semantic::Index>> expr_supports_;

  const std::vector<semantic::Index> &leaf_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (leaf_ready_[pos] == 2U) {
      return leaf_supports_[pos];
    }
    if (leaf_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic leaf dependency in exact support builder");
    }
    leaf_ready_[pos] = 1U;
    std::vector<semantic::Index> support{idx};
    const auto onset_kind =
        static_cast<semantic::OnsetKind>(program_.onset_kind[pos]);
    if (onset_kind != semantic::OnsetKind::Absolute) {
      support = merge_sorted_support(
          std::move(support),
          source_support(
              static_cast<semantic::SourceKind>(program_.onset_source_kind[pos]),
              program_.onset_source_index[pos]));
    }
    leaf_supports_[pos] = std::move(support);
    leaf_ready_[pos] = 2U;
    return leaf_supports_[pos];
  }

  const std::vector<semantic::Index> &pool_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (pool_ready_[pos] == 2U) {
      return pool_supports_[pos];
    }
    if (pool_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic pool dependency in exact support builder");
    }
    pool_ready_[pos] = 1U;
    std::vector<semantic::Index> support;
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      support = merge_sorted_support(
          std::move(support),
          source_support(
              static_cast<semantic::SourceKind>(
                  program_.pool_member_kind[static_cast<std::size_t>(i)]),
              program_.pool_member_indices[static_cast<std::size_t>(i)]));
    }
    pool_supports_[pos] = std::move(support);
    pool_ready_[pos] = 2U;
    return pool_supports_[pos];
  }

  const std::vector<semantic::Index> &expr_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (expr_ready_[pos] == 2U) {
      return expr_supports_[pos];
    }
    if (expr_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic expression dependency in exact support builder");
    }
    expr_ready_[pos] = 1U;
    std::vector<semantic::Index> support;
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    switch (kind) {
    case semantic::ExprKind::Event:
      support = source_support(
          static_cast<semantic::SourceKind>(program_.expr_source_kind[pos]),
          program_.expr_source_index[pos]);
      break;
    case semantic::ExprKind::And:
    case semantic::ExprKind::Or:
    case semantic::ExprKind::Not: {
      const auto begin = program_.expr_arg_offsets[pos];
      const auto end = program_.expr_arg_offsets[pos + 1U];
      for (semantic::Index i = begin; i < end; ++i) {
        support = merge_sorted_support(
            std::move(support),
            expr_support(program_.expr_args[static_cast<std::size_t>(i)]));
      }
      break;
    }
    case semantic::ExprKind::Guard:
      support = merge_sorted_support(
          std::move(support),
          expr_support(program_.expr_ref_child[pos]));
      support = merge_sorted_support(
          std::move(support),
          expr_support(program_.expr_blocker_child[pos]));
      for (semantic::Index i = program_.expr_arg_offsets[pos];
           i < program_.expr_arg_offsets[pos + 1U];
           ++i) {
        support = merge_sorted_support(
            std::move(support),
            expr_support(program_.expr_args[static_cast<std::size_t>(i)]));
      }
      break;
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      break;
    }
    expr_supports_[pos] = std::move(support);
    expr_ready_[pos] = 2U;
    return expr_supports_[pos];
  }

  std::vector<semantic::Index> source_support(const semantic::SourceKind kind,
                                              const semantic::Index index) {
    switch (kind) {
    case semantic::SourceKind::Leaf:
      return leaf_support(index);
    case semantic::SourceKind::Pool:
      return pool_support(index);
    case semantic::SourceKind::Special:
      break;
    }
    return {};
  }
};

inline semantic::Index child_event_source_index(const runtime::ExactProgram &program,
                                                const semantic::Index expr_idx) {
  return program.expr_source_index[static_cast<std::size_t>(expr_idx)];
}

inline semantic::SourceKind child_event_source_kind(
    const runtime::ExactProgram &program,
    const semantic::Index expr_idx) {
  return static_cast<semantic::SourceKind>(
      program.expr_source_kind[static_cast<std::size_t>(expr_idx)]);
}

inline void validate_exact_expr(const runtime::LoweredExactVariant &lowered,
                                const semantic::Index expr_idx) {
  const auto &program = lowered.program;
  const auto kind = static_cast<semantic::ExprKind>(
      program.expr_kind[static_cast<std::size_t>(expr_idx)]);
  if (kind == semantic::ExprKind::Event) {
    return;
  }
  if (kind == semantic::ExprKind::Impossible || kind == semantic::ExprKind::TrueExpr) {
    return;
  }
  if (kind == semantic::ExprKind::Not) {
    validate_exact_expr(
        lowered,
        program.expr_args[static_cast<std::size_t>(
            program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)])]);
    return;
  }
  if (kind == semantic::ExprKind::And || kind == semantic::ExprKind::Or) {
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      validate_exact_expr(lowered, program.expr_args[static_cast<std::size_t>(i)]);
    }
    return;
  }
  if (kind == semantic::ExprKind::Guard) {
    validate_exact_expr(
        lowered,
        program.expr_ref_child[static_cast<std::size_t>(expr_idx)]);
    validate_exact_expr(
        lowered,
        program.expr_blocker_child[static_cast<std::size_t>(expr_idx)]);
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      validate_exact_expr(lowered, program.expr_args[static_cast<std::size_t>(i)]);
    }
  }
}

inline semantic::Index source_ordinal(const ExactVariantPlan &plan,
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

inline double clean_signed_value(const double value,
                                 const double eps = 1e-15) {
  if (!std::isfinite(value)) {
    return 0.0;
  }
  return std::fabs(value) <= eps ? 0.0 : value;
}

} // namespace detail
} // namespace accumulatr::eval
