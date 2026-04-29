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

enum class ExactFactorKind : std::uint8_t {
  AtPdf = 0,
  BeforeCdf = 1,
  AfterSurvival = 2
};

struct ExactIndexSpan {
  semantic::Index offset{0};
  semantic::Index size{0};

  [[nodiscard]] bool empty() const noexcept {
    return size == 0;
  }
};

struct ExactSourceKey {
  semantic::SourceKind kind{semantic::SourceKind::Leaf};
  semantic::Index index{semantic::kInvalidIndex};

  bool operator==(const ExactSourceKey &other) const noexcept {
    return kind == other.kind && index == other.index;
  }
};

struct ExactSourceKeyHash {
  std::size_t operator()(const ExactSourceKey &key) const noexcept {
    return (static_cast<std::size_t>(key.kind) << 32U) ^
           static_cast<std::size_t>(static_cast<std::uint32_t>(key.index));
  }
};

struct ExactSourceConstraint {
  ExactSourceKey key{};
  ExactRelation relation{ExactRelation::Unknown};
};

struct ExactSourceOrderFact {
  semantic::Index before_source_id{semantic::kInvalidIndex};
  semantic::Index after_source_id{semantic::kInvalidIndex};
};

struct ExactGuardUpperBoundFact {
  semantic::Index expr_id{semantic::kInvalidIndex};
  semantic::Index ref_source_id{semantic::kInvalidIndex};
  semantic::Index blocker_source_id{semantic::kInvalidIndex};
};

struct ExactTimedSourceFact {
  semantic::Index source_id{semantic::kInvalidIndex};
  double time{std::numeric_limits<double>::quiet_NaN()};
};

struct ExactTimedExprUpperBound {
  semantic::Index expr_id{semantic::kInvalidIndex};
  double time{std::numeric_limits<double>::quiet_NaN()};
  double normalizer{0.0};
};

struct ExactTimedGuardUpperBound {
  semantic::Index ref_source_id{semantic::kInvalidIndex};
  semantic::Index blocker_source_id{semantic::kInvalidIndex};
  double time{std::numeric_limits<double>::quiet_NaN()};
  double normalizer{0.0};
};

struct ExactRuntimeTermCondition {
  semantic::Index exact_source_id{semantic::kInvalidIndex};
  std::vector<semantic::Index> upper_bound_source_ids;
  std::vector<semantic::Index> lower_bound_source_ids;
  std::vector<semantic::Index> upper_bound_expr_ids;
  std::vector<ExactSourceOrderFact> source_order_facts;
  std::vector<ExactGuardUpperBoundFact> guard_upper_bound_facts;
  semantic::Index observed_condition_id{semantic::kInvalidIndex};
  semantic::Index readiness_condition_id{semantic::kInvalidIndex};
};

inline bool runtime_condition_contains_source(
    const std::vector<semantic::Index> &sources,
    const semantic::Index source_id) {
  return std::find(sources.begin(), sources.end(), source_id) != sources.end();
}

inline void append_runtime_condition_index(
    std::vector<semantic::Index> *sources,
    const semantic::Index index) {
  if (index == semantic::kInvalidIndex) {
    return;
  }
  if (std::find(sources->begin(), sources->end(), index) ==
      sources->end()) {
    sources->push_back(index);
  }
}

inline void append_runtime_source_order_fact(
    std::vector<ExactSourceOrderFact> *facts,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex ||
      before_source_id == after_source_id) {
    return;
  }
  for (const auto &fact : *facts) {
    if (fact.before_source_id == before_source_id &&
        fact.after_source_id == after_source_id) {
      return;
    }
  }
  facts->push_back(ExactSourceOrderFact{before_source_id, after_source_id});
}

inline void append_runtime_guard_upper_bound_fact(
    std::vector<ExactGuardUpperBoundFact> *facts,
    const semantic::Index expr_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (expr_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex ||
      ref_source_id == blocker_source_id) {
    return;
  }
  for (const auto &fact : *facts) {
    if (fact.expr_id == expr_id &&
        fact.ref_source_id == ref_source_id &&
        fact.blocker_source_id == blocker_source_id) {
      return;
    }
  }
  facts->push_back(
      ExactGuardUpperBoundFact{expr_id, ref_source_id, blocker_source_id});
}

inline bool runtime_condition_empty(
    const ExactRuntimeTermCondition &condition) {
  return condition.exact_source_id == semantic::kInvalidIndex &&
         condition.upper_bound_source_ids.empty() &&
         condition.lower_bound_source_ids.empty() &&
         condition.upper_bound_expr_ids.empty() &&
         condition.source_order_facts.empty() &&
         condition.guard_upper_bound_facts.empty();
}

inline bool runtime_condition_equal(const ExactRuntimeTermCondition &lhs,
                                    const ExactRuntimeTermCondition &rhs) {
  return lhs.exact_source_id == rhs.exact_source_id &&
         lhs.upper_bound_source_ids == rhs.upper_bound_source_ids &&
         lhs.lower_bound_source_ids == rhs.lower_bound_source_ids &&
         lhs.upper_bound_expr_ids == rhs.upper_bound_expr_ids &&
         lhs.source_order_facts.size() == rhs.source_order_facts.size() &&
         lhs.guard_upper_bound_facts.size() == rhs.guard_upper_bound_facts.size() &&
         std::equal(
             lhs.source_order_facts.begin(),
             lhs.source_order_facts.end(),
             rhs.source_order_facts.begin(),
             [](const auto &a, const auto &b) {
               return a.before_source_id == b.before_source_id &&
                      a.after_source_id == b.after_source_id;
             }) &&
         std::equal(
             lhs.guard_upper_bound_facts.begin(),
             lhs.guard_upper_bound_facts.end(),
             rhs.guard_upper_bound_facts.begin(),
             [](const auto &a, const auto &b) {
               return a.expr_id == b.expr_id &&
                      a.ref_source_id == b.ref_source_id &&
                      a.blocker_source_id == b.blocker_source_id;
             });
}

inline bool runtime_condition_order_contradiction(
    const ExactRuntimeTermCondition &condition) {
  for (const auto &fact : condition.source_order_facts) {
    if (condition.exact_source_id == fact.before_source_id &&
        runtime_condition_contains_source(
            condition.upper_bound_source_ids,
            fact.after_source_id)) {
      return true;
    }
    if (condition.exact_source_id == fact.after_source_id &&
        runtime_condition_contains_source(
            condition.lower_bound_source_ids,
            fact.before_source_id)) {
      return true;
    }
    if (runtime_condition_contains_source(
            condition.lower_bound_source_ids,
            fact.before_source_id) &&
        runtime_condition_contains_source(
            condition.upper_bound_source_ids,
            fact.after_source_id)) {
      return true;
    }
    for (const auto &guard_fact : condition.guard_upper_bound_facts) {
      if (fact.before_source_id == guard_fact.blocker_source_id &&
          fact.after_source_id == guard_fact.ref_source_id) {
        return true;
      }
    }
  }
  return false;
}

struct ExactRelationTemplate {
  std::vector<semantic::Index> source_ids;
  std::vector<ExactRelation> relations;

  bool empty() const noexcept {
    return source_ids.empty();
  }
};

struct ExactScenarioFactor {
  ExactSourceKey key{};
  ExactFactorKind kind{ExactFactorKind::AtPdf};
};

struct ExactTransitionScenario {
  ExactSourceKey active_key{};
  semantic::Index active_source_id{semantic::kInvalidIndex};
  std::vector<ExactSourceKey> before_keys;
  std::vector<semantic::Index> before_source_ids;
  ExactIndexSpan before_source_span{};
  std::vector<ExactSourceKey> after_keys;
  std::vector<semantic::Index> after_source_ids;
  ExactIndexSpan after_source_span{};
  std::vector<semantic::Index> ready_exprs;
  ExactIndexSpan ready_expr_span{};
  std::vector<semantic::Index> tail_exprs;
  ExactIndexSpan tail_expr_span{};
  std::vector<ExactScenarioFactor> factors;
  std::vector<ExactSourceConstraint> forced;
  std::vector<ExactSourceOrderFact> source_order_facts;
  ExactRelationTemplate relation_template;
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

struct ExactOutcomePlan {
  semantic::Index expr_root{semantic::kInvalidIndex};
  std::vector<ExactTransitionScenario> scenarios;
};

struct ExactCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
  std::vector<semantic::Index> expr_roots;
  semantic::Index singleton_expr_root{semantic::kInvalidIndex};
  int inclusion_sign{1};
  std::vector<ExactTransitionScenario> scenarios;
};

struct ExactCompetitorBlockPlan {
  std::vector<ExactCompetitorSubsetPlan> subsets;
};

struct ExactTargetCompetitorPlan {
  std::vector<ExactCompetitorBlockPlan> blocks;
};

struct ExactRuntimeFactors {
  std::vector<semantic::Index> source_pdf;
  std::vector<semantic::Index> source_cdf;
  std::vector<semantic::Index> source_survival;
  std::vector<semantic::Index> expr_density;
  std::vector<semantic::Index> expr_cdf;
  std::vector<semantic::Index> expr_survival;
};

struct ExactRuntimeProductTerm {
  ExactRuntimeFactors factors;
  ExactRuntimeTermCondition condition;
  ExactRuntimeTermCondition tail_condition;
  ExactRuntimeTermCondition competitor_condition;
  semantic::Index context_group_index{semantic::kInvalidIndex};
  semantic::Index compiled_root_id{semantic::kInvalidIndex};
  bool condition_impossible{false};
};

struct ExactRuntimeCompetitorSubsetMask {
  std::vector<std::vector<std::uint8_t>> affected;
  bool any_affected{false};
};

struct ExactRuntimeConditionGroup {
  ExactRuntimeTermCondition tail_condition;
  ExactRuntimeTermCondition competitor_condition;
  ExactRuntimeCompetitorSubsetMask competitor_subset_mask;
  semantic::Index combined_competitor_condition_id{semantic::kInvalidIndex};
};

struct ExactRuntimeConditionContextPlan {
  semantic::Index condition_id{0};
  semantic::Index after_survival_root_id{semantic::kInvalidIndex};
  bool impossible{false};

  bool operator==(const ExactRuntimeConditionContextPlan &other) const
      noexcept {
    return condition_id == other.condition_id &&
           impossible == other.impossible;
  }
};

struct ExactRuntimeConditionGroupKernel {
  semantic::Index group_index{semantic::kInvalidIndex};
  semantic::Index target_tail_context_index{semantic::kInvalidIndex};
  semantic::Index competitor_context_index{semantic::kInvalidIndex};
  semantic::Index combined_competitor_condition_id{semantic::kInvalidIndex};
  semantic::Index competitor_non_win_plan_index{semantic::kInvalidIndex};
  bool has_tail_condition{false};
  bool has_competitor_condition{false};
};

struct ExactRuntimeReadinessTermKernel {
  semantic::Index term_index{semantic::kInvalidIndex};
  semantic::Index group_kernel_index{semantic::kInvalidIndex};
};

struct ExactRuntimeCompetitorNonWinSubsetPlan {
  semantic::Index block_index{semantic::kInvalidIndex};
  semantic::Index subset_index{semantic::kInvalidIndex};
  int inclusion_sign{1};
  semantic::Index win_cdf_root_id{semantic::kInvalidIndex};
};

struct ExactRuntimeCompetitorNonWinBlockPlan {
  ExactIndexSpan subset_span{};
};

struct ExactRuntimeCompetitorNonWinPlan {
  std::vector<ExactRuntimeCompetitorNonWinBlockPlan> blocks;
  std::vector<ExactRuntimeCompetitorNonWinSubsetPlan> subsets;
  semantic::Index root_id{semantic::kInvalidIndex};
  bool impossible{false};
};

struct ExactRuntimeScenarioExecutionKernel {
  semantic::Index target_scenario_context_index{semantic::kInvalidIndex};
  semantic::Index competitor_scenario_context_index{semantic::kInvalidIndex};
  semantic::Index tail_competitor_context_index{semantic::kInvalidIndex};
  bool has_tail_competitor_condition{false};
  bool has_conditioned_readiness_terms{false};
  bool has_readiness{false};
  std::vector<ExactRuntimeConditionContextPlan> condition_contexts;
  std::vector<ExactRuntimeConditionGroupKernel> condition_groups;
  std::vector<ExactRuntimeReadinessTermKernel> readiness_terms;
  std::vector<ExactRuntimeCompetitorNonWinPlan> competitor_non_win_plans;
  semantic::Index base_competitor_non_win_plan_index{semantic::kInvalidIndex};
  semantic::Index tail_competitor_non_win_plan_index{semantic::kInvalidIndex};
};

struct ExactRuntimeTruthFormula {
  ExactRuntimeFactors product;
  std::vector<ExactRuntimeProductTerm> sum_terms;
  double empty_value{0.0};
  bool sum_of_products{false};
  bool clean_signed{false};
  semantic::Index compiled_root_id{semantic::kInvalidIndex};
};

struct ExactRuntimeScenarioFormula {
  semantic::Index active_source_id{semantic::kInvalidIndex};
  semantic::Index active_observed_condition_id{semantic::kInvalidIndex};
  semantic::Index source_view_id{0};
  ExactRelationTemplate relation_template;
  std::vector<ExactSourceOrderFact> source_order_facts;
  ExactRuntimeTermCondition tail_condition;
  ExactRuntimeTermCondition tail_competitor_condition;
  ExactRuntimeCompetitorSubsetMask tail_competitor_subset_mask;
  std::vector<ExactRuntimeConditionGroup> condition_groups;
  bool has_tail_competitor_condition{false};
  bool has_conditioned_readiness_terms{false};
  bool has_readiness{false};
  ExactRuntimeTruthFormula readiness_cdf;
  ExactRuntimeTruthFormula readiness_density;
  ExactRuntimeTruthFormula after_survival;
  semantic::Index probability_root_id{semantic::kInvalidIndex};
  ExactRuntimeScenarioExecutionKernel execution_kernel;
};

struct ExactConditionedRoot {
  semantic::Index condition_id{semantic::kInvalidIndex};
  semantic::Index root_id{semantic::kInvalidIndex};
};

struct ExactRuntimeCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
  std::vector<semantic::Index> expr_roots;
  int inclusion_sign{1};
  semantic::Index singleton_expr_root{semantic::kInvalidIndex};
  semantic::Index win_cdf_root_id{semantic::kInvalidIndex};
  std::vector<ExactConditionedRoot> conditioned_win_cdf_roots;
  std::vector<ExactRuntimeScenarioFormula> scenarios;
};

struct ExactRuntimeCompetitorBlockPlan {
  std::vector<ExactRuntimeCompetitorSubsetPlan> subsets;
};

struct ExactRuntimeScenarioTransitionPlan {
  semantic::Index probability_root_id{semantic::kInvalidIndex};
  semantic::Index active_source_id{semantic::kInvalidIndex};
  ExactIndexSpan before_source_span{};
  ExactIndexSpan ready_expr_span{};
  bool ranked_supported{false};
};

struct ExactRuntimeSuccessorDistributionPlan {
  std::vector<ExactRuntimeScenarioTransitionPlan> transitions;
  semantic::Index total_probability_root_id{semantic::kInvalidIndex};
};

struct ExactRuntimeOutcomePlan {
  std::vector<ExactRuntimeScenarioFormula> scenarios;
  std::vector<ExactRuntimeCompetitorBlockPlan> competitor_blocks;
  ExactRuntimeSuccessorDistributionPlan successor_distribution;
};

struct ExactRuntimeVariantPlan {
  std::vector<ExactRuntimeOutcomePlan> outcomes;
};

struct ExactNoResponsePlan {
  bool terminal_leaf_survival_product{false};
  std::vector<semantic::Index> source_ids;
};

struct ExactSimpleRacePlan {
  bool terminal_leaf_top1{false};
  std::vector<semantic::Index> source_ids;
  std::vector<semantic::Index> outcome_source_ids;
  std::vector<semantic::Index> source_leaf_indices;
};

enum class ExactProbabilityOpKind : std::uint8_t {
  Constant = 0,
  Top1LeafRaceDensity = 1,
  TerminalNoResponseProbability = 2,
  GenericTransitionDensity = 3,
  GenericTransitionProbability = 4,
  RankedTransitionSequence = 5,
  Integral = 6,
  WeightedTriggerSum = 7,
  Log = 8
};

enum class ExactProbabilityValueKind : std::uint8_t {
  Probability = 0,
  Density = 1,
  Log = 2
};

struct ExactProbabilityOp {
  ExactProbabilityOpKind kind{ExactProbabilityOpKind::Constant};
  ExactProbabilityValueKind value_kind{ExactProbabilityValueKind::Probability};
  semantic::Index outcome_index{semantic::kInvalidIndex};
  semantic::Index target_source_id{semantic::kInvalidIndex};
  semantic::Index time_slot{semantic::kInvalidIndex};
  ExactIndexSpan source_span{};
  ExactIndexSpan trigger_state_span{};
  ExactIndexSpan children{};
  double constant{0.0};
};

struct ExactProbabilityProgram {
  std::vector<ExactProbabilityOp> ops;
  std::vector<semantic::Index> child_ops;
  semantic::Index root{semantic::kInvalidIndex};
  semantic::Index root_child{semantic::kInvalidIndex};
  ExactProbabilityValueKind value_kind{ExactProbabilityValueKind::Probability};
  bool requires_trigger_enumeration{false};

  [[nodiscard]] bool empty() const noexcept {
    return root == semantic::kInvalidIndex || ops.empty();
  }
};

struct ExactProbabilityProgramSet {
  std::vector<ExactProbabilityProgram> programs;
  std::vector<semantic::Index> density_by_outcome;
  std::vector<semantic::Index> finite_probability_by_outcome;
  std::vector<semantic::Index> source_ids;
  semantic::Index no_response_probability{semantic::kInvalidIndex};
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
  std::vector<ExactTargetCompetitorPlan> competitor_plans;
  ExactRuntimeVariantPlan runtime;
  ExactNoResponsePlan no_response;
  ExactSimpleRacePlan simple_race;
  ExactProbabilityProgramSet probability_programs;
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
  std::vector<semantic::Index> scenario_source_ids;
  std::vector<semantic::Index> scenario_expr_ids;
  std::vector<semantic::Index> shared_trigger_indices;
  ExactCompiledTriggerStateTable trigger_state_table;
  std::vector<std::uint8_t> compiled_source_view_relations;
  semantic::Index compiled_source_view_source_count{0};
  semantic::Index expr_union_kernel_cache_slot_count{0};
  bool ranked_supported{true};
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

inline bool runtime_source_relevant_to_competitors(
    const semantic::Index source_id,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  for (const auto &block : runtime_outcome.competitor_blocks) {
    for (const auto &subset : block.subsets) {
      for (const auto expr_root : subset.expr_roots) {
        if (expr_support_contains_source(plan, expr_root, source_id)) {
          return true;
        }
      }
    }
  }
  return false;
}

inline bool runtime_expr_relevant_to_competitors(
    const semantic::Index expr_id,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  for (const auto &block : runtime_outcome.competitor_blocks) {
    for (const auto &subset : block.subsets) {
      for (const auto expr_root : subset.expr_roots) {
        if (expr_supports_overlap(plan, expr_id, expr_root)) {
          return true;
        }
      }
    }
  }
  return false;
}

inline bool runtime_source_order_relevant_to_competitors(
    const ExactSourceOrderFact &fact,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  for (const auto &block : runtime_outcome.competitor_blocks) {
    for (const auto &subset : block.subsets) {
      for (const auto expr_root : subset.expr_roots) {
        if (expr_contains_simple_guard_pair(
                plan,
                expr_root,
                fact.after_source_id,
                fact.before_source_id)) {
          return true;
        }
      }
    }
  }
  return false;
}

inline bool runtime_guard_upper_relevant_to_competitors(
    const ExactGuardUpperBoundFact &fact,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  for (const auto &block : runtime_outcome.competitor_blocks) {
    for (const auto &subset : block.subsets) {
      for (const auto expr_root : subset.expr_roots) {
        if (expr_contains_simple_guard_pair(
                plan,
                expr_root,
                fact.ref_source_id,
                fact.blocker_source_id)) {
          return true;
        }
      }
    }
  }
  return false;
}

inline bool runtime_source_relevant_to_competitor_subset(
    const semantic::Index source_id,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactVariantPlan &plan) {
  for (const auto expr_root : subset.expr_roots) {
    if (expr_support_contains_source(plan, expr_root, source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_expr_relevant_to_competitor_subset(
    const semantic::Index expr_id,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactVariantPlan &plan) {
  for (const auto expr_root : subset.expr_roots) {
    if (expr_supports_overlap(plan, expr_id, expr_root)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_source_order_relevant_to_competitor_subset(
    const ExactSourceOrderFact &fact,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactVariantPlan &plan) {
  for (const auto expr_root : subset.expr_roots) {
    if (expr_contains_simple_guard_pair(
            plan,
            expr_root,
            fact.after_source_id,
            fact.before_source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_guard_upper_relevant_to_competitor_subset(
    const ExactGuardUpperBoundFact &fact,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactVariantPlan &plan) {
  for (const auto expr_root : subset.expr_roots) {
    if (expr_contains_simple_guard_pair(
            plan,
            expr_root,
            fact.ref_source_id,
            fact.blocker_source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_condition_relevant_to_competitor_subset(
    const ExactRuntimeTermCondition &condition,
    const ExactRuntimeCompetitorSubsetPlan &subset,
    const ExactVariantPlan &plan) {
  if (condition.exact_source_id != semantic::kInvalidIndex &&
      runtime_source_relevant_to_competitor_subset(
          condition.exact_source_id, subset, plan)) {
    return true;
  }
  for (const auto source_id : condition.upper_bound_source_ids) {
    if (runtime_source_relevant_to_competitor_subset(source_id, subset, plan)) {
      return true;
    }
  }
  for (const auto source_id : condition.lower_bound_source_ids) {
    if (runtime_source_relevant_to_competitor_subset(source_id, subset, plan)) {
      return true;
    }
  }
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    if (runtime_expr_relevant_to_competitor_subset(expr_id, subset, plan)) {
      return true;
    }
  }
  for (const auto &fact : condition.source_order_facts) {
    if (runtime_source_order_relevant_to_competitor_subset(
            fact, subset, plan)) {
      return true;
    }
  }
  for (const auto &fact : condition.guard_upper_bound_facts) {
    if (runtime_guard_upper_relevant_to_competitor_subset(
            fact, subset, plan)) {
      return true;
    }
  }
  return false;
}

inline ExactRuntimeCompetitorSubsetMask runtime_competitor_subset_mask(
    const ExactRuntimeTermCondition &condition,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  ExactRuntimeCompetitorSubsetMask mask;
  mask.affected.reserve(runtime_outcome.competitor_blocks.size());
  if (runtime_condition_empty(condition)) {
    return mask;
  }
  for (const auto &block : runtime_outcome.competitor_blocks) {
    std::vector<std::uint8_t> block_mask;
    block_mask.reserve(block.subsets.size());
    for (const auto &subset : block.subsets) {
      const bool affected =
          runtime_condition_relevant_to_competitor_subset(
              condition, subset, plan);
      block_mask.push_back(affected ? 1U : 0U);
      mask.any_affected = mask.any_affected || affected;
    }
    mask.affected.push_back(std::move(block_mask));
  }
  if (!mask.any_affected) {
    mask.affected.clear();
  }
  return mask;
}

inline bool runtime_competitor_subset_mask_affected(
    const ExactRuntimeCompetitorSubsetMask *mask,
    const std::size_t block_idx,
    const std::size_t subset_idx) {
  if (mask == nullptr) {
    return false;
  }
  if (!mask->any_affected) {
    return false;
  }
  return mask->affected[block_idx][subset_idx] != 0U;
}

inline ExactRuntimeTermCondition filter_runtime_condition_for_competitors(
    const ExactRuntimeTermCondition &condition,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactVariantPlan &plan) {
  ExactRuntimeTermCondition filtered;
  if (condition.exact_source_id != semantic::kInvalidIndex &&
      runtime_source_relevant_to_competitors(
          condition.exact_source_id, runtime_outcome, plan)) {
    filtered.exact_source_id = condition.exact_source_id;
  }
  for (const auto source_id : condition.upper_bound_source_ids) {
    if (runtime_source_relevant_to_competitors(
            source_id, runtime_outcome, plan)) {
      append_runtime_condition_index(
          &filtered.upper_bound_source_ids, source_id);
    }
  }
  for (const auto source_id : condition.lower_bound_source_ids) {
    if (runtime_source_relevant_to_competitors(
            source_id, runtime_outcome, plan)) {
      append_runtime_condition_index(
          &filtered.lower_bound_source_ids, source_id);
    }
  }
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    if (runtime_expr_relevant_to_competitors(
            expr_id, runtime_outcome, plan)) {
      append_runtime_condition_index(&filtered.upper_bound_expr_ids, expr_id);
    }
  }
  for (const auto &fact : condition.source_order_facts) {
    if (runtime_source_order_relevant_to_competitors(
            fact, runtime_outcome, plan)) {
      append_runtime_source_order_fact(
          &filtered.source_order_facts,
          fact.before_source_id,
          fact.after_source_id);
    }
  }
  for (const auto &fact : condition.guard_upper_bound_facts) {
    if (runtime_guard_upper_relevant_to_competitors(
            fact, runtime_outcome, plan)) {
      append_runtime_guard_upper_bound_fact(
          &filtered.guard_upper_bound_facts,
          fact.expr_id,
          fact.ref_source_id,
          fact.blocker_source_id);
    }
  }
  return filtered;
}

inline bool runtime_factors_source_relevant(const ExactRuntimeFactors &factors,
                                            const ExactVariantPlan &plan,
                                            const semantic::Index source_id) {
  if (runtime_condition_contains_source(factors.source_pdf, source_id) ||
      runtime_condition_contains_source(factors.source_cdf, source_id) ||
      runtime_condition_contains_source(factors.source_survival, source_id)) {
    return true;
  }
  for (const auto expr_id : factors.expr_density) {
    if (expr_support_contains_source(plan, expr_id, source_id)) {
      return true;
    }
  }
  for (const auto expr_id : factors.expr_cdf) {
    if (expr_support_contains_source(plan, expr_id, source_id)) {
      return true;
    }
  }
  for (const auto expr_id : factors.expr_survival) {
    if (expr_support_contains_source(plan, expr_id, source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_truth_source_relevant(const ExactRuntimeTruthFormula &formula,
                                          const ExactVariantPlan &plan,
                                          const semantic::Index source_id) {
  if (runtime_factors_source_relevant(formula.product, plan, source_id)) {
    return true;
  }
  for (const auto &term : formula.sum_terms) {
    if (runtime_factors_source_relevant(term.factors, plan, source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_factors_expr_relevant(const ExactRuntimeFactors &factors,
                                          const ExactVariantPlan &plan,
                                          const semantic::Index expr_id) {
  for (const auto source_id :
       plan.expr_supports[static_cast<std::size_t>(expr_id)]) {
    if (runtime_factors_source_relevant(factors, plan, source_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_truth_expr_relevant(const ExactRuntimeTruthFormula &formula,
                                        const ExactVariantPlan &plan,
                                        const semantic::Index expr_id) {
  if (runtime_factors_expr_relevant(formula.product, plan, expr_id)) {
    return true;
  }
  for (const auto &term : formula.sum_terms) {
    if (runtime_factors_expr_relevant(term.factors, plan, expr_id)) {
      return true;
    }
  }
  return false;
}

inline bool runtime_factors_contain_guard_pair(
    const ExactRuntimeFactors &factors,
    const ExactVariantPlan &plan,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  auto contains_pair = [&](const std::vector<semantic::Index> &expr_ids) {
    for (const auto expr_id : expr_ids) {
      if (expr_contains_simple_guard_pair(
              plan, expr_id, ref_source_id, blocker_source_id)) {
        return true;
      }
    }
    return false;
  };
  return contains_pair(factors.expr_density) ||
         contains_pair(factors.expr_cdf) ||
         contains_pair(factors.expr_survival);
}

inline bool runtime_truth_contains_guard_pair(
    const ExactRuntimeTruthFormula &formula,
    const ExactVariantPlan &plan,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (runtime_factors_contain_guard_pair(
          formula.product, plan, ref_source_id, blocker_source_id)) {
    return true;
  }
  for (const auto &term : formula.sum_terms) {
    if (runtime_factors_contain_guard_pair(
            term.factors, plan, ref_source_id, blocker_source_id)) {
      return true;
    }
  }
  return false;
}

inline ExactRuntimeTermCondition filter_runtime_condition_for_truth(
    const ExactRuntimeTermCondition &condition,
    const ExactRuntimeTruthFormula &formula,
    const ExactVariantPlan &plan) {
  ExactRuntimeTermCondition filtered;
  if (condition.exact_source_id != semantic::kInvalidIndex &&
      runtime_truth_source_relevant(formula, plan, condition.exact_source_id)) {
    filtered.exact_source_id = condition.exact_source_id;
  }
  for (const auto source_id : condition.upper_bound_source_ids) {
    if (runtime_truth_source_relevant(formula, plan, source_id)) {
      append_runtime_condition_index(
          &filtered.upper_bound_source_ids, source_id);
    }
  }
  for (const auto source_id : condition.lower_bound_source_ids) {
    if (runtime_truth_source_relevant(formula, plan, source_id)) {
      append_runtime_condition_index(
          &filtered.lower_bound_source_ids, source_id);
    }
  }
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    if (runtime_truth_expr_relevant(formula, plan, expr_id)) {
      append_runtime_condition_index(&filtered.upper_bound_expr_ids, expr_id);
    }
  }
  for (const auto &fact : condition.source_order_facts) {
    if (runtime_truth_contains_guard_pair(
            formula, plan, fact.after_source_id, fact.before_source_id)) {
      append_runtime_source_order_fact(
          &filtered.source_order_facts,
          fact.before_source_id,
          fact.after_source_id);
    }
  }
  for (const auto &fact : condition.guard_upper_bound_facts) {
    if (runtime_truth_contains_guard_pair(
            formula, plan, fact.ref_source_id, fact.blocker_source_id)) {
      append_runtime_guard_upper_bound_fact(
          &filtered.guard_upper_bound_facts,
          fact.expr_id,
          fact.ref_source_id,
          fact.blocker_source_id);
    }
  }
  return filtered;
}

inline semantic::Index runtime_product_term_source_density_id(
    const ExactRuntimeProductTerm &term) {
  if (term.factors.source_pdf.size() != 1U ||
      !term.factors.expr_density.empty()) {
    return semantic::kInvalidIndex;
  }
  return term.factors.source_pdf.front();
}

inline void append_runtime_factor_source_conditions(
    const ExactRuntimeProductTerm &term,
    ExactRuntimeTermCondition *condition) {
  for (const auto source_id : term.factors.source_cdf) {
    append_runtime_condition_index(
        &condition->upper_bound_source_ids,
        source_id);
  }
  for (const auto source_id : term.factors.source_survival) {
    append_runtime_condition_index(
        &condition->lower_bound_source_ids,
        source_id);
  }
}

inline void append_runtime_factor_expr_conditions(
    const ExactRuntimeProductTerm &term,
    const ExactVariantPlan &plan,
    ExactRuntimeTermCondition *condition) {
  for (const auto expr_id : term.factors.expr_cdf) {
    semantic::Index ref_source_id{semantic::kInvalidIndex};
    semantic::Index blocker_source_id{semantic::kInvalidIndex};
    if (simple_event_guard_sources(
            plan, expr_id, &ref_source_id, &blocker_source_id)) {
      append_runtime_guard_upper_bound_fact(
          &condition->guard_upper_bound_facts,
          expr_id,
          ref_source_id,
          blocker_source_id);
      append_runtime_source_order_fact(
          &condition->source_order_facts,
          ref_source_id,
          blocker_source_id);
      continue;
    }
    append_runtime_condition_index(&condition->upper_bound_expr_ids, expr_id);
  }
}

inline void append_runtime_condition_order_facts(
    ExactRuntimeTermCondition *condition,
    const std::vector<ExactSourceOrderFact> &facts) {
  for (const auto &fact : facts) {
    append_runtime_source_order_fact(
        &condition->source_order_facts,
        fact.before_source_id,
        fact.after_source_id);
  }
}

inline ExactRuntimeTermCondition runtime_product_term_condition(
    const ExactRuntimeProductTerm &term,
    const ExactVariantPlan &plan,
    const std::vector<ExactSourceOrderFact> &source_order_facts = {}) {
  ExactRuntimeTermCondition condition;
  append_runtime_factor_source_conditions(term, &condition);
  append_runtime_factor_expr_conditions(term, plan, &condition);
  append_runtime_condition_order_facts(&condition, source_order_facts);
  const auto source_id = runtime_product_term_source_density_id(term);
  if (source_id != semantic::kInvalidIndex) {
    condition.exact_source_id = source_id;
    return condition;
  }

  if (term.factors.expr_density.size() != 1U ||
      !term.factors.source_pdf.empty()) {
    return condition;
  }
  const auto expr_id = term.factors.expr_density.front();
  semantic::Index ref_source_id{semantic::kInvalidIndex};
  semantic::Index blocker_source_id{semantic::kInvalidIndex};
  if (!simple_event_guard_sources(
          plan, expr_id, &ref_source_id, &blocker_source_id)) {
    append_runtime_condition_index(&condition.upper_bound_expr_ids, expr_id);
    return condition;
  }

  append_runtime_guard_upper_bound_fact(
      &condition.guard_upper_bound_facts,
      expr_id,
      ref_source_id,
      blocker_source_id);
  condition.exact_source_id = ref_source_id;
  append_runtime_condition_index(
      &condition.lower_bound_source_ids,
      blocker_source_id);
  append_runtime_source_order_fact(
      &condition.source_order_facts,
      ref_source_id,
      blocker_source_id);
  return condition;
}

inline ExactRuntimeTermCondition runtime_scenario_tail_condition(
    const ExactRuntimeScenarioFormula &scenario_formula) {
  ExactRuntimeTermCondition condition;
  for (const auto source_id :
       scenario_formula.after_survival.product.source_survival) {
    append_runtime_condition_index(
        &condition.lower_bound_source_ids,
        source_id);
  }
  append_runtime_condition_order_facts(
      &condition,
      scenario_formula.source_order_facts);
  return condition;
}

inline bool has_reason(const std::vector<std::string> &reasons,
                       const std::string &needle) {
  return std::find(reasons.begin(), reasons.end(), needle) != reasons.end();
}

inline std::string source_key_string(const runtime::LoweredExactVariant &variant,
                                     const ExactSourceKey key) {
  if (key.kind == semantic::SourceKind::Leaf &&
      key.index >= 0 &&
      key.index < static_cast<semantic::Index>(variant.leaf_ids.size())) {
    return "leaf '" + variant.leaf_ids[static_cast<std::size_t>(key.index)] + "'";
  }
  if (key.kind == semantic::SourceKind::Pool &&
      key.index >= 0 &&
      key.index < static_cast<semantic::Index>(variant.pool_ids.size())) {
    return "pool '" + variant.pool_ids[static_cast<std::size_t>(key.index)] + "'";
  }
  return "source";
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

inline ExactSourceKey source_key(const semantic::SourceKind kind,
                                 const semantic::Index index) {
  return ExactSourceKey{kind, index};
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

inline semantic::Index source_ordinal(const ExactVariantPlan &plan,
                                      const ExactSourceKey key) {
  return source_ordinal(plan, key.kind, key.index);
}

inline double clean_signed_value(const double value,
                                 const double eps = 1e-15) {
  if (!std::isfinite(value)) {
    return 0.0;
  }
  return std::fabs(value) <= eps ? 0.0 : value;
}

inline bool append_constraint(
    std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> *map,
    const ExactSourceKey key,
    const ExactRelation relation) {
  const auto it = map->find(key);
  if (it == map->end()) {
    map->emplace(key, relation);
    return true;
  }
  return it->second == relation;
}

} // namespace detail
} // namespace accumulatr::eval
