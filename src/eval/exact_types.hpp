#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

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

struct ExactRelationTemplate {
  std::vector<semantic::Index> source_ids;
  std::vector<ExactRelation> relations;

  bool empty() const noexcept {
    return source_ids.empty();
  }

  ExactRelation relation_for(const semantic::Index source_id) const noexcept {
    const auto it =
        std::lower_bound(source_ids.begin(), source_ids.end(), source_id);
    if (it == source_ids.end() || *it != source_id) {
      return ExactRelation::Unknown;
    }
    const auto pos = static_cast<std::size_t>(it - source_ids.begin());
    return relations[pos];
  }
};

class RelationView {
public:
  RelationView() = default;

  RelationView with_overlay(const ExactRelationTemplate *overlay) const noexcept {
    if (overlay == nullptr || overlay->empty()) {
      return *this;
    }
    if (empty()) {
      return RelationView(nullptr, overlay);
    }
    return RelationView(this, overlay);
  }

  ExactRelation relation_for(const semantic::Index source_id) const noexcept {
    if (overlay_ != nullptr) {
      const auto relation = overlay_->relation_for(source_id);
      if (relation != ExactRelation::Unknown) {
        return relation;
      }
    }
    return parent_ != nullptr ? parent_->relation_for(source_id)
                              : ExactRelation::Unknown;
  }

  bool empty() const noexcept {
    return parent_ == nullptr && overlay_ == nullptr;
  }

private:
  RelationView(const RelationView *parent,
               const ExactRelationTemplate *overlay) noexcept
      : parent_(parent), overlay_(overlay) {}

  const RelationView *parent_{nullptr};
  const ExactRelationTemplate *overlay_{nullptr};
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
  std::vector<std::uint8_t> shared_started;
};

struct ExactSequenceState {
  double lower_bound{0.0};
  std::vector<double> exact_times;
};

struct ExactStepBranch {
  double probability{0.0};
  ExactSequenceState next_state;
};

struct ExactStepResult {
  double total_probability{0.0};
  std::vector<ExactStepBranch> branches;
};

struct ExactOutcomePlan {
  semantic::Index expr_root{semantic::kInvalidIndex};
  std::vector<ExactTransitionScenario> scenarios;
};

struct ExactCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
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
};

struct ExactRuntimeTruthFormula {
  ExactRuntimeFactors product;
  std::vector<ExactRuntimeProductTerm> sum_terms;
  double empty_value{0.0};
  bool sum_of_products{false};
  bool clean_signed{false};
  bool requires_scenario{false};
};

struct ExactRuntimeScenarioFormula {
  semantic::Index active_source_id{semantic::kInvalidIndex};
  ExactRelationTemplate relation_template;
  std::vector<ExactSourceOrderFact> source_order_facts;
  bool has_readiness{false};
  ExactRuntimeTruthFormula readiness_cdf;
  ExactRuntimeTruthFormula readiness_density;
  ExactRuntimeTruthFormula after_survival;
};

struct ExactRuntimeCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
  int inclusion_sign{1};
  semantic::Index singleton_expr_root{semantic::kInvalidIndex};
  std::vector<ExactRuntimeScenarioFormula> scenarios;
};

struct ExactRuntimeCompetitorBlockPlan {
  std::vector<ExactRuntimeCompetitorSubsetPlan> subsets;
};

struct ExactRuntimeScenarioSubsetView {
  std::vector<semantic::Index> same_active_scenario_indices;
};

struct ExactRuntimeScenarioBlockView {
  std::vector<ExactRuntimeScenarioSubsetView> subsets;
};

struct ExactRuntimeScenarioCompetitorView {
  std::vector<ExactRuntimeScenarioBlockView> blocks;
};

struct ExactRuntimeOutcomePlan {
  std::vector<ExactRuntimeScenarioFormula> scenarios;
  std::vector<ExactRuntimeCompetitorBlockPlan> competitor_blocks;
  std::vector<ExactRuntimeScenarioCompetitorView> competitor_by_scenario;
};

struct ExactRuntimeVariantPlan {
  std::vector<ExactRuntimeOutcomePlan> outcomes;
};

struct ExactVariantPlan {
  runtime::LoweredExactVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<ExactOutcomePlan> outcomes;
  std::vector<ExactTargetCompetitorPlan> competitor_plans;
  ExactRuntimeVariantPlan runtime;
  std::vector<std::vector<semantic::Index>> leaf_supports;
  std::vector<std::vector<semantic::Index>> pool_supports;
  std::vector<std::vector<semantic::Index>> expr_supports;
  semantic::Index source_count{0};
  std::vector<semantic::Index> leaf_source_ids;
  std::vector<semantic::Index> pool_source_ids;
  std::vector<semantic::Index> scenario_source_ids;
  std::vector<semantic::Index> scenario_expr_ids;
  std::vector<semantic::Index> shared_trigger_indices;
  bool ranked_supported{true};
};

inline ExactSequenceState make_exact_sequence_state(const ExactVariantPlan &plan) {
  ExactSequenceState state;
  state.exact_times.assign(
      static_cast<std::size_t>(plan.source_count),
      std::numeric_limits<double>::quiet_NaN());
  return state;
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
