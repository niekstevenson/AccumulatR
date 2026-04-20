#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "direct_kernel.hpp"
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

struct ExactScenarioFactor {
  ExactSourceKey key{};
  ExactFactorKind kind{ExactFactorKind::AtPdf};
};

struct ExactTransitionScenario {
  ExactSourceKey active_key{};
  std::vector<ExactSourceKey> before_keys;
  std::vector<ExactSourceKey> after_keys;
  std::vector<semantic::Index> ready_exprs;
  std::vector<semantic::Index> tail_exprs;
  std::vector<ExactScenarioFactor> factors;
  std::vector<ExactSourceConstraint> forced;
};

struct ExactTriggerState {
  double weight{1.0};
  std::vector<std::uint8_t> shared_started;
};

struct ExactSequenceState {
  double lower_bound{0.0};
  std::unordered_map<ExactSourceKey, double, ExactSourceKeyHash> exact_times;
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
  std::vector<semantic::Index> support;
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
  std::vector<semantic::Index> outcome_indices;
  std::vector<ExactCompetitorSubsetPlan> subsets;
};

struct ExactTargetCompetitorPlan {
  std::vector<ExactCompetitorBlockPlan> blocks;
};

struct ExactVariantPlan {
  runtime::LoweredExactVariant lowered;
  std::unordered_map<std::string, semantic::Index> outcome_index;
  std::vector<ExactOutcomePlan> outcomes;
  std::vector<ExactTargetCompetitorPlan> competitor_plans;
  std::vector<std::vector<semantic::Index>> leaf_supports;
  std::vector<std::vector<semantic::Index>> pool_supports;
  std::vector<std::vector<semantic::Index>> expr_supports;
  std::vector<semantic::Index> shared_trigger_indices;
  bool ranked_supported{true};
};

struct ExactObservedRank {
  std::string outcome_label;
  double rt{NA_REAL};
};

struct ExactTrialObservation {
  semantic::Index variant_index{semantic::kInvalidIndex};
  std::string component_id;
  std::vector<ExactObservedRank> ranks;
  bool valid{true};
};

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
      : lowered_(lowered),
        program_(lowered.program),
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
  const runtime::LoweredExactVariant &lowered_;
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
    const auto onset_kind = static_cast<semantic::OnsetKind>(
        program_.onset_kind[pos]);
    if (onset_kind != semantic::OnsetKind::Absolute) {
      support = merge_sorted_support(
          std::move(support),
          source_support(static_cast<semantic::SourceKind>(
                             program_.onset_source_kind[pos]),
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
          source_support(static_cast<semantic::SourceKind>(
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
    validate_exact_expr(lowered,
                        program.expr_ref_child[static_cast<std::size_t>(expr_idx)]);
    validate_exact_expr(lowered,
                        program.expr_blocker_child[static_cast<std::size_t>(expr_idx)]);
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      validate_exact_expr(lowered, program.expr_args[static_cast<std::size_t>(i)]);
    }
    return;
  }
}

inline ExactSourceKey source_key(const semantic::SourceKind kind,
                                 const semantic::Index index) {
  return ExactSourceKey{kind, index};
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
  if (it->second == relation) {
    return true;
  }
  return false;
}

inline bool append_factor(ExactTransitionScenario *scenario,
                          const ExactScenarioFactor factor,
                          const ExactVariantPlan &plan) {
  const auto factor_support =
      factor.key.kind == semantic::SourceKind::Leaf
          ? plan.leaf_supports[static_cast<std::size_t>(factor.key.index)]
          : plan.pool_supports[static_cast<std::size_t>(factor.key.index)];
  for (const auto &existing : scenario->factors) {
    if (existing.key == factor.key) {
      continue;
    }
    const auto &existing_support =
        existing.key.kind == semantic::SourceKind::Leaf
            ? plan.leaf_supports[static_cast<std::size_t>(existing.key.index)]
            : plan.pool_supports[static_cast<std::size_t>(existing.key.index)];
    if (supports_overlap(existing_support, factor_support)) {
      return false;
    }
  }
  scenario->factors.push_back(factor);
  return true;
}

inline bool append_key(std::vector<ExactSourceKey> *keys,
                       const ExactSourceKey key) {
  const auto it = std::find_if(
      keys->begin(),
      keys->end(),
      [&](const ExactSourceKey &existing) { return existing == key; });
  if (it == keys->end()) {
    keys->push_back(key);
  }
  return true;
}

inline bool expr_is_runtime_const(const runtime::ExactProgram &program,
                                  const semantic::Index expr_idx) {
  const auto kind = static_cast<semantic::ExprKind>(
      program.expr_kind[static_cast<std::size_t>(expr_idx)]);
  return kind == semantic::ExprKind::Impossible ||
         kind == semantic::ExprKind::TrueExpr;
}

inline bool expr_is_simple_event(const runtime::ExactProgram &program,
                                 const semantic::Index expr_idx) {
  return static_cast<semantic::ExprKind>(
             program.expr_kind[static_cast<std::size_t>(expr_idx)]) ==
         semantic::ExprKind::Event;
}

inline void append_ready_requirement(ExactTransitionScenario *scenario,
                                     const runtime::ExactProgram &program,
                                     const semantic::Index expr_idx) {
  if (expr_is_runtime_const(program, expr_idx)) {
    return;
  }
  if (expr_is_simple_event(program, expr_idx)) {
    append_key(&scenario->before_keys,
               ExactSourceKey{child_event_source_kind(program, expr_idx),
                              child_event_source_index(program, expr_idx)});
    return;
  }
  scenario->ready_exprs.push_back(expr_idx);
}

inline void append_tail_requirement(ExactTransitionScenario *scenario,
                                    const runtime::ExactProgram &program,
                                    const semantic::Index expr_idx) {
  if (expr_is_runtime_const(program, expr_idx)) {
    return;
  }
  if (expr_is_simple_event(program, expr_idx)) {
    append_key(&scenario->after_keys,
               ExactSourceKey{child_event_source_kind(program, expr_idx),
                              child_event_source_index(program, expr_idx)});
    return;
  }
  scenario->tail_exprs.push_back(expr_idx);
}

inline bool append_source_truth_constraints(
    const ExactVariantPlan &plan,
    const ExactSourceKey key,
    const ExactRelation relation,
    std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> *map) {
  if (!append_constraint(map, key, relation)) {
    return false;
  }
  if (key.kind != semantic::SourceKind::Leaf) {
    return true;
  }
  if (relation == ExactRelation::After) {
    return true;
  }
  const auto pos = static_cast<std::size_t>(key.index);
  const auto onset_kind = static_cast<semantic::OnsetKind>(
      plan.lowered.program.onset_kind[pos]);
  if (onset_kind == semantic::OnsetKind::Absolute) {
    return true;
  }
  return append_source_truth_constraints(
      plan,
      ExactSourceKey{
          static_cast<semantic::SourceKind>(
              plan.lowered.program.onset_source_kind[pos]),
          plan.lowered.program.onset_source_index[pos]},
      ExactRelation::Before,
      map);
}

inline std::vector<ExactSourceConstraint> constraints_from_map(
    const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> &map) {
  std::vector<ExactSourceConstraint> out;
  out.reserve(map.size());
  for (const auto &[key, relation] : map) {
    out.push_back(ExactSourceConstraint{key, relation});
  }
  return out;
}

inline bool expand_member_subsets(
    const std::vector<ExactSourceKey> &members,
    const int k,
    const int active_idx,
    const int member_idx,
    const int before_needed,
    std::vector<int> *before_indices,
    std::vector<std::vector<int>> *subsets) {
  if (before_needed < 0) {
    return false;
  }
  if (member_idx == static_cast<int>(members.size())) {
    if (before_needed == 0) {
      subsets->push_back(*before_indices);
      return true;
    }
    return false;
  }
  if (member_idx == active_idx) {
    return expand_member_subsets(
        members, k, active_idx, member_idx + 1, before_needed, before_indices, subsets);
  }
  if (before_needed > 0) {
    before_indices->push_back(member_idx);
    expand_member_subsets(
        members, k, active_idx, member_idx + 1, before_needed - 1, before_indices, subsets);
    before_indices->pop_back();
  }
  expand_member_subsets(
      members, k, active_idx, member_idx + 1, before_needed, before_indices, subsets);
  return true;
}

inline std::vector<ExactTransitionScenario> build_source_transition_scenarios(
    const ExactVariantPlan &plan,
    const ExactSourceKey key) {
  std::vector<ExactTransitionScenario> out;
  if (key.kind == semantic::SourceKind::Leaf) {
    ExactTransitionScenario scenario;
    scenario.active_key = key;
    if (!append_factor(&scenario, ExactScenarioFactor{key, ExactFactorKind::AtPdf}, plan)) {
      return out;
    }
    std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
    if (!append_source_truth_constraints(plan, key, ExactRelation::At, &forced)) {
      return out;
    }
    scenario.forced = constraints_from_map(forced);
    out.push_back(std::move(scenario));
    return out;
  }

  const auto pool_idx = static_cast<std::size_t>(key.index);
  const auto begin = plan.lowered.program.pool_member_offsets[pool_idx];
  const auto end = plan.lowered.program.pool_member_offsets[pool_idx + 1U];
  const auto k = plan.lowered.program.pool_k[pool_idx];
  std::vector<ExactSourceKey> members;
  members.reserve(static_cast<std::size_t>(end - begin));
  for (semantic::Index i = begin; i < end; ++i) {
    members.push_back(ExactSourceKey{
        static_cast<semantic::SourceKind>(
            plan.lowered.program.pool_member_kind[static_cast<std::size_t>(i)]),
        plan.lowered.program.pool_member_indices[static_cast<std::size_t>(i)]});
  }
  if (k < 1 || k > static_cast<int>(members.size())) {
    return out;
  }

  for (int active_idx = 0; active_idx < static_cast<int>(members.size()); ++active_idx) {
    std::vector<std::vector<int>> before_subsets;
    std::vector<int> current;
    expand_member_subsets(
        members, k, active_idx, 0, k - 1, &current, &before_subsets);
    for (const auto &subset : before_subsets) {
      ExactTransitionScenario scenario;
      scenario.active_key = members[static_cast<std::size_t>(active_idx)];
      std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
      if (!append_factor(&scenario,
                         ExactScenarioFactor{members[static_cast<std::size_t>(active_idx)],
                                             ExactFactorKind::AtPdf},
                         plan)) {
        continue;
      }
      if (!append_source_truth_constraints(
              plan, members[static_cast<std::size_t>(active_idx)], ExactRelation::At, &forced)) {
        continue;
      }
      bool ok = true;
      for (int member_idx = 0; member_idx < static_cast<int>(members.size()); ++member_idx) {
        if (member_idx == active_idx) {
          continue;
        }
        const bool is_before =
            std::find(subset.begin(), subset.end(), member_idx) != subset.end();
        const auto relation =
            is_before ? ExactRelation::Before : ExactRelation::After;
        const auto factor_kind =
            is_before ? ExactFactorKind::BeforeCdf : ExactFactorKind::AfterSurvival;
        append_key(
            is_before ? &scenario.before_keys : &scenario.after_keys,
            members[static_cast<std::size_t>(member_idx)]);
        if (!append_factor(&scenario,
                           ExactScenarioFactor{
                               members[static_cast<std::size_t>(member_idx)],
                               factor_kind},
                           plan)) {
          ok = false;
          break;
        }
        if (!append_source_truth_constraints(
                plan,
                members[static_cast<std::size_t>(member_idx)],
                relation,
                &forced)) {
          ok = false;
          break;
        }
      }
      if (!ok) {
        continue;
      }
      scenario.forced = constraints_from_map(forced);
      out.push_back(std::move(scenario));
    }
  }
  return out;
}

inline std::vector<ExactTransitionScenario> build_expr_transition_scenarios(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx);

inline std::vector<ExactTransitionScenario> build_logical_transition_scenarios(
    const ExactVariantPlan &plan,
    const std::vector<semantic::Index> &children,
    const semantic::ExprKind kind) {
  const auto &program = plan.lowered.program;
  std::vector<semantic::Index> normalized_children;
  normalized_children.reserve(children.size());
  for (const auto child : children) {
    const auto child_kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(child)]);
    if (child_kind == semantic::ExprKind::Impossible) {
      return {};
    }
    if (child_kind == semantic::ExprKind::TrueExpr) {
      continue;
    }
    normalized_children.push_back(child);
  }
  std::vector<std::vector<ExactTransitionScenario>> child_scenarios(
      normalized_children.size());
  for (std::size_t i = 0; i < normalized_children.size(); ++i) {
    child_scenarios[i] =
        build_expr_transition_scenarios(plan, normalized_children[i]);
    if (kind == semantic::ExprKind::Or && child_scenarios[i].empty()) {
      throw std::runtime_error(
          "exact kernel does not support logical-not branches inside first_of()/or outcomes");
    }
  }
  std::vector<ExactTransitionScenario> out;
  for (std::size_t active = 0; active < normalized_children.size(); ++active) {
    if (child_scenarios[active].empty()) {
      continue;
    }
    for (const auto &child_scenario : child_scenarios[active]) {
      ExactTransitionScenario scenario = child_scenario;
      std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
      for (const auto &constraint : scenario.forced) {
        if (!append_constraint(&forced, constraint.key, constraint.relation)) {
          forced.clear();
          break;
        }
      }
      if (forced.empty() && !scenario.forced.empty()) {
        continue;
      }
      bool ok = true;
      for (std::size_t other = 0; other < normalized_children.size(); ++other) {
        if (other == active) {
          continue;
        }
        const auto child = normalized_children[other];
        if (kind == semantic::ExprKind::And) {
          append_ready_requirement(&scenario, program, child);
          if (expr_is_simple_event(program, child)) {
            const auto child_key =
                ExactSourceKey{child_event_source_kind(program, child),
                               child_event_source_index(program, child)};
            ok = append_factor(&scenario,
                               ExactScenarioFactor{child_key, ExactFactorKind::BeforeCdf},
                               plan) &&
                 append_source_truth_constraints(
                     plan, child_key, ExactRelation::Before, &forced);
          }
        } else {
          append_tail_requirement(&scenario, program, child);
          if (expr_is_simple_event(program, child)) {
            const auto child_key =
                ExactSourceKey{child_event_source_kind(program, child),
                               child_event_source_index(program, child)};
            ok = append_factor(
                     &scenario,
                     ExactScenarioFactor{child_key, ExactFactorKind::AfterSurvival},
                     plan) &&
                 append_source_truth_constraints(
                     plan, child_key, ExactRelation::After, &forced);
          }
        }
        if (!ok) {
          break;
        }
      }
      if (!ok) {
        continue;
      }
      scenario.forced = constraints_from_map(forced);
      out.push_back(std::move(scenario));
    }
  }
  return out;
}

inline std::vector<ExactTransitionScenario> build_expr_conjunction_scenarios(
    const ExactVariantPlan &plan,
    const std::vector<semantic::Index> &children) {
  if (children.empty()) {
    return {};
  }
  const auto &program = plan.lowered.program;
  std::vector<semantic::Index> flattened;
  flattened.reserve(children.size());
  std::function<void(semantic::Index)> collect = [&](const semantic::Index expr_idx) {
    const auto kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(expr_idx)]);
    if (kind == semantic::ExprKind::And) {
      const auto begin =
          program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
      const auto end =
          program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
      for (semantic::Index i = begin; i < end; ++i) {
        collect(program.expr_args[static_cast<std::size_t>(i)]);
      }
      return;
    }
    if (kind == semantic::ExprKind::TrueExpr) {
      return;
    }
    flattened.push_back(expr_idx);
  };
  for (const auto child : children) {
    collect(child);
  }
  if (flattened.empty()) {
    return {};
  }
  std::vector<semantic::Index> normalized;
  normalized.reserve(flattened.size());
  auto has_equivalent = [&](const semantic::Index expr_idx) {
    const auto kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(expr_idx)]);
    for (const auto existing : normalized) {
      const auto existing_kind = static_cast<semantic::ExprKind>(
          program.expr_kind[static_cast<std::size_t>(existing)]);
      if (kind != existing_kind) {
        continue;
      }
      if (expr_idx == existing) {
        return true;
      }
      if (kind == semantic::ExprKind::Event &&
          program.expr_event_k[static_cast<std::size_t>(expr_idx)] ==
              program.expr_event_k[static_cast<std::size_t>(existing)] &&
          child_event_source_kind(program, expr_idx) ==
              child_event_source_kind(program, existing) &&
          child_event_source_index(program, expr_idx) ==
              child_event_source_index(program, existing)) {
        return true;
      }
    }
    return false;
  };
  for (const auto expr_idx : flattened) {
    if (!has_equivalent(expr_idx)) {
      normalized.push_back(expr_idx);
    }
  }
  if (normalized.size() == 1U) {
    return build_expr_transition_scenarios(plan, normalized.front());
  }
  return build_logical_transition_scenarios(plan, normalized, semantic::ExprKind::And);
}

inline std::vector<ExactTransitionScenario> build_expr_transition_scenarios(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx) {
  const auto &program = plan.lowered.program;
  const auto kind = static_cast<semantic::ExprKind>(
      program.expr_kind[static_cast<std::size_t>(expr_idx)]);

  if (kind == semantic::ExprKind::Impossible || kind == semantic::ExprKind::TrueExpr) {
    return {};
  }

  if (kind == semantic::ExprKind::Event) {
    return build_source_transition_scenarios(
        plan,
        ExactSourceKey{child_event_source_kind(program, expr_idx),
                       child_event_source_index(program, expr_idx)});
  }

  if (kind == semantic::ExprKind::Not) {
    return {};
  }

  if (kind == semantic::ExprKind::And || kind == semantic::ExprKind::Or) {
    std::vector<semantic::Index> children;
    const auto begin = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
    const auto end = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
    children.reserve(static_cast<std::size_t>(end - begin));
    for (semantic::Index i = begin; i < end; ++i) {
      children.push_back(program.expr_args[static_cast<std::size_t>(i)]);
    }
    return build_logical_transition_scenarios(plan, children, kind);
  }

  if (kind == semantic::ExprKind::Guard) {
    const auto ref = program.expr_ref_child[static_cast<std::size_t>(expr_idx)];
    const auto blocker = program.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
    const auto blocker_kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(blocker)]);
    const bool blocker_impossible = blocker_kind == semantic::ExprKind::Impossible;
    const bool blocker_true = blocker_kind == semantic::ExprKind::TrueExpr;
    bool any_unless_true = false;
    std::vector<semantic::Index> unless_events;
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      const auto child = program.expr_args[static_cast<std::size_t>(i)];
      const auto child_kind = static_cast<semantic::ExprKind>(
          program.expr_kind[static_cast<std::size_t>(child)]);
      if (child_kind == semantic::ExprKind::TrueExpr) {
        any_unless_true = true;
        continue;
      }
      if (child_kind == semantic::ExprKind::Impossible) {
        continue;
      }
      unless_events.push_back(child);
    }
    if (any_unless_true || blocker_impossible) {
      return build_expr_transition_scenarios(plan, ref);
    }
    if (blocker_true && unless_events.empty()) {
      return {};
    }

    const auto ref_scenarios = build_expr_transition_scenarios(plan, ref);
    std::vector<ExactTransitionScenario> out;
    const auto blocker_key =
        ExactSourceKey{child_event_source_kind(program, blocker),
                       child_event_source_index(program, blocker)};

    for (const auto &ref_scenario : ref_scenarios) {
      if (!blocker_true) {
        ExactTransitionScenario scenario = ref_scenario;
        std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
        for (const auto &constraint : ref_scenario.forced) {
          if (!append_constraint(&forced, constraint.key, constraint.relation)) {
            forced.clear();
            break;
          }
        }
        if (!(forced.empty() && !ref_scenario.forced.empty())) {
          bool ok = true;
          append_tail_requirement(&scenario, program, blocker);
          if (expr_is_simple_event(program, blocker)) {
            ok = append_factor(
                     &scenario,
                     ExactScenarioFactor{blocker_key, ExactFactorKind::AfterSurvival},
                     plan) &&
                 append_source_truth_constraints(
                     plan, blocker_key, ExactRelation::After, &forced);
          }
          for (const auto child : unless_events) {
            append_tail_requirement(&scenario, program, child);
            if (expr_is_simple_event(program, child)) {
              const auto unless_key =
                  ExactSourceKey{child_event_source_kind(program, child),
                                 child_event_source_index(program, child)};
              ok = ok &&
                   append_factor(
                       &scenario,
                       ExactScenarioFactor{unless_key, ExactFactorKind::AfterSurvival},
                       plan) &&
                   append_source_truth_constraints(
                       plan, unless_key, ExactRelation::After, &forced);
            }
          }
          if (ok) {
            scenario.forced = constraints_from_map(forced);
            out.push_back(std::move(scenario));
          }
        }
      }

      const int n_unless = static_cast<int>(unless_events.size());
      if (n_unless == 0) {
        continue;
      }
      const int max_mask = 1 << n_unless;
      for (int mask = 1; mask < max_mask; ++mask) {
        ExactTransitionScenario scenario = ref_scenario;
        std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
        for (const auto &constraint : ref_scenario.forced) {
          if (!append_constraint(&forced, constraint.key, constraint.relation)) {
            forced.clear();
            break;
          }
        }
        if (forced.empty() && !ref_scenario.forced.empty()) {
          continue;
        }
        bool ok = true;
        for (int bit = 0; bit < n_unless; ++bit) {
          const auto child = unless_events[static_cast<std::size_t>(bit)];
          const bool active = (mask & (1 << bit)) != 0;
          if (active) {
            append_ready_requirement(&scenario, program, child);
          } else {
            append_tail_requirement(&scenario, program, child);
          }
          if (expr_is_simple_event(program, child)) {
            const auto unless_key =
                ExactSourceKey{child_event_source_kind(program, child),
                               child_event_source_index(program, child)};
            ok = ok &&
                 append_factor(
                     &scenario,
                     ExactScenarioFactor{
                         unless_key,
                         active ? ExactFactorKind::BeforeCdf
                                : ExactFactorKind::AfterSurvival},
                     plan) &&
                 append_source_truth_constraints(
                     plan,
                     unless_key,
                     active ? ExactRelation::Before : ExactRelation::After,
                     &forced);
          }
        }
        if (ok) {
          scenario.forced = constraints_from_map(forced);
          out.push_back(std::move(scenario));
        }
      }
    }
    return out;
  }

  throw std::runtime_error("unsupported exact outcome expression");
}

inline void enumerate_competitor_subsets(
    const ExactVariantPlan &plan,
    const std::vector<semantic::Index> &block_outcomes,
    const std::size_t start,
    std::vector<semantic::Index> *current,
    std::vector<ExactCompetitorSubsetPlan> *out) {
  for (std::size_t i = start; i < block_outcomes.size(); ++i) {
    current->push_back(block_outcomes[i]);
    ExactCompetitorSubsetPlan subset;
    subset.outcome_indices = *current;
    subset.inclusion_sign = (current->size() % 2U == 1U) ? 1 : -1;
    subset.expr_roots.reserve(current->size());
    for (const auto outcome_idx : *current) {
      subset.expr_roots.push_back(
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].expr_root);
    }
    if (current->size() == 1U) {
      const auto outcome_idx = current->front();
      subset.singleton_expr_root =
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].expr_root;
      subset.scenarios =
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].scenarios;
    } else {
      subset.scenarios = build_expr_conjunction_scenarios(plan, subset.expr_roots);
    }
    if (!subset.scenarios.empty()) {
      out->push_back(std::move(subset));
    }
    enumerate_competitor_subsets(plan, block_outcomes, i + 1U, current, out);
    current->pop_back();
  }
}

inline std::vector<std::vector<semantic::Index>> build_competitor_overlap_blocks(
    const ExactVariantPlan &plan,
    const semantic::Index target_outcome_idx) {
  std::vector<semantic::Index> competitors;
  competitors.reserve(plan.outcomes.size());
  for (semantic::Index outcome_idx = 0;
       outcome_idx < static_cast<semantic::Index>(plan.outcomes.size());
       ++outcome_idx) {
    if (outcome_idx == target_outcome_idx) {
      continue;
    }
    if (plan.outcomes[static_cast<std::size_t>(outcome_idx)].scenarios.empty()) {
      continue;
    }
    competitors.push_back(outcome_idx);
  }
  std::vector<std::vector<semantic::Index>> blocks;
  std::vector<std::uint8_t> visited(competitors.size(), 0U);
  for (std::size_t start = 0; start < competitors.size(); ++start) {
    if (visited[start] != 0U) {
      continue;
    }
    std::vector<std::size_t> stack{start};
    visited[start] = 1U;
    std::vector<semantic::Index> block;
    while (!stack.empty()) {
      const auto idx = stack.back();
      stack.pop_back();
      const auto outcome_idx = competitors[idx];
      block.push_back(outcome_idx);
      for (std::size_t other = 0; other < competitors.size(); ++other) {
        if (visited[other] != 0U) {
          continue;
        }
        if (!supports_overlap(
                plan.outcomes[static_cast<std::size_t>(outcome_idx)].support,
                plan.outcomes[static_cast<std::size_t>(competitors[other])].support)) {
          continue;
        }
        visited[other] = 1U;
        stack.push_back(other);
      }
    }
    std::sort(block.begin(), block.end());
    blocks.push_back(std::move(block));
  }
  return blocks;
}

inline ExactTargetCompetitorPlan build_target_competitor_plan(
    const ExactVariantPlan &plan,
    const semantic::Index target_outcome_idx) {
  ExactTargetCompetitorPlan target_plan;
  const auto blocks = build_competitor_overlap_blocks(plan, target_outcome_idx);
  target_plan.blocks.reserve(blocks.size());
  for (const auto &block_outcomes : blocks) {
    ExactCompetitorBlockPlan block;
    block.outcome_indices = block_outcomes;
    std::vector<semantic::Index> current;
    enumerate_competitor_subsets(
        plan, block.outcome_indices, 0U, &current, &block.subsets);
    target_plan.blocks.push_back(std::move(block));
  }
  return target_plan;
}

inline bool scenario_supports_ranked_sequence(
    const ExactTransitionScenario &scenario) {
  return scenario.before_keys.empty() && scenario.ready_exprs.empty();
}

inline bool variant_supports_ranked_sequence(const ExactVariantPlan &plan) {
  for (const auto &outcome : plan.outcomes) {
    for (const auto &scenario : outcome.scenarios) {
      if (!scenario_supports_ranked_sequence(scenario)) {
        return false;
      }
    }
  }
  return true;
}

inline ExactVariantPlan make_exact_variant_plan(
    const runtime::LoweredExactVariant &lowered) {
  ExactVariantPlan plan;
  plan.lowered = lowered;
  ExactSupportBuilder builder(plan.lowered);
  plan.leaf_supports = builder.build_leaf_supports();
  plan.pool_supports = builder.build_pool_supports();
  plan.expr_supports = builder.build_expr_supports();

  const auto &program = plan.lowered.program;
  for (int i = 0; i < program.layout.n_triggers; ++i) {
    if (static_cast<semantic::TriggerKind>(
            program.trigger_kind[static_cast<std::size_t>(i)]) ==
            semantic::TriggerKind::Shared &&
        program.trigger_member_offsets[static_cast<std::size_t>(i + 1)] -
                program.trigger_member_offsets[static_cast<std::size_t>(i)] >
            1) {
      plan.shared_trigger_indices.push_back(i);
    }
  }

  plan.outcomes.reserve(plan.lowered.outcome_labels.size());
  for (std::size_t i = 0; i < plan.lowered.outcome_labels.size(); ++i) {
    const auto expr_root = program.outcome_expr_root[i];
    validate_exact_expr(plan.lowered, expr_root);
    ExactOutcomePlan outcome;
    outcome.expr_root = expr_root;
    outcome.support = plan.expr_supports[static_cast<std::size_t>(expr_root)];
    outcome.scenarios = build_expr_transition_scenarios(plan, expr_root);
    plan.outcome_index.emplace(plan.lowered.outcome_labels[i],
                               static_cast<semantic::Index>(i));
    plan.outcomes.push_back(std::move(outcome));
  }

  plan.competitor_plans.reserve(plan.outcomes.size());
  for (semantic::Index target_outcome_idx = 0;
       target_outcome_idx < static_cast<semantic::Index>(plan.outcomes.size());
       ++target_outcome_idx) {
    plan.competitor_plans.push_back(
        build_target_competitor_plan(plan, target_outcome_idx));
  }
  plan.ranked_supported = variant_supports_ranked_sequence(plan);

  return plan;
}

template <typename Fn>
double simpson_integrate(Fn &&fn, double lower, double upper, int n = 256) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return 0.0;
  }
  if (n < 2) {
    n = 2;
  }
  if (n % 2 != 0) {
    ++n;
  }
  const double h = (upper - lower) / static_cast<double>(n);
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double x = lower + h * static_cast<double>(i);
    const double fx = fn(x);
    const double w = (i == 0 || i == n) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
    sum += w * (std::isfinite(fx) ? fx : 0.0);
  }
  return sum * h / 3.0;
}

inline leaf::EventChannels forced_channels(const ExactRelation relation) {
  switch (relation) {
  case ExactRelation::Before:
  case ExactRelation::At:
    return leaf::EventChannels::certain();
  case ExactRelation::After:
    return leaf::EventChannels::impossible();
  case ExactRelation::Unknown:
    break;
  }
  return leaf::EventChannels::impossible();
}

class ExactSourceOracle {
public:
  ExactSourceOracle(const ExactVariantPlan &plan,
                    const ParamView &params,
                    const int first_param_row,
                    const ExactTriggerState &trigger_state,
                    const ExactSequenceState &sequence_state,
                    const double top_time)
      : plan_(plan),
        program_(plan.lowered.program),
        params_(params),
        first_param_row_(first_param_row),
        trigger_state_(trigger_state),
        sequence_state_(sequence_state),
        top_time_(top_time),
        leaf_cache_(static_cast<std::size_t>(program_.layout.n_leaves)),
        leaf_ready_(static_cast<std::size_t>(program_.layout.n_leaves), 0U),
        pool_cache_(static_cast<std::size_t>(program_.layout.n_pools)),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U) {}

  leaf::EventChannels source_channels(const semantic::SourceKind kind,
                                      const semantic::Index index,
                                      const double t) {
    if (std::fabs(t - top_time_) <= 1e-12) {
      if (kind == semantic::SourceKind::Leaf) {
        return cached_leaf(index, t);
      }
      if (kind == semantic::SourceKind::Pool) {
        return cached_pool(index, t);
      }
      throw std::runtime_error("exact kernel does not support special sources");
    }
    return compute_source_channels(kind, index, t);
  }

private:
  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  const ParamView &params_;
  int first_param_row_;
  const ExactTriggerState &trigger_state_;
  const ExactSequenceState &sequence_state_;
  double top_time_;
  std::vector<leaf::EventChannels> leaf_cache_;
  std::vector<std::uint8_t> leaf_ready_;
  std::vector<leaf::EventChannels> pool_cache_;
  std::vector<std::uint8_t> pool_ready_;

  const double *exact_time_for(const ExactSourceKey key) const {
    const auto it = sequence_state_.exact_times.find(key);
    if (it == sequence_state_.exact_times.end()) {
      return nullptr;
    }
    return &it->second;
  }

  leaf::EventChannels conditionalize(const leaf::EventChannels uncond,
                                     const leaf::EventChannels lower) const {
    if (!(sequence_state_.lower_bound > 0.0)) {
      return uncond;
    }
    const double surv_lb = lower.survival;
    if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
      return impossible_channels();
    }
    leaf::EventChannels out;
    out.pdf = safe_density(uncond.pdf / surv_lb);
    out.cdf = clamp_probability((uncond.cdf - lower.cdf) / surv_lb);
    out.survival = clamp_probability(uncond.survival / surv_lb);
    return out;
  }

  double leaf_q(const semantic::Index leaf_index, const int row) const {
    const auto trigger_index =
        program_.leaf_trigger_index[static_cast<std::size_t>(leaf_index)];
    if (trigger_index != semantic::kInvalidIndex &&
        static_cast<semantic::TriggerKind>(
            program_.trigger_kind[static_cast<std::size_t>(trigger_index)]) ==
            semantic::TriggerKind::Shared &&
        trigger_state_.shared_started[static_cast<std::size_t>(trigger_index)] <= 1U) {
      return trigger_state_.shared_started[static_cast<std::size_t>(trigger_index)] == 1U
                 ? 0.0
                 : 1.0;
    }
    return params_.q(row);
  }

  leaf::EventChannels cached_leaf(const semantic::Index index, const double t) {
    const auto pos = static_cast<std::size_t>(index);
    if (leaf_ready_[pos] != 0U) {
      return leaf_cache_[pos];
    }
    leaf_ready_[pos] = 1U;
    leaf_cache_[pos] = compute_source_channels(semantic::SourceKind::Leaf, index, t);
    return leaf_cache_[pos];
  }

  leaf::EventChannels cached_pool(const semantic::Index index, const double t) {
    const auto pos = static_cast<std::size_t>(index);
    if (pool_ready_[pos] != 0U) {
      return pool_cache_[pos];
    }
    pool_ready_[pos] = 1U;
    pool_cache_[pos] = compute_source_channels(semantic::SourceKind::Pool, index, t);
    return pool_cache_[pos];
  }

  leaf::EventChannels compute_source_channels(const semantic::SourceKind kind,
                                              const semantic::Index index,
                                              const double t) {
    const ExactSourceKey key{kind, index};
    if (const double *exact_time = exact_time_for(key)) {
      if (!(t >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    const auto uncond = base_source_channels(kind, index, t);
    if (!(sequence_state_.lower_bound > 0.0)) {
      return uncond;
    }
    const auto lower =
        base_source_channels(kind, index, sequence_state_.lower_bound);
    return conditionalize(uncond, lower);
  }

  leaf::EventChannels base_source_channels(const semantic::SourceKind kind,
                                           const semantic::Index index,
                                           const double t) {
    const ExactSourceKey key{kind, index};
    if (const double *exact_time = exact_time_for(key)) {
      if (!(t >= *exact_time)) {
        return impossible_channels();
      }
      return leaf::EventChannels::certain();
    }
    if (kind == semantic::SourceKind::Leaf) {
      return base_leaf_channels(index, t);
    }
    if (kind == semantic::SourceKind::Pool) {
      return base_pool_channels(index, t);
    }
    throw std::runtime_error("exact kernel does not support special sources");
  }

  leaf::EventChannels base_leaf_channels(const semantic::Index index,
                                         const double t) {
    const auto pos = static_cast<std::size_t>(index);
    const int row = first_param_row_ + index;
    const auto begin =
        program_.parameter_layout.leaf_param_offsets[pos];
    const auto end =
        program_.parameter_layout.leaf_param_offsets[pos + 1U];
    double local_params[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const int n_local = std::min<int>(end - begin, 8);
    for (int j = 0; j < n_local; ++j) {
      local_params[j] = params_.p(row, j);
    }

    const double q = leaf_q(index, row);
    const double t0 = params_.t0(row);
    const auto onset_kind =
        static_cast<semantic::OnsetKind>(program_.onset_kind[pos]);
    if (onset_kind == semantic::OnsetKind::Absolute) {
      return standard_leaf_channels(
          program_.leaf_dist_kind[pos],
          local_params,
          end - begin,
          q,
          t0,
          t - program_.onset_abs_value[pos]);
    }

    const double lag = program_.onset_lag[pos];
    const double upper = t - lag;
    if (!(upper > 0.0)) {
      return impossible_channels();
    }
    const auto source_kind =
        static_cast<semantic::SourceKind>(program_.onset_source_kind[pos]);
    const auto source_index = program_.onset_source_index[pos];
    const ExactSourceKey onset_key{source_kind, source_index};
    auto shifted_channels = [&](const double source_time) {
      return standard_leaf_channels(
          program_.leaf_dist_kind[pos],
          local_params,
          end - begin,
          q,
          t0,
          t - source_time - lag);
    };
    if (const double *exact_time = exact_time_for(onset_key)) {
      return shifted_channels(*exact_time);
    }
    const auto pdf = simpson_integrate(
        [&](const double u) {
          const auto source = base_source_channels(source_kind, source_index, u);
          return source.pdf * shifted_channels(u).pdf;
        },
        0.0,
        upper);
    const auto cdf = simpson_integrate(
        [&](const double u) {
          const auto source = base_source_channels(source_kind, source_index, u);
          return source.pdf * shifted_channels(u).cdf;
        },
        0.0,
        upper);

    leaf::EventChannels out;
    out.pdf = safe_density(pdf);
    out.cdf = clamp_probability(cdf);
    out.survival = clamp_probability(1.0 - out.cdf);
    return out;
  }

  leaf::EventChannels base_pool_channels(const semantic::Index index,
                                         const double t) {
    const auto pos = static_cast<std::size_t>(index);
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    const auto k = program_.pool_k[pos];
    const auto n_members = static_cast<int>(end - begin);
    if (n_members <= 0 || k < 1 || k > n_members) {
      return impossible_channels();
    }

    std::vector<leaf::EventChannels> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(base_source_channels(
          static_cast<semantic::SourceKind>(
              program_.pool_member_kind[static_cast<std::size_t>(i)]),
          program_.pool_member_indices[static_cast<std::size_t>(i)],
          t));
    }

    std::vector<std::vector<double>> prefix(
        static_cast<std::size_t>(n_members + 1),
        std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    std::vector<std::vector<double>> suffix(
        static_cast<std::size_t>(n_members + 1),
        std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    prefix[0][0] = 1.0;
    for (int i = 0; i < n_members; ++i) {
      for (int m = 0; m <= i; ++m) {
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].survival;
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m + 1)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }
    suffix[static_cast<std::size_t>(n_members)][0] = 1.0;
    for (int i = n_members - 1; i >= 0; --i) {
      const int count = n_members - i - 1;
      for (int m = 0; m <= count; ++m) {
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].survival;
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m + 1)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)].cdf;
      }
    }

    double surv = 0.0;
    for (int m = 0; m < k; ++m) {
      surv += prefix[static_cast<std::size_t>(n_members)][static_cast<std::size_t>(m)];
    }
    double pdf = 0.0;
    for (int i = 0; i < n_members; ++i) {
      double others_exact = 0.0;
      for (int left = 0; left < k; ++left) {
        const int right = (k - 1) - left;
        if (right < 0 || right > (n_members - i - 1)) {
          continue;
        }
        others_exact +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(left)] *
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(right)];
      }
      pdf += members[static_cast<std::size_t>(i)].pdf * others_exact;
    }

    leaf::EventChannels out;
    out.pdf = safe_density(pdf);
    out.survival = clamp_probability(surv);
    out.cdf = clamp_probability(1.0 - out.survival);
    return out;
  }
};

inline std::vector<ExactTriggerState> enumerate_trigger_states(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row) {
  std::vector<ExactTriggerState> states(1);
  states.front().weight = 1.0;
  states.front().shared_started.assign(
      static_cast<std::size_t>(plan.lowered.program.layout.n_triggers), 2U);

  for (const auto trigger_index : plan.shared_trigger_indices) {
    const auto pos = static_cast<std::size_t>(trigger_index);
    double q = 0.0;
    if (plan.lowered.program.trigger_has_fixed_q[pos] != 0U) {
      q = plan.lowered.program.trigger_fixed_q[pos];
    } else {
      const auto member_begin = plan.lowered.program.trigger_member_offsets[pos];
      if (member_begin != plan.lowered.program.trigger_member_offsets[pos + 1U]) {
        q = params.q(first_param_row +
                     plan.lowered.program.trigger_member_indices[static_cast<std::size_t>(
                         member_begin)]);
      }
    }
    q = clamp_probability(q);

    std::vector<ExactTriggerState> next;
    next.reserve(states.size() * 2U);
    for (const auto &state : states) {
      if (q > 0.0) {
        auto fail = state;
        fail.weight *= q;
        fail.shared_started[pos] = 0U;
        next.push_back(std::move(fail));
      }
      if (q < 1.0) {
        auto start = state;
        start.weight *= (1.0 - q);
        start.shared_started[pos] = 1U;
        next.push_back(std::move(start));
      }
    }
    states.swap(next);
  }

  return states;
}

class ForcedExprEvaluator {
public:
  ForcedExprEvaluator(const ExactVariantPlan &plan,
                      ExactSourceOracle *oracle,
                      const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash>
                          &forced)
      : plan_(plan),
        program_(plan.lowered.program),
        oracle_(oracle),
        forced_(forced),
        cdf_cache_(program_.expr_kind.size(), 0.0),
        cdf_time_(program_.expr_kind.size(),
                  std::numeric_limits<double>::quiet_NaN()),
        density_cache_(program_.expr_kind.size(), 0.0),
        density_time_(program_.expr_kind.size(),
                      std::numeric_limits<double>::quiet_NaN()) {}

  double expr_cdf(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    if (cdf_time_[pos] == oracle_time_) {
      return cdf_cache_[pos];
    }
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    double value = 0.0;
    switch (kind) {
    case semantic::ExprKind::Impossible:
      value = 0.0;
      break;
    case semantic::ExprKind::TrueExpr:
      value = 1.0;
      break;
    case semantic::ExprKind::Event:
      value = source_channels_forced(
                 ExactSourceKey{child_event_source_kind(program_, expr_idx),
                                child_event_source_index(program_, expr_idx)})
                  .cdf;
      break;
    case semantic::ExprKind::And: {
      double cdf = 1.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        cdf *= expr_cdf(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(cdf);
      break;
    }
    case semantic::ExprKind::Or: {
      double surv = 1.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        surv *= expr_survival(program_.expr_args[static_cast<std::size_t>(i)]);
      }
      value = clamp_probability(1.0 - surv);
      break;
    }
    case semantic::ExprKind::Guard: {
      const double current_time = oracle_time_;
      if (!(current_time > 0.0)) {
        value = 0.0;
        break;
      }
      const double cdf = simpson_integrate(
          [&](const double u) {
            oracle_time_ = u;
            return expr_density(expr_idx);
          },
          0.0,
          current_time);
      oracle_time_ = current_time;
      value = clamp_probability(cdf);
      break;
    }
    case semantic::ExprKind::Not:
      value = clamp_probability(
          1.0 - expr_cdf(program_.expr_args[static_cast<std::size_t>(
                    program_.expr_arg_offsets[pos])]));
      break;
    }
    cdf_time_[pos] = oracle_time_;
    cdf_cache_[pos] = value;
    return value;
  }

  double expr_survival(const semantic::Index expr_idx) {
    return clamp_probability(1.0 - expr_cdf(expr_idx));
  }

  double source_cdf(const ExactSourceKey key) {
    return clamp_probability(source_channels_forced(key).cdf);
  }

  double source_pdf(const ExactSourceKey key) {
    return safe_density(source_channels_forced(key).pdf);
  }

  double source_survival(const ExactSourceKey key) {
    return clamp_probability(source_channels_forced(key).survival);
  }

  ExactSourceOracle *oracle() const {
    return oracle_;
  }

  const ExactVariantPlan &plan() const {
    return plan_;
  }

  const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> &
  forced_map() const {
    return forced_;
  }

  double expr_density(const semantic::Index expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    if (density_time_[pos] == oracle_time_) {
      return density_cache_[pos];
    }
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    double value = 0.0;
    switch (kind) {
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      value = 0.0;
      break;
    case semantic::ExprKind::Event:
      value = source_channels_forced(
                 ExactSourceKey{child_event_source_kind(program_, expr_idx),
                                child_event_source_index(program_, expr_idx)})
                  .pdf;
      break;
    case semantic::ExprKind::And: {
      double total = 0.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        const auto child = program_.expr_args[static_cast<std::size_t>(i)];
        double term = expr_density(child);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
             j < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
             ++j) {
          if (i == j) {
            continue;
          }
          term *= expr_cdf(program_.expr_args[static_cast<std::size_t>(j)]);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        total += term;
      }
      value = clean_signed_value(total);
      break;
    }
    case semantic::ExprKind::Or: {
      double total = 0.0;
      for (semantic::Index i = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
           i < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
           ++i) {
        const auto child = program_.expr_args[static_cast<std::size_t>(i)];
        double term = expr_density(child);
        if (!std::isfinite(term) || term == 0.0) {
          continue;
        }
        for (semantic::Index j = program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
             j < program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
             ++j) {
          if (i == j) {
            continue;
          }
          term *= expr_survival(program_.expr_args[static_cast<std::size_t>(j)]);
          if (!std::isfinite(term) || term == 0.0) {
            break;
          }
        }
        total += term;
      }
      value = clean_signed_value(total);
      break;
    }
    case semantic::ExprKind::Guard: {
      const auto ref = program_.expr_ref_child[static_cast<std::size_t>(expr_idx)];
      const auto blocker = program_.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
      double blocker_survival = 0.0;
      if (program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)] ==
          program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)]) {
        blocker_survival = expr_survival(blocker);
      } else {
        const double current_time = oracle_time_;
        blocker_survival = clamp_probability(
            1.0 - simpson_integrate(
                      [&](const double u) {
                        oracle_time_ = u;
                        double term = expr_density(blocker);
                        if (!std::isfinite(term) || term == 0.0) {
                          return 0.0;
                        }
                        for (semantic::Index i =
                                 program_.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
                             i < program_.expr_arg_offsets[static_cast<std::size_t>(
                                      expr_idx + 1)];
                             ++i) {
                          term *= expr_survival(
                              program_.expr_args[static_cast<std::size_t>(i)]);
                          if (!std::isfinite(term) || term == 0.0) {
                            return 0.0;
                          }
                        }
                        return term;
                      },
                      0.0,
                      current_time));
        oracle_time_ = current_time;
      }
      value = clean_signed_value(expr_density(ref) * blocker_survival);
      break;
    }
    case semantic::ExprKind::Not:
      value = clean_signed_value(
          -expr_density(program_.expr_args[static_cast<std::size_t>(
               program_.expr_arg_offsets[pos])]));
      break;
    }
    density_time_[pos] = oracle_time_;
    density_cache_[pos] = value;
    return value;
  }

  double conjunction_density(const std::vector<semantic::Index> &exprs) {
    if (exprs.empty()) {
      return 0.0;
    }
    if (exprs.size() == 1U) {
      return expr_density(exprs.front());
    }
    double total = 0.0;
    for (std::size_t i = 0; i < exprs.size(); ++i) {
      double term = expr_density(exprs[i]);
      if (!std::isfinite(term) || term == 0.0) {
        continue;
      }
      for (std::size_t j = 0; j < exprs.size(); ++j) {
        if (i == j) {
          continue;
        }
        term *= expr_cdf(exprs[j]);
        if (!std::isfinite(term) || term == 0.0) {
          break;
        }
      }
      total += term;
    }
    return clean_signed_value(total);
  }

  double conjunction_cdf(const std::vector<semantic::Index> &exprs) {
    if (exprs.empty()) {
      return 0.0;
    }
    if (exprs.size() == 1U) {
      return expr_cdf(exprs.front());
    }
    const double current_time = oracle_time_;
    if (!(current_time > 0.0)) {
      return 0.0;
    }
    const double cdf = simpson_integrate(
        [&](const double u) {
          oracle_time_ = u;
          return conjunction_density(exprs);
        },
        0.0,
        current_time);
    oracle_time_ = current_time;
    return clamp_probability(cdf);
  }

private:
  const ExactVariantPlan &plan_;
  const runtime::ExactProgram &program_;
  ExactSourceOracle *oracle_;
  const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> &forced_;
  std::vector<double> cdf_cache_;
  std::vector<double> cdf_time_;
  std::vector<double> density_cache_;
  std::vector<double> density_time_;

  const std::vector<semantic::Index> &support_for(const ExactSourceKey key) const {
    if (key.kind == semantic::SourceKind::Leaf) {
      return plan_.leaf_supports[static_cast<std::size_t>(key.index)];
    }
    return plan_.pool_supports[static_cast<std::size_t>(key.index)];
  }

  leaf::EventChannels source_channels_forced(const ExactSourceKey key) {
    const auto it = forced_.find(key);
    if (it != forced_.end()) {
      return forced_channels(it->second);
    }
    const auto &support = support_for(key);
    for (const auto &[forced_key, relation] : forced_) {
      (void) relation;
      if (forced_key == key) {
        continue;
      }
      if (key.kind == semantic::SourceKind::Leaf &&
          forced_key.kind == semantic::SourceKind::Pool &&
          supports_overlap(support, support_for(forced_key))) {
        throw std::runtime_error(
            "exact kernel cannot condition a leaf on an overlapping forced pool state yet");
      }
      if (key.kind == semantic::SourceKind::Pool &&
          forced_key.kind == semantic::SourceKind::Pool &&
          supports_overlap(support, support_for(forced_key))) {
        throw std::runtime_error(
            "exact kernel cannot condition overlapping pool states yet");
      }
    }
    return oracle_->source_channels(key.kind, key.index, oracle_time_);
  }

public:
  double oracle_time_{0.0};
};

inline double source_cdf_at(ExactSourceOracle *oracle,
                            const ExactSourceKey key,
                            const double t) {
  return clamp_probability(oracle->source_channels(key.kind, key.index, t).cdf);
}

inline double source_pdf_at(ExactSourceOracle *oracle,
                            const ExactSourceKey key,
                            const double t) {
  return safe_density(oracle->source_channels(key.kind, key.index, t).pdf);
}

inline double source_survival_at(ExactSourceOracle *oracle,
                                 const ExactSourceKey key,
                                 const double t) {
  return clamp_probability(oracle->source_channels(key.kind, key.index, t).survival);
}

inline double readiness_cdf(const ExactTransitionScenario &scenario,
                            ForcedExprEvaluator *evaluator,
                            const double t) {
  if (scenario.before_keys.empty() && scenario.ready_exprs.empty()) {
    return 1.0;
  }
  if (t < 0.0) {
    return 0.0;
  }
  double value = 1.0;
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  for (const auto &key : scenario.before_keys) {
    value *= evaluator->source_cdf(key);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> expr_forced =
      evaluator->forced_map();
  for (const auto &constraint : scenario.forced) {
    if (!append_constraint(&expr_forced, constraint.key, constraint.relation)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  ForcedExprEvaluator expr_evaluator(evaluator->plan(), evaluator->oracle(), expr_forced);
  expr_evaluator.oracle_time_ = t;
  for (const auto expr_idx : scenario.ready_exprs) {
    value *= expr_evaluator.expr_cdf(expr_idx);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  evaluator->oracle_time_ = current_time;
  return clamp_probability(value);
}

inline double readiness_density(const ExactTransitionScenario &scenario,
                                ForcedExprEvaluator *evaluator,
                                const double t) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  if (scenario.before_keys.empty() && scenario.ready_exprs.empty()) {
    return 0.0;
  }
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> expr_forced =
      evaluator->forced_map();
  for (const auto &constraint : scenario.forced) {
    if (!append_constraint(&expr_forced, constraint.key, constraint.relation)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  ForcedExprEvaluator expr_evaluator(evaluator->plan(), evaluator->oracle(), expr_forced);
  expr_evaluator.oracle_time_ = t;
  double total = 0.0;
  for (std::size_t i = 0; i < scenario.before_keys.size(); ++i) {
    double term = evaluator->source_pdf(scenario.before_keys[i]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < scenario.before_keys.size(); ++j) {
      if (j == i) {
        continue;
      }
      term *= evaluator->source_cdf(scenario.before_keys[j]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (const auto expr_idx : scenario.ready_exprs) {
      term *= expr_evaluator.expr_cdf(expr_idx);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  for (std::size_t i = 0; i < scenario.ready_exprs.size(); ++i) {
    double term = expr_evaluator.expr_density(scenario.ready_exprs[i]);
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (const auto &key : scenario.before_keys) {
      term *= evaluator->source_cdf(key);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    if (!std::isfinite(term) || term == 0.0) {
      continue;
    }
    for (std::size_t j = 0; j < scenario.ready_exprs.size(); ++j) {
      if (i == j) {
        continue;
      }
      term *= expr_evaluator.expr_cdf(scenario.ready_exprs[j]);
      if (!std::isfinite(term) || term == 0.0) {
        break;
      }
    }
    total += term;
  }
  evaluator->oracle_time_ = current_time;
  return clean_signed_value(total);
}

inline double after_survival(const ExactTransitionScenario &scenario,
                             ForcedExprEvaluator *evaluator,
                             const double t) {
  double value = 1.0;
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  for (const auto &key : scenario.after_keys) {
    value *= evaluator->source_survival(key);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> expr_forced =
      evaluator->forced_map();
  for (const auto &constraint : scenario.forced) {
    if (!append_constraint(&expr_forced, constraint.key, constraint.relation)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  ForcedExprEvaluator expr_evaluator(evaluator->plan(), evaluator->oracle(), expr_forced);
  expr_evaluator.oracle_time_ = t;
  for (const auto expr_idx : scenario.tail_exprs) {
    value *= expr_evaluator.expr_survival(expr_idx);
    if (!(value > 0.0)) {
      evaluator->oracle_time_ = current_time;
      return 0.0;
    }
  }
  evaluator->oracle_time_ = current_time;
  return clamp_probability(value);
}

inline double same_active_win_mass(const ExactTransitionScenario &scenario,
                                   ForcedExprEvaluator *evaluator,
                                   const double t,
                                   const double ready_upper) {
  if (!(ready_upper > 0.0)) {
    return 0.0;
  }
  const double tail = after_survival(scenario, evaluator, t);
  if (!(tail > 0.0)) {
    return 0.0;
  }
  return tail * readiness_cdf(scenario, evaluator, ready_upper);
}

inline double scenario_truth_cdf(const ExactTransitionScenario &scenario,
                                 ForcedExprEvaluator *evaluator,
                                 const double t) {
  if (!(t > 0.0)) {
    return 0.0;
  }
  const double current_time = evaluator->oracle_time_;
  const double integrated = simpson_integrate(
      [&](const double u) {
        evaluator->oracle_time_ = u;
        double value = evaluator->source_pdf(scenario.active_key);
        if (!(value > 0.0)) {
          return 0.0;
        }
        value *= after_survival(scenario, evaluator, u);
        if (!(value > 0.0)) {
          return 0.0;
        }
        if (!scenario.before_keys.empty() || !scenario.ready_exprs.empty()) {
          value *= readiness_cdf(scenario, evaluator, u);
          if (!(value > 0.0)) {
            return 0.0;
          }
        }
        return value;
      },
      0.0,
      t);
  evaluator->oracle_time_ = current_time;
  return clamp_probability(integrated);
}

inline double competitor_subset_win_mass(
    const ExactCompetitorSubsetPlan &subset_plan,
    const ExactSourceKey target_active_key,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
    const double t) {
  const double current_time = evaluator->oracle_time_;
  evaluator->oracle_time_ = t;
  double win_prob = 0.0;
  if (subset_plan.outcome_indices.size() == 1U &&
      subset_plan.singleton_expr_root != semantic::kInvalidIndex) {
    win_prob = evaluator->expr_cdf(subset_plan.singleton_expr_root);
  } else {
    for (const auto &scenario : subset_plan.scenarios) {
      win_prob += scenario_truth_cdf(scenario, evaluator, t);
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
    same_active_t += same_active_win_mass(scenario, evaluator, t, t);
    same_active_ready +=
        same_active_win_mass(scenario, evaluator, t, readiness_upper);
  }
  evaluator->oracle_time_ = current_time;
  return clamp_probability(win_prob - same_active_t + same_active_ready);
}

inline double competitor_non_win_probability(
    const ExactTargetCompetitorPlan &competitor_plan,
    const ExactSourceKey target_active_key,
    const double readiness_upper,
    ForcedExprEvaluator *evaluator,
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
                        subset, target_active_key, readiness_upper, evaluator, t);
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
    const ExactSequenceState &sequence_state,
    const ExactTransitionScenario &scenario,
    const double observed_time) {
  ExactSequenceState next_state = sequence_state;
  next_state.lower_bound = observed_time;
  next_state.exact_times[scenario.active_key] = observed_time;
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
    const bool collect_successors = false) {
  ExactStepResult result;
  const auto target_pos = static_cast<std::size_t>(target_idx);
  const auto &competitor_plan = plan.competitor_plans[target_pos];
  ExactSourceOracle oracle(
      plan, params, first_param_row, trigger_state, sequence_state, observed_time);

  for (const auto &scenario : plan.outcomes[target_pos].scenarios) {
    if (collect_successors && !scenario_supports_ranked_sequence(scenario)) {
      throw std::logic_error(
          "ranked exact step requested for unsupported latent-readiness scenario");
    }

    std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash> forced;
    bool ok = true;
    for (const auto &constraint : scenario.forced) {
      if (!append_constraint(&forced, constraint.key, constraint.relation)) {
        ok = false;
        break;
      }
    }
    if (!ok) {
      continue;
    }

    const std::unordered_map<ExactSourceKey, ExactRelation, ExactSourceKeyHash>
        no_forced;
    ForcedExprEvaluator target_evaluator(plan, &oracle, no_forced);
    target_evaluator.oracle_time_ = observed_time;
    ForcedExprEvaluator competitor_evaluator(plan, &oracle, forced);
    competitor_evaluator.oracle_time_ = observed_time;

    const double active_pdf =
        source_pdf_at(&oracle, scenario.active_key, observed_time);
    const double tail = after_survival(scenario, &target_evaluator, observed_time);
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
              observed_time,
              used_outcomes);
      if (!(scenario_prob > 0.0)) {
        continue;
      }
      result.total_probability += scenario_prob;
      if (collect_successors) {
        result.branches.push_back(ExactStepBranch{
            scenario_prob,
            advance_exact_sequence_state(sequence_state, scenario, observed_time)});
      }
      continue;
    }

    const double initial_ready = readiness_cdf(scenario, &target_evaluator, 0.0);
    if (initial_ready > 0.0) {
      result.total_probability += initial_ready * active_pdf * tail *
                                  competitor_non_win_probability(
                                      competitor_plan,
                                      scenario.active_key,
                                      0.0,
                                      &competitor_evaluator,
                                      observed_time,
                                      used_outcomes);
    }

    result.total_probability += simpson_integrate(
        [&](const double readiness_time) {
          double value =
              readiness_density(scenario, &target_evaluator, readiness_time);
          if (!std::isfinite(value) || value == 0.0) {
            return 0.0;
          }
          value *= active_pdf * tail;
          if (!std::isfinite(value) || value == 0.0) {
            return 0.0;
          }
          value *= competitor_non_win_probability(
              competitor_plan,
              scenario.active_key,
              readiness_time,
              &competitor_evaluator,
              observed_time,
              used_outcomes);
          if (!std::isfinite(value) || value == 0.0) {
            return 0.0;
          }
          return value;
        },
        0.0,
        observed_time);
  }

  return result;
}

inline std::vector<ExactTrialObservation> collapse_exact_observations(
    const Rcpp::DataFrame &data,
    const std::unordered_map<std::string, semantic::Index> &variant_index_by_component) {
  const auto n_rows = data.nrows();
  if (n_rows == 0) {
    return {};
  }
  if (!data.containsElementNamed("R") || !data.containsElementNamed("rt")) {
    throw std::runtime_error("exact evaluator data must include columns R and rt");
  }

  int max_rank = 1;
  for (int rank = 2;; ++rank) {
    const auto r_col = "R" + std::to_string(rank);
    const auto rt_col = "rt" + std::to_string(rank);
    const bool has_r = data.containsElementNamed(r_col.c_str());
    const bool has_rt = data.containsElementNamed(rt_col.c_str());
    if (has_r != has_rt) {
      throw std::runtime_error(
          "exact evaluator ranked columns must appear as matched Rk/rtk pairs");
    }
    if (!has_r) {
      break;
    }
    max_rank = rank;
  }

  std::vector<Rcpp::CharacterVector> labels(static_cast<std::size_t>(max_rank + 1));
  std::vector<Rcpp::NumericVector> times(static_cast<std::size_t>(max_rank + 1));
  labels[1] = Rcpp::as<Rcpp::CharacterVector>(data["R"]);
  times[1] = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  for (int rank = 2; rank <= max_rank; ++rank) {
    labels[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::CharacterVector>(data[("R" + std::to_string(rank)).c_str()]);
    times[static_cast<std::size_t>(rank)] =
        Rcpp::as<Rcpp::NumericVector>(data[("rt" + std::to_string(rank)).c_str()]);
  }

  Rcpp::IntegerVector trial =
      data.containsElementNamed("trial")
          ? Rcpp::as<Rcpp::IntegerVector>(data["trial"])
          : Rcpp::seq_len(n_rows);
  const bool has_component = data.containsElementNamed("component");
  Rcpp::CharacterVector component =
      has_component ? Rcpp::as<Rcpp::CharacterVector>(data["component"])
                    : Rcpp::CharacterVector(n_rows, NA_STRING);

  std::vector<ExactTrialObservation> out;
  out.reserve(static_cast<std::size_t>(n_rows));
  int last_trial = NA_INTEGER;
  for (R_xlen_t i = 0; i < n_rows; ++i) {
    const int trial_id = trial[i];
    if (i > 0 && trial_id == last_trial) {
      continue;
    }
    last_trial = trial_id;

    std::string component_id = "__default__";
    if (has_component && component[i] != NA_STRING) {
      component_id = Rcpp::as<std::string>(component[i]);
    }
    const auto it = variant_index_by_component.find(component_id);
    if (it == variant_index_by_component.end()) {
      throw std::runtime_error(
          "exact evaluator found no exact variant for component '" + component_id + "'");
    }

    ExactTrialObservation obs;
    obs.variant_index = it->second;
    obs.component_id = component_id;
    double prev_rt = -std::numeric_limits<double>::infinity();
    std::unordered_map<std::string, std::uint8_t> seen_labels;
    for (int rank = 1; rank <= max_rank; ++rank) {
      const auto &rank_labels = labels[static_cast<std::size_t>(rank)];
      const auto &rank_times = times[static_cast<std::size_t>(rank)];
      const SEXP label_sexp = rank_labels[i];
      const bool label_missing = label_sexp == NA_STRING;
      const double rt = rank_times[i];
      const bool time_missing = Rcpp::NumericVector::is_na(rt);
      if (label_missing && time_missing) {
        break;
      }
      if (label_missing || time_missing || !std::isfinite(rt) || !(rt > prev_rt)) {
        obs.valid = false;
        break;
      }
      const auto label = Rcpp::as<std::string>(rank_labels[i]);
      if (seen_labels.find(label) != seen_labels.end()) {
        obs.valid = false;
        break;
      }
      seen_labels.emplace(label, 1U);
      obs.ranks.push_back(ExactObservedRank{label, rt});
      prev_rt = rt;
    }
    if (obs.ranks.empty()) {
      obs.valid = false;
    }
    out.push_back(std::move(obs));
  }
  return out;
}

inline std::vector<runtime::TrialBlock> build_trial_blocks(
    const std::vector<ExactTrialObservation> &observations) {
  std::vector<runtime::TrialBlock> blocks;
  if (observations.empty()) {
    return blocks;
  }
  runtime::TrialBlock current;
  current.variant_index = observations.front().variant_index;
  current.start_row = 0;
  current.row_count = 1;
  for (std::size_t i = 1; i < observations.size(); ++i) {
    if (observations[i].variant_index == current.variant_index) {
      ++current.row_count;
      continue;
    }
    blocks.push_back(current);
    current.variant_index = observations[i].variant_index;
    current.start_row = static_cast<int>(i);
    current.row_count = 1;
  }
  blocks.push_back(current);
  return blocks;
}

inline double exact_loglik_for_trial(const ExactVariantPlan &plan,
                                     const ParamView &params,
                                     const int first_param_row,
                                     const std::string &label,
                                     const double rt,
                                     const double min_ll) {
  if (!std::isfinite(rt) || !(rt > 0.0)) {
    return min_ll;
  }
  const auto target_it = plan.outcome_index.find(label);
  if (target_it == plan.outcome_index.end()) {
    return min_ll;
  }
  const auto target_idx = target_it->second;
  const auto trigger_states = enumerate_trigger_states(plan, params, first_param_row);
  double total = 0.0;

  for (const auto &trigger_state : trigger_states) {
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    const ExactStepResult step = evaluate_exact_step(
        plan,
        params,
        first_param_row,
        trigger_state,
        ExactSequenceState{},
        target_idx,
        rt);
    total += trigger_state.weight * step.total_probability;
  }

  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

inline double exact_ranked_sequence_probability(
    const ExactVariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const ExactTriggerState &trigger_state,
    const std::vector<ExactObservedRank> &ranks,
    const std::size_t rank_idx,
    const ExactSequenceState &sequence_state,
    std::vector<std::uint8_t> *used_outcomes) {
  if (rank_idx >= ranks.size()) {
    return 1.0;
  }

  const auto &rank = ranks[rank_idx];
  const auto target_it = plan.outcome_index.find(rank.outcome_label);
  if (target_it == plan.outcome_index.end()) {
    return 0.0;
  }
  const auto target_idx = static_cast<std::size_t>(target_it->second);
  if ((*used_outcomes)[target_idx] != 0U) {
    return 0.0;
  }
  const ExactStepResult step = evaluate_exact_step(
      plan,
      params,
      first_param_row,
      trigger_state,
      sequence_state,
      target_it->second,
      rank.rt,
      used_outcomes,
      true);
  double total = 0.0;
  for (const auto &branch : step.branches) {
    (*used_outcomes)[target_idx] = 1U;
    const double tail_prob = exact_ranked_sequence_probability(
        plan,
        params,
        first_param_row,
        trigger_state,
        ranks,
        rank_idx + 1U,
        branch.next_state,
        used_outcomes);
    (*used_outcomes)[target_idx] = 0U;
    total += branch.probability * tail_prob;
  }
  return total;
}

inline double exact_ranked_loglik_for_trial(const ExactVariantPlan &plan,
                                            const ParamView &params,
                                            const int first_param_row,
                                            const ExactTrialObservation &obs,
                                            const double min_ll) {
  if (!obs.valid || obs.ranks.empty()) {
    return min_ll;
  }
  if (obs.ranks.size() == 1U) {
    return exact_loglik_for_trial(
        plan, params, first_param_row, obs.ranks.front().outcome_label, obs.ranks.front().rt, min_ll);
  }

  const auto trigger_states = enumerate_trigger_states(plan, params, first_param_row);
  double total = 0.0;
  for (const auto &trigger_state : trigger_states) {
    if (!(trigger_state.weight > 0.0)) {
      continue;
    }
    std::vector<std::uint8_t> used_outcomes(plan.outcomes.size(), 0U);
    ExactSequenceState state;
    total += trigger_state.weight *
             exact_ranked_sequence_probability(
                 plan, params, first_param_row, trigger_state, obs.ranks, 0U, state, &used_outcomes);
  }
  if (!std::isfinite(total) || !(total > 0.0)) {
    return min_ll;
  }
  return std::log(total);
}

inline Rcpp::List evaluate_exact_outcome_probabilities(
    const compile::CompiledModel &compiled,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  std::unordered_map<std::string, semantic::Index> variant_index_by_component;
  std::vector<runtime::LoweredExactVariant> lowered_variants;
  lowered_variants.reserve(compiled.variants.size());

  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    variant_index_by_component.emplace(
        variant.component_id,
        static_cast<semantic::Index>(lowered_variants.size()));
    lowered_variants.push_back(runtime::lower_exact_variant(variant));
  }
  if (lowered_variants.empty()) {
    throw std::runtime_error("exact probability evaluator found no exact variants");
  }

  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const auto queries = collapse_probability_queries(data, variant_index_by_component);
  const auto blocks = build_trial_blocks(queries);

  std::vector<ExactVariantPlan> plans;
  plans.reserve(lowered_variants.size());
  for (const auto &variant : lowered_variants) {
    plans.push_back(make_exact_variant_plan(variant));
  }

  Rcpp::NumericVector prob(queries.size());
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    for (int local = 0; local < block.row_count; ++local) {
      const auto query_idx = static_cast<std::size_t>(block.start_row + local);
      if (param_row + leaf_count > static_cast<std::size_t>(params.matrix.nrow())) {
        throw std::runtime_error(
            "parameter rows do not match exact probability query layout");
      }
      prob[static_cast<R_xlen_t>(query_idx)] = integrate_to_infinity(
          [&](const double rt) {
            const double log_density = exact_loglik_for_trial(
                plan,
                params,
                static_cast<int>(param_row),
                queries[query_idx].outcome_label,
                rt,
                -std::numeric_limits<double>::infinity());
            if (!std::isfinite(log_density)) {
              return 0.0;
            }
            const double density = std::exp(log_density);
            return std::isfinite(density) && density > 0.0 ? density : 0.0;
          });
      param_row += leaf_count;
    }
  }

  if (param_row != static_cast<std::size_t>(params.matrix.nrow())) {
    throw std::runtime_error(
        "parameter rows contain unused entries for the exact probability layout");
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = prob,
      Rcpp::Named("n_trials") = static_cast<int>(queries.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()));
}

inline SEXP evaluate_exact_trials(const compile::CompiledModel &compiled,
                                  SEXP paramsSEXP,
                                  SEXP dataSEXP,
                                  const double min_ll) {
  std::unordered_map<std::string, semantic::Index> variant_index_by_component;
  std::vector<runtime::LoweredExactVariant> lowered_variants;
  lowered_variants.reserve(compiled.variants.size());

  for (const auto &variant : compiled.variants) {
    if (variant.backend != compile::BackendKind::Exact) {
      continue;
    }
    variant_index_by_component.emplace(
        variant.component_id,
        static_cast<semantic::Index>(lowered_variants.size()));
    lowered_variants.push_back(runtime::lower_exact_variant(variant));
  }
  if (lowered_variants.empty()) {
    throw std::runtime_error("exact evaluator found no exact variants");
  }

  ParamView params(paramsSEXP);
  Rcpp::DataFrame data(dataSEXP);
  const auto observations = collapse_exact_observations(data, variant_index_by_component);
  const auto blocks = build_trial_blocks(observations);

  std::vector<ExactVariantPlan> plans;
  plans.reserve(lowered_variants.size());
  for (const auto &variant : lowered_variants) {
    plans.push_back(make_exact_variant_plan(variant));
  }

  Rcpp::NumericVector loglik(observations.size(), min_ll);
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans[static_cast<std::size_t>(block.variant_index)];
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    if (!plan.ranked_supported) {
      for (int i = 0; i < block.row_count; ++i) {
        const auto row = static_cast<std::size_t>(block.start_row + i);
        if (observations[row].ranks.size() > 1U) {
          throw std::runtime_error(
              "exact ranked kernel does not support latent prerequisite timings yet");
        }
      }
    }
    for (int i = 0; i < block.row_count; ++i) {
      const auto row = block.start_row + i;
      if (param_row + leaf_count > static_cast<std::size_t>(params.matrix.nrow())) {
        throw std::runtime_error(
            "parameter rows do not match exact trial layout");
      }
      loglik[row] = exact_ranked_loglik_for_trial(
          plan,
          params,
          static_cast<int>(param_row),
          observations[static_cast<std::size_t>(row)],
          min_ll);
      param_row += leaf_count;
    }
  }

  if (param_row != static_cast<std::size_t>(params.matrix.nrow())) {
    throw std::runtime_error(
        "parameter rows contain unused entries for the exact trial layout");
  }

  Rcpp::List block_list(blocks.size());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    const auto &block = blocks[i];
    block_list[i] = Rcpp::List::create(
        Rcpp::Named("component_id") =
            lowered_variants[static_cast<std::size_t>(block.variant_index)].component_id,
        Rcpp::Named("start_row") = block.start_row + 1,
        Rcpp::Named("row_count") = block.row_count);
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    total_loglik += static_cast<double>(loglik[i]);
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("blocks") = block_list);
}

} // namespace detail
} // namespace accumulatr::eval
