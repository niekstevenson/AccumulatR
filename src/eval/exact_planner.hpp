#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

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

inline bool scenario_key_happens_by_transition(
    const ExactTransitionScenario &scenario,
    const ExactSourceKey key) {
  if (scenario.active_key == key) {
    return true;
  }
  return std::find_if(
             scenario.before_keys.begin(),
             scenario.before_keys.end(),
             [&](const ExactSourceKey &existing) { return existing == key; }) !=
         scenario.before_keys.end();
}

inline bool scenario_key_known_before_active(
    const ExactTransitionScenario &scenario,
    const ExactSourceKey key) {
  return std::find_if(
             scenario.before_keys.begin(),
             scenario.before_keys.end(),
             [&](const ExactSourceKey &existing) { return existing == key; }) !=
         scenario.before_keys.end();
}

inline semantic::Index source_ordinal(const ExactVariantPlan &plan,
                                      const ExactSourceKey key);

inline void append_source_order_fact(
    std::vector<ExactSourceOrderFact> *facts,
    const ExactVariantPlan &plan,
    const ExactSourceKey before_key,
    const ExactSourceKey after_key) {
  const auto before_source_id = source_ordinal(plan, before_key);
  const auto after_source_id = source_ordinal(plan, after_key);
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

inline void append_scenario_source_order_fact(
    ExactTransitionScenario *scenario,
    const ExactVariantPlan &plan,
    const ExactSourceKey before_key,
    const ExactSourceKey after_key) {
  append_source_order_fact(
      &scenario->source_order_facts, plan, before_key, after_key);
}

inline bool append_tail_order_requirement(ExactTransitionScenario *scenario,
                                          const ExactVariantPlan &plan,
                                          const semantic::Index expr_idx) {
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind != semantic::ExprKind::Guard ||
      program.expr_arg_offsets[pos] != program.expr_arg_offsets[pos + 1U]) {
    return true;
  }
  const auto ref = program.expr_ref_child[pos];
  const auto blocker = program.expr_blocker_child[pos];
  if (!expr_is_simple_event(program, ref) ||
      !expr_is_simple_event(program, blocker)) {
    return true;
  }
  const auto ref_key =
      ExactSourceKey{child_event_source_kind(program, ref),
                     child_event_source_index(program, ref)};
  const auto blocker_key =
      ExactSourceKey{child_event_source_kind(program, blocker),
                     child_event_source_index(program, blocker)};
  if (!scenario_key_happens_by_transition(*scenario, ref_key) ||
      !scenario_key_happens_by_transition(*scenario, blocker_key)) {
    return true;
  }
  if (scenario->active_key == blocker_key &&
      scenario_key_known_before_active(*scenario, ref_key)) {
    return false;
  }
  const auto before_source_id = source_ordinal(plan, blocker_key);
  const auto after_source_id = source_ordinal(plan, ref_key);
  for (const auto &fact : scenario->source_order_facts) {
    if (fact.before_source_id == before_source_id &&
        fact.after_source_id == after_source_id) {
      return true;
    }
  }
  append_scenario_source_order_fact(scenario, plan, blocker_key, ref_key);
  return true;
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

inline const std::vector<semantic::Index> &planned_source_support_for(
    const ExactVariantPlan &plan,
    const ExactSourceKey key) {
  if (key.kind == semantic::SourceKind::Leaf) {
    return plan.leaf_supports[static_cast<std::size_t>(key.index)];
  }
  return plan.pool_supports[static_cast<std::size_t>(key.index)];
}

inline bool scenario_forced_sources_supported(const ExactVariantPlan &plan,
                                              const ExactTransitionScenario &scenario) {
  std::vector<ExactSourceKey> referenced_keys;
  referenced_keys.reserve(
      1U + scenario.before_keys.size() + scenario.after_keys.size());
  referenced_keys.push_back(scenario.active_key);
  referenced_keys.insert(
      referenced_keys.end(), scenario.before_keys.begin(), scenario.before_keys.end());
  referenced_keys.insert(
      referenced_keys.end(), scenario.after_keys.begin(), scenario.after_keys.end());

  for (const auto &constraint : scenario.forced) {
    if (constraint.key.kind != semantic::SourceKind::Pool) {
      continue;
    }
    const auto &forced_support = planned_source_support_for(plan, constraint.key);
    for (const auto &other : scenario.forced) {
      if (other.key == constraint.key || other.key.kind != semantic::SourceKind::Pool) {
        continue;
      }
      if (supports_overlap(
              forced_support,
              planned_source_support_for(plan, other.key))) {
        return false;
      }
    }
    for (const auto &key : referenced_keys) {
      if (key == constraint.key) {
        continue;
      }
      if (supports_overlap(forced_support, planned_source_support_for(plan, key))) {
        return false;
      }
    }
    for (const auto expr_idx : scenario.ready_exprs) {
      if (supports_overlap(
              forced_support,
              plan.expr_supports[static_cast<std::size_t>(expr_idx)])) {
        return false;
      }
    }
    for (const auto expr_idx : scenario.tail_exprs) {
      if (supports_overlap(
              forced_support,
              plan.expr_supports[static_cast<std::size_t>(expr_idx)])) {
        return false;
      }
    }
  }
  return true;
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
    if (scenario_forced_sources_supported(plan, scenario)) {
      out.push_back(std::move(scenario));
    }
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
      if (scenario_forced_sources_supported(plan, scenario)) {
        out.push_back(std::move(scenario));
      }
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
          append_ready_requirement(&scenario, plan.lowered.program, child);
          if (expr_is_simple_event(plan.lowered.program, child)) {
            const auto child_key =
                ExactSourceKey{child_event_source_kind(plan.lowered.program, child),
                               child_event_source_index(plan.lowered.program, child)};
            ok = append_factor(&scenario,
                               ExactScenarioFactor{child_key, ExactFactorKind::BeforeCdf},
                               plan) &&
                 append_source_truth_constraints(
                     plan, child_key, ExactRelation::Before, &forced);
          }
        } else {
          append_tail_requirement(&scenario, plan.lowered.program, child);
          ok = append_tail_order_requirement(&scenario, plan, child);
          if (expr_is_simple_event(plan.lowered.program, child)) {
            const auto child_key =
                ExactSourceKey{child_event_source_kind(plan.lowered.program, child),
                               child_event_source_index(plan.lowered.program, child)};
            ok = ok &&
                 append_factor(
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
      if (scenario_forced_sources_supported(plan, scenario)) {
        out.push_back(std::move(scenario));
      }
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
    const auto blocker =
        program.expr_blocker_child[static_cast<std::size_t>(expr_idx)];
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
            if (scenario_forced_sources_supported(plan, scenario)) {
              out.push_back(std::move(scenario));
            }
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
          if (scenario_forced_sources_supported(plan, scenario)) {
            out.push_back(std::move(scenario));
          }
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
    std::vector<semantic::Index> expr_roots;
    expr_roots.reserve(current->size());
    for (const auto outcome_idx : *current) {
      expr_roots.push_back(
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].expr_root);
    }
    if (current->size() == 1U) {
      const auto outcome_idx = current->front();
      subset.singleton_expr_root =
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].expr_root;
      subset.scenarios =
          plan.outcomes[static_cast<std::size_t>(outcome_idx)].scenarios;
    } else {
      subset.scenarios = build_expr_conjunction_scenarios(plan, expr_roots);
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
                plan.expr_supports[static_cast<std::size_t>(
                    plan.outcomes[static_cast<std::size_t>(outcome_idx)]
                        .expr_root)],
                plan.expr_supports[static_cast<std::size_t>(
                    plan.outcomes[static_cast<std::size_t>(competitors[other])]
                        .expr_root)])) {
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
    std::vector<semantic::Index> current;
    enumerate_competitor_subsets(plan, block_outcomes, 0U, &current, &block.subsets);
    target_plan.blocks.push_back(std::move(block));
  }
  return target_plan;
}

inline bool scenario_supports_ranked_sequence(
    const ExactTransitionScenario &scenario) {
  return scenario.before_source_span.empty() && scenario.ready_expr_span.empty();
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

inline void compile_index_span(const std::vector<semantic::Index> &values,
                               std::vector<semantic::Index> *arena,
                               ExactIndexSpan *span) {
  span->offset = static_cast<semantic::Index>(arena->size());
  span->size = static_cast<semantic::Index>(values.size());
  arena->insert(arena->end(), values.begin(), values.end());
}

template <typename T>
inline void release_runtime_vector(std::vector<T> *values) {
  std::vector<T>().swap(*values);
}

inline void release_scenario_planning_fields(ExactTransitionScenario *scenario) {
  release_runtime_vector(&scenario->before_keys);
  release_runtime_vector(&scenario->before_source_ids);
  release_runtime_vector(&scenario->after_keys);
  release_runtime_vector(&scenario->after_source_ids);
  release_runtime_vector(&scenario->ready_exprs);
  release_runtime_vector(&scenario->tail_exprs);
  release_runtime_vector(&scenario->factors);
  release_runtime_vector(&scenario->forced);
  release_runtime_vector(&scenario->source_order_facts);
}

inline void compile_scenario_runtime_fields(ExactVariantPlan *plan,
                                            ExactTransitionScenario *scenario) {
  scenario->active_source_id = source_ordinal(*plan, scenario->active_key);
  scenario->before_source_ids.clear();
  scenario->before_source_ids.reserve(scenario->before_keys.size());
  for (const auto &key : scenario->before_keys) {
    scenario->before_source_ids.push_back(source_ordinal(*plan, key));
  }
  scenario->after_source_ids.clear();
  scenario->after_source_ids.reserve(scenario->after_keys.size());
  for (const auto &key : scenario->after_keys) {
    scenario->after_source_ids.push_back(source_ordinal(*plan, key));
  }
  compile_index_span(
      scenario->before_source_ids,
      &plan->scenario_source_ids,
      &scenario->before_source_span);
  compile_index_span(
      scenario->after_source_ids,
      &plan->scenario_source_ids,
      &scenario->after_source_span);
  compile_index_span(
      scenario->ready_exprs,
      &plan->scenario_expr_ids,
      &scenario->ready_expr_span);
  compile_index_span(
      scenario->tail_exprs,
      &plan->scenario_expr_ids,
      &scenario->tail_expr_span);

  std::vector<std::pair<semantic::Index, ExactRelation>> relation_pairs;
  relation_pairs.reserve(scenario->forced.size());
  for (const auto &constraint : scenario->forced) {
    relation_pairs.emplace_back(
        source_ordinal(*plan, constraint.key), constraint.relation);
  }
  std::sort(
      relation_pairs.begin(),
      relation_pairs.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
  scenario->relation_template.source_ids.clear();
  scenario->relation_template.relations.clear();
  scenario->relation_template.source_ids.reserve(relation_pairs.size());
  scenario->relation_template.relations.reserve(relation_pairs.size());
  for (const auto &[source_id, relation] : relation_pairs) {
    if (!scenario->relation_template.source_ids.empty() &&
        scenario->relation_template.source_ids.back() == source_id) {
      scenario->relation_template.relations.back() = relation;
      continue;
    }
    scenario->relation_template.source_ids.push_back(source_id);
    scenario->relation_template.relations.push_back(relation);
  }
}

inline void compile_program_source_runtime_fields(ExactVariantPlan *plan) {
  auto &program = plan->lowered.program;

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    const auto source_id = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.onset_source_kind[pos]),
        program.onset_source_index[pos]);
    program.onset_source_ids[pos] = source_id;
    program.leaf_descriptors[pos].onset_source_id = source_id;
  }

  for (std::size_t i = 0; i < program.pool_member_indices.size(); ++i) {
    program.pool_member_source_ids[i] = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.pool_member_kind[i]),
        program.pool_member_indices[i]);
  }

  for (std::size_t i = 0; i < program.expr_source_index.size(); ++i) {
    program.expr_source_ids[i] = source_ordinal(
        *plan,
        static_cast<semantic::SourceKind>(program.expr_source_kind[i]),
        program.expr_source_index[i]);
  }
}

inline void append_runtime_span_values(
    const std::vector<semantic::Index> &arena,
    const ExactIndexSpan span,
    std::vector<semantic::Index> *values) {
  for (semantic::Index i = 0; i < span.size; ++i) {
    values->push_back(arena[static_cast<std::size_t>(span.offset + i)]);
  }
}

inline ExactRuntimeTruthFormula compile_runtime_readiness_cdf(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeTruthFormula formula;
  formula.empty_value = 1.0;
  append_runtime_span_values(
      plan.scenario_source_ids,
      scenario.before_source_span,
      &formula.product.source_cdf);
  append_runtime_span_values(
      plan.scenario_expr_ids,
      scenario.ready_expr_span,
      &formula.product.expr_cdf);
  formula.requires_scenario = !scenario.ready_expr_span.empty();
  return formula;
}

inline ExactRuntimeTruthFormula compile_runtime_after_survival(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeTruthFormula formula;
  formula.empty_value = 1.0;
  append_runtime_span_values(
      plan.scenario_source_ids,
      scenario.after_source_span,
      &formula.product.source_survival);
  append_runtime_span_values(
      plan.scenario_expr_ids,
      scenario.tail_expr_span,
      &formula.product.expr_survival);
  formula.requires_scenario = !scenario.tail_expr_span.empty();
  return formula;
}

inline ExactRuntimeTruthFormula compile_runtime_readiness_density(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeTruthFormula formula;
  formula.sum_of_products = true;
  formula.clean_signed = true;
  formula.empty_value = 0.0;
  formula.requires_scenario = !scenario.ready_expr_span.empty();

  for (semantic::Index i = 0; i < scenario.before_source_span.size; ++i) {
    const auto source_id = plan.scenario_source_ids[
        static_cast<std::size_t>(scenario.before_source_span.offset + i)];
    ExactRuntimeProductTerm term;
    term.factors.source_pdf.push_back(source_id);
    for (semantic::Index j = 0; j < scenario.before_source_span.size; ++j) {
      if (i == j) {
        continue;
      }
      term.factors.source_cdf.push_back(
          plan.scenario_source_ids[static_cast<std::size_t>(
              scenario.before_source_span.offset + j)]);
    }
    append_runtime_span_values(
        plan.scenario_expr_ids,
        scenario.ready_expr_span,
        &term.factors.expr_cdf);
    term.condition =
        runtime_product_term_condition(term, plan, scenario.source_order_facts);
    formula.sum_terms.push_back(std::move(term));
  }

  for (semantic::Index i = 0; i < scenario.ready_expr_span.size; ++i) {
    const auto expr_id = plan.scenario_expr_ids[
        static_cast<std::size_t>(scenario.ready_expr_span.offset + i)];
    ExactRuntimeProductTerm term;
    term.factors.expr_density.push_back(expr_id);
    append_runtime_span_values(
        plan.scenario_source_ids,
        scenario.before_source_span,
        &term.factors.source_cdf);
    for (semantic::Index j = 0; j < scenario.ready_expr_span.size; ++j) {
      if (i == j) {
        continue;
      }
      term.factors.expr_cdf.push_back(
          plan.scenario_expr_ids[static_cast<std::size_t>(
              scenario.ready_expr_span.offset + j)]);
    }
    term.condition =
        runtime_product_term_condition(term, plan, scenario.source_order_facts);
    formula.sum_terms.push_back(std::move(term));
  }

  return formula;
}

inline ExactRuntimeScenarioFormula compile_runtime_scenario_formula(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeScenarioFormula formula;
  formula.active_source_id = scenario.active_source_id;
  formula.relation_template = scenario.relation_template;
  formula.source_order_facts = scenario.source_order_facts;
  formula.has_readiness =
      !scenario.before_source_span.empty() || !scenario.ready_expr_span.empty();
  formula.readiness_cdf = compile_runtime_readiness_cdf(plan, scenario);
  formula.readiness_density = compile_runtime_readiness_density(plan, scenario);
  formula.after_survival = compile_runtime_after_survival(plan, scenario);
  formula.tail_condition = runtime_scenario_tail_condition(formula);
  return formula;
}

inline void compile_runtime_scenario_contextual_conditions(
    const ExactVariantPlan &plan,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    ExactRuntimeScenarioFormula *formula) {
  formula->tail_competitor_condition =
      filter_runtime_condition_for_competitors(
          formula->tail_condition,
          runtime_outcome,
          nullptr,
          plan);
  formula->has_tail_competitor_condition =
      !runtime_condition_empty(formula->tail_competitor_condition);
  formula->tail_competitor_subset_mask =
      runtime_competitor_subset_mask(
          formula->tail_competitor_condition,
          runtime_outcome,
          plan);
  formula->has_conditioned_readiness_terms = false;
  formula->condition_groups.clear();
  for (auto &term : formula->readiness_density.sum_terms) {
    term.condition_impossible =
        runtime_condition_order_contradiction(term.condition);
    term.tail_condition = filter_runtime_condition_for_truth(
        term.condition,
        formula->after_survival,
        plan);
    term.competitor_condition = filter_runtime_condition_for_competitors(
        term.condition,
        runtime_outcome,
        nullptr,
        plan);
    if (term.condition_impossible ||
        !runtime_condition_empty(term.tail_condition) ||
        !runtime_condition_empty(term.competitor_condition)) {
      formula->has_conditioned_readiness_terms = true;
    }
    term.context_group_index = semantic::kInvalidIndex;
    if (term.condition_impossible) {
      continue;
    }
    for (semantic::Index group_idx = 0;
         group_idx < static_cast<semantic::Index>(formula->condition_groups.size());
         ++group_idx) {
      const auto &group =
          formula->condition_groups[static_cast<std::size_t>(group_idx)];
      if (runtime_condition_equal(group.tail_condition, term.tail_condition) &&
          runtime_condition_equal(
              group.competitor_condition,
              term.competitor_condition)) {
        term.context_group_index = group_idx;
        break;
      }
    }
    if (term.context_group_index == semantic::kInvalidIndex) {
      term.context_group_index =
          static_cast<semantic::Index>(formula->condition_groups.size());
      formula->condition_groups.push_back(
          ExactRuntimeConditionGroup{
              term.tail_condition,
              term.competitor_condition,
              runtime_competitor_subset_mask(
                  term.competitor_condition,
                  runtime_outcome,
                  plan)});
    }
  }
}

inline ExactRuntimeVariantPlan compile_exact_runtime_plan(
    const ExactVariantPlan &plan) {
  ExactRuntimeVariantPlan runtime;
  runtime.outcomes.reserve(plan.outcomes.size());

  for (semantic::Index target_idx = 0;
       target_idx < static_cast<semantic::Index>(plan.outcomes.size());
       ++target_idx) {
    const auto target_pos = static_cast<std::size_t>(target_idx);
    const auto &outcome = plan.outcomes[target_pos];
    const auto &competitor_plan = plan.competitor_plans[target_pos];

    ExactRuntimeOutcomePlan runtime_outcome;
    runtime_outcome.scenarios.reserve(outcome.scenarios.size());
    for (const auto &scenario : outcome.scenarios) {
      runtime_outcome.scenarios.push_back(
          compile_runtime_scenario_formula(plan, scenario));
    }

    runtime_outcome.competitor_blocks.reserve(competitor_plan.blocks.size());
    for (const auto &block : competitor_plan.blocks) {
      ExactRuntimeCompetitorBlockPlan runtime_block;
      runtime_block.subsets.reserve(block.subsets.size());
      for (const auto &subset : block.subsets) {
        ExactRuntimeCompetitorSubsetPlan runtime_subset;
        runtime_subset.outcome_indices = subset.outcome_indices;
        runtime_subset.inclusion_sign = subset.inclusion_sign;
        runtime_subset.singleton_expr_root = subset.singleton_expr_root;
        runtime_subset.scenarios.reserve(subset.scenarios.size());
        for (const auto &scenario : subset.scenarios) {
          runtime_subset.scenarios.push_back(
              compile_runtime_scenario_formula(plan, scenario));
        }
        runtime_block.subsets.push_back(std::move(runtime_subset));
      }
      runtime_outcome.competitor_blocks.push_back(std::move(runtime_block));
    }

    for (auto &scenario : runtime_outcome.scenarios) {
      compile_runtime_scenario_contextual_conditions(
          plan,
          runtime_outcome,
          &scenario);
    }

    runtime_outcome.competitor_by_scenario.reserve(outcome.scenarios.size());
    for (const auto &target_scenario : outcome.scenarios) {
      ExactRuntimeScenarioCompetitorView scenario_view;
      scenario_view.blocks.reserve(competitor_plan.blocks.size());
      for (std::size_t block_idx = 0; block_idx < competitor_plan.blocks.size();
           ++block_idx) {
        const auto &block = competitor_plan.blocks[block_idx];
        ExactRuntimeScenarioBlockView block_view;
        block_view.subsets.reserve(block.subsets.size());
        for (const auto &subset : block.subsets) {
          ExactRuntimeScenarioSubsetView subset_view;
          for (semantic::Index scenario_idx = 0;
               scenario_idx < static_cast<semantic::Index>(subset.scenarios.size());
               ++scenario_idx) {
            const auto &scenario =
                subset.scenarios[static_cast<std::size_t>(scenario_idx)];
            if (scenario.active_source_id == target_scenario.active_source_id) {
              subset_view.same_active_scenario_indices.push_back(scenario_idx);
            }
          }
          block_view.subsets.push_back(std::move(subset_view));
        }
        scenario_view.blocks.push_back(std::move(block_view));
      }
      runtime_outcome.competitor_by_scenario.push_back(std::move(scenario_view));
    }

    runtime.outcomes.push_back(std::move(runtime_outcome));
  }

  return runtime;
}

inline ExactVariantPlan make_exact_variant_plan(
    const runtime::LoweredExactVariant &lowered,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_outcome_codes) {
  ExactVariantPlan plan;
  plan.lowered = lowered;
  plan.outcome_index_by_code.assign(
      n_outcome_codes + 1U,
      semantic::kInvalidIndex);
  ExactSupportBuilder builder(plan.lowered);
  plan.leaf_supports = builder.build_leaf_supports();
  plan.pool_supports = builder.build_pool_supports();
  plan.expr_supports = builder.build_expr_supports();
  plan.source_count = static_cast<semantic::Index>(
      plan.lowered.program.layout.n_leaves + plan.lowered.program.layout.n_pools);
  plan.leaf_source_ids.resize(static_cast<std::size_t>(plan.lowered.program.layout.n_leaves));
  plan.pool_source_ids.resize(static_cast<std::size_t>(plan.lowered.program.layout.n_pools));
  for (semantic::Index i = 0; i < plan.lowered.program.layout.n_leaves; ++i) {
    plan.leaf_source_ids[static_cast<std::size_t>(i)] = i;
  }
  for (semantic::Index i = 0; i < plan.lowered.program.layout.n_pools; ++i) {
    plan.pool_source_ids[static_cast<std::size_t>(i)] =
        static_cast<semantic::Index>(plan.lowered.program.layout.n_leaves + i);
  }
  compile_program_source_runtime_fields(&plan);

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
    outcome.scenarios = build_expr_transition_scenarios(plan, expr_root);
    const auto code_it =
        outcome_code_by_label.find(plan.lowered.outcome_labels[i]);
    if (code_it == outcome_code_by_label.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared outcome code for '" +
          plan.lowered.outcome_labels[i] + "'");
    }
    plan.outcome_index_by_code[static_cast<std::size_t>(code_it->second)] =
        static_cast<semantic::Index>(i);
    plan.outcomes.push_back(std::move(outcome));
  }

  plan.competitor_plans.reserve(plan.outcomes.size());
  for (semantic::Index target_outcome_idx = 0;
       target_outcome_idx < static_cast<semantic::Index>(plan.outcomes.size());
       ++target_outcome_idx) {
    plan.competitor_plans.push_back(
        build_target_competitor_plan(plan, target_outcome_idx));
  }
  plan.scenario_source_ids.clear();
  plan.scenario_expr_ids.clear();
  for (auto &outcome : plan.outcomes) {
    for (auto &scenario : outcome.scenarios) {
      compile_scenario_runtime_fields(&plan, &scenario);
    }
  }
  for (auto &target_plan : plan.competitor_plans) {
    for (auto &block : target_plan.blocks) {
      for (auto &subset : block.subsets) {
        for (auto &scenario : subset.scenarios) {
          compile_scenario_runtime_fields(&plan, &scenario);
        }
      }
    }
  }
  plan.runtime = compile_exact_runtime_plan(plan);
  plan.ranked_supported = variant_supports_ranked_sequence(plan);
  for (auto &outcome : plan.outcomes) {
    for (auto &scenario : outcome.scenarios) {
      release_scenario_planning_fields(&scenario);
    }
  }
  for (auto &target_plan : plan.competitor_plans) {
    for (auto &block : target_plan.blocks) {
      for (auto &subset : block.subsets) {
        for (auto &scenario : subset.scenarios) {
          release_scenario_planning_fields(&scenario);
        }
      }
    }
  }

  return plan;
}

} // namespace detail
} // namespace accumulatr::eval
