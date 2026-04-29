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

inline bool scenario_key_known_after_active(
    const ExactTransitionScenario &scenario,
    const ExactSourceKey key) {
  return std::find_if(
             scenario.after_keys.begin(),
             scenario.after_keys.end(),
             [&](const ExactSourceKey &existing) { return existing == key; }) !=
         scenario.after_keys.end();
}

inline bool scenario_key_known_before_key(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    const ExactSourceKey before_key,
    const ExactSourceKey after_key) {
  if (before_key == after_key) {
    return false;
  }
  if (scenario.active_key == after_key &&
      scenario_key_known_before_active(scenario, before_key)) {
    return true;
  }
  if (scenario.active_key == before_key &&
      scenario_key_known_after_active(scenario, after_key)) {
    return true;
  }
  const auto before_source_id = source_ordinal(plan, before_key);
  const auto after_source_id = source_ordinal(plan, after_key);
  if (before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex) {
    return false;
  }
  for (const auto &fact : scenario.source_order_facts) {
    if (fact.before_source_id == before_source_id &&
        fact.after_source_id == after_source_id) {
      return true;
    }
  }
  return false;
}

inline bool expr_certainly_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    const semantic::Index expr_idx);

inline bool expr_children_any_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    const ExactIndexSpan children) {
  const auto &program = plan.lowered.program;
  for (semantic::Index i = 0; i < children.size; ++i) {
    const auto child = program.expr_args[
        static_cast<std::size_t>(children.offset + i)];
    if (expr_certainly_not_completed_by_active(plan, scenario, child)) {
      return true;
    }
  }
  return false;
}

inline bool expr_children_all_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    const ExactIndexSpan children) {
  if (children.empty()) {
    return false;
  }
  const auto &program = plan.lowered.program;
  for (semantic::Index i = 0; i < children.size; ++i) {
    const auto child = program.expr_args[
        static_cast<std::size_t>(children.offset + i)];
    if (!expr_certainly_not_completed_by_active(plan, scenario, child)) {
      return false;
    }
  }
  return true;
}

inline bool expr_certainly_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    const semantic::Index expr_idx) {
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::Impossible:
    return true;
  case semantic::ExprKind::TrueExpr:
    return false;
  case semantic::ExprKind::Event:
    return scenario_key_known_after_active(
        scenario,
        ExactSourceKey{
            child_event_source_kind(program, expr_idx),
            child_event_source_index(program, expr_idx)});
  case semantic::ExprKind::And:
    return expr_children_any_not_completed_by_active(
        plan,
        scenario,
        ExactIndexSpan{
            program.expr_arg_offsets[pos],
            static_cast<semantic::Index>(
                program.expr_arg_offsets[pos + 1U] -
                program.expr_arg_offsets[pos])});
  case semantic::ExprKind::Or:
    return expr_children_all_not_completed_by_active(
        plan,
        scenario,
        ExactIndexSpan{
            program.expr_arg_offsets[pos],
            static_cast<semantic::Index>(
                program.expr_arg_offsets[pos + 1U] -
                program.expr_arg_offsets[pos])});
  case semantic::ExprKind::Guard: {
    const auto ref = program.expr_ref_child[pos];
    if (expr_certainly_not_completed_by_active(plan, scenario, ref)) {
      return true;
    }
    const auto blocker = program.expr_blocker_child[pos];
    if (!expr_is_simple_event(program, ref) ||
        !expr_is_simple_event(program, blocker)) {
      return false;
    }
    const auto ref_key =
        ExactSourceKey{
            child_event_source_kind(program, ref),
            child_event_source_index(program, ref)};
    const auto blocker_key =
        ExactSourceKey{
            child_event_source_kind(program, blocker),
            child_event_source_index(program, blocker)};
    return scenario_key_known_before_key(
        plan,
        scenario,
        blocker_key,
        ref_key);
  }
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
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
	          if (expr_certainly_not_completed_by_active(
	                  plan,
	                  scenario,
	                  child)) {
	            continue;
	          }
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

struct ExactCompetitorCandidate {
  semantic::Index expr_root{semantic::kInvalidIndex};
  semantic::Index outcome_index{semantic::kInvalidIndex};
};

inline void enumerate_competitor_subsets(
    const ExactVariantPlan &plan,
    const std::vector<ExactCompetitorCandidate> &block_competitors,
    const std::size_t start,
    std::vector<ExactCompetitorCandidate> *current,
    std::vector<ExactCompetitorSubsetPlan> *out) {
  for (std::size_t i = start; i < block_competitors.size(); ++i) {
    current->push_back(block_competitors[i]);
    ExactCompetitorSubsetPlan subset;
    subset.inclusion_sign = (current->size() % 2U == 1U) ? 1 : -1;
    subset.expr_roots.reserve(current->size());
    subset.outcome_indices.reserve(current->size());
    for (const auto &competitor : *current) {
      subset.expr_roots.push_back(competitor.expr_root);
      if (competitor.outcome_index != semantic::kInvalidIndex) {
        subset.outcome_indices.push_back(competitor.outcome_index);
      }
    }
    if (current->size() == 1U) {
      const auto &competitor = current->front();
      subset.singleton_expr_root = competitor.expr_root;
      if (competitor.outcome_index != semantic::kInvalidIndex &&
          competitor.outcome_index <
              static_cast<semantic::Index>(plan.outcomes.size()) &&
          plan.outcomes[static_cast<std::size_t>(competitor.outcome_index)]
                  .expr_root == competitor.expr_root) {
        subset.scenarios =
            plan.outcomes[static_cast<std::size_t>(competitor.outcome_index)]
                .scenarios;
      } else {
        subset.scenarios =
            build_expr_transition_scenarios(plan, competitor.expr_root);
      }
    } else {
      subset.scenarios =
          build_expr_conjunction_scenarios(plan, subset.expr_roots);
    }
    if (!subset.scenarios.empty()) {
      out->push_back(std::move(subset));
    }
    enumerate_competitor_subsets(
        plan, block_competitors, i + 1U, current, out);
    current->pop_back();
  }
}

inline std::vector<ExactCompetitorCandidate> planned_competitor_candidates(
    const ExactVariantPlan &plan,
    const semantic::Index target_outcome_idx) {
  const auto &program = plan.lowered.program;
  if (target_outcome_idx < 0 ||
      target_outcome_idx + 1 >=
          static_cast<semantic::Index>(program.outcome_competitor_offsets.size())) {
    throw std::runtime_error("exact competitor plan has invalid target outcome");
  }

  const auto begin =
      program.outcome_competitor_offsets[static_cast<std::size_t>(target_outcome_idx)];
  const auto end =
      program.outcome_competitor_offsets[static_cast<std::size_t>(target_outcome_idx + 1)];
  std::vector<ExactCompetitorCandidate> competitors;
  competitors.reserve(static_cast<std::size_t>(end - begin));
  for (semantic::Index pos = begin; pos < end; ++pos) {
    const auto idx = static_cast<std::size_t>(pos);
    const auto expr_root = program.outcome_competitor_expr_roots[idx];
    if (expr_root == semantic::kInvalidIndex) {
      continue;
    }
    ExactCompetitorCandidate candidate;
    candidate.expr_root = expr_root;
    candidate.outcome_index =
        idx < program.outcome_competitor_indices.size()
            ? program.outcome_competitor_indices[idx]
            : semantic::kInvalidIndex;
    const bool duplicate = std::any_of(
        competitors.begin(),
        competitors.end(),
        [&](const ExactCompetitorCandidate &existing) {
          return existing.expr_root == candidate.expr_root &&
                 existing.outcome_index == candidate.outcome_index;
        });
    if (!duplicate) {
      competitors.push_back(candidate);
    }
  }
  return competitors;
}

inline std::vector<std::vector<ExactCompetitorCandidate>>
build_competitor_overlap_blocks(
    const ExactVariantPlan &plan,
    const semantic::Index target_outcome_idx) {
  const auto competitors = planned_competitor_candidates(plan, target_outcome_idx);
  std::vector<std::vector<ExactCompetitorCandidate>> blocks;
  std::vector<std::uint8_t> visited(competitors.size(), 0U);
  for (std::size_t start = 0; start < competitors.size(); ++start) {
    if (visited[start] != 0U) {
      continue;
    }
    std::vector<std::size_t> stack{start};
    visited[start] = 1U;
    std::vector<ExactCompetitorCandidate> block;
    while (!stack.empty()) {
      const auto idx = stack.back();
      stack.pop_back();
      const auto competitor = competitors[idx];
      block.push_back(competitor);
      for (std::size_t other = 0; other < competitors.size(); ++other) {
        if (visited[other] != 0U) {
          continue;
        }
        if (!supports_overlap(
                plan.expr_supports[static_cast<std::size_t>(
                    competitor.expr_root)],
                plan.expr_supports[static_cast<std::size_t>(
                    competitors[other].expr_root)])) {
          continue;
        }
        visited[other] = 1U;
        stack.push_back(other);
      }
    }
    std::sort(
        block.begin(),
        block.end(),
        [](const ExactCompetitorCandidate &lhs,
           const ExactCompetitorCandidate &rhs) {
          if (lhs.expr_root != rhs.expr_root) {
            return lhs.expr_root < rhs.expr_root;
          }
          return lhs.outcome_index < rhs.outcome_index;
        });
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
  for (const auto &block_competitors : blocks) {
    ExactCompetitorBlockPlan block;
    std::vector<ExactCompetitorCandidate> current;
    enumerate_competitor_subsets(
        plan, block_competitors, 0U, &current, &block.subsets);
    target_plan.blocks.push_back(std::move(block));
  }
  return target_plan;
}

inline bool scenario_supports_ranked_sequence(
    const ExactTransitionScenario &scenario) {
  (void)scenario;
  return true;
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

inline bool scenario_is_terminal_no_response_source(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario,
    semantic::Index *source_id) {
  if (scenario.active_source_id == semantic::kInvalidIndex ||
      !scenario.before_source_span.empty() ||
      !scenario.after_source_span.empty() ||
      !scenario.ready_expr_span.empty() ||
      !scenario.tail_expr_span.empty() ||
      scenario.factors.size() != 1U ||
      scenario.factors.front().kind != ExactFactorKind::AtPdf) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(scenario.active_source_id);
  if (pos >= plan.source_kernels.size() ||
      plan.source_kernels[pos].kind !=
          CompiledSourceChannelKernelKind::LeafAbsolute) {
    return false;
  }
  if (source_id != nullptr) {
    *source_id = scenario.active_source_id;
  }
  return true;
}

inline ExactNoResponsePlan compile_no_response_plan(
    const ExactVariantPlan &plan) {
  ExactNoResponsePlan no_response;
  if (plan.source_count <= 0 || plan.outcomes.empty()) {
    return no_response;
  }

  std::vector<std::uint8_t> covered(
      static_cast<std::size_t>(plan.source_count), 0U);
  for (const auto &outcome : plan.outcomes) {
    if (outcome.scenarios.size() != 1U) {
      return ExactNoResponsePlan{};
    }
    semantic::Index source_id{semantic::kInvalidIndex};
    if (!scenario_is_terminal_no_response_source(
            plan, outcome.scenarios.front(), &source_id)) {
      return ExactNoResponsePlan{};
    }
    const auto pos = static_cast<std::size_t>(source_id);
    if (pos >= covered.size() || covered[pos] != 0U) {
      return ExactNoResponsePlan{};
    }
    covered[pos] = 1U;
  }

  no_response.source_ids.reserve(covered.size());
  for (std::size_t i = 0; i < covered.size(); ++i) {
    if (covered[i] == 0U) {
      return ExactNoResponsePlan{};
    }
    no_response.source_ids.push_back(static_cast<semantic::Index>(i));
  }
  no_response.terminal_leaf_survival_product = true;
  return no_response;
}

inline ExactSimpleRacePlan compile_simple_race_plan(
    const ExactVariantPlan &plan) {
  ExactSimpleRacePlan simple;
  if (!plan.no_response.terminal_leaf_survival_product ||
      plan.outcomes.empty()) {
    return simple;
  }

  simple.source_ids = plan.no_response.source_ids;
  simple.outcome_source_ids.assign(
      plan.outcomes.size(),
      semantic::kInvalidIndex);
  simple.source_leaf_indices.assign(
      static_cast<std::size_t>(plan.source_count),
      semantic::kInvalidIndex);

  for (const auto source_id : simple.source_ids) {
    const auto source_pos = static_cast<std::size_t>(source_id);
    if (source_id == semantic::kInvalidIndex ||
        source_pos >= plan.source_kernels.size()) {
      return ExactSimpleRacePlan{};
    }
    const auto &kernel = plan.source_kernels[source_pos];
    if (kernel.kind != CompiledSourceChannelKernelKind::LeafAbsolute ||
        kernel.leaf_index == semantic::kInvalidIndex) {
      return ExactSimpleRacePlan{};
    }
    simple.source_leaf_indices[source_pos] = kernel.leaf_index;
  }

  for (std::size_t outcome_index = 0; outcome_index < plan.outcomes.size();
       ++outcome_index) {
    if (plan.outcomes[outcome_index].scenarios.size() != 1U) {
      return ExactSimpleRacePlan{};
    }
    semantic::Index source_id{semantic::kInvalidIndex};
    if (!scenario_is_terminal_no_response_source(
            plan,
            plan.outcomes[outcome_index].scenarios.front(),
            &source_id)) {
      return ExactSimpleRacePlan{};
    }
    const auto source_pos = static_cast<std::size_t>(source_id);
    if (source_pos >= simple.source_leaf_indices.size() ||
        simple.source_leaf_indices[source_pos] == semantic::kInvalidIndex) {
      return ExactSimpleRacePlan{};
    }
    simple.outcome_source_ids[outcome_index] = source_id;
  }

  simple.terminal_leaf_top1 = true;
  return simple;
}

inline semantic::Index append_exact_probability_op(
    ExactProbabilityProgram *program,
    ExactProbabilityOp op,
    const std::vector<semantic::Index> &children = {}) {
  op.children.offset =
      static_cast<semantic::Index>(program->child_ops.size());
  op.children.size = static_cast<semantic::Index>(children.size());
  program->child_ops.insert(
      program->child_ops.end(), children.begin(), children.end());
  const auto op_index = static_cast<semantic::Index>(program->ops.size());
  program->ops.push_back(op);
  return op_index;
}

inline ExactIndexSpan append_exact_probability_sources(
    ExactProbabilityProgramSet *programs,
    const std::vector<semantic::Index> &source_ids) {
  const auto offset =
      static_cast<semantic::Index>(programs->source_ids.size());
  programs->source_ids.insert(
      programs->source_ids.end(), source_ids.begin(), source_ids.end());
  return ExactIndexSpan{
      offset,
      static_cast<semantic::Index>(source_ids.size())};
}

inline semantic::Index append_exact_probability_program(
    ExactProbabilityProgramSet *programs,
    ExactProbabilityProgram program) {
  const auto program_index =
      static_cast<semantic::Index>(programs->programs.size());
  programs->programs.push_back(std::move(program));
  return program_index;
}

inline ExactProbabilityProgram make_top1_leaf_race_density_program(
    ExactProbabilityProgramSet *programs,
    const semantic::Index outcome_index,
    const semantic::Index target_source_id,
    const std::vector<semantic::Index> &source_ids) {
  ExactProbabilityProgram program;
  program.value_kind = ExactProbabilityValueKind::Density;
  program.requires_trigger_enumeration = false;

  ExactProbabilityOp density;
  density.kind = ExactProbabilityOpKind::Top1LeafRaceDensity;
  density.value_kind = ExactProbabilityValueKind::Density;
  density.outcome_index = outcome_index;
  density.target_source_id = target_source_id;
  density.time_slot =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  density.source_span =
      append_exact_probability_sources(programs, source_ids);
  const auto density_op =
      append_exact_probability_op(&program, density);

  ExactProbabilityOp root;
  root.kind = ExactProbabilityOpKind::WeightedTriggerSum;
  root.value_kind = ExactProbabilityValueKind::Density;
  program.root_child = density_op;
  program.root = append_exact_probability_op(&program, root, {density_op});
  return program;
}

inline ExactProbabilityProgram make_terminal_no_response_probability_program(
    ExactProbabilityProgramSet *programs,
    const std::vector<semantic::Index> &source_ids) {
  ExactProbabilityProgram program;
  program.value_kind = ExactProbabilityValueKind::Probability;
  program.requires_trigger_enumeration = false;

  ExactProbabilityOp no_response;
  no_response.kind = ExactProbabilityOpKind::TerminalNoResponseProbability;
  no_response.value_kind = ExactProbabilityValueKind::Probability;
  no_response.source_span =
      append_exact_probability_sources(programs, source_ids);
  const auto no_response_op =
      append_exact_probability_op(&program, no_response);

  ExactProbabilityOp root;
  root.kind = ExactProbabilityOpKind::WeightedTriggerSum;
  root.value_kind = ExactProbabilityValueKind::Probability;
  program.root_child = no_response_op;
  program.root = append_exact_probability_op(&program, root, {no_response_op});
  return program;
}

inline ExactProbabilityProgram make_generic_transition_density_program(
    const semantic::Index outcome_index) {
  ExactProbabilityProgram program;
  program.value_kind = ExactProbabilityValueKind::Density;
  program.requires_trigger_enumeration = false;

  ExactProbabilityOp density;
  density.kind = ExactProbabilityOpKind::GenericTransitionDensity;
  density.value_kind = ExactProbabilityValueKind::Density;
  density.outcome_index = outcome_index;
  density.time_slot =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto density_op =
      append_exact_probability_op(&program, density);

  ExactProbabilityOp root;
  root.kind = ExactProbabilityOpKind::WeightedTriggerSum;
  root.value_kind = ExactProbabilityValueKind::Density;
  program.root_child = density_op;
  program.root = append_exact_probability_op(&program, root, {density_op});
  return program;
}

inline ExactProbabilityProgram make_top1_leaf_race_probability_program(
    ExactProbabilityProgramSet *programs,
    const semantic::Index outcome_index,
    const semantic::Index target_source_id,
    const std::vector<semantic::Index> &source_ids) {
  ExactProbabilityProgram program;
  program.value_kind = ExactProbabilityValueKind::Probability;
  program.requires_trigger_enumeration = false;

  ExactProbabilityOp density;
  density.kind = ExactProbabilityOpKind::Top1LeafRaceDensity;
  density.value_kind = ExactProbabilityValueKind::Density;
  density.outcome_index = outcome_index;
  density.target_source_id = target_source_id;
  density.time_slot =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  density.source_span =
      append_exact_probability_sources(programs, source_ids);
  const auto density_op =
      append_exact_probability_op(&program, density);

  ExactProbabilityOp weighted;
  weighted.kind = ExactProbabilityOpKind::WeightedTriggerSum;
  weighted.value_kind = ExactProbabilityValueKind::Density;
  const auto weighted_op =
      append_exact_probability_op(&program, weighted, {density_op});

  ExactProbabilityOp integral;
  integral.kind = ExactProbabilityOpKind::Integral;
  integral.value_kind = ExactProbabilityValueKind::Probability;
  integral.outcome_index = outcome_index;
  const auto integral_op =
      append_exact_probability_op(&program, integral, {weighted_op});

  program.root_child = density_op;
  program.root = integral_op;
  return program;
}

inline ExactProbabilityProgram make_generic_transition_probability_program(
    const semantic::Index outcome_index) {
  ExactProbabilityProgram program;
  program.value_kind = ExactProbabilityValueKind::Probability;
  program.requires_trigger_enumeration = false;

  ExactProbabilityOp density;
  density.kind = ExactProbabilityOpKind::GenericTransitionDensity;
  density.value_kind = ExactProbabilityValueKind::Density;
  density.outcome_index = outcome_index;
  density.time_slot =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto density_op =
      append_exact_probability_op(&program, density);

  ExactProbabilityOp weighted;
  weighted.kind = ExactProbabilityOpKind::WeightedTriggerSum;
  weighted.value_kind = ExactProbabilityValueKind::Density;
  const auto weighted_op =
      append_exact_probability_op(&program, weighted, {density_op});

  ExactProbabilityOp integral;
  integral.kind = ExactProbabilityOpKind::Integral;
  integral.value_kind = ExactProbabilityValueKind::Probability;
  integral.outcome_index = outcome_index;
  const auto integral_op =
      append_exact_probability_op(&program, integral, {weighted_op});

  program.root_child = density_op;
  program.root = integral_op;
  return program;
}

inline ExactProbabilityProgramSet compile_exact_probability_programs(
    const ExactVariantPlan &plan) {
  ExactProbabilityProgramSet programs;
  programs.density_by_outcome.assign(
      plan.outcomes.size(),
      semantic::kInvalidIndex);
  programs.finite_probability_by_outcome.assign(
      plan.outcomes.size(),
      semantic::kInvalidIndex);

  for (semantic::Index outcome_index = 0;
       outcome_index < static_cast<semantic::Index>(plan.outcomes.size());
       ++outcome_index) {
    if (plan.simple_race.terminal_leaf_top1) {
      const auto target_source_id =
          plan.simple_race.outcome_source_ids[
              static_cast<std::size_t>(outcome_index)];
      ExactProbabilityProgram program =
          make_top1_leaf_race_density_program(
              &programs,
              outcome_index,
              target_source_id,
              plan.simple_race.source_ids);
      program.requires_trigger_enumeration =
          !plan.shared_trigger_indices.empty();
      programs.density_by_outcome[
          static_cast<std::size_t>(outcome_index)] =
              append_exact_probability_program(&programs, std::move(program));
      ExactProbabilityProgram probability_program =
          make_top1_leaf_race_probability_program(
              &programs,
              outcome_index,
              target_source_id,
              plan.simple_race.source_ids);
      probability_program.requires_trigger_enumeration =
          !plan.shared_trigger_indices.empty();
      programs.finite_probability_by_outcome[
          static_cast<std::size_t>(outcome_index)] =
              append_exact_probability_program(
                  &programs,
                  std::move(probability_program));
      continue;
    }
    ExactProbabilityProgram program =
        make_generic_transition_density_program(outcome_index);
    program.requires_trigger_enumeration =
        !plan.shared_trigger_indices.empty();
    programs.density_by_outcome[
        static_cast<std::size_t>(outcome_index)] =
            append_exact_probability_program(&programs, std::move(program));
    ExactProbabilityProgram probability_program =
        make_generic_transition_probability_program(outcome_index);
    probability_program.requires_trigger_enumeration =
        !plan.shared_trigger_indices.empty();
    programs.finite_probability_by_outcome[
        static_cast<std::size_t>(outcome_index)] =
            append_exact_probability_program(
                &programs,
                std::move(probability_program));
  }

  if (plan.no_response.terminal_leaf_survival_product) {
    ExactProbabilityProgram program =
        make_terminal_no_response_probability_program(
            &programs,
            plan.no_response.source_ids);
    program.requires_trigger_enumeration =
        !plan.shared_trigger_indices.empty();
    programs.no_response_probability =
        append_exact_probability_program(&programs, std::move(program));
  }

  return programs;
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

inline void compile_source_kernels(ExactVariantPlan *plan) {
  const auto &program = plan->lowered.program;
  plan->source_kernels.assign(
      static_cast<std::size_t>(plan->source_count), ExactSourceKernel{});

  for (semantic::Index i = 0; i < program.layout.n_leaves; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    auto &kernel = plan->source_kernels[pos];
    kernel.source_id = i;
    kernel.leaf_index = i;
    kernel.onset_source_id = program.onset_source_ids[pos];
    kernel.kind = static_cast<semantic::OnsetKind>(
                      program.onset_kind[pos]) ==
                          semantic::OnsetKind::Absolute
                      ? CompiledSourceChannelKernelKind::LeafAbsolute
                      : CompiledSourceChannelKernelKind::LeafOnsetConvolution;
  }

  for (semantic::Index i = 0; i < program.layout.n_pools; ++i) {
    const auto source_id =
        static_cast<semantic::Index>(program.layout.n_leaves + i);
    const auto pos = static_cast<std::size_t>(i);
    auto &kernel =
        plan->source_kernels[static_cast<std::size_t>(source_id)];
    kernel.kind = CompiledSourceChannelKernelKind::PoolKOfN;
    kernel.source_id = source_id;
    kernel.pool_index = i;
    kernel.pool_member_offset = program.pool_member_offsets[pos];
    kernel.pool_member_count =
        program.pool_member_offsets[pos + 1U] -
        program.pool_member_offsets[pos];
    kernel.pool_k = program.pool_k[pos];
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

inline void append_expr_union_subset(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &children,
    const std::vector<semantic::Index> &child_positions) {
  const auto offset = static_cast<semantic::Index>(
      plan->expr_union_subset_children.size());
  plan->expr_union_subset_children.insert(
      plan->expr_union_subset_children.end(),
      children.begin(),
      children.end());
  const auto position_offset = static_cast<semantic::Index>(
      plan->expr_union_subset_child_positions.size());
  plan->expr_union_subset_child_positions.insert(
      plan->expr_union_subset_child_positions.end(),
      child_positions.begin(),
      child_positions.end());
  plan->expr_union_subsets.push_back(
      ExactExprUnionSubset{
          ExactIndexSpan{
              offset,
              static_cast<semantic::Index>(children.size())},
          ExactIndexSpan{
              position_offset,
              static_cast<semantic::Index>(child_positions.size())},
          children.size() % 2U == 1U ? 1 : -1});
}

inline void enumerate_expr_union_subsets(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &children,
    const std::size_t start,
    std::vector<semantic::Index> *current,
    std::vector<semantic::Index> *current_positions) {
  for (std::size_t i = start; i < children.size(); ++i) {
    current->push_back(children[i]);
    current_positions->push_back(static_cast<semantic::Index>(i));
    append_expr_union_subset(plan, *current, *current_positions);
    enumerate_expr_union_subsets(
        plan, children, i + 1U, current, current_positions);
    current_positions->pop_back();
    current->pop_back();
  }
}

inline ExactIndexSpan compile_expr_union_subsets(
    ExactVariantPlan *plan,
    const ExactIndexSpan children_span) {
  const auto offset =
      static_cast<semantic::Index>(plan->expr_union_subsets.size());
  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(children_span.size));
  const auto &program = plan->lowered.program;
  for (semantic::Index i = 0; i < children_span.size; ++i) {
    children.push_back(
        program.expr_args[static_cast<std::size_t>(children_span.offset + i)]);
  }
  std::vector<semantic::Index> current;
  current.reserve(children.size());
  std::vector<semantic::Index> current_positions;
  current_positions.reserve(children.size());
  enumerate_expr_union_subsets(
      plan, children, 0U, &current, &current_positions);
  return ExactIndexSpan{
      offset,
      static_cast<semantic::Index>(
          plan->expr_union_subsets.size() - static_cast<std::size_t>(offset))};
}

inline bool expr_children_supports_overlap(
    const ExactVariantPlan &plan,
    const ExactIndexSpan children_span) {
  const auto &program = plan.lowered.program;
  for (semantic::Index i = 0; i < children_span.size; ++i) {
    const auto lhs_idx =
        program.expr_args[static_cast<std::size_t>(children_span.offset + i)];
    const auto &lhs = plan.expr_supports[static_cast<std::size_t>(lhs_idx)];
    for (semantic::Index j = i + 1; j < children_span.size; ++j) {
      const auto rhs_idx =
          program.expr_args[static_cast<std::size_t>(children_span.offset + j)];
      const auto &rhs = plan.expr_supports[static_cast<std::size_t>(rhs_idx)];
      if (supports_overlap(lhs, rhs)) {
        return true;
      }
    }
  }
  return false;
}

inline bool expr_static_subset_of_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id);

inline bool expr_static_can_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id);

inline bool expr_span_static_has_event_source_subset(
    const ExactVariantPlan &plan,
    const ExactIndexSpan span,
    const semantic::Index source_id) {
  const auto &program = plan.lowered.program;
  for (semantic::Index i = span.offset; i < span.offset + span.size; ++i) {
    if (expr_static_subset_of_event_source(
            plan,
            program.expr_args[static_cast<std::size_t>(i)],
            source_id)) {
      return true;
    }
  }
  return false;
}

inline bool expr_static_subset_of_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::Event:
    return program.expr_source_ids[pos] == source_id;
  case semantic::ExprKind::And:
    return expr_span_static_has_event_source_subset(
        plan,
        ExactIndexSpan{
            program.expr_arg_offsets[pos],
            static_cast<semantic::Index>(
                program.expr_arg_offsets[pos + 1U] -
                program.expr_arg_offsets[pos])},
        source_id);
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_static_subset_of_event_source(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              source_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Guard:
    return expr_static_subset_of_event_source(
        plan,
        program.expr_ref_child[pos],
        source_id);
  case semantic::ExprKind::Impossible:
    return true;
  case semantic::ExprKind::TrueExpr:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_static_can_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx,
    const semantic::Index source_id) {
  if (expr_idx == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_idx);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    const auto &kernel = plan.expr_kernels[pos];
    return kernel.simple_event_guard &&
           !kernel.has_unless &&
           kernel.guard_ref_source_id == source_id;
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    if (expr_static_can_equiv_event_source(
            plan,
            program.expr_args[static_cast<std::size_t>(i)],
            source_id)) {
      return true;
    }
  }
  return false;
}

inline ExactIndexSpan compile_expr_union_absorption_candidates(
    ExactVariantPlan *plan,
    const semantic::Index expr_idx,
    const ExactExprKernel &kernel) {
  const auto offset = static_cast<semantic::Index>(
      plan->expr_union_absorption_sources.size());
  const auto &program = plan->lowered.program;
  const auto &support = plan->expr_supports[static_cast<std::size_t>(expr_idx)];
  for (const auto source_id : support) {
    bool all_subset = true;
    bool any_equiv = false;
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      const auto child = program.expr_args[
          static_cast<std::size_t>(kernel.children.offset + i)];
      if (!expr_static_subset_of_event_source(*plan, child, source_id)) {
        all_subset = false;
        break;
      }
      any_equiv =
          any_equiv ||
          expr_static_can_equiv_event_source(*plan, child, source_id);
    }
    if (all_subset && any_equiv) {
      plan->expr_union_absorption_sources.push_back(source_id);
    }
  }
  return ExactIndexSpan{
      offset,
      static_cast<semantic::Index>(
          plan->expr_union_absorption_sources.size() -
          static_cast<std::size_t>(offset))};
}

inline void compile_exact_expr_kernels(ExactVariantPlan *plan) {
  const auto &program = plan->lowered.program;
  plan->expr_kernels.assign(program.expr_kind.size(), ExactExprKernel{});
  plan->expr_union_subsets.clear();
  plan->expr_union_subset_children.clear();
  plan->expr_union_subset_child_positions.clear();
  plan->expr_union_absorption_sources.clear();
  plan->expr_union_kernel_cache_slot_count = 0;

  for (semantic::Index expr_idx = 0;
       expr_idx < static_cast<semantic::Index>(program.expr_kind.size());
       ++expr_idx) {
    const auto pos = static_cast<std::size_t>(expr_idx);
    auto &kernel = plan->expr_kernels[pos];
    kernel.kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
    kernel.children = ExactIndexSpan{
        program.expr_arg_offsets[pos],
        static_cast<semantic::Index>(
            program.expr_arg_offsets[pos + 1U] -
            program.expr_arg_offsets[pos])};

    if (kernel.kind == semantic::ExprKind::Event) {
      kernel.event_source_id = program.expr_source_ids[pos];
      continue;
    }

    if (kernel.kind == semantic::ExprKind::Or) {
      kernel.children_overlap =
          expr_children_supports_overlap(*plan, kernel.children);
      if (kernel.children_overlap) {
        kernel.union_subset_span =
            compile_expr_union_subsets(plan, kernel.children);
        kernel.union_absorption_candidate_span =
            compile_expr_union_absorption_candidates(
                plan, expr_idx, kernel);
        kernel.union_kernel_cache_slot =
            plan->expr_union_kernel_cache_slot_count++;
      }
      continue;
    }

    if (kernel.kind != semantic::ExprKind::Guard) {
      continue;
    }

    kernel.guard_ref_expr_id = program.expr_ref_child[pos];
    kernel.guard_blocker_expr_id = program.expr_blocker_child[pos];
    kernel.has_unless = !kernel.children.empty();
    const auto ref_pos =
        static_cast<std::size_t>(kernel.guard_ref_expr_id);
    const auto blocker_pos =
        static_cast<std::size_t>(kernel.guard_blocker_expr_id);
    kernel.guard_ref_kind =
        static_cast<semantic::ExprKind>(program.expr_kind[ref_pos]);
    kernel.guard_blocker_kind =
        static_cast<semantic::ExprKind>(program.expr_kind[blocker_pos]);
    kernel.simple_event_guard =
        kernel.guard_ref_kind == semantic::ExprKind::Event &&
        kernel.guard_blocker_kind == semantic::ExprKind::Event;
    if (kernel.simple_event_guard) {
      kernel.guard_ref_source_id = program.expr_source_ids[ref_pos];
      kernel.guard_blocker_source_id = program.expr_source_ids[blocker_pos];
    }
    if (kernel.guard_ref_kind == semantic::ExprKind::And) {
      kernel.guard_ref_child_span = ExactIndexSpan{
          program.expr_arg_offsets[ref_pos],
          static_cast<semantic::Index>(
              program.expr_arg_offsets[ref_pos + 1U] -
              program.expr_arg_offsets[ref_pos])};
    }
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
  return formula;
}

inline ExactRuntimeTruthFormula compile_runtime_readiness_density(
    const ExactVariantPlan &plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeTruthFormula formula;
  formula.sum_of_products = true;
  formula.clean_signed = true;
  formula.empty_value = 0.0;

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

inline bool compile_runtime_factors_empty(const ExactRuntimeFactors &factors) {
  return factors.source_pdf.empty() &&
         factors.source_cdf.empty() &&
         factors.source_survival.empty() &&
         factors.expr_density.empty() &&
         factors.expr_cdf.empty() &&
         factors.expr_survival.empty();
}

inline semantic::Index compile_relation_condition_id(
    ExactVariantPlan *plan,
    const ExactRelationTemplate &relation_template) {
  CompiledMathConditionKey key;
  key.source_ids = relation_template.source_ids;
  key.relations.reserve(relation_template.relations.size());
  for (const auto relation : relation_template.relations) {
    key.relations.push_back(static_cast<std::uint8_t>(relation));
  }
  return compiled_math_intern_condition(&plan->compiled_math, std::move(key));
}

inline bool relation_template_equal(const ExactRelationTemplate &lhs,
                                    const ExactRelationTemplate &rhs) {
  return lhs.source_ids == rhs.source_ids && lhs.relations == rhs.relations;
}

inline semantic::Index compile_source_view_id(
    ExactVariantPlan *plan,
    const ExactRelationTemplate &relation_template) {
  if (relation_template.empty()) {
    return 0;
  }
  for (semantic::Index i = 0;
       i < static_cast<semantic::Index>(plan->compiled_source_views.size());
       ++i) {
    if (relation_template_equal(
            plan->compiled_source_views[static_cast<std::size_t>(i)],
            relation_template)) {
      return i + 1U;
    }
  }
  plan->compiled_source_views.push_back(relation_template);
  return static_cast<semantic::Index>(plan->compiled_source_views.size());
}

inline semantic::Index compile_expr_value_node(
    ExactVariantPlan *plan,
    semantic::Index expr_id,
    CompiledMathNodeKind value_kind,
    semantic::Index condition_id,
    semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    semantic::Index source_view_id = 0);

struct CompiledConditionNormalizer {
  semantic::Index node_id{semantic::kInvalidIndex};
  semantic::Index root_id{semantic::kInvalidIndex};
};

inline CompiledConditionNormalizer compile_condition_normalizer(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathTimeSlot time_slot) {
  if (expr_id == semantic::kInvalidIndex) {
    return {};
  }
  const auto normalizer_node =
      compile_expr_value_node(
          plan,
          expr_id,
          CompiledMathNodeKind::ExprCdf,
          0,
          static_cast<semantic::Index>(time_slot));
  return CompiledConditionNormalizer{
      normalizer_node,
      compiled_math_make_root(&plan->compiled_math, normalizer_node)};
}

inline void append_compiled_condition_fact(
    CompiledMathConditionKey *key,
    const CompiledMathConditionFactKind kind,
    const semantic::Index subject_id,
    const semantic::Index aux_id,
    const semantic::Index aux2_id,
    const CompiledMathTimeSlot time_slot,
    const CompiledConditionNormalizer normalizer = {}) {
  if (subject_id == semantic::kInvalidIndex) {
    return;
  }
  key->fact_kinds.push_back(static_cast<std::uint8_t>(kind));
  key->fact_subject_ids.push_back(subject_id);
  key->fact_aux_ids.push_back(aux_id);
  key->fact_aux2_ids.push_back(aux2_id);
  key->fact_time_ids.push_back(static_cast<semantic::Index>(time_slot));
  key->fact_normalizer_node_ids.push_back(normalizer.node_id);
  key->fact_normalizer_root_ids.push_back(normalizer.root_id);
}

inline void append_runtime_condition_to_compiled_key(
    ExactVariantPlan *plan,
    CompiledMathConditionKey *key,
    const ExactRuntimeTermCondition &condition,
    const CompiledMathTimeSlot time_slot) {
  if (condition.exact_source_id != semantic::kInvalidIndex) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::SourceExact,
        condition.exact_source_id,
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        time_slot);
  }
  for (const auto source_id : condition.upper_bound_source_ids) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::SourceUpperBound,
        source_id,
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        time_slot);
  }
  for (const auto source_id : condition.lower_bound_source_ids) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::SourceLowerBound,
        source_id,
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        time_slot);
  }
  for (const auto expr_id : condition.upper_bound_expr_ids) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::ExprUpperBound,
        expr_id,
        semantic::kInvalidIndex,
        semantic::kInvalidIndex,
        time_slot,
        compile_condition_normalizer(plan, expr_id, time_slot));
  }
  for (const auto &fact : condition.source_order_facts) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::SourceOrder,
        fact.before_source_id,
        fact.after_source_id,
        semantic::kInvalidIndex,
        time_slot);
  }
  for (const auto &fact : condition.guard_upper_bound_facts) {
    append_compiled_condition_fact(
        key,
        CompiledMathConditionFactKind::GuardUpperBound,
        fact.expr_id,
        fact.ref_source_id,
        fact.blocker_source_id,
        time_slot,
        compile_condition_normalizer(plan, fact.expr_id, time_slot));
  }
}

inline semantic::Index compile_runtime_condition_id(
    ExactVariantPlan *plan,
    const ExactRuntimeTermCondition &condition,
    const CompiledMathTimeSlot time_slot) {
  if (runtime_condition_empty(condition)) {
    return 0;
  }
  CompiledMathConditionKey key;
  key.impossible = runtime_condition_order_contradiction(condition);
  append_runtime_condition_to_compiled_key(plan, &key, condition, time_slot);
  return compiled_math_intern_condition(&plan->compiled_math, std::move(key));
}

inline semantic::Index compile_combined_runtime_condition_id(
    ExactVariantPlan *plan,
    const ExactRuntimeTermCondition &readiness_condition,
    const ExactRuntimeTermCondition &observed_condition) {
  if (runtime_condition_empty(readiness_condition) &&
      runtime_condition_empty(observed_condition)) {
    return 0;
  }
  CompiledMathConditionKey key;
  key.impossible =
      runtime_condition_order_contradiction(readiness_condition) ||
      runtime_condition_order_contradiction(observed_condition);
  append_runtime_condition_to_compiled_key(
      plan, &key, readiness_condition, CompiledMathTimeSlot::Readiness);
  append_runtime_condition_to_compiled_key(
      plan, &key, observed_condition, CompiledMathTimeSlot::Observed);
  return compiled_math_intern_condition(&plan->compiled_math, std::move(key));
}

inline void compile_runtime_condition_ids(
    ExactVariantPlan *plan,
    ExactRuntimeTermCondition *condition) {
  condition->observed_condition_id =
      compile_runtime_condition_id(
          plan, *condition, CompiledMathTimeSlot::Observed);
  condition->readiness_condition_id =
      compile_runtime_condition_id(
          plan, *condition, CompiledMathTimeSlot::Readiness);
}

inline semantic::Index compile_source_exact_condition_id(
    ExactVariantPlan *plan,
    const semantic::Index source_id,
    const CompiledMathTimeSlot time_slot) {
  if (source_id == semantic::kInvalidIndex) {
    return 0;
  }
  CompiledMathConditionKey key;
  append_compiled_condition_fact(
      &key,
      CompiledMathConditionFactKind::SourceExact,
      source_id,
      semantic::kInvalidIndex,
      semantic::kInvalidIndex,
      time_slot);
  return compiled_math_intern_condition(&plan->compiled_math, std::move(key));
}

inline semantic::Index compile_expr_source_node(
    ExactVariantPlan *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0) {
  return compiled_math_source_node(
      &plan->compiled_math,
      kind,
      source_id,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_expr_child_value_node(
    ExactVariantPlan *plan,
    const ExactIndexSpan children,
    const semantic::Index index,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0) {
  const auto &program = plan->lowered.program;
  return compile_expr_value_node(
      plan,
      program.expr_args[static_cast<std::size_t>(children.offset + index)],
      value_kind,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_expr_unsupported_node(
    ExactVariantPlan *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  (void)plan;
  (void)kind;
  (void)condition_id;
  throw std::runtime_error(
      "exact expression compilation reached unsupported expression " +
      std::to_string(expr_id) +
      "; runtime expression interpretation is disabled");
}

inline semantic::Index compile_guard_unless_allowed_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  const auto &program = plan->lowered.program;
  std::vector<semantic::Index> blocked_factors;
  blocked_factors.push_back(
      compile_expr_value_node(
          plan,
          kernel.guard_blocker_expr_id,
          CompiledMathNodeKind::ExprCdf,
          condition_id,
          time_id,
          source_view_id));
  bool any_unless_true = false;
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    const auto child = program.expr_args[
        static_cast<std::size_t>(kernel.children.offset + i)];
    const auto child_kind = static_cast<semantic::ExprKind>(
        program.expr_kind[static_cast<std::size_t>(child)]);
    if (child_kind == semantic::ExprKind::TrueExpr) {
      any_unless_true = true;
      continue;
    }
    if (child_kind == semantic::ExprKind::Impossible) {
      continue;
    }
    blocked_factors.push_back(
        compile_expr_value_node(
            plan,
            child,
            CompiledMathNodeKind::ExprSurvival,
            condition_id,
            time_id,
            source_view_id));
  }
  if (any_unless_true) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  const auto blocked_node =
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Product,
          std::move(blocked_factors),
          CompiledMathValueKind::Cdf);
  return compiled_math_unary_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Complement,
      blocked_node,
      CompiledMathValueKind::Cdf);
}

inline semantic::Index compile_guard_unless_density_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  const auto allowed_node =
      compile_guard_unless_allowed_node(
          plan,
          kernel,
          condition_id,
          time_id,
          source_view_id);
  const auto ref_density =
      compile_expr_value_node(
          plan,
          kernel.guard_ref_expr_id,
          CompiledMathNodeKind::ExprDensity,
          condition_id,
          time_id,
          source_view_id);
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Product,
      std::vector<semantic::Index>{ref_density, allowed_node},
      CompiledMathValueKind::Density);
}

inline bool compiled_condition_has_union_conditioning(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (const auto kind_value : condition.fact_kinds) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(kind_value);
    if (kind == CompiledMathConditionFactKind::SourceExact ||
        kind == CompiledMathConditionFactKind::SourceUpperBound ||
        kind == CompiledMathConditionFactKind::SourceLowerBound ||
        kind == CompiledMathConditionFactKind::ExprUpperBound ||
        kind == CompiledMathConditionFactKind::SourceOrder ||
        kind == CompiledMathConditionFactKind::GuardUpperBound) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_has_source_order(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]);
    if (kind == CompiledMathConditionFactKind::SourceOrder &&
        condition.fact_subject_ids[i] == before_source_id &&
        condition.fact_aux_ids[i] == after_source_id) {
      return true;
    }
  }
  return false;
}

inline semantic::Index compiled_condition_source_exact_time_id(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::SourceExact &&
        condition.fact_subject_ids[i] == source_id) {
      return condition.fact_time_ids[i];
    }
  }
  return semantic::kInvalidIndex;
}

inline bool compiled_condition_has_expr_upper_bound(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index expr_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]);
    if (kind == CompiledMathConditionFactKind::ExprUpperBound &&
        condition.fact_subject_ids[i] == expr_id) {
      return true;
    }
  }
  return false;
}

inline CompiledMathIndexSpan compiled_condition_expr_upper_bound_fact_span(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index expr_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      expr_id == semantic::kInvalidIndex) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return {};
  }
  const auto &condition = program.conditions[condition_pos];
  const auto expr_pos = static_cast<std::size_t>(expr_id);
  if (expr_pos >= condition.expr_upper_fact_spans.size()) {
    return {};
  }
  return condition.expr_upper_fact_spans[expr_pos];
}

inline CompiledMathIndexSpan compiled_condition_guard_upper_bound_fact_span(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return {};
  }
  const auto &lookup =
      program.conditions[condition_pos].guard_upper_fact_lookup;
  const auto found = std::lower_bound(
      lookup.entries.begin(),
      lookup.entries.end(),
      std::pair<semantic::Index, semantic::Index>{
          ref_source_id,
          blocker_source_id},
      [](const CompiledMathPairFactEntry &entry, const auto &target) {
        return entry.first < target.first ||
               (entry.first == target.first && entry.second < target.second);
      });
  if (found == lookup.entries.end() ||
      found->first != ref_source_id ||
      found->second != blocker_source_id) {
    return {};
  }
  return found->facts;
}

inline CompiledMathIndexSpan compile_timed_upper_bound_terms(
    CompiledMathProgram *program,
    const semantic::Index condition_id,
    const std::vector<semantic::Index> &fact_indices,
    const CompiledMathIndexSpan fact_span) {
  if (condition_id == 0 ||
      condition_id == semantic::kInvalidIndex ||
      fact_span.empty()) {
    return {};
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program->conditions.size()) {
    return {};
  }
  const auto &condition = program->conditions[condition_pos];
  const auto offset =
      static_cast<semantic::Index>(
          program->timed_upper_bound_terms.size());
  for (semantic::Index i = 0; i < fact_span.size; ++i) {
    const auto fact_pos = static_cast<std::size_t>(
        fact_indices[static_cast<std::size_t>(fact_span.offset + i)]);
    const auto time_id =
        fact_pos < condition.fact_time_ids.size()
            ? condition.fact_time_ids[fact_pos]
            : semantic::kInvalidIndex;
    const auto normalizer_node_id =
        fact_pos < condition.fact_normalizer_node_ids.size()
            ? condition.fact_normalizer_node_ids[fact_pos]
            : semantic::kInvalidIndex;
    program->timed_upper_bound_terms.push_back(
        CompiledMathTimedUpperBoundTerm{time_id, normalizer_node_id});
  }
  return CompiledMathIndexSpan{
      offset,
      static_cast<semantic::Index>(
          program->timed_upper_bound_terms.size() -
          static_cast<std::size_t>(offset))};
}

inline bool compiled_condition_has_source_relation(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id,
    const ExactRelation relation) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.source_ids.size(); ++i) {
    if (condition.source_ids[i] == source_id &&
        static_cast<ExactRelation>(condition.relations[i]) == relation) {
      return true;
    }
  }
  return false;
}

inline bool compiled_source_view_knows_before(
    const ExactVariantPlan &plan,
    const semantic::Index source_view_id,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (source_view_id == 0 ||
      source_view_id == semantic::kInvalidIndex ||
      before_source_id == semantic::kInvalidIndex ||
      after_source_id == semantic::kInvalidIndex ||
      before_source_id == after_source_id) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(source_view_id - 1U);
  if (pos >= plan.compiled_source_views.size()) {
    return false;
  }
  const auto &source_view = plan.compiled_source_views[pos];
  auto relation_for = [&](const semantic::Index source_id) {
    for (std::size_t i = 0; i < source_view.source_ids.size(); ++i) {
      if (source_view.source_ids[i] == source_id) {
        return source_view.relations[i];
      }
    }
    return ExactRelation::Unknown;
  };
  const auto before_relation = relation_for(before_source_id);
  const auto after_relation = relation_for(after_source_id);
  return before_relation != ExactRelation::Unknown &&
         after_relation != ExactRelation::Unknown &&
         before_relation < after_relation;
}

inline bool compiled_guard_order_blocks(
    const ExactVariantPlan &plan,
    const semantic::Index condition_id,
    const semantic::Index source_view_id,
    const semantic::Index blocker_source_id,
    const semantic::Index ref_source_id) {
  return compiled_condition_has_source_order(
             plan.compiled_math,
             condition_id,
             blocker_source_id,
             ref_source_id) ||
         compiled_source_view_knows_before(
             plan,
             source_view_id,
             blocker_source_id,
             ref_source_id);
}

inline bool compiled_condition_has_source_fact(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const CompiledMathConditionFactKind fact_kind,
    const semantic::Index source_id,
    const CompiledMathTimeSlot time_slot) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  const auto target_time = static_cast<semantic::Index>(time_slot);
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            fact_kind &&
        condition.fact_subject_ids[i] == source_id &&
        condition.fact_time_ids[i] == target_time) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_forces_source_after_observed(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (compiled_condition_has_source_relation(
          program, condition_id, source_id, ExactRelation::After)) {
    return true;
  }
  return compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Observed) ||
         compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Readiness) ||
         compiled_condition_has_source_fact(
             program,
             condition_id,
             CompiledMathConditionFactKind::SourceLowerBound,
             source_id,
             CompiledMathTimeSlot::Active);
}

inline bool compiled_condition_forces_source_certain(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  return compiled_condition_has_source_relation(
             program, condition_id, source_id, ExactRelation::Before) ||
         compiled_condition_has_source_relation(
             program, condition_id, source_id, ExactRelation::At);
}

inline bool expr_condition_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id);

inline bool compiled_condition_has_guard_upper_bound_expr(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index expr_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::GuardUpperBound &&
        condition.fact_subject_ids[i] == expr_id) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_has_guard_upper_bound_fact(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (const auto kind_value : condition.fact_kinds) {
    if (static_cast<CompiledMathConditionFactKind>(kind_value) ==
        CompiledMathConditionFactKind::GuardUpperBound) {
      return true;
    }
  }
  return false;
}

inline bool compiled_condition_has_guard_upper_bound_signature(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index ref_source_id,
    const semantic::Index blocker_source_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex ||
      ref_source_id == semantic::kInvalidIndex ||
      blocker_source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= program.conditions.size()) {
    return false;
  }
  const auto &condition = program.conditions[condition_pos];
  for (std::size_t i = 0; i < condition.fact_kinds.size(); ++i) {
    if (static_cast<CompiledMathConditionFactKind>(condition.fact_kinds[i]) ==
            CompiledMathConditionFactKind::GuardUpperBound &&
        condition.fact_aux_ids[i] == ref_source_id &&
        condition.fact_aux2_ids[i] == blocker_source_id) {
      return true;
    }
  }
  return false;
}

inline bool expr_condition_certain(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::TrueExpr:
    return true;
  case semantic::ExprKind::Event:
    return compiled_condition_forces_source_certain(
        plan.compiled_math,
        condition_id,
        program.expr_source_ids[pos]);
  case semantic::ExprKind::And: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_condition_certain(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      if (expr_condition_certain(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return true;
      }
    }
    return false;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::Guard:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_condition_equiv_event_source(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    const auto &kernel = plan.expr_kernels[pos];
    return kernel.simple_event_guard &&
           !kernel.has_unless &&
           kernel.guard_ref_source_id == source_id &&
           compiled_condition_forces_source_after_observed(
               plan.compiled_math,
               condition_id,
               kernel.guard_blocker_source_id);
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  bool found_candidate = false;
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    const auto child = program.expr_args[static_cast<std::size_t>(i)];
    if (expr_condition_equiv_event_source(
            plan, child, source_id, condition_id)) {
      if (found_candidate) {
        return false;
      }
      found_candidate = true;
      continue;
    }
    if (!expr_condition_certain(plan, child, condition_id)) {
      return false;
    }
  }
  return found_candidate;
}

inline bool expr_condition_certain_at_observed_cdf(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  switch (kind) {
  case semantic::ExprKind::TrueExpr:
    return true;
  case semantic::ExprKind::Event: {
    const auto source_id = program.expr_source_ids[pos];
    return compiled_condition_forces_source_certain(
               plan.compiled_math,
               condition_id,
               source_id) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Observed) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Readiness) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceExact,
               source_id,
               CompiledMathTimeSlot::Active) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Observed) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Readiness) ||
           compiled_condition_has_source_fact(
               plan.compiled_math,
               condition_id,
               CompiledMathConditionFactKind::SourceUpperBound,
               source_id,
               CompiledMathTimeSlot::Active);
    }
  case semantic::ExprKind::Guard: {
    const auto &kernel = plan.expr_kernels[pos];
    return compiled_condition_has_guard_upper_bound_expr(
               plan.compiled_math,
               condition_id,
               expr_id) ||
           (kernel.simple_event_guard &&
            !kernel.has_unless &&
            compiled_condition_has_guard_upper_bound_signature(
                plan.compiled_math,
                condition_id,
                kernel.guard_ref_source_id,
                kernel.guard_blocker_source_id));
  }
  case semantic::ExprKind::And: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    if (begin == end) {
      return false;
    }
    for (semantic::Index i = begin; i < end; ++i) {
      if (!expr_condition_certain_at_observed_cdf(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return false;
      }
    }
    return true;
  }
  case semantic::ExprKind::Or: {
    const auto begin = program.expr_arg_offsets[pos];
    const auto end = program.expr_arg_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      if (expr_condition_certain_at_observed_cdf(
              plan,
              program.expr_args[static_cast<std::size_t>(i)],
              condition_id)) {
        return true;
      }
    }
    return false;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool expr_condition_equiv_event_source_at_observed_cdf(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id,
    const semantic::Index source_id,
    const semantic::Index condition_id) {
  if (expr_id == semantic::kInvalidIndex ||
      source_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto &program = plan.lowered.program;
  const auto pos = static_cast<std::size_t>(expr_id);
  const auto kind = static_cast<semantic::ExprKind>(program.expr_kind[pos]);
  if (kind == semantic::ExprKind::Event) {
    return program.expr_source_ids[pos] == source_id;
  }
  if (kind == semantic::ExprKind::Guard) {
    return expr_condition_equiv_event_source(
        plan, expr_id, source_id, condition_id);
  }
  if (kind != semantic::ExprKind::And) {
    return false;
  }
  bool found_candidate = false;
  const auto begin = program.expr_arg_offsets[pos];
  const auto end = program.expr_arg_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    const auto child = program.expr_args[static_cast<std::size_t>(i)];
    if (expr_condition_equiv_event_source_at_observed_cdf(
            plan, child, source_id, condition_id)) {
      if (found_candidate) {
        return false;
      }
      found_candidate = true;
      continue;
    }
    if (!expr_condition_certain_at_observed_cdf(
            plan, child, condition_id)) {
      return false;
    }
  }
  return found_candidate;
}

inline bool compiled_condition_absorbs_union_to_source(
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const bool observed_cdf_context,
    semantic::Index *source_id) {
  for (semantic::Index i = 0;
       i < kernel.union_absorption_candidate_span.size;
       ++i) {
    const auto candidate = plan.expr_union_absorption_sources[
        static_cast<std::size_t>(
            kernel.union_absorption_candidate_span.offset + i)];
    for (semantic::Index j = 0; j < kernel.children.size; ++j) {
      const auto child = plan.lowered.program.expr_args[
          static_cast<std::size_t>(kernel.children.offset + j)];
      const bool child_equiv =
          observed_cdf_context
              ? expr_condition_equiv_event_source_at_observed_cdf(
                    plan, child, candidate, condition_id)
              : expr_condition_equiv_event_source(
                    plan, child, candidate, condition_id);
      if (child_equiv) {
        if (source_id != nullptr) {
          *source_id = candidate;
        }
        return true;
      }
    }
  }
  return false;
}

inline semantic::Index compile_integral_zero_to_current_node(
    ExactVariantPlan *plan,
    const semantic::Index integrand_node,
    const semantic::Index condition_id,
    const semantic::Index time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto integrand_root =
      compiled_math_make_root(&plan->compiled_math, integrand_node);
  return compiled_math_integral_zero_to_current_node(
      &plan->compiled_math,
      integrand_root,
      condition_id,
      time_id,
      source_view_id,
      bind_time_id);
}

inline semantic::Index compile_simple_guard_raw_density_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  std::vector<semantic::Index> factors;
  factors.reserve(2U);
  factors.push_back(
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourcePdf,
          kernel.guard_ref_source_id,
          condition_id,
          time_id,
          source_view_id));
  factors.push_back(
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourceSurvival,
          kernel.guard_blocker_source_id,
          condition_id,
          time_id,
          source_view_id));
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Product,
      std::move(factors),
      CompiledMathValueKind::Density);
}

inline semantic::Index compile_simple_guard_upper_wrapper_node(
    ExactVariantPlan *plan,
    const CompiledMathNodeKind kind,
    const semantic::Index expr_id,
    const semantic::Index raw_node,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = kind == CompiledMathNodeKind::SimpleGuardDensity
                       ? CompiledMathValueKind::Density
                       : CompiledMathValueKind::Cdf;
  key.subject_id = expr_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.source_view_id = source_view_id;
  if (expr_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(expr_id) >= plan->expr_kernels.size()) {
    throw std::runtime_error(
        "compiled simple guard upper wrapper has no expression kernel");
  }
  const auto &kernel = plan->expr_kernels[static_cast<std::size_t>(expr_id)];
  const auto upper_span =
      compiled_condition_guard_upper_bound_fact_span(
          plan->compiled_math,
          condition_id,
          kernel.guard_ref_source_id,
          kernel.guard_blocker_source_id);
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (upper_span.empty() ||
      condition_id == 0 ||
      condition_id == semantic::kInvalidIndex ||
      condition_pos >= plan->compiled_math.conditions.size()) {
    throw std::runtime_error(
        "compiled simple guard upper wrapper has no planned fact span");
  }
  const auto term_span =
      compile_timed_upper_bound_terms(
          &plan->compiled_math,
          condition_id,
          plan->compiled_math.conditions[condition_pos]
              .guard_upper_fact_lookup.fact_indices,
          upper_span);
  key.aux_id =
      term_span.empty() ? semantic::kInvalidIndex : term_span.offset;
  key.aux2_id = term_span.size;
  key.children.push_back(raw_node);
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline bool compile_simple_guard_needs_upper_wrapper(
    const ExactVariantPlan &plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id) {
  return compiled_condition_has_guard_upper_bound_signature(
             plan.compiled_math,
             condition_id,
             kernel.guard_ref_source_id,
             kernel.guard_blocker_source_id);
}

inline semantic::Index compile_simple_guard_density_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (compiled_guard_order_blocks(
          *plan,
          condition_id,
          source_view_id,
          kernel.guard_blocker_source_id,
          kernel.guard_ref_source_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  const auto raw_node =
      compile_simple_guard_raw_density_node(
          plan, kernel, condition_id, time_id, source_view_id);
  if (!compile_simple_guard_needs_upper_wrapper(
          *plan, kernel, condition_id)) {
    return raw_node;
  }
  return compile_simple_guard_upper_wrapper_node(
      plan,
      CompiledMathNodeKind::SimpleGuardDensity,
      expr_id,
      raw_node,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_simple_guard_raw_cdf_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (compiled_guard_order_blocks(
          *plan,
          condition_id,
          source_view_id,
          kernel.guard_blocker_source_id,
          kernel.guard_ref_source_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }

  const auto ref_exact_time_id =
      compiled_condition_source_exact_time_id(
          plan->compiled_math,
          condition_id,
          kernel.guard_ref_source_id);
  if (ref_exact_time_id != semantic::kInvalidIndex) {
    const auto blocker_survival_at_ref =
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceSurvival,
            kernel.guard_blocker_source_id,
            condition_id,
            ref_exact_time_id,
            source_view_id);
    return compiled_math_time_gate_node(
        &plan->compiled_math,
        blocker_survival_at_ref,
        time_id,
        ref_exact_time_id,
        CompiledMathValueKind::Cdf);
  }

  const auto blocker_exact_time_id =
      compiled_condition_source_exact_time_id(
          plan->compiled_math,
          condition_id,
          kernel.guard_blocker_source_id);
  if (blocker_exact_time_id != semantic::kInvalidIndex) {
    return compiled_math_source_node(
        &plan->compiled_math,
        CompiledMathNodeKind::SourceCdf,
        kernel.guard_ref_source_id,
        condition_id,
        time_id,
        source_view_id,
        blocker_exact_time_id);
  }

  (void)expr_id;
  const auto bind_time_id =
      time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
          ? time_id
          : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
  const auto raw_density_node =
      compile_simple_guard_raw_density_node(
          plan, kernel, condition_id, bind_time_id, source_view_id);
  return compile_integral_zero_to_current_node(
      plan,
      raw_density_node,
      condition_id,
      time_id,
      source_view_id,
      bind_time_id);
}

inline semantic::Index compile_simple_guard_cdf_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  const auto raw_node =
      compile_simple_guard_raw_cdf_node(
          plan, expr_id, kernel, condition_id, time_id, source_view_id);
  if (!compile_simple_guard_needs_upper_wrapper(
          *plan, kernel, condition_id)) {
    return raw_node;
  }
  return compile_simple_guard_upper_wrapper_node(
      plan,
      CompiledMathNodeKind::SimpleGuardCdf,
      expr_id,
      raw_node,
      condition_id,
      time_id,
      source_view_id);
}

inline semantic::Index compile_signed_term_node(
    ExactVariantPlan *plan,
    const semantic::Index node_id,
    const int sign) {
  if (sign >= 0) {
    return node_id;
  }
  return compiled_math_unary_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Negate,
      node_id);
}

inline semantic::Index compile_outcome_subset_unused_node(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &outcome_indices) {
  if (outcome_indices.empty()) {
    return compiled_math_constant(&plan->compiled_math, 1.0);
  }
  const auto offset =
      static_cast<semantic::Index>(
          plan->compiled_outcome_gate_indices.size());
  plan->compiled_outcome_gate_indices.insert(
      plan->compiled_outcome_gate_indices.end(),
      outcome_indices.begin(),
      outcome_indices.end());
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::OutcomeSubsetUnused;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.subject_id = offset;
  key.aux_id = static_cast<semantic::Index>(outcome_indices.size());
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline semantic::Index compile_expr_span_conjunction_density_node(
    ExactVariantPlan *plan,
    const std::vector<semantic::Index> &arena,
    const ExactIndexSpan span,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (span.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if (span.size == 1U) {
    return compile_expr_value_node(
        plan,
        arena[static_cast<std::size_t>(span.offset)],
        CompiledMathNodeKind::ExprDensity,
        condition_id,
        time_id,
        source_view_id);
  }

  std::vector<semantic::Index> terms;
  terms.reserve(static_cast<std::size_t>(span.size));
  for (semantic::Index i = 0; i < span.size; ++i) {
    std::vector<semantic::Index> factors;
    factors.reserve(static_cast<std::size_t>(span.size));
    factors.push_back(
        compile_expr_value_node(
            plan,
            arena[static_cast<std::size_t>(span.offset + i)],
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
    for (semantic::Index j = 0; j < span.size; ++j) {
      if (i == j) {
        continue;
      }
      factors.push_back(
          compile_expr_value_node(
              plan,
              arena[static_cast<std::size_t>(span.offset + j)],
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    terms.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors)));
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline bool exact_union_has_multi_subset(const ExactVariantPlan &plan,
                                         const ExactExprKernel &kernel) {
  for (semantic::Index i = 0; i < kernel.union_subset_span.size; ++i) {
    const auto &subset = plan.expr_union_subsets[
        static_cast<std::size_t>(kernel.union_subset_span.offset + i)];
    if (subset.child_positions.size > 1U) {
      return true;
    }
  }
  return false;
}

inline semantic::Index compile_union_subset_density_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const ExactExprUnionSubset &subset,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (subset.child_positions.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if (subset.child_positions.size == 1U) {
    const auto active_pos =
        plan->expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset)];
    return compile_expr_child_value_node(
        plan,
        kernel.children,
        active_pos,
        CompiledMathNodeKind::ExprDensity,
        condition_id,
        time_id,
        source_view_id);
  }

  std::vector<semantic::Index> terms;
  terms.reserve(static_cast<std::size_t>(subset.child_positions.size));
  for (semantic::Index i = 0; i < subset.child_positions.size; ++i) {
    const auto active_pos =
        plan->expr_union_subset_child_positions[
            static_cast<std::size_t>(subset.child_positions.offset + i)];
    std::vector<semantic::Index> factors;
    factors.reserve(static_cast<std::size_t>(subset.child_positions.size));
    factors.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            active_pos,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
    for (semantic::Index j = 0; j < subset.child_positions.size; ++j) {
      if (i == j) {
        continue;
      }
      const auto other_pos =
          plan->expr_union_subset_child_positions[
              static_cast<std::size_t>(subset.child_positions.offset + j)];
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              other_pos,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    terms.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors),
            CompiledMathValueKind::Density));
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline semantic::Index compile_union_density_sum_node(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id,
    const bool multi_subset_only) {
  std::vector<semantic::Index> terms;
  for (semantic::Index i = 0; i < kernel.union_subset_span.size; ++i) {
    const auto &subset =
        plan->expr_union_subsets[
            static_cast<std::size_t>(kernel.union_subset_span.offset + i)];
    if (multi_subset_only && subset.child_positions.size <= 1U) {
      continue;
    }
    const auto subset_node =
        compile_union_subset_density_node(
            plan,
            kernel,
            subset,
            condition_id,
            time_id,
            source_view_id);
    terms.push_back(
        compile_signed_term_node(plan, subset_node, subset.sign));
  }
  if (terms.empty()) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::CleanSignedSum,
      std::move(terms),
      CompiledMathValueKind::Density);
}

inline std::vector<semantic::Index> compile_union_child_density_cdf_nodes(
    ExactVariantPlan *plan,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(kernel.children.size) * 2U);
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id));
  }
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id));
  }
  return children;
}

inline semantic::Index compile_union_multi_subset_density_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto density_time_id =
      bind_time_id == semantic::kInvalidIndex ? time_id : bind_time_id;
  (void)expr_id;
  return compile_union_density_sum_node(
      plan,
      kernel,
      condition_id,
      density_time_id,
      source_view_id,
      true);
}

inline semantic::Index compile_overlapping_union_kernel_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (value_kind == CompiledMathNodeKind::ExprSurvival) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        compile_overlapping_union_kernel_node(
            plan,
            expr_id,
            kernel,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Survival);
  }

  auto *program = &plan->compiled_math;
  if (value_kind == CompiledMathNodeKind::ExprDensity) {
    (void)expr_id;
    return compile_union_density_sum_node(
        plan,
        kernel,
        condition_id,
        time_id,
        source_view_id,
        false);
  }

  std::vector<semantic::Index> children;
  children.reserve(static_cast<std::size_t>(kernel.children.size) + 1U);
  for (semantic::Index i = 0; i < kernel.children.size; ++i) {
    children.push_back(
        compile_expr_child_value_node(
            plan,
            kernel.children,
            i,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id));
  }
  if (exact_union_has_multi_subset(*plan, kernel)) {
    const auto bind_time_id =
        time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
            ? time_id
            : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
    const auto multi_density_node =
        compile_union_multi_subset_density_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id,
            bind_time_id);
    const auto multi_density_root =
        compiled_math_make_root(program, multi_density_node);
    children.push_back(
        compiled_math_union_kernel_node(
            program,
            CompiledMathNodeKind::UnionKernelMultiSubsetCdf,
            expr_id,
            condition_id,
            {},
            multi_density_root,
            time_id,
            source_view_id,
            bind_time_id));
  }
  return compiled_math_union_kernel_node(
      program,
      CompiledMathNodeKind::UnionKernelCdf,
      expr_id,
      condition_id,
      std::move(children),
      semantic::kInvalidIndex,
      time_id,
      source_view_id);
}

inline semantic::Index compile_overlapping_union_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const ExactExprKernel &kernel,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (value_kind == CompiledMathNodeKind::ExprSurvival) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        compile_overlapping_union_node(
            plan,
            expr_id,
            kernel,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Survival);
  }

  if (!kernel.union_absorption_candidate_span.empty() &&
      compiled_condition_has_union_conditioning(plan->compiled_math, condition_id)) {
    semantic::Index absorbed_source_id{semantic::kInvalidIndex};
    const bool absorbed = compiled_condition_absorbs_union_to_source(
        *plan,
        kernel,
        condition_id,
        value_kind == CompiledMathNodeKind::ExprCdf,
        &absorbed_source_id);
    if (absorbed) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_expr_source_node(
            plan,
            CompiledMathNodeKind::SourcePdf,
            absorbed_source_id,
            condition_id,
            time_id,
            source_view_id);
      }
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceCdf,
          absorbed_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
  }

  return compile_overlapping_union_kernel_node(
      plan, expr_id, kernel, value_kind, condition_id, time_id, source_view_id);
}

inline semantic::Index compile_expr_value_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id);

inline semantic::Index compile_expr_value_node_raw(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  const auto &program = plan->lowered.program;
  const auto &kernel = plan->expr_kernels[static_cast<std::size_t>(expr_id)];
  const auto constant = [&](const double value) {
    return compiled_math_constant(&plan->compiled_math, value);
  };
  const auto complement = [&](const semantic::Index child) {
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Complement,
        child,
        CompiledMathValueKind::Cdf);
  };
  const auto unsupported = [&]() {
    return compile_expr_unsupported_node(plan, value_kind, expr_id, condition_id);
  };

  switch (kernel.kind) {
  case semantic::ExprKind::Impossible:
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return constant(1.0);
    }
    return constant(0.0);

  case semantic::ExprKind::TrueExpr:
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return constant(0.0);
    }
    return constant(1.0);

  case semantic::ExprKind::Event:
    if (value_kind == CompiledMathNodeKind::ExprDensity) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourcePdf,
          kernel.event_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return compile_expr_source_node(
          plan,
          CompiledMathNodeKind::SourceCdf,
          kernel.event_source_id,
          condition_id,
          time_id,
          source_view_id);
    }
    return compile_expr_source_node(
        plan,
        CompiledMathNodeKind::SourceSurvival,
        kernel.event_source_id,
        condition_id,
        time_id,
        source_view_id);

  case semantic::ExprKind::And: {
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      std::vector<semantic::Index> children;
      children.reserve(static_cast<std::size_t>(kernel.children.size));
      for (semantic::Index i = 0; i < kernel.children.size; ++i) {
        children.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                i,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id));
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(children)),
          CompiledMathValueKind::Cdf);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return complement(
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    std::vector<semantic::Index> terms;
    terms.reserve(static_cast<std::size_t>(kernel.children.size));
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      std::vector<semantic::Index> factors;
      factors.reserve(static_cast<std::size_t>(kernel.children.size));
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              i,
              CompiledMathNodeKind::ExprDensity,
              condition_id,
              time_id,
              source_view_id));
      for (semantic::Index j = 0; j < kernel.children.size; ++j) {
        if (i == j) {
          continue;
        }
        factors.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                j,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id));
      }
      terms.push_back(
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(factors)));
    }
    return compiled_math_algebra_node(
        &plan->compiled_math,
        CompiledMathNodeKind::CleanSignedSum,
        std::move(terms),
        CompiledMathValueKind::Density);
  }

	  case semantic::ExprKind::Or: {
	    if (kernel.children_overlap) {
	      return compile_overlapping_union_node(
	          plan,
	          expr_id,
          kernel,
          value_kind,
          condition_id,
          time_id,
          source_view_id);
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      std::vector<semantic::Index> children;
      children.reserve(static_cast<std::size_t>(kernel.children.size));
      for (semantic::Index i = 0; i < kernel.children.size; ++i) {
        children.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                i,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(children)),
          CompiledMathValueKind::Survival);
    }
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return complement(
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprSurvival,
              condition_id,
              time_id,
              source_view_id));
    }
    std::vector<semantic::Index> terms;
    terms.reserve(static_cast<std::size_t>(kernel.children.size));
    for (semantic::Index i = 0; i < kernel.children.size; ++i) {
      std::vector<semantic::Index> factors;
      factors.reserve(static_cast<std::size_t>(kernel.children.size));
      factors.push_back(
          compile_expr_child_value_node(
              plan,
              kernel.children,
              i,
              CompiledMathNodeKind::ExprDensity,
              condition_id,
              time_id,
              source_view_id));
      for (semantic::Index j = 0; j < kernel.children.size; ++j) {
        if (i == j) {
          continue;
        }
        factors.push_back(
            compile_expr_child_value_node(
                plan,
                kernel.children,
                j,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
      }
      terms.push_back(
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(factors)));
    }
    return compiled_math_algebra_node(
        &plan->compiled_math,
        CompiledMathNodeKind::CleanSignedSum,
        std::move(terms),
        CompiledMathValueKind::Density);
  }

  case semantic::ExprKind::Not: {
    const auto child =
        program.expr_args[static_cast<std::size_t>(kernel.children.offset)];
    if (value_kind == CompiledMathNodeKind::ExprCdf) {
      return complement(
          compile_expr_value_node(
              plan,
              child,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id));
    }
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return compile_expr_value_node(
          plan,
          child,
          CompiledMathNodeKind::ExprCdf,
          condition_id,
          time_id,
          source_view_id);
    }
    return compiled_math_unary_node(
        &plan->compiled_math,
        CompiledMathNodeKind::Negate,
        compile_expr_value_node(
            plan,
            child,
            CompiledMathNodeKind::ExprDensity,
            condition_id,
            time_id,
            source_view_id),
        CompiledMathValueKind::Density);
  }

  case semantic::ExprKind::Guard:
    if (kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_guard_unless_density_node(
            plan,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      } else if (value_kind == CompiledMathNodeKind::ExprCdf) {
        const auto bind_time_id =
            time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
                ? time_id
                : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
        const auto density_node =
            compile_guard_unless_density_node(
                plan,
                kernel,
                condition_id,
                bind_time_id,
                source_view_id);
        return compile_integral_zero_to_current_node(
            plan,
            density_node,
            condition_id,
            time_id,
            source_view_id,
            bind_time_id);
      } else {
        return compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Complement,
            compile_expr_value_node(
                plan,
                expr_id,
                CompiledMathNodeKind::ExprCdf,
                condition_id,
                time_id,
                source_view_id),
            CompiledMathValueKind::Survival);
      }
    }
    if (kernel.simple_event_guard && !kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        return compile_simple_guard_density_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      }
      if (value_kind == CompiledMathNodeKind::ExprCdf) {
        return compile_simple_guard_cdf_node(
            plan,
            expr_id,
            kernel,
            condition_id,
            time_id,
            source_view_id);
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Complement,
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id),
          CompiledMathValueKind::Survival);
    }
    if (!kernel.has_unless) {
      if (value_kind == CompiledMathNodeKind::ExprDensity) {
        std::vector<semantic::Index> factors;
        factors.reserve(2U);
        factors.push_back(
            compile_expr_value_node(
                plan,
                kernel.guard_ref_expr_id,
                CompiledMathNodeKind::ExprDensity,
                condition_id,
                time_id,
                source_view_id));
        factors.push_back(
            compile_expr_value_node(
                plan,
                kernel.guard_blocker_expr_id,
                CompiledMathNodeKind::ExprSurvival,
                condition_id,
                time_id,
                source_view_id));
        return compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::move(factors),
            CompiledMathValueKind::Density);
      }
      if (value_kind == CompiledMathNodeKind::ExprCdf) {
        const auto bind_time_id =
            time_id == static_cast<semantic::Index>(CompiledMathTimeSlot::Active)
                ? time_id
                : static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
        const auto density_node =
            compile_expr_value_node(
                plan,
                expr_id,
                CompiledMathNodeKind::ExprDensity,
                condition_id,
                bind_time_id,
                source_view_id);
        return compile_integral_zero_to_current_node(
            plan,
            density_node,
            condition_id,
            time_id,
            source_view_id,
            bind_time_id);
      }
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Complement,
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id),
          CompiledMathValueKind::Survival);
    }
    return unsupported();
  }

  return unsupported();
}

inline semantic::Index compile_expr_upper_bound_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const semantic::Index child_node,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  CompiledMathNodeKey key;
  key.kind = value_kind == CompiledMathNodeKind::ExprDensity
                 ? CompiledMathNodeKind::ExprUpperBoundDensity
                 : CompiledMathNodeKind::ExprUpperBoundCdf;
  key.value_kind = value_kind == CompiledMathNodeKind::ExprDensity
                       ? CompiledMathValueKind::Density
                       : CompiledMathValueKind::Cdf;
  key.subject_id = expr_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.source_view_id = source_view_id;
  const auto upper_span =
      compiled_condition_expr_upper_bound_fact_span(
          plan->compiled_math, condition_id, expr_id);
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  const auto term_span =
      condition_id == 0 ||
              condition_id == semantic::kInvalidIndex ||
              condition_pos >= plan->compiled_math.conditions.size()
          ? CompiledMathIndexSpan{}
          : compile_timed_upper_bound_terms(
                &plan->compiled_math,
                condition_id,
                plan->compiled_math.conditions[condition_pos]
                    .expr_upper_fact_indices,
                upper_span);
  key.aux_id =
      term_span.empty() ? semantic::kInvalidIndex : term_span.offset;
  key.aux2_id = term_span.size;
  key.children.push_back(child_node);
  return compiled_math_intern_node(&plan->compiled_math, std::move(key));
}

inline bool sequence_expr_upper_bound_used(
    const ExactVariantPlan &plan,
    const semantic::Index expr_id) {
  return expr_id != semantic::kInvalidIndex &&
         static_cast<std::size_t>(expr_id) <
             plan.sequence_expr_upper_bound_used.size() &&
         plan.sequence_expr_upper_bound_used[
             static_cast<std::size_t>(expr_id)] != 0U;
}

inline semantic::Index compile_expr_value_node(
    ExactVariantPlan *plan,
    const semantic::Index expr_id,
    const CompiledMathNodeKind value_kind,
    const semantic::Index condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id) {
  if (compiled_condition_has_expr_upper_bound(
          plan->compiled_math,
          condition_id,
          expr_id) ||
      sequence_expr_upper_bound_used(*plan, expr_id)) {
    if (value_kind == CompiledMathNodeKind::ExprSurvival) {
      return compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Complement,
          compile_expr_value_node(
              plan,
              expr_id,
              CompiledMathNodeKind::ExprCdf,
              condition_id,
              time_id,
              source_view_id),
          CompiledMathValueKind::Survival);
    }
    const auto raw_node =
        compile_expr_value_node_raw(
            plan, expr_id, value_kind, condition_id, time_id, source_view_id);
    if (value_kind == CompiledMathNodeKind::ExprDensity ||
        value_kind == CompiledMathNodeKind::ExprCdf) {
      return compile_expr_upper_bound_node(
          plan,
          expr_id,
          raw_node,
          value_kind,
          condition_id,
          time_id,
          source_view_id);
    }
    return raw_node;
  }
  return compile_expr_value_node_raw(
      plan, expr_id, value_kind, condition_id, time_id, source_view_id);
}

inline semantic::Index compile_runtime_factors_node(
    ExactVariantPlan *plan,
    const ExactRuntimeFactors &factors,
    const semantic::Index source_condition_id,
    const semantic::Index expr_condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  std::vector<semantic::Index> children;
  children.reserve(
      factors.source_pdf.size() +
      factors.source_cdf.size() +
      factors.source_survival.size() +
      factors.expr_density.size() +
      factors.expr_cdf.size() +
      factors.expr_survival.size());

  for (const auto source_id : factors.source_pdf) {
    children.push_back(
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourcePdf,
            source_id,
            source_condition_id,
            time_id));
  }
  for (const auto source_id : factors.source_cdf) {
    children.push_back(
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceCdf,
            source_id,
            source_condition_id,
            time_id));
  }
  for (const auto source_id : factors.source_survival) {
    children.push_back(
        compiled_math_source_node(
            &plan->compiled_math,
            CompiledMathNodeKind::SourceSurvival,
            source_id,
            source_condition_id,
            time_id));
  }
  for (const auto expr_id : factors.expr_density) {
    children.push_back(
        compile_expr_value_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprDensity,
            expr_condition_id,
            time_id,
            source_view_id));
  }
  for (const auto expr_id : factors.expr_cdf) {
    children.push_back(
        compile_expr_value_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprCdf,
            expr_condition_id,
            time_id,
            source_view_id));
  }
  for (const auto expr_id : factors.expr_survival) {
    children.push_back(
        compile_expr_value_node(
            plan,
            expr_id,
            CompiledMathNodeKind::ExprSurvival,
            expr_condition_id,
            time_id,
            source_view_id));
  }

  return compiled_math_algebra_node(
      &plan->compiled_math,
      CompiledMathNodeKind::Product,
      std::move(children));
}

inline void compile_runtime_truth_formula_math(
    ExactVariantPlan *plan,
    ExactRuntimeTruthFormula *formula,
    const semantic::Index source_condition_id,
    const semantic::Index expr_condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (!formula->sum_of_products) {
    const auto product_node =
        compile_runtime_factors_empty(formula->product)
            ? compiled_math_constant(&plan->compiled_math, formula->empty_value)
            : compile_runtime_factors_node(
                  plan,
                  formula->product,
                  source_condition_id,
                  expr_condition_id,
                  time_id,
                  source_view_id);
    const auto root_node =
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::ClampProbability,
            product_node,
            CompiledMathValueKind::Cdf);
    formula->compiled_root_id =
        compiled_math_make_root(&plan->compiled_math, root_node);
    return;
  }

  std::vector<semantic::Index> term_nodes;
  term_nodes.reserve(formula->sum_terms.size());
  for (auto &term : formula->sum_terms) {
    const auto term_node =
        compile_runtime_factors_node(
            plan,
            term.factors,
            source_condition_id,
            expr_condition_id,
            time_id,
            source_view_id);
    term.compiled_root_id =
        compiled_math_make_root(&plan->compiled_math, term_node);
    term_nodes.push_back(term_node);
  }
  const auto sum_node =
      term_nodes.empty()
          ? compiled_math_constant(&plan->compiled_math, formula->empty_value)
          : compiled_math_algebra_node(
                &plan->compiled_math,
                formula->clean_signed
                    ? CompiledMathNodeKind::CleanSignedSum
                    : CompiledMathNodeKind::Sum,
                std::move(term_nodes));
  formula->compiled_root_id =
      compiled_math_make_root(&plan->compiled_math, sum_node);
}

inline semantic::Index compile_runtime_truth_formula_root(
    ExactVariantPlan *plan,
    const ExactRuntimeTruthFormula &formula,
    const semantic::Index source_condition_id,
    const semantic::Index expr_condition_id,
    const semantic::Index time_id,
    const semantic::Index source_view_id = 0) {
  if (!formula.sum_of_products) {
    const auto product_node =
        compile_runtime_factors_empty(formula.product)
            ? compiled_math_constant(&plan->compiled_math, formula.empty_value)
            : compile_runtime_factors_node(
                  plan,
                  formula.product,
                  source_condition_id,
                  expr_condition_id,
                  time_id,
                  source_view_id);
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::ClampProbability,
            product_node,
            CompiledMathValueKind::Cdf));
  }

  std::vector<semantic::Index> term_nodes;
  term_nodes.reserve(formula.sum_terms.size());
  for (const auto &term : formula.sum_terms) {
    term_nodes.push_back(
        compile_runtime_factors_node(
            plan,
            term.factors,
            source_condition_id,
            expr_condition_id,
            time_id,
            source_view_id));
  }
  const auto sum_node =
      term_nodes.empty()
          ? compiled_math_constant(&plan->compiled_math, formula.empty_value)
          : compiled_math_algebra_node(
                &plan->compiled_math,
                formula.clean_signed
                    ? CompiledMathNodeKind::CleanSignedSum
                    : CompiledMathNodeKind::Sum,
                std::move(term_nodes));
  return compiled_math_make_root(&plan->compiled_math, sum_node);
}

inline semantic::Index compiled_math_root_node_id(
    const CompiledMathProgram &program,
    const semantic::Index root_id) {
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  return program.roots[static_cast<std::size_t>(root_id)].node_id;
}

inline semantic::Index compile_runtime_scenario_same_active_root(
    ExactVariantPlan *plan,
    const ExactRuntimeScenarioFormula &formula,
    const semantic::Index extra_condition_id,
    const semantic::Index ready_upper_time_id =
        static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
    const ExactRuntimeTruthFormula *covered_after_survival = nullptr);

inline void erase_runtime_factor_once(std::vector<semantic::Index> *factors,
                                      const semantic::Index value) {
  const auto it = std::find(factors->begin(), factors->end(), value);
  if (it != factors->end()) {
    factors->erase(it);
  }
}

inline ExactRuntimeTruthFormula runtime_after_survival_excluding_covered(
    ExactRuntimeTruthFormula after_survival,
    const ExactRuntimeTruthFormula *covered_after_survival) {
  if (covered_after_survival == nullptr) {
    return after_survival;
  }
  for (const auto source_id : covered_after_survival->product.source_survival) {
    erase_runtime_factor_once(&after_survival.product.source_survival, source_id);
  }
  for (const auto expr_id : covered_after_survival->product.expr_survival) {
    erase_runtime_factor_once(&after_survival.product.expr_survival, expr_id);
  }
  return after_survival;
}

inline semantic::Index compile_runtime_scenario_truth_cdf_root(
    ExactVariantPlan *plan,
    const ExactRuntimeScenarioFormula &formula,
    const semantic::Index extra_condition_id) {
  const auto active_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Active);
  const auto upper_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto active_condition_id =
      compile_source_exact_condition_id(
          plan,
          formula.active_source_id,
          CompiledMathTimeSlot::Active);
  const auto relation_condition_id =
      compile_relation_condition_id(plan, formula.relation_template);
  const auto source_condition_id =
      compiled_math_merge_conditions(
          &plan->compiled_math,
          {active_condition_id, extra_condition_id});
  const auto expr_condition_id =
      compiled_math_merge_conditions(
          &plan->compiled_math,
          {relation_condition_id, source_condition_id});

  std::vector<semantic::Index> factors;
  factors.push_back(
      compiled_math_source_node(
          &plan->compiled_math,
          CompiledMathNodeKind::SourcePdf,
          formula.active_source_id,
          0,
          active_time_id));

  const auto after_node =
      compiled_math_root_node_id(
          plan->compiled_math,
          compile_runtime_truth_formula_root(
              plan,
              formula.after_survival,
              source_condition_id,
              expr_condition_id,
              active_time_id,
              formula.source_view_id));
  if (after_node != semantic::kInvalidIndex) {
    factors.push_back(after_node);
  }
  if (formula.has_readiness) {
    const auto readiness_node =
        compile_runtime_factors_empty(formula.readiness_cdf.product)
            ? compiled_math_constant(
                  &plan->compiled_math,
                  formula.readiness_cdf.empty_value)
            : compile_runtime_factors_node(
                  plan,
                  formula.readiness_cdf.product,
                  source_condition_id,
                  expr_condition_id,
                  active_time_id,
                  formula.source_view_id);
    factors.push_back(
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::ClampProbability,
            readiness_node,
            CompiledMathValueKind::Cdf));
  }

  const auto integrand_root =
      compiled_math_make_root(
          &plan->compiled_math,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(factors)));
  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_integral_zero_to_current_node(
          &plan->compiled_math,
          integrand_root,
          source_condition_id,
          upper_time_id,
          formula.source_view_id,
          active_time_id));
}

inline semantic::Index compile_runtime_scenario_same_active_root(
    ExactVariantPlan *plan,
    const ExactRuntimeScenarioFormula &formula,
    const semantic::Index extra_condition_id,
    const semantic::Index ready_upper_time_id,
    const ExactRuntimeTruthFormula *covered_after_survival) {
  const auto relation_condition_id =
      compile_relation_condition_id(plan, formula.relation_template);
  const auto base_expr_condition_id =
      compiled_math_merge_conditions(
          &plan->compiled_math,
          {relation_condition_id, extra_condition_id});
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);

  const auto after_survival =
      runtime_after_survival_excluding_covered(
          formula.after_survival,
          covered_after_survival);

  auto compile_tail_node = [&](const semantic::Index condition_id) {
    const auto expr_condition_id =
        compiled_math_merge_conditions(
            &plan->compiled_math,
            {relation_condition_id, condition_id});
    return compiled_math_root_node_id(
        plan->compiled_math,
        compile_runtime_truth_formula_root(
            plan,
            after_survival,
            condition_id,
            expr_condition_id,
            observed_time_id,
            formula.source_view_id));
  };

  if (!formula.has_readiness) {
    return compile_runtime_truth_formula_root(
        plan,
        after_survival,
        extra_condition_id,
        base_expr_condition_id,
        observed_time_id,
        formula.source_view_id);
  }

  std::vector<semantic::Index> branches;
  const auto base_tail_node = compile_tail_node(extra_condition_id);
  const auto initial_ready_node =
      compiled_math_root_node_id(
          plan->compiled_math,
          compile_runtime_truth_formula_root(
              plan,
              formula.readiness_cdf,
              extra_condition_id,
              base_expr_condition_id,
              zero_time_id,
              formula.source_view_id));
  if (initial_ready_node != semantic::kInvalidIndex &&
      base_tail_node != semantic::kInvalidIndex) {
    branches.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{initial_ready_node, base_tail_node},
            CompiledMathValueKind::Scalar));
  }

  std::vector<semantic::Index> term_nodes;
  term_nodes.reserve(formula.readiness_density.sum_terms.size());
  for (const auto &term : formula.readiness_density.sum_terms) {
    if (runtime_condition_order_contradiction(term.condition)) {
      continue;
    }
    const auto tail_condition =
        filter_runtime_condition_for_truth(
            term.condition,
            after_survival,
            *plan);
    ExactRuntimeTermCondition compiled_tail_condition = tail_condition;
    compile_runtime_condition_ids(plan, &compiled_tail_condition);
    const auto tail_context_id =
        compiled_math_merge_conditions(
            &plan->compiled_math,
            {extra_condition_id,
             compiled_tail_condition.readiness_condition_id});
    if (compiled_condition_impossible(plan->compiled_math, tail_context_id)) {
      continue;
    }
    const auto term_node =
        compile_runtime_factors_node(
            plan,
            term.factors,
            extra_condition_id,
            base_expr_condition_id,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
            formula.source_view_id);
    const auto tail_node = compile_tail_node(tail_context_id);
    if (term_node == semantic::kInvalidIndex ||
        tail_node == semantic::kInvalidIndex) {
      continue;
    }
    term_nodes.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{term_node, tail_node},
            CompiledMathValueKind::Scalar));
  }

  const auto integrand_node =
      term_nodes.empty()
          ? compiled_math_constant(&plan->compiled_math, 0.0)
          : compiled_math_algebra_node(
                &plan->compiled_math,
                CompiledMathNodeKind::CleanSignedSum,
                std::move(term_nodes),
                CompiledMathValueKind::Scalar);
  const auto integrand_root =
      compiled_math_make_root(&plan->compiled_math, integrand_node);
  branches.push_back(
      compiled_math_raw_integral_zero_to_current_node(
          &plan->compiled_math,
          integrand_root,
          0,
          ready_upper_time_id,
          0,
          static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness)));

  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::CleanSignedSum,
              std::move(branches),
              CompiledMathValueKind::Scalar),
          CompiledMathValueKind::Cdf));
}

inline semantic::Index compile_competitor_subset_win_cdf_root(
    ExactVariantPlan *plan,
    ExactRuntimeCompetitorSubsetPlan *subset,
    const semantic::Index condition_id) {
  if (condition_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  if (condition_id == 0 &&
      subset->win_cdf_root_id != semantic::kInvalidIndex) {
    return subset->win_cdf_root_id;
  }
  for (const auto &root : subset->conditioned_win_cdf_roots) {
    if (root.condition_id == condition_id) {
      return root.root_id;
    }
  }
  semantic::Index root_id = semantic::kInvalidIndex;
  if (compiled_condition_impossible(plan->compiled_math, condition_id)) {
    root_id =
        compiled_math_make_root(
            &plan->compiled_math,
            compiled_math_constant(&plan->compiled_math, 0.0));
  } else if (subset->singleton_expr_root != semantic::kInvalidIndex) {
    const auto node =
        compile_expr_value_node(
            plan,
            subset->singleton_expr_root,
            CompiledMathNodeKind::ExprCdf,
            condition_id,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Observed));
    root_id = compiled_math_make_root(&plan->compiled_math, node);
	  } else if (subset->scenarios.empty()) {
	    root_id =
	        compiled_math_make_root(
	            &plan->compiled_math,
	            compiled_math_constant(&plan->compiled_math, 0.0));
	  } else {
    std::vector<semantic::Index> scenario_nodes;
    scenario_nodes.reserve(subset->scenarios.size());
    for (const auto &scenario : subset->scenarios) {
      const auto scenario_root =
          compile_runtime_scenario_truth_cdf_root(
              plan,
              scenario,
              condition_id);
      scenario_nodes.push_back(
          compiled_math_root_node_id(plan->compiled_math, scenario_root));
    }
    const auto sum_node =
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Sum,
            std::move(scenario_nodes),
            CompiledMathValueKind::Cdf);
    root_id =
        compiled_math_make_root(
            &plan->compiled_math,
            compiled_math_unary_node(
                &plan->compiled_math,
                CompiledMathNodeKind::ClampProbability,
                sum_node,
                CompiledMathValueKind::Cdf));
  }
  subset->conditioned_win_cdf_roots.push_back(
      ExactConditionedRoot{condition_id, root_id});
  return root_id;
}

inline void compile_competitor_conditioned_roots(
    ExactVariantPlan *plan,
    ExactRuntimeOutcomePlan *runtime_outcome,
    const semantic::Index condition_id) {
  if (condition_id == semantic::kInvalidIndex || condition_id == 0) {
    return;
  }
  for (auto &block : runtime_outcome->competitor_blocks) {
    for (auto &subset : block.subsets) {
      compile_competitor_subset_win_cdf_root(plan, &subset, condition_id);
    }
  }
}

inline ExactRuntimeScenarioFormula compile_runtime_scenario_formula(
    ExactVariantPlan *plan,
    const ExactTransitionScenario &scenario) {
  ExactRuntimeScenarioFormula formula;
  formula.active_source_id = scenario.active_source_id;
  formula.active_observed_condition_id =
      compile_source_exact_condition_id(
          plan,
          scenario.active_source_id,
          CompiledMathTimeSlot::Observed);
  formula.relation_template = scenario.relation_template;
  formula.source_view_id =
      compile_source_view_id(plan, scenario.relation_template);
  formula.source_order_facts = scenario.source_order_facts;
  formula.has_readiness =
      !scenario.before_source_span.empty() || !scenario.ready_expr_span.empty();
  const auto relation_condition_id =
      compile_relation_condition_id(plan, scenario.relation_template);
  const auto source_condition_id = formula.active_observed_condition_id;
  const auto expr_condition_id =
      compiled_math_merge_conditions(
          &plan->compiled_math,
          {relation_condition_id, source_condition_id});
  formula.readiness_cdf = compile_runtime_readiness_cdf(*plan, scenario);
  formula.readiness_density = compile_runtime_readiness_density(*plan, scenario);
  formula.after_survival = compile_runtime_after_survival(*plan, scenario);
  formula.tail_condition = runtime_scenario_tail_condition(formula);
  compile_runtime_condition_ids(plan, &formula.tail_condition);
  compile_runtime_truth_formula_math(
      plan,
      &formula.readiness_cdf,
      source_condition_id,
      expr_condition_id,
      static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
      formula.source_view_id);
  compile_runtime_truth_formula_math(
      plan,
      &formula.readiness_density,
      source_condition_id,
      expr_condition_id,
      static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
      formula.source_view_id);
  compile_runtime_truth_formula_math(
      plan,
      &formula.after_survival,
      source_condition_id,
      expr_condition_id,
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
      formula.source_view_id);
  return formula;
}

inline void compile_runtime_scenario_contextual_conditions(
    ExactVariantPlan *plan,
    ExactRuntimeOutcomePlan *runtime_outcome,
    ExactRuntimeScenarioFormula *formula) {
  formula->tail_competitor_condition =
      filter_runtime_condition_for_competitors(
          formula->tail_condition,
          *runtime_outcome,
          *plan);
  compile_runtime_condition_ids(plan, &formula->tail_competitor_condition);
  compile_competitor_conditioned_roots(
      plan,
      runtime_outcome,
      formula->tail_competitor_condition.observed_condition_id);
  formula->has_tail_competitor_condition =
      !runtime_condition_empty(formula->tail_competitor_condition);
  formula->tail_competitor_subset_mask =
      runtime_competitor_subset_mask(
          formula->tail_competitor_condition,
          *runtime_outcome,
          *plan);
  formula->has_conditioned_readiness_terms = false;
  formula->condition_groups.clear();
  for (auto &term : formula->readiness_density.sum_terms) {
    term.condition_impossible =
        runtime_condition_order_contradiction(term.condition);
    compile_runtime_condition_ids(plan, &term.condition);
    term.tail_condition = filter_runtime_condition_for_truth(
        term.condition,
        formula->after_survival,
        *plan);
    compile_runtime_condition_ids(plan, &term.tail_condition);
    term.competitor_condition = filter_runtime_condition_for_competitors(
        term.condition,
        *runtime_outcome,
        *plan);
    compile_runtime_condition_ids(plan, &term.competitor_condition);
    const auto term_combined_competitor_condition_id =
        compile_combined_runtime_condition_id(
            plan,
            term.competitor_condition,
            formula->tail_competitor_condition);
    compile_competitor_conditioned_roots(
        plan,
        runtime_outcome,
        term_combined_competitor_condition_id);
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
                  *runtime_outcome,
                  *plan),
              term_combined_competitor_condition_id});
      compile_runtime_condition_ids(
          plan,
          &formula->condition_groups.back().tail_condition);
      compile_runtime_condition_ids(
          plan,
          &formula->condition_groups.back().competitor_condition);
      compile_competitor_conditioned_roots(
          plan,
          runtime_outcome,
          formula->condition_groups.back().combined_competitor_condition_id);
    }
  }
}

inline semantic::Index append_runtime_kernel_context(
    ExactVariantPlan *plan,
    ExactRuntimeScenarioExecutionKernel *kernel,
    const std::initializer_list<semantic::Index> condition_ids) {
  ExactRuntimeConditionContextPlan context;
  context.condition_id =
      compiled_math_merge_conditions(&plan->compiled_math, condition_ids);
  context.impossible =
      compiled_condition_impossible(plan->compiled_math, context.condition_id);
  for (semantic::Index i = 0;
       i < static_cast<semantic::Index>(kernel->condition_contexts.size());
       ++i) {
    if (kernel->condition_contexts[static_cast<std::size_t>(i)] == context) {
      return i;
    }
  }
  const auto index =
      static_cast<semantic::Index>(kernel->condition_contexts.size());
  kernel->condition_contexts.push_back(std::move(context));
  return index;
}

inline void compile_runtime_kernel_context_roots(
    ExactVariantPlan *plan,
    ExactRuntimeScenarioExecutionKernel *kernel,
    const ExactRuntimeScenarioFormula &formula) {
  const auto relation_condition_id =
      compile_relation_condition_id(plan, formula.relation_template);
  for (auto &context : kernel->condition_contexts) {
    if (context.impossible) {
      continue;
    }
    const auto expr_condition_id =
        compiled_math_merge_conditions(
            &plan->compiled_math,
            {relation_condition_id, context.condition_id});
    context.after_survival_root_id =
        compile_runtime_truth_formula_root(
            plan,
            formula.after_survival,
            context.condition_id,
            expr_condition_id,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Observed),
            formula.source_view_id);
  }
}

inline void compile_runtime_scenario_execution_kernel(
    ExactVariantPlan *plan,
    ExactRuntimeScenarioFormula *formula) {
  auto &kernel = formula->execution_kernel;
  kernel = ExactRuntimeScenarioExecutionKernel{};
  kernel.has_readiness = formula->has_readiness;
  kernel.has_tail_competitor_condition =
      formula->has_tail_competitor_condition;
  kernel.has_conditioned_readiness_terms =
      formula->has_conditioned_readiness_terms;
  kernel.target_scenario_context_index =
      append_runtime_kernel_context(
          plan,
          &kernel,
          {formula->active_observed_condition_id});
  kernel.competitor_scenario_context_index =
      append_runtime_kernel_context(
          plan,
          &kernel,
          {formula->active_observed_condition_id});
  kernel.tail_competitor_context_index =
      formula->has_tail_competitor_condition
          ? append_runtime_kernel_context(
                plan,
                &kernel,
                {formula->active_observed_condition_id,
                 formula->tail_competitor_condition.observed_condition_id})
          : kernel.competitor_scenario_context_index;

  kernel.condition_groups.reserve(formula->condition_groups.size());
  for (semantic::Index group_idx = 0;
       group_idx < static_cast<semantic::Index>(
                       formula->condition_groups.size());
       ++group_idx) {
    const auto &group =
        formula->condition_groups[static_cast<std::size_t>(group_idx)];
    const bool has_tail_condition =
        !runtime_condition_empty(group.tail_condition);
    const bool has_competitor_condition =
        !runtime_condition_empty(group.competitor_condition);
    const auto target_tail_context =
        has_tail_condition
            ? append_runtime_kernel_context(
                  plan,
                  &kernel,
                  {formula->active_observed_condition_id,
                   group.tail_condition.readiness_condition_id})
            : kernel.target_scenario_context_index;
    const auto competitor_context =
        (has_competitor_condition ||
         formula->has_tail_competitor_condition)
            ? append_runtime_kernel_context(
                  plan,
                  &kernel,
                  {formula->active_observed_condition_id,
                   has_competitor_condition
                       ? group.competitor_condition.readiness_condition_id
                       : semantic::kInvalidIndex,
                   formula->has_tail_competitor_condition
                       ? formula->tail_competitor_condition.observed_condition_id
                       : semantic::kInvalidIndex})
            : kernel.competitor_scenario_context_index;
	    kernel.condition_groups.push_back(
	        ExactRuntimeConditionGroupKernel{
	            group_idx,
	            target_tail_context,
	            competitor_context,
	            group.combined_competitor_condition_id,
	            semantic::kInvalidIndex,
	            has_tail_condition,
	            has_competitor_condition});
  }

  for (semantic::Index term_idx = 0;
       term_idx < static_cast<semantic::Index>(
                      formula->readiness_density.sum_terms.size());
       ++term_idx) {
    const auto &term =
        formula->readiness_density.sum_terms[static_cast<std::size_t>(term_idx)];
    if (!term.condition_impossible) {
      kernel.readiness_terms.push_back(
          ExactRuntimeReadinessTermKernel{
              term_idx,
              term.context_group_index});
    }
  }
  compile_runtime_kernel_context_roots(plan, &kernel, *formula);
}

inline semantic::Index compile_runtime_competitor_subset_precedence_root(
    ExactVariantPlan *plan,
    ExactRuntimeCompetitorSubsetPlan *subset,
    const semantic::Index condition_id,
    const ExactRuntimeScenarioFormula &target_formula,
    const semantic::Index ready_upper_time_id) {
  if (condition_id == semantic::kInvalidIndex ||
      compiled_condition_impossible(plan->compiled_math, condition_id)) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
  }
  if (subset->scenarios.empty()) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
  }

  const bool has_same_active =
      std::any_of(
          subset->scenarios.begin(),
          subset->scenarios.end(),
          [&](const ExactRuntimeScenarioFormula &scenario) {
            return scenario.active_source_id == target_formula.active_source_id;
          });
  if (!has_same_active) {
    return compile_competitor_subset_win_cdf_root(
        plan,
        subset,
        condition_id);
  }

  std::vector<semantic::Index> scenario_nodes;
  scenario_nodes.reserve(subset->scenarios.size());
  for (const auto &scenario : subset->scenarios) {
    const auto root_id =
        scenario.active_source_id == target_formula.active_source_id
            ? compile_runtime_scenario_same_active_root(
                  plan,
                  scenario,
                  condition_id,
                  ready_upper_time_id,
                  &target_formula.after_survival)
            : compile_runtime_scenario_truth_cdf_root(
                  plan,
                  scenario,
                  condition_id);
    const auto node_id = compiled_math_root_node_id(plan->compiled_math, root_id);
    if (node_id != semantic::kInvalidIndex) {
      scenario_nodes.push_back(node_id);
    }
  }

  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::CleanSignedSum,
              std::move(scenario_nodes),
              CompiledMathValueKind::Cdf),
          CompiledMathValueKind::Cdf));
}

inline semantic::Index compile_runtime_competitor_non_win_root(
    ExactVariantPlan *plan,
    const ExactRuntimeOutcomePlan &runtime_outcome,
    const ExactRuntimeCompetitorNonWinPlan &non_win_plan) {
  if (non_win_plan.impossible) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 0.0));
  }
  if (non_win_plan.blocks.empty()) {
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_constant(&plan->compiled_math, 1.0));
  }

  std::vector<semantic::Index> block_non_win_nodes;
  block_non_win_nodes.reserve(non_win_plan.blocks.size());
  for (const auto &block : non_win_plan.blocks) {
    std::vector<semantic::Index> subset_terms;
    subset_terms.reserve(static_cast<std::size_t>(block.subset_span.size));
    for (semantic::Index i = 0; i < block.subset_span.size; ++i) {
      const auto &subset =
          non_win_plan.subsets[
              static_cast<std::size_t>(block.subset_span.offset + i)];
      if (subset.block_index == semantic::kInvalidIndex ||
          subset.subset_index == semantic::kInvalidIndex ||
          subset.win_cdf_root_id == semantic::kInvalidIndex) {
        continue;
      }
      const auto &runtime_subset =
          runtime_outcome
              .competitor_blocks[static_cast<std::size_t>(subset.block_index)]
              .subsets[static_cast<std::size_t>(subset.subset_index)];
      const auto subset_mass =
          compiled_math_unary_node(
              &plan->compiled_math,
              CompiledMathNodeKind::ClampProbability,
              compiled_math_root_node_id(
                  plan->compiled_math,
                  subset.win_cdf_root_id),
              CompiledMathValueKind::Cdf);
      const auto gate =
          compile_outcome_subset_unused_node(
              plan,
              runtime_subset.outcome_indices);
      const auto gated_mass =
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::vector<semantic::Index>{gate, subset_mass},
              CompiledMathValueKind::Cdf);
      subset_terms.push_back(
          compile_signed_term_node(
              plan,
              gated_mass,
              subset.inclusion_sign));
    }
    const auto union_node =
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::ClampProbability,
            compiled_math_algebra_node(
                &plan->compiled_math,
                CompiledMathNodeKind::CleanSignedSum,
                std::move(subset_terms),
                CompiledMathValueKind::Cdf),
            CompiledMathValueKind::Cdf);
    block_non_win_nodes.push_back(
        compiled_math_unary_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Complement,
            union_node,
            CompiledMathValueKind::Cdf));
  }

  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_unary_node(
          &plan->compiled_math,
          CompiledMathNodeKind::ClampProbability,
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::move(block_non_win_nodes)),
          CompiledMathValueKind::Cdf));
}

inline semantic::Index append_runtime_competitor_non_win_plan(
    ExactVariantPlan *plan,
    ExactRuntimeOutcomePlan *runtime_outcome,
    ExactRuntimeScenarioExecutionKernel *kernel,
    const ExactRuntimeScenarioFormula &target_formula,
    const semantic::Index context_index,
    const semantic::Index ready_upper_time_id,
    const ExactRuntimeCompetitorSubsetMask *mask_a,
    const ExactRuntimeCompetitorSubsetMask *mask_b) {
  const auto &context_plan =
      kernel->condition_contexts[static_cast<std::size_t>(context_index)];
  ExactRuntimeCompetitorNonWinPlan non_win_plan;
  non_win_plan.impossible = context_plan.impossible;
  const auto context_condition_id = context_plan.condition_id;
  const auto base_condition_id =
      kernel->condition_contexts[
          static_cast<std::size_t>(
              kernel->competitor_scenario_context_index)].condition_id;
  const bool has_subset_mask = mask_a != nullptr || mask_b != nullptr;

  non_win_plan.blocks.reserve(runtime_outcome->competitor_blocks.size());
  for (std::size_t block_idx = 0;
       block_idx < runtime_outcome->competitor_blocks.size();
       ++block_idx) {
    auto &block = runtime_outcome->competitor_blocks[block_idx];
    const auto subset_offset =
        static_cast<semantic::Index>(non_win_plan.subsets.size());
    for (std::size_t subset_idx = 0; subset_idx < block.subsets.size();
         ++subset_idx) {
      auto &subset = block.subsets[subset_idx];
      const bool affected =
          !has_subset_mask ||
          runtime_competitor_subset_mask_affected(mask_a, block_idx, subset_idx) ||
          runtime_competitor_subset_mask_affected(mask_b, block_idx, subset_idx);
      const auto subset_condition_id =
          affected ? context_condition_id : base_condition_id;
      non_win_plan.subsets.push_back(
          ExactRuntimeCompetitorNonWinSubsetPlan{
              static_cast<semantic::Index>(block_idx),
              static_cast<semantic::Index>(subset_idx),
              subset.inclusion_sign,
              compile_runtime_competitor_subset_precedence_root(
                  plan,
                  &subset,
                  subset_condition_id,
                  target_formula,
                  ready_upper_time_id)});
    }
    non_win_plan.blocks.push_back(
        ExactRuntimeCompetitorNonWinBlockPlan{
            ExactIndexSpan{
                subset_offset,
                static_cast<semantic::Index>(
                    non_win_plan.subsets.size()) -
                    subset_offset}});
  }

  const auto index =
      static_cast<semantic::Index>(kernel->competitor_non_win_plans.size());
  non_win_plan.root_id =
      compile_runtime_competitor_non_win_root(
          plan,
          *runtime_outcome,
          non_win_plan);
  kernel->competitor_non_win_plans.push_back(std::move(non_win_plan));
  return index;
}

inline void compile_runtime_scenario_competitor_non_win_plans(
    ExactVariantPlan *plan,
    ExactRuntimeOutcomePlan *runtime_outcome,
    ExactRuntimeScenarioFormula *formula) {
  auto &kernel = formula->execution_kernel;
  kernel.base_competitor_non_win_plan_index =
      append_runtime_competitor_non_win_plan(
          plan,
          runtime_outcome,
          &kernel,
          *formula,
          kernel.competitor_scenario_context_index,
          static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
          nullptr,
          nullptr);
  const auto *tail_mask =
      kernel.has_tail_competitor_condition
          ? &formula->tail_competitor_subset_mask
          : nullptr;
  kernel.tail_competitor_non_win_plan_index =
      append_runtime_competitor_non_win_plan(
          plan,
          runtime_outcome,
          &kernel,
          *formula,
          kernel.tail_competitor_context_index,
          static_cast<semantic::Index>(CompiledMathTimeSlot::Zero),
          nullptr,
          tail_mask);
  for (auto &group_kernel : kernel.condition_groups) {
    const auto &group = formula->condition_groups[
        static_cast<std::size_t>(group_kernel.group_index)];
    group_kernel.competitor_non_win_plan_index =
        append_runtime_competitor_non_win_plan(
            plan,
            runtime_outcome,
            &kernel,
            *formula,
            group_kernel.competitor_context_index,
            static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness),
            group_kernel.has_competitor_condition
                ? &group.competitor_subset_mask
                : nullptr,
            tail_mask);
  }
}

struct CompiledMathSourceViewCloneMemo {
  struct Entry {
    semantic::Index original_id{semantic::kInvalidIndex};
    semantic::Index source_view_id{0};
    semantic::Index cloned_id{semantic::kInvalidIndex};
  };

  std::vector<Entry> nodes;
  std::vector<Entry> roots;
  std::vector<Entry> conditions;
};

inline bool compose_source_view_template(
    const ExactRelationTemplate &base,
    const ExactRelationTemplate &overlay,
    ExactRelationTemplate *out) {
  *out = base;
  for (std::size_t i = 0; i < overlay.source_ids.size(); ++i) {
    const auto source_id = overlay.source_ids[i];
    const auto relation = overlay.relations[i];
    const auto it =
        std::lower_bound(out->source_ids.begin(), out->source_ids.end(), source_id);
    const auto pos = static_cast<std::size_t>(it - out->source_ids.begin());
    if (it != out->source_ids.end() && *it == source_id) {
      if (out->relations[pos] != relation) {
        return false;
      }
      continue;
    }
    out->source_ids.insert(it, source_id);
    out->relations.insert(out->relations.begin() + static_cast<std::ptrdiff_t>(pos), relation);
  }
  return true;
}

inline bool composed_source_view_id(
    ExactVariantPlan *plan,
    const semantic::Index base_source_view_id,
    const semantic::Index overlay_source_view_id,
    semantic::Index *out_source_view_id) {
  if (overlay_source_view_id == 0 ||
      overlay_source_view_id == semantic::kInvalidIndex) {
    *out_source_view_id =
        base_source_view_id == semantic::kInvalidIndex ? 0 : base_source_view_id;
    return true;
  }
  if (base_source_view_id == 0 ||
      base_source_view_id == semantic::kInvalidIndex) {
    *out_source_view_id = overlay_source_view_id;
    return true;
  }
  const auto base_pos = static_cast<std::size_t>(base_source_view_id - 1U);
  const auto overlay_pos = static_cast<std::size_t>(overlay_source_view_id - 1U);
  if (base_pos >= plan->compiled_source_views.size() ||
      overlay_pos >= plan->compiled_source_views.size()) {
    return false;
  }
  ExactRelationTemplate composed;
  if (!compose_source_view_template(
          plan->compiled_source_views[base_pos],
          plan->compiled_source_views[overlay_pos],
          &composed)) {
    return false;
  }
  *out_source_view_id = compile_source_view_id(plan, composed);
  return true;
}

inline semantic::Index clone_compiled_math_node_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index node_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo);

inline semantic::Index clone_compiled_math_root_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index root_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (root_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->roots) {
    if (entry.original_id == root_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto root_pos = static_cast<std::size_t>(root_id);
  if (root_pos >= plan->compiled_math.roots.size()) {
    return semantic::kInvalidIndex;
  }
  const auto cloned_node_id =
      clone_compiled_math_node_for_source_view(
          plan,
          plan->compiled_math.roots[root_pos].node_id,
          normalized_source_view_id,
          memo);
  const auto cloned_root_id =
      compiled_math_make_root(&plan->compiled_math, cloned_node_id);
  memo->roots.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          root_id, normalized_source_view_id, cloned_root_id});
  return cloned_root_id;
}

inline semantic::Index clone_compiled_math_condition_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index condition_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return 0;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->conditions) {
    if (entry.original_id == condition_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto condition_pos = static_cast<std::size_t>(condition_id - 1U);
  if (condition_pos >= plan->compiled_math.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &condition = plan->compiled_math.conditions[condition_pos];
  CompiledMathConditionKey key;
  key.impossible = condition.impossible;
  key.source_ids = condition.source_ids;
  key.relations = condition.relations;
  key.fact_kinds = condition.fact_kinds;
  key.fact_subject_ids = condition.fact_subject_ids;
  key.fact_aux_ids = condition.fact_aux_ids;
  key.fact_aux2_ids = condition.fact_aux2_ids;
  key.fact_time_ids = condition.fact_time_ids;
  key.fact_normalizer_node_ids.reserve(
      condition.fact_normalizer_node_ids.size());
  for (const auto normalizer_node_id : condition.fact_normalizer_node_ids) {
    key.fact_normalizer_node_ids.push_back(
        normalizer_node_id == semantic::kInvalidIndex
            ? semantic::kInvalidIndex
            : clone_compiled_math_node_for_source_view(
                  plan,
                  normalizer_node_id,
                  normalized_source_view_id,
                  memo));
  }
  key.fact_normalizer_root_ids.reserve(
      condition.fact_normalizer_root_ids.size());
  for (const auto normalizer_root_id : condition.fact_normalizer_root_ids) {
    key.fact_normalizer_root_ids.push_back(
        normalizer_root_id == semantic::kInvalidIndex
            ? semantic::kInvalidIndex
            : clone_compiled_math_root_for_source_view(
                  plan,
                  normalizer_root_id,
                  normalized_source_view_id,
                  memo));
  }
  const auto cloned_condition_id =
      compiled_math_intern_condition(&plan->compiled_math, std::move(key));
  memo->conditions.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          condition_id, normalized_source_view_id, cloned_condition_id});
  return cloned_condition_id;
}

inline semantic::Index clone_compiled_math_node_for_source_view(
    ExactVariantPlan *plan,
    const semantic::Index node_id,
    const semantic::Index source_view_id,
    CompiledMathSourceViewCloneMemo *memo) {
  if (node_id == semantic::kInvalidIndex) {
    return semantic::kInvalidIndex;
  }
  const auto normalized_source_view_id =
      source_view_id == semantic::kInvalidIndex ? 0 : source_view_id;
  for (const auto &entry : memo->nodes) {
    if (entry.original_id == node_id &&
        entry.source_view_id == normalized_source_view_id) {
      return entry.cloned_id;
    }
  }
  const auto node_pos = static_cast<std::size_t>(node_id);
  if (node_pos >= plan->compiled_math.nodes.size()) {
    return semantic::kInvalidIndex;
  }
  const auto &node = plan->compiled_math.nodes[node_pos];
  semantic::Index composed_view_id{0};
  if (!composed_source_view_id(
          plan,
          normalized_source_view_id,
          node.source_view_id,
          &composed_view_id)) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }

  CompiledMathNodeKey key;
  key.kind = node.kind;
  key.value_kind = node.value_kind;
  key.subject_id = node.subject_id;
  key.condition_id =
      clone_compiled_math_condition_for_source_view(
          plan,
          node.condition_id,
          normalized_source_view_id,
          memo);
  key.time_id = node.time_id;
  key.aux_id = node.aux_id;
  key.aux2_id = node.aux2_id;
  key.source_view_id = composed_view_id;
  key.constant = node.constant;
  key.children.reserve(static_cast<std::size_t>(node.children.size));
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    key.children.push_back(
        clone_compiled_math_node_for_source_view(
            plan,
            plan->compiled_math.child_nodes[
                static_cast<std::size_t>(node.children.offset + i)],
            normalized_source_view_id,
            memo));
  }
  if (node.kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
      node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw) {
    key.subject_id =
        clone_compiled_math_root_for_source_view(
            plan,
            node.subject_id,
            normalized_source_view_id,
            memo);
  }
  const auto cloned_node_id =
      compiled_math_intern_node(&plan->compiled_math, std::move(key));
  memo->nodes.push_back(
      CompiledMathSourceViewCloneMemo::Entry{
          node_id, normalized_source_view_id, cloned_node_id});
  return cloned_node_id;
}

inline semantic::Index compile_runtime_root_value_node(
    ExactVariantPlan *plan,
    const semantic::Index root_id,
    const semantic::Index source_view_id = 0,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  if (root_id == semantic::kInvalidIndex) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if ((source_view_id == 0 || source_view_id == semantic::kInvalidIndex) &&
      static_cast<std::size_t>(root_id) < plan->compiled_math.roots.size()) {
    return plan->compiled_math.roots[static_cast<std::size_t>(root_id)].node_id;
  }
  (void)value_kind;
  CompiledMathSourceViewCloneMemo memo;
  const auto cloned_root_id =
      clone_compiled_math_root_for_source_view(
          plan,
          root_id,
          source_view_id,
          &memo);
  if (cloned_root_id == semantic::kInvalidIndex) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  return plan->compiled_math.roots[
      static_cast<std::size_t>(cloned_root_id)].node_id;
}

inline semantic::Index compile_runtime_competitor_non_win_value_node(
    ExactVariantPlan *plan,
    const ExactRuntimeScenarioExecutionKernel &kernel,
    const semantic::Index plan_index,
    const semantic::Index target_source_view_id) {
  if (plan_index == semantic::kInvalidIndex ||
      static_cast<std::size_t>(plan_index) >=
          kernel.competitor_non_win_plans.size()) {
    throw std::runtime_error(
        "scenario probability root is missing a competitor non-win plan");
  }
  const auto &non_win_plan =
      kernel.competitor_non_win_plans[static_cast<std::size_t>(plan_index)];
  if (non_win_plan.impossible) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  if (non_win_plan.root_id == semantic::kInvalidIndex) {
    throw std::runtime_error(
        "scenario probability root reached an uncompiled competitor plan");
  }
  return compile_runtime_root_value_node(
      plan,
      non_win_plan.root_id,
      target_source_view_id,
      CompiledMathValueKind::Cdf);
}

inline semantic::Index compile_runtime_context_tail_node(
    ExactVariantPlan *plan,
    const ExactRuntimeScenarioFormula &formula,
    const semantic::Index context_index) {
  const auto &kernel = formula.execution_kernel;
  if (context_index == semantic::kInvalidIndex ||
      static_cast<std::size_t>(context_index) >=
          kernel.condition_contexts.size()) {
    throw std::runtime_error(
        "scenario probability root is missing a tail context");
  }
  const auto &context =
      kernel.condition_contexts[static_cast<std::size_t>(context_index)];
  if (context.impossible ||
      context.after_survival_root_id == semantic::kInvalidIndex) {
    return compiled_math_constant(&plan->compiled_math, 0.0);
  }
  return compile_runtime_root_value_node(
      plan,
      context.after_survival_root_id,
      0,
      CompiledMathValueKind::Survival);
}

inline semantic::Index compile_runtime_scenario_probability_root(
    ExactVariantPlan *plan,
    ExactRuntimeScenarioFormula *formula) {
  const auto &kernel = formula->execution_kernel;
  const auto observed_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Observed);
  const auto readiness_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Readiness);
  const auto zero_time_id =
      static_cast<semantic::Index>(CompiledMathTimeSlot::Zero);

  const auto active_pdf_root =
      compiled_math_make_root(
          &plan->compiled_math,
          compiled_math_source_node(
              &plan->compiled_math,
              CompiledMathNodeKind::SourcePdf,
              formula->active_source_id,
              0,
              observed_time_id));
  const auto active_pdf_node =
      compile_runtime_root_value_node(
          plan,
          active_pdf_root,
          0,
          CompiledMathValueKind::Pdf);
  const auto base_tail_node =
      compile_runtime_context_tail_node(
          plan,
          *formula,
          kernel.target_scenario_context_index);

  if (!kernel.has_readiness) {
    const auto non_win_node =
        compile_runtime_competitor_non_win_value_node(
            plan,
            kernel,
            kernel.tail_competitor_non_win_plan_index,
            formula->source_view_id);
    return compiled_math_make_root(
        &plan->compiled_math,
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{
                active_pdf_node,
                base_tail_node,
                non_win_node},
            CompiledMathValueKind::Scalar));
  }

  std::vector<semantic::Index> branches;
  branches.reserve(2U);

  const auto relation_condition_id =
      compile_relation_condition_id(plan, formula->relation_template);
  const auto expr_condition_id =
      compiled_math_merge_conditions(
          &plan->compiled_math,
          {relation_condition_id, formula->active_observed_condition_id});
  const auto initial_ready_root =
      compile_runtime_truth_formula_root(
          plan,
          formula->readiness_cdf,
          formula->active_observed_condition_id,
          expr_condition_id,
          zero_time_id,
          formula->source_view_id);
  const auto initial_ready_node =
      compile_runtime_root_value_node(
          plan,
          initial_ready_root,
          0,
          CompiledMathValueKind::Cdf);
  const auto initial_non_win_node =
      compile_runtime_competitor_non_win_value_node(
          plan,
          kernel,
          kernel.tail_competitor_non_win_plan_index,
          formula->source_view_id);
  branches.push_back(
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::Product,
          std::vector<semantic::Index>{
              initial_ready_node,
              active_pdf_node,
              base_tail_node,
              initial_non_win_node},
          CompiledMathValueKind::Scalar));

  semantic::Index integrand_node{semantic::kInvalidIndex};
  if (!kernel.has_conditioned_readiness_terms &&
      !kernel.has_tail_competitor_condition) {
    const auto readiness_density_node =
        compile_runtime_root_value_node(
            plan,
            formula->readiness_density.compiled_root_id,
            0,
            CompiledMathValueKind::Density);
    const auto non_win_node =
        compile_runtime_competitor_non_win_value_node(
            plan,
            kernel,
            kernel.base_competitor_non_win_plan_index,
            formula->source_view_id);
    integrand_node =
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{
                readiness_density_node,
                non_win_node},
            CompiledMathValueKind::Scalar);
    const auto integrand_root =
        compiled_math_make_root(&plan->compiled_math, integrand_node);
    const auto integral_node =
        compiled_math_raw_integral_zero_to_current_node(
            &plan->compiled_math,
            integrand_root,
            0,
            observed_time_id,
            0,
            readiness_time_id);
    const auto integral_value_node =
        compile_runtime_root_value_node(
            plan,
            compiled_math_make_root(&plan->compiled_math, integral_node));
    branches.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{
                active_pdf_node,
                base_tail_node,
                integral_value_node},
            CompiledMathValueKind::Scalar));
  } else {
    std::vector<semantic::Index> term_nodes;
    term_nodes.reserve(kernel.readiness_terms.size());
    const auto &readiness_terms = formula->readiness_density.sum_terms;
    for (const auto &term_kernel : kernel.readiness_terms) {
      if (term_kernel.group_kernel_index == semantic::kInvalidIndex) {
        continue;
      }
      const auto term_pos =
          static_cast<std::size_t>(term_kernel.term_index);
      if (term_pos >= readiness_terms.size()) {
        throw std::runtime_error(
            "scenario probability root points outside readiness terms");
      }
      const auto &group_kernel =
          kernel.condition_groups[
              static_cast<std::size_t>(term_kernel.group_kernel_index)];
      const auto term_node =
          compile_runtime_root_value_node(
              plan,
              readiness_terms[term_pos].compiled_root_id,
              0,
              CompiledMathValueKind::Density);
      const auto tail_node =
          group_kernel.has_tail_condition
              ? compile_runtime_context_tail_node(
                    plan,
                    *formula,
                    group_kernel.target_tail_context_index)
              : base_tail_node;
      const auto non_win_node =
          compile_runtime_competitor_non_win_value_node(
              plan,
              kernel,
              group_kernel.competitor_non_win_plan_index,
              formula->source_view_id);
      term_nodes.push_back(
          compiled_math_algebra_node(
              &plan->compiled_math,
              CompiledMathNodeKind::Product,
              std::vector<semantic::Index>{
                  term_node,
                  tail_node,
                  non_win_node},
              CompiledMathValueKind::Scalar));
    }
    integrand_node =
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::CleanSignedSum,
            std::move(term_nodes),
            CompiledMathValueKind::Scalar);
    const auto integrand_root =
        compiled_math_make_root(&plan->compiled_math, integrand_node);
    const auto integral_node =
        compiled_math_raw_integral_zero_to_current_node(
            &plan->compiled_math,
            integrand_root,
            0,
            observed_time_id,
            0,
            readiness_time_id);
    const auto integral_value_node =
        compile_runtime_root_value_node(
            plan,
            compiled_math_make_root(&plan->compiled_math, integral_node));
    branches.push_back(
        compiled_math_algebra_node(
            &plan->compiled_math,
            CompiledMathNodeKind::Product,
            std::vector<semantic::Index>{
                active_pdf_node,
                integral_value_node},
            CompiledMathValueKind::Scalar));
  }

  return compiled_math_make_root(
      &plan->compiled_math,
      compiled_math_algebra_node(
          &plan->compiled_math,
          CompiledMathNodeKind::CleanSignedSum,
          std::move(branches),
          CompiledMathValueKind::Scalar));
}

inline semantic::Index compile_runtime_outcome_probability_root(
    ExactVariantPlan *plan,
    const ExactRuntimeOutcomePlan &runtime_outcome) {
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
    const ExactTransitionScenario &scenario) {
  for (semantic::Index i = 0; i < scenario.ready_expr_span.size; ++i) {
    const auto expr_id =
        plan->scenario_expr_ids[
            static_cast<std::size_t>(scenario.ready_expr_span.offset + i)];
    if (expr_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(expr_id) <
            plan->sequence_expr_upper_bound_used.size()) {
      plan->sequence_expr_upper_bound_used[
          static_cast<std::size_t>(expr_id)] = 1U;
    }
  }
}

inline void compile_sequence_expr_upper_bound_roots(ExactVariantPlan *plan) {
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
  for (const auto &target_plan : plan->competitor_plans) {
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
    ExactVariantPlan *plan) {
  const ExactVariantPlan &plan_ref = *plan;
  ExactRuntimeVariantPlan runtime;
  runtime.outcomes.reserve(plan_ref.outcomes.size());

  for (semantic::Index target_idx = 0;
       target_idx < static_cast<semantic::Index>(plan_ref.outcomes.size());
       ++target_idx) {
    const auto target_pos = static_cast<std::size_t>(target_idx);
    const auto &outcome = plan_ref.outcomes[target_pos];
    const auto &competitor_plan = plan_ref.competitor_plans[target_pos];

    ExactRuntimeOutcomePlan runtime_outcome;
    runtime_outcome.scenarios.reserve(outcome.scenarios.size());
    for (const auto &scenario : outcome.scenarios) {
      runtime_outcome.scenarios.push_back(
          compile_runtime_scenario_formula(
              plan,
              scenario));
    }

    runtime_outcome.competitor_blocks.reserve(competitor_plan.blocks.size());
    for (const auto &block : competitor_plan.blocks) {
      ExactRuntimeCompetitorBlockPlan runtime_block;
      runtime_block.subsets.reserve(block.subsets.size());
      for (const auto &subset : block.subsets) {
        ExactRuntimeCompetitorSubsetPlan runtime_subset;
        runtime_subset.outcome_indices = subset.outcome_indices;
        runtime_subset.expr_roots = subset.expr_roots;
        runtime_subset.inclusion_sign = subset.inclusion_sign;
        runtime_subset.singleton_expr_root = subset.singleton_expr_root;
        runtime_subset.scenarios.reserve(subset.scenarios.size());
        for (const auto &scenario : subset.scenarios) {
          runtime_subset.scenarios.push_back(
              compile_runtime_scenario_formula(
                  plan,
                  scenario));
        }
        runtime_subset.win_cdf_root_id =
            compile_competitor_subset_win_cdf_root(
                plan,
                &runtime_subset,
                0);
        runtime_block.subsets.push_back(std::move(runtime_subset));
      }
      runtime_outcome.competitor_blocks.push_back(std::move(runtime_block));
    }

    for (auto &scenario : runtime_outcome.scenarios) {
      compile_runtime_scenario_contextual_conditions(
          plan,
          &runtime_outcome,
          &scenario);
      compile_runtime_scenario_execution_kernel(
          plan,
          &scenario);
    }

    for (std::size_t scenario_idx = 0;
         scenario_idx < runtime_outcome.scenarios.size();
         ++scenario_idx) {
      compile_runtime_scenario_competitor_non_win_plans(
          plan,
          &runtime_outcome,
          &runtime_outcome.scenarios[scenario_idx]);
      runtime_outcome.scenarios[scenario_idx].probability_root_id =
          compile_runtime_scenario_probability_root(
              plan,
              &runtime_outcome.scenarios[scenario_idx]);
    }
    runtime_outcome.successor_distribution.total_probability_root_id =
        compile_runtime_outcome_probability_root(plan, runtime_outcome);
    runtime_outcome.successor_distribution.transitions.reserve(
        runtime_outcome.scenarios.size());
    for (std::size_t scenario_idx = 0;
         scenario_idx < runtime_outcome.scenarios.size();
         ++scenario_idx) {
      runtime_outcome.successor_distribution.transitions.push_back(
          ExactRuntimeScenarioTransitionPlan{
              runtime_outcome.scenarios[scenario_idx].probability_root_id,
              runtime_outcome.scenarios[scenario_idx].active_source_id,
              outcome.scenarios[scenario_idx].before_source_span,
              outcome.scenarios[scenario_idx].ready_expr_span,
              scenario_idx < outcome.scenarios.size() &&
                  scenario_supports_ranked_sequence(
                      outcome.scenarios[scenario_idx])});
    }

    runtime.outcomes.push_back(std::move(runtime_outcome));
  }

  return runtime;
}

inline void validate_compiled_math_has_no_interpreter_expr_nodes(
    const CompiledMathProgram &program) {
  for (std::size_t i = 0; i < program.nodes.size(); ++i) {
    const auto kind = program.nodes[i].kind;
    if (kind == CompiledMathNodeKind::ExprDensity ||
        kind == CompiledMathNodeKind::ExprCdf ||
        kind == CompiledMathNodeKind::ExprSurvival) {
      throw std::runtime_error(
          "exact compiled math contains an unlowered semantic node at " +
          std::to_string(i));
    }
  }
}

inline semantic::Index compiled_source_bound_plan_slot(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    const semantic::Index source_id) {
  if (source_id == semantic::kInvalidIndex ||
      program.source_condition_bound_source_count <= 0) {
    return semantic::kInvalidIndex;
  }
  const auto source_count =
      static_cast<std::size_t>(program.source_condition_bound_source_count);
  const auto source_pos = static_cast<std::size_t>(source_id);
  if (source_pos >= source_count) {
    return semantic::kInvalidIndex;
  }
  const auto condition_slot =
      condition_id == semantic::kInvalidIndex ? 0 : condition_id;
  const auto condition_pos = static_cast<std::size_t>(condition_slot);
  if (condition_pos > program.conditions.size()) {
    return semantic::kInvalidIndex;
  }
  const auto slot = condition_pos * source_count + source_pos;
  if (slot >= program.source_condition_bound_plans.size()) {
    return semantic::kInvalidIndex;
  }
  return static_cast<semantic::Index>(slot);
}

inline void compile_source_condition_bound_plans(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  const auto source_count = static_cast<std::size_t>(plan->source_count);
  program.source_condition_bound_source_count = plan->source_count;
  program.condition_source_relation_source_count = plan->source_count;
  program.source_condition_bound_plans.assign(
      (program.conditions.size() + 1U) * source_count,
      CompiledSourceBoundPlan{});
  program.source_condition_bound_terms.clear();
  program.condition_source_relations.assign(
      (program.conditions.size() + 1U) * source_count,
      static_cast<std::uint8_t>(ExactRelation::Unknown));
  if (source_count == 0U) {
    return;
  }
  auto append_bound_terms = [&program](
                                const CompiledMathCondition &condition,
                                const std::vector<semantic::Index> &fact_indices,
                                const CompiledMathIndexSpan fact_span) {
    const auto offset =
        static_cast<semantic::Index>(
            program.source_condition_bound_terms.size());
    for (semantic::Index i = 0; i < fact_span.size; ++i) {
      const auto fact_pos = static_cast<std::size_t>(
          fact_indices[
              static_cast<std::size_t>(fact_span.offset + i)]);
      const auto time_id =
          fact_pos < condition.fact_time_ids.size()
              ? condition.fact_time_ids[fact_pos]
              : static_cast<semantic::Index>(
                    CompiledMathTimeSlot::Observed);
      program.source_condition_bound_terms.push_back(
          CompiledSourceBoundTerm{time_id});
    }
    return CompiledMathIndexSpan{
        offset,
        static_cast<semantic::Index>(
            program.source_condition_bound_terms.size() -
            static_cast<std::size_t>(offset))};
  };
  for (std::size_t condition_pos = 0;
       condition_pos < program.conditions.size();
       ++condition_pos) {
    const auto &condition = program.conditions[condition_pos];
    const auto condition_id =
        static_cast<semantic::Index>(condition_pos + 1U);
    const auto relation_offset = (condition_pos + 1U) * source_count;
    for (std::size_t i = 0; i < condition.source_ids.size(); ++i) {
      const auto source_id = condition.source_ids[i];
      if (source_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(source_id) >= source_count) {
        continue;
      }
      program.condition_source_relations[
          relation_offset + static_cast<std::size_t>(source_id)] =
          condition.relations[i];
    }
    for (const auto &entry : condition.source_order_fact_lookup.entries) {
      const auto before_pos = static_cast<std::size_t>(entry.first);
      const auto after_pos = static_cast<std::size_t>(entry.second);
      if (entry.first == semantic::kInvalidIndex ||
          entry.second == semantic::kInvalidIndex ||
          before_pos >= source_count ||
          after_pos >= source_count) {
        continue;
      }
      program.condition_source_relations[relation_offset + before_pos] =
          static_cast<std::uint8_t>(ExactRelation::Before);
      program.condition_source_relations[relation_offset + after_pos] =
          static_cast<std::uint8_t>(ExactRelation::After);
    }
    for (std::size_t source_pos = 0; source_pos < source_count; ++source_pos) {
      const auto slot_id =
          compiled_source_bound_plan_slot(
              program,
              condition_id,
              static_cast<semantic::Index>(source_pos));
      if (slot_id == semantic::kInvalidIndex) {
        continue;
      }
      const auto slot = static_cast<std::size_t>(slot_id);
      auto &bounds = program.source_condition_bound_plans[slot];
      if (source_pos < condition.source_exact_fact_spans.size()) {
        bounds.exact =
            append_bound_terms(
                condition,
                condition.source_exact_fact_indices,
                condition.source_exact_fact_spans[source_pos]);
        bounds.has_condition_exact = !bounds.exact.empty();
      }
      if (source_pos < condition.source_lower_fact_spans.size()) {
        bounds.lower =
            append_bound_terms(
                condition,
                condition.source_lower_fact_indices,
                condition.source_lower_fact_spans[source_pos]);
        bounds.has_condition_lower = !bounds.lower.empty();
      }
      if (source_pos < condition.source_upper_fact_spans.size()) {
        bounds.upper =
            append_bound_terms(
                condition,
                condition.source_upper_fact_indices,
                condition.source_upper_fact_spans[source_pos]);
        bounds.has_condition_upper = !bounds.upper.empty();
      }
    }
  }
}

inline void compile_condition_cache_plans(CompiledMathProgram *program) {
  program->condition_cache_plans.clear();
  program->condition_cache_time_dependencies.clear();
  program->condition_cache_plans.reserve(program->conditions.size() + 1U);
  program->condition_cache_plans.push_back(CompiledConditionCachePlan{});
  for (std::size_t condition_pos = 0;
       condition_pos < program->conditions.size();
       ++condition_pos) {
    const auto condition_id =
        static_cast<semantic::Index>(condition_pos + 1U);
    const auto &condition = program->conditions[condition_pos];
    const auto offset =
        static_cast<semantic::Index>(
            program->condition_cache_time_dependencies.size());
    program->condition_cache_time_dependencies.insert(
        program->condition_cache_time_dependencies.end(),
        condition.time_dependency_ids.begin(),
        condition.time_dependency_ids.end());
    const auto size =
        static_cast<semantic::Index>(
            program->condition_cache_time_dependencies.size() -
            static_cast<std::size_t>(offset));
    program->condition_cache_plans.push_back(
        CompiledConditionCachePlan{
            CompiledMathIndexSpan{offset, size},
            compiled_math_static_condition_cache_id(condition_id),
            size > 0});
  }
}

inline void compile_source_channel_plans(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  compile_condition_cache_plans(&program);
  compile_source_condition_bound_plans(plan);
  program.source_channel_plans.clear();
  program.source_channel_plans.reserve(program.source_channel_keys.size());
  for (std::size_t request_pos = 0;
       request_pos < program.source_channel_keys.size();
       ++request_pos) {
    const auto &request = program.source_channel_keys[request_pos];
    CompiledSourceChannelPlan channel_plan;
    channel_plan.request = request;
    if (request_pos < program.source_channel_required_channels.size() &&
        program.source_channel_required_channels[request_pos] != 0U) {
      channel_plan.required_channels =
          program.source_channel_required_channels[request_pos];
    }
    if (request.source_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(request.source_id) <
            plan->source_kernels.size()) {
      channel_plan.source_kernel_slot = request.source_id;
      channel_plan.kernel =
          plan->source_kernels[static_cast<std::size_t>(request.source_id)]
              .kind;
    }
    channel_plan.bound_plan_slot =
        compiled_source_bound_plan_slot(
            program, request.condition_id, request.source_id);
    if (channel_plan.bound_plan_slot != semantic::kInvalidIndex) {
      channel_plan.bounds =
          program.source_condition_bound_plans[
              static_cast<std::size_t>(channel_plan.bound_plan_slot)];
    }
    channel_plan.has_source_condition_overlay =
        channel_plan.bounds.has_condition_exact ||
        channel_plan.bounds.has_condition_lower ||
        channel_plan.bounds.has_condition_upper;
    channel_plan.direct_leaf_absolute_candidate =
        request.source_id != semantic::kInvalidIndex &&
        channel_plan.source_kernel_slot != semantic::kInvalidIndex &&
        channel_plan.kernel == CompiledSourceChannelKernelKind::LeafAbsolute &&
        !channel_plan.has_source_condition_overlay;
    const auto condition_pos =
        request.condition_id == semantic::kInvalidIndex
            ? static_cast<std::size_t>(0U)
            : static_cast<std::size_t>(request.condition_id);
    if (condition_pos < program.condition_cache_plans.size()) {
      const auto &cache_plan = program.condition_cache_plans[condition_pos];
      channel_plan.condition_time_dependencies =
          cache_plan.time_dependencies;
      channel_plan.condition_static_cache_id =
          cache_plan.static_cache_id;
      channel_plan.condition_cache_dynamic =
          cache_plan.dynamic;
    }
    channel_plan.source_view_condition_static_cache_id =
        channel_plan.direct_leaf_absolute_candidate
            ? 0
            : compiled_math_static_source_view_condition_cache_id(
                  request.source_view_id,
                  request.condition_id);
    channel_plan.source_view_condition_cache_dynamic =
        !channel_plan.direct_leaf_absolute_candidate &&
        channel_plan.condition_cache_dynamic;
    program.source_channel_plans.push_back(channel_plan);
  }
}

inline void compile_source_product_channel_programs(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  for (auto &channel : program.integral_kernel_source_product_channels) {
    const auto slot = static_cast<std::size_t>(channel.source_channel_slot);
    if (channel.source_channel_slot == semantic::kInvalidIndex ||
        slot >= program.source_channel_plans.size()) {
      continue;
    }
    const auto &channel_plan = program.source_channel_plans[slot];
    channel.source_id = channel_plan.request.source_id;
    channel.condition_id = channel_plan.request.condition_id;
    channel.source_kernel_slot = channel_plan.source_kernel_slot;
    channel.kernel = channel_plan.kernel;
    if (channel.source_view_id != 0 &&
        channel.source_view_id != semantic::kInvalidIndex &&
        channel.source_id != semantic::kInvalidIndex) {
      channel.static_source_view_relation = static_cast<std::uint8_t>(
          exact_compiled_source_view_relation(
              *plan, channel.source_view_id, channel.source_id));
      channel.has_static_source_view_relation = true;
    }
    if (channel.source_kernel_slot != semantic::kInvalidIndex &&
        static_cast<std::size_t>(channel.source_kernel_slot) <
            plan->source_kernels.size()) {
      const auto &source_kernel =
          plan->source_kernels[
              static_cast<std::size_t>(channel.source_kernel_slot)];
      channel.leaf_index = source_kernel.leaf_index;
      if (channel.leaf_index != semantic::kInvalidIndex &&
          static_cast<std::size_t>(channel.leaf_index) <
              plan->lowered.program.leaf_descriptors.size()) {
        const auto &leaf =
            plan->lowered.program.leaf_descriptors[
                static_cast<std::size_t>(channel.leaf_index)];
        channel.leaf_dist_kind = leaf.dist_kind;
        channel.leaf_param_count = leaf.param_count;
        channel.leaf_onset_abs_value = leaf.onset_abs_value;
      }
    }
    channel.bounds = channel_plan.bounds;
    channel.condition_time_dependencies =
        channel_plan.condition_time_dependencies;
    channel.condition_static_cache_id = channel_plan.condition_static_cache_id;
    channel.direct_leaf_absolute_candidate =
        channel_plan.direct_leaf_absolute_candidate;
    channel.has_source_condition_overlay =
        channel_plan.has_source_condition_overlay;
    channel.condition_cache_dynamic = channel_plan.condition_cache_dynamic;
  }
}

inline CompiledMathSourceProductOpKind source_product_generic_op_kind(
    const CompiledMathNodeKind factor_kind) noexcept {
  switch (factor_kind) {
  case CompiledMathNodeKind::SourcePdf:
    return CompiledMathSourceProductOpKind::GenericPdf;
  case CompiledMathNodeKind::SourceCdf:
    return CompiledMathSourceProductOpKind::GenericCdf;
  case CompiledMathNodeKind::SourceSurvival:
    return CompiledMathSourceProductOpKind::GenericSurvival;
  default:
    break;
  }
  return CompiledMathSourceProductOpKind::ConstantZero;
}

inline CompiledMathSourceProductOpKind source_product_forced_relation_op_kind(
    const ExactRelation relation,
    const CompiledMathNodeKind factor_kind) noexcept {
  if (relation == ExactRelation::Before) {
    return factor_kind == CompiledMathNodeKind::SourceCdf
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  if (relation == ExactRelation::At) {
    if (factor_kind == CompiledMathNodeKind::SourcePdf) {
      return source_product_generic_op_kind(factor_kind);
    }
    return factor_kind == CompiledMathNodeKind::SourceCdf
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  if (relation == ExactRelation::After) {
    return factor_kind == CompiledMathNodeKind::SourceSurvival
               ? CompiledMathSourceProductOpKind::ConstantOne
               : CompiledMathSourceProductOpKind::ConstantZero;
  }
  return source_product_generic_op_kind(factor_kind);
}

inline CompiledMathSourceProductOpKind source_product_leaf_op_kind(
    const std::uint8_t leaf_dist_kind,
    const CompiledMathNodeKind factor_kind) noexcept {
  const auto dist_kind = static_cast<leaf::DistKind>(leaf_dist_kind);
  switch (dist_kind) {
  case leaf::DistKind::Lognormal:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafLognormalPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafLognormalCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafLognormalSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::Gamma:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafGammaPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafGammaCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafGammaSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::Exgauss:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafExgaussPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafExgaussCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafExgaussSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::LBA:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafLbaPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafLbaCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafLbaSurvival;
    default:
      break;
    }
    break;
  case leaf::DistKind::RDM:
    switch (factor_kind) {
    case CompiledMathNodeKind::SourcePdf:
      return CompiledMathSourceProductOpKind::LeafRdmPdf;
    case CompiledMathNodeKind::SourceCdf:
      return CompiledMathSourceProductOpKind::LeafRdmCdf;
    case CompiledMathNodeKind::SourceSurvival:
      return CompiledMathSourceProductOpKind::LeafRdmSurvival;
    default:
      break;
    }
    break;
  }
  return source_product_generic_op_kind(factor_kind);
}

inline CompiledMathSourceProductOpKind source_product_op_kind_for_factor(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathNodeKind factor_kind) noexcept {
  if (channel.has_static_source_view_relation &&
      !channel.has_source_condition_overlay) {
    const auto relation =
        static_cast<ExactRelation>(channel.static_source_view_relation);
    if (relation != ExactRelation::Unknown) {
      if (relation != ExactRelation::At ||
          factor_kind != CompiledMathNodeKind::SourcePdf) {
        return source_product_forced_relation_op_kind(relation, factor_kind);
      }
    }
  }
  if (channel.direct_leaf_absolute_candidate) {
    return source_product_leaf_op_kind(channel.leaf_dist_kind, factor_kind);
  }
  return source_product_generic_op_kind(factor_kind);
}

inline std::uint8_t source_product_op_fill_mask(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathSourceProductOpKind op_kind) noexcept {
  const auto requested =
      compiled_math_source_product_op_channel_mask(op_kind);
  if ((requested & (2U | 4U)) != 0U) {
    const auto paired =
        channel.required_channels & (2U | 4U);
    return paired == 0U ? requested : paired;
  }
  return requested;
}

inline CompiledMathIndexSpan compile_source_product_ops_for_factor_span(
    CompiledMathProgram *program,
    const CompiledMathIndexSpan factors) {
  const auto offset = static_cast<semantic::Index>(
      program->integral_kernel_source_product_ops.size());
  for (semantic::Index i = 0; i < factors.size; ++i) {
    const auto factor_pos =
        static_cast<std::size_t>(factors.offset + i);
    const auto &factor =
        program->integral_kernel_source_value_factors[factor_pos];
    const auto channel_pos =
        static_cast<std::size_t>(factor.source_product_channel_id);
    const auto &channel =
        program->integral_kernel_source_product_channels[channel_pos];
    const auto kind = source_product_op_kind_for_factor(channel, factor.kind);
    if (kind == CompiledMathSourceProductOpKind::ConstantZero) {
      program->integral_kernel_source_product_ops.resize(
          static_cast<std::size_t>(offset));
      program->integral_kernel_source_product_ops.push_back(
          CompiledMathSourceProductOp{
              kind,
              semantic::kInvalidIndex,
              0U,
              0U,
              0.0});
      return CompiledMathIndexSpan{offset, 1};
    }
    if (kind == CompiledMathSourceProductOpKind::ConstantOne) {
      continue;
    }
    const auto value_mask =
        compiled_math_source_product_op_channel_mask(kind);
    program->integral_kernel_source_product_channels[channel_pos]
        .scalar_op_count++;
    program->integral_kernel_source_product_ops.push_back(
        CompiledMathSourceProductOp{
            kind,
            factor.source_product_channel_id,
            value_mask,
            static_cast<std::uint8_t>(
                value_mask == 0U
                    ? 0U
                    : source_product_op_fill_mask(channel, kind)),
            kind == CompiledMathSourceProductOpKind::ConstantOne ? 1.0 : 0.0});
  }
  return CompiledMathIndexSpan{
      offset,
      static_cast<semantic::Index>(
          program->integral_kernel_source_product_ops.size() -
          static_cast<std::size_t>(offset))};
}

inline void compile_source_product_scalar_ops(ExactVariantPlan *plan) {
  auto &program = plan->compiled_math;
  program.integral_kernel_source_product_ops.clear();
  for (auto &channel : program.integral_kernel_source_product_channels) {
    channel.scalar_op_count = 0;
  }
  for (auto &kernel : program.integral_kernels) {
    if (kernel.kind == CompiledMathIntegralKernelKind::SourceProduct) {
      kernel.source_product_ops =
          compile_source_product_ops_for_factor_span(
              &program, kernel.source_value_factors);
      continue;
    }
    for (semantic::Index i = 0; i < kernel.source_product_terms.size; ++i) {
      auto &term =
          program.integral_kernel_source_product_terms[
              static_cast<std::size_t>(
                  kernel.source_product_terms.offset + i)];
      term.source_product_ops =
          compile_source_product_ops_for_factor_span(
              &program, term.source_value_factors);
    }
  }
  for (auto &op : program.integral_kernel_source_product_ops) {
    if (op.value_channel_mask == 0U) {
      op.cache_result = false;
      continue;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    op.cache_result = channel.scalar_op_count > 1 ||
                      op.fill_channel_mask != op.value_channel_mask;
  }
}

inline void compile_trigger_state_table(ExactVariantPlan *plan) {
  struct TriggerStateBuilder {
    double fixed_weight{1.0};
    std::vector<ExactCompiledTriggerWeightTerm> weight_terms;
    std::vector<std::uint8_t> shared_started;
  };

  auto &table = plan->trigger_state_table;
  table.states.clear();
  table.weight_terms.clear();
  table.shared_started_values.clear();
  table.trigger_count = plan->lowered.program.layout.n_triggers;

  std::vector<TriggerStateBuilder> builders;
  builders.push_back(TriggerStateBuilder{});
  builders.front().shared_started.assign(
      static_cast<std::size_t>(table.trigger_count), 2U);

  const auto &program = plan->lowered.program;
  for (const auto trigger_index : plan->shared_trigger_indices) {
    const auto trigger_pos = static_cast<std::size_t>(trigger_index);
    const bool fixed_q = program.trigger_has_fixed_q[trigger_pos] != 0U;
    const double raw_fixed_q =
        fixed_q ? program.trigger_fixed_q[trigger_pos] : 0.0;
    const double q_fixed =
        fixed_q
            ? (!std::isfinite(raw_fixed_q)
                   ? 0.0
                   : (raw_fixed_q <= 0.0
                          ? 0.0
                          : (raw_fixed_q >= 1.0 ? 1.0 : raw_fixed_q)))
            : 0.0;
    semantic::Index q_leaf_index{semantic::kInvalidIndex};
    if (!fixed_q) {
      const auto member_begin = program.trigger_member_offsets[trigger_pos];
      const auto member_end = program.trigger_member_offsets[trigger_pos + 1U];
      if (member_begin != member_end) {
        q_leaf_index =
            program.trigger_member_indices[
                static_cast<std::size_t>(member_begin)];
      }
    }

    std::vector<TriggerStateBuilder> next;
    next.reserve(builders.size() * 2U);
    for (const auto &builder : builders) {
      auto append_state =
          [&](const std::uint8_t shared_started, const double fixed_weight) {
            auto out = builder;
            out.fixed_weight *= fixed_weight;
            out.shared_started[trigger_pos] = shared_started;
            next.push_back(std::move(out));
          };
      auto append_variable_state =
          [&](const std::uint8_t shared_started) {
            auto out = builder;
            out.shared_started[trigger_pos] = shared_started;
            out.weight_terms.push_back(
                ExactCompiledTriggerWeightTerm{
                    q_leaf_index, shared_started});
            next.push_back(std::move(out));
          };

      if (fixed_q) {
        if (q_fixed > 0.0) {
          append_state(0U, q_fixed);
        }
        if (q_fixed < 1.0) {
          append_state(1U, 1.0 - q_fixed);
        }
      } else {
        append_variable_state(0U);
        append_variable_state(1U);
      }
    }
    builders.swap(next);
  }

  table.states.reserve(builders.size());
  table.shared_started_values.reserve(
      builders.size() * static_cast<std::size_t>(table.trigger_count));
  for (const auto &builder : builders) {
    const auto shared_offset =
        static_cast<semantic::Index>(table.shared_started_values.size());
    table.shared_started_values.insert(
        table.shared_started_values.end(),
        builder.shared_started.begin(),
        builder.shared_started.end());
    const auto weight_offset =
        static_cast<semantic::Index>(table.weight_terms.size());
    table.weight_terms.insert(
        table.weight_terms.end(),
        builder.weight_terms.begin(),
        builder.weight_terms.end());
    table.states.push_back(
        ExactCompiledTriggerState{
            builder.fixed_weight,
            ExactIndexSpan{
                weight_offset,
                static_cast<semantic::Index>(
                    table.weight_terms.size() -
                    static_cast<std::size_t>(weight_offset))},
            shared_offset});
  }
}

inline void compile_source_view_relation_tables(ExactVariantPlan *plan) {
  const auto source_count = static_cast<std::size_t>(plan->source_count);
  plan->compiled_source_view_source_count = plan->source_count;
  plan->compiled_source_view_relations.assign(
      plan->compiled_source_views.size() * source_count,
      static_cast<std::uint8_t>(ExactRelation::Unknown));
  if (source_count == 0U) {
    return;
  }
  for (std::size_t view_pos = 0;
       view_pos < plan->compiled_source_views.size();
       ++view_pos) {
    const auto &view = plan->compiled_source_views[view_pos];
    const auto view_offset = view_pos * source_count;
    for (std::size_t i = 0; i < view.source_ids.size(); ++i) {
      const auto source_id = view.source_ids[i];
      if (source_id == semantic::kInvalidIndex ||
          static_cast<std::size_t>(source_id) >= source_count) {
        continue;
      }
      plan->compiled_source_view_relations[
          view_offset + static_cast<std::size_t>(source_id)] =
          static_cast<std::uint8_t>(view.relations[i]);
    }
  }
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
  compile_source_kernels(&plan);
  compile_exact_expr_kernels(&plan);

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
  compile_trigger_state_table(&plan);

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
  plan.no_response = compile_no_response_plan(plan);
  for (auto &target_plan : plan.competitor_plans) {
    for (auto &block : target_plan.blocks) {
      for (auto &subset : block.subsets) {
        for (auto &scenario : subset.scenarios) {
          compile_scenario_runtime_fields(&plan, &scenario);
        }
      }
    }
  }
  compile_sequence_expr_upper_bound_roots(&plan);
  plan.runtime = compile_exact_runtime_plan(&plan);
  compile_source_channel_plans(&plan);
  compile_source_view_relation_tables(&plan);
  compile_source_product_channel_programs(&plan);
  compile_source_product_scalar_ops(&plan);
  plan.simple_race = compile_simple_race_plan(plan);
  plan.probability_programs = compile_exact_probability_programs(plan);
  validate_compiled_math_has_no_interpreter_expr_nodes(plan.compiled_math);
  compiled_math_release_planning_fields(&plan.compiled_math);
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
