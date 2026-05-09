#pragma once

#include "exact_types.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index compile_source_view_id(
    ExactVariantPlan *plan,
    const ExactRelationTemplate &relation_template);

inline const std::vector<semantic::Index> &planned_source_support_for_id(
    const ExactVariantPlan &plan,
    const semantic::Index source_id) {
  const auto leaf_count = plan.lowered.program.layout.n_leaves;
  if (source_id < leaf_count) {
    return plan.leaf_supports[static_cast<std::size_t>(source_id)];
  }
  return plan.pool_supports[static_cast<std::size_t>(source_id - leaf_count)];
}

inline bool source_id_is_pool(const ExactVariantPlan &plan,
                              const semantic::Index source_id) {
  const auto leaf_count = plan.lowered.program.layout.n_leaves;
  return source_id >= leaf_count && source_id < plan.source_count;
}

inline bool append_transition_guard(
    ExactTransitionGuardSet *guards,
    const ExactTransitionGuardKind kind,
    const semantic::Index subject_id) {
  if (guards == nullptr || subject_id == semantic::kInvalidIndex) {
    return false;
  }
  for (const auto &guard : guards->guards) {
    if (guard.kind == kind && guard.subject_id == subject_id) {
      return true;
    }
  }
  guards->guards.push_back(ExactTransitionGuard{kind, subject_id});
  return true;
}

inline bool transition_has_source_guard(
    const ExactTransitionGuardSet &guards,
    const ExactTransitionGuardKind kind,
    const semantic::Index source_id) {
  return std::find_if(
             guards.guards.begin(),
             guards.guards.end(),
             [&](const ExactTransitionGuard &guard) {
               return guard.kind == kind && guard.subject_id == source_id;
             }) != guards.guards.end();
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

inline bool append_ready_requirement(ExactSymbolicTransitionScenario *scenario,
                                     const ExactVariantPlan &plan,
                                     const semantic::Index expr_idx) {
  const auto &program = plan.lowered.program;
  if (expr_is_runtime_const(program, expr_idx)) {
    return true;
  }
  if (expr_is_simple_event(program, expr_idx)) {
    return append_transition_guard(
        &scenario->transition.readiness_time_expr.requirements,
        ExactTransitionGuardKind::SourceBefore,
        program.expr_source_ids[static_cast<std::size_t>(expr_idx)]);
  }
  return append_transition_guard(
      &scenario->transition.readiness_time_expr.requirements,
      ExactTransitionGuardKind::ExprBefore,
      expr_idx);
}

inline bool append_tail_requirement(ExactSymbolicTransitionScenario *scenario,
                                    const ExactVariantPlan &plan,
                                    const semantic::Index expr_idx) {
  const auto &program = plan.lowered.program;
  if (expr_is_runtime_const(program, expr_idx)) {
    return true;
  }
  if (expr_is_simple_event(program, expr_idx)) {
    return append_transition_guard(
        &scenario->transition.guards,
        ExactTransitionGuardKind::SourceAfter,
        program.expr_source_ids[static_cast<std::size_t>(expr_idx)]);
  }
  return append_transition_guard(
      &scenario->transition.guards,
      ExactTransitionGuardKind::ExprAfter,
      expr_idx);
}

inline bool scenario_source_happens_by_transition(
    const ExactSymbolicTransitionScenario &scenario,
    const semantic::Index source_id) {
  if (exact_symbolic_transition_release_source_id(scenario.transition) ==
      source_id) {
    return true;
  }
  return transition_has_source_guard(
      scenario.transition.readiness_time_expr.requirements,
      ExactTransitionGuardKind::SourceBefore,
      source_id);
}

inline bool scenario_source_known_before_active(
    const ExactSymbolicTransitionScenario &scenario,
    const semantic::Index source_id) {
  return transition_has_source_guard(
      scenario.transition.readiness_time_expr.requirements,
      ExactTransitionGuardKind::SourceBefore,
      source_id);
}

inline void append_source_order_fact(
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

inline void append_scenario_source_order_fact(
    ExactSymbolicTransitionScenario *scenario,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  append_source_order_fact(
      &scenario->transition.order_region.source_order_facts,
      before_source_id,
      after_source_id);
}

inline bool append_tail_order_requirement(ExactSymbolicTransitionScenario *scenario,
                                          const ExactVariantPlan &plan,
                                          const semantic::Index expr_idx,
                                          bool *handled) {
  *handled = false;
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
  const auto ref_source_id =
      program.expr_source_ids[static_cast<std::size_t>(ref)];
  const auto blocker_source_id =
      program.expr_source_ids[static_cast<std::size_t>(blocker)];
  if (!scenario_source_happens_by_transition(*scenario, ref_source_id) ||
      !scenario_source_happens_by_transition(*scenario, blocker_source_id)) {
    return true;
  }
  if (exact_symbolic_transition_release_source_id(scenario->transition) ==
          blocker_source_id &&
      scenario_source_known_before_active(*scenario, ref_source_id)) {
    return false;
  }
  for (const auto &fact :
       scenario->transition.order_region.source_order_facts) {
    if (fact.before_source_id == blocker_source_id &&
        fact.after_source_id == ref_source_id) {
      *handled = true;
      return true;
    }
  }
  append_scenario_source_order_fact(scenario, blocker_source_id, ref_source_id);
  *handled = true;
  return true;
}

inline bool scenario_source_known_after_active(
    const ExactSymbolicTransitionScenario &scenario,
    const semantic::Index source_id) {
  return transition_has_source_guard(
      scenario.transition.guards,
      ExactTransitionGuardKind::SourceAfter,
      source_id);
}

inline bool scenario_source_known_before_source(
    const ExactSymbolicTransitionScenario &scenario,
    const semantic::Index before_source_id,
    const semantic::Index after_source_id) {
  if (before_source_id == after_source_id) {
    return false;
  }
  const auto active_source_id =
      exact_symbolic_transition_release_source_id(scenario.transition);
  if (active_source_id == after_source_id &&
      scenario_source_known_before_active(scenario, before_source_id)) {
    return true;
  }
  if (active_source_id == before_source_id &&
      scenario_source_known_after_active(scenario, after_source_id)) {
    return true;
  }
  for (const auto &fact :
       scenario.transition.order_region.source_order_facts) {
    if (fact.before_source_id == before_source_id &&
        fact.after_source_id == after_source_id) {
      return true;
    }
  }
  return false;
}

inline bool expr_certainly_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionScenario &scenario,
    const semantic::Index expr_idx);

inline bool expr_children_any_not_completed_by_active(
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionScenario &scenario,
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
    const ExactSymbolicTransitionScenario &scenario,
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
    const ExactSymbolicTransitionScenario &scenario,
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
    return scenario_source_known_after_active(
        scenario,
        program.expr_source_ids[static_cast<std::size_t>(expr_idx)]);
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
    return scenario_source_known_before_source(
        scenario,
        program.expr_source_ids[static_cast<std::size_t>(blocker)],
        program.expr_source_ids[static_cast<std::size_t>(ref)]);
  }
  case semantic::ExprKind::Not:
    return false;
  }
  return false;
}

inline bool append_transition_relation(
    ExactSymbolicTransitionScenario *scenario,
    const semantic::Index source_id,
    const ExactRelation relation) {
  if (scenario == nullptr ||
      source_id == semantic::kInvalidIndex ||
      relation == ExactRelation::Unknown) {
    return false;
  }
  auto &relation_template = scenario->transition.relation_template;
  for (std::size_t i = 0; i < relation_template.source_ids.size(); ++i) {
    if (relation_template.source_ids[i] == source_id) {
      return relation_template.relations[i] == relation;
    }
  }
  relation_template.source_ids.push_back(source_id);
  relation_template.relations.push_back(relation);
  return true;
}

inline bool append_source_truth_constraints(
    const ExactVariantPlan &plan,
    const semantic::Index source_id,
    const ExactRelation relation,
    ExactSymbolicTransitionScenario *scenario) {
  if (!append_transition_relation(scenario, source_id, relation)) {
    return false;
  }
  const auto leaf_count = plan.lowered.program.layout.n_leaves;
  if (source_id >= leaf_count) {
    return true;
  }
  if (relation == ExactRelation::After) {
    return true;
  }
  const auto pos = static_cast<std::size_t>(source_id);
  const auto onset_kind = static_cast<semantic::OnsetKind>(
      plan.lowered.program.onset_kind[pos]);
  if (onset_kind == semantic::OnsetKind::Absolute) {
    return true;
  }
  return append_source_truth_constraints(
      plan,
      plan.lowered.program.onset_source_ids[pos],
      ExactRelation::Before,
      scenario);
}

inline void append_unique_source_id(std::vector<semantic::Index> *source_ids,
                                    const semantic::Index source_id) {
  if (source_id == semantic::kInvalidIndex ||
      std::find(source_ids->begin(), source_ids->end(), source_id) !=
          source_ids->end()) {
    return;
  }
  source_ids->push_back(source_id);
}

inline bool scenario_sources_supported(
    const ExactVariantPlan &plan,
    const ExactSymbolicTransitionScenario &scenario) {
  std::vector<semantic::Index> referenced_sources;
  referenced_sources.reserve(
      1U +
      scenario.transition.readiness_time_expr.requirements.guards.size() +
      scenario.transition.guards.guards.size());
  append_unique_source_id(
      &referenced_sources,
      exact_symbolic_transition_release_source_id(scenario.transition));
  for (const auto &guard :
       scenario.transition.readiness_time_expr.requirements.guards) {
    if (guard.kind == ExactTransitionGuardKind::SourceBefore) {
      append_unique_source_id(&referenced_sources, guard.subject_id);
    }
  }
  for (const auto &guard : scenario.transition.guards.guards) {
    if (guard.kind == ExactTransitionGuardKind::SourceAfter) {
      append_unique_source_id(&referenced_sources, guard.subject_id);
    }
  }

  for (std::size_t i = 0; i < referenced_sources.size(); ++i) {
    const auto lhs = referenced_sources[i];
    const auto &lhs_support = planned_source_support_for_id(plan, lhs);
    for (std::size_t j = i + 1U; j < referenced_sources.size(); ++j) {
      const auto rhs = referenced_sources[j];
      if (supports_overlap(lhs_support, planned_source_support_for_id(plan, rhs))) {
        return false;
      }
    }
  }

  const auto &relation_template = scenario.transition.relation_template;
  for (std::size_t constraint_idx = 0;
       constraint_idx < relation_template.source_ids.size();
       ++constraint_idx) {
    const auto source_id = relation_template.source_ids[constraint_idx];
    if (!source_id_is_pool(plan, source_id)) {
      continue;
    }
    const auto &forced_support = planned_source_support_for_id(plan, source_id);
    for (std::size_t other_idx = constraint_idx + 1U;
         other_idx < relation_template.source_ids.size();
         ++other_idx) {
      const auto other_source_id = relation_template.source_ids[other_idx];
      if (other_source_id == source_id || !source_id_is_pool(plan, other_source_id)) {
        continue;
      }
      if (supports_overlap(
              forced_support,
              planned_source_support_for_id(plan, other_source_id))) {
        return false;
      }
    }
    for (const auto referenced_source_id : referenced_sources) {
      if (referenced_source_id == source_id) {
        continue;
      }
      if (supports_overlap(
              forced_support,
              planned_source_support_for_id(plan, referenced_source_id))) {
        return false;
      }
    }
    for (const auto &guard :
         scenario.transition.readiness_time_expr.requirements.guards) {
      if (guard.kind == ExactTransitionGuardKind::ExprBefore &&
          supports_overlap(
              forced_support,
              plan.expr_supports[static_cast<std::size_t>(guard.subject_id)])) {
          return false;
      }
    }
    for (const auto &guard : scenario.transition.guards.guards) {
      if (guard.kind == ExactTransitionGuardKind::ExprAfter &&
          supports_overlap(
              forced_support,
              plan.expr_supports[static_cast<std::size_t>(guard.subject_id)])) {
          return false;
      }
    }
  }
  return true;
}

inline bool expand_member_subsets(
    const std::vector<semantic::Index> &members,
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

inline ExactSymbolicTransitionScenario make_source_release_transition_scenario(
    const semantic::Index source_id) {
  ExactSymbolicTransitionScenario scenario;
  scenario.transition.transition_time_expr =
      ExactSymbolicTransitionTimeExpr{
          ExactSymbolicTransitionTimeKind::SourceRelease,
          source_id};
  scenario.transition.active_sources.push_back(source_id);
  scenario.transition.readiness_time_expr.requirements.empty_value = 1.0;
  scenario.transition.guards.empty_value = 1.0;
  return scenario;
}

inline std::vector<ExactSymbolicTransitionScenario> build_source_transition_scenarios(
    const ExactVariantPlan &plan,
    const semantic::Index source_id) {
  std::vector<ExactSymbolicTransitionScenario> out;
  const auto leaf_count = plan.lowered.program.layout.n_leaves;
  if (source_id < leaf_count) {
    auto scenario = make_source_release_transition_scenario(source_id);
    if (!append_source_truth_constraints(
            plan, source_id, ExactRelation::At, &scenario)) {
      return out;
    }
    if (scenario_sources_supported(plan, scenario)) {
      out.push_back(std::move(scenario));
    }
    return out;
  }

  const auto pool_idx = static_cast<std::size_t>(source_id - leaf_count);
  const auto begin = plan.lowered.program.pool_member_offsets[pool_idx];
  const auto end = plan.lowered.program.pool_member_offsets[pool_idx + 1U];
  const auto k = plan.lowered.program.pool_k[pool_idx];
  std::vector<semantic::Index> members;
  members.reserve(static_cast<std::size_t>(end - begin));
  for (semantic::Index i = begin; i < end; ++i) {
    members.push_back(
        plan.lowered.program.pool_member_source_ids[static_cast<std::size_t>(i)]);
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
      auto scenario = make_source_release_transition_scenario(
          members[static_cast<std::size_t>(active_idx)]);
      if (!append_source_truth_constraints(
              plan,
              members[static_cast<std::size_t>(active_idx)],
              ExactRelation::At,
              &scenario)) {
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
        const auto member_source_id = members[static_cast<std::size_t>(member_idx)];
        ok = append_transition_guard(
                 is_before
                     ? &scenario.transition.readiness_time_expr.requirements
                     : &scenario.transition.guards,
                 is_before ? ExactTransitionGuardKind::SourceBefore
                           : ExactTransitionGuardKind::SourceAfter,
                 member_source_id) &&
             append_source_truth_constraints(
                 plan,
                 member_source_id,
                 relation,
                 &scenario);
        if (!ok) {
          ok = false;
          break;
        }
      }
      if (!ok) {
        continue;
      }
      if (scenario_sources_supported(plan, scenario)) {
        out.push_back(std::move(scenario));
      }
    }
  }
  return out;
}

inline std::vector<ExactSymbolicTransitionScenario> build_expr_transition_scenarios(
    const ExactVariantPlan &plan,
    const semantic::Index expr_idx);

inline std::vector<ExactSymbolicTransitionScenario> build_logical_transition_scenarios(
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
  std::vector<std::vector<ExactSymbolicTransitionScenario>> child_scenarios(
      normalized_children.size());
  for (std::size_t i = 0; i < normalized_children.size(); ++i) {
    child_scenarios[i] =
        build_expr_transition_scenarios(plan, normalized_children[i]);
    if (kind == semantic::ExprKind::Or && child_scenarios[i].empty()) {
      throw std::runtime_error(
          "exact kernel does not support logical-not branches inside first_of()/or outcomes");
    }
  }
  std::vector<ExactSymbolicTransitionScenario> out;
  for (std::size_t active = 0; active < normalized_children.size(); ++active) {
    if (child_scenarios[active].empty()) {
      continue;
    }
    for (const auto &child_scenario : child_scenarios[active]) {
      auto scenario = child_scenario;
      bool ok = true;
      for (std::size_t other = 0; other < normalized_children.size(); ++other) {
        if (other == active) {
          continue;
        }
        const auto child = normalized_children[other];
        if (kind == semantic::ExprKind::And) {
          ok = append_ready_requirement(&scenario, plan, child);
          if (ok && expr_is_simple_event(plan.lowered.program, child)) {
            ok = append_source_truth_constraints(
                plan,
                plan.lowered.program.expr_source_ids[
                    static_cast<std::size_t>(child)],
                ExactRelation::Before,
                &scenario);
          }
        } else {
          if (expr_certainly_not_completed_by_active(
                  plan,
                  scenario,
                  child)) {
            continue;
          }
          bool tail_order_handled = false;
          ok = append_tail_order_requirement(
              &scenario, plan, child, &tail_order_handled);
          if (ok && !tail_order_handled) {
            ok = append_tail_requirement(&scenario, plan, child);
          }
          if (ok && expr_is_simple_event(plan.lowered.program, child)) {
            ok = append_source_truth_constraints(
                plan,
                plan.lowered.program.expr_source_ids[
                    static_cast<std::size_t>(child)],
                ExactRelation::After,
                &scenario);
          }
        }
        if (!ok) {
          break;
        }
      }
      if (!ok) {
        continue;
      }
      if (scenario_sources_supported(plan, scenario)) {
        out.push_back(std::move(scenario));
      }
    }
  }
  return out;
}

inline std::vector<ExactSymbolicTransitionScenario> build_expr_conjunction_transition_scenarios(
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

inline std::vector<ExactSymbolicTransitionScenario> build_expr_transition_scenarios(
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
        program.expr_source_ids[static_cast<std::size_t>(expr_idx)]);
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
    std::vector<ExactSymbolicTransitionScenario> out;

    for (const auto &ref_scenario : ref_scenarios) {
      if (!blocker_true) {
        auto scenario = ref_scenario;
        bool ok = append_tail_requirement(&scenario, plan, blocker);
        if (ok && expr_is_simple_event(program, blocker)) {
          ok = append_source_truth_constraints(
              plan,
              program.expr_source_ids[static_cast<std::size_t>(blocker)],
              ExactRelation::After,
              &scenario);
        }
        for (const auto child : unless_events) {
          ok = ok && append_tail_requirement(&scenario, plan, child);
          if (ok && expr_is_simple_event(program, child)) {
            ok = append_source_truth_constraints(
                plan,
                program.expr_source_ids[static_cast<std::size_t>(child)],
                ExactRelation::After,
                &scenario);
          }
        }
        if (ok && scenario_sources_supported(plan, scenario)) {
          out.push_back(std::move(scenario));
        }
      }

      const int n_unless = static_cast<int>(unless_events.size());
      if (n_unless == 0) {
        continue;
      }
      const int max_mask = 1 << n_unless;
      for (int mask = 1; mask < max_mask; ++mask) {
        auto scenario = ref_scenario;
        bool ok = true;
        for (int bit = 0; bit < n_unless; ++bit) {
          const auto child = unless_events[static_cast<std::size_t>(bit)];
          const bool active = (mask & (1 << bit)) != 0;
          ok = active ? append_ready_requirement(&scenario, plan, child)
                      : append_tail_requirement(&scenario, plan, child);
          if (ok && expr_is_simple_event(program, child)) {
            ok = append_source_truth_constraints(
                plan,
                program.expr_source_ids[static_cast<std::size_t>(child)],
                active ? ExactRelation::Before : ExactRelation::After,
                &scenario);
          }
          if (!ok) {
            break;
          }
        }
        if (ok && scenario_sources_supported(plan, scenario)) {
          out.push_back(std::move(scenario));
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

struct ExactCompetitorSubsetPlan {
  std::vector<semantic::Index> outcome_indices;
  std::vector<semantic::Index> expr_roots;
  semantic::Index singleton_expr_root{semantic::kInvalidIndex};
  int inclusion_sign{1};
  std::vector<ExactSymbolicTransitionScenario> scenarios;
};

struct ExactCompetitorBlockPlan {
  std::vector<ExactCompetitorSubsetPlan> subsets;
};

struct ExactTargetCompetitorPlan {
  std::vector<ExactCompetitorBlockPlan> blocks;
};

inline void canonicalize_transition_relation_template(
    ExactRelationTemplate *relation_template) {
  std::vector<std::pair<semantic::Index, ExactRelation>> relation_pairs;
  relation_pairs.reserve(relation_template->source_ids.size());
  for (std::size_t i = 0; i < relation_template->source_ids.size(); ++i) {
    const auto source_id = relation_template->source_ids[i];
    if (source_id == semantic::kInvalidIndex) {
      continue;
    }
    relation_pairs.emplace_back(source_id, relation_template->relations[i]);
  }
  std::sort(
      relation_pairs.begin(),
      relation_pairs.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
  relation_template->source_ids.clear();
  relation_template->relations.clear();
  relation_template->source_ids.reserve(relation_pairs.size());
  relation_template->relations.reserve(relation_pairs.size());
  for (const auto &[source_id, relation] : relation_pairs) {
    if (!relation_template->source_ids.empty() &&
        relation_template->source_ids.back() == source_id) {
      if (relation_template->relations.back() != relation) {
        throw std::runtime_error(
            "symbolic transition relation template has conflicting source relation");
      }
      continue;
    }
    relation_template->source_ids.push_back(source_id);
    relation_template->relations.push_back(relation);
  }
}

inline void finalize_symbolic_transition_scenario(
    ExactVariantPlan *plan,
    ExactSymbolicTransitionScenario *scenario,
    const semantic::Index visible_outcome) {
  scenario->visible_outcome = visible_outcome;
  canonicalize_transition_relation_template(
      &scenario->transition.relation_template);
  scenario->transition.source_view_id =
      compile_source_view_id(plan, scenario->transition.relation_template);
}

inline std::vector<ExactSymbolicTransitionScenario>
finalize_symbolic_transition_scenarios(
    ExactVariantPlan *plan,
    std::vector<ExactSymbolicTransitionScenario> scenarios,
    const semantic::Index visible_outcome = semantic::kInvalidIndex) {
  for (auto &scenario : scenarios) {
    finalize_symbolic_transition_scenario(plan, &scenario, visible_outcome);
  }
  return scenarios;
}

inline void enumerate_competitor_subsets(
    ExactVariantPlan *plan,
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
              static_cast<semantic::Index>(plan->outcomes.size()) &&
          plan->outcomes[static_cast<std::size_t>(competitor.outcome_index)]
                  .expr_root == competitor.expr_root) {
        subset.scenarios =
            plan->outcomes[static_cast<std::size_t>(competitor.outcome_index)]
                .scenarios;
      } else {
        subset.scenarios =
            finalize_symbolic_transition_scenarios(
                plan,
                build_expr_transition_scenarios(*plan, competitor.expr_root),
                competitor.outcome_index);
      }
    } else {
      subset.scenarios =
          finalize_symbolic_transition_scenarios(
              plan,
              build_expr_conjunction_transition_scenarios(
                  *plan, subset.expr_roots));
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
    ExactVariantPlan *plan,
    const semantic::Index target_outcome_idx) {
  ExactTargetCompetitorPlan target_plan;
  const auto blocks = build_competitor_overlap_blocks(*plan, target_outcome_idx);
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

inline std::vector<ExactTargetCompetitorPlan> compile_target_competitor_plans(
    ExactVariantPlan *plan) {
  std::vector<ExactTargetCompetitorPlan> competitor_plans;
  competitor_plans.reserve(plan->outcomes.size());
  for (semantic::Index target_outcome_idx = 0;
       target_outcome_idx < static_cast<semantic::Index>(plan->outcomes.size());
       ++target_outcome_idx) {
    competitor_plans.push_back(
        build_target_competitor_plan(plan, target_outcome_idx));
  }
  return competitor_plans;
}

inline void compile_exact_outcome_transition_scenarios(
    ExactVariantPlan *plan,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_outcome_codes) {
  plan->outcome_index_by_code.assign(
      n_outcome_codes + 1U,
      semantic::kInvalidIndex);
  const auto &program = plan->lowered.program;
  plan->outcomes.clear();
  plan->outcomes.reserve(plan->lowered.outcome_labels.size());
  for (std::size_t i = 0; i < plan->lowered.outcome_labels.size(); ++i) {
    const auto expr_root = program.outcome_expr_root[i];
    validate_exact_expr(plan->lowered, expr_root);
    ExactOutcomePlan outcome;
    outcome.expr_root = expr_root;
    outcome.scenarios = finalize_symbolic_transition_scenarios(
        plan,
        build_expr_transition_scenarios(*plan, expr_root),
        static_cast<semantic::Index>(i));
    const auto code_it =
        outcome_code_by_label.find(plan->lowered.outcome_labels[i]);
    if (code_it == outcome_code_by_label.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared outcome code for '" +
          plan->lowered.outcome_labels[i] + "'");
    }
    plan->outcome_index_by_code[static_cast<std::size_t>(code_it->second)] =
        static_cast<semantic::Index>(i);
    plan->outcomes.push_back(std::move(outcome));
  }
}


} // namespace detail
} // namespace accumulatr::eval
