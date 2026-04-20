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
          if (expr_is_simple_event(plan.lowered.program, child)) {
            const auto child_key =
                ExactSourceKey{child_event_source_kind(plan.lowered.program, child),
                               child_event_source_index(plan.lowered.program, child)};
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

} // namespace detail
} // namespace accumulatr::eval
