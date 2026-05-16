#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "../runtime/exact_evaluation_program.hpp"
#include "project_semantic.hpp"

namespace accumulatr::compile {

inline runtime::ExactEvaluationProgram lower_exact_evaluation_program(
    const CompiledVariant &variant,
    const std::unordered_map<std::string, semantic::Index>
        &outcome_code_by_label) {
  const auto &model = variant.model;
  runtime::ExactEvaluationProgram program;

  program.layout.n_leaves = static_cast<int>(model.leaves.size());
  program.layout.n_pools = static_cast<int>(model.pools.size());
  program.layout.n_outcomes = static_cast<int>(model.outcomes.size());
  program.layout.n_triggers = static_cast<int>(model.triggers.size());

  program.leaf_dist_kind.reserve(model.leaves.size());
  program.leaf_descriptors.reserve(model.leaves.size());
  program.onset_kind.reserve(model.leaves.size());
  program.onset_source_kind.reserve(model.leaves.size());
  program.onset_source_index.reserve(model.leaves.size());
  program.onset_source_ids.reserve(model.leaves.size());
  program.onset_lag.reserve(model.leaves.size());
  program.onset_abs_value.reserve(model.leaves.size());
  program.leaf_trigger_index.reserve(model.leaves.size());

  for (const auto &leaf : model.leaves) {
    program.leaf_dist_kind.push_back(static_cast<std::uint8_t>(leaf.dist));
    program.onset_kind.push_back(static_cast<std::uint8_t>(leaf.onset.kind));
    program.onset_source_kind.push_back(
        static_cast<std::uint8_t>(leaf.onset.source.kind));
    program.onset_source_index.push_back(leaf.onset.source.index);
    program.onset_source_ids.push_back(semantic::kInvalidIndex);
    program.onset_lag.push_back(leaf.onset.lag);
    program.onset_abs_value.push_back(leaf.onset.absolute_value);
    program.leaf_trigger_index.push_back(leaf.trigger_index);

    program.leaf_descriptors.push_back(runtime::LeafRuntimeDescriptor{
        static_cast<std::uint8_t>(leaf.dist),
        static_cast<std::uint8_t>(leaf.onset.kind),
        static_cast<std::uint8_t>(leaf.onset.source.kind),
        leaf.onset.source.index,
        semantic::kInvalidIndex,
        leaf.onset.lag,
        leaf.onset.absolute_value,
        leaf.trigger_index,
        static_cast<int>(leaf.params.dist_param_names.size())});
  }

  program.trigger_member_offsets.reserve(model.triggers.size() + 1U);
  program.trigger_member_offsets.push_back(0);
  for (const auto &trigger : model.triggers) {
    for (const auto leaf_index : trigger.leaf_indices) {
      program.trigger_member_indices.push_back(leaf_index);
    }
    program.trigger_member_offsets.push_back(
        static_cast<semantic::Index>(program.trigger_member_indices.size()));
  }

  program.pool_k.reserve(model.pools.size());
  program.pool_member_offsets.reserve(model.pools.size() + 1U);
  program.pool_member_source_ids.reserve(model.pools.size());
  program.pool_member_offsets.push_back(0);
  for (const auto &pool : model.pools) {
    program.pool_k.push_back(pool.k);
    for (const auto &member : pool.members) {
      program.pool_member_indices.push_back(member.index);
      program.pool_member_kind.push_back(
          static_cast<std::uint8_t>(member.kind));
      program.pool_member_source_ids.push_back(semantic::kInvalidIndex);
    }
    program.pool_member_offsets.push_back(
        static_cast<semantic::Index>(program.pool_member_indices.size()));
  }

  program.expr_kind.reserve(model.expr_nodes.size());
  program.expr_arg_offsets.reserve(model.expr_nodes.size() + 1U);
  program.expr_ref_child.reserve(model.expr_nodes.size());
  program.expr_blocker_child.reserve(model.expr_nodes.size());
  program.expr_source_index.reserve(model.expr_nodes.size());
  program.expr_source_kind.reserve(model.expr_nodes.size());
  program.expr_source_ids.reserve(model.expr_nodes.size());
  program.expr_event_k.reserve(model.expr_nodes.size());
  program.expr_arg_offsets.push_back(0);
  for (const auto &expr : model.expr_nodes) {
    program.expr_kind.push_back(static_cast<std::uint8_t>(expr.kind));
    for (const auto child : expr.children) {
      program.expr_args.push_back(child);
    }
    for (const auto child : expr.unless_children) {
      program.expr_args.push_back(child);
    }
    program.expr_arg_offsets.push_back(
        static_cast<semantic::Index>(program.expr_args.size()));
    program.expr_ref_child.push_back(expr.reference_child);
    program.expr_blocker_child.push_back(expr.blocker_child);
    program.expr_source_index.push_back(expr.source.index);
    program.expr_source_kind.push_back(
        static_cast<std::uint8_t>(expr.source.kind));
    program.expr_source_ids.push_back(semantic::kInvalidIndex);
    program.expr_event_k.push_back(expr.event_k);
  }

  program.outcome_expr_root.reserve(model.outcomes.size());
  program.outcome_codes.reserve(model.outcomes.size());
  program.outcome_competitor_offsets.reserve(model.outcomes.size() + 1U);
  program.outcome_competitor_offsets.push_back(0);
  for (const auto &outcome : model.outcomes) {
    const auto code_it = outcome_code_by_label.find(outcome.label);
    if (code_it == outcome_code_by_label.end()) {
      throw std::runtime_error(
          "exact evaluator found no prepared outcome code for '" +
          outcome.label + "'");
    }
    program.outcome_expr_root.push_back(outcome.expr_root);
    program.outcome_codes.push_back(code_it->second);
    for (std::size_t j = 0; j < outcome.competitor_expr_roots.size(); ++j) {
      program.outcome_competitor_expr_roots.push_back(
          outcome.competitor_expr_roots[j]);
      program.outcome_competitor_indices.push_back(
          j < outcome.competitor_outcome_indices.size()
              ? outcome.competitor_outcome_indices[j]
              : semantic::kInvalidIndex);
    }
    program.outcome_competitor_offsets.push_back(
        static_cast<semantic::Index>(
            program.outcome_competitor_expr_roots.size()));
  }

  return program;
}

} // namespace accumulatr::compile
