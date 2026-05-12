#pragma once

#include <algorithm>
#include <initializer_list>
#include <vector>

#include "compiled_math_kernel_planning.hpp"
#include "compiled_math_workspace.hpp"

namespace accumulatr::eval {
namespace detail {

inline semantic::Index compiled_math_intern_node(
    CompiledMathProgram *program,
    CompiledMathNodeKey key) {
  const auto found = program->node_index.find(key);
  if (found != program->node_index.end()) {
    return found->second;
  }

  const auto node_id = static_cast<semantic::Index>(program->nodes.size());
  const auto child_offset =
      static_cast<semantic::Index>(program->child_nodes.size());
  program->child_nodes.insert(
      program->child_nodes.end(), key.children.begin(), key.children.end());

  CompiledMathNode node;
  node.kind = key.kind;
  node.value_kind = key.value_kind;
  node.subject_id = key.subject_id;
  node.condition_id = key.condition_id;
  node.time_id = key.time_id;
  node.aux_id = key.aux_id;
  node.aux2_id = key.aux2_id;
  node.source_view_id = key.source_view_id;
  node.children = CompiledMathIndexSpan{
      child_offset,
      static_cast<semantic::Index>(key.children.size())};
  node.cache_slot = node_id;
  node.integral_kernel_slot = compiled_math_integral_kernel_slot(program, node);
  node.constant = key.constant;
  program->nodes.push_back(node);
  program->node_index.emplace(std::move(key), node_id);
  return node_id;
}

inline void compiled_math_append_fact_index(
    std::vector<std::vector<semantic::Index>> *by_key,
    const semantic::Index key,
    const semantic::Index fact_index) {
  if (key == semantic::kInvalidIndex) {
    return;
  }
  const auto pos = static_cast<std::size_t>(key);
  if (by_key->size() <= pos) {
    by_key->resize(pos + 1U);
  }
  (*by_key)[pos].push_back(fact_index);
}

inline void compiled_math_finish_fact_spans(
    const std::vector<std::vector<semantic::Index>> &by_key,
    std::vector<CompiledMathIndexSpan> *spans,
    std::vector<semantic::Index> *indices) {
  spans->assign(by_key.size(), CompiledMathIndexSpan{});
  indices->clear();
  for (std::size_t i = 0; i < by_key.size(); ++i) {
    const auto offset = static_cast<semantic::Index>(indices->size());
    indices->insert(indices->end(), by_key[i].begin(), by_key[i].end());
    (*spans)[i] = CompiledMathIndexSpan{
        offset,
        static_cast<semantic::Index>(by_key[i].size())};
  }
}

inline void compiled_math_append_pair_fact_index(
    std::vector<CompiledMathPairFactEntry> *entries,
    std::vector<semantic::Index> *indices,
    const semantic::Index first,
    const semantic::Index second,
    const semantic::Index fact_index) {
  if (first == semantic::kInvalidIndex ||
      second == semantic::kInvalidIndex) {
    return;
  }
  entries->push_back(CompiledMathPairFactEntry{
      first,
      second,
      CompiledMathIndexSpan{
          static_cast<semantic::Index>(indices->size()),
          1}});
  indices->push_back(fact_index);
}

inline void compiled_math_finish_pair_lookup(
    CompiledMathPairFactLookup *lookup) {
  if (lookup->entries.empty()) {
    return;
  }
  std::vector<CompiledMathPairFactEntry> old_entries =
      std::move(lookup->entries);
  std::vector<semantic::Index> old_indices = std::move(lookup->fact_indices);
  std::vector<std::size_t> order(old_entries.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  std::sort(
      order.begin(),
      order.end(),
      [&](const auto lhs, const auto rhs) {
        const auto &a = old_entries[lhs];
        const auto &b = old_entries[rhs];
        return a.first < b.first ||
               (a.first == b.first && a.second < b.second);
      });
  lookup->entries.clear();
  lookup->fact_indices.clear();
  std::size_t order_pos = 0;
  while (order_pos < order.size()) {
    const auto first = old_entries[order[order_pos]].first;
    const auto second = old_entries[order[order_pos]].second;
    const auto offset =
        static_cast<semantic::Index>(lookup->fact_indices.size());
    while (order_pos < order.size() &&
           old_entries[order[order_pos]].first == first &&
           old_entries[order[order_pos]].second == second) {
      const auto old_offset = static_cast<std::size_t>(
          old_entries[order[order_pos]].facts.offset);
      lookup->fact_indices.push_back(old_indices[old_offset]);
      ++order_pos;
    }
    lookup->entries.push_back(CompiledMathPairFactEntry{
        first,
        second,
        CompiledMathIndexSpan{
            offset,
            static_cast<semantic::Index>(
                lookup->fact_indices.size() -
                static_cast<std::size_t>(offset))}});
  }
}

inline void compiled_math_build_condition_access(
    CompiledMathCondition *condition) {
  std::vector<std::vector<semantic::Index>> source_exact;
  std::vector<std::vector<semantic::Index>> source_lower;
  std::vector<std::vector<semantic::Index>> source_upper;
  std::vector<std::vector<semantic::Index>> expr_upper;
  condition->time_dependency_ids.clear();
  condition->source_order_fact_lookup.entries.clear();
  condition->source_order_fact_lookup.fact_indices.clear();

  for (std::size_t i = 0; i < condition->fact_kinds.size(); ++i) {
    const auto time_id = condition->fact_time_ids[i];
    if (time_id != semantic::kInvalidIndex &&
        std::find(
            condition->time_dependency_ids.begin(),
            condition->time_dependency_ids.end(),
            time_id) == condition->time_dependency_ids.end()) {
      condition->time_dependency_ids.push_back(time_id);
    }
    const auto kind =
        static_cast<CompiledMathConditionFactKind>(condition->fact_kinds[i]);
    const auto fact_index = static_cast<semantic::Index>(i);
    switch (kind) {
    case CompiledMathConditionFactKind::SourceExact:
      compiled_math_append_fact_index(
          &source_exact, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceLowerBound:
      compiled_math_append_fact_index(
          &source_lower, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceUpperBound:
      compiled_math_append_fact_index(
          &source_upper, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::ExprUpperBound:
      compiled_math_append_fact_index(
          &expr_upper, condition->fact_subject_ids[i], fact_index);
      break;
    case CompiledMathConditionFactKind::SourceOrder:
      compiled_math_append_pair_fact_index(
          &condition->source_order_fact_lookup.entries,
          &condition->source_order_fact_lookup.fact_indices,
          condition->fact_subject_ids[i],
          condition->fact_aux_ids[i],
          fact_index);
      break;
    }
  }

  compiled_math_finish_fact_spans(
      source_exact,
      &condition->source_exact_fact_spans,
      &condition->source_exact_fact_indices);
  compiled_math_finish_fact_spans(
      source_lower,
      &condition->source_lower_fact_spans,
      &condition->source_lower_fact_indices);
  compiled_math_finish_fact_spans(
      source_upper,
      &condition->source_upper_fact_spans,
      &condition->source_upper_fact_indices);
  compiled_math_finish_fact_spans(
      expr_upper,
      &condition->expr_upper_fact_spans,
      &condition->expr_upper_fact_indices);
  compiled_math_finish_pair_lookup(&condition->source_order_fact_lookup);
}

inline semantic::Index compiled_math_constant(CompiledMathProgram *program,
                                             const double value) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::Constant;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.constant = value;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_intern_condition(
    CompiledMathProgram *program,
    CompiledMathConditionKey key) {
  if (key.source_ids.empty() && key.fact_kinds.empty()) {
    return 0;
  }
  const auto found = program->condition_index.find(key);
  if (found != program->condition_index.end()) {
    return found->second;
  }
  const auto condition_id =
      static_cast<semantic::Index>(program->conditions.size() + 1U);
  program->conditions.push_back(
      CompiledMathCondition{
          key.impossible,
          key.source_ids,
          key.relations,
          key.fact_kinds,
          key.fact_subject_ids,
          key.fact_aux_ids,
          key.fact_aux2_ids,
          key.fact_time_ids,
          key.fact_normalizer_node_ids,
          key.fact_normalizer_root_ids});
  compiled_math_build_condition_access(&program->conditions.back());
  program->condition_index.emplace(std::move(key), condition_id);
  return condition_id;
}

inline bool compiled_condition_impossible(
    const CompiledMathProgram &program,
    const semantic::Index condition_id) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return false;
  }
  const auto pos = static_cast<std::size_t>(condition_id - 1U);
  return pos >= program.conditions.size() || program.conditions[pos].impossible;
}

inline void compiled_math_append_condition_to_key(
    const CompiledMathProgram &program,
    const semantic::Index condition_id,
    CompiledMathConditionKey *key) {
  if (condition_id == 0 || condition_id == semantic::kInvalidIndex) {
    return;
  }
  const auto pos = static_cast<std::size_t>(condition_id - 1U);
  if (pos >= program.conditions.size()) {
    key->impossible = true;
    return;
  }
  const auto &condition = program.conditions[pos];
  key->impossible = key->impossible || condition.impossible;
  key->source_ids.insert(
      key->source_ids.end(),
      condition.source_ids.begin(),
      condition.source_ids.end());
  key->relations.insert(
      key->relations.end(),
      condition.relations.begin(),
      condition.relations.end());
  key->fact_kinds.insert(
      key->fact_kinds.end(),
      condition.fact_kinds.begin(),
      condition.fact_kinds.end());
  key->fact_subject_ids.insert(
      key->fact_subject_ids.end(),
      condition.fact_subject_ids.begin(),
      condition.fact_subject_ids.end());
  key->fact_aux_ids.insert(
      key->fact_aux_ids.end(),
      condition.fact_aux_ids.begin(),
      condition.fact_aux_ids.end());
  key->fact_aux2_ids.insert(
      key->fact_aux2_ids.end(),
      condition.fact_aux2_ids.begin(),
      condition.fact_aux2_ids.end());
  key->fact_time_ids.insert(
      key->fact_time_ids.end(),
      condition.fact_time_ids.begin(),
      condition.fact_time_ids.end());
  key->fact_normalizer_node_ids.insert(
      key->fact_normalizer_node_ids.end(),
      condition.fact_normalizer_node_ids.begin(),
      condition.fact_normalizer_node_ids.end());
  key->fact_normalizer_root_ids.insert(
      key->fact_normalizer_root_ids.end(),
      condition.fact_normalizer_root_ids.begin(),
      condition.fact_normalizer_root_ids.end());
}

inline semantic::Index compiled_math_merge_conditions(
    CompiledMathProgram *program,
    const std::initializer_list<semantic::Index> condition_ids) {
  CompiledMathConditionKey key;
  for (const auto condition_id : condition_ids) {
    compiled_math_append_condition_to_key(*program, condition_id, &key);
  }
  return compiled_math_intern_condition(program, std::move(key));
}

inline semantic::Index compiled_math_source_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index time_cap_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.subject_id = source_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux_id = time_cap_id;
  key.source_view_id = source_view_id;
  if (kind == CompiledMathNodeKind::SourcePdf) {
    key.value_kind = CompiledMathValueKind::Pdf;
  } else if (kind == CompiledMathNodeKind::SourceCdf) {
    key.value_kind = CompiledMathValueKind::Cdf;
  } else {
    key.value_kind = CompiledMathValueKind::Survival;
  }
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_algebra_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    std::vector<semantic::Index> children,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  if (children.empty()) {
    if (kind == CompiledMathNodeKind::Product) {
      return compiled_math_constant(program, 1.0);
    }
    return compiled_math_constant(program, 0.0);
  }
  if (children.size() == 1U &&
      (kind == CompiledMathNodeKind::Product ||
       kind == CompiledMathNodeKind::Sum ||
       kind == CompiledMathNodeKind::CleanSignedSum)) {
    return children.front();
  }
  if (kind == CompiledMathNodeKind::Product ||
      kind == CompiledMathNodeKind::Sum) {
    std::sort(children.begin(), children.end());
  }
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = value_kind;
  key.children = std::move(children);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_unary_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index child,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.value_kind = value_kind;
  key.children.push_back(child);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_time_gate_node(
    CompiledMathProgram *program,
    const semantic::Index child,
    const semantic::Index current_time_id,
    const semantic::Index gate_time_id,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::TimeGate;
  key.value_kind = value_kind;
  key.time_id = current_time_id;
  key.aux_id = gate_time_id;
  key.children.push_back(child);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_strict_time_gate_node(
    CompiledMathProgram *program,
    const semantic::Index child,
    const semantic::Index current_time_id,
    const semantic::Index gate_time_id,
    const CompiledMathValueKind value_kind = CompiledMathValueKind::Scalar) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::StrictTimeGate;
  key.value_kind = value_kind;
  key.time_id = current_time_id;
  key.aux_id = gate_time_id;
  key.children.push_back(child);
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_integral_zero_to_current_node(
    CompiledMathProgram *program,
    const semantic::Index integrand_root_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::IntegralZeroToCurrent;
  key.value_kind = CompiledMathValueKind::Cdf;
  key.subject_id = integrand_root_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  return compiled_math_intern_node(program, std::move(key));
}

inline semantic::Index compiled_math_raw_integral_zero_to_current_node(
    CompiledMathProgram *program,
    const semantic::Index integrand_root_id,
    const semantic::Index condition_id = 0,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = CompiledMathNodeKind::IntegralZeroToCurrentRaw;
  key.value_kind = CompiledMathValueKind::Scalar;
  key.subject_id = integrand_root_id;
  key.condition_id = condition_id;
  key.time_id = time_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  return compiled_math_intern_node(program, std::move(key));
}

inline void compiled_math_append_schedule_node(
    const CompiledMathProgram &program,
    const semantic::Index node_id,
    std::vector<std::uint8_t> *visited,
    std::vector<semantic::Index> *schedule) {
  const auto pos = static_cast<std::size_t>(node_id);
  if ((*visited)[pos] != 0U) {
    return;
  }
  (*visited)[pos] = 1U;
  const auto &node = program.nodes[pos];
  if (node.condition_id != 0 &&
      node.condition_id != semantic::kInvalidIndex) {
    const auto condition_pos = static_cast<std::size_t>(node.condition_id - 1U);
    if (condition_pos < program.conditions.size()) {
      const auto &condition = program.conditions[condition_pos];
      for (const auto normalizer_node_id :
           condition.fact_normalizer_node_ids) {
        if (normalizer_node_id != semantic::kInvalidIndex) {
          compiled_math_append_schedule_node(
              program, normalizer_node_id, visited, schedule);
        }
      }
    }
  }
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_id = program.child_nodes[
        static_cast<std::size_t>(node.children.offset + i)];
    compiled_math_append_schedule_node(program, child_id, visited, schedule);
  }
  schedule->push_back(node_id);
}

inline semantic::Index compiled_math_make_root(CompiledMathProgram *program,
                                              const semantic::Index node_id) {
  for (semantic::Index root_id = 0;
       root_id < static_cast<semantic::Index>(program->roots.size());
       ++root_id) {
    if (program->roots[static_cast<std::size_t>(root_id)].node_id == node_id) {
      return root_id;
    }
  }
  const auto root_id = static_cast<semantic::Index>(program->roots.size());
  const auto offset =
      static_cast<semantic::Index>(program->root_schedule_nodes.size());
  std::vector<std::uint8_t> visited(program->nodes.size(), 0U);
  std::vector<semantic::Index> schedule;
  schedule.reserve(program->nodes.size());
  compiled_math_append_schedule_node(*program, node_id, &visited, &schedule);
  program->root_schedule_nodes.insert(
      program->root_schedule_nodes.end(), schedule.begin(), schedule.end());
  program->roots.push_back(
      CompiledMathRoot{
          node_id,
          CompiledMathIndexSpan{
              offset,
              static_cast<semantic::Index>(schedule.size())}});
  return root_id;
}

inline void compiled_math_release_planning_fields(
    CompiledMathProgram *program) {
  decltype(program->node_index)().swap(program->node_index);
  decltype(program->condition_index)().swap(program->condition_index);
}

} // namespace detail
} // namespace accumulatr::eval
