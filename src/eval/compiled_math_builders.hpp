#pragma once

#include <algorithm>
#include <initializer_list>
#include <vector>

#include "compiled_math_kernel_planning.hpp"

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
  node.subject_id = key.subject_id;
  node.condition_id = key.condition_id;
  node.time_id = key.time_id;
  node.aux_id = key.aux_id;
  node.aux2_id = key.aux2_id;
  node.source_view_id = key.source_view_id;
  node.children = CompiledMathIndexSpan{
      child_offset,
      static_cast<semantic::Index>(key.children.size())};
  node.integral_kernel_slot = compiled_math_integral_kernel_slot(program, node);
  node.constant = key.constant;
  program->nodes.push_back(node);
  program->node_index.emplace(std::move(key), node_id);
  return node_id;
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
  if (!key.impossible && key.source_ids.empty()) {
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
          key.relations});
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
  for (auto &kernel : program->integral_kernels) {
    kernel.execution.source_value_factors = CompiledMathIndexSpan{};
    kernel.initial_execution.source_value_factors = CompiledMathIndexSpan{};
  }
  for (auto &root : program->roots) {
    root.execution.source_value_factors = CompiledMathIndexSpan{};
    root.initial_execution.source_value_factors = CompiledMathIndexSpan{};
  }
  for (auto &term : program->source_product_terms) {
    term.source_value_factors = CompiledMathIndexSpan{};
  }
  decltype(program->conditions)().swap(program->conditions);
  decltype(program->source_value_factors)().swap(
      program->source_value_factors);
  decltype(program->source_product_channels)().swap(
      program->source_product_channels);
  decltype(program->node_index)().swap(program->node_index);
  decltype(program->condition_index)().swap(program->condition_index);
}

} // namespace detail
} // namespace accumulatr::eval
