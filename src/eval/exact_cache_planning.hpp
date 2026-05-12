#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "compiled_math_kernel_planning.hpp"
#include "compiled_math_workspace.hpp"
#include "exact_source_view.hpp"

namespace accumulatr::eval {
namespace detail {

inline std::size_t exact_condition_cache_id(
    const CompiledSourceView *evaluator) {
  return evaluator == nullptr
             ? 0U
             : static_cast<std::size_t>(evaluator->condition_cache_id());
}

inline void exact_hash_condition_part(std::size_t *seed,
                                      const std::size_t value) noexcept {
  *seed ^= value + 0x9e3779b97f4a7c15ULL + (*seed << 6U) + (*seed >> 2U);
}

inline std::size_t compiled_node_condition_cache_id(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace) {
  const auto evaluator_id = exact_condition_cache_id(evaluator);
  if (node.condition_id == 0 ||
      node.condition_id == semantic::kInvalidIndex) {
    return evaluator_id;
  }
  const auto condition_pos = static_cast<std::size_t>(node.condition_id);
  const auto &cache_plan = program.condition_cache_plans[condition_pos];
  if (!cache_plan.dynamic) {
    return static_cast<std::size_t>(
        compiled_math_static_source_view_condition_cache_id(
            static_cast<semantic::Index>(evaluator_id),
            node.condition_id));
  }
  std::size_t seed = 0x9e3779b97f4a7c15ULL;
  exact_hash_condition_part(&seed, static_cast<std::size_t>(evaluator_id));
  exact_hash_condition_part(&seed, static_cast<std::size_t>(node.condition_id));
  for (semantic::Index i = 0; i < cache_plan.time_dependencies.size; ++i) {
    const auto time_id =
        program.condition_cache_time_dependencies[
            static_cast<std::size_t>(
                cache_plan.time_dependencies.offset + i)];
    const auto time_pos = static_cast<std::size_t>(time_id);
    exact_hash_condition_part(&seed, static_cast<std::size_t>(time_id) + 1U);
    exact_hash_condition_part(
        &seed, std::hash<double>{}(workspace.time_values[time_pos]));
  }
  return seed;
}

inline bool compiled_math_outcome_gate_open_for_node(
    const CompiledMathNode &node,
    const CompiledMathWorkspace &workspace,
    const CompiledSourceView *parent);

inline bool compiled_math_node_cacheable(
    const CompiledMathProgram &,
    const CompiledMathNode &node) noexcept {
  return node.kind == CompiledMathNodeKind::IntegralZeroToCurrent ||
         node.kind == CompiledMathNodeKind::IntegralZeroToCurrentRaw;
}

inline bool compiled_math_cache_time_excluded(
    const semantic::Index time_id,
    const std::vector<semantic::Index> &excluded_time_ids) noexcept {
  return std::find(
             excluded_time_ids.begin(),
             excluded_time_ids.end(),
             time_id) != excluded_time_ids.end();
}

inline void compiled_math_hash_time_cache_dependency(
    std::size_t *seed,
    const semantic::Index time_id,
    const CompiledMathWorkspace &workspace,
    const std::vector<semantic::Index> &excluded_time_ids) {
  if (time_id == semantic::kInvalidIndex ||
      compiled_math_cache_time_excluded(time_id, excluded_time_ids)) {
    return;
  }
  exact_hash_condition_part(seed, static_cast<std::size_t>(time_id) + 1U);
  const auto time_pos = static_cast<std::size_t>(time_id);
  if (time_pos < workspace.time_valid.size() &&
      workspace.time_valid[time_pos] != 0U) {
    exact_hash_condition_part(
        seed, std::hash<double>{}(workspace.time_values[time_pos]));
  } else {
    exact_hash_condition_part(seed, 0x7fbadbadULL);
  }
}

inline void compiled_math_hash_node_cache_dependencies(
    const CompiledMathProgram &program,
    semantic::Index node_id,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    std::vector<semantic::Index> *excluded_time_ids,
    std::size_t *seed,
    std::vector<std::uint8_t> *visiting) {
  if (node_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(node_id) >= program.nodes.size()) {
    return;
  }
  const auto node_pos = static_cast<std::size_t>(node_id);
  if ((*visiting)[node_pos] != 0U) {
    return;
  }
  (*visiting)[node_pos] = 1U;
  const auto &node = program.nodes[node_pos];
  exact_hash_condition_part(
      seed,
      compiled_node_condition_cache_id(program, node, evaluator, workspace));

  switch (node.kind) {
  case CompiledMathNodeKind::SourcePdf:
  case CompiledMathNodeKind::SourceCdf:
  case CompiledMathNodeKind::SourceSurvival:
    compiled_math_hash_time_cache_dependency(
        seed, node.time_id, workspace, *excluded_time_ids);
    compiled_math_hash_time_cache_dependency(
        seed, node.aux_id, workspace, *excluded_time_ids);
    break;
  case CompiledMathNodeKind::TimeGate:
  case CompiledMathNodeKind::StrictTimeGate:
    compiled_math_hash_time_cache_dependency(
        seed, node.time_id, workspace, *excluded_time_ids);
    compiled_math_hash_time_cache_dependency(
        seed, node.aux_id, workspace, *excluded_time_ids);
    break;
  case CompiledMathNodeKind::ExprUpperBoundDensity:
  case CompiledMathNodeKind::ExprUpperBoundCdf:
    compiled_math_hash_time_cache_dependency(
        seed, node.time_id, workspace, *excluded_time_ids);
    if (node.aux_id != semantic::kInvalidIndex &&
        node.aux2_id != semantic::kInvalidIndex) {
      const auto offset = static_cast<std::size_t>(node.aux_id);
      const auto size = static_cast<std::size_t>(node.aux2_id);
      for (std::size_t i = 0; i < size &&
                              offset + i < program.timed_upper_bound_terms.size();
           ++i) {
        const auto &term = program.timed_upper_bound_terms[offset + i];
        compiled_math_hash_time_cache_dependency(
            seed, term.time_id, workspace, *excluded_time_ids);
        compiled_math_hash_node_cache_dependencies(
            program,
            term.normalizer_node_id,
            evaluator,
            workspace,
            excluded_time_ids,
            seed,
            visiting);
      }
    }
    break;
  case CompiledMathNodeKind::OutcomeSubsetUnused:
  case CompiledMathNodeKind::OutcomeSubsetUsed:
    exact_hash_condition_part(
        seed,
        compiled_math_outcome_gate_open_for_node(node, workspace, evaluator)
            ? 1U
            : 0U);
    break;
  case CompiledMathNodeKind::IntegralZeroToCurrent:
  case CompiledMathNodeKind::IntegralZeroToCurrentRaw: {
    compiled_math_hash_time_cache_dependency(
        seed, node.time_id, workspace, *excluded_time_ids);
    const auto root_id = compiled_math_integral_root_id(node);
    if (root_id != semantic::kInvalidIndex &&
        static_cast<std::size_t>(root_id) < program.roots.size()) {
      excluded_time_ids->push_back(compiled_math_integral_bind_time_id(node));
      compiled_math_hash_node_cache_dependencies(
          program,
          program.roots[static_cast<std::size_t>(root_id)].node_id,
          evaluator,
          workspace,
          excluded_time_ids,
          seed,
          visiting);
      excluded_time_ids->pop_back();
    }
    break;
  }
  default:
    break;
  }

  for (semantic::Index i = 0; i < node.children.size; ++i) {
    const auto child_pos =
        static_cast<std::size_t>(node.children.offset + i);
    if (child_pos < program.child_nodes.size()) {
      compiled_math_hash_node_cache_dependencies(
          program,
          program.child_nodes[child_pos],
          evaluator,
          workspace,
          excluded_time_ids,
          seed,
          visiting);
    }
  }
  (*visiting)[node_pos] = 0U;
}

inline std::size_t compiled_math_node_cache_dependency_id(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace) {
  std::size_t seed = 0x243f6a8885a308d3ULL;
  exact_hash_condition_part(
      &seed,
      compiled_node_condition_cache_id(program, node, evaluator, workspace));
  if (!compiled_math_is_integral_node(node.kind)) {
    return seed;
  }
  const auto root_id = compiled_math_integral_root_id(node);
  if (root_id == semantic::kInvalidIndex ||
      static_cast<std::size_t>(root_id) >= program.roots.size()) {
    return seed;
  }
  std::vector<semantic::Index> excluded_time_ids{
      compiled_math_integral_bind_time_id(node)};
  std::vector<std::uint8_t> visiting(program.nodes.size(), 0U);
  compiled_math_hash_node_cache_dependencies(
      program,
      program.roots[static_cast<std::size_t>(root_id)].node_id,
      evaluator,
      workspace,
      &excluded_time_ids,
      &seed,
      &visiting);
  return seed;
}

inline bool compiled_math_load_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    const CompiledMathWorkspace &workspace,
    double *value);

inline void compiled_math_store_node_cache_entry(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    const CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace,
    double value);

} // namespace detail
} // namespace accumulatr::eval
