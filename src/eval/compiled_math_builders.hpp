#pragma once

#include <algorithm>
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

inline semantic::Index compiled_math_source_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index source_id,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index time_cap_id = semantic::kInvalidIndex) {
  CompiledMathNodeKey key;
  key.kind = kind;
  key.subject_id = source_id;
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

inline semantic::Index compiled_math_make_root(CompiledMathProgram *program,
                                              semantic::Index node_id);

// Free-time dependence respects the binding scope of nested integrals.
inline bool compiled_math_depends_on_integral_scope(
    const CompiledMathProgram &program, const semantic::Index node_id,
    const semantic::Index time_id, const semantic::Index source_view_id = 0) {
  const auto &node = program.nodes[node_id];
  if (compiled_math_is_integral_node(node.kind)) {
    return node.time_id == time_id ||
        (compiled_math_integral_bind_time_id(node) != time_id &&
         compiled_math_depends_on_integral_scope(program,
             program.roots[node.subject_id].node_id, time_id));
  }
  if (compiled_math_is_source_value_node(node.kind) &&
      source_view_id != 0 && node.source_view_id == 0) return true;
  if (compiled_math_is_source_value_node(node.kind) ||
      node.kind == CompiledMathNodeKind::TimeGate ||
      node.kind == CompiledMathNodeKind::StrictTimeGate ||
      node.kind == CompiledMathNodeKind::ExprUpperBoundDensity ||
      node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
    if (node.time_id == time_id || node.aux_id == time_id) return true;
  }
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    if (compiled_math_depends_on_integral_scope(program,
            program.child_nodes[node.children.offset + i], time_id, source_view_id)) return true;
  }
  return false;
}

inline void compiled_math_partition_integrand(
    const CompiledMathProgram &program, const semantic::Index node_id,
    const semantic::Index upper_time, const semantic::Index bind_time,
    std::vector<semantic::Index> *inside,
    std::vector<semantic::Index> *outside, bool *stripped_gate,
    const semantic::Index source_view_id) {
  const auto &node = program.nodes[node_id];
  // The quadrature domain already enforces 0 < s < upper.
  if ((node.kind == CompiledMathNodeKind::TimeGate ||
       node.kind == CompiledMathNodeKind::StrictTimeGate) &&
      upper_time != bind_time &&
      node.time_id == upper_time && node.aux_id == bind_time) {
    *stripped_gate = true;
    compiled_math_partition_integrand(program,
        program.child_nodes[node.children.offset], upper_time, bind_time,
        inside, outside, stripped_gate, source_view_id);
  } else if (node.kind == CompiledMathNodeKind::Product) {
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      compiled_math_partition_integrand(program,
          program.child_nodes[node.children.offset + i], upper_time, bind_time,
          inside, outside, stripped_gate, source_view_id);
    }
  } else {
    (compiled_math_depends_on_integral_scope(program, node_id, bind_time, source_view_id)
         ? inside : outside)->push_back(node_id);
  }
}

inline semantic::Index compiled_math_integral_node(
    CompiledMathProgram *program,
    const CompiledMathNodeKind kind,
    const semantic::Index integrand_root_id,
    const semantic::Index time_id = 0,
    const semantic::Index source_view_id = 0,
    const semantic::Index bind_time_id = semantic::kInvalidIndex) {
  const auto bind = bind_time_id == semantic::kInvalidIndex ? time_id : bind_time_id;
  std::vector<semantic::Index> inside, outside;
  bool stripped_gate = false;
  compiled_math_partition_integrand(*program,
      program->roots[integrand_root_id].node_id, time_id, bind, &inside, &outside,
      &stripped_gate, source_view_id);
  const auto root = outside.empty() && !stripped_gate ? integrand_root_id
      : compiled_math_make_root(program, compiled_math_algebra_node(
            program, CompiledMathNodeKind::Product, std::move(inside)));
  CompiledMathNodeKey key;
  key.kind = outside.empty() ? kind : CompiledMathNodeKind::IntegralZeroToCurrentRaw;
  key.value_kind = key.kind == CompiledMathNodeKind::IntegralZeroToCurrent
      ? CompiledMathValueKind::Cdf : CompiledMathValueKind::Scalar;
  key.subject_id = root;
  key.time_id = time_id;
  key.aux2_id = bind_time_id;
  key.source_view_id = source_view_id;
  const auto integral = compiled_math_intern_node(program, std::move(key));
  if (outside.empty()) return integral;
  outside.push_back(integral);
  const auto product = compiled_math_algebra_node(
      program, CompiledMathNodeKind::Product, std::move(outside));
  return kind == CompiledMathNodeKind::IntegralZeroToCurrent
      ? compiled_math_unary_node(program, CompiledMathNodeKind::ClampProbability,
                                product, CompiledMathValueKind::Cdf)
      : product;
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
  if (node.kind != CompiledMathNodeKind::OutcomeSelect) {
    for (semantic::Index i = 0; i < node.children.size; ++i) {
      const auto child_id = program.child_nodes[node.children.offset + i];
      compiled_math_append_schedule_node(program, child_id, visited, schedule);
    }
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
  decltype(program->source_value_factors)().swap(
      program->source_value_factors);
  decltype(program->source_product_channels)().swap(
      program->source_product_channels);
  decltype(program->node_index)().swap(program->node_index);
}

} // namespace detail
} // namespace accumulatr::eval
