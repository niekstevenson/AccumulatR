#pragma once

#include "compiled_math_kernel_planning.hpp"

namespace accumulatr::eval::detail {

inline void compiled_integral_source_leaves(
    const CompiledMathProgram &program, const semantic::Index id,
    std::vector<semantic::Index> *leaves) {
  if (id < 0) return;
  const auto &source = program.source_programs[id];
  if (source.leaf_index >= 0) leaves->push_back(source.leaf_index);
  compiled_integral_source_leaves(program, source.child_program_id, leaves);
  compiled_integral_source_leaves(program, source.onset_source_program_id, leaves);
  for (semantic::Index i = 0; i < source.member_programs.size; ++i) {
    compiled_integral_source_leaves(program,
        program.source_program_members[source.member_programs.offset + i], leaves);
  }
}

inline bool compiled_integral_dependencies(
    const CompiledMathProgram &program, const semantic::Index id,
    std::vector<semantic::Index> *bound,
    std::vector<semantic::Index> *leaves) {
  const auto &node = program.nodes[id];
  const auto local_time = [&](semantic::Index time) {
    return time == semantic::kInvalidIndex ||
        time == static_cast<semantic::Index>(CompiledMathTimeSlot::Zero) ||
        std::find(bound->begin(), bound->end(), time) != bound->end();
  };
  if (compiled_math_is_integral_node(node.kind)) {
    if (!local_time(node.time_id)) return false;
    bound->push_back(compiled_math_integral_bind_time_id(node));
    const bool valid = compiled_integral_dependencies(program,
        program.roots[node.subject_id].node_id, bound, leaves);
    bound->pop_back();
    return valid;
  }
  if (compiled_math_is_source_value_node(node.kind)) {
    if (!local_time(node.time_id) || !local_time(node.aux_id)) return false;
    compiled_integral_source_leaves(program, node.source_program_id, leaves);
  } else if (node.kind == CompiledMathNodeKind::TimeGate ||
             node.kind == CompiledMathNodeKind::StrictTimeGate) {
    if (!local_time(node.time_id) || !local_time(node.aux_id)) return false;
  } else if (node.kind == CompiledMathNodeKind::ExprUpperBoundDensity ||
             node.kind == CompiledMathNodeKind::ExprUpperBoundCdf) {
    return false;
  }
  for (semantic::Index i = 0; i < node.children.size; ++i) {
    if (!compiled_integral_dependencies(program,
            program.child_nodes[node.children.offset + i], bound, leaves)) return false;
  }
  return true;
}

inline void compile_cumulative_integral_dependencies(CompiledMathProgram *program) {
  for (auto &kernel : program->integral_kernels) {
    std::vector<semantic::Index> bound{kernel.bind_time_id};
    auto &leaves = kernel.cumulative_leaves;
    if (!compiled_integral_dependencies(*program,
            program->roots[kernel.root_id].node_id, &bound, &leaves)) {
      leaves.clear();
    }
    std::sort(leaves.begin(), leaves.end());
    leaves.erase(std::unique(leaves.begin(), leaves.end()), leaves.end());
    const auto support = [&](const CompiledMathSourceProductOps ops) {
      for (semantic::Index i = 0; i < ops.size; ++i) {
        const auto &op = program->source_product_ops[ops.offset + i];
        if (!(op.value_channel_mask & 3U) || op.time_id != kernel.bind_time_id) continue;
        const auto &source = program->source_programs[op.source_product_program_id];
        const auto initial = op.value_channel_mask & 1U
            ? source.initial_with_pdf_program_id : source.initial_without_pdf_program_id;
        if (initial >= 0 && program->source_programs[initial].kind ==
                CompiledMathSourceProductProgramKind::LeafAbsolute) {
          kernel.support_leaves.push_back(program->source_programs[initial].leaf_index);
        }
      }
    };
    support(kernel.initial_execution.source_product_ops);
    if (kernel.initial_execution.kind == CompiledMathExecutionKind::SourceProductSum &&
        kernel.initial_execution.source_product_terms.size == 1) {
      support(program->source_product_terms[
          kernel.initial_execution.source_product_terms.offset].source_product_ops);
    }
    auto &bounds = kernel.support_leaves;
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
  }
}

} // namespace accumulatr::eval::detail
