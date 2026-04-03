#include "vector_runtime.h"

#include <cstddef>
#include <vector>

namespace uuber {
namespace {

inline VectorOpCode tree_vector_code_from_ir(IrNodeOp op) {
  switch (op) {
  case IrNodeOp::EventAcc:
  case IrNodeOp::EventPool:
    return VectorOpCode::TreeEvent;
  case IrNodeOp::And:
    return VectorOpCode::TreeAnd;
  case IrNodeOp::Or:
    return VectorOpCode::TreeOr;
  case IrNodeOp::Not:
    return VectorOpCode::TreeNot;
  case IrNodeOp::Guard:
    return VectorOpCode::TreeGuard;
  default:
    return VectorOpCode::TreeEvent;
  }
}

void dfs_topo(const IrContext &ir, int node_idx, std::vector<int> &marks,
              std::vector<int> &topo, bool &has_cycle) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ir.nodes.size())) {
    return;
  }
  if (marks[static_cast<std::size_t>(node_idx)] == 2) {
    return;
  }
  if (marks[static_cast<std::size_t>(node_idx)] == 1) {
    has_cycle = true;
    return;
  }

  marks[static_cast<std::size_t>(node_idx)] = 1;
  const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];

  if (node.child_begin >= 0 && node.child_count > 0 &&
      node.child_begin + node.child_count <=
          static_cast<int>(ir.node_children.size())) {
    for (int i = 0; i < node.child_count; ++i) {
      const int child =
          ir.node_children[static_cast<std::size_t>(node.child_begin + i)];
      dfs_topo(ir, child, marks, topo, has_cycle);
    }
  }
  if (node.reference_idx >= 0) {
    dfs_topo(ir, node.reference_idx, marks, topo, has_cycle);
  }
  if (node.blocker_idx >= 0) {
    dfs_topo(ir, node.blocker_idx, marks, topo, has_cycle);
  }

  marks[static_cast<std::size_t>(node_idx)] = 2;
  topo.push_back(node_idx);
}

} // namespace

VectorProgram compile_tree_vector_program(const IrContext &ir) {
  VectorProgram program;
  if (!ir.valid || ir.nodes.empty()) {
    return program;
  }

  std::vector<int> marks(ir.nodes.size(), 0);
  std::vector<int> topo;
  topo.reserve(ir.nodes.size());
  bool has_cycle = false;
  for (int node_idx = 0; node_idx < static_cast<int>(ir.nodes.size()); ++node_idx) {
    if (marks[static_cast<std::size_t>(node_idx)] == 0) {
      dfs_topo(ir, node_idx, marks, topo, has_cycle);
    }
  }
  if (has_cycle) {
    return program;
  }

  program.outputs.node_idx_to_slot.assign(ir.nodes.size(), -1);
  program.ops.reserve(topo.size());

  for (std::size_t slot = 0; slot < topo.size(); ++slot) {
    const int node_idx = topo[slot];
    program.outputs.node_idx_to_slot[static_cast<std::size_t>(node_idx)] =
        static_cast<int>(slot);
  }

  for (std::size_t slot = 0; slot < topo.size(); ++slot) {
    const int node_idx = topo[slot];
    const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];

    VectorOp op;
    op.code = tree_vector_code_from_ir(node.op);
    op.dst_slot = static_cast<int>(slot);
    op.node_idx = node_idx;
    op.event_idx = node.event_idx;
    op.flags = node.flags;
    op.ref_slot =
        (node.reference_idx >= 0 &&
         node.reference_idx <
             static_cast<int>(program.outputs.node_idx_to_slot.size()))
            ? program.outputs.node_idx_to_slot[static_cast<std::size_t>(
                  node.reference_idx)]
            : -1;
    op.blocker_slot =
        (node.blocker_idx >= 0 &&
         node.blocker_idx <
             static_cast<int>(program.outputs.node_idx_to_slot.size()))
            ? program.outputs.node_idx_to_slot[static_cast<std::size_t>(
                  node.blocker_idx)]
            : -1;
    if (node.child_begin >= 0 && node.child_count > 0 &&
        node.child_begin + node.child_count <=
            static_cast<int>(ir.node_children.size())) {
      op.child_begin = static_cast<int>(program.children.size());
      op.child_count = node.child_count;
      if (op.child_count > program.max_child_count) {
        program.max_child_count = op.child_count;
      }
      for (int i = 0; i < node.child_count; ++i) {
        const int child_idx =
            ir.node_children[static_cast<std::size_t>(node.child_begin + i)];
        if (child_idx >= 0 &&
            child_idx <
                static_cast<int>(program.outputs.node_idx_to_slot.size())) {
          program.children.push_back(
              program.outputs.node_idx_to_slot[static_cast<std::size_t>(
                  child_idx)]);
        } else {
          program.children.push_back(-1);
        }
      }
    }
    program.ops.push_back(op);
  }

  program.valid = !program.ops.empty();
  return program;
}

} // namespace uuber
