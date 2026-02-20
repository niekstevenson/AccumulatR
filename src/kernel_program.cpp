#include "kernel_program.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace uuber {
namespace {

inline KernelOpCode kernel_code_from_ir(IrNodeOp op) {
  switch (op) {
  case IrNodeOp::EventAcc:
  case IrNodeOp::EventPool:
    return KernelOpCode::Event;
  case IrNodeOp::And:
    return KernelOpCode::And;
  case IrNodeOp::Or:
    return KernelOpCode::Or;
  case IrNodeOp::Not:
    return KernelOpCode::Not;
  case IrNodeOp::Guard:
    return KernelOpCode::Guard;
  default:
    return KernelOpCode::Event;
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
      int child = ir.node_children[static_cast<std::size_t>(node.child_begin + i)];
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

KernelProgram compile_kernel_program(const IrContext &ir) {
  KernelProgram program;
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
  program.outputs.slot_to_node_idx.reserve(topo.size());
  program.outputs.outcome_idx_to_slot.assign(ir.outcomes.size(), -1);
  program.ops.reserve(topo.size());

  for (std::size_t slot = 0; slot < topo.size(); ++slot) {
    int node_idx = topo[slot];
    program.outputs.node_idx_to_slot[static_cast<std::size_t>(node_idx)] =
        static_cast<int>(slot);
    program.outputs.slot_to_node_idx.push_back(node_idx);
  }
  for (std::size_t oi = 0; oi < ir.outcomes.size(); ++oi) {
    int node_idx = ir.outcomes[oi].node_idx;
    if (node_idx < 0 ||
        node_idx >= static_cast<int>(program.outputs.node_idx_to_slot.size())) {
      continue;
    }
    program.outputs.outcome_idx_to_slot[oi] =
        program.outputs.node_idx_to_slot[static_cast<std::size_t>(node_idx)];
  }

  for (std::size_t slot = 0; slot < topo.size(); ++slot) {
    int node_idx = topo[slot];
    const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];

    KernelOp op;
    op.code = kernel_code_from_ir(node.op);
    op.node_idx = node_idx;
    op.out_slot = static_cast<int>(slot);
    op.event_idx = node.event_idx;
    op.flags = node.flags;
    op.ref_slot = (node.reference_idx >= 0 &&
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
        int child_idx =
            ir.node_children[static_cast<std::size_t>(node.child_begin + i)];
        if (child_idx >= 0 &&
            child_idx <
                static_cast<int>(program.outputs.node_idx_to_slot.size())) {
          int child_slot = program.outputs.node_idx_to_slot[static_cast<std::size_t>(
              child_idx)];
          program.children.push_back(child_slot);
        } else {
          program.children.push_back(-1);
        }
      }
    }
    if (op.code == KernelOpCode::Guard) {
      program.has_guard = true;
    }
    program.ops.push_back(op);
  }

  program.valid = !program.ops.empty();
  return program;
}

} // namespace uuber
