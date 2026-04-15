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

inline std::uint8_t need_mask(EvalNeed need) {
  return static_cast<std::uint8_t>(need);
}

inline bool has_any_need_mask(std::uint8_t mask) { return mask != 0u; }

inline bool mask_needs_density(std::uint8_t mask) {
  return (mask & need_mask(EvalNeed::kDensity)) != 0u;
}

inline bool mask_needs_survival(std::uint8_t mask) {
  return (mask & need_mask(EvalNeed::kSurvival)) != 0u;
}

inline bool mask_needs_cdf(std::uint8_t mask) {
  return (mask & need_mask(EvalNeed::kCDF)) != 0u;
}

inline void propagate_need(std::vector<std::uint8_t> &slot_need_masks, int slot,
                           std::uint8_t mask) {
  if (slot < 0 || !has_any_need_mask(mask)) {
    return;
  }
  slot_need_masks[static_cast<std::size_t>(slot)] |= mask;
}

inline std::size_t single_query_plan_offset(int target_slot, std::uint8_t need) {
  return static_cast<std::size_t>(target_slot) * 8u +
         static_cast<std::size_t>(need);
}

KernelQueryPlan compile_kernel_query_plan_from_slots(
    const KernelProgram &program, const std::vector<int> &target_slots,
    const std::vector<std::uint8_t> &target_need_masks) {
  KernelQueryPlan plan;
  if (target_slots.empty() || target_slots.size() != target_need_masks.size() ||
      !program.valid) {
    return plan;
  }

  plan.target_slots = target_slots;
  plan.slot_need_masks.assign(program.ops.size(), 0u);
  plan.max_slot = -1;

  for (std::size_t i = 0; i < target_slots.size(); ++i) {
    const int slot = target_slots[i];
    if (slot < 0 || slot >= static_cast<int>(program.ops.size())) {
      return KernelQueryPlan{};
    }
    plan.slot_need_masks[static_cast<std::size_t>(slot)] |= target_need_masks[i];
    if (slot > plan.max_slot) {
      plan.max_slot = slot;
    }
  }

  for (int op_idx = plan.max_slot; op_idx >= 0; --op_idx) {
    const std::uint8_t slot_need =
        plan.slot_need_masks[static_cast<std::size_t>(op_idx)];
    if (!has_any_need_mask(slot_need)) {
      continue;
    }
    const KernelOp &op = program.ops[static_cast<std::size_t>(op_idx)];
    switch (op.code) {
    case KernelOpCode::Event:
      break;
    case KernelOpCode::And: {
      std::uint8_t child_need = 0u;
      if (mask_needs_density(slot_need)) {
        child_need |= need_mask(EvalNeed::kDensity);
        child_need |= need_mask(EvalNeed::kCDF);
      }
      if (mask_needs_survival(slot_need) || mask_needs_cdf(slot_need)) {
        child_need |= need_mask(EvalNeed::kCDF);
      }
      if (!has_any_need_mask(child_need) || op.child_begin < 0 ||
          op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        break;
      }
      for (int i = 0; i < op.child_count; ++i) {
        propagate_need(
            plan.slot_need_masks,
            program.children[static_cast<std::size_t>(op.child_begin + i)],
            child_need);
      }
      break;
    }
    case KernelOpCode::Or: {
      std::uint8_t child_need = 0u;
      if (mask_needs_density(slot_need)) {
        child_need |= need_mask(EvalNeed::kDensity);
        child_need |= need_mask(EvalNeed::kSurvival);
      }
      if (mask_needs_survival(slot_need) || mask_needs_cdf(slot_need)) {
        child_need |= need_mask(EvalNeed::kSurvival);
      }
      if (!has_any_need_mask(child_need) || op.child_begin < 0 ||
          op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        break;
      }
      for (int i = 0; i < op.child_count; ++i) {
        propagate_need(
            plan.slot_need_masks,
            program.children[static_cast<std::size_t>(op.child_begin + i)],
            child_need);
      }
      break;
    }
    case KernelOpCode::Not:
      if (!(mask_needs_survival(slot_need) || mask_needs_cdf(slot_need)) ||
          op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin >= static_cast<int>(program.children.size())) {
        break;
      }
      propagate_need(
          plan.slot_need_masks,
          program.children[static_cast<std::size_t>(op.child_begin)],
          need_mask(EvalNeed::kCDF));
      break;
    case KernelOpCode::Guard:
      if (mask_needs_density(slot_need)) {
        propagate_need(plan.slot_need_masks, op.ref_slot,
                       need_mask(EvalNeed::kDensity));
        propagate_need(plan.slot_need_masks, op.blocker_slot,
                       need_mask(EvalNeed::kSurvival));
      }
      break;
    }
  }

  plan.valid = plan.max_slot >= 0;
  return plan;
}

const KernelQueryPlan *lookup_single_kernel_query_plan(const KernelProgram &program,
                                                       int target_slot,
                                                       std::uint8_t need) {
  if (!program.valid || target_slot < 0 ||
      target_slot >= static_cast<int>(program.ops.size()) || need == 0u) {
    return nullptr;
  }
  const KernelQueryPlan &plan =
      program.single_query_plans[single_query_plan_offset(target_slot, need)];
  if (plan.valid) {
    return &plan;
  }
  return nullptr;
}

const KernelQueryPlan *ensure_single_kernel_query_plan(KernelProgram &program,
                                                       int target_slot,
                                                       std::uint8_t need) {
  const KernelQueryPlan *existing =
      lookup_single_kernel_query_plan(program, target_slot, need);
  if (existing != nullptr) {
    return existing;
  }
  KernelQueryPlan plan = compile_kernel_query_plan_from_slots(
      program, std::vector<int>{target_slot}, std::vector<std::uint8_t>{need});
  if (!plan.valid) {
    return nullptr;
  }
  KernelQueryPlan &stored =
      program.single_query_plans[single_query_plan_offset(target_slot, need)];
  stored = std::move(plan);
  return &stored;
}

KernelMultiQuerySpec make_multi_query_spec(
    const KernelProgram &program, const std::vector<int> &target_node_indices,
    const std::vector<EvalNeed> &target_needs) {
  KernelMultiQuerySpec spec;
  if (!program.valid || target_node_indices.size() != target_needs.size()) {
    return spec;
  }
  spec.target_slots.reserve(target_node_indices.size());
  spec.target_need_masks.reserve(target_needs.size());
  for (std::size_t i = 0; i < target_node_indices.size(); ++i) {
    const int node_idx = target_node_indices[i];
    if (node_idx < 0 ||
        node_idx >= static_cast<int>(program.outputs.node_idx_to_slot.size()) ||
        !has_any_need(target_needs[i])) {
      spec.target_slots.clear();
      spec.target_need_masks.clear();
      return spec;
    }
    const int slot = program.outputs.node_idx_to_slot[static_cast<std::size_t>(
        target_node_indices[i])];
    if (slot < 0 || slot >= static_cast<int>(program.ops.size())) {
      spec.target_slots.clear();
      spec.target_need_masks.clear();
      return spec;
    }
    spec.target_slots.push_back(slot);
    spec.target_need_masks.push_back(need_mask(target_needs[i]));
  }
  return spec;
}

int lookup_precomputed_multi_query_plan_index(const KernelProgram &program,
                                              const KernelMultiQuerySpec &spec) {
  if (spec.target_slots.empty() ||
      spec.target_slots.size() != spec.target_need_masks.size()) {
    return -1;
  }
  for (int i = 0; i < static_cast<int>(program.multi_query_plans.size()); ++i) {
    const KernelMultiQueryPlan &entry =
        program.multi_query_plans[static_cast<std::size_t>(i)];
    if (entry.spec == spec) {
      return i;
    }
  }
  return -1;
}

int precompute_multi_query_plan(KernelProgram &program,
                                const KernelMultiQuerySpec &spec) {
  const int existing_index =
      lookup_precomputed_multi_query_plan_index(program, spec);
  if (existing_index >= 0) {
    return existing_index;
  }
  KernelQueryPlan plan = compile_kernel_query_plan_from_slots(
      program, spec.target_slots, spec.target_need_masks);
  if (!plan.valid) {
    return -1;
  }
  program.multi_query_plans.push_back(KernelMultiQueryPlan{spec, std::move(plan)});
  return static_cast<int>(program.multi_query_plans.size()) - 1;
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
  if (program.valid) {
    program.single_query_plans.assign(program.ops.size() * 8u,
                                      KernelQueryPlan{});
  }
  return program;
}

const KernelQueryPlan *precompute_kernel_query_plan(KernelProgram &program,
                                                    int target_node_idx,
                                                    EvalNeed need) {
  if (!program.valid || target_node_idx < 0 || !has_any_need(need) ||
      target_node_idx >=
          static_cast<int>(program.outputs.node_idx_to_slot.size())) {
    return nullptr;
  }
  const int target_slot = program.outputs.node_idx_to_slot
      [static_cast<std::size_t>(target_node_idx)];
  return ensure_single_kernel_query_plan(program, target_slot, need_mask(need));
}

int precompute_kernel_multi_query_plan(KernelProgram &program,
                                       const std::vector<int> &target_node_indices,
                                       const std::vector<EvalNeed> &target_needs) {
  if (!program.valid || target_node_indices.empty() ||
      target_node_indices.size() != target_needs.size()) {
    return -1;
  }

  if (target_node_indices.size() == 1u) {
    return -1;
  }

  KernelMultiQuerySpec spec =
      make_multi_query_spec(program, target_node_indices, target_needs);
  if (spec.target_slots.empty()) {
    return -1;
  }
  return precompute_multi_query_plan(program, spec);
}

const KernelQueryPlan *get_kernel_query_plan(const KernelProgram &program,
                                             int target_node_idx,
                                             EvalNeed need) {
  if (!program.valid || target_node_idx < 0 || !has_any_need(need) ||
      target_node_idx >=
          static_cast<int>(program.outputs.node_idx_to_slot.size())) {
    return nullptr;
  }
  const int target_slot = program.outputs.node_idx_to_slot
      [static_cast<std::size_t>(target_node_idx)];
  return lookup_single_kernel_query_plan(program, target_slot, need_mask(need));
}

const KernelQueryPlan *get_kernel_multi_query_plan(const KernelProgram &program,
                                                   int query_plan_index) {
  if (!program.valid || query_plan_index < 0 ||
      query_plan_index >= static_cast<int>(program.multi_query_plans.size())) {
    return nullptr;
  }
  const KernelQueryPlan &plan =
      program.multi_query_plans[static_cast<std::size_t>(query_plan_index)].plan;
  return plan.valid ? &plan : nullptr;
}

} // namespace uuber
