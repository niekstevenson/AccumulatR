#include "coupling_program.h"

#include <Rcpp.h>

#include "native_utils.h"

namespace {

inline std::uint64_t ir_mix_lookup_hash(std::uint64_t hash, int value) {
  const std::uint32_t v = static_cast<std::uint32_t>(value);
  for (int i = 0; i < 4; ++i) {
    const std::uint8_t b =
        static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
    hash ^= static_cast<std::uint64_t>(b);
    hash *= kFNV64Prime;
  }
  return hash;
}

inline std::uint64_t ir_outcome_coupling_lookup_key(
    int node_idx, const std::vector<int> &competitors) {
  std::uint64_t hash = kFNV64Offset;
  hash = ir_mix_lookup_hash(hash, node_idx);
  hash = ir_mix_lookup_hash(hash, static_cast<int>(competitors.size()));
  for (int comp : competitors) {
    hash = ir_mix_lookup_hash(hash, comp);
  }
  return hash;
}

inline int ir_resolve_node_idx(const uuber::NativeContext &ctx, int node_id) {
  if (!ctx.ir.valid) {
    return -1;
  }
  auto it = ctx.ir.id_to_node_idx.find(node_id);
  if (it != ctx.ir.id_to_node_idx.end()) {
    return it->second;
  }
  if (node_id >= 0 && node_id < static_cast<int>(ctx.ir.nodes.size())) {
    return node_id;
  }
  return -1;
}

inline int ir_dense_idx_to_node_id(const uuber::NativeContext &ctx,
                                   int node_idx) {
  if (!ctx.ir.valid || node_idx < 0 ||
      node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return node_idx;
  }
  const int node_id = ctx.ir.nodes[static_cast<std::size_t>(node_idx)].node_id;
  return node_id >= 0 ? node_id : node_idx;
}

inline bool ir_node_source_masks_overlap(const uuber::NativeContext &ctx,
                                         const uuber::IrNode &a,
                                         const uuber::IrNode &b) {
  if (!ctx.ir.valid || a.source_mask_begin < 0 || a.source_mask_count <= 0 ||
      b.source_mask_begin < 0 || b.source_mask_count <= 0 ||
      a.source_mask_count != b.source_mask_count) {
    return false;
  }
  for (int i = 0; i < a.source_mask_count; ++i) {
    const std::uint64_t wa = ctx.ir.node_source_masks[static_cast<std::size_t>(
        a.source_mask_begin + i)];
    const std::uint64_t wb = ctx.ir.node_source_masks[static_cast<std::size_t>(
        b.source_mask_begin + i)];
    if ((wa & wb) != 0u) {
      return true;
    }
  }
  return false;
}

inline bool ir_generic_requires_exact_scenario_eval(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids) {
  if (!ctx.ir.valid) {
    return false;
  }
  const int node_idx = ir_resolve_node_idx(ctx, node_id);
  if (node_idx < 0) {
    return false;
  }
  const uuber::IrNode &target = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  for (int comp_node_id : competitor_node_ids) {
    const int comp_idx = ir_resolve_node_idx(ctx, comp_node_id);
    if (comp_idx < 0 ||
        comp_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      continue;
    }
    const uuber::IrNode &competitor =
        ctx.ir.nodes[static_cast<std::size_t>(comp_idx)];
    // Exact scenario evaluation is only needed when the target scenario can
    // change a competitor through shared sources. Disjoint guard/composite
    // competitors already evaluate correctly on the dense path.
    if (ir_node_source_masks_overlap(ctx, target, competitor)) {
      return true;
    }
  }
  return false;
}

bool ir_outcome_coupling_pair_lookup(const uuber::NativeContext &ctx,
                                     int node_id, int competitor_node_id,
                                     OutcomeCouplingPairPayload &out) {
  if (!ctx.ir.valid) {
    return false;
  }
  const int node_idx = ir_resolve_node_idx(ctx, node_id);
  const int comp_idx = ir_resolve_node_idx(ctx, competitor_node_id);
  if (node_idx < 0 || comp_idx < 0) {
    return false;
  }
  const std::vector<int> competitor_dense{comp_idx};
  const std::uint64_t key =
      ir_outcome_coupling_lookup_key(node_idx, competitor_dense);
  auto it = ctx.ir.outcome_coupling_lookup.find(key);
  if (it == ctx.ir.outcome_coupling_lookup.end()) {
    return false;
  }
  const int spec_idx = it->second;
  if (spec_idx < 0 ||
      spec_idx >= static_cast<int>(ctx.ir.outcome_coupling_ops.size())) {
    return false;
  }
  const uuber::IrOutcomeCouplingOp &spec =
      ctx.ir.outcome_coupling_ops[static_cast<std::size_t>(spec_idx)];
  if (spec.kind != uuber::IrOutcomeCouplingKind::Pair) {
    return false;
  }
  if (spec.node_idx != node_idx || spec.competitor_count != 1) {
    return false;
  }
  if (spec.competitor_begin < 0 ||
      spec.competitor_begin >=
          static_cast<int>(ctx.ir.outcome_competitors.size())) {
    return false;
  }
  if (ctx.ir.outcome_competitors[static_cast<std::size_t>(spec.competitor_begin)] !=
      comp_idx) {
    return false;
  }
  auto fill_from_event = [&](int event_idx,
                             OutcomeCouplingEventRefPayload &payload) -> bool {
    if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
      return false;
    }
    const uuber::IrEvent &event =
        ctx.ir.events[static_cast<std::size_t>(event_idx)];
    if (event.node_idx < 0 ||
        event.node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      return false;
    }
    payload.ref.label_id = event.label_id;
    payload.ref.acc_idx = event.acc_idx;
    payload.ref.pool_idx = event.pool_idx;
    payload.ref.outcome_idx = event.outcome_idx;
    const uuber::IrNode &node =
        ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)];
    payload.node_id = ir_dense_idx_to_node_id(ctx, event.node_idx);
    payload.node_flags = node.flags;
    return true;
  };
  if (!fill_from_event(spec.target_event_idx, out.x_ref)) {
    return false;
  }
  if (!fill_from_event(spec.gate_event_idx, out.c_ref)) {
    return false;
  }
  if (spec.aux_event_count != 1 || spec.aux_event_begin < 0) {
    return false;
  }
  if (spec.aux_event_begin >=
      static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }
  const int y_event_idx = ctx.ir.outcome_coupling_event_indices
      [static_cast<std::size_t>(spec.aux_event_begin)];
  if (!fill_from_event(y_event_idx, out.y_ref)) {
    return false;
  }
  return true;
}

bool ir_outcome_coupling_nway_lookup(const uuber::NativeContext &ctx,
                                     int node_id,
                                     const std::vector<int> &competitor_node_ids,
                                     OutcomeCouplingNWayPayload &out) {
  if (!ctx.ir.valid || competitor_node_ids.empty()) {
    return false;
  }
  const int node_idx = ir_resolve_node_idx(ctx, node_id);
  if (node_idx < 0) {
    return false;
  }
  std::vector<int> competitor_dense;
  competitor_dense.reserve(competitor_node_ids.size());
  for (int comp_id : competitor_node_ids) {
    const int comp_idx = ir_resolve_node_idx(ctx, comp_id);
    if (comp_idx < 0) {
      return false;
    }
    competitor_dense.push_back(comp_idx);
  }
  const std::uint64_t key =
      ir_outcome_coupling_lookup_key(node_idx, competitor_dense);
  auto it = ctx.ir.outcome_coupling_lookup.find(key);
  if (it == ctx.ir.outcome_coupling_lookup.end()) {
    return false;
  }
  const int spec_idx = it->second;
  if (spec_idx < 0 ||
      spec_idx >= static_cast<int>(ctx.ir.outcome_coupling_ops.size())) {
    return false;
  }
  const uuber::IrOutcomeCouplingOp &spec =
      ctx.ir.outcome_coupling_ops[static_cast<std::size_t>(spec_idx)];
  if (spec.kind != uuber::IrOutcomeCouplingKind::NWay) {
    return false;
  }
  if (spec.node_idx != node_idx ||
      spec.competitor_count != static_cast<int>(competitor_dense.size())) {
    return false;
  }
  for (int k = 0; k < spec.competitor_count; ++k) {
    const int idx = spec.competitor_begin + k;
    if (idx < 0 || idx >= static_cast<int>(ctx.ir.outcome_competitors.size())) {
      return false;
    }
    if (ctx.ir.outcome_competitors[static_cast<std::size_t>(idx)] !=
        competitor_dense[static_cast<std::size_t>(k)]) {
      return false;
    }
  }

  auto ref_from_event = [&](int event_idx,
                            OutcomeCouplingEventRefPayload &payload) -> bool {
    if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
      return false;
    }
    const uuber::IrEvent &event =
        ctx.ir.events[static_cast<std::size_t>(event_idx)];
    if (event.node_idx < 0 ||
        event.node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      return false;
    }
    payload.ref.label_id = event.label_id;
    payload.ref.acc_idx = event.acc_idx;
    payload.ref.pool_idx = event.pool_idx;
    payload.ref.outcome_idx = event.outcome_idx;
    const uuber::IrNode &node =
        ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)];
    payload.node_id = ir_dense_idx_to_node_id(ctx, event.node_idx);
    payload.node_flags = node.flags;
    return true;
  };

  if (!ref_from_event(spec.gate_event_idx, out.gate_ref)) {
    return false;
  }
  if (!ref_from_event(spec.target_event_idx, out.target_ref)) {
    return false;
  }
  out.competitor_refs.clear();
  if (spec.aux_event_count <= 0 || spec.aux_event_begin < 0) {
    return false;
  }
  out.competitor_refs.reserve(static_cast<std::size_t>(spec.aux_event_count));
  for (int i = 0; i < spec.aux_event_count; ++i) {
    const int idx = spec.aux_event_begin + i;
    if (idx < 0 ||
        idx >= static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
      return false;
    }
    const int event_idx =
        ctx.ir.outcome_coupling_event_indices[static_cast<std::size_t>(idx)];
    OutcomeCouplingEventRefPayload comp_ref;
    if (!ref_from_event(event_idx, comp_ref)) {
      return false;
    }
    out.competitor_refs.push_back(comp_ref);
  }
  return true;
}

bool ir_outcome_coupling_generic_lookup(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids,
    OutcomeCouplingGenericNodeIntegralPayload &out) {
  if (!ctx.ir.valid || competitor_node_ids.empty()) {
    return false;
  }
  const int node_idx = ir_resolve_node_idx(ctx, node_id);
  if (node_idx < 0) {
    return false;
  }
  std::vector<int> competitor_dense;
  competitor_dense.reserve(competitor_node_ids.size());
  for (int comp_id : competitor_node_ids) {
    const int comp_idx = ir_resolve_node_idx(ctx, comp_id);
    if (comp_idx < 0) {
      return false;
    }
    competitor_dense.push_back(comp_idx);
  }
  const std::uint64_t key =
      ir_outcome_coupling_lookup_key(node_idx, competitor_dense);
  auto it = ctx.ir.outcome_coupling_lookup.find(key);
  if (it == ctx.ir.outcome_coupling_lookup.end()) {
    return false;
  }
  const int spec_idx = it->second;
  if (spec_idx < 0 ||
      spec_idx >= static_cast<int>(ctx.ir.outcome_coupling_ops.size())) {
    return false;
  }
  const uuber::IrOutcomeCouplingOp &spec =
      ctx.ir.outcome_coupling_ops[static_cast<std::size_t>(spec_idx)];
  if (spec.kind != uuber::IrOutcomeCouplingKind::GenericNodeIntegral) {
    return false;
  }
  if (spec.node_idx != node_idx ||
      spec.competitor_count != static_cast<int>(competitor_dense.size())) {
    return false;
  }
  for (int k = 0; k < spec.competitor_count; ++k) {
    const int idx = spec.competitor_begin + k;
    if (idx < 0 || idx >= static_cast<int>(ctx.ir.outcome_competitors.size())) {
      return false;
    }
    if (ctx.ir.outcome_competitors[static_cast<std::size_t>(idx)] !=
        competitor_dense[static_cast<std::size_t>(k)]) {
      return false;
    }
  }
  out = OutcomeCouplingGenericNodeIntegralPayload{};
  out.node_id = ir_dense_idx_to_node_id(ctx, node_idx);
  out.requires_exact_scenario_eval = spec.requires_exact_scenario_eval;
  out.competitor_node_ids.reserve(competitor_dense.size());
  for (int comp_idx : competitor_dense) {
    out.competitor_node_ids.push_back(ir_dense_idx_to_node_id(ctx, comp_idx));
  }
  return true;
}

OutcomeCouplingProgram resolve_outcome_coupling_program_impl(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty,
    bool include_generic_runtime) {
  OutcomeCouplingProgram program;
  if (node_id < 0) {
    return program;
  }
  if (forced_empty && competitor_node_ids.size() == 1 &&
      competitor_node_ids[0] != NA_INTEGER &&
      ir_outcome_coupling_pair_lookup(ctx, node_id, competitor_node_ids[0],
                                      program.pair)) {
    program.kind = OutcomeCouplingOpKind::Pair;
    program.valid = true;
    return program;
  }
  if (forced_empty && competitor_node_ids.size() >= 2 &&
      ir_outcome_coupling_nway_lookup(ctx, node_id, competitor_node_ids,
                                      program.nway)) {
    program.kind = OutcomeCouplingOpKind::NWay;
    program.valid = true;
    return program;
  }
  if (include_generic_runtime && !competitor_node_ids.empty() &&
      ir_outcome_coupling_generic_lookup(ctx, node_id, competitor_node_ids,
                                         program.generic)) {
    program.kind = OutcomeCouplingOpKind::GenericNodeIntegral;
    program.valid = true;
    return program;
  }
  if (include_generic_runtime) {
    program.kind = OutcomeCouplingOpKind::GenericNodeIntegral;
    program.valid = true;
    program.generic.node_id = node_id;
    program.generic.competitor_node_ids = competitor_node_ids;
    program.generic.requires_exact_scenario_eval =
        ir_generic_requires_exact_scenario_eval(ctx, node_id,
                                                competitor_node_ids);
    return program;
  }
  return program;
}

} // namespace

OutcomeCouplingProgram resolve_outcome_coupling_program_with_generic(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty) {
  return resolve_outcome_coupling_program_impl(
      ctx, node_id, competitor_node_ids, forced_empty, true);
}
