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
    if (ir_node_source_masks_overlap(ctx, target, competitor)) {
      return true;
    }
  }
  return false;
}

inline bool fill_event_payload_from_event(
    const uuber::NativeContext &ctx, int event_idx, uuber::VectorPayloadRole role,
    uuber::VectorEventRefPayload &payload) {
  if (event_idx < 0 || event_idx >= static_cast<int>(ctx.ir.events.size())) {
    return false;
  }
  const uuber::IrEvent &event =
      ctx.ir.events[static_cast<std::size_t>(event_idx)];
  if (event.node_idx < 0 ||
      event.node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return false;
  }
  payload.role = role;
  payload.ref.label_id = event.label_id;
  payload.ref.acc_idx = event.acc_idx;
  payload.ref.pool_idx = event.pool_idx;
  payload.ref.outcome_idx = event.outcome_idx;
  const uuber::IrNode &node =
      ctx.ir.nodes[static_cast<std::size_t>(event.node_idx)];
  payload.node_id = ir_dense_idx_to_node_id(ctx, event.node_idx);
  payload.node_flags = node.flags;
  return true;
}

inline int append_event_payload(uuber::VectorProgram &program,
                                const uuber::VectorEventRefPayload &payload) {
  program.event_payloads.push_back(payload);
  return static_cast<int>(program.event_payloads.size()) - 1;
}

inline int append_generic_payload(uuber::VectorProgram &program, int node_id,
                                  const std::vector<int> &competitor_node_ids,
                                  bool requires_exact_scenario_eval) {
  uuber::VectorGenericNodePayload payload;
  payload.node_id = node_id;
  payload.competitor_node_ids = competitor_node_ids;
  payload.requires_exact_scenario_eval = requires_exact_scenario_eval;
  program.generic_payloads.push_back(std::move(payload));
  return static_cast<int>(program.generic_payloads.size()) - 1;
}

inline int append_children(uuber::VectorProgram &program,
                           const std::vector<int> &child_slots) {
  const int begin = static_cast<int>(program.children.size());
  program.children.insert(program.children.end(), child_slots.begin(),
                          child_slots.end());
  return begin;
}

inline int append_op(uuber::VectorProgram &program, uuber::VectorOpCode code,
                     int payload_idx = -1,
                     const std::vector<int> &child_slots = {},
                     double scalar = 1.0,
                     std::uint32_t flags = 0u) {
  const int dst_slot = static_cast<int>(program.ops.size());
  uuber::VectorOp op;
  op.code = code;
  op.dst_slot = dst_slot;
  op.payload_idx = payload_idx;
  op.scalar = scalar;
  op.flags = flags;
  if (!child_slots.empty()) {
    op.child_begin = append_children(program, child_slots);
    op.child_count = static_cast<int>(child_slots.size());
  }
  program.ops.push_back(op);
  return dst_slot;
}

inline uuber::VectorProgram make_empty_coupling_program() {
  uuber::VectorProgram program;
  program.domain = uuber::VectorProgramDomain::OutcomeCoupling;
  return program;
}

bool lower_pair_program(const uuber::NativeContext &ctx, int node_id,
                        int competitor_node_id, uuber::VectorProgram &out) {
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

  out = make_empty_coupling_program();
  uuber::VectorEventRefPayload x_payload;
  uuber::VectorEventRefPayload y_payload;
  uuber::VectorEventRefPayload c_payload;
  if (!fill_event_payload_from_event(ctx, spec.target_event_idx,
                                     uuber::VectorPayloadRole::PairX,
                                     x_payload) ||
      !fill_event_payload_from_event(ctx, spec.gate_event_idx,
                                     uuber::VectorPayloadRole::PairC,
                                     c_payload)) {
    return false;
  }
  if (spec.aux_event_count != 1 || spec.aux_event_begin < 0 ||
      spec.aux_event_begin >=
          static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
    return false;
  }
  const int y_event_idx = ctx.ir.outcome_coupling_event_indices
      [static_cast<std::size_t>(spec.aux_event_begin)];
  if (!fill_event_payload_from_event(ctx, y_event_idx,
                                     uuber::VectorPayloadRole::PairY,
                                     y_payload)) {
    return false;
  }

  const int x_idx = append_event_payload(out, x_payload);
  const int y_idx = append_event_payload(out, y_payload);
  const int c_idx = append_event_payload(out, c_payload);

  const int c_cdf = append_op(out, uuber::VectorOpCode::EventCDF, c_idx);
  const int c_density =
      append_op(out, uuber::VectorOpCode::EventDensity, c_idx);
  const int x_survival =
      append_op(out, uuber::VectorOpCode::EventSurvival, x_idx);
  const int x_density =
      append_op(out, uuber::VectorOpCode::EventDensity, x_idx);
  const int y_survival =
      append_op(out, uuber::VectorOpCode::EventSurvival, y_idx);
  const int c_survival =
      append_op(out, uuber::VectorOpCode::EventSurvival, c_idx);
  const int fx =
      append_op(out, uuber::VectorOpCode::Complement, -1, {x_survival});
  const int fy =
      append_op(out, uuber::VectorOpCode::Complement, -1, {y_survival});
  const int fc =
      append_op(out, uuber::VectorOpCode::Complement, -1, {c_survival});
  const int term_fc_fx =
      append_op(out, uuber::VectorOpCode::Multiply, -1, {c_density, fx});
  const int term_fx_fc =
      append_op(out, uuber::VectorOpCode::Multiply, -1, {x_density, fc});
  const int term_fx_fy =
      append_op(out, uuber::VectorOpCode::Multiply, -1, {x_density, fy});
  const int int_fc_fx =
      append_op(out, uuber::VectorOpCode::Integral, -1, {term_fc_fx});
  const int int_fx_fc =
      append_op(out, uuber::VectorOpCode::Integral, -1, {term_fx_fc});
  const int int_fx_fy =
      append_op(out, uuber::VectorOpCode::Integral, -1, {term_fx_fy});
  const int coeff_term =
      append_op(out, uuber::VectorOpCode::Multiply, -1, {c_cdf, int_fx_fy});
  const int partial_sum =
      append_op(out, uuber::VectorOpCode::Add, -1, {int_fc_fx, int_fx_fc});
  const int result =
      append_op(out, uuber::VectorOpCode::Subtract, -1,
                {partial_sum, coeff_term});

  out.outputs.coeff_slot = c_cdf;
  out.outputs.aux_slots = {int_fc_fx, int_fx_fc, int_fx_fy};
  out.outputs.result_slot = result;
  out.valid = true;
  return true;
}

bool lower_nway_program(const uuber::NativeContext &ctx, int node_id,
                        const std::vector<int> &competitor_node_ids,
                        uuber::VectorProgram &out) {
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

  out = make_empty_coupling_program();
  uuber::VectorEventRefPayload gate_payload;
  uuber::VectorEventRefPayload target_payload;
  if (!fill_event_payload_from_event(ctx, spec.gate_event_idx,
                                     uuber::VectorPayloadRole::NWayGate,
                                     gate_payload) ||
      !fill_event_payload_from_event(ctx, spec.target_event_idx,
                                     uuber::VectorPayloadRole::NWayTarget,
                                     target_payload)) {
    return false;
  }
  const int gate_idx = append_event_payload(out, gate_payload);
  const int target_idx = append_event_payload(out, target_payload);
  const int gate_cdf = append_op(out, uuber::VectorOpCode::EventCDF, gate_idx);
  const int target_density =
      append_op(out, uuber::VectorOpCode::EventDensity, target_idx);

  int running_product = target_density;
  if (spec.aux_event_count <= 0 || spec.aux_event_begin < 0) {
    return false;
  }
  for (int i = 0; i < spec.aux_event_count; ++i) {
    const int idx = spec.aux_event_begin + i;
    if (idx < 0 ||
        idx >= static_cast<int>(ctx.ir.outcome_coupling_event_indices.size())) {
      return false;
    }
    const int event_idx =
        ctx.ir.outcome_coupling_event_indices[static_cast<std::size_t>(idx)];
    uuber::VectorEventRefPayload comp_payload;
    if (!fill_event_payload_from_event(ctx, event_idx,
                                       uuber::VectorPayloadRole::NWayCompetitor,
                                       comp_payload)) {
      return false;
    }
    const int comp_idx = append_event_payload(out, comp_payload);
    const int comp_survival =
        append_op(out, uuber::VectorOpCode::EventSurvival, comp_idx);
    running_product = append_op(out, uuber::VectorOpCode::Multiply, -1,
                                {running_product, comp_survival});
  }

  const int order_mass =
      append_op(out, uuber::VectorOpCode::Integral, -1, {running_product});
  const int result =
      append_op(out, uuber::VectorOpCode::Multiply, -1, {gate_cdf, order_mass});
  out.outputs.coeff_slot = gate_cdf;
  out.outputs.aux_slots = {order_mass};
  out.outputs.result_slot = result;
  out.valid = true;
  return true;
}

bool lower_generic_program(const uuber::NativeContext &ctx, int node_id,
                           const std::vector<int> &competitor_node_ids,
                           uuber::VectorProgram &out) {
  out = make_empty_coupling_program();
  const bool requires_exact =
      ir_generic_requires_exact_scenario_eval(ctx, node_id, competitor_node_ids);
  const int payload_idx = append_generic_payload(
      out, node_id, competitor_node_ids, requires_exact);
  out.requires_exact_scenario_eval = requires_exact;
  if (competitor_node_ids.empty()) {
    out.outputs.result_slot =
        append_op(out, uuber::VectorOpCode::GenericDirectCDF, payload_idx);
  } else {
    const int term =
        append_op(out, uuber::VectorOpCode::GenericIntegrand, payload_idx);
    out.outputs.aux_slots = {term};
    out.outputs.result_slot =
        append_op(out, uuber::VectorOpCode::Integral, -1, {term});
  }
  out.valid = true;
  return true;
}

uuber::VectorProgram resolve_outcome_coupling_program_impl(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty,
    bool include_generic_runtime) {
  uuber::VectorProgram program = make_empty_coupling_program();
  if (node_id < 0) {
    return program;
  }
  if (forced_empty && competitor_node_ids.size() == 1 &&
      competitor_node_ids[0] != NA_INTEGER &&
      lower_pair_program(ctx, node_id, competitor_node_ids[0], program)) {
    return program;
  }
  if (forced_empty && competitor_node_ids.size() >= 2 &&
      lower_nway_program(ctx, node_id, competitor_node_ids, program)) {
    return program;
  }
  if (include_generic_runtime) {
    lower_generic_program(ctx, node_id, competitor_node_ids, program);
  }
  return program;
}

} // namespace

uuber::VectorProgram resolve_outcome_coupling_program_with_generic(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty) {
  return resolve_outcome_coupling_program_impl(
      ctx, node_id, competitor_node_ids, forced_empty, true);
}
