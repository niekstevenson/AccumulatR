#pragma once

#include <cstdint>
#include <vector>

#include "context.h"

struct OutcomeCouplingEventRefPayload {
  uuber::LabelRef ref;
  int node_id{-1};
  std::uint32_t node_flags{0u};
};

struct OutcomeCouplingPairPayload {
  OutcomeCouplingEventRefPayload x_ref;
  OutcomeCouplingEventRefPayload y_ref;
  OutcomeCouplingEventRefPayload c_ref;
};

struct OutcomeCouplingNWayPayload {
  OutcomeCouplingEventRefPayload gate_ref;
  OutcomeCouplingEventRefPayload target_ref;
  std::vector<OutcomeCouplingEventRefPayload> competitor_refs;
};

struct OutcomeCouplingGenericNodeIntegralPayload {
  int node_id{-1};
  std::vector<int> competitor_node_ids;
  bool requires_exact_scenario_eval{false};
};

enum class OutcomeCouplingOpKind : std::uint8_t {
  None = 0,
  Pair = 1,
  NWay = 2,
  GenericNodeIntegral = 3
};

struct OutcomeCouplingProgram {
  OutcomeCouplingOpKind kind{OutcomeCouplingOpKind::None};
  OutcomeCouplingPairPayload pair{};
  OutcomeCouplingNWayPayload nway{};
  OutcomeCouplingGenericNodeIntegralPayload generic{};
  bool valid{false};
};

OutcomeCouplingProgram resolve_outcome_coupling_program_with_generic(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty);
