#pragma once

#include <cstdint>
#include <vector>

#include "context.h"

struct OutcomeCouplingPairPayload {
  uuber::LabelRef x_ref;
  uuber::LabelRef y_ref;
  uuber::LabelRef c_ref;
};

struct OutcomeCouplingNWayPayload {
  uuber::LabelRef gate_ref;
  uuber::LabelRef target_ref;
  std::vector<uuber::LabelRef> competitor_refs;
};

struct OutcomeCouplingGenericNodeIntegralPayload {
  int node_id{-1};
  std::vector<int> competitor_node_ids;
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
