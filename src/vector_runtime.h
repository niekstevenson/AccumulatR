#pragma once

#include <cstdint>
#include <vector>

#include "bitset_state.h"
#include "context.h"
#include "time_constraints.h"

namespace uuber {

enum class VectorProgramDomain : std::uint8_t {
  None = 0,
  OutcomeCoupling = 1
};

// Transitional metadata for VT1/VT2 while the runtime still has shape-specific
// evaluator bodies. The end-state is one evaluator over these ops.
enum class VectorProgramPattern : std::uint8_t {
  None = 0,
  CouplingPair = 1,
  CouplingNWay = 2,
  CouplingGenericDirectCdf = 3,
  CouplingGenericIntegral = 4
};

enum class VectorPayloadRole : std::uint8_t {
  None = 0,
  PairX = 1,
  PairY = 2,
  PairC = 3,
  NWayGate = 4,
  NWayTarget = 5,
  NWayCompetitor = 6,
  GenericNode = 7
};

enum class VectorOpCode : std::uint8_t {
  EventDensity = 0,
  EventCDF = 1,
  EventSurvival = 2,
  GenericDirectCDF = 3,
  GenericIntegrand = 4,
  Complement = 5,
  Multiply = 6,
  Add = 7,
  Subtract = 8,
  Integral = 9
};

struct VectorEventRefPayload {
  VectorPayloadRole role{VectorPayloadRole::None};
  LabelRef ref;
  int node_id{-1};
  std::uint32_t node_flags{0u};
};

struct VectorGenericNodePayload {
  VectorPayloadRole role{VectorPayloadRole::GenericNode};
  int node_id{-1};
  std::vector<int> competitor_node_ids;
  bool requires_exact_scenario_eval{false};
};

struct VectorOp {
  VectorOpCode code{VectorOpCode::EventDensity};
  int dst_slot{-1};
  int child_begin{-1};
  int child_count{0};
  int payload_idx{-1};
  double scalar{1.0};
  std::uint32_t flags{0u};
};

struct VectorProgramOutputs {
  int result_slot{-1};
  int coeff_slot{-1};
  std::vector<int> aux_slots;
};

struct VectorProgram {
  VectorProgramDomain domain{VectorProgramDomain::None};
  VectorProgramPattern pattern{VectorProgramPattern::None};
  std::vector<VectorOp> ops;
  std::vector<int> children;
  std::vector<VectorEventRefPayload> event_payloads;
  std::vector<VectorGenericNodePayload> generic_payloads;
  VectorProgramOutputs outputs;
  bool requires_exact_scenario_eval{false};
  bool valid{false};
};

inline const VectorEventRefPayload *
find_vector_event_payload(const VectorProgram &program, VectorPayloadRole role,
                          std::size_t ordinal = 0u) {
  std::size_t seen = 0u;
  for (const auto &payload : program.event_payloads) {
    if (payload.role != role) {
      continue;
    }
    if (seen == ordinal) {
      return &payload;
    }
    ++seen;
  }
  return nullptr;
}

inline std::vector<const VectorEventRefPayload *>
collect_vector_event_payloads(const VectorProgram &program,
                              VectorPayloadRole role) {
  std::vector<const VectorEventRefPayload *> out;
  for (const auto &payload : program.event_payloads) {
    if (payload.role == role) {
      out.push_back(&payload);
    }
  }
  return out;
}

inline const VectorGenericNodePayload *
find_vector_generic_payload(const VectorProgram &program,
                            VectorPayloadRole role =
                                VectorPayloadRole::GenericNode) {
  for (const auto &payload : program.generic_payloads) {
    if (payload.role == role) {
      return &payload;
    }
  }
  return nullptr;
}

struct VectorLaneBatch {
  std::vector<double> t;
  std::vector<std::size_t> density_index;
  std::vector<double> upper;
  std::vector<double> weight;
  std::vector<const TrialParamsSoA *> trial_params_soa;
  std::vector<BitsetState> forced_complete_bits;
  std::vector<BitsetState> forced_survive_bits;
  std::vector<std::uint8_t> forced_complete_bits_valid;
  std::vector<std::uint8_t> forced_survive_bits_valid;
  std::vector<TimeConstraintMap> time_constraints;
  std::vector<int> group_id;
  std::vector<std::uint8_t> active_mask;

  void clear() {
    t.clear();
    density_index.clear();
    upper.clear();
    weight.clear();
    trial_params_soa.clear();
    forced_complete_bits.clear();
    forced_survive_bits.clear();
    forced_complete_bits_valid.clear();
    forced_survive_bits_valid.clear();
    time_constraints.clear();
    group_id.clear();
    active_mask.clear();
  }

  void reserve(std::size_t n) {
    t.reserve(n);
    density_index.reserve(n);
    upper.reserve(n);
    weight.reserve(n);
    trial_params_soa.reserve(n);
    forced_complete_bits.reserve(n);
    forced_survive_bits.reserve(n);
    forced_complete_bits_valid.reserve(n);
    forced_survive_bits_valid.reserve(n);
    time_constraints.reserve(n);
    group_id.reserve(n);
    active_mask.reserve(n);
  }

  std::size_t size() const noexcept { return t.size(); }
  bool empty() const noexcept { return t.empty(); }

  void swap(VectorLaneBatch &other) noexcept {
    t.swap(other.t);
    density_index.swap(other.density_index);
    upper.swap(other.upper);
    weight.swap(other.weight);
    trial_params_soa.swap(other.trial_params_soa);
    forced_complete_bits.swap(other.forced_complete_bits);
    forced_survive_bits.swap(other.forced_survive_bits);
    forced_complete_bits_valid.swap(other.forced_complete_bits_valid);
    forced_survive_bits_valid.swap(other.forced_survive_bits_valid);
    time_constraints.swap(other.time_constraints);
    group_id.swap(other.group_id);
    active_mask.swap(other.active_mask);
  }
};

struct VectorLaneRef {
  double t{0.0};
  std::size_t density_index{0u};
  double weight{0.0};
  const TrialParamsSoA *trial_params_soa{nullptr};
  const BitsetState *forced_complete_bits{nullptr};
  const BitsetState *forced_survive_bits{nullptr};
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  const TimeConstraintMap *time_constraints{nullptr};
};

inline VectorLaneRef vector_lane_ref(const VectorLaneBatch &batch,
                                     std::size_t idx) {
  return VectorLaneRef{
      batch.t[idx],
      batch.density_index[idx],
      batch.weight[idx],
      idx < batch.trial_params_soa.size() ? batch.trial_params_soa[idx]
                                          : nullptr,
      &batch.forced_complete_bits[idx],
      &batch.forced_survive_bits[idx],
      idx < batch.forced_complete_bits_valid.size() &&
          batch.forced_complete_bits_valid[idx] != 0u,
      idx < batch.forced_survive_bits_valid.size() &&
          batch.forced_survive_bits_valid[idx] != 0u,
      &batch.time_constraints[idx]};
}

struct VectorEvalScratch {
  std::vector<double> slot_values;
  std::vector<double> tmp_values;
  std::vector<int> partition_offsets;
  std::vector<int> partition_indices;
  std::vector<std::uint8_t> partition_mask;

  void clear() {
    slot_values.clear();
    tmp_values.clear();
    partition_offsets.clear();
    partition_indices.clear();
    partition_mask.clear();
  }
};

} // namespace uuber
