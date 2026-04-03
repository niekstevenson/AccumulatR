#pragma once

#include <cstdint>
#include <vector>

#include "bitset_state.h"
#include "context.h"
#include "time_constraints.h"

namespace uuber {

enum class VectorOpCode : std::uint8_t {
  TreeEvent = 0,
  TreeAnd = 1,
  TreeOr = 2,
  TreeNot = 3,
  TreeGuard = 4
};

struct VectorOp {
  VectorOpCode code{VectorOpCode::TreeEvent};
  int dst_slot{-1};
  int node_idx{-1};
  int event_idx{-1};
  int child_begin{-1};
  int child_count{0};
  int ref_slot{-1};
  int blocker_slot{-1};
  std::uint32_t flags{0u};
};

struct VectorProgramOutputs {
  std::vector<int> node_idx_to_slot;
};

struct VectorProgram {
  std::vector<VectorOp> ops;
  std::vector<int> children;
  VectorProgramOutputs outputs;
  int max_child_count{0};
  bool valid{false};
};

VectorProgram compile_tree_vector_program(const IrContext &ir);

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
