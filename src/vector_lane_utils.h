#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "forced_state.h"
#include "trial_params.h"
#include "vector_runtime.h"

namespace uuber {

struct VectorRelevantStateStorage {
  TimeConstraintMap time_constraints;
  BitsetState forced_complete_bits;
  BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
};

inline ForcedStateView make_vector_lane_forced_state_view(
    const NativeContext &ctx, const VectorLaneRef &lane) {
  return make_forced_state_view(
      nullptr, lane.forced_complete_bits_valid ? lane.forced_complete_bits
                                               : nullptr,
      &lane.forced_complete_bits_valid,
      lane.forced_survive_bits_valid ? lane.forced_survive_bits : nullptr,
      &lane.forced_survive_bits_valid,
      ctx.ir.label_id_to_bit_idx.empty() ? nullptr : &ctx.ir.label_id_to_bit_idx);
}

inline const TrialParamsSoA *resolve_vector_lane_trial_params_soa(
    const NativeContext &ctx, const VectorLaneRef &lane,
    const TrialParamSet *trial_params,
    const TrialParamsSoA *fallback_trial_params_soa = nullptr) {
  if (lane.trial_params_soa != nullptr) {
    return lane.trial_params_soa;
  }
  if (fallback_trial_params_soa != nullptr) {
    return fallback_trial_params_soa;
  }
  return resolve_trial_params_soa(ctx, trial_params);
}

inline void build_vector_lane_filtered_time_constraints(
    const TimeConstraintMap &source,
    const std::vector<int> &relevant_source_ids,
    TimeConstraintMap &filtered) {
  filtered.clear();
  for (int source_id : relevant_source_ids) {
    auto it = source.find(source_id);
    if (it != source.end()) {
      filtered.emplace(it->first, it->second);
    }
  }
}

inline void build_vector_lane_filtered_forced_bits(
    const NativeContext &ctx, const BitsetState &source_bits,
    bool source_bits_valid, const std::vector<int> &relevant_source_ids,
    BitsetState &filtered_bits, bool &filtered_bits_valid) {
  filtered_bits_valid = false;
  if (!source_bits_valid) {
    return;
  }
  (void)ensure_forced_bitset_capacity(ctx, filtered_bits, filtered_bits_valid);
  filtered_bits.clear_all();
  for (int source_id : relevant_source_ids) {
    if (forced_bits_contains_label_id_strict(ctx, source_id, source_bits,
                                             source_bits_valid)) {
      set_forced_id_bit_strict(ctx, source_id, filtered_bits,
                               filtered_bits_valid);
    }
  }
}

template <typename PointVec>
inline bool vector_lanes_match_relevant_state(
    const PointVec &points, std::size_t lhs_idx, std::size_t rhs_idx,
    const std::vector<int> &relevant_source_ids, const NativeContext &ctx) {
  const VectorLaneRef lhs = vector_lane_ref(points, lhs_idx);
  const VectorLaneRef rhs = vector_lane_ref(points, rhs_idx);
  for (int source_id : relevant_source_ids) {
    if (forced_bits_contains_label_id_strict(
            ctx, source_id, *lhs.forced_complete_bits,
            lhs.forced_complete_bits_valid) !=
            forced_bits_contains_label_id_strict(
                ctx, source_id, *rhs.forced_complete_bits,
                rhs.forced_complete_bits_valid) ||
        forced_bits_contains_label_id_strict(
            ctx, source_id, *lhs.forced_survive_bits,
            lhs.forced_survive_bits_valid) !=
            forced_bits_contains_label_id_strict(
                ctx, source_id, *rhs.forced_survive_bits,
                rhs.forced_survive_bits_valid)) {
      return false;
    }
    const SourceTimeConstraint *lhs_constraint =
        time_constraints_find(lhs.time_constraints, source_id);
    const SourceTimeConstraint *rhs_constraint =
        time_constraints_find(rhs.time_constraints, source_id);
    const SourceTimeConstraint lhs_value =
        lhs_constraint ? *lhs_constraint : SourceTimeConstraint{};
    const SourceTimeConstraint rhs_value =
        rhs_constraint ? *rhs_constraint : SourceTimeConstraint{};
    if (!time_constraint_equal(lhs_value, rhs_value)) {
      return false;
    }
  }
  return true;
}

template <typename PointVec>
inline bool vector_lanes_share_relevant_state(
    const PointVec &points, const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const NativeContext &ctx,
    std::size_t &anchor_idx) {
  anchor_idx = points.size();
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    anchor_idx = i;
    break;
  }
  if (anchor_idx >= points.size()) {
    return false;
  }
  for (std::size_t i = anchor_idx + 1u; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    if (!vector_lanes_match_relevant_state(points, anchor_idx, i,
                                           relevant_source_ids, ctx)) {
      return false;
    }
  }
  return true;
}

template <typename PointVec>
inline bool prepare_vector_lane_shared_state_storage(
    const NativeContext &ctx, const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, std::size_t &anchor_idx,
    VectorRelevantStateStorage &storage) {
  anchor_idx = points.size();
  if (!vector_lanes_share_relevant_state(points, active_mask, relevant_source_ids,
                                         ctx, anchor_idx)) {
    return false;
  }
  const VectorLaneRef anchor = vector_lane_ref(points, anchor_idx);
  build_vector_lane_filtered_time_constraints(*anchor.time_constraints,
                                              relevant_source_ids,
                                              storage.time_constraints);
  build_vector_lane_filtered_forced_bits(
      ctx, *anchor.forced_complete_bits, anchor.forced_complete_bits_valid,
      relevant_source_ids, storage.forced_complete_bits,
      storage.forced_complete_bits_valid);
  build_vector_lane_filtered_forced_bits(
      ctx, *anchor.forced_survive_bits, anchor.forced_survive_bits_valid,
      relevant_source_ids, storage.forced_survive_bits,
      storage.forced_survive_bits_valid);
  return true;
}

template <typename PointVec>
inline bool partition_vector_lanes_by_relevant_state(
    const PointVec &points, const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, const NativeContext &ctx,
    std::vector<std::vector<std::uint8_t>> &group_masks) {
  group_masks.clear();
  std::vector<std::size_t> anchor_indices;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (active_mask && (i >= active_mask->size() || (*active_mask)[i] == 0u)) {
      continue;
    }
    bool assigned = false;
    for (std::size_t group_idx = 0; group_idx < anchor_indices.size();
         ++group_idx) {
      if (vector_lanes_match_relevant_state(points, anchor_indices[group_idx], i,
                                            relevant_source_ids, ctx)) {
        group_masks[group_idx][i] = 1u;
        assigned = true;
        break;
      }
    }
    if (assigned) {
      continue;
    }
    anchor_indices.push_back(i);
    group_masks.emplace_back(points.size(), 0u);
    group_masks.back()[i] = 1u;
  }
  return !group_masks.empty();
}

inline bool vector_lane_active(const std::vector<std::uint8_t> *active_mask,
                               std::size_t idx) {
  return !active_mask ||
         (idx < active_mask->size() && (*active_mask)[idx] != 0u);
}

inline std::size_t count_active_vector_lanes(
    const std::vector<std::uint8_t> *active_mask, std::size_t point_count) {
  if (active_mask == nullptr) {
    return point_count;
  }
  std::size_t active_count = 0u;
  const std::size_t limit = std::min(point_count, active_mask->size());
  for (std::size_t i = 0; i < limit; ++i) {
    active_count += vector_lane_active(active_mask, i) ? 1u : 0u;
  }
  return active_count;
}

template <typename PointVec, typename EvalGroupFn, typename FallbackFn,
          typename RecordPartitionFn>
inline bool eval_vector_lanes_with_shared_state(
    const NativeContext &ctx, const PointVec &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<int> &relevant_source_ids, EvalGroupFn &&eval_group,
    FallbackFn &&fallback_eval, RecordPartitionFn &&record_partition,
    bool fallback_on_singleton_partition = false) {
  const std::size_t active_point_count =
      count_active_vector_lanes(active_mask, points.size());
  if (active_point_count == 0u) {
    return true;
  }

  std::size_t anchor_idx = points.size();
  VectorRelevantStateStorage storage;
  if (!prepare_vector_lane_shared_state_storage(
          ctx, points, active_mask, relevant_source_ids, anchor_idx, storage)) {
    std::vector<std::vector<std::uint8_t>> group_masks;
    if (!partition_vector_lanes_by_relevant_state(
            points, active_mask, relevant_source_ids, ctx, group_masks) ||
        group_masks.size() <= 1u) {
      return std::forward<FallbackFn>(fallback_eval)();
    }
    if (fallback_on_singleton_partition &&
        group_masks.size() >= active_point_count) {
      return std::forward<FallbackFn>(fallback_eval)();
    }
    std::forward<RecordPartitionFn>(record_partition)(
        static_cast<std::uint64_t>(active_point_count),
        static_cast<std::uint64_t>(group_masks.size()));
    for (const std::vector<std::uint8_t> &group_mask : group_masks) {
      std::size_t group_anchor_idx = points.size();
      VectorRelevantStateStorage group_storage;
      if (!prepare_vector_lane_shared_state_storage(
              ctx, points, &group_mask, relevant_source_ids, group_anchor_idx,
              group_storage) ||
          !std::forward<EvalGroupFn>(eval_group)(group_anchor_idx,
                                                 group_storage, &group_mask)) {
        return std::forward<FallbackFn>(fallback_eval)();
      }
    }
    return true;
  }
  return std::forward<EvalGroupFn>(eval_group)(anchor_idx, storage, active_mask);
}

} // namespace uuber
