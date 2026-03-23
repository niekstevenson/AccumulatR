#pragma once

#include <cstdint>

struct UnifiedOutcomeRuntimeStats {
  std::uint64_t coupling_eval_calls{0};
  std::uint64_t node_batch_calls{0};
  std::uint64_t node_batch_active_lanes_total{0};
  std::uint64_t node_batch_active_lanes_max{0};
  std::uint64_t node_batch_recursive_calls{0};
  std::uint64_t node_batch_recursive_active_lanes_total{0};
  std::uint64_t lane_partition_calls{0};
  std::uint64_t lane_partition_groups_total{0};
  std::uint64_t lane_partition_active_lanes_total{0};
  std::uint64_t leaf_batch_calls{0};
  std::uint64_t leaf_batch_lanes_total{0};
  std::uint64_t leaf_batch_param_matrix_calls{0};
  std::uint64_t leaf_batch_param_matrix_lanes_total{0};
  std::uint64_t adaptive_segments_accepted{0};
  std::uint64_t adaptive_segments_split{0};
};

extern UnifiedOutcomeRuntimeStats g_unified_outcome_runtime_stats;
extern bool g_unified_outcome_runtime_stats_enabled;

void reset_unified_outcome_runtime_stats();
void record_unified_outcome_coupling_eval_call();
void record_unified_outcome_node_batch_call(std::uint64_t active_lanes,
                                            bool recursive_call);
void record_unified_outcome_lane_partition(
    std::uint64_t active_lanes, std::uint64_t group_count);
void record_unified_outcome_leaf_batch_call(std::uint64_t lane_count,
                                            bool uses_param_matrix);
void record_unified_outcome_adaptive_segment_accept();
void record_unified_outcome_adaptive_segment_split();
