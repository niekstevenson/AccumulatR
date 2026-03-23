#pragma once

#include <cstdint>

struct UnifiedOutcomeRuntimeStats {
  std::uint64_t evaluate_outcome_coupling_unified_calls{0};
  std::uint64_t exact_node_batch_calls{0};
  std::uint64_t exact_node_batch_active_points_total{0};
  std::uint64_t exact_node_batch_active_points_max{0};
  std::uint64_t exact_node_batch_recursive_calls{0};
  std::uint64_t exact_node_batch_recursive_active_points_total{0};
  std::uint64_t exact_shared_state_partition_calls{0};
  std::uint64_t exact_shared_state_partition_groups_total{0};
  std::uint64_t exact_shared_state_partition_active_points_total{0};
  std::uint64_t exact_competitor_batch_calls{0};
  std::uint64_t exact_competitor_batch_active_points_total{0};
  std::uint64_t exact_competitor_batch_active_points_max{0};
  std::uint64_t exact_competitor_batch_compiled_ops_total{0};
  std::uint64_t exact_competitor_batch_compiled_ops_max{0};
  std::uint64_t exact_guard_prepared_batch_calls{0};
  std::uint64_t exact_guard_prepared_batch_points_total{0};
  std::uint64_t exact_guard_prepared_cache_hits{0};
  std::uint64_t exact_guard_prepared_cache_misses{0};
  std::uint64_t exact_competitor_guard_fastpath_calls{0};
  std::uint64_t adaptive_segments_accepted{0};
  std::uint64_t adaptive_segments_split{0};
  std::uint64_t shared_trigger_mask_batch_calls{0};
  std::uint64_t shared_trigger_mask_batch_masks_total{0};
  std::uint64_t shared_trigger_mask_batch_points_total{0};
  std::uint64_t generic_labelref_batch_fastpath_calls{0};
  std::uint64_t generic_poolref_batch_fastpath_calls{0};
  std::uint64_t generic_noderef_batch_calls{0};
  std::uint64_t direct_batch_spec_attempts{0};
  std::uint64_t direct_batch_spec_reject_input{0};
  std::uint64_t direct_batch_spec_reject_contribution{0};
  std::uint64_t direct_batch_spec_reject_shape{0};
  std::uint64_t direct_batch_spec_reject_node{0};
  std::uint64_t kernel_event_batch_calls{0};
  std::uint64_t kernel_event_batch_points_total{0};
  std::uint64_t kernel_event_batch_param_matrix_calls{0};
  std::uint64_t kernel_event_batch_param_matrix_points_total{0};
  std::uint64_t direct_node_density_batch_calls{0};
  std::uint64_t direct_node_density_batch_points_total{0};
  std::uint64_t direct_node_density_batch_param_matrix_calls{0};
  std::uint64_t direct_node_density_batch_param_matrix_points_total{0};
  std::uint64_t direct_node_density_batch_simple_event_fastpath_calls{0};
  std::uint64_t direct_node_density_batch_simple_competing_fastpath_calls{0};
  std::uint64_t direct_node_density_batch_kernel_eval_calls{0};
  std::uint64_t direct_node_density_batch_competitor_calls{0};
};

extern UnifiedOutcomeRuntimeStats g_unified_outcome_runtime_stats;
extern bool g_unified_outcome_runtime_stats_enabled;

void reset_unified_outcome_runtime_stats();
void record_unified_outcome_coupling_eval_call();
void record_unified_outcome_exact_node_batch_call(std::uint64_t active_points,
                                                  bool recursive_call);
void record_unified_outcome_exact_shared_state_partition(
    std::uint64_t active_points, std::uint64_t group_count);
void record_unified_outcome_exact_competitor_batch_call(
    std::uint64_t active_points, std::uint64_t compiled_op_count);
void record_unified_outcome_exact_guard_prepared_batch_call(
    std::uint64_t point_count);
void record_unified_outcome_exact_guard_prepared_cache_hit();
void record_unified_outcome_exact_guard_prepared_cache_miss();
void record_unified_outcome_exact_competitor_guard_fastpath_call();
void record_unified_outcome_adaptive_segment_accept();
void record_unified_outcome_adaptive_segment_split();
void record_unified_outcome_shared_trigger_mask_batch_call(
    std::uint64_t mask_count, std::uint64_t point_count);
void record_unified_outcome_generic_labelref_batch_fastpath_call();
void record_unified_outcome_generic_poolref_batch_fastpath_call();
void record_unified_outcome_generic_noderef_batch_call();
void record_unified_outcome_direct_batch_spec_attempt();
void record_unified_outcome_direct_batch_spec_reject_input();
void record_unified_outcome_direct_batch_spec_reject_contribution();
void record_unified_outcome_direct_batch_spec_reject_shape();
void record_unified_outcome_direct_batch_spec_reject_node();
void record_unified_outcome_kernel_event_batch_call(
    std::uint64_t point_count, bool uses_param_matrix);
void record_unified_outcome_direct_node_density_batch_call(
    std::uint64_t point_count, bool uses_param_matrix);
void record_unified_outcome_direct_node_density_simple_event_fastpath_call();
void record_unified_outcome_direct_node_density_simple_competing_fastpath_call();
void record_unified_outcome_direct_node_density_kernel_eval_call();
void record_unified_outcome_direct_node_density_competitor_call();
