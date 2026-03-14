#pragma once

#include <cstdint>

struct UnifiedOutcomeRuntimeStats {
  std::uint64_t evaluate_outcome_coupling_unified_calls{0};
  std::uint64_t kernel_node_density_entry_idx_calls{0};
  std::uint64_t exact_scalar_density_calls{0};
  std::uint64_t exact_batch_density_calls{0};
  std::uint64_t exact_node_batch_calls{0};
  std::uint64_t exact_node_batch_active_points_total{0};
  std::uint64_t exact_node_batch_active_points_max{0};
  std::uint64_t exact_node_batch_recursive_calls{0};
  std::uint64_t exact_node_batch_recursive_active_points_total{0};
  std::uint64_t exact_node_batch_pointwise_fallback_calls{0};
  std::uint64_t exact_node_batch_pointwise_fallback_active_points_total{0};
  std::uint64_t exact_competitor_batch_calls{0};
  std::uint64_t exact_competitor_batch_active_points_total{0};
  std::uint64_t exact_competitor_batch_active_points_max{0};
  std::uint64_t exact_competitor_batch_compiled_ops_total{0};
  std::uint64_t exact_competitor_batch_compiled_ops_max{0};
  std::uint64_t adaptive_segments_accepted{0};
  std::uint64_t adaptive_segments_split{0};
  std::uint64_t generic_labelref_batch_fastpath_calls{0};
  std::uint64_t generic_noderef_batch_calls{0};
  std::uint64_t generic_scalar_fallback_calls{0};
};

extern UnifiedOutcomeRuntimeStats g_unified_outcome_runtime_stats;
extern bool g_unified_outcome_runtime_stats_enabled;

void reset_unified_outcome_runtime_stats();
void record_unified_outcome_coupling_eval_call();
void record_unified_outcome_kernel_density_call();
void record_unified_outcome_exact_scalar_density_call();
void record_unified_outcome_exact_batch_density_call();
void record_unified_outcome_exact_node_batch_call(std::uint64_t active_points,
                                                  bool recursive_call);
void record_unified_outcome_exact_node_batch_pointwise_fallback(
    std::uint64_t active_points);
void record_unified_outcome_exact_competitor_batch_call(
    std::uint64_t active_points, std::uint64_t compiled_op_count);
void record_unified_outcome_adaptive_segment_accept();
void record_unified_outcome_adaptive_segment_split();
void record_unified_outcome_generic_labelref_batch_fastpath_call();
void record_unified_outcome_generic_noderef_batch_call();
