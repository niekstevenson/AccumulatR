#pragma once

#include <cstdint>

struct UnifiedOutcomeRuntimeStats {
  std::uint64_t evaluate_outcome_coupling_unified_calls{0};
  std::uint64_t kernel_node_density_entry_idx_calls{0};
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
void record_unified_outcome_adaptive_segment_accept();
void record_unified_outcome_adaptive_segment_split();
void record_unified_outcome_generic_labelref_batch_fastpath_call();
void record_unified_outcome_generic_noderef_batch_call();
