// [[Rcpp::depends(Rcpp)]]

#include "runtime_stats.h"

#include <Rcpp.h>

UnifiedOutcomeRuntimeStats g_unified_outcome_runtime_stats;
bool g_unified_outcome_runtime_stats_enabled{false};

void reset_unified_outcome_runtime_stats() {
  g_unified_outcome_runtime_stats = UnifiedOutcomeRuntimeStats{};
}

void record_unified_outcome_coupling_eval_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.evaluate_outcome_coupling_unified_calls;
}

void record_unified_outcome_kernel_density_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.kernel_node_density_entry_idx_calls;
}

void record_unified_outcome_exact_scalar_density_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_scalar_density_calls;
}

void record_unified_outcome_exact_batch_density_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_batch_density_calls;
}

void record_unified_outcome_adaptive_segment_accept() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.adaptive_segments_accepted;
}

void record_unified_outcome_adaptive_segment_split() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.adaptive_segments_split;
}

void record_unified_outcome_generic_labelref_batch_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.generic_labelref_batch_fastpath_calls;
}

void record_unified_outcome_generic_noderef_batch_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.generic_noderef_batch_calls;
}

// [[Rcpp::export]]
void unified_outcome_stats_reset_cpp() {
  g_unified_outcome_runtime_stats_enabled = true;
  reset_unified_outcome_runtime_stats();
}

// [[Rcpp::export]]
Rcpp::List unified_outcome_stats_cpp() {
  return Rcpp::List::create(
      Rcpp::Named("evaluate_outcome_coupling_unified_calls") =
          static_cast<double>(
              g_unified_outcome_runtime_stats
                  .evaluate_outcome_coupling_unified_calls),
      Rcpp::Named("kernel_node_density_entry_idx_calls") =
          static_cast<double>(
              g_unified_outcome_runtime_stats
                  .kernel_node_density_entry_idx_calls),
      Rcpp::Named("exact_scalar_density_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_scalar_density_calls),
      Rcpp::Named("exact_batch_density_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_batch_density_calls),
      Rcpp::Named("adaptive_segments_accepted") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_accepted),
      Rcpp::Named("adaptive_segments_split") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_split),
      Rcpp::Named("generic_labelref_batch_fastpath_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats
              .generic_labelref_batch_fastpath_calls),
      Rcpp::Named("generic_noderef_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.generic_noderef_batch_calls),
      Rcpp::Named("generic_scalar_fallback_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.generic_scalar_fallback_calls));
}
