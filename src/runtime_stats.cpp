// [[Rcpp::depends(Rcpp)]]

#include "runtime_stats.h"

#include <algorithm>

#include <Rcpp.h>

UnifiedOutcomeRuntimeStats g_unified_outcome_runtime_stats;
bool g_unified_outcome_runtime_stats_enabled{false};

void reset_unified_outcome_runtime_stats() {
  g_unified_outcome_runtime_stats = UnifiedOutcomeRuntimeStats{};
}

void record_unified_outcome_node_batch_call(std::uint64_t active_lanes,
                                            bool recursive_call) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.node_batch_calls;
  g_unified_outcome_runtime_stats.node_batch_active_lanes_total += active_lanes;
  g_unified_outcome_runtime_stats.node_batch_active_lanes_max =
      std::max(g_unified_outcome_runtime_stats.node_batch_active_lanes_max,
               active_lanes);
  if (recursive_call) {
    ++g_unified_outcome_runtime_stats.node_batch_recursive_calls;
    g_unified_outcome_runtime_stats.node_batch_recursive_active_lanes_total +=
        active_lanes;
  }
}

void record_unified_outcome_lane_partition(std::uint64_t active_lanes,
                                           std::uint64_t group_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.lane_partition_calls;
  g_unified_outcome_runtime_stats.lane_partition_groups_total += group_count;
  g_unified_outcome_runtime_stats.lane_partition_active_lanes_total +=
      active_lanes;
}

void record_unified_outcome_leaf_batch_call(std::uint64_t lane_count,
                                            bool uses_param_matrix) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.leaf_batch_calls;
  g_unified_outcome_runtime_stats.leaf_batch_lanes_total += lane_count;
  if (uses_param_matrix) {
    ++g_unified_outcome_runtime_stats.leaf_batch_param_matrix_calls;
    g_unified_outcome_runtime_stats.leaf_batch_param_matrix_lanes_total +=
        lane_count;
  }
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

// [[Rcpp::export]]
void unified_outcome_stats_reset_cpp() {
  g_unified_outcome_runtime_stats_enabled = true;
  reset_unified_outcome_runtime_stats();
}

// [[Rcpp::export]]
Rcpp::List unified_outcome_stats_cpp() {
  return Rcpp::List::create(
      Rcpp::Named("node_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.node_batch_calls),
      Rcpp::Named("node_batch_active_lanes_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.node_batch_active_lanes_total),
      Rcpp::Named("node_batch_active_lanes_max") = static_cast<double>(
          g_unified_outcome_runtime_stats.node_batch_active_lanes_max),
      Rcpp::Named("node_batch_recursive_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.node_batch_recursive_calls),
      Rcpp::Named("node_batch_recursive_active_lanes_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .node_batch_recursive_active_lanes_total),
      Rcpp::Named("lane_partition_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.lane_partition_calls),
      Rcpp::Named("lane_partition_groups_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.lane_partition_groups_total),
      Rcpp::Named("lane_partition_active_lanes_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.lane_partition_active_lanes_total),
      Rcpp::Named("leaf_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.leaf_batch_calls),
      Rcpp::Named("leaf_batch_lanes_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.leaf_batch_lanes_total),
      Rcpp::Named("leaf_batch_param_matrix_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.leaf_batch_param_matrix_calls),
      Rcpp::Named("leaf_batch_param_matrix_lanes_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.leaf_batch_param_matrix_lanes_total),
      Rcpp::Named("adaptive_segments_accepted") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_accepted),
      Rcpp::Named("adaptive_segments_split") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_split));
}
