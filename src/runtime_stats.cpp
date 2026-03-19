// [[Rcpp::depends(Rcpp)]]

#include "runtime_stats.h"

#include <algorithm>

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

void record_unified_outcome_exact_batch_density_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_batch_density_calls;
}

void record_unified_outcome_exact_node_batch_call(std::uint64_t active_points,
                                                  bool recursive_call) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_node_batch_calls;
  g_unified_outcome_runtime_stats.exact_node_batch_active_points_total +=
      active_points;
  g_unified_outcome_runtime_stats.exact_node_batch_active_points_max =
      std::max(g_unified_outcome_runtime_stats.exact_node_batch_active_points_max,
               active_points);
  if (recursive_call) {
    ++g_unified_outcome_runtime_stats.exact_node_batch_recursive_calls;
    g_unified_outcome_runtime_stats
        .exact_node_batch_recursive_active_points_total += active_points;
  }
}

void record_unified_outcome_exact_shared_state_partition(
    std::uint64_t active_points, std::uint64_t group_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_shared_state_partition_calls;
  g_unified_outcome_runtime_stats.exact_shared_state_partition_groups_total +=
      group_count;
  g_unified_outcome_runtime_stats
      .exact_shared_state_partition_active_points_total += active_points;
}

void record_unified_outcome_exact_competitor_batch_call(
    std::uint64_t active_points, std::uint64_t compiled_op_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_competitor_batch_calls;
  g_unified_outcome_runtime_stats.exact_competitor_batch_active_points_total +=
      active_points;
  g_unified_outcome_runtime_stats.exact_competitor_batch_active_points_max =
      std::max(g_unified_outcome_runtime_stats
                   .exact_competitor_batch_active_points_max,
               active_points);
  g_unified_outcome_runtime_stats.exact_competitor_batch_compiled_ops_total +=
      compiled_op_count;
  g_unified_outcome_runtime_stats.exact_competitor_batch_compiled_ops_max =
      std::max(g_unified_outcome_runtime_stats
                   .exact_competitor_batch_compiled_ops_max,
               compiled_op_count);
}

void record_unified_outcome_exact_guard_prepared_batch_call(
    std::uint64_t point_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_guard_prepared_batch_calls;
  g_unified_outcome_runtime_stats.exact_guard_prepared_batch_points_total +=
      point_count;
}

void record_unified_outcome_exact_guard_prepared_cache_hit() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_guard_prepared_cache_hits;
}

void record_unified_outcome_exact_guard_prepared_cache_miss() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_guard_prepared_cache_misses;
}

void record_unified_outcome_exact_competitor_guard_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.exact_competitor_guard_fastpath_calls;
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

void record_unified_outcome_shared_trigger_mask_batch_call(
    std::uint64_t mask_count, std::uint64_t point_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.shared_trigger_mask_batch_calls;
  g_unified_outcome_runtime_stats.shared_trigger_mask_batch_masks_total +=
      mask_count;
  g_unified_outcome_runtime_stats.shared_trigger_mask_batch_points_total +=
      point_count;
}

void record_unified_outcome_generic_labelref_batch_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.generic_labelref_batch_fastpath_calls;
}

void record_unified_outcome_generic_poolref_batch_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.generic_poolref_batch_fastpath_calls;
}

void record_unified_outcome_generic_noderef_batch_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.generic_noderef_batch_calls;
}

void record_unified_outcome_direct_batch_spec_attempt() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_batch_spec_attempts;
}

void record_unified_outcome_direct_batch_spec_reject_input() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_batch_spec_reject_input;
}

void record_unified_outcome_direct_batch_spec_reject_contribution() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_batch_spec_reject_contribution;
}

void record_unified_outcome_direct_batch_spec_reject_shape() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_batch_spec_reject_shape;
}

void record_unified_outcome_direct_batch_spec_reject_node() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_batch_spec_reject_node;
}

void record_unified_outcome_ranked_batch_spec_attempt() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_spec_attempts;
}

void record_unified_outcome_ranked_batch_spec_reject_input() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_spec_reject_input;
}

void record_unified_outcome_ranked_batch_spec_reject_contribution() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_spec_reject_contribution;
}

void record_unified_outcome_ranked_batch_spec_reject_shape() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_spec_reject_shape;
}

void record_unified_outcome_ranked_batch_template_cache_hit() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_template_cache_hits;
}

void record_unified_outcome_ranked_batch_template_cache_miss() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.ranked_batch_template_cache_misses;
}

void record_unified_outcome_cpp_loglik_fallback_group_call(
    std::uint64_t trial_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.cpp_loglik_fallback_group_calls;
  g_unified_outcome_runtime_stats.cpp_loglik_fallback_group_trials_total +=
      trial_count;
}

void record_unified_outcome_cpp_loglik_outcome_mass_group_batch_call(
    std::uint64_t trial_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.cpp_loglik_outcome_mass_group_batch_calls;
  g_unified_outcome_runtime_stats
      .cpp_loglik_outcome_mass_group_batch_trials_total += trial_count;
}

void record_unified_outcome_cpp_loglik_outcome_mass_group_batch_exec_failure() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats
       .cpp_loglik_outcome_mass_group_batch_exec_failures;
}

void record_unified_outcome_cpp_loglik_direct_group_batch_call(
    std::uint64_t trial_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.cpp_loglik_direct_group_batch_calls;
  g_unified_outcome_runtime_stats.cpp_loglik_direct_group_batch_trials_total +=
      trial_count;
}

void record_unified_outcome_cpp_loglik_direct_group_batch_exec_failure() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.cpp_loglik_direct_group_batch_exec_failures;
}

void record_unified_outcome_cpp_loglik_direct_shared_trigger_group_batch_call(
    std::uint64_t expanded_point_count, std::uint64_t compressed_point_count) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats
       .cpp_loglik_direct_shared_trigger_group_batch_calls;
  g_unified_outcome_runtime_stats
      .cpp_loglik_direct_shared_trigger_expanded_points_total +=
      expanded_point_count;
  g_unified_outcome_runtime_stats
      .cpp_loglik_direct_shared_trigger_compressed_points_total +=
      compressed_point_count;
}

void record_unified_outcome_kernel_event_batch_call(
    std::uint64_t point_count, bool uses_param_matrix) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.kernel_event_batch_calls;
  g_unified_outcome_runtime_stats.kernel_event_batch_points_total += point_count;
  if (uses_param_matrix) {
    ++g_unified_outcome_runtime_stats.kernel_event_batch_param_matrix_calls;
    g_unified_outcome_runtime_stats
        .kernel_event_batch_param_matrix_points_total += point_count;
  }
}

void record_unified_outcome_direct_node_density_batch_call(
    std::uint64_t point_count, bool uses_param_matrix) {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_node_density_batch_calls;
  g_unified_outcome_runtime_stats.direct_node_density_batch_points_total +=
      point_count;
  if (uses_param_matrix) {
    ++g_unified_outcome_runtime_stats.direct_node_density_batch_param_matrix_calls;
    g_unified_outcome_runtime_stats
        .direct_node_density_batch_param_matrix_points_total += point_count;
  }
}

void record_unified_outcome_direct_node_density_simple_event_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats
       .direct_node_density_batch_simple_event_fastpath_calls;
}

void record_unified_outcome_direct_node_density_simple_competing_fastpath_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats
       .direct_node_density_batch_simple_competing_fastpath_calls;
}

void record_unified_outcome_direct_node_density_kernel_eval_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_node_density_batch_kernel_eval_calls;
}

void record_unified_outcome_direct_node_density_competitor_call() {
  if (!g_unified_outcome_runtime_stats_enabled) {
    return;
  }
  ++g_unified_outcome_runtime_stats.direct_node_density_batch_competitor_calls;
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
      Rcpp::Named("exact_batch_density_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_batch_density_calls),
      Rcpp::Named("exact_node_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_node_batch_calls),
      Rcpp::Named("exact_node_batch_active_points_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_node_batch_active_points_total),
      Rcpp::Named("exact_node_batch_active_points_max") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_node_batch_active_points_max),
      Rcpp::Named("exact_node_batch_recursive_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_node_batch_recursive_calls),
      Rcpp::Named("exact_node_batch_recursive_active_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_node_batch_recursive_active_points_total),
      Rcpp::Named("exact_shared_state_partition_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_shared_state_partition_calls),
      Rcpp::Named("exact_shared_state_partition_groups_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_shared_state_partition_groups_total),
      Rcpp::Named("exact_shared_state_partition_active_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_shared_state_partition_active_points_total),
      Rcpp::Named("exact_competitor_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_competitor_batch_calls),
      Rcpp::Named("exact_competitor_batch_active_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_competitor_batch_active_points_total),
      Rcpp::Named("exact_competitor_batch_active_points_max") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_competitor_batch_active_points_max),
      Rcpp::Named("exact_competitor_batch_compiled_ops_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_competitor_batch_compiled_ops_total),
      Rcpp::Named("exact_competitor_batch_compiled_ops_max") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_competitor_batch_compiled_ops_max),
      Rcpp::Named("exact_guard_prepared_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_guard_prepared_batch_calls),
      Rcpp::Named("exact_guard_prepared_batch_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_guard_prepared_batch_points_total),
      Rcpp::Named("exact_guard_prepared_cache_hits") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_guard_prepared_cache_hits),
      Rcpp::Named("exact_guard_prepared_cache_misses") = static_cast<double>(
          g_unified_outcome_runtime_stats.exact_guard_prepared_cache_misses),
      Rcpp::Named("exact_competitor_guard_fastpath_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .exact_competitor_guard_fastpath_calls),
      Rcpp::Named("adaptive_segments_accepted") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_accepted),
      Rcpp::Named("adaptive_segments_split") = static_cast<double>(
          g_unified_outcome_runtime_stats.adaptive_segments_split),
      Rcpp::Named("shared_trigger_mask_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.shared_trigger_mask_batch_calls),
      Rcpp::Named("shared_trigger_mask_batch_masks_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .shared_trigger_mask_batch_masks_total),
      Rcpp::Named("shared_trigger_mask_batch_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .shared_trigger_mask_batch_points_total),
      Rcpp::Named("generic_labelref_batch_fastpath_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats
              .generic_labelref_batch_fastpath_calls),
      Rcpp::Named("generic_poolref_batch_fastpath_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats
              .generic_poolref_batch_fastpath_calls),
      Rcpp::Named("generic_noderef_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.generic_noderef_batch_calls),
      Rcpp::Named("direct_batch_spec_attempts") = static_cast<double>(
          g_unified_outcome_runtime_stats.direct_batch_spec_attempts),
      Rcpp::Named("direct_batch_spec_reject_input") = static_cast<double>(
          g_unified_outcome_runtime_stats.direct_batch_spec_reject_input),
      Rcpp::Named("direct_batch_spec_reject_contribution") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_batch_spec_reject_contribution),
      Rcpp::Named("direct_batch_spec_reject_shape") = static_cast<double>(
          g_unified_outcome_runtime_stats.direct_batch_spec_reject_shape),
      Rcpp::Named("direct_batch_spec_reject_node") = static_cast<double>(
          g_unified_outcome_runtime_stats.direct_batch_spec_reject_node),
      Rcpp::Named("ranked_batch_spec_attempts") = static_cast<double>(
          g_unified_outcome_runtime_stats.ranked_batch_spec_attempts),
      Rcpp::Named("ranked_batch_spec_reject_input") = static_cast<double>(
          g_unified_outcome_runtime_stats.ranked_batch_spec_reject_input),
      Rcpp::Named("ranked_batch_spec_reject_contribution") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .ranked_batch_spec_reject_contribution),
      Rcpp::Named("ranked_batch_spec_reject_shape") = static_cast<double>(
          g_unified_outcome_runtime_stats.ranked_batch_spec_reject_shape),
      Rcpp::Named("ranked_batch_template_cache_hits") = static_cast<double>(
          g_unified_outcome_runtime_stats.ranked_batch_template_cache_hits),
      Rcpp::Named("ranked_batch_template_cache_misses") = static_cast<double>(
          g_unified_outcome_runtime_stats.ranked_batch_template_cache_misses),
      Rcpp::Named("cpp_loglik_fallback_group_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_fallback_group_calls),
      Rcpp::Named("cpp_loglik_fallback_group_trials_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_fallback_group_trials_total),
      Rcpp::Named("cpp_loglik_outcome_mass_group_batch_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_outcome_mass_group_batch_calls),
      Rcpp::Named("cpp_loglik_outcome_mass_group_batch_trials_total") =
          static_cast<double>(
              g_unified_outcome_runtime_stats
                  .cpp_loglik_outcome_mass_group_batch_trials_total),
      Rcpp::Named("cpp_loglik_outcome_mass_group_batch_exec_failures") =
          static_cast<double>(
              g_unified_outcome_runtime_stats
                  .cpp_loglik_outcome_mass_group_batch_exec_failures),
      Rcpp::Named("cpp_loglik_direct_group_batch_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_group_batch_calls),
      Rcpp::Named("cpp_loglik_direct_group_batch_trials_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_group_batch_trials_total),
      Rcpp::Named("cpp_loglik_direct_group_batch_exec_failures") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_group_batch_exec_failures),
      Rcpp::Named("cpp_loglik_direct_shared_trigger_group_batch_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_shared_trigger_group_batch_calls),
      Rcpp::Named("cpp_loglik_direct_shared_trigger_expanded_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_shared_trigger_expanded_points_total),
      Rcpp::Named("cpp_loglik_direct_shared_trigger_compressed_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .cpp_loglik_direct_shared_trigger_compressed_points_total),
      Rcpp::Named("kernel_event_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.kernel_event_batch_calls),
      Rcpp::Named("kernel_event_batch_points_total") = static_cast<double>(
          g_unified_outcome_runtime_stats.kernel_event_batch_points_total),
      Rcpp::Named("kernel_event_batch_param_matrix_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .kernel_event_batch_param_matrix_calls),
      Rcpp::Named("kernel_event_batch_param_matrix_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .kernel_event_batch_param_matrix_points_total),
      Rcpp::Named("direct_node_density_batch_calls") = static_cast<double>(
          g_unified_outcome_runtime_stats.direct_node_density_batch_calls),
      Rcpp::Named("direct_node_density_batch_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_points_total),
      Rcpp::Named("direct_node_density_batch_param_matrix_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_param_matrix_calls),
      Rcpp::Named("direct_node_density_batch_param_matrix_points_total") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_param_matrix_points_total),
      Rcpp::Named("direct_node_density_batch_simple_event_fastpath_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_simple_event_fastpath_calls),
      Rcpp::Named("direct_node_density_batch_simple_competing_fastpath_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_simple_competing_fastpath_calls),
      Rcpp::Named("direct_node_density_batch_kernel_eval_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_kernel_eval_calls),
      Rcpp::Named("direct_node_density_batch_competitor_calls") =
          static_cast<double>(g_unified_outcome_runtime_stats
                                  .direct_node_density_batch_competitor_calls));
}
