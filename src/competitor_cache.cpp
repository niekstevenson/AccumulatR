#include "competitor_cache.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <utility>
#include <vector>

#include "evaluator_internal.h"

namespace {

struct CompetitorMeta {
  int node_id{-1};
  int node_idx{-1};
  int guard_transition_idx{-1};
  int transition_mask_begin{-1};
  int transition_mask_count{0};
  int invalidate_slot{0};
  bool transition_required{false};
  std::vector<int> sources;
  bool scenario_sensitive{false};
};

std::uint64_t competitor_cache_key_hash(const std::vector<int> &ids) {
  std::uint64_t hash = kFNV64Offset;
  const std::int32_t count = static_cast<std::int32_t>(ids.size());
  hash_append_bytes(hash, &count, sizeof(count));
  for (int id : ids) {
    const std::int32_t value = static_cast<std::int32_t>(id);
    hash_append_bytes(hash, &value, sizeof(value));
  }
  return hash;
}

inline void validate_competitor_transition_meta(const uuber::NativeContext &ctx,
                                                const CompetitorMeta &meta) {
  if (meta.node_idx < 0) {
    Rcpp::stop("IR competitor node index missing for node %d", meta.node_id);
  }
  if (!meta.transition_required) {
    return;
  }
  if (meta.transition_mask_begin < 0 || meta.transition_mask_count <= 0 ||
      meta.transition_mask_begin + meta.transition_mask_count >
          static_cast<int>(ctx.ir.node_source_masks.size())) {
    Rcpp::stop("IR competitor transition mask invalid for node %d",
               meta.node_id);
  }
  if (ctx.kernel_program.valid &&
      (meta.invalidate_slot < 0 ||
       meta.invalidate_slot >= static_cast<int>(ctx.kernel_program.ops.size()))) {
    Rcpp::stop(
        "IR competitor transition invalidate slot invalid for node %d",
        meta.node_id);
  }
  if (meta.guard_transition_idx >= 0 &&
      meta.guard_transition_idx >=
          static_cast<int>(ctx.kernel_state_graph.guard_transitions.size())) {
    Rcpp::stop("IR competitor guard-transition index invalid for node %d",
               meta.node_id);
  }
}

bool share_sources(const std::vector<int> &a, const std::vector<int> &b) {
  if (a.empty() || b.empty()) {
    return false;
  }
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j]) {
      return true;
    }
    if (a[i] < b[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return false;
}

std::vector<std::vector<int>>
build_competitor_clusters(const std::vector<CompetitorMeta> &metas) {
  std::vector<std::vector<int>> clusters;
  const std::size_t n = metas.size();
  if (n == 0) {
    return clusters;
  }
  std::vector<bool> visited(n, false);
  for (std::size_t i = 0; i < n; ++i) {
    if (visited[i]) {
      continue;
    }
    std::vector<int> cluster;
    std::queue<std::size_t> q;
    q.push(i);
    visited[i] = true;
    while (!q.empty()) {
      std::size_t idx = q.front();
      q.pop();
      cluster.push_back(static_cast<int>(idx));
      for (std::size_t j = 0; j < n; ++j) {
        if (visited[j]) {
          continue;
        }
        if (share_sources(metas[idx].sources, metas[j].sources)) {
          visited[j] = true;
          q.push(j);
        }
      }
    }
    clusters.push_back(std::move(cluster));
  }
  return clusters;
}

inline bool competitor_kernel_runtime_usable(const NodeEvalState &state) {
  return state.kernel_runtime_ready && state.kernel_runtime_ptr &&
         kernel_runtime_cache_safe(state);
}

double run_competitor_batch_survival_op(
    const uuber::CompetitorCompiledOp &op, NodeEvalState &state,
    const uuber::KernelEventEvalFn &event_eval_cb,
    const uuber::KernelGuardEvalFn &guard_eval_cb) {
  if (op.target_node_indices.empty()) {
    return 1.0;
  }

  uuber::KernelEvalNeed kernel_need;
  kernel_need.density = false;
  kernel_need.survival = true;
  kernel_need.cdf = true;
  std::vector<uuber::KernelNodeValues> kernel_values;
  const bool use_state_runtime = competitor_kernel_runtime_usable(state);
  thread_local uuber::KernelRuntimeState fallback_runtime;
  uuber::KernelRuntimeState *runtime_ptr =
      use_state_runtime ? state.kernel_runtime_ptr : &fallback_runtime;
  if (!use_state_runtime) {
    uuber::invalidate_kernel_runtime_from_slot(*runtime_ptr, 0);
  }
  if (!uuber::eval_kernel_nodes_incremental(
          state.ctx.kernel_program, *runtime_ptr, op.target_node_indices,
          kernel_need, event_eval_cb, guard_eval_cb, kernel_values) ||
      kernel_values.size() != op.target_node_indices.size()) {
    Rcpp::stop("IR kernel execution failed for competitor op");
  }
  if (!use_state_runtime) {
    uuber::invalidate_kernel_runtime_from_slot(*runtime_ptr, 0);
  }

  double product = 1.0;
  for (const auto &vals : kernel_values) {
    const double surv = clamp_probability(vals.survival);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    product *= surv;
    if (!std::isfinite(product) || product <= 0.0) {
      return 0.0;
    }
  }
  if (op.transition_mask_begin >= 0 && op.transition_mask_count > 0) {
    apply_transition_mask_words(
        state.ctx, op.transition_mask_begin, op.transition_mask_count,
        op.transition_invalidate_slot, state.forced_survive_bits,
        state.forced_survive_bits_valid,
        (state.kernel_runtime_ready && state.kernel_runtime_ptr)
            ? state.kernel_runtime_ptr
            : nullptr);
  }
  return clamp_probability(product);
}

double run_competitor_compiled_ops_product(
    const std::vector<uuber::CompetitorCompiledOp> &compiled_ops,
    NodeEvalState &state, const uuber::KernelEventEvalFn &event_eval_cb,
    const uuber::KernelGuardEvalFn &guard_eval_cb) {
  if (compiled_ops.empty()) {
    return 1.0;
  }
  double product = 1.0;
  for (const uuber::CompetitorCompiledOp &op : compiled_ops) {
    const double surv =
        run_competitor_batch_survival_op(op, state, event_eval_cb, guard_eval_cb);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    product *= surv;
    if (!std::isfinite(product) || product <= 0.0) {
      return 0.0;
    }
  }
  return clamp_probability(product);
}

} // namespace

const uuber::CompetitorClusterCacheEntry &
fetch_competitor_cluster_cache(const uuber::NativeContext &ctx,
                               const std::vector<int> &competitor_ids) {
  uuber::CompetitorCacheMap &cache = ctx.competitor_cache;
  const std::uint64_t key_hash = competitor_cache_key_hash(competitor_ids);
  auto it = cache.find(key_hash);
  if (it != cache.end()) {
    for (const uuber::CompetitorCacheRecord &record : it->second) {
      if (record.competitor_ids == competitor_ids) {
        return record.entry;
      }
    }
  }

  uuber::CompetitorClusterCacheEntry entry;
  std::vector<CompetitorMeta> metas;
  metas.reserve(competitor_ids.size());
  if (ctx.kernel_state_graph.node_contains_guard.size() != ctx.ir.nodes.size() ||
      ctx.kernel_state_graph.node_competitor_guard_transition_idx.size() !=
          ctx.ir.nodes.size() ||
      ctx.kernel_state_graph.node_competitor_transition_mask_begin.size() !=
          ctx.ir.nodes.size() ||
      ctx.kernel_state_graph.node_competitor_transition_mask_count.size() !=
          ctx.ir.nodes.size() ||
      ctx.kernel_state_graph.node_competitor_transition_invalidate_slot.size() !=
          ctx.ir.nodes.size()) {
    Rcpp::stop("IR competitor transition metadata unavailable");
  }
  for (int node_id : competitor_ids) {
    int node_idx = resolve_dense_node_idx_required(ctx, node_id);
    const uuber::IrNode &node = ir_node_required(ctx, node_idx);
    CompetitorMeta meta;
    meta.node_id = node_id;
    meta.node_idx = node_idx;
    meta.sources = ensure_source_ids(ctx, node);
    sort_unique(meta.sources);
    const bool contains_guard =
        ctx.kernel_state_graph
            .node_contains_guard[static_cast<std::size_t>(node_idx)] != 0u;
    meta.scenario_sensitive =
        (node.flags & uuber::IR_NODE_FLAG_SCENARIO_SENSITIVE) != 0u;
    meta.transition_required = contains_guard && !meta.sources.empty();
    if (meta.transition_required) {
      meta.guard_transition_idx = ctx.kernel_state_graph
          .node_competitor_guard_transition_idx[static_cast<std::size_t>(
              node_idx)];
      meta.transition_mask_begin = ctx.kernel_state_graph
          .node_competitor_transition_mask_begin[static_cast<std::size_t>(
              node_idx)];
      meta.transition_mask_count = ctx.kernel_state_graph
          .node_competitor_transition_mask_count[static_cast<std::size_t>(
              node_idx)];
      meta.invalidate_slot = ctx.kernel_state_graph
          .node_competitor_transition_invalidate_slot[static_cast<std::size_t>(
              node_idx)];
    }
    validate_competitor_transition_meta(ctx, meta);
    metas.push_back(std::move(meta));
  }

  const std::vector<std::vector<int>> clusters =
      build_competitor_clusters(metas);
  std::size_t op_estimate = 0;
  for (const std::vector<int> &cluster : clusters) {
    op_estimate += cluster.size();
  }
  entry.compiled_ops.reserve(op_estimate == 0 ? clusters.size() : op_estimate);
  for (const std::vector<int> &cluster : clusters) {
    bool requires_transition = false;
    for (int idx : cluster) {
      if (metas[static_cast<std::size_t>(idx)].transition_required) {
        requires_transition = true;
        break;
      }
    }
    if (!requires_transition) {
      uuber::CompetitorCompiledOp op;
      op.transition_guard_idx = -1;
      op.transition_mask_begin = -1;
      op.transition_mask_count = 0;
      op.transition_invalidate_slot = 0;
      op.target_node_indices.reserve(cluster.size());
      for (int idx : cluster) {
        if (idx < 0 || idx >= static_cast<int>(metas.size())) {
          Rcpp::stop("IR competitor batch index out of range: %d", idx);
        }
        const CompetitorMeta &meta = metas[static_cast<std::size_t>(idx)];
        if (meta.node_idx < 0) {
          Rcpp::stop("IR competitor node index missing for node %d",
                     meta.node_id);
        }
        op.target_node_indices.push_back(meta.node_idx);
      }
      entry.compiled_ops.push_back(std::move(op));
      continue;
    }

    struct OrderingMeta {
      int index;
      std::size_t source_count;
      bool scenario_sensitive;
    };
    std::vector<OrderingMeta> order;
    order.reserve(cluster.size());
    for (int idx : cluster) {
      const CompetitorMeta &meta = metas[static_cast<std::size_t>(idx)];
      order.push_back({idx, meta.sources.size(), meta.scenario_sensitive});
    }
    std::sort(order.begin(), order.end(),
              [&](const OrderingMeta &a, const OrderingMeta &b) {
                if (a.scenario_sensitive != b.scenario_sensitive) {
                  return !a.scenario_sensitive && b.scenario_sensitive;
                }
                if (a.source_count != b.source_count) {
                  return a.source_count < b.source_count;
                }
                return a.index < b.index;
              });
    for (const OrderingMeta &rec : order) {
      uuber::CompetitorCompiledOp op;
      const CompetitorMeta &meta = metas[static_cast<std::size_t>(rec.index)];
      op.target_node_indices.push_back(meta.node_idx);
      op.transition_guard_idx =
          meta.transition_required ? meta.guard_transition_idx : -1;
      op.transition_mask_begin =
          meta.transition_required ? meta.transition_mask_begin : -1;
      op.transition_mask_count =
          meta.transition_required ? meta.transition_mask_count : 0;
      op.transition_invalidate_slot =
          meta.transition_required ? meta.invalidate_slot : 0;
      if (op.transition_mask_count > 0 && op.transition_mask_begin >= 0) {
        entry.mutates_forced_survive = true;
      }
      entry.compiled_ops.push_back(std::move(op));
    }
  }

  uuber::CompetitorCacheRecord record;
  record.competitor_ids = competitor_ids;
  record.entry = std::move(entry);
  auto &bucket = cache[key_hash];
  bucket.push_back(std::move(record));
  return bucket.back().entry;
}

double competitor_survival_from_state_compiled_ops(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache) {
  if (competitor_ids.empty()) {
    return 1.0;
  }
  const uuber::CompetitorClusterCacheEntry &cache =
      competitor_cache ? *competitor_cache
                       : fetch_competitor_cluster_cache(ctx, competitor_ids);
  if (cache.compiled_ops.empty()) {
    return 1.0;
  }

  const bool include_na_seed = state.include_na_donors;
  const int outcome_idx_seed = state.outcome_idx;
  const bool mutates_forced_survive = cache.mutates_forced_survive;
  uuber::BitsetState forced_survive_seed;
  bool forced_survive_seed_valid = false;
  if (mutates_forced_survive) {
    forced_survive_seed = state.forced_survive_bits;
    forced_survive_seed_valid = state.forced_survive_bits_valid;
  }

  state.include_na_donors = false;
  state.outcome_idx = -1;
  uuber::KernelEventEvalFn event_eval_cb =
      evaluator_make_kernel_event_eval(state);
  uuber::KernelGuardEvalFn guard_eval_cb =
      evaluator_make_kernel_guard_eval(state);
  const double out = run_competitor_compiled_ops_product(
      cache.compiled_ops, state, event_eval_cb, guard_eval_cb);

  state.include_na_donors = include_na_seed;
  state.outcome_idx = outcome_idx_seed;
  if (mutates_forced_survive) {
    state.forced_survive_bits = forced_survive_seed;
    state.forced_survive_bits_valid = forced_survive_seed_valid;
  }
  return out;
}

double competitor_survival_internal(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    double t, int component_idx,
    const uuber::BitsetState *forced_complete_bits_in,
    bool forced_complete_bits_in_valid,
    const uuber::BitsetState *forced_survive_bits_in,
    bool forced_survive_bits_in_valid, const std::string &trial_type_key,
    const TrialParamSet *trial_params,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds,
    uuber::KernelRuntimeState *kernel_runtime) {
  if (competitor_ids.empty()) {
    return 1.0;
  }
  const uuber::CompetitorClusterCacheEntry &cache =
      fetch_competitor_cluster_cache(ctx, competitor_ids);
  if (cache.compiled_ops.empty()) {
    return 1.0;
  }

  uuber::BitsetState forced_complete_bits_seed;
  uuber::BitsetState forced_survive_bits_seed;
  bool forced_complete_bits_seed_valid = false;
  bool forced_survive_bits_seed_valid = false;
  if (forced_complete_bits_in && forced_complete_bits_in_valid) {
    forced_complete_bits_seed = *forced_complete_bits_in;
    forced_complete_bits_seed_valid = true;
  }
  if (forced_survive_bits_in && forced_survive_bits_in_valid) {
    forced_survive_bits_seed = *forced_survive_bits_in;
    forced_survive_bits_seed_valid = true;
  }

  const bool kernel_runtime_usable =
      kernel_runtime &&
      !forced_bits_any(forced_complete_bits_in, forced_complete_bits_in_valid) &&
      !forced_bits_any(forced_survive_bits_in, forced_survive_bits_in_valid) &&
      exact_source_times == nullptr && source_time_bounds == nullptr;
  NodeEvalState state(
      ctx, t, component_idx, trial_params, trial_type_key, false, -1,
      exact_source_times, source_time_bounds, nullptr,
      forced_complete_bits_seed_valid ? &forced_complete_bits_seed : nullptr,
      forced_complete_bits_seed_valid,
      forced_survive_bits_seed_valid ? &forced_survive_bits_seed : nullptr,
      forced_survive_bits_seed_valid,
      kernel_runtime_usable ? kernel_runtime : nullptr);
  uuber::KernelEventEvalFn event_eval_cb =
      evaluator_make_kernel_event_eval(state);
  uuber::KernelGuardEvalFn guard_eval_cb =
      evaluator_make_kernel_guard_eval(state);
  return run_competitor_compiled_ops_product(cache.compiled_ops, state,
                                             event_eval_cb, guard_eval_cb);
}

void invalidate_kernel_runtime_root(NodeEvalState &state) {
  if (state.kernel_runtime_ready && state.kernel_runtime_ptr) {
    uuber::invalidate_kernel_runtime_from_slot(*state.kernel_runtime_ptr, 0);
  }
}
