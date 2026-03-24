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
  if (ctx.tree_program && ctx.tree_program->valid &&
      (meta.invalidate_slot < 0 ||
       meta.invalidate_slot >= static_cast<int>(ctx.tree_program->ops.size()))) {
    Rcpp::stop(
        "IR competitor transition invalidate slot invalid for node %d",
        meta.node_id);
  }
  if (meta.guard_transition_idx >= 0 &&
      meta.guard_transition_idx >=
          static_cast<int>(ctx.tree_runtime.guard_transitions.size())) {
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

[[noreturn]] inline void competitor_batch_invariant_failure(
    const char *reason, std::size_t point_count, std::size_t compiled_op_count) {
  Rcpp::stop("Competitor batch invariant failed: %s (points=%d compiled_ops=%d)",
             reason, static_cast<int>(point_count),
             static_cast<int>(compiled_op_count));
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
  if (ctx.tree_runtime.node_contains_guard.size() != ctx.ir.nodes.size() ||
      ctx.tree_runtime.node_competitor_guard_transition_idx.size() !=
          ctx.ir.nodes.size() ||
      ctx.tree_runtime.node_competitor_transition_mask_begin.size() !=
          ctx.ir.nodes.size() ||
      ctx.tree_runtime.node_competitor_transition_mask_count.size() !=
          ctx.ir.nodes.size() ||
      ctx.tree_runtime.node_competitor_transition_invalidate_slot.size() !=
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
        ctx.tree_runtime
            .node_contains_guard[static_cast<std::size_t>(node_idx)] != 0u;
    meta.scenario_sensitive =
        (node.flags & uuber::IR_NODE_FLAG_SCENARIO_SENSITIVE) != 0u;
    meta.transition_required = contains_guard && !meta.sources.empty();
    if (meta.transition_required) {
      meta.guard_transition_idx = ctx.tree_runtime
          .node_competitor_guard_transition_idx[static_cast<std::size_t>(
              node_idx)];
      meta.transition_mask_begin = ctx.tree_runtime
          .node_competitor_transition_mask_begin[static_cast<std::size_t>(
              node_idx)];
      meta.transition_mask_count = ctx.tree_runtime
          .node_competitor_transition_mask_count[static_cast<std::size_t>(
              node_idx)];
      meta.invalidate_slot = ctx.tree_runtime
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
      uuber::TreeCompetitorOp op;
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
      uuber::TreeCompetitorOp op;
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

void competitor_survival_batch_from_state_compiled_ops(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    NodeEvalState &state,
    const std::vector<double> &times, std::vector<double> &survival_out,
    const uuber::CompetitorClusterCacheEntry *competitor_cache) {
  survival_out.assign(times.size(), 1.0);
  if (competitor_ids.empty()) {
    return;
  }
  const uuber::CompetitorClusterCacheEntry &cache =
      competitor_cache ? *competitor_cache
                       : fetch_competitor_cluster_cache(ctx, competitor_ids);
  if (cache.compiled_ops.empty()) {
    return;
  }
  if (state.trial_params_soa_batch != nullptr &&
      state.trial_params_soa_batch->size() != times.size()) {
    competitor_batch_invariant_failure(
        "trial_params_soa_batch size mismatch", times.size(),
        cache.compiled_ops.size());
  }

  struct StateRestoreGuard {
    NodeEvalState &state;
    bool include_na_seed;
    int outcome_idx_seed;
    bool mutates_forced_survive;
    uuber::BitsetState forced_survive_seed;
    bool forced_survive_seed_valid;

    ~StateRestoreGuard() {
      state.include_na_donors = include_na_seed;
      state.outcome_idx = outcome_idx_seed;
      if (mutates_forced_survive) {
        state.forced_survive_bits = forced_survive_seed;
        state.forced_survive_bits_valid = forced_survive_seed_valid;
      }
    }
  };

  StateRestoreGuard restore{state,
                            state.include_na_donors,
                            state.outcome_idx,
                            cache.mutates_forced_survive,
                            state.forced_survive_bits,
                            state.forced_survive_bits_valid};

  state.include_na_donors = false;
  state.outcome_idx = -1;

  std::vector<std::uint8_t> active(times.size(), 0u);
  for (std::size_t i = 0; i < times.size(); ++i) {
    active[i] =
        (std::isfinite(times[i]) && times[i] >= 0.0 && survival_out[i] > 0.0)
            ? 1u
            : 0u;
  }

  uuber::TreeNodeBatchValues node_values;
  for (const uuber::TreeCompetitorOp &op : cache.compiled_ops) {
    for (int target_node_idx : op.target_node_indices) {
      if (!evaluator_eval_node_batch_with_state_dense(
              target_node_idx, times, state, EvalNeed::kSurvival, node_values) ||
          node_values.survival.size() != survival_out.size()) {
        competitor_batch_invariant_failure(
            "tree vector competitor product batch execution failed",
            times.size(), cache.compiled_ops.size());
      }
      for (std::size_t i = 0; i < survival_out.size(); ++i) {
        if (active[i] == 0u) {
          continue;
        }
        const double surv = clamp_probability(node_values.survival[i]);
        if (!std::isfinite(surv) || surv <= 0.0) {
          survival_out[i] = 0.0;
          active[i] = 0u;
          continue;
        }
        survival_out[i] *= surv;
        if (!std::isfinite(survival_out[i]) || survival_out[i] <= 0.0) {
          survival_out[i] = 0.0;
          active[i] = 0u;
        }
      }
    }
    if (op.transition_mask_begin < 0 || op.transition_mask_count <= 0) {
      continue;
    }
    apply_transition_mask_words(state.ctx, op.transition_mask_begin,
                                op.transition_mask_count,
                                op.transition_invalidate_slot,
                                state.forced_survive_bits,
                                state.forced_survive_bits_valid);
  }
}
