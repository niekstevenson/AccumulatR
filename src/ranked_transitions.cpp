#include "ranked_transitions.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "evaluator_internal.h"
#include "exact_outcome_density.h"

namespace {

struct SequenceStateKey {
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  int forced_complete_bit_count{0};
  int forced_survive_bit_count{0};
  std::uint64_t forced_complete_hash{0};
  std::uint64_t forced_survive_hash{0};
  std::uint64_t time_constraint_hash{0};
  std::size_t time_constraint_size{0};

  bool operator==(const SequenceStateKey &other) const noexcept {
    return forced_complete_bits_valid == other.forced_complete_bits_valid &&
           forced_survive_bits_valid == other.forced_survive_bits_valid &&
           forced_complete_bit_count == other.forced_complete_bit_count &&
           forced_survive_bit_count == other.forced_survive_bit_count &&
           forced_complete_hash == other.forced_complete_hash &&
           forced_survive_hash == other.forced_survive_hash &&
           time_constraint_hash == other.time_constraint_hash &&
           time_constraint_size == other.time_constraint_size;
  }
};

struct SequenceStateKeyHash {
  std::size_t operator()(const SequenceStateKey &key) const noexcept {
    std::size_t seed = static_cast<std::size_t>(key.forced_complete_hash);
    auto combine = [](std::size_t &seed_ref, std::size_t value) {
      seed_ref ^= value + 0x9e3779b97f4a7c15ULL + (seed_ref << 6) +
                  (seed_ref >> 2);
    };
    combine(seed, static_cast<std::size_t>(key.forced_survive_hash));
    combine(seed, static_cast<std::size_t>(key.time_constraint_hash));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bits_valid));
    combine(seed, static_cast<std::size_t>(key.forced_complete_bit_count));
    combine(seed, static_cast<std::size_t>(key.forced_survive_bit_count));
    combine(seed, key.time_constraint_size);
    return seed;
  }
};

struct RankedFrontierStorage {
  std::vector<int> trial_slot;
  std::vector<double> weight;
  std::vector<uuber::BitsetState> forced_complete_bits;
  std::vector<uuber::BitsetState> forced_survive_bits;
  std::vector<std::uint8_t> forced_complete_bits_valid;
  std::vector<std::uint8_t> forced_survive_bits_valid;
  std::vector<TimeConstraintMap> time_constraints;

  void clear() {
    trial_slot.clear();
    weight.clear();
    forced_complete_bits.clear();
    forced_survive_bits.clear();
    forced_complete_bits_valid.clear();
    forced_survive_bits_valid.clear();
    time_constraints.clear();
  }

  void reserve(std::size_t n) {
    trial_slot.reserve(n);
    weight.reserve(n);
    forced_complete_bits.reserve(n);
    forced_survive_bits.reserve(n);
    forced_complete_bits_valid.reserve(n);
    forced_survive_bits_valid.reserve(n);
    time_constraints.reserve(n);
  }

  std::size_t size() const noexcept { return trial_slot.size(); }
  bool empty() const noexcept { return trial_slot.empty(); }

  void swap(RankedFrontierStorage &other) noexcept {
    trial_slot.swap(other.trial_slot);
    weight.swap(other.weight);
    forced_complete_bits.swap(other.forced_complete_bits);
    forced_survive_bits.swap(other.forced_survive_bits);
    forced_complete_bits_valid.swap(other.forced_complete_bits_valid);
    forced_survive_bits_valid.swap(other.forced_survive_bits_valid);
    time_constraints.swap(other.time_constraints);
  }
};

struct RankedFrontierStateRef {
  std::size_t state_index{0};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  SequenceStateKey state_key;
};

struct RankedStateCollapseIndex {
  std::size_t candidate_index{0};
  int trial_slot{-1};
  SequenceStateKey state_key;
};

struct RankedFrontierGroup {
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  std::size_t begin{0};
  std::size_t end{0};
};

struct RankedLaneGroupBatch {
  ExactScenarioBatch rank_batch;
  ExactScenarioBatch denom_batch;
  std::vector<std::size_t> source_state_indices;

  void clear() {
    rank_batch.clear();
    denom_batch.clear();
    source_state_indices.clear();
  }

  void reserve(std::size_t n, bool need_denom) {
    rank_batch.reserve(n);
    if (need_denom) {
      denom_batch.reserve(n);
    }
    source_state_indices.reserve(n);
  }
};

inline std::uint64_t ranked_hash_bitset(const uuber::BitsetState &bits,
                                        bool bits_valid) {
  if (!bits_valid) {
    return 0ULL;
  }
  std::uint64_t hash = kFNV64Offset;
  const std::int32_t bit_count = static_cast<std::int32_t>(bits.bit_count());
  hash_append_bytes(hash, &bit_count, sizeof(bit_count));
  if (bits.words().empty()) {
    const std::uint64_t word = bits.small_word();
    hash_append_bytes(hash, &word, sizeof(word));
  } else {
    for (std::uint64_t word : bits.words()) {
      hash_append_bytes(hash, &word, sizeof(word));
    }
  }
  return hash;
}

inline std::uint64_t ranked_hash_time_constraints(
    const TimeConstraintMap &time_constraints) {
  std::uint64_t hash = kFNV64Offset;
  for (const auto &kv : time_constraints) {
    const std::int32_t id = static_cast<std::int32_t>(kv.first);
    hash_append_bytes(hash, &id, sizeof(id));
    const std::uint8_t has_exact = kv.second.has_exact ? 1u : 0u;
    const std::uint8_t has_lower = kv.second.has_lower ? 1u : 0u;
    const std::uint8_t has_upper = kv.second.has_upper ? 1u : 0u;
    hash_append_bytes(hash, &has_exact, sizeof(has_exact));
    hash_append_bytes(hash, &has_lower, sizeof(has_lower));
    hash_append_bytes(hash, &has_upper, sizeof(has_upper));
    if (kv.second.has_exact) {
      const std::uint64_t exact_bits = canonical_double_bits(kv.second.exact_time);
      hash_append_bytes(hash, &exact_bits, sizeof(exact_bits));
    }
    if (kv.second.has_lower) {
      const std::uint64_t lower_bits = canonical_double_bits(kv.second.lower);
      hash_append_bytes(hash, &lower_bits, sizeof(lower_bits));
    }
    if (kv.second.has_upper) {
      const std::uint64_t upper_bits = canonical_double_bits(kv.second.upper);
      hash_append_bytes(hash, &upper_bits, sizeof(upper_bits));
    }
  }
  return hash;
}

inline SequenceStateKey sequence_state_key(
    const uuber::BitsetState &forced_complete_bits,
    bool forced_complete_bits_valid,
    const uuber::BitsetState &forced_survive_bits,
    bool forced_survive_bits_valid,
    const TimeConstraintMap &time_constraints) {
  SequenceStateKey key;
  key.forced_complete_bits_valid = forced_complete_bits_valid;
  key.forced_survive_bits_valid = forced_survive_bits_valid;
  key.forced_complete_bit_count =
      forced_complete_bits_valid ? forced_complete_bits.bit_count() : 0;
  key.forced_survive_bit_count =
      forced_survive_bits_valid ? forced_survive_bits.bit_count() : 0;
  key.forced_complete_hash =
      ranked_hash_bitset(forced_complete_bits, forced_complete_bits_valid);
  key.forced_survive_hash =
      ranked_hash_bitset(forced_survive_bits, forced_survive_bits_valid);
  key.time_constraint_hash = ranked_hash_time_constraints(time_constraints);
  key.time_constraint_size = time_constraints.size();
  return key;
}

inline bool sequence_state_key_less(const SequenceStateKey &lhs,
                                    const SequenceStateKey &rhs) noexcept {
  if (lhs.forced_complete_bits_valid != rhs.forced_complete_bits_valid) {
    return lhs.forced_complete_bits_valid < rhs.forced_complete_bits_valid;
  }
  if (lhs.forced_survive_bits_valid != rhs.forced_survive_bits_valid) {
    return lhs.forced_survive_bits_valid < rhs.forced_survive_bits_valid;
  }
  if (lhs.forced_complete_bit_count != rhs.forced_complete_bit_count) {
    return lhs.forced_complete_bit_count < rhs.forced_complete_bit_count;
  }
  if (lhs.forced_survive_bit_count != rhs.forced_survive_bit_count) {
    return lhs.forced_survive_bit_count < rhs.forced_survive_bit_count;
  }
  if (lhs.forced_complete_hash != rhs.forced_complete_hash) {
    return lhs.forced_complete_hash < rhs.forced_complete_hash;
  }
  if (lhs.forced_survive_hash != rhs.forced_survive_hash) {
    return lhs.forced_survive_hash < rhs.forced_survive_hash;
  }
  if (lhs.time_constraint_hash != rhs.time_constraint_hash) {
    return lhs.time_constraint_hash < rhs.time_constraint_hash;
  }
  return lhs.time_constraint_size < rhs.time_constraint_size;
}

inline bool ranked_frontier_state_ref_less(const RankedFrontierStateRef &lhs,
                                           const RankedFrontierStateRef &rhs) noexcept {
  if (lhs.trial_params_soa != rhs.trial_params_soa) {
    return std::less<const uuber::TrialParamsSoA *>{}(lhs.trial_params_soa,
                                                      rhs.trial_params_soa);
  }
  return sequence_state_key_less(lhs.state_key, rhs.state_key);
}

inline bool ranked_frontier_state_ref_same_group(
    const RankedFrontierStateRef &lhs,
    const RankedFrontierStateRef &rhs) noexcept {
  return lhs.trial_params_soa == rhs.trial_params_soa &&
         lhs.state_key == rhs.state_key;
}

inline bool ranked_state_collapse_index_less(
    const RankedStateCollapseIndex &lhs,
    const RankedStateCollapseIndex &rhs) noexcept {
  if (lhs.trial_slot != rhs.trial_slot) {
    return lhs.trial_slot < rhs.trial_slot;
  }
  return sequence_state_key_less(lhs.state_key, rhs.state_key);
}

inline bool ranked_state_collapse_index_same_key(
    const RankedStateCollapseIndex &lhs,
    const RankedStateCollapseIndex &rhs) noexcept {
  return lhs.trial_slot == rhs.trial_slot && lhs.state_key == rhs.state_key;
}

inline SequenceStateKey ranked_frontier_state_key(
    const RankedFrontierStorage &states, std::size_t state_index) {
  return sequence_state_key(
      states.forced_complete_bits[state_index],
      state_index < states.forced_complete_bits_valid.size() &&
          states.forced_complete_bits_valid[state_index] != 0u,
      states.forced_survive_bits[state_index],
      state_index < states.forced_survive_bits_valid.size() &&
          states.forced_survive_bits_valid[state_index] != 0u,
      states.time_constraints[state_index]);
}

inline void ranked_frontier_append_state(
    RankedFrontierStorage &states, int trial_slot, double weight,
    uuber::BitsetState &&forced_complete_bits, bool forced_complete_bits_valid,
    uuber::BitsetState &&forced_survive_bits, bool forced_survive_bits_valid,
    TimeConstraintMap &&time_constraints) {
  states.trial_slot.push_back(trial_slot);
  states.weight.push_back(weight);
  states.forced_complete_bits.push_back(std::move(forced_complete_bits));
  states.forced_survive_bits.push_back(std::move(forced_survive_bits));
  states.forced_complete_bits_valid.push_back(forced_complete_bits_valid ? 1u
                                                                         : 0u);
  states.forced_survive_bits_valid.push_back(forced_survive_bits_valid ? 1u
                                                                       : 0u);
  states.time_constraints.push_back(std::move(time_constraints));
}

inline void collapse_ranked_state_candidates(
    RankedFrontierStorage &states,
    std::vector<RankedStateCollapseIndex> &collapse_order,
    double branch_eps) {
  if (states.empty() || collapse_order.empty()) {
    states.clear();
    collapse_order.clear();
    return;
  }

  std::sort(collapse_order.begin(), collapse_order.end(),
            ranked_state_collapse_index_less);

  RankedFrontierStorage collapsed_states;
  std::vector<RankedStateCollapseIndex> collapsed_order;
  collapsed_states.reserve(collapse_order.size());
  collapsed_order.reserve(collapse_order.size());

  bool have_prev = false;
  RankedStateCollapseIndex prev_key;
  for (const RankedStateCollapseIndex &entry : collapse_order) {
    if (entry.candidate_index >= states.size()) {
      continue;
    }
    const std::size_t idx = entry.candidate_index;
    const double candidate_weight = states.weight[idx];
    if (!std::isfinite(candidate_weight) || candidate_weight <= branch_eps) {
      continue;
    }
    if (have_prev && ranked_state_collapse_index_same_key(prev_key, entry)) {
      collapsed_states.weight.back() += candidate_weight;
      continue;
    }
    const std::size_t collapsed_idx = collapsed_states.size();
    ranked_frontier_append_state(
        collapsed_states, states.trial_slot[idx], candidate_weight,
        std::move(states.forced_complete_bits[idx]),
        states.forced_complete_bits_valid[idx] != 0u,
        std::move(states.forced_survive_bits[idx]),
        states.forced_survive_bits_valid[idx] != 0u,
        std::move(states.time_constraints[idx]));
    RankedStateCollapseIndex collapsed_entry = entry;
    collapsed_entry.candidate_index = collapsed_idx;
    collapsed_order.push_back(std::move(collapsed_entry));
    prev_key = entry;
    have_prev = true;
  }

  states.swap(collapsed_states);
  collapse_order.swap(collapsed_order);
}

inline bool ranked_apply_persistent_survive_sources(
    const uuber::NativeContext &ctx,
    const std::vector<int> &persistent_sources, double t,
    uuber::BitsetState &forced_survive_bits, bool &forced_survive_bits_valid,
    TimeConstraintMap &time_constraints) {
  for (int src_id : persistent_sources) {
    set_forced_id_bit_strict(ctx, src_id, forced_survive_bits,
                             forced_survive_bits_valid);
    if (!time_constraints_mark_survive(src_id, t, time_constraints)) {
      return false;
    }
  }
  return true;
}

inline bool ranked_apply_persistent_survive_sources(
    const uuber::NativeContext &ctx,
    const std::vector<int> &persistent_sources, ExactScenarioBatch &points,
    std::size_t idx) {
  if (idx >= points.size()) {
    return false;
  }
  const bool valid =
      idx < points.forced_survive_bits_valid.size() &&
      points.forced_survive_bits_valid[idx] != 0u;
  bool mutable_valid = valid;
  const bool ok = ranked_apply_persistent_survive_sources(
      ctx, persistent_sources, points.t[idx], points.forced_survive_bits[idx],
      mutable_valid, points.time_constraints[idx]);
  if (idx < points.forced_survive_bits_valid.size()) {
    points.forced_survive_bits_valid[idx] = mutable_valid ? 1u : 0u;
  }
  return ok;
}

inline void ranked_lane_batch_append_frontier_state(
    ExactScenarioBatch &batch, double t, std::size_t density_index,
    const RankedFrontierStorage &states, std::size_t state_index,
    const uuber::TrialParamsSoA *trial_params_soa) {
  batch.t.push_back(t);
  batch.density_index.push_back(density_index);
  batch.weight.push_back(1.0);
  batch.trial_params_soa.push_back(trial_params_soa);
  batch.forced_complete_bits.push_back(states.forced_complete_bits[state_index]);
  batch.forced_survive_bits.push_back(states.forced_survive_bits[state_index]);
  batch.forced_complete_bits_valid.push_back(
      states.forced_complete_bits_valid[state_index]);
  batch.forced_survive_bits_valid.push_back(
      states.forced_survive_bits_valid[state_index]);
  batch.time_constraints.push_back(states.time_constraints[state_index]);
}

inline void ranked_materialize_lane_group_batch(
    const RankedFrontierStorage &states,
    const std::vector<RankedFrontierStateRef> &active_state_refs,
    const std::vector<const double *> &times_by_trial, std::size_t rank_idx,
    std::size_t group_begin, std::size_t group_end,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    RankedLaneGroupBatch &lane_batch) {
  lane_batch.clear();
  const bool need_denom = rank_idx > 0u;
  lane_batch.reserve(group_end - group_begin, need_denom);

  std::size_t local_idx = 0u;
  for (std::size_t ref_idx = group_begin; ref_idx < group_end;
       ++ref_idx, ++local_idx) {
    const RankedFrontierStateRef &state_ref = active_state_refs[ref_idx];
    const std::size_t state_index = state_ref.state_index;
    const std::size_t trial_slot =
        static_cast<std::size_t>(states.trial_slot[state_index]);
    const uuber::TrialParamsSoA *trial_params_soa =
        uniform_trial_params_soa ? uniform_trial_params_soa
                                 : state_ref.trial_params_soa;
    ranked_lane_batch_append_frontier_state(
        lane_batch.rank_batch, times_by_trial[trial_slot][rank_idx], local_idx,
        states, state_index, trial_params_soa);
    if (need_denom) {
      ranked_lane_batch_append_frontier_state(
          lane_batch.denom_batch, times_by_trial[trial_slot][rank_idx - 1u],
          local_idx, states, state_index, trial_params_soa);
    }
    lane_batch.source_state_indices.push_back(state_index);
  }
}

inline SequenceStateKey ranked_scenario_state_key(const ExactScenarioBatch &points,
                                                  std::size_t idx) {
  return sequence_state_key(
      points.forced_complete_bits[idx],
      idx < points.forced_complete_bits_valid.size() &&
          points.forced_complete_bits_valid[idx] != 0u,
      points.forced_survive_bits[idx],
      idx < points.forced_survive_bits_valid.size() &&
          points.forced_survive_bits_valid[idx] != 0u,
      points.time_constraints[idx]);
}

inline void collapse_ranked_scenario_candidates_into_frontier(
    ExactScenarioBatch &points, std::vector<int> &trial_slots,
    RankedFrontierStorage &next_states,
    std::vector<RankedStateCollapseIndex> &next_state_collapse_order,
    std::vector<RankedStateCollapseIndex> &collapse_order,
    double branch_eps);

inline bool ranked_reduce_scenario_batch_into_frontier(
    const uuber::NativeContext &ctx, const RankedFrontierStorage &states,
    const std::vector<std::size_t> &source_state_indices,
    const std::vector<double> &denom,
    const std::vector<int> &persistent_sources,
    bool apply_persistent_sources,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> *scenario_survival, ExactScenarioBatch &points,
    RankedFrontierStorage &next_states,
    std::vector<RankedStateCollapseIndex> &next_state_collapse_order,
    std::vector<int> &scenario_trial_slots,
    std::vector<RankedStateCollapseIndex> &scenario_collapse_order,
    double branch_eps) {
  scenario_trial_slots.assign(points.size(), -1);
  scenario_collapse_order.clear();
  scenario_collapse_order.reserve(points.size());

  for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
    if (active_mask != nullptr &&
        (point_idx >= active_mask->size() || (*active_mask)[point_idx] == 0u)) {
      continue;
    }
    if (points.density_index[point_idx] >= source_state_indices.size()) {
      return false;
    }
    const std::size_t lane_idx = points.density_index[point_idx];
    const std::size_t source_state_index = source_state_indices[lane_idx];
    const double denom_val = (lane_idx < denom.size()) ? denom[lane_idx] : 0.0;
    if (!std::isfinite(denom_val) || denom_val <= branch_eps) {
      continue;
    }
    double weight =
        states.weight[source_state_index] * (points.weight[point_idx] / denom_val);
    if (!std::isfinite(weight) || weight <= branch_eps) {
      continue;
    }
    if (scenario_survival != nullptr) {
      if (point_idx >= scenario_survival->size()) {
        return false;
      }
      const double surv = clamp_probability((*scenario_survival)[point_idx]);
      if (!std::isfinite(surv) || surv <= branch_eps) {
        continue;
      }
      weight *= surv;
      if (!std::isfinite(weight) || weight <= branch_eps) {
        continue;
      }
    }
    if (apply_persistent_sources && !persistent_sources.empty()) {
      if (!ranked_apply_persistent_survive_sources(ctx, persistent_sources,
                                                   points, point_idx)) {
        continue;
      }
    }

    points.weight[point_idx] = weight;
    scenario_trial_slots[point_idx] = states.trial_slot[source_state_index];
    RankedStateCollapseIndex collapse_index;
    collapse_index.candidate_index = point_idx;
    collapse_index.trial_slot = states.trial_slot[source_state_index];
    collapse_index.state_key = ranked_scenario_state_key(points, point_idx);
    scenario_collapse_order.push_back(std::move(collapse_index));
  }

  collapse_ranked_scenario_candidates_into_frontier(
      points, scenario_trial_slots, next_states, next_state_collapse_order,
      scenario_collapse_order, branch_eps);
  return true;
}

inline void collapse_ranked_scenario_candidates_into_frontier(
    ExactScenarioBatch &points, std::vector<int> &trial_slots,
    RankedFrontierStorage &next_states,
    std::vector<RankedStateCollapseIndex> &next_state_collapse_order,
    std::vector<RankedStateCollapseIndex> &collapse_order,
    double branch_eps) {
  if (points.empty() || trial_slots.size() != points.size() ||
      collapse_order.empty()) {
    points.clear();
    trial_slots.clear();
    collapse_order.clear();
    return;
  }

  std::sort(collapse_order.begin(), collapse_order.end(),
            ranked_state_collapse_index_less);

  next_states.reserve(next_states.size() + collapse_order.size());
  next_state_collapse_order.reserve(next_state_collapse_order.size() +
                                    collapse_order.size());

  bool have_prev = false;
  RankedStateCollapseIndex prev_key;
  std::size_t appended_index = 0u;
  for (const RankedStateCollapseIndex &entry : collapse_order) {
    if (entry.candidate_index >= points.size() ||
        entry.candidate_index >= trial_slots.size()) {
      continue;
    }
    const std::size_t point_idx = entry.candidate_index;
    const int trial_slot = trial_slots[point_idx];
    const double point_weight = points.weight[point_idx];
    if (trial_slot < 0 || !std::isfinite(point_weight) ||
        point_weight <= branch_eps) {
      continue;
    }
    if (have_prev && ranked_state_collapse_index_same_key(prev_key, entry)) {
      next_states.weight[appended_index] += point_weight;
      continue;
    }

    appended_index = next_states.size();
    ranked_frontier_append_state(
        next_states, trial_slot, point_weight,
        std::move(points.forced_complete_bits[point_idx]),
        point_idx < points.forced_complete_bits_valid.size() &&
            points.forced_complete_bits_valid[point_idx] != 0u,
        std::move(points.forced_survive_bits[point_idx]),
        point_idx < points.forced_survive_bits_valid.size() &&
            points.forced_survive_bits_valid[point_idx] != 0u,
        std::move(points.time_constraints[point_idx]));
    RankedStateCollapseIndex collapse_index;
    collapse_index.candidate_index = appended_index;
    collapse_index.trial_slot = trial_slot;
    collapse_index.state_key = entry.state_key;
    next_state_collapse_order.push_back(std::move(collapse_index));
    prev_key = entry;
    have_prev = true;
  }

  points.clear();
  trial_slots.clear();
  collapse_order.clear();
}

inline std::vector<int>
ranked_plan_child_ids(const uuber::NativeContext &ctx,
                      const uuber::IrNode &node) {
  std::vector<int> out;
  if (node.child_begin >= 0 && node.child_count > 0 &&
      node.child_begin + node.child_count <=
          static_cast<int>(ctx.ir.node_children.size())) {
    out.assign(ctx.ir.node_children.begin() + node.child_begin,
               ctx.ir.node_children.begin() + node.child_begin +
                   node.child_count);
  }
  return out;
}

inline std::vector<int>
ranked_source_bits_for_ids(const uuber::NativeContext &ctx,
                           const std::vector<int> &source_ids) {
  std::vector<int> out;
  out.reserve(source_ids.size());
  for (int id : source_ids) {
    int bit_idx = -1;
    auto bit_it = ctx.ir.label_id_to_bit_idx.find(id);
    if (bit_it != ctx.ir.label_id_to_bit_idx.end()) {
      bit_idx = bit_it->second;
    }
    out.push_back(bit_idx);
  }
  return out;
}

inline int ranked_program_slot_for_node(const uuber::NativeContext &ctx,
                                        int node_idx) {
  if (!ctx.tree_program || !ctx.tree_program->valid) {
    return 0;
  }
  if (node_idx < 0 || node_idx >=
                          static_cast<int>(
                              ctx.tree_program->outputs.node_idx_to_slot.size())) {
    return 0;
  }
  const int out_slot = ctx.tree_program->outputs.node_idx_to_slot
      [static_cast<std::size_t>(node_idx)];
  if (out_slot < 0 || out_slot >= static_cast<int>(ctx.tree_program->ops.size())) {
    return 0;
  }
  return out_slot;
}

inline RankedProgramEvalRef ranked_make_eval_ref(
    const uuber::NativeContext &ctx, RankedProgramEvalKind kind, int node_idx) {
  RankedProgramEvalRef ref;
  ref.kind = kind;
  ref.node_idx = node_idx;
  ref.slot = ranked_program_slot_for_node(ctx, node_idx);
  return ref;
}

inline void ranked_attach_delta_for_node(const uuber::NativeContext &ctx,
                                         int node_idx,
                                         RankedStateDelta &delta) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    return;
  }
  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  if (node.source_mask_begin < 0 || node.source_mask_count <= 0) {
    return;
  }
  if (node.source_mask_begin + node.source_mask_count >
      static_cast<int>(ctx.ir.node_source_masks.size())) {
    return;
  }
  delta.source_mask_begin = node.source_mask_begin;
  delta.source_mask_count = node.source_mask_count;
  delta.invalidate_slot = ranked_program_slot_for_node(ctx, node_idx);
}

inline bool ranked_batch_plan_is_deterministic(const RankedBatchPlan &plan) {
  if (!plan.valid || plan.slice.evals.size() != 1u ||
      plan.slice.evals.front().kind != RankedProgramEvalKind::Density) {
    return false;
  }
  for (const RankedStateDelta &delta : plan.deltas) {
    if (delta.kind == RankedStateDeltaKind::OrWitnessSurvive) {
      return false;
    }
  }
  return true;
}

inline void ranked_finalize_batch_plan(RankedBatchPlan &plan) {
  plan.valid = !plan.slice.empty();
  plan.deterministic = ranked_batch_plan_is_deterministic(plan);
}

} // namespace

RankedBatchPlanner::RankedBatchPlanner(const uuber::NativeContext &ctx_)
    : ctx(ctx_), plans(ctx_.ir.nodes.size()) {}

const RankedNodeBatchPlan &RankedBatchPlanner::plan_for_node(int node_idx) {
  static const RankedNodeBatchPlan kInvalidPlan;
  if (node_idx < 0 || node_idx >= static_cast<int>(plans.size())) {
    return kInvalidPlan;
  }
  RankedNodeBatchPlan &plan = plans[static_cast<std::size_t>(node_idx)];
  if (!plan.compiled) {
    compile_node(node_idx);
  }
  return plan;
}

void RankedBatchPlanner::compile_node(int node_idx) {
  RankedNodeBatchPlan &plan = plans[static_cast<std::size_t>(node_idx)];
  if (plan.compiled || plan.compiling) {
    return;
  }
  plan.compiling = true;
  plan.batch_plans.clear();
  plan.valid = false;

  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
    plan.compiled = true;
    plan.compiling = false;
    return;
  }

  const uuber::IrNode &node = ctx.ir.nodes[static_cast<std::size_t>(node_idx)];
  switch (node.op) {
  case uuber::IrNodeOp::EventAcc:
  case uuber::IrNodeOp::EventPool: {
    RankedBatchPlan batch_plan;
    batch_plan.slice.evals.push_back(
        ranked_make_eval_ref(ctx, RankedProgramEvalKind::Density, node_idx));
    std::vector<int> source_ids = ensure_source_ids(ctx, node);
    if (!source_ids.empty()) {
      RankedStateDelta delta;
      delta.kind = RankedStateDeltaKind::CompleteSources;
      delta.bind_exact_current_time =
          (node.op == uuber::IrNodeOp::EventAcc);
      delta.source_ids = std::move(source_ids);
      delta.source_bits = ranked_source_bits_for_ids(ctx, delta.source_ids);
      ranked_attach_delta_for_node(ctx, node_idx, delta);
      batch_plan.deltas.push_back(std::move(delta));
    }
    ranked_finalize_batch_plan(batch_plan);
    plan.batch_plans.push_back(std::move(batch_plan));
    break;
  }
  case uuber::IrNodeOp::And: {
    std::vector<int> child_ids = ranked_plan_child_ids(ctx, node);
    if (!child_ids.empty()) {
      std::vector<std::vector<int>> child_sources;
      child_sources.reserve(child_ids.size());
      for (int child_idx : child_ids) {
        child_sources.push_back(
            ensure_source_ids(ctx, ir_node_required(ctx, child_idx)));
      }
      for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
        const RankedNodeBatchPlan &child_plan =
            plan_for_node(child_ids[idx]);
        if (!child_plan.valid || child_plan.batch_plans.empty()) {
          continue;
        }
        for (const RankedBatchPlan &child_batch_plan : child_plan.batch_plans) {
          RankedBatchPlan batch_plan = child_batch_plan;
          for (std::size_t other_idx = 0; other_idx < child_ids.size();
               ++other_idx) {
            if (other_idx == idx) {
              continue;
            }
            batch_plan.slice.evals.push_back(ranked_make_eval_ref(
                ctx, RankedProgramEvalKind::CDF, child_ids[other_idx]));
            if (!child_sources[other_idx].empty()) {
              RankedStateDelta delta;
              delta.kind = RankedStateDeltaKind::CompleteSources;
              delta.source_ids = child_sources[other_idx];
              delta.source_bits =
                  ranked_source_bits_for_ids(ctx, delta.source_ids);
              ranked_attach_delta_for_node(ctx, child_ids[other_idx], delta);
              batch_plan.deltas.push_back(std::move(delta));
            }
          }
          ranked_finalize_batch_plan(batch_plan);
          plan.batch_plans.push_back(std::move(batch_plan));
        }
      }
    }
    break;
  }
  case uuber::IrNodeOp::Or: {
    std::vector<int> child_ids = ranked_plan_child_ids(ctx, node);
    if (!child_ids.empty()) {
      std::vector<std::vector<int>> child_sources;
      child_sources.reserve(child_ids.size());
      for (int child_idx : child_ids) {
        child_sources.push_back(
            ensure_source_ids(ctx, ir_node_required(ctx, child_idx)));
      }
      for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
        const RankedNodeBatchPlan &child_plan =
            plan_for_node(child_ids[idx]);
        if (!child_plan.valid || child_plan.batch_plans.empty()) {
          continue;
        }
        for (const RankedBatchPlan &child_batch_plan : child_plan.batch_plans) {
          RankedBatchPlan batch_plan = child_batch_plan;
          for (std::size_t other_idx = 0; other_idx < child_ids.size();
               ++other_idx) {
            if (other_idx == idx || child_sources[other_idx].empty()) {
              continue;
            }
            RankedStateDelta delta;
            delta.kind = RankedStateDeltaKind::OrWitnessSurvive;
            delta.source_ids = child_sources[other_idx];
            delta.source_bits = ranked_source_bits_for_ids(ctx, delta.source_ids);
            batch_plan.deltas.push_back(std::move(delta));
          }
          ranked_finalize_batch_plan(batch_plan);
          plan.batch_plans.push_back(std::move(batch_plan));
        }
      }
    }
    break;
  }
  case uuber::IrNodeOp::Guard: {
    if (node.reference_idx >= 0 && node.blocker_idx >= 0) {
      const RankedNodeBatchPlan &ref_plan = plan_for_node(node.reference_idx);
      if (ref_plan.valid && !ref_plan.batch_plans.empty()) {
        const std::vector<int> blocker_sources =
            evaluator_gather_blocker_sources(ctx, node_idx);
        for (const RankedBatchPlan &ref_batch_plan : ref_plan.batch_plans) {
          RankedBatchPlan batch_plan = ref_batch_plan;
          batch_plan.slice.evals.push_back(ranked_make_eval_ref(
              ctx, RankedProgramEvalKind::Survival, node.blocker_idx));
          if (!blocker_sources.empty()) {
            RankedStateDelta delta;
            delta.kind = RankedStateDeltaKind::SurviveSources;
            delta.source_ids = blocker_sources;
            delta.source_bits = ranked_source_bits_for_ids(ctx, delta.source_ids);
            ranked_attach_delta_for_node(ctx, node.blocker_idx, delta);
            batch_plan.deltas.push_back(std::move(delta));
          }
          ranked_finalize_batch_plan(batch_plan);
          plan.batch_plans.push_back(std::move(batch_plan));
        }
      }
    }
    break;
  }
  default:
    break;
  }

  plan.valid = !plan.batch_plans.empty();
  plan.compiled = true;
  plan.compiling = false;
}

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids) {
  if (competitor_ids.empty()) {
    return {};
  }
  std::vector<int> out;
  for (int node_id : competitor_ids) {
    if (node_id == NA_INTEGER) {
      continue;
    }
    std::vector<int> source_ids = ensure_source_ids(ctx, node_id);
    out.insert(out.end(), source_ids.begin(), source_ids.end());
  }
  sort_unique(out);
  return out;
}

bool sequence_prefix_density_batch_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices,
    const std::vector<const double *> &times_by_trial, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_by_trial,
    RankedBatchPlanner *transition_planner,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs,
    const std::vector<std::vector<int>> *step_persistent_sources) {
  density_out.assign(times_by_trial.size(), 0.0);
  const std::size_t trial_count = times_by_trial.size();
  const std::size_t rank_count = outcome_indices.size();
  if (trial_count == 0u) {
    return true;
  }
  if (outcome_indices.empty() || node_indices.size() != rank_count ||
      step_competitor_ids_ptrs == nullptr ||
      step_persistent_sources == nullptr ||
      step_competitor_ids_ptrs->size() != rank_count ||
      step_persistent_sources->size() != rank_count) {
    return false;
  }
  for (const double *times_ptr : times_by_trial) {
    if (times_ptr == nullptr) {
      return false;
    }
  }

  RankedBatchPlanner local_transition_planner(ctx);
  RankedBatchPlanner *planner =
      transition_planner ? transition_planner : &local_transition_planner;
  const uuber::TrialParamsSoA *default_trial_params_soa =
      resolve_trial_params_soa(ctx, trial_params);

  constexpr double kBranchEps = 1e-18;
  RankedFrontierStorage states;
  states.reserve(trial_count);
  for (std::size_t trial_slot = 0; trial_slot < trial_count; ++trial_slot) {
    uuber::BitsetState forced_complete_bits;
    uuber::BitsetState forced_survive_bits;
    bool forced_complete_bits_valid = false;
    bool forced_survive_bits_valid = false;
    (void)ensure_forced_bitset_capacity(ctx, forced_complete_bits,
                                        forced_complete_bits_valid);
    (void)ensure_forced_bitset_capacity(ctx, forced_survive_bits,
                                        forced_survive_bits_valid);
    ranked_frontier_append_state(states, static_cast<int>(trial_slot), 1.0,
                                 std::move(forced_complete_bits),
                                 forced_complete_bits_valid,
                                 std::move(forced_survive_bits),
                                 forced_survive_bits_valid,
                                 TimeConstraintMap{});
  }

  for (std::size_t rank_idx = 0; rank_idx < rank_count; ++rank_idx) {
    const int outcome_idx = outcome_indices[rank_idx];
    if (outcome_idx < 0 ||
        outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      return false;
    }
    const int outcome_node_idx = node_indices[rank_idx];
    if (outcome_node_idx < 0 ||
        outcome_node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
      return false;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    const int info_node_idx = resolve_dense_node_idx_required(ctx, info.node_id);
    static const std::vector<int> kEmptyCompetitorIds;
    const std::vector<int> &competitors =
        ((*step_competitor_ids_ptrs)[rank_idx] != nullptr)
            ? *(*step_competitor_ids_ptrs)[rank_idx]
            : kEmptyCompetitorIds;
    const std::vector<int> &persistent_sources =
        (*step_persistent_sources)[rank_idx];

    RankedFrontierStorage next_states_collapsed;
    next_states_collapsed.reserve(states.size() * 2u);
    std::vector<RankedStateCollapseIndex> next_state_collapse_order;
    next_state_collapse_order.reserve(states.size() * 2u);

    const bool use_frontier_groups = (trial_params_soa_by_trial != nullptr);
    std::vector<RankedFrontierStateRef> active_state_refs;
    active_state_refs.reserve(states.size());
    std::vector<RankedFrontierGroup> frontier_groups;
    if (use_frontier_groups) {
      frontier_groups.reserve(states.size());
    }

    bool valid_rank = true;
    const uuber::TrialParamsSoA *frontier_uniform_trial_params_soa =
        default_trial_params_soa;
    bool frontier_uniform_trial_params_known = false;
    for (std::size_t state_idx = 0; state_idx < states.size(); ++state_idx) {
      if (!std::isfinite(states.weight[state_idx]) ||
          states.weight[state_idx] <= kBranchEps) {
        continue;
      }
      const int trial_slot = states.trial_slot[state_idx];
      if (trial_slot < 0 ||
          static_cast<std::size_t>(trial_slot) >= times_by_trial.size()) {
        valid_rank = false;
        break;
      }
      const double t =
          times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx];
      if (!std::isfinite(t) || t < 0.0) {
        valid_rank = false;
        break;
      }
      if (rank_idx > 0u) {
        const double lower_t =
            times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx - 1u];
        if (!std::isfinite(lower_t) || lower_t < 0.0) {
          valid_rank = false;
          break;
        }
      }

      RankedFrontierStateRef state_ref;
      state_ref.state_index = state_idx;
      state_ref.trial_params_soa = default_trial_params_soa;
      if (trial_params_soa_by_trial != nullptr &&
          static_cast<std::size_t>(trial_slot) < trial_params_soa_by_trial->size()) {
        state_ref.trial_params_soa =
            (*trial_params_soa_by_trial)[static_cast<std::size_t>(trial_slot)];
      }
      state_ref.state_key = ranked_frontier_state_key(states, state_idx);
      active_state_refs.push_back(std::move(state_ref));

      const uuber::TrialParamsSoA *state_trial_params_soa =
          active_state_refs.back().trial_params_soa;
      if (!frontier_uniform_trial_params_known) {
        frontier_uniform_trial_params_soa = state_trial_params_soa;
        frontier_uniform_trial_params_known = true;
      } else if (frontier_uniform_trial_params_soa != state_trial_params_soa) {
        frontier_uniform_trial_params_soa = nullptr;
      }
    }
    if (!valid_rank) {
      return false;
    }
    if (active_state_refs.empty()) {
      return true;
    }
    if (use_frontier_groups && active_state_refs.size() > 1u) {
      std::sort(active_state_refs.begin(), active_state_refs.end(),
                ranked_frontier_state_ref_less);
      std::size_t group_begin = 0u;
      while (group_begin < active_state_refs.size()) {
        std::size_t group_end = group_begin + 1u;
        while (group_end < active_state_refs.size() &&
               ranked_frontier_state_ref_same_group(
                   active_state_refs[group_begin],
                   active_state_refs[group_end])) {
          ++group_end;
        }
        RankedFrontierGroup group;
        group.trial_params_soa = active_state_refs[group_begin].trial_params_soa;
        group.begin = group_begin;
        group.end = group_end;
        frontier_groups.push_back(group);
        group_begin = group_end;
      }
    }

    RankedLaneGroupBatch lane_group_batch;
    lane_group_batch.reserve(active_state_refs.size(), rank_idx > 0u);
    std::vector<double> denom;
    denom.reserve(active_state_refs.size());
    ExactScenarioBatch scenario_points;
    scenario_points.reserve(active_state_refs.size() * 2u);
    std::vector<double> scenario_survival;
    scenario_survival.reserve(active_state_refs.size() * 2u);
    std::vector<int> scenario_trial_slots;
    scenario_trial_slots.reserve(active_state_refs.size() * 2u);
    std::vector<RankedStateCollapseIndex> scenario_collapse_order;
    scenario_collapse_order.reserve(active_state_refs.size() * 2u);
    ExactScenarioBatch deterministic_scenario_points;
    deterministic_scenario_points.reserve(active_state_refs.size());
    std::vector<std::uint8_t> deterministic_scenario_active;
    deterministic_scenario_active.reserve(active_state_refs.size());
    std::vector<double> deterministic_competitor_survival;
    deterministic_competitor_survival.reserve(active_state_refs.size());
    const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr = nullptr;
    if (!competitors.empty()) {
      competitor_cache_ptr = &fetch_competitor_cluster_cache(ctx, competitors);
    }

    const bool split_frontier_groups =
        use_frontier_groups && frontier_groups.size() > 1u &&
        frontier_groups.size() * 2u < active_state_refs.size();
    const std::size_t frontier_begin = 0u;
    const std::size_t frontier_end = active_state_refs.size();

    auto process_state_group =
        [&](std::size_t group_begin, std::size_t group_end,
            const uuber::TrialParamsSoA *uniform_trial_params_soa) -> bool {
      lane_group_batch.clear();
      denom.clear();
      ranked_materialize_lane_group_batch(
          states, active_state_refs, times_by_trial, rank_idx, group_begin,
          group_end, uniform_trial_params_soa, lane_group_batch);

      const std::size_t group_size = group_end - group_begin;
      denom.assign(group_size, 1.0);
      if (rank_idx > 0u) {
        uuber::KernelNodeBatchValues denom_values;
        if (!exact_eval_node_batch_from_batch(
                ctx, info_node_idx, lane_group_batch.denom_batch,
                component_idx, trial_params, trial_type_key,
                EvalNeed::kSurvival, denom_values, uniform_trial_params_soa) ||
            denom_values.survival.size() != denom.size()) {
          return false;
        }
        for (std::size_t i = 0; i < denom.size(); ++i) {
          denom[i] = denom_values.survival[i];
        }
      }

      scenario_points.clear();
      scenario_trial_slots.clear();
      scenario_collapse_order.clear();

      deterministic_scenario_points.clear();
      deterministic_scenario_active.clear();
      if (exact_collect_deterministic_scenarios_batch_aligned_from_batch(
              *planner, ctx, outcome_node_idx, lane_group_batch.rank_batch,
              component_idx, trial_params, trial_type_key,
              deterministic_scenario_points, deterministic_scenario_active,
              uniform_trial_params_soa)) {
        deterministic_competitor_survival.clear();
        if (competitor_cache_ptr != nullptr) {
          deterministic_competitor_survival.assign(
              deterministic_scenario_points.size(), 1.0);
          if (!deterministic_scenario_points.empty()) {
            exact_competitor_survival_batch(
                ctx, *competitor_cache_ptr, component_idx, trial_params,
                trial_type_key, deterministic_scenario_points,
                deterministic_competitor_survival, uniform_trial_params_soa);
          }
        }

        return ranked_reduce_scenario_batch_into_frontier(
            ctx, states, lane_group_batch.source_state_indices, denom,
            persistent_sources, competitor_cache_ptr != nullptr,
            &deterministic_scenario_active,
            competitor_cache_ptr != nullptr ? &deterministic_competitor_survival
                                            : nullptr,
            deterministic_scenario_points, next_states_collapsed,
            next_state_collapse_order, scenario_trial_slots,
            scenario_collapse_order, kBranchEps);
      }

      if (!exact_collect_scenarios_batch_aligned_from_batch(
              *planner, ctx, outcome_node_idx, lane_group_batch.rank_batch,
              component_idx, trial_params, trial_type_key, scenario_points,
              uniform_trial_params_soa)) {
        return true;
      }

      scenario_survival.assign(scenario_points.size(), 1.0);
      if (competitor_cache_ptr != nullptr) {
        exact_competitor_survival_batch(
            ctx, *competitor_cache_ptr, component_idx, trial_params,
            trial_type_key, scenario_points, scenario_survival,
            uniform_trial_params_soa);
      }

      return ranked_reduce_scenario_batch_into_frontier(
          ctx, states, lane_group_batch.source_state_indices, denom,
          persistent_sources, competitor_cache_ptr != nullptr, nullptr,
          competitor_cache_ptr != nullptr ? &scenario_survival : nullptr,
          scenario_points, next_states_collapsed, next_state_collapse_order,
          scenario_trial_slots, scenario_collapse_order, kBranchEps);
    };

    if (split_frontier_groups) {
      for (const RankedFrontierGroup &group : frontier_groups) {
        if (!process_state_group(group.begin, group.end, group.trial_params_soa)) {
          return false;
        }
      }
    } else if (!process_state_group(frontier_begin, frontier_end,
                                    frontier_uniform_trial_params_soa)) {
      return false;
    }

    collapse_ranked_state_candidates(next_states_collapsed, next_state_collapse_order,
                                     kBranchEps);
    states.swap(next_states_collapsed);
    if (states.empty()) {
      return true;
    }
  }

  for (std::size_t state_idx = 0; state_idx < states.size(); ++state_idx) {
    const int trial_slot = states.trial_slot[state_idx];
    const double weight = states.weight[state_idx];
    if (!std::isfinite(weight) || weight <= 0.0 || trial_slot < 0 ||
        static_cast<std::size_t>(trial_slot) >= density_out.size()) {
      continue;
    }
    density_out[static_cast<std::size_t>(trial_slot)] += weight;
  }
  for (double &value : density_out) {
    if (!std::isfinite(value) || value <= 0.0) {
      value = 0.0;
    }
  }

  return true;
}
