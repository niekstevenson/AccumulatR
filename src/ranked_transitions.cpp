#include "ranked_transitions.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "competitor_cache.h"
#include "evaluator_internal.h"
#include "exact_outcome_density.h"

namespace {

struct RankedFrontierCollapseKey {
  int trial_slot{-1};
  SequenceStateKey state_key;

  bool operator==(const RankedFrontierCollapseKey &other) const noexcept {
    return trial_slot == other.trial_slot && state_key == other.state_key;
  }
};

struct RankedFrontierCollapseKeyHash {
  std::size_t operator()(const RankedFrontierCollapseKey &key) const noexcept {
    std::size_t seed = std::hash<int>{}(key.trial_slot);
    seed ^= SequenceStateKeyHash{}(key.state_key) + 0x9e3779b97f4a7c15ULL +
            (seed << 6) + (seed >> 2);
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
  std::vector<SequenceStateKey> state_keys;

  void clear() {
    trial_slot.clear();
    weight.clear();
    forced_complete_bits.clear();
    forced_survive_bits.clear();
    forced_complete_bits_valid.clear();
    forced_survive_bits_valid.clear();
    time_constraints.clear();
    state_keys.clear();
  }

  void reserve(std::size_t n) {
    trial_slot.reserve(n);
    weight.reserve(n);
    forced_complete_bits.reserve(n);
    forced_survive_bits.reserve(n);
    forced_complete_bits_valid.reserve(n);
    forced_survive_bits_valid.reserve(n);
    time_constraints.reserve(n);
    state_keys.reserve(n);
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
    state_keys.swap(other.state_keys);
  }
};

struct RankedStateCollapseIndex {
  std::size_t candidate_index{0};
  int trial_slot{-1};
  SequenceStateKey state_key;
};

struct RankedFrontierStateRef {
  std::size_t state_index{0};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
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

struct RankedFrontierExecutionScratch {
  RankedLaneGroupBatch lane_group_batch;
  std::vector<double> denom;
  std::vector<double> scenario_survival;

  void reserve(std::size_t n) {
    lane_group_batch.reserve(n, true);
    denom.reserve(n);
    scenario_survival.reserve(n * 2u);
  }

  void clear_materialized() {
    lane_group_batch.clear();
    denom.clear();
    scenario_survival.clear();
  }
};

template <typename PointVec>
inline bool ranked_reduce_exact_batch_into_frontier(
    const uuber::NativeContext &ctx, const PointVec &points,
    const RankedFrontierStorage &states,
    const std::vector<std::size_t> &source_state_indices,
    const std::vector<double> &scenario_weights, const std::vector<double> &denom,
    const std::vector<int> &persistent_sources,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> *scenario_survival,
    RankedFrontierStorage &next_states,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order,
    double branch_eps);

class RankedFrontierReducer final : public ExactTransitionReducer {
public:
  RankedFrontierReducer(
      const uuber::NativeContext &ctx_, const RankedFrontierStorage &states_,
      const std::vector<std::size_t> &source_state_indices_,
      const std::vector<double> &denom_,
      const std::vector<int> &persistent_sources_,
      const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr_,
      int component_idx_, const TrialParamSet *trial_params_,
      const std::string &trial_type_key_,
      const uuber::TrialParamsSoA *uniform_trial_params_soa_,
      std::vector<double> &scenario_survival_,
      RankedFrontierStorage &next_states_,
      std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                         RankedFrontierCollapseKeyHash> *next_state_lookup_,
      std::vector<RankedStateCollapseIndex> *next_state_collapse_order_,
      double branch_eps_)
      : ctx(ctx_), states(states_), source_state_indices(source_state_indices_),
        denom(denom_), persistent_sources(persistent_sources_),
        competitor_cache_ptr(competitor_cache_ptr_),
        component_idx(component_idx_), trial_params(trial_params_),
        trial_type_key(trial_type_key_),
        uniform_trial_params_soa(uniform_trial_params_soa_),
        scenario_survival(scenario_survival_), next_states(next_states_),
        next_state_lookup(next_state_lookup_),
        next_state_collapse_order(next_state_collapse_order_),
        branch_eps(branch_eps_) {}

  bool consume(const ExactScenarioBatch &points,
               const std::vector<std::uint8_t> *active_mask,
               const std::vector<double> &weights) override {
    return consume_impl(points, active_mask, weights);
  }

  bool consume(const uuber::ExactScenarioLaneViewBatch &points,
               const std::vector<std::uint8_t> *active_mask,
               const std::vector<double> &weights) override {
    return consume_impl(points, active_mask, weights);
  }

private:
  template <typename PointVec>
  bool consume_impl(const PointVec &points,
                    const std::vector<std::uint8_t> *active_mask,
                    const std::vector<double> &weights) {
    scenario_survival.clear();
    if (competitor_cache_ptr != nullptr) {
      scenario_survival.assign(points.size(), 1.0);
      exact_competitor_survival_batch(
          ctx, *competitor_cache_ptr, component_idx, trial_params,
          trial_type_key, points, scenario_survival, uniform_trial_params_soa);
    }
    return ranked_reduce_exact_batch_into_frontier(
        ctx, points, states, source_state_indices, weights, denom,
        persistent_sources, active_mask,
        competitor_cache_ptr != nullptr ? &scenario_survival : nullptr,
        next_states, next_state_lookup, next_state_collapse_order, branch_eps);
  }

  const uuber::NativeContext &ctx;
  const RankedFrontierStorage &states;
  const std::vector<std::size_t> &source_state_indices;
  const std::vector<double> &denom;
  const std::vector<int> &persistent_sources;
  const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr;
  int component_idx;
  const TrialParamSet *trial_params;
  const std::string &trial_type_key;
  const uuber::TrialParamsSoA *uniform_trial_params_soa;
  std::vector<double> &scenario_survival;
  RankedFrontierStorage &next_states;
  std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                     RankedFrontierCollapseKeyHash> *next_state_lookup;
  std::vector<RankedStateCollapseIndex> *next_state_collapse_order;
  double branch_eps;
};

struct RankedActiveFrontier {
  bool valid{true};
  const uuber::TrialParamsSoA *uniform_trial_params_soa{nullptr};
  std::vector<RankedFrontierStateRef> refs;
  std::vector<RankedFrontierGroup> groups;
  bool split_groups{false};
};

struct RankedStepExecutionContext {
  const RankedNodeBatchPlan *transition_plan{nullptr};
  int info_node_idx{-1};
  const std::vector<int> *persistent_sources{nullptr};
  const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr{nullptr};
};

inline SequenceStateKey ranked_frontier_state_key(
    const RankedFrontierStorage &states, std::size_t state_index) {
  return state_index < states.state_keys.size() ? states.state_keys[state_index]
                                                : SequenceStateKey{};
}

inline bool ranked_frontier_state_ref_less(
    const RankedFrontierStateRef &lhs,
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

inline void ranked_frontier_append_state(
    RankedFrontierStorage &states, int trial_slot, double weight,
    uuber::BitsetState &&forced_complete_bits, bool forced_complete_bits_valid,
    uuber::BitsetState &&forced_survive_bits, bool forced_survive_bits_valid,
    TimeConstraintMap &&time_constraints,
    const SequenceStateKey *state_key = nullptr) {
  const SequenceStateKey computed_state_key =
      state_key != nullptr
          ? *state_key
          : sequence_state_key(forced_complete_bits, forced_complete_bits_valid,
                               forced_survive_bits, forced_survive_bits_valid,
                               time_constraints);
  states.trial_slot.push_back(trial_slot);
  states.weight.push_back(weight);
  states.forced_complete_bits.push_back(std::move(forced_complete_bits));
  states.forced_survive_bits.push_back(std::move(forced_survive_bits));
  states.forced_complete_bits_valid.push_back(forced_complete_bits_valid ? 1u
                                                                         : 0u);
  states.forced_survive_bits_valid.push_back(forced_survive_bits_valid ? 1u
                                                                       : 0u);
  states.time_constraints.push_back(std::move(time_constraints));
  states.state_keys.push_back(computed_state_key);
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
        std::move(states.time_constraints[idx]), &entry.state_key);
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

inline bool ranked_append_frontier_candidate(
    int trial_slot, double weight, uuber::BitsetState &&forced_complete_bits,
    bool forced_complete_bits_valid, uuber::BitsetState &&forced_survive_bits,
    bool forced_survive_bits_valid, TimeConstraintMap &&time_constraints,
    RankedFrontierStorage &next_states,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order) {
  if (next_state_lookup != nullptr) {
    RankedFrontierCollapseKey collapse_key;
    collapse_key.trial_slot = trial_slot;
    collapse_key.state_key =
        sequence_state_key(forced_complete_bits, forced_complete_bits_valid,
                           forced_survive_bits, forced_survive_bits_valid,
                           time_constraints);
    const auto existing_it = next_state_lookup->find(collapse_key);
    if (existing_it != next_state_lookup->end()) {
      next_states.weight[existing_it->second] += weight;
      return true;
    }
    const std::size_t next_state_index = next_states.size();
    ranked_frontier_append_state(
        next_states, trial_slot, weight, std::move(forced_complete_bits),
        forced_complete_bits_valid, std::move(forced_survive_bits),
        forced_survive_bits_valid, std::move(time_constraints),
        &collapse_key.state_key);
    next_state_lookup->emplace(std::move(collapse_key), next_state_index);
    return true;
  }

  if (next_state_collapse_order == nullptr) {
    return false;
  }
  const std::size_t next_state_index = next_states.size();
  RankedStateCollapseIndex collapse_index;
  collapse_index.candidate_index = next_state_index;
  collapse_index.trial_slot = trial_slot;
  collapse_index.state_key =
      sequence_state_key(forced_complete_bits, forced_complete_bits_valid,
                         forced_survive_bits, forced_survive_bits_valid,
                         time_constraints);
  ranked_frontier_append_state(
      next_states, trial_slot, weight, std::move(forced_complete_bits),
      forced_complete_bits_valid, std::move(forced_survive_bits),
      forced_survive_bits_valid, std::move(time_constraints),
      &collapse_index.state_key);
  next_state_collapse_order->push_back(std::move(collapse_index));
  return true;
}

template <typename PointVec>
inline bool ranked_reduce_scenario_weight_view_into_frontier(
    const uuber::NativeContext &ctx, const PointVec &points,
    const RankedFrontierStorage &states,
    const std::vector<std::size_t> &source_state_indices,
    const std::vector<double> &scenario_weights,
    const std::vector<int> &persistent_sources,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> *scenario_survival,
    RankedFrontierStorage &next_states,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order,
    double branch_eps);

template <typename PointVec>
inline bool ranked_reduce_exact_batch_into_frontier(
    const uuber::NativeContext &ctx, const PointVec &points,
    const RankedFrontierStorage &states,
    const std::vector<std::size_t> &source_state_indices,
    const std::vector<double> &scenario_weights, const std::vector<double> &denom,
    const std::vector<int> &persistent_sources,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> *scenario_survival,
    RankedFrontierStorage &next_states,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order,
    double branch_eps) {
  if (points.size() != scenario_weights.size()) {
    return false;
  }

  std::vector<double> normalized_weights(points.size(), 0.0);
  for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
    if (!uuber::vector_lane_active(active_mask, point_idx) ||
        point_idx >= points.density_index.size()) {
      continue;
    }
    const std::size_t lane_idx = points.density_index[point_idx];
    const double denom_val = (lane_idx < denom.size()) ? denom[lane_idx] : 0.0;
    if (!std::isfinite(denom_val) || denom_val <= branch_eps) {
      continue;
    }
    const double normalized_weight = scenario_weights[point_idx] / denom_val;
    if (!std::isfinite(normalized_weight) || normalized_weight <= branch_eps) {
      continue;
    }
    normalized_weights[point_idx] = normalized_weight;
  }

  return ranked_reduce_scenario_weight_view_into_frontier(
      ctx, points, states, source_state_indices, normalized_weights,
      persistent_sources, active_mask, scenario_survival, next_states,
      next_state_lookup, next_state_collapse_order, branch_eps);
}

template <typename PointVec>
inline bool ranked_reduce_scenario_weight_view_into_frontier(
    const uuber::NativeContext &ctx, const PointVec &points,
    const RankedFrontierStorage &states,
    const std::vector<std::size_t> &source_state_indices,
    const std::vector<double> &scenario_weights,
    const std::vector<int> &persistent_sources,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> *scenario_survival,
    RankedFrontierStorage &next_states,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order,
    double branch_eps) {
  if (points.size() != scenario_weights.size()) {
    return false;
  }

  for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
    if (!uuber::vector_lane_active(active_mask, point_idx)) {
      continue;
    }
    const uuber::VectorLaneRef point = uuber::vector_lane_ref(points, point_idx);
    if (point.density_index >= source_state_indices.size()) {
      return false;
    }
    const std::size_t lane_idx = point.density_index;
    const std::size_t source_state_index = source_state_indices[lane_idx];
    const double source_weight = scenario_weights[point_idx];
    if (!std::isfinite(source_weight) || source_weight <= branch_eps) {
      continue;
    }
    double weight = states.weight[source_state_index] * source_weight;
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

    uuber::BitsetState forced_complete_bits;
    if (point.forced_complete_bits != nullptr) {
      forced_complete_bits = *point.forced_complete_bits;
    }
    const bool forced_complete_bits_valid = point.forced_complete_bits_valid;
    uuber::BitsetState forced_survive_bits;
    if (point.forced_survive_bits != nullptr) {
      forced_survive_bits = *point.forced_survive_bits;
    }
    bool forced_survive_bits_valid = point.forced_survive_bits_valid;
    TimeConstraintMap time_constraints;
    if (point.time_constraints != nullptr) {
      time_constraints = *point.time_constraints;
    }
    const double t = point.t;

    if (!persistent_sources.empty() &&
        !ranked_apply_persistent_survive_sources(
            ctx, persistent_sources, t, forced_survive_bits,
            forced_survive_bits_valid, time_constraints)) {
      continue;
    }

    const int trial_slot = states.trial_slot[source_state_index];
    if (trial_slot < 0) {
      return false;
    }
    if (!ranked_append_frontier_candidate(
            trial_slot, weight, std::move(forced_complete_bits),
            forced_complete_bits_valid, std::move(forced_survive_bits),
            forced_survive_bits_valid, std::move(time_constraints), next_states,
            next_state_lookup, next_state_collapse_order)) {
      return false;
    }
  }

  return true;
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

inline bool ranked_batch_plan_supports_deterministic_fastpath_impl(
    const RankedBatchPlan &plan) {
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

inline void ranked_finalize_batch_plan(const uuber::NativeContext &ctx,
                                       RankedBatchPlan &plan) {
  plan.slice.node_indices.clear();
  plan.slice.relevant_source_ids.clear();
  for (const RankedProgramEvalRef &eval_ref : plan.slice.evals) {
    plan.slice.node_indices.push_back(eval_ref.node_idx);
    std::vector<int> eval_sources =
        relevant_source_ids_for_batch_node(ctx, eval_ref.node_idx);
    plan.slice.relevant_source_ids.insert(plan.slice.relevant_source_ids.end(),
                                          eval_sources.begin(),
                                          eval_sources.end());
  }
  sort_unique(plan.slice.relevant_source_ids);
  plan.valid = !plan.slice.empty();
  plan.deterministic =
      ranked_batch_plan_supports_deterministic_fastpath_impl(plan);
}

inline bool ranked_state_delta_same(const RankedStateDelta &lhs,
                                    const RankedStateDelta &rhs) noexcept {
  return lhs.kind == rhs.kind && lhs.trigger_bit == rhs.trigger_bit &&
         lhs.source_mask_begin == rhs.source_mask_begin &&
         lhs.source_mask_count == rhs.source_mask_count &&
         lhs.invalidate_slot == rhs.invalidate_slot &&
         lhs.bind_exact_current_time == rhs.bind_exact_current_time &&
         lhs.source_ids == rhs.source_ids && lhs.source_bits == rhs.source_bits;
}

inline bool ranked_batch_plan_same(const RankedBatchPlan &lhs,
                                   const RankedBatchPlan &rhs) noexcept {
  if (lhs.valid != rhs.valid || lhs.deterministic != rhs.deterministic ||
      !ranked_program_slice_same(lhs.slice, rhs.slice) ||
      lhs.deltas.size() != rhs.deltas.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.deltas.size(); ++i) {
    if (!ranked_state_delta_same(lhs.deltas[i], rhs.deltas[i])) {
      return false;
    }
  }
  return true;
}

inline void ranked_dedupe_batch_plans(std::vector<RankedBatchPlan> &batch_plans) {
  if (batch_plans.size() < 2u) {
    return;
  }
  std::vector<RankedBatchPlan> unique_plans;
  unique_plans.reserve(batch_plans.size());
  for (RankedBatchPlan &plan : batch_plans) {
    bool duplicate = false;
    for (const RankedBatchPlan &existing : unique_plans) {
      if (ranked_batch_plan_same(existing, plan)) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) {
      unique_plans.push_back(std::move(plan));
    }
  }
  batch_plans.swap(unique_plans);
}

inline void ranked_group_batch_plans(RankedNodeBatchPlan &plan) {
  plan.batch_plan_groups.clear();
  plan.has_shared_slice_groups = false;
  if (plan.batch_plans.empty()) {
    return;
  }
  plan.batch_plan_groups.reserve(plan.batch_plans.size());
  for (std::size_t plan_idx = 0; plan_idx < plan.batch_plans.size(); ++plan_idx) {
    const RankedBatchPlan &batch_plan = plan.batch_plans[plan_idx];
    bool grouped = false;
    for (RankedBatchPlanGroup &group : plan.batch_plan_groups) {
      const RankedBatchPlan &leader_plan = plan.batch_plans[group.leader_index];
      if (!ranked_program_slice_same(leader_plan.slice, batch_plan.slice)) {
        continue;
      }
      group.plan_indices.push_back(plan_idx);
      group.has_deltas = group.has_deltas || !batch_plan.deltas.empty();
      plan.has_shared_slice_groups = true;
      grouped = true;
      break;
    }
    if (grouped) {
      continue;
    }
    RankedBatchPlanGroup group;
    group.leader_index = plan_idx;
    group.plan_indices.push_back(plan_idx);
    group.has_deltas = !batch_plan.deltas.empty();
    plan.batch_plan_groups.push_back(std::move(group));
  }
}

} // namespace

bool ranked_batch_plan_supports_deterministic_fastpath(
    const RankedBatchPlan &plan) {
  return ranked_batch_plan_supports_deterministic_fastpath_impl(plan);
}

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
  compile_ranked_node_batch_plan(
      ctx, node_idx,
      [this](int child_node_idx) -> const RankedNodeBatchPlan & {
        return plan_for_node(child_node_idx);
      },
      plan);
  plan.compiled = true;
  plan.compiling = false;
}

void compile_ranked_node_batch_plan(
    const uuber::NativeContext &ctx, int node_idx,
    const std::function<const RankedNodeBatchPlan &(int)> &plan_lookup,
    RankedNodeBatchPlan &plan) {
  plan.batch_plans.clear();
  plan.batch_plan_groups.clear();
  plan.valid = false;
  plan.has_shared_slice_groups = false;

  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.ir.nodes.size())) {
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
    ranked_finalize_batch_plan(ctx, batch_plan);
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
        const RankedNodeBatchPlan &child_plan = plan_lookup(child_ids[idx]);
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
          ranked_finalize_batch_plan(ctx, batch_plan);
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
        const RankedNodeBatchPlan &child_plan = plan_lookup(child_ids[idx]);
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
            delta.source_bits =
                ranked_source_bits_for_ids(ctx, delta.source_ids);
            batch_plan.deltas.push_back(std::move(delta));
          }
          ranked_finalize_batch_plan(ctx, batch_plan);
          plan.batch_plans.push_back(std::move(batch_plan));
        }
      }
    }
    break;
  }
  case uuber::IrNodeOp::Guard: {
    if (node.reference_idx >= 0 && node.blocker_idx >= 0) {
      const RankedNodeBatchPlan &ref_plan = plan_lookup(node.reference_idx);
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
            delta.source_bits =
                ranked_source_bits_for_ids(ctx, delta.source_ids);
            ranked_attach_delta_for_node(ctx, node.blocker_idx, delta);
            batch_plan.deltas.push_back(std::move(delta));
          }
          ranked_finalize_batch_plan(ctx, batch_plan);
          plan.batch_plans.push_back(std::move(batch_plan));
        }
      }
    }
    break;
  }
  default:
    break;
  }

  ranked_dedupe_batch_plans(plan.batch_plans);
  ranked_group_batch_plans(plan);
  plan.valid = !plan.batch_plans.empty();
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

inline bool ranked_build_active_frontier(
    const RankedFrontierStorage &states,
    const std::vector<const double *> &times_by_trial, std::size_t rank_idx,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_by_trial,
    const uuber::TrialParamsSoA *default_trial_params_soa, double branch_eps,
    RankedActiveFrontier &out) {
  out = RankedActiveFrontier{};
  const bool use_frontier_groups = (trial_params_soa_by_trial != nullptr);
  out.refs.reserve(states.size());

  bool frontier_uniform_trial_params_known = false;
  out.uniform_trial_params_soa = default_trial_params_soa;
  for (std::size_t state_idx = 0; state_idx < states.size(); ++state_idx) {
    if (!std::isfinite(states.weight[state_idx]) ||
        states.weight[state_idx] <= branch_eps) {
      continue;
    }
    const int trial_slot = states.trial_slot[state_idx];
    if (trial_slot < 0 ||
        static_cast<std::size_t>(trial_slot) >= times_by_trial.size()) {
      return false;
    }
    const double t = times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx];
    if (!std::isfinite(t) || t < 0.0) {
      return false;
    }
    if (rank_idx > 0u) {
      const double lower_t =
          times_by_trial[static_cast<std::size_t>(trial_slot)][rank_idx - 1u];
      if (!std::isfinite(lower_t) || lower_t < 0.0) {
        return false;
      }
    }

    const uuber::TrialParamsSoA *state_trial_params_soa =
        default_trial_params_soa;
    if (trial_params_soa_by_trial != nullptr &&
        static_cast<std::size_t>(trial_slot) < trial_params_soa_by_trial->size()) {
      state_trial_params_soa =
          (*trial_params_soa_by_trial)[static_cast<std::size_t>(trial_slot)];
    }
    RankedFrontierStateRef state_ref;
    state_ref.state_index = state_idx;
    state_ref.trial_params_soa = state_trial_params_soa;
    state_ref.state_key = ranked_frontier_state_key(states, state_idx);
    out.refs.push_back(std::move(state_ref));

    if (!frontier_uniform_trial_params_known) {
      out.uniform_trial_params_soa = state_trial_params_soa;
      frontier_uniform_trial_params_known = true;
    } else if (out.uniform_trial_params_soa != state_trial_params_soa) {
      out.uniform_trial_params_soa = nullptr;
    }
  }

  if (use_frontier_groups && out.refs.size() > 1u) {
    std::sort(out.refs.begin(), out.refs.end(), ranked_frontier_state_ref_less);
    out.groups.reserve(out.refs.size());
    std::size_t group_begin = 0u;
    while (group_begin < out.refs.size()) {
      std::size_t group_end = group_begin + 1u;
      while (group_end < out.refs.size() &&
             ranked_frontier_state_ref_same_group(out.refs[group_begin],
                                                  out.refs[group_end])) {
        ++group_end;
      }
      RankedFrontierGroup group;
      group.trial_params_soa = out.refs[group_begin].trial_params_soa;
      group.begin = group_begin;
      group.end = group_end;
      out.groups.push_back(std::move(group));
      group_begin = group_end;
    }
  }

  out.split_groups = use_frontier_groups && out.groups.size() > 1u &&
                     out.groups.size() * 2u < out.refs.size();
  return true;
}

inline bool ranked_prepare_group_denom(
    const uuber::NativeContext &ctx, int info_node_idx, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::size_t rank_idx, const uuber::TrialParamsSoA *uniform_trial_params_soa,
    RankedFrontierExecutionScratch &scratch) {
  RankedLaneGroupBatch &lane_group_batch = scratch.lane_group_batch;
  std::vector<double> &denom = scratch.denom;
  denom.assign(lane_group_batch.rank_batch.size(), 1.0);
  if (rank_idx == 0u) {
    return true;
  }
  uuber::TreeNodeBatchValues denom_values;
  if (!exact_eval_node_batch_from_batch(
          ctx, info_node_idx, lane_group_batch.denom_batch, component_idx,
          trial_params, trial_type_key, EvalNeed::kSurvival, denom_values,
          uniform_trial_params_soa) ||
      denom_values.survival.size() != denom.size()) {
    return false;
  }
  for (std::size_t i = 0; i < denom.size(); ++i) {
    denom[i] = denom_values.survival[i];
  }
  return true;
}

inline bool ranked_execute_frontier_group(
    const RankedNodeBatchPlan &transition_plan, const uuber::NativeContext &ctx,
    int info_node_idx, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const RankedFrontierStorage &states,
    const std::vector<RankedFrontierStateRef> &active_state_refs,
    const std::vector<const double *> &times_by_trial, std::size_t rank_idx,
    std::size_t group_begin, std::size_t group_end,
    const uuber::TrialParamsSoA *uniform_trial_params_soa,
    const std::vector<int> &persistent_sources,
    const uuber::CompetitorClusterCacheEntry *competitor_cache_ptr,
    RankedFrontierExecutionScratch &scratch,
    RankedFrontierStorage &next_states_collapsed,
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash> *next_state_lookup,
    std::vector<RankedStateCollapseIndex> *next_state_collapse_order,
    double branch_eps) {
  scratch.clear_materialized();
  ranked_materialize_lane_group_batch(
      states, active_state_refs, times_by_trial, rank_idx, group_begin,
      group_end, uniform_trial_params_soa, scratch.lane_group_batch);
  if (!ranked_prepare_group_denom(
          ctx, info_node_idx, component_idx, trial_params, trial_type_key,
          rank_idx, uniform_trial_params_soa, scratch)) {
    return false;
  }
  RankedFrontierReducer reducer(
      ctx, states, scratch.lane_group_batch.source_state_indices, scratch.denom,
      persistent_sources, competitor_cache_ptr, component_idx, trial_params,
      trial_type_key, uniform_trial_params_soa, scratch.scenario_survival,
      next_states_collapsed, next_state_lookup, next_state_collapse_order,
      branch_eps);
  return exact_execute_compiled_transition_plan_from_batch(
             transition_plan, ctx, scratch.lane_group_batch.rank_batch,
             component_idx, trial_params, trial_type_key, reducer,
             uniform_trial_params_soa) != ExactTransitionExecutionResult::Error;
}

inline bool ranked_build_step_execution_contexts(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices, int component_idx,
    RankedBatchPlanner &planner,
    const std::vector<const std::vector<int> *> &step_competitor_ids_ptrs,
    const std::vector<std::vector<int>> &step_persistent_sources,
    std::vector<RankedStepExecutionContext> &step_contexts) {
  const std::size_t rank_count = outcome_indices.size();
  if (node_indices.size() != rank_count ||
      step_competitor_ids_ptrs.size() != rank_count ||
      step_persistent_sources.size() != rank_count) {
    return false;
  }
  step_contexts.assign(rank_count, RankedStepExecutionContext{});
  static const std::vector<int> kEmptyCompetitorIds;
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
    RankedStepExecutionContext &step_ctx = step_contexts[rank_idx];
    const int transition_node_idx =
        exact_resolve_transition_node_idx(ctx, outcome_node_idx, component_idx);
    if (transition_node_idx != NA_INTEGER) {
      step_ctx.transition_plan = &planner.plan_for_node(transition_node_idx);
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    step_ctx.info_node_idx = resolve_dense_node_idx_required(ctx, info.node_id);
    step_ctx.persistent_sources = &step_persistent_sources[rank_idx];
    const std::vector<int> &competitors =
        step_competitor_ids_ptrs[rank_idx] != nullptr
            ? *step_competitor_ids_ptrs[rank_idx]
            : kEmptyCompetitorIds;
    if (!competitors.empty()) {
      step_ctx.competitor_cache_ptr =
          &fetch_competitor_cluster_cache(ctx, competitors);
    }
  }
  return true;
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
  if (outcome_indices.empty() || step_competitor_ids_ptrs == nullptr ||
      step_persistent_sources == nullptr) {
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
  std::vector<RankedStepExecutionContext> step_contexts;
  if (!ranked_build_step_execution_contexts(
          ctx, outcome_indices, node_indices, component_idx, *planner,
          *step_competitor_ids_ptrs, *step_persistent_sources, step_contexts)) {
    return false;
  }
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
    const RankedStepExecutionContext &step_ctx = step_contexts[rank_idx];
    if (step_ctx.transition_plan == nullptr) {
      continue;
    }

    RankedFrontierStorage next_states_collapsed;
    next_states_collapsed.reserve(states.size() * 2u);
    std::unordered_map<RankedFrontierCollapseKey, std::size_t,
                       RankedFrontierCollapseKeyHash>
        next_state_lookup;
    std::vector<RankedStateCollapseIndex> next_state_collapse_order;
    next_state_collapse_order.reserve(states.size() * 2u);

    RankedActiveFrontier active_frontier;
    if (!ranked_build_active_frontier(states, times_by_trial, rank_idx,
                                      trial_params_soa_by_trial,
                                      default_trial_params_soa, kBranchEps,
                                      active_frontier)) {
      return false;
    }
    if (active_frontier.refs.empty()) {
      return true;
    }
    RankedFrontierExecutionScratch scratch;
    scratch.reserve(active_frontier.refs.size());

    const bool use_hash_frontier_collapse = active_frontier.split_groups;
    if (use_hash_frontier_collapse) {
      next_state_lookup.reserve(states.size() * 2u);
    }
    if (active_frontier.split_groups) {
      for (const RankedFrontierGroup &group : active_frontier.groups) {
        if (group.begin >= group.end) {
          continue;
        }
        if (!ranked_execute_frontier_group(
                *step_ctx.transition_plan, ctx, step_ctx.info_node_idx, component_idx,
                trial_params, trial_type_key, states, active_frontier.refs,
                times_by_trial, rank_idx, group.begin, group.end,
                group.trial_params_soa, *step_ctx.persistent_sources,
                step_ctx.competitor_cache_ptr,
                scratch, next_states_collapsed,
                use_hash_frontier_collapse ? &next_state_lookup : nullptr,
                use_hash_frontier_collapse ? nullptr
                                           : &next_state_collapse_order,
                kBranchEps)) {
          return false;
        }
      }
    } else {
      if (!ranked_execute_frontier_group(
              *step_ctx.transition_plan, ctx, step_ctx.info_node_idx, component_idx,
              trial_params, trial_type_key, states, active_frontier.refs,
              times_by_trial, rank_idx, 0u, active_frontier.refs.size(),
              active_frontier.uniform_trial_params_soa,
              *step_ctx.persistent_sources, step_ctx.competitor_cache_ptr,
              scratch, next_states_collapsed,
              use_hash_frontier_collapse ? &next_state_lookup : nullptr,
              use_hash_frontier_collapse ? nullptr
                                         : &next_state_collapse_order,
              kBranchEps)) {
        return false;
      }
    }

    if (!use_hash_frontier_collapse) {
      collapse_ranked_state_candidates(next_states_collapsed,
                                       next_state_collapse_order, kBranchEps);
    }
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
