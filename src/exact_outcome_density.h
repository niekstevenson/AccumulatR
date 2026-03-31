#pragma once

#include <functional>
#include <string>
#include <vector>

#include "vector_runtime.h"

namespace uuber {

struct ExactScenarioLaneViewBatch {
  std::vector<double> t;
  std::vector<std::size_t> density_index;
  std::vector<double> weight;
  std::vector<const TrialParamsSoA *> trial_params_soa;
  std::vector<const BitsetState *> forced_complete_bits;
  std::vector<const BitsetState *> forced_survive_bits;
  std::vector<std::uint8_t> forced_complete_bits_valid;
  std::vector<std::uint8_t> forced_survive_bits_valid;
  std::vector<const TimeConstraintMap *> time_constraints;
  std::vector<TimeConstraintMap> mutable_time_constraints;
  std::vector<BitsetState> mutable_forced_survive_bits;

  void clear() {
    t.clear();
    density_index.clear();
    weight.clear();
    trial_params_soa.clear();
    forced_complete_bits.clear();
    forced_survive_bits.clear();
    forced_complete_bits_valid.clear();
    forced_survive_bits_valid.clear();
    time_constraints.clear();
    mutable_time_constraints.clear();
    mutable_forced_survive_bits.clear();
  }

  void reserve(std::size_t n) {
    t.reserve(n);
    density_index.reserve(n);
    weight.reserve(n);
    trial_params_soa.reserve(n);
    forced_complete_bits.reserve(n);
    forced_survive_bits.reserve(n);
    forced_complete_bits_valid.reserve(n);
    forced_survive_bits_valid.reserve(n);
    time_constraints.reserve(n);
    mutable_time_constraints.reserve(n);
    mutable_forced_survive_bits.reserve(n);
  }

  std::size_t size() const noexcept { return t.size(); }
  bool empty() const noexcept { return t.empty(); }
};

inline VectorLaneRef vector_lane_ref(const ExactScenarioLaneViewBatch &batch,
                                     std::size_t idx) {
  return VectorLaneRef{
      batch.t[idx],
      batch.density_index[idx],
      batch.weight[idx],
      idx < batch.trial_params_soa.size() ? batch.trial_params_soa[idx]
                                          : nullptr,
      idx < batch.forced_complete_bits.size() ? batch.forced_complete_bits[idx]
                                              : nullptr,
      idx < batch.forced_survive_bits.size() ? batch.forced_survive_bits[idx]
                                             : nullptr,
      idx < batch.forced_complete_bits_valid.size() &&
          batch.forced_complete_bits_valid[idx] != 0u,
      idx < batch.forced_survive_bits_valid.size() &&
          batch.forced_survive_bits_valid[idx] != 0u,
      idx < batch.time_constraints.size() ? batch.time_constraints[idx]
                                          : nullptr};
}

} // namespace uuber

#include "evaluator_internal.h"

class RankedBatchPlanner;

using ExactScenarioBatch = uuber::VectorLaneBatch;
using ExactScenarioBatchProcessFn = std::function<bool(
    const ExactScenarioBatch &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> &weights)>;
using ExactScenarioProcessFn = std::function<bool(
    const uuber::ExactScenarioLaneViewBatch &points,
    const std::vector<std::uint8_t> *active_mask,
    const std::vector<double> &weights)>;

bool exact_eval_node_batch_from_batch(
    const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::TreeNodeBatchValues &out_values,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_outcome_density_batch_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const std::vector<double> &times, std::vector<double> &density_out,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);

bool exact_collect_scenarios_batch_process_from_batch(
    RankedBatchPlanner &planner, const uuber::NativeContext &ctx, int node_idx,
    const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    const ExactScenarioBatchProcessFn &process_exact_batch,
    const ExactScenarioProcessFn &process_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const ExactScenarioBatch &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const uuber::ExactScenarioLaneViewBatch &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);
