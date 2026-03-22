#pragma once

#include <string>
#include <vector>

#include "evaluator_internal.h"

class RankedTransitionCompiler;

struct ExactScenarioPoint {
  double t{0.0};
  std::size_t density_index{0u};
  double weight{0.0};
  const uuber::TrialParamsSoA *trial_params_soa{nullptr};
  uuber::BitsetState forced_complete_bits;
  uuber::BitsetState forced_survive_bits;
  bool forced_complete_bits_valid{false};
  bool forced_survive_bits_valid{false};
  TimeConstraintMap time_constraints;
};

struct ExactScenarioBatch {
  std::vector<double> t;
  std::vector<std::size_t> density_index;
  std::vector<double> weight;
  std::vector<const uuber::TrialParamsSoA *> trial_params_soa;
  std::vector<uuber::BitsetState> forced_complete_bits;
  std::vector<uuber::BitsetState> forced_survive_bits;
  std::vector<std::uint8_t> forced_complete_bits_valid;
  std::vector<std::uint8_t> forced_survive_bits_valid;
  std::vector<TimeConstraintMap> time_constraints;

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
  }

  std::size_t size() const noexcept { return t.size(); }
  bool empty() const noexcept { return t.empty(); }

  void swap(ExactScenarioBatch &other) noexcept {
    t.swap(other.t);
    density_index.swap(other.density_index);
    weight.swap(other.weight);
    trial_params_soa.swap(other.trial_params_soa);
    forced_complete_bits.swap(other.forced_complete_bits);
    forced_survive_bits.swap(other.forced_survive_bits);
    forced_complete_bits_valid.swap(other.forced_complete_bits_valid);
    forced_survive_bits_valid.swap(other.forced_survive_bits_valid);
    time_constraints.swap(other.time_constraints);
  }
};

bool exact_eval_node_batch_from_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_outcome_density_batch_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const std::vector<double> &times, std::vector<double> &density_out,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);

bool exact_collect_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &scenario_points,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_collect_scenarios_batch_aligned_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, ExactScenarioBatch &scenario_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_collect_deterministic_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &aligned_points,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_collect_deterministic_scenarios_batch_aligned_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, ExactScenarioBatch &aligned_batch,
    std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const std::vector<ExactScenarioPoint> &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const ExactScenarioBatch &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);
