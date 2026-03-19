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

double exact_outcome_density_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);

bool exact_eval_node_batch_from_points(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<ExactScenarioPoint> &points, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    EvalNeed need, uuber::KernelNodeBatchValues &out_values);

bool exact_outcome_density_batch_from_state(
    const uuber::NativeContext &ctx, int node_idx,
    const std::vector<int> &competitor_ids, NodeEvalState &state,
    const std::vector<double> &times, std::vector<double> &density_out,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);

bool exact_collect_scenarios_batch(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<double> &times, NodeEvalState &state,
    std::vector<ExactScenarioPoint> &scenario_points);

bool exact_collect_scenarios_batch_from_points(
    RankedTransitionCompiler &compiler, const uuber::NativeContext &ctx,
    int node_idx, const std::vector<ExactScenarioPoint> &seed_points,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    std::vector<ExactScenarioPoint> &scenario_points);

bool exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    const std::vector<ExactScenarioPoint> &points,
    std::vector<double> &survival_out);
