#pragma once

#include <string>
#include <vector>

#include "evaluator_internal.h"
#include "vector_runtime.h"

class RankedBatchPlanner;

using ExactScenarioBatch = uuber::VectorLaneBatch;

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

bool exact_collect_scenarios_batch_aligned_from_batch(
    RankedBatchPlanner &planner, const uuber::NativeContext &ctx,
    int node_idx, const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ExactScenarioBatch &scenario_batch,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

bool exact_collect_deterministic_scenarios_batch_aligned_from_batch(
    RankedBatchPlanner &planner, const uuber::NativeContext &ctx,
    int node_idx, const ExactScenarioBatch &seed_batch, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    ExactScenarioBatch &aligned_batch, std::vector<std::uint8_t> &active_mask,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);

void exact_competitor_survival_batch(
    const uuber::NativeContext &ctx,
    const uuber::CompetitorClusterCacheEntry &competitor_cache,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const ExactScenarioBatch &points,
    std::vector<double> &survival_out,
    const uuber::TrialParamsSoA *uniform_trial_params_soa = nullptr);
