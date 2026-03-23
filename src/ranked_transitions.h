#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "context.h"
#include "evaluator_internal.h"
#include "trial_params.h"

enum class RankedProgramEvalKind : std::uint8_t {
  Density = 0,
  CDF = 1,
  Survival = 2
};

struct RankedProgramEvalRef {
  RankedProgramEvalKind kind{RankedProgramEvalKind::Density};
  int node_idx{-1};
  int slot{-1};
};

struct RankedProgramSliceRef {
  std::vector<RankedProgramEvalRef> evals;

  bool empty() const noexcept { return evals.empty(); }
};

enum class RankedStateDeltaKind : std::uint8_t {
  CompleteSources = 0,
  SurviveSources = 1,
  OrWitnessSurvive = 2
};

struct RankedStateDelta {
  RankedStateDeltaKind kind{RankedStateDeltaKind::CompleteSources};
  int trigger_bit{-1};
  int source_mask_begin{-1};
  int source_mask_count{0};
  int invalidate_slot{0};
  bool bind_exact_current_time{false};
  std::vector<int> source_ids;
  std::vector<int> source_bits;
};

struct RankedBatchPlan {
  std::vector<RankedStateDelta> deltas;
  RankedProgramSliceRef slice;
  bool deterministic{false};
  bool valid{false};
};

struct RankedNodeBatchPlan {
  bool compiled{false};
  bool compiling{false};
  bool valid{false};
  std::vector<RankedBatchPlan> batch_plans;
};

class RankedBatchPlanner {
public:
  explicit RankedBatchPlanner(const uuber::NativeContext &ctx);

  const RankedNodeBatchPlan &plan_for_node(int node_idx);

private:
  void compile_node(int node_idx);

  const uuber::NativeContext &ctx;
  std::vector<RankedNodeBatchPlan> plans;
};

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids);

bool sequence_prefix_density_batch_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices,
    const std::vector<const double *> &times_by_trial, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    std::vector<double> &density_out,
    const std::vector<const uuber::TrialParamsSoA *> *trial_params_soa_by_trial =
        nullptr,
    RankedBatchPlanner *transition_planner = nullptr,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs =
        nullptr,
    const std::vector<std::vector<int>> *step_persistent_sources = nullptr);
