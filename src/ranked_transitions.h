#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "context.h"
#include "trial_params.h"

enum class RankedTransitionStepKind : std::uint8_t {
  EvalDensityNode = 0,
  EvalCDFNode = 1,
  EvalGuardEffective = 2,
  AddCompleteSources = 3,
  AddSurviveSources = 4,
  AddOrWitnessFromSources = 5
};

struct RankedTransitionStep {
  RankedTransitionStepKind kind{RankedTransitionStepKind::EvalDensityNode};
  int node_idx{-1};
  int source_mask_begin{-1};
  int source_mask_count{0};
  int invalidate_slot{0};
  bool source_mask_covers_ids{false};
  std::vector<int> source_ids;
  std::vector<int> source_bits;
};

struct RankedTransitionTemplate {
  std::vector<RankedTransitionStep> steps;
};

struct RankedNodeTransitionPlan {
  bool compiled{false};
  bool compiling{false};
  bool valid{false};
  std::vector<RankedTransitionTemplate> transitions;
};

class RankedTransitionCompiler {
public:
  explicit RankedTransitionCompiler(const uuber::NativeContext &ctx);

  const RankedNodeTransitionPlan &plan_for_node(int node_idx);

private:
  void compile_node(int node_idx);

  const uuber::NativeContext &ctx;
  std::vector<RankedNodeTransitionPlan> plans;
};

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids);

double sequence_prefix_density_resolved(
    const uuber::NativeContext &ctx, const std::vector<int> &outcome_indices,
    const std::vector<int> &node_indices, const double *times,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key,
    uuber::KernelRuntimeState *kernel_runtime = nullptr,
    RankedTransitionCompiler *transition_compiler = nullptr,
    const std::vector<const std::vector<int> *> *step_competitor_ids_ptrs =
        nullptr,
    const std::vector<std::vector<int>> *step_persistent_sources = nullptr);
