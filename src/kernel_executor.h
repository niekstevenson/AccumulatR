#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "context.h"

namespace uuber {

struct KernelEvalNeed {
  bool density{false};
  bool survival{false};
  bool cdf{false};
};

struct KernelNodeValues {
  double density{0.0};
  double survival{1.0};
  double cdf{0.0};
};

struct KernelNodeBatchValues {
  std::vector<double> density;
  std::vector<double> survival;
  std::vector<double> cdf;
};

struct KernelRuntimeState {
  const KernelProgram *program{nullptr};
  std::vector<KernelNodeValues> slots;
  std::vector<double> child_primary;
  std::vector<double> child_density;
  std::vector<double> prefix;
  std::vector<double> suffix;
  int computed_upto{-1};
  bool initialized{false};
};

struct KernelBatchRuntimeState {
  const KernelProgram *program{nullptr};
  std::vector<KernelNodeBatchValues> slots;
  std::vector<double> child_primary;
  std::vector<double> child_density;
  std::vector<double> prefix;
  std::vector<double> suffix;
  std::size_t point_count{0};
  int computed_upto{-1};
  bool initialized{false};
};

using KernelEventEvalFn = std::function<KernelNodeValues(int event_idx)>;
using KernelGuardEvalFn =
    std::function<KernelNodeValues(const KernelOp &op,
                                   const KernelNodeValues &reference_value,
                                   const KernelNodeValues &blocker_value,
                                   const KernelEvalNeed &need)>;
using KernelEventBatchEvalFn =
    std::function<bool(int event_idx, const std::vector<double> &times,
                       const KernelEvalNeed &need,
                       KernelNodeBatchValues &out_values)>;
using KernelGuardBatchEvalFn = std::function<bool(
    const KernelOp &op, const std::vector<double> &times,
    const KernelNodeBatchValues &reference_values,
    const KernelNodeBatchValues &blocker_values, const KernelEvalNeed &need,
    KernelNodeBatchValues &out_values)>;
using KernelBatchTransitionApplyFn =
    std::function<void(const CompetitorCompiledOp &op,
                       KernelBatchRuntimeState &runtime)>;

bool eval_kernel_node(const KernelProgram &program, int target_node_idx,
                      const KernelEvalNeed &need,
                      const KernelEventEvalFn &event_eval,
                      const KernelGuardEvalFn &guard_eval,
                      KernelNodeValues &out_values);

void reset_kernel_runtime(const KernelProgram &program,
                          KernelRuntimeState &runtime);

void reset_kernel_batch_runtime(const KernelProgram &program,
                                KernelBatchRuntimeState &runtime,
                                std::size_t point_count);

void invalidate_kernel_runtime_from_slot(KernelRuntimeState &runtime,
                                         int slot_begin);

void invalidate_kernel_batch_runtime_from_slot(KernelBatchRuntimeState &runtime,
                                               int slot_begin);

bool eval_kernel_node_incremental(const KernelProgram &program,
                                  KernelRuntimeState &runtime,
                                  int target_node_idx,
                                  const KernelEvalNeed &need,
                                  const KernelEventEvalFn &event_eval,
                                  const KernelGuardEvalFn &guard_eval,
                                  KernelNodeValues &out_values);

bool eval_kernel_nodes_incremental(const KernelProgram &program,
                                   KernelRuntimeState &runtime,
                                   const std::vector<int> &target_node_indices,
                                   const KernelEvalNeed &need,
                                   const KernelEventEvalFn &event_eval,
                                   const KernelGuardEvalFn &guard_eval,
                                   std::vector<KernelNodeValues> &out_values);

bool eval_kernel_node_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    int target_node_idx, const std::vector<double> &times,
    const KernelEvalNeed &need, const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    KernelNodeBatchValues &out_values);

bool eval_kernel_nodes_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    const std::vector<int> &target_node_indices, const std::vector<double> &times,
    const KernelEvalNeed &need, const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    std::vector<KernelNodeBatchValues> &out_values);

bool eval_kernel_competitor_product_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    const std::vector<CompetitorCompiledOp> &compiled_ops,
    const std::vector<double> &times,
    const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    const KernelBatchTransitionApplyFn &apply_transition_batch,
    std::vector<double> &survival_product_out);

TrialParamsSoA build_base_trial_params_soa(const NativeContext &ctx);

} // namespace uuber
