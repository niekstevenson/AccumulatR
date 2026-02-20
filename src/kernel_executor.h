#pragma once

#include <functional>

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

using KernelEventEvalFn = std::function<KernelNodeValues(int event_idx)>;
using KernelGuardEvalFn =
    std::function<KernelNodeValues(const KernelOp &op,
                                   const KernelNodeValues &reference_value,
                                   const KernelNodeValues &blocker_value,
                                   const KernelEvalNeed &need)>;

bool eval_kernel_node(const KernelProgram &program, int target_node_idx,
                      const KernelEvalNeed &need,
                      const KernelEventEvalFn &event_eval,
                      const KernelGuardEvalFn &guard_eval,
                      KernelNodeValues &out_values);

void reset_kernel_runtime(const KernelProgram &program,
                          KernelRuntimeState &runtime);

void invalidate_kernel_runtime_from_slot(KernelRuntimeState &runtime,
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

TrialParamsSoA build_base_trial_params_soa(const NativeContext &ctx);

} // namespace uuber
