#pragma once

#include "context.h"

namespace uuber {

enum class EvalNeed : std::uint8_t {
  kNone = 0,
  kDensity = 1 << 0,
  kSurvival = 1 << 1,
  kCDF = 1 << 2
};

constexpr EvalNeed kEvalAll =
    static_cast<EvalNeed>(static_cast<std::uint8_t>(EvalNeed::kDensity) |
                          static_cast<std::uint8_t>(EvalNeed::kSurvival) |
                          static_cast<std::uint8_t>(EvalNeed::kCDF));

inline EvalNeed operator|(EvalNeed lhs, EvalNeed rhs) {
  return static_cast<EvalNeed>(static_cast<std::uint8_t>(lhs) |
                               static_cast<std::uint8_t>(rhs));
}

inline EvalNeed operator&(EvalNeed lhs, EvalNeed rhs) {
  return static_cast<EvalNeed>(static_cast<std::uint8_t>(lhs) &
                               static_cast<std::uint8_t>(rhs));
}

inline EvalNeed &operator|=(EvalNeed &lhs, EvalNeed rhs) {
  lhs = lhs | rhs;
  return lhs;
}

inline bool needs_density(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kDensity) != 0;
}

inline bool needs_survival(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kSurvival) != 0;
}

inline bool needs_cdf(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kCDF) != 0;
}

inline bool has_any_need(EvalNeed need) {
  return static_cast<std::uint8_t>(need) != 0;
}

struct KernelNodeValues {
  double density{0.0};
  double survival{1.0};
  double cdf{0.0};
};

struct KernelRuntimeState {
  const KernelProgram *program{nullptr};
  const KernelQueryPlan *active_plan{nullptr};
  std::vector<KernelNodeValues> slots;
  std::vector<double> child_primary;
  std::vector<double> child_density;
  std::vector<double> prefix;
  std::vector<double> suffix;
  int computed_upto{-1};
  bool initialized{false};
};

using KernelEventEvalFnPtr = KernelNodeValues (*)(const void *context,
                                                  int event_idx,
                                                  EvalNeed need);
using KernelGuardEvalFnPtr =
    KernelNodeValues (*)(const void *context, const KernelOp &op,
                         const KernelNodeValues &reference_value,
                         const KernelNodeValues &blocker_value,
                         EvalNeed need);

struct KernelEventEvaluator {
  const void *context{nullptr};
  KernelEventEvalFnPtr fn{nullptr};

  explicit operator bool() const { return fn != nullptr; }

  KernelNodeValues operator()(int event_idx, EvalNeed need) const {
    return fn ? fn(context, event_idx, need) : KernelNodeValues{};
  }
};

struct KernelGuardEvaluator {
  const void *context{nullptr};
  KernelGuardEvalFnPtr fn{nullptr};

  explicit operator bool() const { return fn != nullptr; }

  KernelNodeValues operator()(const KernelOp &op,
                              const KernelNodeValues &reference_value,
                              const KernelNodeValues &blocker_value,
                              EvalNeed need) const {
    return fn ? fn(context, op, reference_value, blocker_value, need)
              : KernelNodeValues{};
  }
};

bool eval_kernel_query_plan(const KernelProgram &program,
                            const KernelQueryPlan &plan,
                            const KernelEventEvaluator &event_eval,
                            const KernelGuardEvaluator &guard_eval,
                            std::vector<KernelNodeValues> &out_values);

bool eval_kernel_single_query_plan(const KernelProgram &program,
                                   const KernelQueryPlan &plan,
                                   const KernelEventEvaluator &event_eval,
                                   const KernelGuardEvaluator &guard_eval,
                                   KernelNodeValues &out_values);

void reset_kernel_runtime(const KernelProgram &program,
                          KernelRuntimeState &runtime);

void invalidate_kernel_runtime_from_slot(KernelRuntimeState &runtime,
                                         int slot_begin);

bool eval_kernel_query_plan_incremental(const KernelProgram &program,
                                        KernelRuntimeState &runtime,
                                        const KernelQueryPlan &plan,
                                        const KernelEventEvaluator &event_eval,
                                        const KernelGuardEvaluator &guard_eval,
                                        std::vector<KernelNodeValues> &out_values);

bool eval_kernel_single_query_plan_incremental(
    const KernelProgram &program, KernelRuntimeState &runtime,
    const KernelQueryPlan &plan, const KernelEventEvaluator &event_eval,
    const KernelGuardEvaluator &guard_eval, KernelNodeValues &out_values);

TrialParamsSoA build_base_trial_params_soa(const NativeContext &ctx);

} // namespace uuber
