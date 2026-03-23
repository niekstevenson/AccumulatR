#pragma once

#include <functional>
#include <vector>

#include "context.h"

namespace uuber {

struct KernelEvalNeed {
  bool density{false};
  bool survival{false};
  bool cdf{false};
};

struct KernelNodeBatchValues {
  std::vector<double> density;
  std::vector<double> survival;
  std::vector<double> cdf;
};

using KernelEventBatchEvalFn =
    std::function<bool(int event_idx, const std::vector<double> &times,
                       const KernelEvalNeed &need,
                       KernelNodeBatchValues &out_values)>;

using KernelGuardBatchEvalFn = std::function<bool(
    const KernelOp &op, const std::vector<double> &times,
    const KernelNodeBatchValues &reference_values,
    const KernelNodeBatchValues &blocker_values, const KernelEvalNeed &need,
    KernelNodeBatchValues &out_values)>;

} // namespace uuber
