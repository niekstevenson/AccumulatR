#pragma once

#include <functional>
#include <vector>

#include "context.h"

namespace uuber {

struct TreeEvalNeed {
  bool density{false};
  bool survival{false};
  bool cdf{false};
};

struct TreeNodeBatchValues {
  std::vector<double> density;
  std::vector<double> survival;
  std::vector<double> cdf;
};

using TreeEventBatchEvalFn =
    std::function<bool(int event_idx, const std::vector<double> &times,
                       const TreeEvalNeed &need,
                       TreeNodeBatchValues &out_values)>;

using TreeGuardBatchEvalFn = std::function<bool(
    const TreeEvalOp &op, const std::vector<double> &times,
    const TreeNodeBatchValues &reference_values,
    const TreeNodeBatchValues &blocker_values, const TreeEvalNeed &need,
    TreeNodeBatchValues &out_values)>;

} // namespace uuber
