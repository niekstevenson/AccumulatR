#pragma once

#include <vector>

namespace uuber {

struct TimeBatch {
  std::vector<double> nodes;
  std::vector<double> weights;
  double lower{0.0};
  double upper{0.0};
};

TimeBatch build_time_batch(double lower, double upper);
TimeBatch build_time_batch_0_to_upper(double upper);

} // namespace uuber
