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
TimeBatch build_time_batch_with_observed(double lower, double upper,
                                         const std::vector<double> &observed_times);
double integrate_time_batch(const TimeBatch &batch,
                            const std::vector<double> &values);

} // namespace uuber

