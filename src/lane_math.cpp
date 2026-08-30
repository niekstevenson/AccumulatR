#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#undef COMPLEX
#endif

#include "eval/lane_math.hpp"

#include <cmath>

namespace accumulatr::eval::detail {

#ifdef __APPLE__
constexpr std::size_t kVForceMinimumValues = 8U;
#endif

void exp_lanes(const double *input,
               double *output,
               const std::size_t size) noexcept {
#ifdef __APPLE__
  if (size >= kVForceMinimumValues) {
    const int lane_count = static_cast<int>(size);
    vvexp(output, input, &lane_count);
    return;
  }
#endif
  for (std::size_t lane = 0; lane < size; ++lane) {
    output[lane] = std::exp(input[lane]);
  }
}

void log_lanes(const double *input,
               double *output,
               const std::size_t size) noexcept {
#ifdef __APPLE__
  if (size >= kVForceMinimumValues) {
    const int lane_count = static_cast<int>(size);
    vvlog(output, input, &lane_count);
    return;
  }
#endif
  for (std::size_t lane = 0; lane < size; ++lane) {
    output[lane] = std::log(input[lane]);
  }
}

} // namespace accumulatr::eval::detail
