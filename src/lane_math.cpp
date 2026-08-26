#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#undef COMPLEX
#endif

#include "eval/lane_math.hpp"

#include <cmath>

namespace accumulatr::eval::detail {

void log_lanes(const double *input,
               double *output,
               const std::size_t size) noexcept {
#ifdef __APPLE__
  const int lane_count = static_cast<int>(size);
  vvlog(output, input, &lane_count);
#else
  for (std::size_t lane = 0; lane < size; ++lane) {
    output[lane] = std::log(input[lane]);
  }
#endif
}

} // namespace accumulatr::eval::detail
