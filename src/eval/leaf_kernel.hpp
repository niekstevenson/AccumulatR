#pragma once

#include <cmath>
#include <cstdint>

namespace accumulatr::eval {
namespace detail {

constexpr double kInverseSqrtTwoPi = 0.39894228040143267794;
constexpr double kSqrtTwoPi = 2.5066282746310005024;

inline double clamp_probability(double x) noexcept {
  if (!std::isfinite(x)) {
    return 0.0;
  }
  if (x <= 0.0) {
    return 0.0;
  }
  if (x >= 1.0) {
    return 1.0;
  }
  return x;
}

inline double safe_density(double value) noexcept {
  return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

constexpr double kNormalHartSplit = 7.07106781186547;

inline double normal_hart_ratio(const double z) noexcept {
  const double numerator =
      (((((3.52624965998911e-02 * z + 0.700383064443688) * z +
          6.37396220353165) * z + 33.912866078383) * z +
        112.079291497871) * z + 221.213596169931) * z +
      220.206867912376;
  const double denominator =
      ((((((8.83883476483184e-02 * z + 1.75566716318264) * z +
           16.064177579207) * z + 86.7807322029461) * z +
         296.564248779674) * z + 637.333633378831) * z +
       793.826512519948) * z + 440.413735824752;
  return numerator / denominator;
}

inline double normal_hart_fraction(const double z) noexcept {
  return z +
         1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
}

inline double normal_tail_factor(const double z) noexcept {
  return z < kNormalHartSplit
             ? normal_hart_ratio(z)
             : 1.0 / (kSqrtTwoPi * normal_hart_fraction(z));
}

inline double prepare_normal_cdf(const double value) noexcept {
  const double z = std::fabs(value);
  const double factor = z > 37.0 ? 0.0 : normal_tail_factor(z);
  return std::copysign(factor, value <= 0.0 ? 1.0 : -1.0);
}

inline double finish_normal_cdf(const double prepared,
                                const double exponential) noexcept {
  const double tail = std::fabs(prepared) * exponential;
  return std::signbit(prepared) ? 1.0 - tail : tail;
}

constexpr double kRdmAEpsilon = 1e-4;
constexpr double kRdmLEpsilon = 1e-4;
constexpr double kRdmKMaximum = 1e6;

inline double rdm_clamp_drift(const double value) noexcept {
  return value > -kRdmLEpsilon && value < kRdmLEpsilon
             ? (value >= 0.0 ? kRdmLEpsilon : -kRdmLEpsilon)
             : value;
}

constexpr std::uint8_t kLeafChannelPdf = 1U;
constexpr std::uint8_t kLeafChannelCdf = 2U;
constexpr std::uint8_t kLeafChannelSurvival = 4U;
constexpr std::uint8_t kLeafChannelAll =
    kLeafChannelPdf | kLeafChannelCdf | kLeafChannelSurvival;

} // namespace detail
} // namespace accumulatr::eval
