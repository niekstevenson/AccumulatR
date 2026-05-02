#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace accumulatr::eval {
namespace detail {

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

constexpr double kInvSqrtTwo = 0.7071067811865475244008443621048490;
constexpr double kInvSqrtTwoPi = 0.3989422804014326779399460599343819;
constexpr double kLogSqrtTwoPi = 0.9189385332046727417803297364056176;

inline double normal_cdf_polevl(double x,
                                const double *coefficients,
                                int degree) noexcept {
  double value = coefficients[0];
  for (int i = 1; i <= degree; ++i) {
    value = value * x + coefficients[i];
  }
  return value;
}

inline double normal_cdf_p1evl(double x,
                               const double *coefficients,
                               int degree) noexcept {
  double value = x + coefficients[0];
  for (int i = 1; i < degree; ++i) {
    value = value * x + coefficients[i];
  }
  return value;
}

inline double normal_erf_small(double x) noexcept {
  static constexpr double p[] = {
      9.60497373987051638749E0,
      9.00260197203842689217E1,
      2.23200534594684319226E3,
      7.00332514112805075473E3,
      5.55923013010394962768E4};
  static constexpr double q[] = {
      3.35617141647503099647E1,
      5.21357949780152679795E2,
      4.59432382970980127987E3,
      2.26290000613890934246E4,
      4.92673942608635921086E4};
  const double z = x * x;
  return x * normal_cdf_polevl(z, p, 4) /
         normal_cdf_p1evl(z, q, 5);
}

inline double normal_erfc_positive_from_exp(double x,
                                            double exp_neg_x_sq) noexcept {
  static constexpr double p[] = {
      2.46196981473530512524E-10,
      5.64189564831068821977E-1,
      7.46321056442269912687E0,
      4.86371970985681366614E1,
      1.96520832956077098242E2,
      5.26445194995477358631E2,
      9.34528527171957607540E2,
      1.02755188689515710272E3,
      5.57535335369399327526E2};
  static constexpr double q[] = {
      1.32281951154744992508E1,
      8.67072140885989742329E1,
      3.54937778887819891062E2,
      9.75708501743205489753E2,
      1.82390916687909736289E3,
      2.24633760818710981792E3,
      1.65666309194161350182E3,
      5.57535340817727675546E2};
  static constexpr double r[] = {
      5.64189583547755073984E-1,
      1.27536670759978104416E0,
      5.01905042251180477414E0,
      6.16021097993053585195E0,
      7.40974269950448939160E0,
      2.97886665372100240670E0};
  static constexpr double s[] = {
      2.26052863220117276590E0,
      9.39603524938001434673E0,
      1.20489539808096656605E1,
      1.70814450747565897222E1,
      9.60896809063285878198E0,
      3.36907645100081516050E0};
  if (x < 1.0) {
    return 1.0 - normal_erf_small(x);
  }
  const double numerator =
      x < 8.0 ? normal_cdf_polevl(x, p, 8)
              : normal_cdf_polevl(x, r, 5);
  const double denominator =
      x < 8.0 ? normal_cdf_p1evl(x, q, 8)
              : normal_cdf_p1evl(x, s, 6);
  return exp_neg_x_sq * numerator / denominator;
}

inline double normal_cdf_rational_fast(double z) noexcept {
  if (!std::isfinite(z)) {
    return z > 0.0 ? 1.0 : 0.0;
  }
  const double x = z * kInvSqrtTwo;
  const double abs_x = std::abs(x);
  if (abs_x < kInvSqrtTwo) {
    return clamp_probability(0.5 + 0.5 * normal_erf_small(x));
  }
  double tail =
      0.5 * normal_erfc_positive_from_exp(abs_x, std::exp(-abs_x * abs_x));
  if (x > 0.0) {
    tail = 1.0 - tail;
  }
  return clamp_probability(tail);
}

constexpr double kLogPi = 1.1447298858494001741434;

constexpr std::uint8_t kLeafChannelPdf = 1U;
constexpr std::uint8_t kLeafChannelCdf = 2U;
constexpr std::uint8_t kLeafChannelSurvival = 4U;
constexpr std::uint8_t kLeafChannelAll =
    kLeafChannelPdf | kLeafChannelCdf | kLeafChannelSurvival;

} // namespace detail
} // namespace accumulatr::eval
