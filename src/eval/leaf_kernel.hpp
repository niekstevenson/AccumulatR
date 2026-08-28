#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>

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

inline double normal_hart_tail(const double z) noexcept {
  if (z > 37.0) {
    return 0.0;
  }
  const double factor = z < kNormalHartSplit
                            ? normal_hart_ratio(z)
                            : 1.0 / (kSqrtTwoPi * normal_hart_fraction(z));
  return std::exp(-0.5 * z * z) * factor;
}

// Hart's near-double rational CDF approximation, matching EMC2.
[[gnu::always_inline]] inline double normal_cdf_fast(
    const double x) noexcept {
  if (std::isnan(x)) {
    return x;
  }
  const double z = std::fabs(x);
  const double tail = normal_hart_tail(z);
  return x <= 0.0 ? tail : 1.0 - tail;
}

[[gnu::always_inline]] inline double normal_survival_fast(
    const double x) noexcept {
  return normal_cdf_fast(-x);
}

[[gnu::always_inline]] inline double normal_log_cdf_fast(
    const double x) noexcept {
  const double cdf = normal_cdf_fast(x);
  return cdf > 0.0 ? std::log(cdf) : -1e30;
}

[[gnu::always_inline]] inline double normal_pdf_fast(
    const double x) noexcept {
  return kInverseSqrtTwoPi * std::exp(-0.5 * x * x);
}

inline double exgauss_raw_pdf(double x, double mu, double sigma,
                              double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  const double tail = normal_cdf_fast(z - sigma_over_tau);
  return inv_tau * std::exp(exponent) * tail;
}

inline double exgauss_raw_cdf(double x, double mu, double sigma,
                              double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  const double tail = normal_cdf_fast(z - sigma_over_tau);
  const double base = normal_cdf_fast(z);
  const double exp_term = std::exp(exponent);
  return clamp_probability(base - exp_term * tail);
}

inline double lba_denom(double v, double sv) noexcept {
  double denom = normal_cdf_fast(v / sv);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

[[gnu::always_inline]] inline double lba_pdf_fast(
    double x, double v, double B, double A, double sv) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) ||
      !std::isfinite(sv) || sv <= 0.0 || !std::isfinite(B) ||
      !std::isfinite(A)) {
    return 0.0;
  }
  const double denom = lba_denom(v, sv);
  if (!(denom > 0.0)) {
    return 0.0;
  }
  double pdf = 0.0;
  if (A > 1e-10) {
    const double zs = x * sv;
    if (!(zs > 0.0) || !std::isfinite(zs)) {
      return 0.0;
    }
    const double cmz = B - x * v;
    const double cz = cmz / zs;
    const double cz_max = (cmz - A) / zs;
    pdf = (v * (normal_cdf_fast(cz) - normal_cdf_fast(cz_max)) +
           sv * (normal_pdf_fast(cz_max) - normal_pdf_fast(cz))) /
          (A * denom);
  } else {
    pdf = normal_pdf_fast((B / x - v) / sv) * B /
          (sv * x * x * denom);
  }
  return pdf;
}

[[gnu::always_inline]] inline double lba_cdf_fast(
    double x, double v, double B, double A, double sv) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) ||
      !std::isfinite(sv) || sv <= 0.0 || !std::isfinite(B) ||
      !std::isfinite(A)) {
    return 0.0;
  }
  const double denom = lba_denom(v, sv);
  if (!(denom > 0.0)) {
    return 0.0;
  }
  double cdf = 0.0;
  if (A > 1e-10) {
    const double zs = x * sv;
    if (!(zs > 0.0) || !std::isfinite(zs)) {
      return 0.0;
    }
    const double cmz = B - x * v;
    const double xx = cmz - A;
    const double cz = cmz / zs;
    const double cz_max = xx / zs;
    cdf = (1.0 + (zs * (normal_pdf_fast(cz_max) -
                        normal_pdf_fast(cz)) +
                  xx * normal_cdf_fast(cz_max) -
                  cmz * normal_cdf_fast(cz)) /
                     A) /
          denom;
  } else {
    cdf = (1.0 - normal_cdf_fast((B / x - v) / sv)) / denom;
  }
  return clamp_probability(cdf);
}

inline double rdm_pigt0(double x, double k, double l) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  if (std::fabs(l) < 1e-12) {
    const double z = k / std::sqrt(x);
    return clamp_probability(2.0 * normal_survival_fast(z));
  }
  if (k == 0.0) {
    return 0.0;
  }
  const double lambda = k * k;
  const double mu = k / l;
  const double scale = std::sqrt(lambda / x);
  const double time_ratio = x / mu;
  const double z1 = scale * (1.0 + time_ratio);
  const double z2 = scale * (1.0 - time_ratio);
  return clamp_probability(
      normal_survival_fast(z1) * std::exp(2.0 * lambda / mu) +
      (1.0 - normal_cdf_fast(z2)));
}

inline double rdm_digt0(double x, double k, double l) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  const double delta = x * l - k;
  return std::exp(-0.5 * delta * delta / x) *
         (std::fabs(k) * kInverseSqrtTwoPi) / x / std::sqrt(x);
}

inline double rdm_pigt(double x, double k, double l, double a,
                       double threshold = 1e-4) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return rdm_pigt0(x, k, l);
  }
  const double sqt = std::sqrt(x);
  const double inv_sqt = 1.0 / sqt;
  const double inv_x = 1.0 / x;
  double cdf = 0.0;
  if (l < threshold) {
    const double plus = k + a;
    const double minus = k - a;
    const double t5a = 2.0 * normal_cdf_fast(plus * inv_sqt) - 1.0;
    const double t5b = 2.0 * normal_cdf_fast(-plus * inv_sqt) - 1.0;
    const double scale = kSqrtTwoPi * inv_sqt / a;
    const double t6a = scale * std::exp(-0.5 * plus * plus * inv_x);
    const double t6b = scale * std::exp(-0.5 * minus * minus * inv_x);
    cdf = 1.0 + t6a - t6b +
          ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
  } else {
    const double delta_a = k - a - x * l;
    const double delta_b = a + k - x * l;
    const double t1a = std::exp(-0.5 * delta_a * delta_a * inv_x);
    const double t1b = std::exp(-0.5 * delta_b * delta_b * inv_x);
    const double t1 = sqt * kInverseSqrtTwoPi * (t1a - t1b);
    const double t2a =
        std::exp(2.0 * l * (k - a) +
                 normal_log_cdf_fast(-(k - a + x * l) * inv_sqt));
    const double t2b =
        std::exp(2.0 * l * (k + a) +
                 normal_log_cdf_fast(-(k + a + x * l) * inv_sqt));
    const double t2 = a + (t2b - t2a) / (2.0 * l);
    const double t4a =
        2.0 * normal_cdf_fast((k + a) * inv_sqt - sqt * l) - 1.0;
    const double t4b =
        2.0 * normal_cdf_fast((k - a) * inv_sqt - sqt * l) - 1.0;
    const double t4 = 0.5 * (x * l - a - k + 0.5 / l) * t4a +
                      0.5 * (k - a - x * l - 0.5 / l) * t4b;
    cdf = 0.5 * (t4 + t2 + t1) / a;
  }
  return clamp_probability(cdf);
}

inline double rdm_digt(double x, double k, double l, double a,
                       double threshold = 1e-4) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return rdm_digt0(x, k, l);
  }
  double pdf = 0.0;
  if (l < threshold) {
    const double term = std::exp(-(k - a) * (k - a) / (2.0 * x)) -
                        std::exp(-(k + a) * (k + a) / (2.0 * x));
    pdf = term * kInverseSqrtTwoPi / (2.0 * a * std::sqrt(x));
  } else {
    const double sqt = std::sqrt(x);
    const double inv_sqt = 1.0 / sqt;
    const double inv_x = 1.0 / x;
    const double delta_a = a - k + x * l;
    const double delta_b = a + k - x * l;
    const double t1a = -0.5 * delta_a * delta_a * inv_x;
    const double t1b = -0.5 * delta_b * delta_b * inv_x;
    const double t1 = kInverseSqrtTwoPi *
                      (std::exp(t1a) - std::exp(t1b)) * inv_sqt;
    const double t2a =
        2.0 * normal_cdf_fast((-k + a) * inv_sqt + sqt * l) - 1.0;
    const double t2b =
        2.0 * normal_cdf_fast((k + a) * inv_sqt - sqt * l) - 1.0;
    const double t2 = 0.5 * l * (t2a + t2b);
    pdf = (t1 + t2) / (2.0 * a);
  }
  return safe_density(pdf);
}

constexpr double kRdmAEpsilon = 1e-4;
constexpr double kRdmLEpsilon = 1e-4;
constexpr double kRdmKMaximum = 1e6;

inline double rdm_clamp_drift(const double value) noexcept {
  return value > -kRdmLEpsilon && value < kRdmLEpsilon
             ? (value >= 0.0 ? kRdmLEpsilon : -kRdmLEpsilon)
             : value;
}

inline double rdm_pdf_fast(double x, double v, double B, double A,
                           double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double inv_s = 1.0 / s;
  const double v_sc = rdm_clamp_drift(v * inv_s);
  if (!std::isfinite(v_sc) || v_sc < 0.0) {
    return 0.0;
  }
  const double B_sc = B * inv_s;
  if (A < kRdmAEpsilon) {
    return B_sc < 0.0 || B_sc > kRdmKMaximum
               ? 0.0
               : safe_density(rdm_digt0(x, B_sc, v_sc));
  }
  const double a = std::max(kRdmAEpsilon, 0.5 * A * inv_s);
  return rdm_digt(x, B_sc + a, v_sc, a);
}

inline double rdm_cdf_fast(double x, double v, double B, double A,
                           double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double inv_s = 1.0 / s;
  const double v_sc = rdm_clamp_drift(v * inv_s);
  if (!std::isfinite(v_sc) || v_sc < 0.0) {
    return 0.0;
  }
  const double B_sc = B * inv_s;
  if (A < kRdmAEpsilon) {
    return B_sc < 0.0 || B_sc > kRdmKMaximum
               ? 0.0
               : clamp_probability(rdm_pigt0(x, B_sc, v_sc));
  }
  const double a = std::max(kRdmAEpsilon, 0.5 * A * inv_s);
  return rdm_pigt(x, B_sc + a, v_sc, a);
}

constexpr std::uint8_t kLeafChannelPdf = 1U;
constexpr std::uint8_t kLeafChannelCdf = 2U;
constexpr std::uint8_t kLeafChannelSurvival = 4U;
constexpr std::uint8_t kLeafChannelAll =
    kLeafChannelPdf | kLeafChannelCdf | kLeafChannelSurvival;

} // namespace detail
} // namespace accumulatr::eval
