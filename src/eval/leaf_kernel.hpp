#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "../leaf/dist_kind.hpp"
#include "quadrature.hpp"

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

[[gnu::always_inline]] inline double normal_cdf_fast(double z) noexcept {
  const double arg = -z * 0.7071067811865475244;
  return clamp_probability(0.5 * std::erfc(arg));
}

inline double log_diff_exp(double log_a, double log_b) noexcept {
  if (!std::isfinite(log_a)) {
    return R_NegInf;
  }
  if (!std::isfinite(log_b)) {
    return log_a;
  }
  if (!(log_a > log_b)) {
    return R_NegInf;
  }
  return log_a + std::log1p(-std::exp(log_b - log_a));
}

inline double log1m_exp(double log_x) noexcept {
  if (!(log_x < 0.0)) {
    return R_NegInf;
  }
  return log_x < -0.69314718055994530942
             ? std::log1p(-std::exp(log_x))
             : std::log(-std::expm1(log_x));
}

inline double exgauss_raw_log_pdf_stable(double x, double mu, double sigma,
                                         double tau) noexcept {
  const double tau_p = std::max(tau, 1e-12);
  const double sigma_p = std::max(sigma, 1e-12);
  const double y = (x - mu) / sigma_p;
  const double sigma_over_tau = sigma_p / tau_p;
  const double z = y - sigma_over_tau;
  return -std::log(tau_p) + (mu - x) / tau_p +
         (sigma_p * sigma_p) / (2.0 * tau_p * tau_p) +
         R::pnorm(z, 0.0, 1.0, 1, 1);
}

inline double exgauss_raw_log_cdf_stable(double x, double mu, double sigma,
                                         double tau) noexcept {
  const double tau_p = std::max(tau, 1e-12);
  const double sigma_p = std::max(sigma, 1e-12);
  const double y = (x - mu) / sigma_p;
  const double sigma_over_tau = sigma_p / tau_p;
  const double log_phi_1 = R::pnorm(y, 0.0, 1.0, 1, 1);
  const double log_phi_2 =
      R::pnorm(y - sigma_over_tau, 0.0, 1.0, 1, 1);
  const double log_second =
      (mu - x) / tau_p + (sigma_p * sigma_p) /
                              (2.0 * tau_p * tau_p) +
      log_phi_2;
  return log_diff_exp(log_phi_1, log_second);
}

[[gnu::always_inline]] inline bool exgauss_raw_pdf_direct(double x, double mu,
                                                     double sigma, double tau,
                                                     double *out) noexcept {
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double y = (x - mu) / sigma;
  const double z = y - sigma_over_tau;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  if (!(z > -8.0) || !(exponent > -745.0 && exponent < 700.0)) {
    return false;
  }
  const double value = inv_tau * std::exp(exponent) * normal_cdf_fast(z);
  if (!std::isfinite(value) || !(value > 0.0)) {
    return false;
  }
  *out = value;
  return true;
}

[[gnu::always_inline]] inline bool exgauss_raw_cdf_direct(double x, double mu,
                                                     double sigma, double tau,
                                                     double *out) noexcept {
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double y = (x - mu) / sigma;
  const double z = y - sigma_over_tau;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  if (!(z > -8.0) || !(exponent > -745.0 && exponent < 700.0)) {
    return false;
  }
  const double base = normal_cdf_fast(y);
  const double second = std::exp(exponent) * normal_cdf_fast(z);
  const double value = base - second;
  if (!std::isfinite(value) || value < 1e-12 || value > 1.0 - 1e-12) {
    return false;
  }
  *out = clamp_probability(value);
  return true;
}

[[gnu::always_inline]] inline bool exgauss_raw_survival_direct(double x, double mu,
                                                         double sigma,
                                                         double tau,
                                                         double *out) noexcept {
  double cdf = 0.0;
  if (!exgauss_raw_cdf_direct(x, mu, sigma, tau, &cdf)) {
    return false;
  }
  *out = clamp_probability(1.0 - cdf);
  return true;
}

[[gnu::always_inline]] inline double exgauss_raw_pdf(double x, double mu,
                                                double sigma,
                                                double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  double direct = 0.0;
  if (exgauss_raw_pdf_direct(x, mu, sigma, tau, &direct)) {
    return direct;
  }
  return safe_density(std::exp(exgauss_raw_log_pdf_stable(x, mu, sigma, tau)));
}

[[gnu::always_inline]] inline double exgauss_raw_cdf(double x, double mu,
                                                double sigma,
                                                double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  double direct = 0.0;
  if (exgauss_raw_cdf_direct(x, mu, sigma, tau, &direct)) {
    return direct;
  }
  return clamp_probability(
      std::exp(exgauss_raw_log_cdf_stable(x, mu, sigma, tau)));
}

[[gnu::always_inline]] inline double exgauss_raw_survival(double x, double mu,
                                                     double sigma,
                                                     double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 1.0;
  }
  double direct = 0.0;
  if (exgauss_raw_survival_direct(x, mu, sigma, tau, &direct)) {
    return direct;
  }
  const double log_cdf = exgauss_raw_log_cdf_stable(x, mu, sigma, tau);
  if (log_cdf == R_NegInf) {
    return 1.0;
  }
  return clamp_probability(std::exp(log1m_exp(log_cdf)));
}

constexpr double kLogPi = 1.1447298858494001741434;

inline double lba_denom(double v, double sv) noexcept {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = R::pnorm(v / sv, 0.0, 1.0, 1, 0);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

inline double lba_pdf_fast(double x, double v, double B, double A,
                           double sv) noexcept {
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
    pdf = (v * (R::pnorm(cz, 0.0, 1.0, 1, 0) -
                R::pnorm(cz_max, 0.0, 1.0, 1, 0)) +
           sv * (R::dnorm(cz_max, 0.0, 1.0, 0) -
                 R::dnorm(cz, 0.0, 1.0, 0))) /
          (A * denom);
  } else {
    pdf = R::dnorm(B / x, v, sv, 0) * B / (x * x * denom);
  }
  return safe_density(pdf);
}

inline double lba_cdf_fast(double x, double v, double B, double A,
                           double sv) noexcept {
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
    cdf = (1.0 + (zs * (R::dnorm(cz_max, 0.0, 1.0, 0) -
                        R::dnorm(cz, 0.0, 1.0, 0)) +
                  xx * R::pnorm(cz_max, 0.0, 1.0, 1, 0) -
                  cmz * R::pnorm(cz, 0.0, 1.0, 1, 0)) /
                     A) /
          denom;
  } else {
    cdf = R::pnorm(B / x, v, sv, 0, 0) / denom;
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
    return clamp_probability(2.0 * (1.0 - R::pnorm(z, 0.0, 1.0, 1, 0)));
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double p1 = 1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 + x / mu),
                                   0.0, 1.0, 1, 0);
  const double p2 = 1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 - x / mu),
                                   0.0, 1.0, 1, 0);
  const double part =
      std::exp(std::exp(std::log(2.0 * lambda) - std::log(mu)) +
               std::log(std::max(1e-300, p1)));
  return clamp_probability(part + p2);
}

inline double rdm_digt0(double x, double k, double l) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  const double lambda = k * k;
  double exponent = 0.0;
  if (l == 0.0) {
    exponent = -0.5 * lambda / x;
  } else {
    const double mu = k / l;
    exponent = -(lambda / (2.0 * x)) *
               ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
  }
  return std::exp(exponent + 0.5 * std::log(lambda) -
                  0.5 * std::log(2.0 * x * x * x * M_PI));
}

inline double rdm_pigt(double x, double k, double l, double a,
                       double threshold = 1e-10) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return rdm_pigt0(x, k, l);
  }
  const double sqt = std::sqrt(x);
  const double lgt = std::log(x);
  double cdf = 0.0;
  if (l < threshold) {
    const double t5a = 2.0 * R::pnorm((k + a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
    const double t5b = 2.0 * R::pnorm((-k - a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
    const double t6a =
        -0.5 * ((k + a) * (k + a) / x - M_LN2 - kLogPi + lgt) - std::log(a);
    const double t6b =
        -0.5 * ((k - a) * (k - a) / x - M_LN2 - kLogPi + lgt) - std::log(a);
    cdf = 1.0 + std::exp(t6a) - std::exp(t6b) +
          ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
  } else {
    const double t1a = std::exp(-0.5 * std::pow(k - a - x * l, 2.0) / x);
    const double t1b = std::exp(-0.5 * std::pow(a + k - x * l, 2.0) / x);
    const double t1 =
        std::exp(0.5 * (lgt - M_LN2 - kLogPi)) * (t1a - t1b);
    const double t2a =
        std::exp(2.0 * l * (k - a) +
                 R::pnorm(-(k - a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2b =
        std::exp(2.0 * l * (k + a) +
                 R::pnorm(-(k + a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2 = a + (t2b - t2a) / (2.0 * l);
    const double t4a =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t4b =
        2.0 * R::pnorm((k - a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t4 = 0.5 * (x * l - a - k + 0.5 / l) * t4a +
                      0.5 * (k - a - x * l - 0.5 / l) * t4b;
    cdf = 0.5 * (t4 + t2 + t1) / a;
  }
  if (!std::isfinite(cdf) || cdf < 0.0) {
    return 0.0;
  }
  return clamp_probability(cdf);
}

inline double rdm_digt(double x, double k, double l, double a,
                       double threshold = 1e-10) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return safe_density(rdm_digt0(x, k, l));
  }
  double pdf = 0.0;
  if (l < threshold) {
    const double term = std::exp(-(k - a) * (k - a) / (2.0 * x)) -
                        std::exp(-(k + a) * (k + a) / (2.0 * x));
    pdf = std::exp(-0.5 * (M_LN2 + kLogPi + std::log(x)) +
                   std::log(std::max(1e-300, term)) - M_LN2 - std::log(a));
  } else {
    const double sqt = std::sqrt(x);
    const double t1a = -std::pow(a - k + x * l, 2.0) / (2.0 * x);
    const double t1b = -std::pow(a + k - x * l, 2.0) / (2.0 * x);
    const double t1 =
        0.7071067811865475244 *
        (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
    const double t2a =
        2.0 * R::pnorm((-k + a) / sqt + sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2b =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2 = std::exp(std::log(0.5) + std::log(l)) * (t2a + t2b);
    pdf = std::exp(std::log(std::max(1e-300, t1 + t2)) - M_LN2 - std::log(a));
  }
  return safe_density(pdf);
}

inline double rdm_pdf_fast(double x, double v, double B, double A,
                           double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  if (!std::isfinite(v_sc) || v_sc < 0.0) {
    return 0.0;
  }
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_digt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double rdm_cdf_fast(double x, double v, double B, double A,
                           double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  if (!std::isfinite(v_sc) || v_sc < 0.0) {
    return 0.0;
  }
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

constexpr std::uint8_t kLeafChannelPdf = 1U;
constexpr std::uint8_t kLeafChannelCdf = 2U;
constexpr std::uint8_t kLeafChannelSurvival = 4U;
constexpr std::uint8_t kLeafChannelAll =
    kLeafChannelPdf | kLeafChannelCdf | kLeafChannelSurvival;

struct LeafChannelFill {
  std::uint8_t mask{0U};
  double pdf{0.0};
  double cdf{0.0};
  double survival{1.0};
};

inline LeafChannelFill leaf_channel_zero_fill(
    const std::uint8_t mask) noexcept {
  LeafChannelFill fill;
  fill.mask = mask;
  return fill;
}

[[gnu::always_inline]] inline LeafChannelFill leaf_channel_from_pdf_cdf(
    const double pdf,
    const double cdf,
    const std::uint8_t mask) noexcept {
  LeafChannelFill fill;
  fill.mask = mask;
  if ((mask & kLeafChannelPdf) != 0U) {
    fill.pdf = safe_density(pdf);
  }
  if ((mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    const double clamped = clamp_probability(cdf);
    if ((mask & kLeafChannelCdf) != 0U) {
      fill.cdf = clamped;
    }
    if ((mask & kLeafChannelSurvival) != 0U) {
      fill.survival = clamp_probability(1.0 - clamped);
    }
  }
  return fill;
}

[[gnu::always_inline]] inline LeafChannelFill lognormal_leaf_kernel_fill(
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  const bool need_pdf = (mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return leaf_channel_from_pdf_cdf(
      need_pdf ? R::dlnorm(x, params[0], params[1], 0) : 0.0,
      need_cdf ? R::plnorm(x, params[0], params[1], 1, 0) : 0.0,
      mask);
}

[[gnu::always_inline]] inline LeafChannelFill gamma_leaf_kernel_fill(
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  const bool need_pdf = (mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double shape = params[0];
  const double scale = 1.0 / params[1];
  return leaf_channel_from_pdf_cdf(
      need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
      need_cdf ? R::pgamma(x, shape, scale, 1, 0) : 0.0,
      mask);
}

[[gnu::always_inline]] inline LeafChannelFill exgauss_leaf_kernel_fill(
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  const bool need_pdf = (mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (mask & kLeafChannelCdf) != 0U;
  const bool need_survival = (mask & kLeafChannelSurvival) != 0U;
  const double mu = params[0];
  const double sigma = params[1];
  const double tau = params[2];
  const double lower_survival = exgauss_raw_survival(0.0, mu, sigma, tau);
  if (!(lower_survival > 0.0)) {
    return leaf_channel_zero_fill(mask);
  }

  LeafChannelFill fill;
  fill.mask = mask;
  if (need_pdf) {
    fill.pdf = safe_density(
        exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival);
  }
  if (need_cdf && need_survival) {
    double raw_cdf = 0.0;
    if (exgauss_raw_cdf_direct(x, mu, sigma, tau, &raw_cdf)) {
      const double cdf = clamp_probability(
          (raw_cdf + lower_survival - 1.0) / lower_survival);
      fill.cdf = cdf;
      fill.survival = clamp_probability(1.0 - cdf);
    } else if (x >= mu) {
      const double survival = clamp_probability(
          exgauss_raw_survival(x, mu, sigma, tau) / lower_survival);
      fill.survival = survival;
      fill.cdf = clamp_probability(1.0 - survival);
    } else {
      const double cdf = clamp_probability(
          (exgauss_raw_cdf(x, mu, sigma, tau) + lower_survival - 1.0) /
          lower_survival);
      fill.cdf = cdf;
      fill.survival = clamp_probability(1.0 - cdf);
    }
    return fill;
  }
  if (need_cdf) {
    const double raw_cdf = exgauss_raw_cdf(x, mu, sigma, tau);
    fill.cdf = clamp_probability(
        (raw_cdf + lower_survival - 1.0) / lower_survival);
  }
  if (need_survival) {
    fill.survival = clamp_probability(
        exgauss_raw_survival(x, mu, sigma, tau) / lower_survival);
  }
  return fill;
}

[[gnu::always_inline]] inline LeafChannelFill lba_leaf_kernel_fill(
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  const bool need_pdf = (mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return leaf_channel_from_pdf_cdf(
      need_pdf ? lba_pdf_fast(x, params[0], params[1], params[2], params[3])
               : 0.0,
      need_cdf ? lba_cdf_fast(x, params[0], params[1], params[2], params[3])
               : 0.0,
      mask);
}

[[gnu::always_inline]] inline LeafChannelFill rdm_leaf_kernel_fill(
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  const bool need_pdf = (mask & kLeafChannelPdf) != 0U;
  const bool need_cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return leaf_channel_from_pdf_cdf(
      need_pdf ? rdm_pdf_fast(x, params[0], params[1], params[2], params[3])
               : 0.0,
      need_cdf ? rdm_cdf_fast(x, params[0], params[1], params[2], params[3])
               : 0.0,
      mask);
}

[[gnu::always_inline]] inline LeafChannelFill leaf_kernel_fill(
    const std::uint8_t leaf_dist_kind,
    const double *params,
    const double x,
    const std::uint8_t mask) noexcept {
  if (!(x > 0.0)) {
    return leaf_channel_zero_fill(mask);
  }
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return lognormal_leaf_kernel_fill(params, x, mask);
  case leaf::DistKind::Gamma:
    return gamma_leaf_kernel_fill(params, x, mask);
  case leaf::DistKind::Exgauss:
    return exgauss_leaf_kernel_fill(params, x, mask);
  case leaf::DistKind::LBA:
    return lba_leaf_kernel_fill(params, x, mask);
  case leaf::DistKind::RDM:
    return rdm_leaf_kernel_fill(params, x, mask);
  }
  return LeafChannelFill{};
}

[[gnu::always_inline]] inline double leaf_kernel_value(
    const std::uint8_t leaf_dist_kind,
    const double *params,
    const double x,
    const std::uint8_t channel_mask) noexcept {
  if (!(x > 0.0)) {
    return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
  }
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(R::dlnorm(x, params[0], params[1], 0));
    }
    if (channel_mask == kLeafChannelCdf) {
      return clamp_probability(R::plnorm(x, params[0], params[1], 1, 0));
    }
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(
                     1.0 - R::plnorm(x, params[0], params[1], 1, 0))
               : 0.0;
  case leaf::DistKind::Gamma: {
    const double shape = params[0];
    const double scale = 1.0 / params[1];
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(R::dgamma(x, shape, scale, 0));
    }
    if (channel_mask == kLeafChannelCdf) {
      return clamp_probability(R::pgamma(x, shape, scale, 1, 0));
    }
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(1.0 - R::pgamma(x, shape, scale, 1, 0))
               : 0.0;
  }
  case leaf::DistKind::Exgauss: {
    const double mu = params[0];
    const double sigma = params[1];
    const double tau = params[2];
    const double lower_survival = exgauss_raw_survival(0.0, mu, sigma, tau);
    if (!(lower_survival > 0.0)) {
      return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
    }
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(
          exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival);
    }
    if (channel_mask == kLeafChannelCdf) {
      return clamp_probability(
          (exgauss_raw_cdf(x, mu, sigma, tau) + lower_survival - 1.0) /
          lower_survival);
    }
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(
                     exgauss_raw_survival(x, mu, sigma, tau) / lower_survival)
               : 0.0;
  }
  case leaf::DistKind::LBA:
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(
          lba_pdf_fast(x, params[0], params[1], params[2], params[3]));
    }
    if (channel_mask == kLeafChannelCdf) {
      return lba_cdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(
                     1.0 -
                     lba_cdf_fast(x, params[0], params[1], params[2],
                                  params[3]))
               : 0.0;
  case leaf::DistKind::RDM:
    if (channel_mask == kLeafChannelPdf) {
      return safe_density(
          rdm_pdf_fast(x, params[0], params[1], params[2], params[3]));
    }
    if (channel_mask == kLeafChannelCdf) {
      return rdm_cdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    return channel_mask == kLeafChannelSurvival
               ? clamp_probability(
                     1.0 -
                     rdm_cdf_fast(x, params[0], params[1], params[2],
                                  params[3]))
               : 0.0;
  }
  return 0.0;
}

template <typename Fn>
double integrate_to_infinity(Fn &&density_fn) {
  return quadrature::integrate_tail_default(
      [&](const double t) {
        const double density = density_fn(t);
        return std::isfinite(density) && density > 0.0 ? density : 0.0;
      });
}

} // namespace detail
} // namespace accumulatr::eval
