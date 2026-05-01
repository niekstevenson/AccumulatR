#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "../leaf/channels.hpp"
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

inline leaf::EventChannels impossible_channels() noexcept {
  return leaf::EventChannels::impossible();
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

inline double normal_pdf_fast(double z) noexcept {
  if (!std::isfinite(z)) {
    return 0.0;
  }
  return kInvSqrtTwoPi * std::exp(-0.5 * z * z);
}

inline double normal_pdf_fast(double x, double mean, double sd) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mean) || !std::isfinite(sd) ||
      sd <= 0.0) {
    return 0.0;
  }
  return normal_pdf_fast((x - mean) / sd) / sd;
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

inline double normal_cdf_fast(double z) noexcept {
  if (!std::isfinite(z)) {
    return z > 0.0 ? 1.0 : 0.0;
  }
  const double arg = -z * kInvSqrtTwo;
  return clamp_probability(0.5 * std::erfc(arg));
}

inline double normal_upper_cdf_fast(double z) noexcept {
  if (!std::isfinite(z)) {
    return z > 0.0 ? 0.0 : 1.0;
  }
  return clamp_probability(0.5 * std::erfc(z * kInvSqrtTwo));
}

inline double normal_log_cdf_fast(double z) noexcept {
  if (!std::isfinite(z)) {
    return z > 0.0 ? 0.0 : -std::numeric_limits<double>::infinity();
  }
  if (z > -10.0) {
    return std::log(std::max(1e-300, normal_cdf_fast(z)));
  }
  const double x = -z;
  const double inv_x2 = 1.0 / (x * x);
  const double series =
      std::max(1e-12, 1.0 - inv_x2 + 3.0 * inv_x2 * inv_x2);
  return -0.5 * z * z - std::log(x) - kLogSqrtTwoPi + std::log(series);
}

inline double lognormal_pdf_fast(double x, double m, double s) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(m) ||
      !std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double log_x = std::log(x);
  const double z = (log_x - m) / s;
  return safe_density(normal_pdf_fast(z) / (x * s));
}

inline double lognormal_cdf_fast(double x, double m, double s) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(m) ||
      !std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  return normal_cdf_fast((std::log(x) - m) / s);
}

inline double gamma_pdf_fast_rate(double x, double shape,
                                  double rate) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(shape) ||
      shape <= 0.0 || !std::isfinite(rate) || rate <= 0.0) {
    return 0.0;
  }
  return safe_density(
      std::exp(shape * std::log(rate) + (shape - 1.0) * std::log(x) -
               rate * x - std::lgamma(shape)));
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
  return safe_density(inv_tau * std::exp(exponent) * tail);
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

constexpr double kLogPi = 1.1447298858494001741434;

inline double lba_denom(double v, double sv) noexcept {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = normal_cdf_fast(v / sv);
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
    pdf = (v * (normal_cdf_fast(cz) - normal_cdf_fast(cz_max)) +
           sv * (normal_pdf_fast(cz_max) - normal_pdf_fast(cz))) /
          (A * denom);
  } else {
    pdf = normal_pdf_fast(B / x, v, sv) * B / (x * x * denom);
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
    cdf = (1.0 + (zs * (normal_pdf_fast(cz_max) -
                        normal_pdf_fast(cz)) +
                  xx * normal_cdf_fast(cz_max) -
                  cmz * normal_cdf_fast(cz)) /
                     A) /
          denom;
  } else {
    cdf = normal_upper_cdf_fast((B / x - v) / sv) / denom;
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
    return clamp_probability(2.0 * normal_upper_cdf_fast(z));
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double p1 =
      normal_upper_cdf_fast(std::sqrt(lambda / x) * (1.0 + x / mu));
  const double p2 =
      normal_upper_cdf_fast(std::sqrt(lambda / x) * (1.0 - x / mu));
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
    const double t5a = 2.0 * normal_cdf_fast((k + a) / sqt) - 1.0;
    const double t5b = 2.0 * normal_cdf_fast((-k - a) / sqt) - 1.0;
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
                 normal_log_cdf_fast(-(k - a + x * l) / sqt));
    const double t2b =
        std::exp(2.0 * l * (k + a) +
                 normal_log_cdf_fast(-(k + a + x * l) / sqt));
    const double t2 = a + (t2b - t2a) / (2.0 * l);
    const double t4a =
        2.0 * normal_cdf_fast((k + a) / sqt - sqt * l) - 1.0;
    const double t4b =
        2.0 * normal_cdf_fast((k - a) / sqt - sqt * l) - 1.0;
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
        kInvSqrtTwo *
        (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
    const double t2a =
        2.0 * normal_cdf_fast((-k + a) / sqt + sqt * l) - 1.0;
    const double t2b =
        2.0 * normal_cdf_fast((k + a) / sqt - sqt * l) - 1.0;
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
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

constexpr std::uint8_t kLeafChannelPdf = 1U;
constexpr std::uint8_t kLeafChannelCdf = 2U;
constexpr std::uint8_t kLeafChannelSurvival = 4U;
constexpr std::uint8_t kLeafChannelAll =
    kLeafChannelPdf | kLeafChannelCdf | kLeafChannelSurvival;

inline leaf::EventChannels standard_leaf_channels_mask(
    const std::uint8_t dist_kind,
    const double *params,
    const int n_params,
    const double q,
    const double t0,
    const double t,
    const std::uint8_t channel_mask) {
  (void)n_params;
  if (!std::isfinite(t) || !std::isfinite(q) || q < 0.0 || q > 1.0 ||
      !std::isfinite(t0)) {
    return impossible_channels();
  }

  const double x = t - t0;
  if (!(x > 0.0)) {
    return impossible_channels();
  }

  if (channel_mask == kLeafChannelAll) {
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    double lower_cdf = 0.0;

    switch (static_cast<leaf::DistKind>(dist_kind)) {
    case leaf::DistKind::Lognormal: {
      const double m = params[0];
      const double s = params[1];
      if (!std::isfinite(m) || !std::isfinite(s) || s <= 0.0) {
        return impossible_channels();
      }
      base_pdf = lognormal_pdf_fast(x, m, s);
      base_cdf = lognormal_cdf_fast(x, m, s);
      break;
    }
    case leaf::DistKind::Gamma: {
      const double shape = params[0];
      const double rate = params[1];
      if (!std::isfinite(shape) || shape <= 0.0 || !std::isfinite(rate) ||
          rate <= 0.0) {
        return impossible_channels();
      }
      const double scale = 1.0 / rate;
      base_pdf = gamma_pdf_fast_rate(x, shape, rate);
      base_cdf = R::pgamma(x, shape, scale, 1, 0);
      break;
    }
    case leaf::DistKind::Exgauss: {
      const double mu = params[0];
      const double sigma = params[1];
      const double tau = params[2];
      if (!std::isfinite(mu) || !std::isfinite(sigma) || sigma <= 0.0 ||
          !std::isfinite(tau) || tau <= 0.0) {
        return impossible_channels();
      }
      base_pdf = exgauss_raw_pdf(x, mu, sigma, tau);
      base_cdf = exgauss_raw_cdf(x, mu, sigma, tau);
      lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
      const double lower_survival = 1.0 - lower_cdf;
      if (!(lower_survival > 0.0)) {
        return impossible_channels();
      }
      base_pdf = safe_density(base_pdf / lower_survival);
      base_cdf = clamp_probability((base_cdf - lower_cdf) / lower_survival);
      break;
    }
    case leaf::DistKind::LBA: {
      base_pdf = lba_pdf_fast(x, params[0], params[1], params[2], params[3]);
      base_cdf = lba_cdf_fast(x, params[0], params[1], params[2], params[3]);
      break;
    }
    case leaf::DistKind::RDM: {
      base_pdf = rdm_pdf_fast(x, params[0], params[1], params[2], params[3]);
      base_cdf = rdm_cdf_fast(x, params[0], params[1], params[2], params[3]);
      break;
    }
    }

    base_pdf = std::isfinite(base_pdf) && base_pdf > 0.0 ? base_pdf : 0.0;
    base_cdf = clamp_probability(base_cdf);

    leaf::EventChannels out;
    const double start_prob = 1.0 - q;
    out.pdf = start_prob * base_pdf;
    out.cdf = clamp_probability(start_prob * base_cdf);
    out.survival = clamp_probability(1.0 - out.cdf);
    return out;
  }

  double base_pdf = 0.0;
  double base_cdf = 0.0;
  double lower_cdf = 0.0;
  const bool need_pdf = (channel_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (channel_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;

  switch (static_cast<leaf::DistKind>(dist_kind)) {
  case leaf::DistKind::Lognormal: {
    const double m = params[0];
    const double s = params[1];
    if (!std::isfinite(m) || !std::isfinite(s) || s <= 0.0) {
      return impossible_channels();
    }
    if (need_pdf) {
      base_pdf = lognormal_pdf_fast(x, m, s);
    }
    if (need_cdf) {
      base_cdf = lognormal_cdf_fast(x, m, s);
    }
    break;
  }
  case leaf::DistKind::Gamma: {
    const double shape = params[0];
    const double rate = params[1];
    if (!std::isfinite(shape) || shape <= 0.0 || !std::isfinite(rate) ||
        rate <= 0.0) {
      return impossible_channels();
    }
    const double scale = 1.0 / rate;
    if (need_pdf) {
      base_pdf = gamma_pdf_fast_rate(x, shape, rate);
    }
    if (need_cdf) {
      base_cdf = R::pgamma(x, shape, scale, 1, 0);
    }
    break;
  }
  case leaf::DistKind::Exgauss: {
    const double mu = params[0];
    const double sigma = params[1];
    const double tau = params[2];
    if (!std::isfinite(mu) || !std::isfinite(sigma) || sigma <= 0.0 ||
        !std::isfinite(tau) || tau <= 0.0) {
      return impossible_channels();
    }
    if (need_pdf) {
      base_pdf = exgauss_raw_pdf(x, mu, sigma, tau);
    }
    if (need_cdf) {
      base_cdf = exgauss_raw_cdf(x, mu, sigma, tau);
    }
    lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      return impossible_channels();
    }
    if (need_pdf) {
      base_pdf = safe_density(base_pdf / lower_survival);
    }
    if (need_cdf) {
      base_cdf = clamp_probability((base_cdf - lower_cdf) / lower_survival);
    }
    break;
  }
  case leaf::DistKind::LBA: {
    if (need_pdf) {
      base_pdf = lba_pdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    if (need_cdf) {
      base_cdf = lba_cdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    break;
  }
  case leaf::DistKind::RDM: {
    if (need_pdf) {
      base_pdf = rdm_pdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    if (need_cdf) {
      base_cdf = rdm_cdf_fast(x, params[0], params[1], params[2], params[3]);
    }
    break;
  }
  }

  leaf::EventChannels out;
  const double start_prob = 1.0 - q;
  if (need_pdf) {
    base_pdf = std::isfinite(base_pdf) && base_pdf > 0.0 ? base_pdf : 0.0;
    out.pdf = start_prob * base_pdf;
  }
  if (need_cdf) {
    base_cdf = clamp_probability(base_cdf);
    const double cdf = clamp_probability(start_prob * base_cdf);
    if ((channel_mask & kLeafChannelCdf) != 0U) {
      out.cdf = cdf;
    }
    if ((channel_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(1.0 - cdf);
    }
  }
  return out;
}

inline leaf::EventChannels standard_leaf_channels(const std::uint8_t dist_kind,
                                                  const double *params,
                                                  const int n_params,
                                                  const double q,
                                                  const double t0,
                                                  const double t) {
  return standard_leaf_channels_mask(
      dist_kind,
      params,
      n_params,
      q,
      t0,
      t,
      kLeafChannelAll);
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
