#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>

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

inline double normal_cdf_fast(double z) noexcept {
  const double arg = -z * 0.7071067811865475244;
  return clamp_probability(0.5 * std::erfc(arg));
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

inline leaf::EventChannels standard_leaf_channels(const std::uint8_t dist_kind,
                                                  const double *params,
                                                  const int n_params,
                                                  const double q,
                                                  const double t0,
                                                  const double t) {
  (void)n_params;
  if (!std::isfinite(t) || !std::isfinite(q) || q < 0.0 || q > 1.0 ||
      !std::isfinite(t0)) {
    return impossible_channels();
  }

  const double x = t - t0;
  if (!(x > 0.0)) {
    return impossible_channels();
  }

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
    base_pdf = R::dlnorm(x, m, s, 0);
    base_cdf = R::plnorm(x, m, s, 1, 0);
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
    base_pdf = R::dgamma(x, shape, scale, 0);
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
