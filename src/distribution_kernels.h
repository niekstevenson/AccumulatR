#pragma once

#include "native_utils.h"

#include <algorithm>
#include <cmath>

namespace uuber::distkernels {

inline constexpr double kLogPi = 1.1447298858494001741434;
inline constexpr double kInvSqrt2 = 0.7071067811865475;
inline constexpr double kInvSqrt2Pi = 0.3989422804014327;

inline double normal_cdf(double z) {
  const double arg = -z * kInvSqrt2;
  const double val = 0.5 * std::erfc(arg);
  return clamp(val, 0.0, 1.0);
}

inline double lognormal_pdf(double x, double meanlog, double sdlog) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) ||
      !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  const double norm = kInvSqrt2Pi * inv_sigma;
  const double val = (norm / x) * std::exp(-0.5 * z * z);
  return std::isfinite(val) && val > 0.0 ? val : 0.0;
}

inline double lognormal_cdf(double x, double meanlog, double sdlog) {
  if (!std::isfinite(meanlog) || !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  if (!std::isfinite(x)) {
    if (std::isnan(x)) {
      return 0.0;
    }
    return x < 0.0 ? 0.0 : 1.0;
  }
  if (x <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  return clamp_probability(normal_cdf(z));
}

inline double gamma_pdf(double x, double shape, double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    return 0.0;
  }
  if (!std::isfinite(x)) {
    return 0.0;
  }
  const double scale = 1.0 / rate;
  const double val = R::dgamma(x, shape, scale, 0);
  return (!std::isfinite(val) || val < 0.0) ? 0.0 : val;
}

inline double gamma_cdf(double x, double shape, double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    return 0.0;
  }
  if (!std::isfinite(x)) {
    if (std::isnan(x)) {
      return 0.0;
    }
    return x < 0.0 ? 0.0 : 1.0;
  }
  const double scale = 1.0 / rate;
  const double val = R::pgamma(x, shape, scale, 1, 0);
  return std::isfinite(val) ? clamp(val, 0.0, 1.0) : 0.0;
}

inline double exgauss_pdf(double x, double mu, double sigma, double tau) {
  if (!std::isfinite(sigma) || sigma <= 0.0 || !std::isfinite(tau) ||
      tau <= 0.0) {
    return 0.0;
  }
  if (!std::isfinite(x)) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = half_sigma_sq_over_tau_sq - (x - mu) * inv_tau;
  const double gaussian_tail = normal_cdf(z - sigma_over_tau);
  const double val = inv_tau * std::exp(exponent) * gaussian_tail;
  return (!std::isfinite(val) || val < 0.0) ? 0.0 : val;
}

inline double exgauss_cdf(double x, double mu, double sigma, double tau) {
  if (!std::isfinite(sigma) || sigma <= 0.0 || !std::isfinite(tau) ||
      tau <= 0.0) {
    return 0.0;
  }
  if (!std::isfinite(x)) {
    if (std::isnan(x)) {
      return 0.0;
    }
    return x < 0.0 ? 0.0 : 1.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double base_cdf = normal_cdf(z);
  const double tail = normal_cdf(z - sigma_over_tau);
  const double exp_term =
      std::exp(half_sigma_sq_over_tau_sq - (x - mu) * inv_tau);
  const double val = base_cdf - exp_term * tail;
  return std::isfinite(val) ? clamp(val, 0.0, 1.0) : 0.0;
}

inline double lba_denom(double v, double sv) {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = R::pnorm(v / sv, 0.0, 1.0, 1, 0);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

inline double lba_pdf(double x, double v, double sv, double B, double A) {
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

inline double lba_cdf(double x, double v, double sv, double B, double A) {
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
    cdf = (1.0 +
           (zs * (R::dnorm(cz_max, 0.0, 1.0, 0) -
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

inline double rdm_pigt0(double x, double k, double l) {
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
  const double p1 =
      1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 + x / mu), 0.0, 1.0, 1, 0);
  const double p2 =
      1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 - x / mu), 0.0, 1.0, 1, 0);
  const double part = std::exp(std::exp(std::log(2.0 * lambda) - std::log(mu)) +
                               std::log(std::max(1e-300, p1)));
  return clamp_probability(part + p2);
}

inline double rdm_digt0(double x, double k, double l) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  const double lambda = k * k;
  double e = 0.0;
  if (l == 0.0) {
    e = -0.5 * lambda / x;
  } else {
    const double mu = k / l;
    e = -(lambda / (2.0 * x)) * ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
  }
  return std::exp(e + 0.5 * std::log(lambda) -
                  0.5 * std::log(2.0 * x * x * x * M_PI));
}

inline double rdm_pigt(double x, double k, double l, double a,
                       double threshold = 1e-10) {
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
    const double t1 = std::exp(0.5 * (lgt - M_LN2 - kLogPi)) * (t1a - t1b);
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
                       double threshold = 1e-10) {
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
        M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
    const double t2a =
        2.0 * R::pnorm((-k + a) / sqt + sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2b =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2 = std::exp(std::log(0.5) + std::log(l)) * (t2a + t2b);
    pdf = std::exp(std::log(std::max(1e-300, t1 + t2)) - M_LN2 - std::log(a));
  }
  return safe_density(pdf);
}

inline double rdm_pdf(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_digt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double rdm_cdf(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

} // namespace uuber::distkernels
