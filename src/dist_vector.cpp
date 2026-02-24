#include "dist_vector.h"

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "accumulator.h"

namespace uuber {
namespace {

constexpr double kLogPi = 1.1447298858494001741434; // log(pi)

inline double clamp_probability(double x) {
  if (!std::isfinite(x)) {
    return 0.0;
  }
  if (x < 0.0) {
    return 0.0;
  }
  if (x > 1.0) {
    return 1.0;
  }
  return x;
}

inline double safe_density(double x) {
  if (!std::isfinite(x) || x < 0.0) {
    return 0.0;
  }
  return x;
}

inline double lba_denom(double v, double sv) {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = ::pnorm(v / sv, 0.0, 1.0, 1, 0);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

inline double lba_pdf_single(double x, double v, double sv, double B, double A) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) || !std::isfinite(sv) ||
      sv <= 0.0 || !std::isfinite(B) || !std::isfinite(A)) {
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
    pdf = (v * (::pnorm(cz, 0.0, 1.0, 1, 0) - ::pnorm(cz_max, 0.0, 1.0, 1, 0)) +
           sv * (::dnorm(cz_max, 0.0, 1.0, 0) - ::dnorm(cz, 0.0, 1.0, 0))) /
          (A * denom);
  } else {
    pdf = ::dnorm(B / x, v, sv, 0) * B / (x * x * denom);
  }
  return safe_density(pdf);
}

inline double lba_cdf_single(double x, double v, double sv, double B, double A) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) || !std::isfinite(sv) ||
      sv <= 0.0 || !std::isfinite(B) || !std::isfinite(A)) {
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
    cdf = (1.0 + (zs * (::dnorm(cz_max, 0.0, 1.0, 0) - ::dnorm(cz, 0.0, 1.0, 0)) +
                  xx * ::pnorm(cz_max, 0.0, 1.0, 1, 0) -
                  cmz * ::pnorm(cz, 0.0, 1.0, 1, 0)) /
                     A) /
          denom;
  } else {
    cdf = ::pnorm(B / x, v, sv, 0, 0) / denom;
  }
  return clamp_probability(cdf);
}

inline double rdm_pigt0(double x, double k, double l) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l)) {
    return 0.0;
  }
  if (std::fabs(l) < 1e-12) {
    // Levy limit used by the simulator when drift approaches zero.
    const double z = k / std::sqrt(x);
    return clamp_probability(2.0 * (1.0 - ::pnorm(z, 0.0, 1.0, 1, 0)));
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double p1 = 1.0 - ::pnorm(std::sqrt(lambda / x) * (1.0 + x / mu), 0.0,
                                    1.0, 1, 0);
  const double p2 = 1.0 - ::pnorm(std::sqrt(lambda / x) * (1.0 - x / mu), 0.0,
                                    1.0, 1, 0);
  const double part = std::exp(std::exp(std::log(2.0 * lambda) - std::log(mu)) +
                               std::log(std::max(1e-300, p1)));
  return clamp_probability(part + p2);
}

inline double rdm_digt0(double x, double k, double l) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l)) {
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
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l) ||
      !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return rdm_pigt0(x, k, l);
  }
  const double sqt = std::sqrt(x);
  const double lgt = std::log(x);
  double cdf = 0.0;
  if (l < threshold) {
    const double t5a = 2.0 * ::pnorm((k + a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
    const double t5b = 2.0 * ::pnorm((-k - a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
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
    const double t2a = std::exp(2.0 * l * (k - a) +
                                ::pnorm(-(k - a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2b = std::exp(2.0 * l * (k + a) +
                                ::pnorm(-(k + a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2 = a + (t2b - t2a) / (2.0 * l);
    const double t4a =
        2.0 * ::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t4b =
        2.0 * ::pnorm((k - a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
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
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l) ||
      !std::isfinite(a)) {
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
        2.0 * ::pnorm((-k + a) / sqt + sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2b =
        2.0 * ::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2 = std::exp(std::log(0.5) + std::log(l)) * (t2a + t2b);
    pdf = std::exp(std::log(std::max(1e-300, t1 + t2)) - M_LN2 - std::log(a));
  }
  return safe_density(pdf);
}

inline double rdm_pdf_single(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_digt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double rdm_cdf_single(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double eval_pdf_single(int dist_code, double p1, double p2, double p3,
                              double p4, double p5, double p6, double p7,
                              double p8, double x) {
  (void)p5;
  (void)p6;
  (void)p7;
  (void)p8;
  switch (dist_code) {
  case ACC_DIST_LOGNORMAL:
    return ::dlnorm(x, p1, p2, 0);
  case ACC_DIST_GAMMA:
    return ::dgamma(x, p1, 1.0 / p2, 0);
  case ACC_DIST_EXGAUSS: {
    if (!std::isfinite(p2) || p2 <= 0.0 || !std::isfinite(p3) || p3 <= 0.0) {
      return 0.0;
    }
    const double inv_tau = 1.0 / p3;
    const double sigma_sq = p2 * p2;
    const double tau_sq = p3 * p3;
    const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
    const double sigma_over_tau = p2 * inv_tau;
    const double z = (x - p1) / p2;
    const double exponent = half_sigma_sq_over_tau_sq - (x - p1) * inv_tau;
    const double gaussian_tail = ::pnorm(z - sigma_over_tau, 0.0, 1.0, 1, 0);
    return inv_tau * std::exp(exponent) * gaussian_tail;
  }
  case ACC_DIST_LBA:
    return lba_pdf_single(x, p1, p2, p3, p4);
  case ACC_DIST_RDM:
    return rdm_pdf_single(x, p1, p2, p3, p4);
  default:
    return 0.0;
  }
}

inline double eval_cdf_single(int dist_code, double p1, double p2, double p3,
                              double p4, double p5, double p6, double p7,
                              double p8, double x) {
  (void)p5;
  (void)p6;
  (void)p7;
  (void)p8;
  switch (dist_code) {
  case ACC_DIST_LOGNORMAL:
    return ::plnorm(x, p1, p2, 1, 0);
  case ACC_DIST_GAMMA:
    return ::pgamma(x, p1, 1.0 / p2, 1, 0);
  case ACC_DIST_EXGAUSS: {
    if (!std::isfinite(p2) || p2 <= 0.0 || !std::isfinite(p3) || p3 <= 0.0) {
      return 0.0;
    }
    const double inv_tau = 1.0 / p3;
    const double sigma_sq = p2 * p2;
    const double tau_sq = p3 * p3;
    const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
    const double sigma_over_tau = p2 * inv_tau;
    const double z = (x - p1) / p2;
    const double base_cdf = ::pnorm(z, 0.0, 1.0, 1, 0);
    const double tail = ::pnorm(z - sigma_over_tau, 0.0, 1.0, 1, 0);
    const double exp_term =
        std::exp(half_sigma_sq_over_tau_sq - (x - p1) * inv_tau);
    return base_cdf - exp_term * tail;
  }
  case ACC_DIST_LBA:
    return lba_cdf_single(x, p1, p2, p3, p4);
  case ACC_DIST_RDM:
    return rdm_cdf_single(x, p1, p2, p3, p4);
  default:
    return 0.0;
  }
}

} // namespace

void eval_pdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = safe_density(
        eval_pdf_single(dist_code, p1, p2, p3, p4, p5, p6, p7, p8, x[i]));
  }
}

void eval_cdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = clamp_probability(
        eval_cdf_single(dist_code, p1, p2, p3, p4, p5, p6, p7, p8, x[i]));
  }
}

} // namespace uuber
