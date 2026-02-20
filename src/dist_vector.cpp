#include "dist_vector.h"

#include <Rmath.h>

#include <cmath>
#include <cstddef>

#include "accumulator.h"

namespace uuber {
namespace {

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

inline double eval_pdf_single(int dist_code, double p1, double p2, double p3,
                              double x) {
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
  default:
    return 0.0;
  }
}

inline double eval_cdf_single(int dist_code, double p1, double p2, double p3,
                              double x) {
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
  default:
    return 0.0;
  }
}

} // namespace

void eval_pdf_vec(int dist_code, double p1, double p2, double p3,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = safe_density(eval_pdf_single(dist_code, p1, p2, p3, x[i]));
  }
}

void eval_cdf_vec(int dist_code, double p1, double p2, double p3,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = clamp_probability(eval_cdf_single(dist_code, p1, p2, p3, x[i]));
  }
}

} // namespace uuber
