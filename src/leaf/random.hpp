#pragma once

#include <Rcpp.h>
#include <cmath>
#include <cstddef>

#include "dist_kind.hpp"

namespace accumulatr::leaf {

inline double positive_normal(const double mean, const double sd) {
  if (mean > 0.0) {
    double value;
    do value = R::rnorm(mean, sd); while (value <= 0.0);
    return value;
  }
  // Exponential rejection in the normal tail avoids inverse-CDF cancellation.
  // Its acceptance probability is at least 0.76, including arbitrarily far tails.
  const double lower = -mean / sd;
  const double rate = 0.5 * lower + 0.5 * std::hypot(lower, 2.0);
  double excess, residual;
  do {
    excess = R::rexp(1.0 / rate);
    residual = excess - 1.0 / rate;
  } while (R::runif(0.0, 1.0) > std::exp(-0.5 * residual * residual));
  return sd * excess;
}

// Parameters are checked once when simulation inputs are prepared. The caller
// owns the RNGScope and adds nondecision time, onset and trigger failures.
inline double sample_time(const DistKind kind, const double *row,
                          const std::ptrdiff_t stride) {
  const double p1 = row[0], p2 = row[stride];
  switch (kind) {
  case DistKind::Lognormal:
    return R::rlnorm(p1, p2);
  case DistKind::Gamma:
    return R::rgamma(p1, 1.0 / p2);
  case DistKind::Exgauss: {
    const double tau = row[2 * stride];
    const double log_positive = R::pnorm(p1 / p2, 0.0, 1.0, true, true);
    const double ratio = p2 / tau;
    const double log_negative = p1 / tau + 0.5 * ratio * ratio +
        R::pnorm(-p1 / p2 - ratio, 0.0, 1.0, true, true);
    const bool positive = R::runif(0.0, 1.0) <
        R::plogis(log_positive - log_negative, 0.0, 1.0, true, false);
    // When N <= 0, conditioning N + Exp > 0 leaves an exponential residual.
    return R::rexp(tau) + (positive ? positive_normal(p1, p2) : 0.0);
  }
  case DistKind::LBA: {
    const double distance = p2 + row[2 * stride] * R::runif(0.0, 1.0);
    return distance / positive_normal(p1, row[3 * stride]);
  }
  case DistKind::RDM: {
    const double scale = row[3 * stride];
    const double distance =
        (p2 + row[2 * stride] * R::runif(0.0, 1.0)) / scale;
    const double drift = p1 / scale;
    const double z = R::rnorm(0.0, 1.0);
    const double half_square = (0.5 * z / distance) * z;
    // Rationalized inverse-Gaussian root; no subtraction of near-equal terms.
    // At zero drift this is exactly the Levy hitting-time distribution.
    const double root = drift + half_square +
        std::sqrt(half_square) * std::sqrt(drift + 0.5 * half_square) *
            std::sqrt(2.0);
    if (drift == 0.0 || R::runif(0.0, 1.0) <= 1.0 / (1.0 + drift / root)) {
      return distance / root;
    }
    return (distance / drift) * (root / drift);
  }
  }
  return NA_REAL;
}

} // namespace accumulatr::leaf
