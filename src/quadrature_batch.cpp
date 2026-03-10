#include "quadrature_batch.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace uuber {
namespace {

constexpr int kGauss15N = 15;
constexpr double kGauss15X[kGauss15N] = {
    -0.9879925180204854, -0.9372733924007060, -0.8482065834104272,
    -0.7244177313601700, -0.5709721726085388, -0.3941513470775634,
    -0.2011940939974345, 0.0, 0.2011940939974345, 0.3941513470775634,
    0.5709721726085388,  0.7244177313601700,  0.8482065834104272,
    0.9372733924007060,  0.9879925180204854};
constexpr double kGauss15W[kGauss15N] = {
    0.0307532419961173, 0.0703660474881081, 0.1071592204671719,
    0.1395706779261543, 0.1662692058169939, 0.1861610000155622,
    0.1984314853271116, 0.2025782419255613, 0.1984314853271116,
    0.1861610000155622, 0.1662692058169939, 0.1395706779261543,
    0.1071592204671719, 0.0703660474881081, 0.0307532419961173};
constexpr int kFiniteSegments = 16;
constexpr int kInfiniteUniformSegments = 24;
constexpr int kInfiniteTailSegments = 20;
constexpr double kInfiniteUniformUpper = 0.9;

} // namespace

TimeBatch build_time_batch(double lower, double upper) {
  TimeBatch batch;
  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    return batch;
  }
  if (upper < lower) {
    std::swap(lower, upper);
  }
  if (!(upper > lower)) {
    return batch;
  }
  batch.lower = lower;
  batch.upper = upper;
  batch.nodes.resize(kGauss15N);
  batch.weights.resize(kGauss15N);
  const double a = 0.5 * (upper - lower);
  const double b = 0.5 * (upper + lower);
  for (int i = 0; i < kGauss15N; ++i) {
    batch.nodes[static_cast<std::size_t>(i)] = b + a * kGauss15X[i];
    batch.weights[static_cast<std::size_t>(i)] = a * kGauss15W[i];
  }
  return batch;
}

TimeBatch build_time_batch_0_to_upper(double upper) {
  if (std::isfinite(upper)) {
    TimeBatch batch;
    if (upper <= 0.0) {
      return batch;
    }

    batch.lower = 0.0;
    batch.upper = upper;
    const double width = upper;
    if (!std::isfinite(width) || width <= 0.0) {
      return {};
    }
    batch.nodes.reserve(static_cast<std::size_t>(kFiniteSegments * kGauss15N));
    batch.weights.reserve(static_cast<std::size_t>(kFiniteSegments * kGauss15N));
    for (int seg = 0; seg < kFiniteSegments; ++seg) {
      const double seg_lo = width * static_cast<double>(seg) /
                            static_cast<double>(kFiniteSegments);
      const double seg_hi = width * static_cast<double>(seg + 1) /
                            static_cast<double>(kFiniteSegments);
      TimeBatch seg_batch = build_time_batch(seg_lo, seg_hi);
      if (seg_batch.nodes.empty() ||
          seg_batch.nodes.size() != seg_batch.weights.size()) {
        return {};
      }
      batch.nodes.insert(batch.nodes.end(), seg_batch.nodes.begin(),
                         seg_batch.nodes.end());
      batch.weights.insert(batch.weights.end(), seg_batch.weights.begin(),
                           seg_batch.weights.end());
    }
    return batch;
  }

  TimeBatch batch;
  if (!(upper > 0.0)) {
    return batch;
  }
  batch.lower = 0.0;
  batch.upper = std::numeric_limits<double>::infinity();
  std::vector<double> x_edges;
  x_edges.reserve(static_cast<std::size_t>(kInfiniteUniformSegments +
                                           kInfiniteTailSegments + 3));
  for (int i = 0; i <= kInfiniteUniformSegments; ++i) {
    x_edges.push_back(kInfiniteUniformUpper * static_cast<double>(i) /
                      static_cast<double>(kInfiniteUniformSegments));
  }
  const double tail_width = 1.0 - kInfiniteUniformUpper;
  for (int i = 1; i <= kInfiniteTailSegments; ++i) {
    const double x =
        1.0 - tail_width * std::pow(0.5, static_cast<double>(i));
    if (x > x_edges.back()) {
      x_edges.push_back(x);
    }
  }
  if (x_edges.back() < 1.0) {
    x_edges.push_back(1.0);
  }

  batch.nodes.reserve((x_edges.size() - 1) * kGauss15N);
  batch.weights.reserve((x_edges.size() - 1) * kGauss15N);
  for (std::size_t seg = 0; seg + 1 < x_edges.size(); ++seg) {
    const double x_lo = x_edges[seg];
    const double x_hi = x_edges[seg + 1];
    if (!(x_hi > x_lo) || x_lo < 0.0 || x_hi > 1.0) {
      continue;
    }
    TimeBatch x_batch = build_time_batch(x_lo, x_hi);
    if (x_batch.nodes.empty() || x_batch.nodes.size() != x_batch.weights.size()) {
      continue;
    }
    for (std::size_t i = 0; i < x_batch.nodes.size(); ++i) {
      const double x = x_batch.nodes[i];
      const double wx = x_batch.weights[i];
      if (!std::isfinite(x) || !std::isfinite(wx) || wx <= 0.0 || x < 0.0) {
        continue;
      }
      double one_minus_x = 1.0 - x;
      if (!(one_minus_x > 0.0)) {
        one_minus_x = std::nextafter(1.0, 0.0) - x;
      }
      if (!(one_minus_x > 0.0)) {
        continue;
      }
      const double t = x / one_minus_x;
      const double jac = 1.0 / (one_minus_x * one_minus_x);
      const double wt = wx * jac;
      if (!std::isfinite(t) || !std::isfinite(wt) || wt <= 0.0) {
        continue;
      }
      batch.nodes.push_back(t);
      batch.weights.push_back(wt);
    }
  }
  return batch;
}

} // namespace uuber
