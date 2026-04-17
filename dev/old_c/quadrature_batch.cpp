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

TimeBatch build_time_batch_with_observed(double lower, double upper,
                                         const std::vector<double> &observed_times) {
  TimeBatch batch = build_time_batch(lower, upper);
  if (batch.nodes.empty()) {
    return batch;
  }
  for (double t : observed_times) {
    if (!std::isfinite(t) || t < lower || t > upper) {
      continue;
    }
    batch.nodes.push_back(t);
    batch.weights.push_back(0.0);
  }
  std::vector<std::size_t> order(batch.nodes.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return batch.nodes[a] < batch.nodes[b];
  });
  std::vector<double> nodes_sorted;
  std::vector<double> weights_sorted;
  nodes_sorted.reserve(order.size());
  weights_sorted.reserve(order.size());
  for (std::size_t idx : order) {
    nodes_sorted.push_back(batch.nodes[idx]);
    weights_sorted.push_back(batch.weights[idx]);
  }
  batch.nodes.swap(nodes_sorted);
  batch.weights.swap(weights_sorted);
  return batch;
}

double integrate_time_batch(const TimeBatch &batch,
                            const std::vector<double> &values) {
  if (batch.nodes.empty() || batch.weights.empty() ||
      values.size() != batch.nodes.size()) {
    return 0.0;
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    double v = values[i];
    if (!std::isfinite(v) || v <= 0.0) {
      continue;
    }
    double w = batch.weights[i];
    if (!std::isfinite(w) || w <= 0.0) {
      continue;
    }
    sum += w * v;
  }
  if (!std::isfinite(sum)) {
    return 0.0;
  }
  return sum;
}

} // namespace uuber

