#pragma once

#include "src/eval/exact_adaptive.hpp"
#include "src/eval/quadrature.hpp"

// Benchmark-only dispatch. The production executor calls this in a temporary copy.
namespace accumulatr::eval::detail {
inline int benchmark_method = 0;
inline std::size_t benchmark_evaluations = 0;

template <std::size_t N, typename Evaluate>
void benchmark_fixed(std::size_t n, const double *lower, const double *upper,
                     Evaluate &evaluate, AdaptiveLaneWorkspace *w,
                     std::vector<double> *out) {
  const auto &rule = quadrature::gauss_legendre_rule<N>();
  w->mapped_times.resize(kAdaptiveNodeBatchSize);
  w->mapped_lanes.resize(kAdaptiveNodeBatchSize);
  w->values.resize(kAdaptiveNodeBatchSize);
  out->assign(n, 0.0);
  w->active.clear();
  for (std::size_t lane = 0; lane < n; ++lane)
    if (upper[lane] > lower[lane]) w->active.push_back(lane);
  const auto total = w->active.size() * N;
  for (std::size_t begin = 0; begin < total; begin += kAdaptiveNodeBatchSize) {
    const auto count = std::min(kAdaptiveNodeBatchSize, total - begin);
    for (std::size_t j = 0; j < count; ++j) {
      const auto k = begin + j;
      const auto lane = w->active[k / N];
      const double half = (upper[lane] - lower[lane]) / 2;
      w->mapped_lanes[j] = lane;
      w->mapped_times[j] = lower[lane] + half * (1 + rule.nodes[k % N]);
    }
    evaluate(w->mapped_times.data(), w->mapped_lanes.data(), count, w->values.data());
    for (std::size_t j = 0; j < count; ++j) {
      const auto k = begin + j;
      const auto lane = w->mapped_lanes[j];
      (*out)[lane] += w->values[j] * rule.weights[k % N] *
        (upper[lane] - lower[lane]) / 2;
    }
  }
}

template <typename Evaluate>
void benchmark_integrate_lane_batch(std::size_t n, const double *lower,
    const double *upper, Evaluate &&evaluate, AdaptiveLaneWorkspace *w,
    std::vector<double> *out) {
  auto counted = [&](const double *t, const std::size_t *lanes, std::size_t count, double *values) {
    benchmark_evaluations += count;
    evaluate(t, lanes, count, values);
  };
  if (benchmark_method == 31)
    benchmark_fixed<31>(n, lower, upper, counted, w, out);
  else
    adaptive_integrate_lane_batch(n, lower, upper, counted, w, out);
}
}
