#pragma once

#include <cstring>
#include "exact_lane_source_state.hpp"
#include "exact_adaptive.hpp"

namespace accumulatr::eval::detail {

inline constexpr std::size_t kCumulativeMinimumRequests = 4U;

struct CumulativeLaneWorkspace {
  std::vector<std::uint64_t> keys, hashes;
  std::vector<std::size_t> order, previous;
  std::vector<double> errors, values, lower, refine_upper, tolerance_divisors;
};

// Exact parameter equality only; no sharing of numerical results across particles.
inline bool prepare_cumulative_lanes(
    const ExactVariantPlan &plan, const CompiledMathIntegralKernel &kernel,
    const ExactLaneSourceState &sources, const CompiledLaneFrame &frame,
    const semantic::Index *lanes, const std::size_t count,
    double *lower, const double *upper, CumulativeLaneWorkspace *work) {
  work->order.clear();
  for (std::size_t i = 0; i < count; ++i)
    if (upper[i] > lower[i] && std::isfinite(upper[i])) work->order.push_back(i);
  if (work->order.size() < kCumulativeMinimumRequests) return false;
  std::size_t width = 0;
  for (const auto leaf : kernel.cumulative_leaves)
    width += 3 + leaf::dist_param_count(
        static_cast<leaf::DistKind>(plan.leaf_descriptors[leaf].dist_kind));
  work->keys.resize(width * count);
  work->hashes.assign(count, 14695981039346656037ULL);
  work->previous.assign(count, count);
  const auto set = [&](std::size_t i, std::size_t field, double value) {
    std::uint64_t bits = 0;
    if (value != 0) std::memcpy(&bits, &value, sizeof(bits));
    work->keys[i * width + field] = bits;
    work->hashes[i] = (work->hashes[i] ^ bits) * 1099511628211ULL;
  };
  std::size_t field = 0;
  for (const auto leaf : kernel.cumulative_leaves) {
    const auto input = sources.leaf_batch(leaf);
    const auto parameters = leaf::dist_param_count(
        static_cast<leaf::DistKind>(plan.leaf_descriptors[leaf].dist_kind));
    for (const auto i : work->order) {
      const auto lane = static_cast<std::size_t>(compiled_frame_lane(lanes, i));
      const auto source = frame.source_lanes_identity ? lane : frame.source_lanes[lane];
      const auto row = input.physical_row(source);
      set(i, field, input.q(source, row));
      set(i, field + 1, input.t0(row));
      set(i, field + 2, input.onset(row));
      for (int p = 0; p < parameters; ++p) set(i, field + 3 + p, input.param(p, row));
    }
    field += 3 + parameters;
  }
  const auto compare_keys = [&](std::size_t a, std::size_t b) {
    return std::memcmp(work->keys.data() + a * width,
                       work->keys.data() + b * width, width * sizeof(std::uint64_t));
  };
  std::sort(work->order.begin(), work->order.end(), [&](auto a, auto b) {
    if (work->hashes[a] != work->hashes[b]) return work->hashes[a] < work->hashes[b];
    const auto comparison = compare_keys(a, b);
    return comparison != 0 ? comparison < 0 : upper[a] < upper[b];
  });
  bool shared = false;
  for (std::size_t begin = 0; begin < work->order.size();) {
    std::size_t end = begin + 1;
    while (end < work->order.size() &&
           work->hashes[work->order[begin]] == work->hashes[work->order[end]] &&
           compare_keys(work->order[begin], work->order[end]) == 0) ++end;
    if (end - begin >= kCumulativeMinimumRequests) {
      if (!shared) {
        work->lower.assign(lower, lower + count);
        work->tolerance_divisors.assign(count, 1.0);
        shared = true;
      }
      std::size_t intervals = 1;
      for (std::size_t i = begin + 1; i < end; ++i)
        intervals += upper[work->order[i]] > upper[work->order[i - 1]];
      for (std::size_t i = begin; i < end; ++i)
        work->tolerance_divisors[work->order[i]] = intervals;
      for (std::size_t i = begin + 1; i < end; ++i) {
        const auto current = work->order[i], previous = work->order[i-1];
        lower[current] = std::max(lower[current], upper[previous]);
        work->previous[current] = previous;
      }
    }
    begin = end;
  }
  return shared;
}

template <typename Evaluate>
inline void finish_cumulative_lanes(
    const std::size_t count, const double *upper,
    Evaluate &evaluate, AdaptiveLaneWorkspace *adaptive,
    CumulativeLaneWorkspace *work, std::vector<double> *out,
    const double *absolute_scales) {
  work->errors.resize(count);
  for (std::size_t i = 0; i < count; ++i) work->errors[i] = adaptive->lanes[i].error;
  work->refine_upper.clear();
  for (const auto i : work->order) {
    const auto previous = work->previous[i];
    if (previous != count) {
      (*out)[i] += (*out)[previous];
      work->errors[i] += work->errors[previous];
    }
    // Signed integrands may cancel: enforce the error budget on every prefix.
    if (!adaptive_integral_converged((*out)[i], work->errors[i],
            kAdaptiveAbsoluteTolerance / (absolute_scales == nullptr ? 1.0 : absolute_scales[i]),
            kAdaptiveRelativeTolerance)) {
      if (work->refine_upper.empty()) work->refine_upper.assign(count, 0.0);
      work->refine_upper[i] = upper[i];
    }
  }
  if (!work->refine_upper.empty()) {
    adaptive_integrate_lane_batch(count, work->lower.data(), work->refine_upper.data(),
        evaluate, adaptive, &work->values,
        kAdaptiveAbsoluteTolerance, kAdaptiveRelativeTolerance, absolute_scales);
    for (std::size_t i = 0; i < count; ++i) {
      if (work->refine_upper[i] > 0) {
        (*out)[i] = work->values[i];
        work->errors[i] = adaptive->lanes[i].error;
      }
    }
  }
}

} // namespace accumulatr::eval::detail
