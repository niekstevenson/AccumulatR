#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <R_ext/Error.h>

#include "exact_source_channels.hpp"

namespace accumulatr::eval {
namespace detail {

struct BatchActiveLaneSpan {
  const semantic::Index *lanes{nullptr};
  std::size_t size{0};
};

struct BatchTimeSlotView {
  const double *values{nullptr};
  std::size_t lane_stride{0};
  const std::uint8_t *valid{nullptr};

  bool has(const semantic::Index time_id,
           const semantic::Index lane) const noexcept {
    if (valid == nullptr) {
      return true;
    }
    return valid[
        static_cast<std::size_t>(time_id) * lane_stride +
        static_cast<std::size_t>(lane)] != 0U;
  }

  double get(const semantic::Index time_id,
             const semantic::Index lane) const noexcept {
    return values[
        static_cast<std::size_t>(time_id) * lane_stride +
        static_cast<std::size_t>(lane)];
  }
};

inline double batch_source_product_channel_time(
    const CompiledMathSourceProductChannel &channel,
    const BatchTimeSlotView &time_slots,
    const semantic::Index lane) noexcept {
  double time = time_slots.get(channel.time_id, lane);
  if (channel.time_cap_id != semantic::kInvalidIndex) {
    time = std::min(time, time_slots.get(channel.time_cap_id, lane));
  }
  return time;
}

inline double batch_source_product_finish_base_value(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t channel_mask) noexcept {
  const double start_prob = 1.0 - q;
  if (channel_mask == kLeafChannelPdf) {
    return start_prob * safe_density(base_pdf);
  }
  const double cdf = clamp_probability(start_prob * clamp_probability(base_cdf));
  if (channel_mask == kLeafChannelCdf) {
    return cdf;
  }
  return channel_mask == kLeafChannelSurvival
             ? clamp_probability(1.0 - cdf)
             : 0.0;
}

inline double batch_source_product_before_onset_value(
    const std::uint8_t channel_mask) noexcept {
  return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
}

inline bool batch_source_product_ops_are_direct_leaf_only(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops) noexcept {
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    const auto channel_pos =
        static_cast<std::size_t>(op.source_product_channel_id);
    if (channel_pos >= program.integral_kernel_source_product_channels.size()) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[channel_pos];
    if (!channel.direct_leaf_absolute_candidate ||
        channel.leaf_index == semantic::kInvalidIndex) {
      return false;
    }
  }
  return true;
}

inline bool batch_source_product_ops_have_available_direct_leaf_inputs(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    ExactSourceChannels *source_channels) noexcept {
  if (source_channels == nullptr) {
    return false;
  }
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    const auto channel_pos =
        static_cast<std::size_t>(op.source_product_channel_id);
    if (channel_pos >= program.integral_kernel_source_product_channels.size()) {
      return false;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[channel_pos];
    if (!channel.direct_leaf_absolute_candidate ||
        channel.leaf_index == semantic::kInvalidIndex ||
        !source_channels->source_product_direct_leaf_available(
            op.source_product_channel_id)) {
      return false;
    }
  }
  return true;
}

inline semantic::Index batch_source_product_max_leaf_index(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops) noexcept {
  semantic::Index max_leaf = semantic::kInvalidIndex;
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    if (op.value_channel_mask == 0U) {
      continue;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    if (channel.leaf_index != semantic::kInvalidIndex &&
        (max_leaf == semantic::kInvalidIndex ||
         channel.leaf_index > max_leaf)) {
      max_leaf = channel.leaf_index;
    }
  }
  return max_leaf;
}

inline bool batch_finite_integral_kernel_direct_leaf_supported(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel) noexcept {
  return kernel.kind == CompiledMathIntegralKernelKind::SourceProduct &&
         kernel.bind_time_id != semantic::kInvalidIndex &&
         batch_source_product_ops_are_direct_leaf_only(
             program, kernel.source_product_ops);
}

} // namespace detail
} // namespace accumulatr::eval
