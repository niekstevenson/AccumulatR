#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <R_ext/Error.h>
#include <R_ext/Print.h>

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

  double get(const semantic::Index time_id,
             const semantic::Index lane) const noexcept {
    return values[
        static_cast<std::size_t>(time_id) * lane_stride +
        static_cast<std::size_t>(lane)];
  }
};

struct BatchLeafInputView {
  const ExactLoadedLeafInput *inputs{nullptr};
  std::size_t lane_stride{0};

  const ExactLoadedLeafInput &get(
      const semantic::Index leaf_index,
      const semantic::Index lane) const noexcept {
    return inputs[
        static_cast<std::size_t>(leaf_index) * lane_stride +
        static_cast<std::size_t>(lane)];
  }
};

struct BatchSourceProductContext {
  std::size_t lane_count{0};
  BatchTimeSlotView time_slots{};
  BatchLeafInputView leaf_inputs{};
  semantic::Index cache_scope_id{0};
};

struct BatchSourceProductWorkspace {
  std::vector<double> source_value_by_lane;
  std::vector<semantic::Index> active_a;
  std::vector<semantic::Index> active_b;

  void ensure_size(const std::size_t lane_count,
                   const std::size_t active_count) {
    if (source_value_by_lane.size() < lane_count) {
      source_value_by_lane.resize(lane_count, 0.0);
    }
    if (active_a.size() < active_count) {
      active_a.resize(active_count);
    }
    if (active_b.size() < active_count) {
      active_b.resize(active_count);
    }
  }
};

struct BatchSourceProductComparison {
  bool supported{false};
  bool matched{false};
  std::size_t compared_lanes{0};
  double max_abs_error{0.0};
  semantic::Index max_error_lane{semantic::kInvalidIndex};
  double scalar_value{0.0};
  double batch_value{0.0};
};

struct BatchSourceProductDebugCounters {
  std::uint64_t checked_spans{0};
  std::uint64_t supported_spans{0};
  std::uint64_t unsupported_spans{0};
  std::uint64_t compared_lanes{0};
  std::uint64_t mismatched_lanes{0};
  double max_abs_error{0.0};
  double scalar_value_at_max{0.0};
  double batch_value_at_max{0.0};
};

struct BatchFiniteIntegralDebugCounters {
  std::uint64_t checked_kernels{0};
  std::uint64_t supported_kernels{0};
  std::uint64_t unsupported_kernels{0};
  std::uint64_t compared_lanes{0};
  std::uint64_t mismatched_lanes{0};
  std::uint64_t child_lanes{0};
  double max_abs_error{0.0};
  double scalar_value_at_max{0.0};
  double batch_value_at_max{0.0};
};

struct BatchFiniteIntegralWorkspace {
  std::vector<double> child_time_values;
  std::vector<ExactLoadedLeafInput> child_leaf_inputs;
  std::vector<semantic::Index> child_active;
  std::vector<semantic::Index> child_parent_lanes;
  std::vector<double> child_weights;
  std::vector<double> child_values;
  BatchSourceProductWorkspace source_product_workspace;

  void ensure_size(const std::size_t time_slot_count,
                   const std::size_t leaf_count,
                   const std::size_t child_capacity) {
    const auto time_value_count = time_slot_count * child_capacity;
    if (child_time_values.size() < time_value_count) {
      child_time_values.resize(time_value_count, 0.0);
    }
    const auto leaf_value_count = leaf_count * child_capacity;
    if (child_leaf_inputs.size() < leaf_value_count) {
      child_leaf_inputs.resize(leaf_value_count);
    }
    if (child_active.size() < child_capacity) {
      child_active.resize(child_capacity);
    }
    if (child_parent_lanes.size() < child_capacity) {
      child_parent_lanes.resize(child_capacity);
    }
    if (child_weights.size() < child_capacity) {
      child_weights.resize(child_capacity, 0.0);
    }
    if (child_values.size() < child_capacity) {
      child_values.resize(child_capacity, 0.0);
    }
  }
};

inline bool batch_source_product_env_enabled(const char *name) {
  const char *value = std::getenv(name);
  return value != nullptr && value[0] != '\0' &&
         std::strcmp(value, "0") != 0 &&
         std::strcmp(value, "false") != 0 &&
         std::strcmp(value, "FALSE") != 0;
}

inline bool batch_source_product_debug_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled("ACCUMULATR_BATCH_SP_COMPARE");
  return enabled;
}

inline bool batch_source_product_debug_report_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled("ACCUMULATR_BATCH_SP_COMPARE_REPORT");
  return enabled;
}

inline bool batch_source_product_debug_strict_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled("ACCUMULATR_BATCH_SP_COMPARE_STRICT");
  return enabled;
}

inline bool batch_finite_integral_debug_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled("ACCUMULATR_BATCH_INTEGRAL_COMPARE");
  return enabled;
}

inline bool batch_finite_integral_debug_report_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled(
          "ACCUMULATR_BATCH_INTEGRAL_COMPARE_REPORT");
  return enabled;
}

inline bool batch_finite_integral_debug_strict_enabled() {
  static const bool enabled =
      batch_source_product_env_enabled(
          "ACCUMULATR_BATCH_INTEGRAL_COMPARE_STRICT");
  return enabled;
}

inline BatchSourceProductDebugCounters &
batch_source_product_debug_counters() {
  static BatchSourceProductDebugCounters counters;
  return counters;
}

inline BatchFiniteIntegralDebugCounters &
batch_finite_integral_debug_counters() {
  static BatchFiniteIntegralDebugCounters counters;
  return counters;
}

inline void batch_source_product_debug_reset() {
  batch_source_product_debug_counters() = BatchSourceProductDebugCounters{};
}

inline void batch_finite_integral_debug_reset() {
  batch_finite_integral_debug_counters() =
      BatchFiniteIntegralDebugCounters{};
}

inline void batch_source_product_debug_report() {
  if (!batch_source_product_debug_enabled() ||
      !batch_source_product_debug_report_enabled()) {
    return;
  }
  const auto &counters = batch_source_product_debug_counters();
  Rprintf(
      "AccumulatR batch source-product compare: checked=%llu, "
      "supported=%llu, unsupported=%llu, compared_lanes=%llu, "
      "mismatched_lanes=%llu, max_abs_error=%.17g, scalar_at_max=%.17g, "
      "batch_at_max=%.17g\n",
      static_cast<unsigned long long>(counters.checked_spans),
      static_cast<unsigned long long>(counters.supported_spans),
      static_cast<unsigned long long>(counters.unsupported_spans),
      static_cast<unsigned long long>(counters.compared_lanes),
      static_cast<unsigned long long>(counters.mismatched_lanes),
      counters.max_abs_error,
      counters.scalar_value_at_max,
      counters.batch_value_at_max);
}

inline void batch_finite_integral_debug_report() {
  if (!batch_finite_integral_debug_enabled() ||
      !batch_finite_integral_debug_report_enabled()) {
    return;
  }
  const auto &counters = batch_finite_integral_debug_counters();
  Rprintf(
      "AccumulatR batch finite-integral compare: checked=%llu, "
      "supported=%llu, unsupported=%llu, compared_lanes=%llu, "
      "child_lanes=%llu, mismatched_lanes=%llu, max_abs_error=%.17g, "
      "scalar_at_max=%.17g, batch_at_max=%.17g\n",
      static_cast<unsigned long long>(counters.checked_kernels),
      static_cast<unsigned long long>(counters.supported_kernels),
      static_cast<unsigned long long>(counters.unsupported_kernels),
      static_cast<unsigned long long>(counters.compared_lanes),
      static_cast<unsigned long long>(counters.child_lanes),
      static_cast<unsigned long long>(counters.mismatched_lanes),
      counters.max_abs_error,
      counters.scalar_value_at_max,
      counters.batch_value_at_max);
}

class BatchSourceProductDebugScope {
public:
  BatchSourceProductDebugScope()
      : enabled_(batch_source_product_debug_enabled()) {
    if (enabled_) {
      batch_source_product_debug_reset();
    }
  }

  BatchSourceProductDebugScope(const BatchSourceProductDebugScope &) = delete;
  BatchSourceProductDebugScope &operator=(
      const BatchSourceProductDebugScope &) = delete;

  ~BatchSourceProductDebugScope() {
    if (enabled_) {
      batch_source_product_debug_report();
    }
  }

private:
  bool enabled_{false};
};

class BatchFiniteIntegralDebugScope {
public:
  BatchFiniteIntegralDebugScope()
      : enabled_(batch_finite_integral_debug_enabled()) {
    if (batch_finite_integral_debug_enabled()) {
      batch_finite_integral_debug_reset();
    }
  }

  BatchFiniteIntegralDebugScope(const BatchFiniteIntegralDebugScope &) = delete;
  BatchFiniteIntegralDebugScope &operator=(
      const BatchFiniteIntegralDebugScope &) = delete;

  ~BatchFiniteIntegralDebugScope() {
    if (!enabled_) {
      return;
    }
    if (batch_finite_integral_debug_enabled()) {
      batch_finite_integral_debug_report();
    }
  }

private:
  bool enabled_{false};
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

inline double batch_source_product_lognormal_leaf_value(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t channel_mask) {
  if (!(x > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  const double m = loaded.params[0];
  const double s = loaded.params[1];
  return batch_source_product_finish_base_value(
      need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
      need_pdf ? 0.0 : R::plnorm(x, m, s, 1, 0),
      loaded.q,
      channel_mask);
}

inline double batch_source_product_gamma_leaf_value(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t channel_mask) {
  if (!(x > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  const double shape = loaded.params[0];
  const double scale = 1.0 / loaded.params[1];
  return batch_source_product_finish_base_value(
      need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
      need_pdf ? 0.0 : R::pgamma(x, shape, scale, 1, 0),
      loaded.q,
      channel_mask);
}

inline double batch_source_product_exgauss_leaf_value(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t channel_mask) {
  if (!(x > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  const double mu = loaded.params[0];
  const double sigma = loaded.params[1];
  const double tau = loaded.params[2];
  const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
  const double lower_survival = 1.0 - lower_cdf;
  if (!(lower_survival > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  return batch_source_product_finish_base_value(
      need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
      need_pdf
          ? 0.0
          : (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                lower_survival,
      loaded.q,
      channel_mask);
}

inline double batch_source_product_lba_leaf_value(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t channel_mask) {
  if (!(x > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  return batch_source_product_finish_base_value(
      need_pdf
          ? lba_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_pdf
          ? 0.0
          : lba_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3]),
      loaded.q,
      channel_mask);
}

inline double batch_source_product_rdm_leaf_value(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t channel_mask) {
  if (!(x > 0.0)) {
    return batch_source_product_before_onset_value(channel_mask);
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  return batch_source_product_finish_base_value(
      need_pdf
          ? rdm_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_pdf
          ? 0.0
          : rdm_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3]),
      loaded.q,
      channel_mask);
}

inline double batch_source_product_direct_leaf_value_for_lane(
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t channel_mask) {
  const double x = time - channel.leaf_onset_abs_value - loaded.t0;
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return batch_source_product_lognormal_leaf_value(
        loaded, x, channel_mask);
  case leaf::DistKind::Gamma:
    return batch_source_product_gamma_leaf_value(loaded, x, channel_mask);
  case leaf::DistKind::Exgauss:
    return batch_source_product_exgauss_leaf_value(
        loaded, x, channel_mask);
  case leaf::DistKind::LBA:
    return batch_source_product_lba_leaf_value(loaded, x, channel_mask);
  case leaf::DistKind::RDM:
    return batch_source_product_rdm_leaf_value(loaded, x, channel_mask);
  }
  return 0.0;
}

template <typename LeafValueFn>
inline void batch_source_product_direct_leaf_values_for_kind(
    const CompiledMathSourceProductChannel &channel,
    const BatchSourceProductContext &context,
    const BatchActiveLaneSpan active,
    const std::uint8_t channel_mask,
    double *values_by_lane,
    LeafValueFn leaf_value_fn) {
  const auto leaf_index = channel.leaf_index;
  for (std::size_t i = 0; i < active.size; ++i) {
    const auto lane = active.lanes[i];
    const auto &loaded = context.leaf_inputs.get(leaf_index, lane);
    const double x =
        batch_source_product_channel_time(channel, context.time_slots, lane) -
        channel.leaf_onset_abs_value - loaded.t0;
    values_by_lane[static_cast<std::size_t>(lane)] =
        leaf_value_fn(loaded, x, channel_mask);
  }
}

inline void batch_source_product_direct_leaf_values(
    const CompiledMathSourceProductChannel &channel,
    const BatchSourceProductContext &context,
    const BatchActiveLaneSpan active,
    const std::uint8_t channel_mask,
    double *values_by_lane) {
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    batch_source_product_direct_leaf_values_for_kind(
        channel,
        context,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_lognormal_leaf_value);
    break;
  case leaf::DistKind::Gamma:
    batch_source_product_direct_leaf_values_for_kind(
        channel,
        context,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_gamma_leaf_value);
    break;
  case leaf::DistKind::Exgauss:
    batch_source_product_direct_leaf_values_for_kind(
        channel,
        context,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_exgauss_leaf_value);
    break;
  case leaf::DistKind::LBA:
    batch_source_product_direct_leaf_values_for_kind(
        channel,
        context,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_lba_leaf_value);
    break;
  case leaf::DistKind::RDM:
    batch_source_product_direct_leaf_values_for_kind(
        channel,
        context,
        active,
        channel_mask,
        values_by_lane,
        batch_source_product_rdm_leaf_value);
    break;
  }
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

inline bool batch_source_product_value_for_ops_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const BatchSourceProductContext &context,
    const BatchActiveLaneSpan active,
    BatchSourceProductWorkspace *workspace,
    double *out_by_lane) {
  if (workspace == nullptr || out_by_lane == nullptr ||
      !batch_source_product_ops_are_direct_leaf_only(
          program, source_product_ops)) {
    return false;
  }

  workspace->ensure_size(context.lane_count, active.size);
  if (active.size > 0U) {
    std::copy(
        active.lanes,
        active.lanes + active.size,
        workspace->active_a.begin());
  }
  std::size_t current_count = active.size;
  auto *current_active = &workspace->active_a;
  auto *next_active = &workspace->active_b;

  for (std::size_t i = 0; i < current_count; ++i) {
    out_by_lane[static_cast<std::size_t>((*current_active)[i])] = 1.0;
  }

  for (semantic::Index op_idx = 0; op_idx < source_product_ops.size; ++op_idx) {
    if (current_count == 0U) {
      break;
    }
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + op_idx)];
    const auto channel_mask = op.value_channel_mask;
    if (channel_mask == 0U) {
      if (op.constant_value == 0.0) {
        for (std::size_t i = 0; i < current_count; ++i) {
          out_by_lane[static_cast<std::size_t>((*current_active)[i])] = 0.0;
        }
        current_count = 0U;
        break;
      }
      for (std::size_t i = 0; i < current_count; ++i) {
        out_by_lane[static_cast<std::size_t>((*current_active)[i])] *=
            op.constant_value;
      }
      continue;
    }

    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    const BatchActiveLaneSpan current_span{
        current_active->data(),
        current_count};
    batch_source_product_direct_leaf_values(
        channel,
        context,
        current_span,
        channel_mask,
        workspace->source_value_by_lane.data());

    std::size_t next_count = 0U;
    for (std::size_t i = 0; i < current_count; ++i) {
      const auto lane = (*current_active)[i];
      const auto lane_pos = static_cast<std::size_t>(lane);
      const double product =
          out_by_lane[lane_pos] * workspace->source_value_by_lane[lane_pos];
      out_by_lane[lane_pos] =
          std::isfinite(product) && product > 0.0 ? product : 0.0;
      if (out_by_lane[lane_pos] > 0.0) {
        (*next_active)[next_count++] = lane;
      }
    }
    std::swap(current_active, next_active);
    current_count = next_count;
  }

  return true;
}

inline double batch_source_product_scalar_reference_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const BatchSourceProductContext &context,
    const semantic::Index lane) {
  double product = 1.0;
  for (semantic::Index op_idx = 0; op_idx < source_product_ops.size; ++op_idx) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + op_idx)];
    if (op.value_channel_mask == 0U) {
      if (op.constant_value == 0.0) {
        return 0.0;
      }
      product *= op.constant_value;
      continue;
    }
    const auto &channel =
        program.integral_kernel_source_product_channels[
            static_cast<std::size_t>(op.source_product_channel_id)];
    product *= batch_source_product_direct_leaf_value_for_lane(
        channel,
        context.leaf_inputs.get(channel.leaf_index, lane),
        batch_source_product_channel_time(channel, context.time_slots, lane),
        op.value_channel_mask);
    if (!std::isfinite(product) || product == 0.0) {
      return 0.0;
    }
  }
  return product;
}

inline BatchSourceProductComparison
batch_compare_source_product_direct_leaf_to_scalar(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    const BatchSourceProductContext &context,
    const BatchActiveLaneSpan active,
    BatchSourceProductWorkspace *workspace,
    double *batch_out_by_lane,
    const double tolerance = 1e-12) {
  BatchSourceProductComparison comparison;
  comparison.supported = batch_source_product_value_for_ops_direct_leaf(
      program,
      source_product_ops,
      context,
      active,
      workspace,
      batch_out_by_lane);
  if (!comparison.supported) {
    return comparison;
  }

  comparison.matched = true;
  comparison.compared_lanes = active.size;
  for (std::size_t i = 0; i < active.size; ++i) {
    const auto lane = active.lanes[i];
    const double scalar = batch_source_product_scalar_reference_direct_leaf(
        program, source_product_ops, context, lane);
    const double batch = batch_out_by_lane[static_cast<std::size_t>(lane)];
    const double error = std::fabs(scalar - batch);
    if (error > comparison.max_abs_error) {
      comparison.max_abs_error = error;
      comparison.max_error_lane = lane;
      comparison.scalar_value = scalar;
      comparison.batch_value = batch;
    }
    if (!(error <= tolerance)) {
      comparison.matched = false;
    }
  }
  return comparison;
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

inline semantic::Index batch_source_product_max_time_id_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops) noexcept {
  semantic::Index max_time = semantic::kInvalidIndex;
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
    if (channel.time_id != semantic::kInvalidIndex &&
        (max_time == semantic::kInvalidIndex ||
         channel.time_id > max_time)) {
      max_time = channel.time_id;
    }
    if (channel.time_cap_id != semantic::kInvalidIndex &&
        (max_time == semantic::kInvalidIndex ||
         channel.time_cap_id > max_time)) {
      max_time = channel.time_cap_id;
    }
  }
  return max_time;
}

inline bool batch_finite_integral_kernel_direct_leaf_supported(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel) noexcept {
  return kernel.kind == CompiledMathIntegralKernelKind::SourceProduct &&
         kernel.bind_time_id != semantic::kInvalidIndex &&
         batch_source_product_ops_are_direct_leaf_only(
             program, kernel.source_product_ops);
}

inline bool batch_finite_integral_source_product_direct_leaf(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    const BatchSourceProductContext &parent_context,
    const BatchActiveLaneSpan parent_active,
    const double *lower_by_lane,
    const double *upper_by_lane,
    const std::size_t time_slot_count,
    const std::size_t leaf_count,
    BatchFiniteIntegralWorkspace *workspace,
    double *out_by_parent_lane,
    std::size_t *child_lane_count_out = nullptr) {
  if (child_lane_count_out != nullptr) {
    *child_lane_count_out = 0U;
  }
  if (workspace == nullptr || out_by_parent_lane == nullptr ||
      lower_by_lane == nullptr || upper_by_lane == nullptr ||
      !batch_finite_integral_kernel_direct_leaf_supported(program, kernel)) {
    return false;
  }
  const auto max_time =
      batch_source_product_max_time_id_for_ops(program, kernel.source_product_ops);
  const auto required_time_count =
      static_cast<std::size_t>(
          std::max(max_time, kernel.bind_time_id)) +
      1U;
  if (time_slot_count < required_time_count) {
    return false;
  }

  const std::size_t child_capacity =
      parent_active.size * quadrature::kDefaultFiniteOrder;
  workspace->ensure_size(time_slot_count, leaf_count, child_capacity);
  for (std::size_t i = 0; i < parent_active.size; ++i) {
    out_by_parent_lane[
        static_cast<std::size_t>(parent_active.lanes[i])] = 0.0;
  }
  if (child_capacity == 0U) {
    return true;
  }

  const auto &rule =
      quadrature::gauss_legendre_rule<quadrature::kDefaultFiniteOrder>();
  std::size_t child_count = 0U;
  for (std::size_t parent_pos = 0; parent_pos < parent_active.size;
       ++parent_pos) {
    const auto parent_lane = parent_active.lanes[parent_pos];
    const auto parent_lane_pos = static_cast<std::size_t>(parent_lane);
    const double lower = lower_by_lane[parent_lane_pos];
    const double upper = upper_by_lane[parent_lane_pos];
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      continue;
    }
    const double scale = 0.5 * (upper - lower);
    const double shift = 0.5 * (upper + lower);
    for (std::size_t q = 0; q < quadrature::kDefaultFiniteOrder; ++q) {
      const auto child_lane = static_cast<semantic::Index>(child_count);
      workspace->child_active[child_count] = child_lane;
      workspace->child_parent_lanes[child_count] = parent_lane;
      workspace->child_weights[child_count] = scale * rule.weights[q];
      for (std::size_t time_slot = 0; time_slot < time_slot_count;
           ++time_slot) {
        workspace->child_time_values[
            time_slot * child_capacity + child_count] =
            parent_context.time_slots.get(
                static_cast<semantic::Index>(time_slot),
                parent_lane);
      }
      workspace->child_time_values[
          static_cast<std::size_t>(kernel.bind_time_id) * child_capacity +
          child_count] = shift + scale * rule.nodes[q];
      for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
        workspace->child_leaf_inputs[leaf * child_capacity + child_count] =
            parent_context.leaf_inputs.get(
                static_cast<semantic::Index>(leaf),
                parent_lane);
      }
      ++child_count;
    }
  }
  if (child_lane_count_out != nullptr) {
    *child_lane_count_out = child_count;
  }
  if (child_count == 0U) {
    return true;
  }

  const BatchSourceProductContext child_context{
      child_capacity,
      BatchTimeSlotView{workspace->child_time_values.data(), child_capacity},
      BatchLeafInputView{workspace->child_leaf_inputs.data(), child_capacity},
      static_cast<semantic::Index>(parent_context.cache_scope_id + 1)};
  const BatchActiveLaneSpan child_active{
      workspace->child_active.data(),
      child_count};
  if (!batch_source_product_value_for_ops_direct_leaf(
          program,
          kernel.source_product_ops,
          child_context,
          child_active,
          &workspace->source_product_workspace,
          workspace->child_values.data())) {
    return false;
  }
  for (std::size_t child_pos = 0; child_pos < child_count; ++child_pos) {
    const double value = workspace->child_values[child_pos];
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    out_by_parent_lane[
        static_cast<std::size_t>(
            workspace->child_parent_lanes[child_pos])] +=
        workspace->child_weights[child_pos] * value;
  }
  return true;
}

inline void batch_debug_compare_source_product_ops_to_scalar(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double scalar_value,
    const double tolerance = 1e-12) {
  if (!batch_source_product_debug_enabled()) {
    return;
  }

  auto &counters = batch_source_product_debug_counters();
  ++counters.checked_spans;

  if (workspace == nullptr) {
    ++counters.unsupported_spans;
    return;
  }

  if (!batch_source_product_ops_have_available_direct_leaf_inputs(
          program,
          source_product_ops,
          source_channels)) {
    ++counters.unsupported_spans;
    return;
  }

  const auto max_leaf =
      batch_source_product_max_leaf_index(program, source_product_ops);
  if (max_leaf == semantic::kInvalidIndex) {
    ++counters.unsupported_spans;
    return;
  }

  std::vector<ExactLoadedLeafInput> leaf_inputs(
      static_cast<std::size_t>(max_leaf) + 1U);
  for (semantic::Index leaf = 0; leaf <= max_leaf; ++leaf) {
    leaf_inputs[static_cast<std::size_t>(leaf)] =
        source_channels->source_product_leaf_input(leaf);
  }

  std::vector<double> time_values = workspace->time_values;
  const semantic::Index lane = 0;
  double batch_value = 0.0;
  BatchSourceProductWorkspace batch_workspace;
  const BatchSourceProductContext context{
      1U,
      BatchTimeSlotView{time_values.data(), 1U},
      BatchLeafInputView{leaf_inputs.data(), 1U},
      static_cast<semantic::Index>(
          workspace->source_product_program_current_epoch)};
  const BatchActiveLaneSpan active{&lane, 1U};

  const bool supported = batch_source_product_value_for_ops_direct_leaf(
      program,
      source_product_ops,
      context,
      active,
      &batch_workspace,
      &batch_value);
  if (!supported) {
    ++counters.unsupported_spans;
    return;
  }

  ++counters.supported_spans;
  ++counters.compared_lanes;
  const double error = std::fabs(scalar_value - batch_value);
  if (error > counters.max_abs_error) {
    counters.max_abs_error = error;
    counters.scalar_value_at_max = scalar_value;
    counters.batch_value_at_max = batch_value;
  }
  if (!(error <= tolerance)) {
    ++counters.mismatched_lanes;
    if (batch_source_product_debug_strict_enabled()) {
      Rf_error(
          "batch source-product comparison failed: scalar=%.17g, "
          "batch=%.17g, abs_error=%.17g, tolerance=%.17g",
          scalar_value,
          batch_value,
          error,
          tolerance);
    }
  }
}

inline void batch_debug_compare_finite_integral_to_scalar(
    const CompiledMathProgram &program,
    const CompiledMathIntegralKernel &kernel,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double lower,
    const double upper,
    const double scalar_value,
    const double tolerance = 1e-12) {
  if (!batch_finite_integral_debug_enabled()) {
    return;
  }

  auto &counters = batch_finite_integral_debug_counters();
  ++counters.checked_kernels;

  if (workspace == nullptr || source_channels == nullptr ||
      !batch_finite_integral_kernel_direct_leaf_supported(program, kernel) ||
      !batch_source_product_ops_have_available_direct_leaf_inputs(
          program,
          kernel.source_product_ops,
          source_channels)) {
    ++counters.unsupported_kernels;
    return;
  }

  const auto max_leaf =
      batch_source_product_max_leaf_index(program, kernel.source_product_ops);
  const std::size_t leaf_count =
      max_leaf == semantic::kInvalidIndex
          ? 0U
          : static_cast<std::size_t>(max_leaf) + 1U;
  std::vector<ExactLoadedLeafInput> leaf_inputs(leaf_count);
  for (semantic::Index leaf = 0; leaf <= max_leaf; ++leaf) {
    leaf_inputs[static_cast<std::size_t>(leaf)] =
        source_channels->source_product_leaf_input(leaf);
  }

  std::vector<double> parent_time_values = workspace->time_values;
  const auto max_time =
      batch_source_product_max_time_id_for_ops(program, kernel.source_product_ops);
  const auto required_time_count =
      static_cast<std::size_t>(
          std::max(max_time, kernel.bind_time_id)) +
      1U;
  if (parent_time_values.size() < required_time_count) {
    parent_time_values.resize(required_time_count, 0.0);
  }

  const semantic::Index parent_lane = 0;
  const double lower_by_lane[1] = {lower};
  const double upper_by_lane[1] = {upper};
  double batch_value_by_lane[1] = {0.0};
  BatchFiniteIntegralWorkspace batch_workspace;
  const BatchSourceProductContext parent_context{
      1U,
      BatchTimeSlotView{parent_time_values.data(), 1U},
      BatchLeafInputView{leaf_inputs.data(), 1U},
      static_cast<semantic::Index>(
          workspace->source_product_program_current_epoch)};
  const BatchActiveLaneSpan parent_active{&parent_lane, 1U};
  std::size_t child_lanes = 0U;
  if (!batch_finite_integral_source_product_direct_leaf(
          program,
          kernel,
          parent_context,
          parent_active,
          lower_by_lane,
          upper_by_lane,
          parent_time_values.size(),
          leaf_count,
          &batch_workspace,
          batch_value_by_lane,
          &child_lanes)) {
    ++counters.unsupported_kernels;
    return;
  }

  ++counters.supported_kernels;
  ++counters.compared_lanes;
  counters.child_lanes += static_cast<std::uint64_t>(child_lanes);
  const double batch_value = batch_value_by_lane[0];
  const double error = std::fabs(scalar_value - batch_value);
  if (error > counters.max_abs_error) {
    counters.max_abs_error = error;
    counters.scalar_value_at_max = scalar_value;
    counters.batch_value_at_max = batch_value;
  }
  if (!(error <= tolerance)) {
    ++counters.mismatched_lanes;
    if (batch_finite_integral_debug_strict_enabled()) {
      Rf_error(
          "batch finite-integral comparison failed: scalar=%.17g, "
          "batch=%.17g, abs_error=%.17g, tolerance=%.17g",
          scalar_value,
          batch_value,
          error,
          tolerance);
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
