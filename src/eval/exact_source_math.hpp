#pragma once

#include <cmath>
#include <cstdint>

#include "exact_lane_source_state.hpp"
#include "leaf_kernel.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactSourceFill {
  double pdf{0.0};
  double cdf{0.0};
  double survival{1.0};
};

inline std::uint8_t exact_source_node_fill_mask(
    const CompiledMathNodeKind kind) noexcept {
  const auto value_mask = compiled_math_source_factor_channel_mask(kind);
  return (value_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U
             ? kLeafChannelCdf | kLeafChannelSurvival
             : value_mask;
}

template <std::uint8_t Mask>
[[gnu::always_inline]] inline ExactSourceFill exact_source_finish_base_fill(
    const double base_pdf,
    const double base_cdf,
    const double q) {
  const double start_probability = 1.0 - q;
  ExactSourceFill fill;
  if constexpr ((Mask & kLeafChannelPdf) != 0U) {
    fill.pdf = safe_density(start_probability * base_pdf);
  }
  if constexpr ((Mask &
                 (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    const double cdf = clamp_probability(start_probability * base_cdf);
    if constexpr ((Mask & kLeafChannelCdf) != 0U) {
      fill.cdf = cdf;
    }
    if constexpr ((Mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0 - cdf;
    }
  }
  return fill;
}

inline ExactSourceFill exact_source_impossible_fill() {
  return {};
}

inline ExactSourceFill exact_source_certain_fill(const std::uint8_t mask) {
  ExactSourceFill fill;
  if ((mask & kLeafChannelCdf) != 0U) {
    fill.cdf = 1.0;
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    fill.survival = 0.0;
  }
  return fill;
}

inline ExactSourceFill exact_source_forced_fill(
    const ExactRelation relation,
    const std::uint8_t mask) {
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    return exact_source_certain_fill(mask);
  }
  return exact_source_impossible_fill();
}

inline bool exact_source_relation_forces_fill(
    const ExactRelation relation,
    const std::uint8_t mask) noexcept {
  return relation != ExactRelation::Unknown &&
         !(relation == ExactRelation::At &&
           (mask & kLeafChannelPdf) != 0U);
}

inline ExactRelation exact_source_program_relation(
    const CompiledMathSourceProductProgram &program) noexcept {
  return static_cast<ExactRelation>(program.static_source_view_relation);
}

inline ExactSourceFill exact_source_conditionalize(
    const ExactSourceFill unconditioned,
    const ExactSourceFill lower,
    const std::uint8_t mask) {
  if (!std::isfinite(lower.survival) || !(lower.survival > 0.0)) {
    return exact_source_impossible_fill();
  }
  ExactSourceFill out;
  if ((mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(unconditioned.pdf / lower.survival);
  }
  if ((mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability(
        (unconditioned.cdf + lower.survival - 1.0) / lower.survival);
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability(
        unconditioned.survival / lower.survival);
  }
  return out;
}

inline ExactSourceFill exact_source_conditionalize_between(
    const ExactSourceFill unconditioned,
    const ExactSourceFill lower,
    const ExactSourceFill upper,
    const std::uint8_t mask) {
  const double mass = upper.cdf - lower.cdf;
  if (!std::isfinite(mass) || !(mass > 0.0)) {
    return exact_source_impossible_fill();
  }
  ExactSourceFill out;
  if ((mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(unconditioned.pdf / mass);
  }
  if ((mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability((unconditioned.cdf - lower.cdf) / mass);
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability((upper.cdf - unconditioned.cdf) / mass);
  }
  return out;
}

} // namespace detail
} // namespace accumulatr::eval
