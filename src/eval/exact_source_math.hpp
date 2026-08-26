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

inline ExactSourceFill exact_source_finish_base_fill(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t mask) {
  const double start_probability = 1.0 - q;
  ExactSourceFill fill;
  if ((mask & kLeafChannelPdf) != 0U) {
    fill.pdf = safe_density(start_probability * base_pdf);
  }
  if ((mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    const double cdf = clamp_probability(start_probability * base_cdf);
    if ((mask & kLeafChannelCdf) != 0U) {
      fill.cdf = cdf;
    }
    if ((mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0 - cdf;
    }
  }
  return fill;
}

inline ExactSourceFill exact_source_lognormal_leaf_fill(
    const double meanlog,
    const double sdlog,
    const double q,
    const double x,
    const std::uint8_t mask) {
  if (!std::isfinite(meanlog) || !std::isfinite(sdlog) ||
      !(sdlog > 0.0)) {
    return {};
  }
  const bool pdf = (mask & kLeafChannelPdf) != 0U;
  const bool cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double z = (std::log(x) - meanlog) / sdlog;
  return exact_source_finish_base_fill(
      pdf ? normal_pdf_fast(z) / (x * sdlog) : 0.0,
      cdf ? normal_cdf_fast(z) : 0.0,
      q,
      mask);
}

inline ExactSourceFill exact_source_gamma_leaf_fill(
    const double shape,
    const double rate,
    const double q,
    const double x,
    const std::uint8_t mask) {
  const bool pdf = (mask & kLeafChannelPdf) != 0U;
  const bool cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double scale = 1.0 / rate;
  return exact_source_finish_base_fill(
      pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
      cdf ? R::pgamma(x, shape, scale, 1, 0) : 0.0,
      q,
      mask);
}

inline ExactSourceFill exact_source_exgauss_leaf_fill(
    const double mu,
    const double sigma,
    const double tau,
    const double q,
    const double x,
    const std::uint8_t mask) {
  const bool pdf = (mask & kLeafChannelPdf) != 0U;
  const bool cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
  const double lower_survival = 1.0 - lower_cdf;
  if (!(lower_survival > 0.0)) {
    ExactSourceFill fill;
    return fill;
  }
  return exact_source_finish_base_fill(
      pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
      cdf ? clamp_probability(
                (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                lower_survival)
          : 0.0,
      q,
      mask);
}

inline ExactSourceFill exact_source_lba_leaf_fill(
    const double v,
    const double B,
    const double A,
    const double sv,
    const double q,
    const double x,
    const std::uint8_t mask) {
  const bool pdf = (mask & kLeafChannelPdf) != 0U;
  const bool cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return exact_source_finish_base_fill(
      pdf ? lba_pdf_fast(x, v, B, A, sv) : 0.0,
      cdf ? lba_cdf_fast(x, v, B, A, sv) : 0.0,
      q,
      mask);
}

inline ExactSourceFill exact_source_rdm_leaf_fill(
    const double v,
    const double B,
    const double A,
    const double s,
    const double q,
    const double x,
    const std::uint8_t mask) {
  const bool pdf = (mask & kLeafChannelPdf) != 0U;
  const bool cdf = (mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return exact_source_finish_base_fill(
      pdf ? rdm_pdf_fast(x, v, B, A, s) : 0.0,
      cdf ? rdm_cdf_fast(x, v, B, A, s) : 0.0,
      q,
      mask);
}

inline ExactSourceFill exact_source_impossible_fill(
    const std::uint8_t mask) {
  ExactSourceFill fill;
  return fill;
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
  return exact_source_impossible_fill(mask);
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
  return program.has_static_source_view_relation
             ? static_cast<ExactRelation>(program.static_source_view_relation)
             : ExactRelation::Unknown;
}

inline ExactSourceFill exact_source_conditionalize(
    const ExactSourceFill unconditioned,
    const ExactSourceFill lower,
    const std::uint8_t mask) {
  if (!std::isfinite(lower.survival) || !(lower.survival > 0.0)) {
    return exact_source_impossible_fill(mask);
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
    return exact_source_impossible_fill(mask);
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
