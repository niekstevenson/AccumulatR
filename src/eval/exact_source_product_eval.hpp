#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "../leaf/dist_kind.hpp"
#include "exact_eval_workspace.hpp"
#include "leaf_kernel.hpp"
#include "quadrature.hpp"

namespace accumulatr::eval {
namespace detail {

inline double compiled_math_source_product_channel_time(
    const CompiledMathSourceProductChannel &channel,
    const CompiledMathWorkspace &workspace) {
  const double time =
      workspace.time_values[static_cast<std::size_t>(channel.time_id)];
  if (channel.time_cap_id == semantic::kInvalidIndex) {
    return time;
  }
  return std::min(
      time,
      workspace.time_values[static_cast<std::size_t>(channel.time_cap_id)]);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_finish_base_fill(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t fill_mask) {
  const double start_prob = 1.0 - q;
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    fill.pdf = start_prob * safe_density(base_pdf);
  }
  if ((fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    const double cdf = clamp_probability(start_prob * clamp_probability(base_cdf));
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      fill.cdf = cdf;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = clamp_probability(1.0 - cdf);
    }
  }
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_lognormal_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double m = loaded.params[0];
  const double s = loaded.params[1];
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
      need_cdf ? R::plnorm(x, m, s, 1, 0) : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_gamma_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double shape = loaded.params[0];
  const double scale = 1.0 / loaded.params[1];
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
      need_cdf ? R::pgamma(x, shape, scale, 1, 0) : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_exgauss_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const double mu = loaded.params[0];
  const double sigma = loaded.params[1];
  const double tau = loaded.params[2];
  const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
  const double lower_survival = 1.0 - lower_cdf;
  if (!(lower_survival > 0.0)) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  return compiled_math_source_product_finish_base_fill(
      need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
      need_cdf
          ? (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                lower_survival
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_lba_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return compiled_math_source_product_finish_base_fill(
      need_pdf
          ? lba_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_cdf
          ? lba_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_rdm_leaf_fill(
    const ExactLoadedLeafInput &loaded,
    const double x,
    const std::uint8_t fill_mask) {
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  return compiled_math_source_product_finish_base_fill(
      need_pdf
          ? rdm_pdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      need_cdf
          ? rdm_cdf_fast(x, loaded.params[0], loaded.params[1],
                         loaded.params[2], loaded.params[3])
          : 0.0,
      loaded.q,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_direct_leaf_fill(
    const std::uint8_t leaf_dist_kind,
    const ExactLoadedLeafInput &loaded,
    const double time_since_onset,
    const std::uint8_t fill_mask) {
  const double x = time_since_onset - loaded.t0;
  if (!(x > 0.0)) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = fill_mask;
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  switch (static_cast<leaf::DistKind>(leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    return compiled_math_source_product_lognormal_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::Gamma:
    return compiled_math_source_product_gamma_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::Exgauss:
    return compiled_math_source_product_exgauss_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::LBA:
    return compiled_math_source_product_lba_leaf_fill(
        loaded, x, fill_mask);
  case leaf::DistKind::RDM:
    return compiled_math_source_product_rdm_leaf_fill(
        loaded, x, fill_mask);
  }
  return ExactSourceChannels::SourceProductScalarFill{};
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_direct_leaf_fill(
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t fill_mask) {
  return compiled_math_source_product_direct_leaf_fill(
      channel.leaf_dist_kind,
      loaded,
      time - loaded.onset_abs,
      fill_mask);
}

inline double compiled_math_source_product_finish_base_scalar(
    const double base_pdf,
    const double base_cdf,
    const double q,
    const std::uint8_t channel_mask) {
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

inline double compiled_math_source_product_direct_leaf_scalar(
    const CompiledMathSourceProductChannel &channel,
    const ExactLoadedLeafInput &loaded,
    const double time,
    const std::uint8_t channel_mask) {
  const double x = time - loaded.onset_abs - loaded.t0;
  if (!(x > 0.0)) {
    return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
  }
  const bool need_pdf = channel_mask == kLeafChannelPdf;
  switch (static_cast<leaf::DistKind>(channel.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal: {
    const double m = loaded.params[0];
    const double s = loaded.params[1];
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? R::dlnorm(x, m, s, 0) : 0.0,
        need_pdf ? 0.0 : R::plnorm(x, m, s, 1, 0),
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::Gamma: {
    const double shape = loaded.params[0];
    const double scale = 1.0 / loaded.params[1];
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? R::dgamma(x, shape, scale, 0) : 0.0,
        need_pdf ? 0.0 : R::pgamma(x, shape, scale, 1, 0),
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::Exgauss: {
    const double mu = loaded.params[0];
    const double sigma = loaded.params[1];
    const double tau = loaded.params[2];
    const double lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      return channel_mask == kLeafChannelSurvival ? 1.0 : 0.0;
    }
    return compiled_math_source_product_finish_base_scalar(
        need_pdf ? exgauss_raw_pdf(x, mu, sigma, tau) / lower_survival : 0.0,
        need_pdf
            ? 0.0
            : (exgauss_raw_cdf(x, mu, sigma, tau) - lower_cdf) /
                  lower_survival,
        loaded.q,
        channel_mask);
  }
  case leaf::DistKind::LBA:
    return compiled_math_source_product_finish_base_scalar(
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
  case leaf::DistKind::RDM:
    return compiled_math_source_product_finish_base_scalar(
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
  return 0.0;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_impossible_fill(const std::uint8_t mask) {
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = mask;
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_certain_fill(const std::uint8_t mask) {
  ExactSourceChannels::SourceProductScalarFill fill;
  fill.mask = mask;
  if ((mask & kLeafChannelCdf) != 0U) {
    fill.cdf = 1.0;
  }
  if ((mask & kLeafChannelSurvival) != 0U) {
    fill.survival = 0.0;
  }
  return fill;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_forced_fill(
    const ExactRelation relation,
    const std::uint8_t mask) {
  if (relation == ExactRelation::Before || relation == ExactRelation::At) {
    return compiled_math_source_product_certain_fill(mask);
  }
  if (relation == ExactRelation::After) {
    ExactSourceChannels::SourceProductScalarFill fill;
    fill.mask = mask;
    if ((mask & kLeafChannelSurvival) != 0U) {
      fill.survival = 1.0;
    }
    return fill;
  }
  return compiled_math_source_product_impossible_fill(mask);
}

inline bool compiled_math_source_product_bounds_have_overlay(
    const CompiledSourceBoundPlan &bounds) noexcept {
  return bounds.has_condition_exact ||
         bounds.has_condition_lower ||
         bounds.has_condition_upper;
}

inline bool compiled_math_source_product_relation_forces_fill(
    const ExactRelation relation,
    const std::uint8_t fill_mask) noexcept {
  return relation != ExactRelation::Unknown &&
         !(relation == ExactRelation::At &&
           (fill_mask & kLeafChannelPdf) != 0U);
}

inline ExactRelation compiled_math_source_product_program_relation(
    const CompiledMathSourceProductProgram &source_program) noexcept {
  return source_program.has_static_source_view_relation
             ? static_cast<ExactRelation>(
                   source_program.static_source_view_relation)
             : ExactRelation::Unknown;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_conditionalize(
    const ExactSourceChannels::SourceProductScalarFill uncond,
    const ExactSourceChannels::SourceProductScalarFill lower,
    const std::uint8_t fill_mask) {
  const double surv_lb = lower.survival;
  if (!std::isfinite(surv_lb) || !(surv_lb > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(uncond.pdf / surv_lb);
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability((uncond.cdf + surv_lb - 1.0) / surv_lb);
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability(uncond.survival / surv_lb);
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_conditionalize_between(
    const ExactSourceChannels::SourceProductScalarFill uncond,
    const ExactSourceChannels::SourceProductScalarFill lower,
    const ExactSourceChannels::SourceProductScalarFill upper,
    const std::uint8_t fill_mask) {
  const double mass = upper.cdf - lower.cdf;
  if (!std::isfinite(mass) || !(mass > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if ((fill_mask & kLeafChannelPdf) != 0U) {
    out.pdf = safe_density(uncond.pdf / mass);
  }
  if ((fill_mask & kLeafChannelCdf) != 0U) {
    out.cdf = clamp_probability((uncond.cdf - lower.cdf) / mass);
  }
  if ((fill_mask & kLeafChannelSurvival) != 0U) {
    out.survival = clamp_probability((upper.cdf - uncond.cdf) / mass);
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask);

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_leaf_fill(
    const CompiledMathSourceProductProgram &source_program,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const auto &loaded =
      source_channels->source_product_leaf_input(source_program.leaf_index);
  return compiled_math_source_product_direct_leaf_fill(
      source_program.leaf_dist_kind,
      loaded,
      current_time - loaded.onset_abs,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_exact_gate_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !compiled_math_source_product_bounds_have_overlay(
          source_program.bounds)) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }
  const auto bounds =
      source_channels->source_product_resolved_bounds(
          source_program.source_id,
          source_program.bounds,
          workspace);
  if (const double *exact_time = bounds.exact_time) {
    if (!(current_time >= *exact_time)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    return compiled_math_source_product_certain_fill(fill_mask);
  }
  return compiled_math_source_product_program_fill(
      program,
      source_program.child_program_id,
      workspace,
      source_channels,
      current_time,
      fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_conditioned_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const auto relation =
      compiled_math_source_product_program_relation(source_program);
  const bool has_condition_overlay =
      compiled_math_source_product_bounds_have_overlay(
          source_program.bounds);
  if (compiled_math_source_product_relation_forces_fill(relation, fill_mask) &&
      !has_condition_overlay) {
    return compiled_math_source_product_forced_fill(relation, fill_mask);
  }

  const auto bounds =
      source_channels->source_product_resolved_bounds(
          source_program.source_id,
          source_program.bounds,
          workspace);
  const double lower_bound = bounds.lower;
  const double upper_bound = bounds.upper;
  if (std::isfinite(upper_bound) && !(upper_bound > lower_bound)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  if (const double *exact_time = bounds.exact_time) {
    if (lower_bound > 0.0 && !(*exact_time > lower_bound)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    if (std::isfinite(upper_bound) && !(*exact_time < upper_bound)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    if (!(current_time >= *exact_time)) {
      return compiled_math_source_product_impossible_fill(fill_mask);
    }
    return compiled_math_source_product_certain_fill(fill_mask);
  }

  std::uint8_t uncond_mask = fill_mask;
  if (std::isfinite(upper_bound) &&
      (fill_mask & kLeafChannelSurvival) != 0U) {
    uncond_mask |= kLeafChannelCdf;
  }
  const auto uncond =
      compiled_math_source_product_program_fill(
          program,
          source_program.child_program_id,
          workspace,
          source_channels,
          current_time,
          uncond_mask);
  if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
    return uncond;
  }

  std::uint8_t lower_mask = kLeafChannelSurvival;
  if (std::isfinite(upper_bound)) {
    lower_mask = kLeafChannelCdf;
  }
  auto lower = compiled_math_source_product_impossible_fill(lower_mask);
  if (lower_bound > 0.0) {
    lower =
        compiled_math_source_product_program_fill(
            program,
            source_program.child_program_id,
            workspace,
            source_channels,
            lower_bound,
            lower_mask);
  }
  if (!std::isfinite(upper_bound)) {
    return compiled_math_source_product_conditionalize(
        uncond, lower, fill_mask);
  }
  const auto upper =
      compiled_math_source_product_program_fill(
          program,
          source_program.child_program_id,
          workspace,
          source_channels,
          upper_bound,
          kLeafChannelCdf);
  if (!std::isfinite(upper.cdf - lower.cdf) ||
      !(upper.cdf - lower.cdf > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  if (current_time >= upper_bound) {
    return compiled_math_source_product_certain_fill(fill_mask);
  }
  if (current_time <= lower_bound) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  return compiled_math_source_product_conditionalize_between(
      uncond, lower, upper, fill_mask);
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_onset_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const double upper = current_time - source_program.leaf_onset_lag;
  if (!(upper > 0.0)) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }

  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t shifted_mask =
      (need_pdf ? kLeafChannelPdf : 0U) |
      (need_cdf ? kLeafChannelCdf : 0U);
  const auto &loaded =
      source_channels->source_product_leaf_input(source_program.leaf_index);
  auto shifted_fill = [&](const double source_time) {
    return compiled_math_source_product_direct_leaf_fill(
        source_program.leaf_dist_kind,
        loaded,
        current_time - source_time - source_program.leaf_onset_lag,
        shifted_mask);
  };

  const auto &onset_program =
      program.integral_kernel_source_product_programs[
          static_cast<std::size_t>(
              source_program.onset_source_program_id)];
  const auto onset_bounds =
      source_channels->source_product_resolved_bounds(
          onset_program.source_id,
          source_program.onset_bounds,
          workspace);
  if (const double *exact_time = onset_bounds.exact_time) {
    const auto shifted = shifted_fill(*exact_time);
    ExactSourceChannels::SourceProductScalarFill out;
    out.mask = fill_mask;
    if (need_pdf) {
      out.pdf = shifted.pdf;
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = shifted.cdf;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(1.0 - shifted.cdf);
    }
    return out;
  }

  const auto batch = quadrature::build_finite_batch(0.0, upper);
  double pdf = 0.0;
  double cdf = 0.0;
  for (std::size_t i = 0; i < batch.nodes.nodes.size(); ++i) {
    const double u = batch.nodes.nodes[i];
    const auto source =
        compiled_math_source_product_program_fill(
            program,
            source_program.onset_source_program_id,
            workspace,
            source_channels,
            u,
            kLeafChannelPdf);
    if (!(source.pdf > 0.0)) {
      continue;
    }
    const auto shifted = shifted_fill(u);
    const double weight = batch.nodes.weights[i] * source.pdf;
    if (need_pdf) {
      pdf += weight * shifted.pdf;
    }
    if (need_cdf) {
      cdf += weight * shifted.cdf;
    }
  }

  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if (need_pdf) {
    out.pdf = safe_density(pdf);
  }
  if (need_cdf) {
    const double clamped = clamp_probability(cdf);
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamped;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamp_probability(1.0 - clamped);
    }
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_pool_fill(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  const int n_members = static_cast<int>(source_program.member_programs.size);
  const int k = static_cast<int>(source_program.pool_k);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const std::uint8_t member_mask =
      kLeafChannelCdf | kLeafChannelSurvival |
      (need_pdf ? kLeafChannelPdf : 0U);

  double *scratch =
      workspace->source_product_scratch.data() +
      static_cast<std::size_t>(
          source_program.source_product_scratch_offset);
  double *member_pdf = scratch;
  double *member_cdf = member_pdf + n_members;
  double *member_survival = member_cdf + n_members;
  const int width = n_members + 1;
  const std::size_t table_size =
      static_cast<std::size_t>(width) * static_cast<std::size_t>(width);
  double *prefix = member_survival + n_members;
  double *suffix = prefix + table_size;
  std::fill(prefix, prefix + table_size, 0.0);
  std::fill(suffix, suffix + table_size, 0.0);

  for (semantic::Index i = 0; i < source_program.member_programs.size; ++i) {
    const auto member_program_id =
        program.integral_kernel_source_product_program_members[
            static_cast<std::size_t>(
                source_program.member_programs.offset + i)];
    const auto member =
        compiled_math_source_product_program_fill(
            program,
            member_program_id,
            workspace,
            source_channels,
            current_time,
            member_mask);
    const auto pos = static_cast<std::size_t>(i);
    member_pdf[pos] = member.pdf;
    member_cdf[pos] = member.cdf;
    member_survival[pos] = member.survival;
  }

  const auto idx = [width](const int row, const int col) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(col);
  };
  prefix[idx(0, 0)] = 1.0;
  for (int i = 0; i < n_members; ++i) {
    for (int m = 0; m <= i; ++m) {
      prefix[idx(i + 1, m)] +=
          prefix[idx(i, m)] * member_survival[static_cast<std::size_t>(i)];
      prefix[idx(i + 1, m + 1)] +=
          prefix[idx(i, m)] * member_cdf[static_cast<std::size_t>(i)];
    }
  }
  suffix[idx(n_members, 0)] = 1.0;
  for (int i = n_members - 1; i >= 0; --i) {
    const int count = n_members - i - 1;
    for (int m = 0; m <= count; ++m) {
      suffix[idx(i, m)] +=
          suffix[idx(i + 1, m)] *
          member_survival[static_cast<std::size_t>(i)];
      suffix[idx(i, m + 1)] +=
          suffix[idx(i + 1, m)] *
          member_cdf[static_cast<std::size_t>(i)];
    }
  }

  double pdf = 0.0;
  if (need_pdf) {
    for (int i = 0; i < n_members; ++i) {
      double others_exact = 0.0;
      for (int left = 0; left < k; ++left) {
        const int right = (k - 1) - left;
        if (right < 0 || right > (n_members - i - 1)) {
          continue;
        }
        others_exact +=
            prefix[idx(i, left)] * suffix[idx(i + 1, right)];
      }
      pdf += member_pdf[static_cast<std::size_t>(i)] * others_exact;
    }
  }
  double surv = 0.0;
  if (need_cdf) {
    for (int m = 0; m < k; ++m) {
      surv += prefix[idx(n_members, m)];
    }
  }

  ExactSourceChannels::SourceProductScalarFill out;
  out.mask = fill_mask;
  if (need_pdf) {
    out.pdf = safe_density(pdf);
  }
  if (need_cdf) {
    const double clamped_survival = clamp_probability(surv);
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      out.survival = clamped_survival;
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      out.cdf = clamp_probability(1.0 - clamped_survival);
    }
  }
  return out;
}

inline ExactSourceChannels::SourceProductScalarFill
compiled_math_source_product_program_fill(
    const CompiledMathProgram &program,
    const semantic::Index source_product_program_id,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels,
    const double current_time,
    const std::uint8_t fill_mask) {
  if (source_product_program_id == semantic::kInvalidIndex) {
    return compiled_math_source_product_impossible_fill(fill_mask);
  }
  const auto &source_program =
      program.integral_kernel_source_product_programs[
          static_cast<std::size_t>(source_product_program_id)];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    return compiled_math_source_product_impossible_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::ConstantOne:
    return compiled_math_source_product_certain_fill(fill_mask);
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    return compiled_math_source_product_program_leaf_fill(
        source_program, source_channels, current_time, fill_mask);
  case CompiledMathSourceProductProgramKind::ExactGate:
    return compiled_math_source_product_program_exact_gate_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::Conditioned:
    return compiled_math_source_product_program_conditioned_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return compiled_math_source_product_program_onset_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return compiled_math_source_product_program_pool_fill(
        program,
        source_program,
        workspace,
        source_channels,
        current_time,
        fill_mask);
  }
  return compiled_math_source_product_impossible_fill(fill_mask);
}

inline void compiled_math_store_source_product_program_fill(
    CompiledMathWorkspace *workspace,
    const std::size_t program_pos,
    const ExactSourceChannels::SourceProductScalarFill &fill) {
  if (workspace->source_product_program_epoch[program_pos] !=
      workspace->source_product_program_current_epoch) {
    workspace->source_product_program_epoch[program_pos] =
        workspace->source_product_program_current_epoch;
    workspace->source_product_program_valid_mask[program_pos] = 0U;
  }
  workspace->source_product_program_valid_mask[program_pos] |= fill.mask;
  if ((fill.mask & kLeafChannelPdf) != 0U) {
    workspace->source_product_program_pdf[program_pos] = fill.pdf;
  }
  if ((fill.mask & kLeafChannelCdf) != 0U) {
    workspace->source_product_program_cdf[program_pos] = fill.cdf;
  }
  if ((fill.mask & kLeafChannelSurvival) != 0U) {
    workspace->source_product_program_survival[program_pos] = fill.survival;
  }
}

inline double compiled_math_source_product_fill_value(
    const ExactSourceChannels::SourceProductScalarFill &fill,
    const std::uint8_t channel_mask) {
  if (channel_mask == kLeafChannelPdf) {
    return fill.pdf;
  }
  if (channel_mask == kLeafChannelCdf) {
    return fill.cdf;
  }
  if (channel_mask == kLeafChannelSurvival) {
    return fill.survival;
  }
  return 0.0;
}

inline std::uint8_t compiled_math_source_node_fill_mask(
    const CompiledMathNodeKind kind) noexcept {
  const auto value_mask = compiled_math_source_factor_channel_mask(kind);
  if ((value_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U) {
    return kLeafChannelCdf | kLeafChannelSurvival;
  }
  return value_mask;
}

inline double compiled_math_source_node_value(
    const CompiledMathProgram &program,
    const CompiledMathNode &node,
    CompiledSourceView *evaluator,
    CompiledMathWorkspace *workspace) {
  auto *source_channels = evaluator->source_channels();
  const auto program_pos =
      static_cast<std::size_t>(node.source_program_id);
  const auto value_mask = compiled_math_source_factor_channel_mask(node.kind);
  const auto fill_mask = compiled_math_source_node_fill_mask(node.kind);
  if (workspace->source_product_program_epoch[program_pos] !=
          workspace->source_product_program_current_epoch ||
      (workspace->source_product_program_valid_mask[program_pos] &
       fill_mask) != fill_mask) {
    const auto fill =
        compiled_math_source_product_program_fill(
            program,
            node.source_program_id,
            workspace,
            source_channels,
            compiled_math_source_node_time(node, *workspace),
            fill_mask);
    compiled_math_store_source_product_program_fill(
        workspace,
        program_pos,
        fill);
  }
  if (value_mask == kLeafChannelPdf) {
    return safe_density(workspace->source_product_program_pdf[program_pos]);
  }
  if (value_mask == kLeafChannelCdf) {
    return clamp_probability(
        workspace->source_product_program_cdf[program_pos]);
  }
  if (value_mask == kLeafChannelSurvival) {
    return clamp_probability(
        workspace->source_product_program_survival[program_pos]);
  }
  return 0.0;
}

inline double compiled_math_source_product_value_for_ops(
    const CompiledMathProgram &program,
    const CompiledMathIndexSpan source_product_ops,
    CompiledMathWorkspace *workspace,
    ExactSourceChannels *source_channels) {
  double product = 1.0;
  for (semantic::Index i = 0; i < source_product_ops.size; ++i) {
    const auto &op =
        program.integral_kernel_source_product_ops[
            static_cast<std::size_t>(source_product_ops.offset + i)];
    double value = op.constant_value;
    const auto channel_mask = op.value_channel_mask;
    if (channel_mask == 0U) {
      if (value == 0.0) {
        return 0.0;
      }
    } else {
      const auto channel_pos =
          static_cast<std::size_t>(op.source_product_channel_id);
      const auto &channel =
          program.integral_kernel_source_product_channels[channel_pos];
      if (!op.cache_result) {
        const double time =
            compiled_math_source_product_channel_time(channel, *workspace);
        if (source_channels->source_product_direct_leaf_available(
                op.source_product_channel_id)) {
          value = compiled_math_source_product_direct_leaf_scalar(
              channel,
              source_channels->source_product_leaf_input(channel.leaf_index),
              time,
              channel_mask);
        } else {
          const auto fill =
              compiled_math_source_product_program_fill(
                  program,
                  op.source_product_program_id,
                  workspace,
                  source_channels,
                  time,
                  channel_mask);
          value = compiled_math_source_product_fill_value(fill, channel_mask);
        }
      } else {
        const auto program_pos =
            static_cast<std::size_t>(op.source_product_program_id);
        if (workspace->source_product_program_epoch[program_pos] !=
                workspace->source_product_program_current_epoch ||
            (workspace->source_product_program_valid_mask[program_pos] &
             channel_mask) == 0U) {
          const double time =
              compiled_math_source_product_channel_time(channel, *workspace);
          const auto fill_mask = op.fill_channel_mask;
          const auto fill =
              source_channels->source_product_direct_leaf_available(
                  op.source_product_channel_id)
                  ? compiled_math_source_product_direct_leaf_fill(
                        channel,
                        source_channels->source_product_leaf_input(
                            channel.leaf_index),
                        time,
                        fill_mask)
                  : compiled_math_source_product_program_fill(
                        program,
                        op.source_product_program_id,
                        workspace,
                        source_channels,
                        time,
                        fill_mask);
          compiled_math_store_source_product_program_fill(
              workspace,
              program_pos,
              fill);
        }
        if (channel_mask == kLeafChannelPdf) {
          value = workspace->source_product_program_pdf[program_pos];
        } else if (channel_mask == kLeafChannelCdf) {
          value = workspace->source_product_program_cdf[program_pos];
        } else {
          value = workspace->source_product_program_survival[program_pos];
        }
      }
    }
    product *= value;
    if (product == 0.0) {
      return 0.0;
    }
  }
  return product;
}

} // namespace detail
} // namespace accumulatr::eval
