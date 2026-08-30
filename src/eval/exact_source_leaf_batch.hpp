#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "exact_source_math.hpp"
#include "lane_math.hpp"

namespace accumulatr::eval {
namespace detail {

struct SourceLaneFill {
  void resize(const std::size_t count, const std::uint8_t requested_mask) {
    mask = requested_mask;
    if ((mask & kLeafChannelPdf) != 0U) {
      pdf.resize(count);
    }
    if ((mask & kLeafChannelCdf) != 0U) {
      cdf.resize(count);
    }
    if ((mask & kLeafChannelSurvival) != 0U) {
      survival.resize(count);
    }
  }

  void assign(const std::size_t count, const std::uint8_t requested_mask) {
    resize(count, requested_mask);
    if ((mask & kLeafChannelPdf) != 0U) {
      std::fill(pdf.begin(), pdf.end(), 0.0);
    }
    if ((mask & kLeafChannelCdf) != 0U) {
      std::fill(cdf.begin(), cdf.end(), 0.0);
    }
    if ((mask & kLeafChannelSurvival) != 0U) {
      std::fill(survival.begin(), survival.end(), 1.0);
    }
  }

  std::uint8_t mask{0U};
  std::vector<double> pdf;
  std::vector<double> cdf;
  std::vector<double> survival;
};

struct PreparedSourceLeafBatch {
  void ensure_lanes(const std::size_t lane_count) {
    count = lane_count;
    const auto value_count =
        (static_cast<std::size_t>(leaf::kMaxDistParamCount) + 2U) * lane_count;
    if (values.size() < value_count) {
      values.resize(value_count);
    }
    if (positions.size() < lane_count) {
      positions.resize(lane_count);
    }
  }

  double *elapsed() noexcept {
    return values.data();
  }

  const double *elapsed() const noexcept {
    return values.data();
  }

  double *q() noexcept {
    return values.data() + count;
  }

  const double *q() const noexcept {
    return values.data() + count;
  }

  double *parameter(const std::size_t slot) noexcept {
    return values.data() + (slot + 2U) * count;
  }

  const double *parameter(const std::size_t slot) const noexcept {
    return values.data() + (slot + 2U) * count;
  }

  void ensure_work(const std::size_t lane_count,
                   const std::size_t plane_count) {
    work_stride = lane_count;
    const auto value_count = plane_count * lane_count;
    if (work.size() < value_count) {
      work.resize(value_count);
    }
  }

  double *work_plane(const std::size_t plane) noexcept {
    return work.data() + plane * work_stride;
  }

  std::size_t count{0U};
  std::size_t work_stride{0U};
  std::vector<double> values;
  std::vector<double> work;
  std::vector<std::size_t> positions;
};

inline void source_lane_store_fill(
    SourceLaneFill *out,
    const std::size_t position,
    const ExactSourceFill &fill) {
  if ((out->mask & kLeafChannelPdf) != 0U) {
    out->pdf[position] = fill.pdf;
  }
  if ((out->mask & kLeafChannelCdf) != 0U) {
    out->cdf[position] = fill.cdf;
  }
  if ((out->mask & kLeafChannelSurvival) != 0U) {
    out->survival[position] = fill.survival;
  }
}

template <std::uint8_t Mask>
inline void source_lane_store_fill_masked(
    SourceLaneFill *out,
    const std::size_t position,
    const ExactSourceFill &fill) {
  if constexpr ((Mask & kLeafChannelPdf) != 0U) {
    out->pdf[position] = fill.pdf;
  }
  if constexpr ((Mask & kLeafChannelCdf) != 0U) {
    out->cdf[position] = fill.cdf;
  }
  if constexpr ((Mask & kLeafChannelSurvival) != 0U) {
    out->survival[position] = fill.survival;
  }
}

template <bool NeedCdf>
inline void prepare_normal_lanes(
    const double *z,
    const std::size_t count,
    double *factors,
    double *exponent_arguments) {
  for (std::size_t i = 0U; i < count; ++i) {
    const double value = z[i];
    const double abs_z = std::fabs(value);
    if constexpr (NeedCdf) {
      double factor = abs_z > 37.0
                          ? 0.0
                          : (abs_z < kNormalHartSplit
                                 ? normal_hart_ratio(abs_z)
                                 : 1.0 / (kSqrtTwoPi *
                                          normal_hart_fraction(abs_z)));
      factors[i] = std::copysign(factor, value <= 0.0 ? 1.0 : -1.0);
    }
    exponent_arguments[i] = -0.5 * abs_z * abs_z;
  }
}

inline double finish_normal_cdf(const double factor,
                                const double exponential) noexcept {
  const double tail = std::fabs(factor) * exponential;
  return std::signbit(factor) ? 1.0 - tail : tail;
}

inline void finish_normal_cdf_lanes(
    const std::size_t count,
    const double *factors,
    const double *exponentials,
    double *cdf) {
  for (std::size_t i = 0U; i < count; ++i) {
    cdf[i] = finish_normal_cdf(factors[i], exponentials[i]);
  }
}

template <std::uint8_t Mask>
inline void evaluate_lognormal_leaf_batch(
    PreparedSourceLeafBatch *input,
    SourceLaneFill *out) {
  constexpr bool need_pdf = (Mask & kLeafChannelPdf) != 0U;
  constexpr bool need_cdf =
      (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  input->ensure_work(input->count, 5U);
  auto *z = input->work_plane(0U);
  auto *factors = input->work_plane(1U);
  auto *exponent_arguments = input->work_plane(2U);
  auto *exponentials = input->work_plane(3U);
  auto *cdf = input->work_plane(4U);
  log_lanes(input->elapsed(), z, input->count);
  const auto *meanlog = input->parameter(0U);
  const auto *sdlog = input->parameter(1U);
  for (std::size_t i = 0U; i < input->count; ++i) {
    const bool valid = input->elapsed()[i] > 0.0 &&
                       std::isfinite(meanlog[i]) &&
                       std::isfinite(sdlog[i]) && sdlog[i] > 0.0;
    z[i] = valid ? (z[i] - meanlog[i]) / sdlog[i] : 0.0;
  }
  prepare_normal_lanes<need_cdf>(
      z, input->count, factors, exponent_arguments);
  exp_lanes(exponent_arguments, exponentials, input->count);
  if constexpr (need_cdf) {
    finish_normal_cdf_lanes(
        input->count, factors, exponentials, cdf);
  }
  const auto *q = input->q();
  for (std::size_t i = 0U; i < input->count; ++i) {
    auto fill = exact_source_impossible_fill(Mask);
    const double x = input->elapsed()[i];
    if (x > 0.0 && std::isfinite(meanlog[i]) &&
        std::isfinite(sdlog[i]) && sdlog[i] > 0.0) {
      fill = exact_source_finish_base_fill<Mask>(
          need_pdf ? kInverseSqrtTwoPi * exponentials[i] /
                         (x * sdlog[i])
                   : 0.0,
          need_cdf ? cdf[i] : 0.0,
          q[i]);
    }
    source_lane_store_fill_masked<Mask>(out, i, fill);
  }
}

template <std::uint8_t Mask, bool WideStart>
inline void evaluate_lba_leaf_group(
    PreparedSourceLeafBatch *input,
    const std::size_t *positions,
    const std::size_t count,
    SourceLaneFill *out) {
  constexpr bool need_pdf = (Mask & kLeafChannelPdf) != 0U;
  constexpr bool need_cdf =
      (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  constexpr std::size_t normal_count = WideStart ? 3U : 2U;
  input->ensure_work(count, 4U * normal_count);
  auto *z = input->work_plane(0U);
  auto *factors = input->work_plane(normal_count);
  auto *exponentials = input->work_plane(2U * normal_count);
  auto *cdf = input->work_plane(3U * normal_count);
  const auto *x = input->elapsed();
  const auto *v = input->parameter(0U);
  const auto *B = input->parameter(1U);
  const auto *A = input->parameter(2U);
  const auto *sv = input->parameter(3U);
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    z[lane] = v[i] / sv[i];
    if constexpr (WideStart) {
      const double zs = x[i] * sv[i];
      const double cmz = B[i] - x[i] * v[i];
      z[count + lane] = cmz / zs;
      z[2U * count + lane] = (cmz - A[i]) / zs;
    } else {
      z[count + lane] = (B[i] / x[i] - v[i]) / sv[i];
    }
  }
  prepare_normal_lanes<true>(z, count, factors, z);
  if constexpr (WideStart) {
    prepare_normal_lanes<true>(
        z + count, 2U * count, factors + count, z + count);
  } else {
    prepare_normal_lanes<need_cdf>(
        z + count, count, factors + count, z + count);
  }
  exp_lanes(z, exponentials, normal_count * count);
  finish_normal_cdf_lanes(count, factors, exponentials, cdf);
  if constexpr (WideStart) {
    finish_normal_cdf_lanes(
        2U * count, factors + count, exponentials + count,
        cdf + count);
  } else if constexpr (need_cdf) {
    finish_normal_cdf_lanes(
        count, factors + count, exponentials + count,
        cdf + count);
  }
  const auto *q = input->q();
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    double denom = cdf[lane];
    if (!std::isfinite(denom) || denom < 1e-10) {
      denom = 1e-10;
    }
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (WideStart) {
      const double zs = x[i] * sv[i];
      const double cmz = B[i] - x[i] * v[i];
      const double xx = cmz - A[i];
      const double cdf_z = cdf[count + lane];
      const double cdf_z_max = cdf[2U * count + lane];
      const double pdf_z =
          kInverseSqrtTwoPi * exponentials[count + lane];
      const double pdf_z_max =
          kInverseSqrtTwoPi * exponentials[2U * count + lane];
      if constexpr (need_pdf) {
        base_pdf =
            (v[i] * (cdf_z - cdf_z_max) +
             sv[i] * (pdf_z_max - pdf_z)) /
            (A[i] * denom);
      }
      if constexpr (need_cdf) {
        base_cdf = clamp_probability(
            (1.0 + (zs * (pdf_z_max - pdf_z) +
                    xx * cdf_z_max - cmz * cdf_z) /
                       A[i]) /
            denom);
      }
    } else {
      if constexpr (need_pdf) {
        const double normal_pdf =
            kInverseSqrtTwoPi * exponentials[count + lane];
        base_pdf = normal_pdf * B[i] /
                   (sv[i] * x[i] * x[i] * denom);
      }
      if constexpr (need_cdf) {
        base_cdf = clamp_probability(
            (1.0 - cdf[count + lane]) / denom);
      }
    }
    source_lane_store_fill_masked<Mask>(
        out, i,
        exact_source_finish_base_fill<Mask>(base_pdf, base_cdf, q[i]));
  }
}

template <std::uint8_t Mask>
inline void evaluate_lba_leaf_batch(
    PreparedSourceLeafBatch *input,
    SourceLaneFill *out) {
  const auto impossible = exact_source_impossible_fill(Mask);
  const auto *x = input->elapsed();
  const auto *v = input->parameter(0U);
  const auto *B = input->parameter(1U);
  const auto *A = input->parameter(2U);
  const auto *sv = input->parameter(3U);
  std::size_t wide_count = 0U;
  std::size_t point_begin = input->count;
  for (std::size_t i = 0U; i < input->count; ++i) {
    source_lane_store_fill_masked<Mask>(out, i, impossible);
    const bool valid = std::isfinite(x[i]) && x[i] > 0.0 &&
                       std::isfinite(v[i]) && std::isfinite(B[i]) &&
                       std::isfinite(A[i]) && std::isfinite(sv[i]) &&
                       sv[i] > 0.0;
    if (!valid) {
      continue;
    }
    if (A[i] > 1e-10) {
      const double zs = x[i] * sv[i];
      if (!(std::isfinite(zs) && zs > 0.0)) {
        continue;
      }
      input->positions[wide_count++] = i;
    } else {
      input->positions[--point_begin] = i;
    }
  }
  if (wide_count != 0U) {
    evaluate_lba_leaf_group<Mask, true>(
        input, input->positions.data(), wide_count, out);
  }
  const std::size_t point_count = input->count - point_begin;
  if (point_count != 0U) {
    evaluate_lba_leaf_group<Mask, false>(
        input, input->positions.data() + point_begin, point_count, out);
  }
}

template <std::uint8_t Mask>
inline void evaluate_exgauss_leaf_batch(
    PreparedSourceLeafBatch *input,
    SourceLaneFill *out) {
  constexpr bool need_pdf = (Mask & kLeafChannelPdf) != 0U;
  constexpr bool need_cdf =
      (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  constexpr std::size_t normal_count = need_cdf ? 4U : 3U;
  constexpr std::size_t exponent_count = normal_count + 2U;
  const auto impossible = exact_source_impossible_fill(Mask);
  const auto *x = input->elapsed();
  const auto *mu = input->parameter(0U);
  const auto *sigma = input->parameter(1U);
  const auto *tau = input->parameter(2U);
  std::size_t count = 0U;
  for (std::size_t i = 0U; i < input->count; ++i) {
    source_lane_store_fill_masked<Mask>(out, i, impossible);
    if (std::isfinite(x[i]) && x[i] > 0.0 &&
        std::isfinite(mu[i]) && std::isfinite(sigma[i]) &&
        sigma[i] > 0.0 && std::isfinite(tau[i]) && tau[i] > 0.0) {
      input->positions[count++] = i;
    }
  }
  if (count == 0U) {
    return;
  }

  input->ensure_work(count, 2U * exponent_count +
                                2U * normal_count);
  auto *exponent_arguments = input->work_plane(0U);
  auto *factors = input->work_plane(exponent_count);
  auto *exponentials =
      input->work_plane(exponent_count + normal_count);
  auto *normal_cdf = input->work_plane(
      2U * exponent_count + normal_count);
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = input->positions[lane];
    const double inv_tau = 1.0 / tau[i];
    const double sigma_sq = sigma[i] * sigma[i];
    const double tau_sq = tau[i] * tau[i];
    const double sigma_over_tau = sigma[i] * inv_tau;
    const double lower_z = -mu[i] / sigma[i];
    const double z = (x[i] - mu[i]) / sigma[i];
    exponent_arguments[lane] = lower_z;
    exponent_arguments[count + lane] = lower_z - sigma_over_tau;
    exponent_arguments[2U * count + lane] = z - sigma_over_tau;
    if constexpr (need_cdf) {
      exponent_arguments[3U * count + lane] = z;
    }
    exponent_arguments[normal_count * count + lane] =
        sigma_sq / (2.0 * tau_sq) - (0.0 - mu[i]) * inv_tau;
    exponent_arguments[(normal_count + 1U) * count + lane] =
        sigma_sq / (2.0 * tau_sq) - (x[i] - mu[i]) * inv_tau;
  }
  prepare_normal_lanes<true>(
      exponent_arguments, normal_count * count,
      factors, exponent_arguments);
  exp_lanes(
      exponent_arguments, exponentials, exponent_count * count);
  finish_normal_cdf_lanes(
      normal_count * count, factors, exponentials, normal_cdf);

  const auto *q = input->q();
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = input->positions[lane];
    const double lower_cdf = clamp_probability(
        normal_cdf[lane] -
        exponentials[normal_count * count + lane] *
            normal_cdf[count + lane]);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      continue;
    }
    const double x_tail = normal_cdf[2U * count + lane];
    const double x_exponential =
        exponentials[(normal_count + 1U) * count + lane];
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (need_pdf) {
      base_pdf = (1.0 / tau[i]) * x_exponential * x_tail /
                 lower_survival;
    }
    if constexpr (need_cdf) {
      const double raw_cdf = clamp_probability(
          normal_cdf[3U * count + lane] - x_exponential * x_tail);
      base_cdf = clamp_probability(
          (raw_cdf - lower_cdf) / lower_survival);
    }
    source_lane_store_fill_masked<Mask>(
        out, i,
        exact_source_finish_base_fill<Mask>(base_pdf, base_cdf, q[i]));
  }
}

template <std::uint8_t Mask>
inline void evaluate_rdm_degenerate_group(
    PreparedSourceLeafBatch *input,
    const std::size_t *positions,
    const std::size_t count,
    SourceLaneFill *out) {
  constexpr bool need_pdf = (Mask & kLeafChannelPdf) != 0U;
  constexpr bool need_cdf =
      (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  constexpr std::size_t normal_count = need_cdf ? 2U : 0U;
  constexpr std::size_t pdf_exponent_count = need_pdf ? 1U : 0U;
  constexpr std::size_t cdf_exponent_count = need_cdf ? 1U : 0U;
  constexpr std::size_t exponent_count =
      normal_count + pdf_exponent_count + cdf_exponent_count;
  constexpr std::size_t operand_count = need_pdf ? 2U : 0U;
  input->ensure_work(
      count, 2U * exponent_count + 2U * normal_count + operand_count);
  auto *exponent_arguments = input->work_plane(0U);
  auto *factors = input->work_plane(exponent_count);
  auto *exponentials =
      input->work_plane(exponent_count + normal_count);
  auto *normal_cdf = input->work_plane(
      2U * exponent_count + normal_count);
  auto *k_values = need_pdf
      ? input->work_plane(2U * exponent_count + 2U * normal_count)
      : nullptr;
  auto *sqrt_x = need_pdf ? k_values + count : nullptr;
  const auto *x = input->elapsed();
  const auto *v = input->parameter(0U);
  const auto *B = input->parameter(1U);
  const auto *s = input->parameter(3U);
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    const double inv_s = 1.0 / s[i];
    const double l = rdm_clamp_drift(v[i] * inv_s);
    const double k = B[i] * inv_s;
    if constexpr (need_pdf) {
      k_values[lane] = k;
      sqrt_x[lane] = std::sqrt(x[i]);
    }
    if constexpr (need_cdf) {
      const double lambda = k * k;
      const double mu = k / l;
      const double scale = std::sqrt(lambda / x[i]);
      const double time_ratio = x[i] / mu;
      exponent_arguments[lane] = -scale * (1.0 + time_ratio);
      exponent_arguments[count + lane] = scale * (1.0 - time_ratio);
      exponent_arguments[
          (normal_count + pdf_exponent_count) * count + lane] =
          2.0 * lambda / mu;
    }
    if constexpr (need_pdf) {
      const double delta = x[i] * l - k;
      exponent_arguments[normal_count * count + lane] =
          -0.5 * delta * delta / x[i];
    }
  }
  if constexpr (need_cdf) {
    prepare_normal_lanes<true>(
        exponent_arguments, normal_count * count,
        factors, exponent_arguments);
  }
  exp_lanes(exponent_arguments, exponentials, exponent_count * count);
  if constexpr (need_cdf) {
    finish_normal_cdf_lanes(
        normal_count * count, factors, exponentials, normal_cdf);
  }

  const auto *q = input->q();
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (need_pdf) {
      const double k = k_values[lane];
      base_pdf = safe_density(
          exponentials[normal_count * count + lane] *
          (std::fabs(k) * kInverseSqrtTwoPi) / x[i] /
          sqrt_x[lane]);
    }
    if constexpr (need_cdf) {
      const double cdf =
          normal_cdf[lane] *
              exponentials[
                  (normal_count + pdf_exponent_count) * count + lane] +
          (1.0 - normal_cdf[count + lane]);
      base_cdf = clamp_probability(cdf);
    }
    source_lane_store_fill_masked<Mask>(
        out, i,
        exact_source_finish_base_fill<Mask>(base_pdf, base_cdf, q[i]));
  }
}

template <std::uint8_t Mask>
inline void evaluate_rdm_regular_group(
    PreparedSourceLeafBatch *input,
    const std::size_t *positions,
    const std::size_t count,
    SourceLaneFill *out) {
  constexpr bool need_pdf = (Mask & kLeafChannelPdf) != 0U;
  constexpr bool need_cdf =
      (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  constexpr bool combined = need_pdf && need_cdf;
  constexpr std::size_t normal_count = combined ? 5U
      : (need_pdf ? 2U : 4U);
  constexpr std::size_t exponent_count = normal_count + 2U;
  constexpr std::size_t second_stage_planes = need_cdf ? 4U : 0U;
  constexpr std::size_t operand_count = 2U +
      (need_cdf ? 2U : 0U) + (need_pdf ? 1U : 0U);
  input->ensure_work(
      count, 2U * exponent_count + 2U * normal_count +
                 second_stage_planes + operand_count);
  auto *exponent_arguments = input->work_plane(0U);
  auto *factors = input->work_plane(exponent_count);
  auto *exponentials =
      input->work_plane(exponent_count + normal_count);
  auto *normal_cdf = input->work_plane(
      2U * exponent_count + normal_count);
  auto *outer_arguments = input->work_plane(
      2U * exponent_count + 2U * normal_count);
  auto *outer_exponentials = need_cdf
      ? input->work_plane(
            2U * exponent_count + 2U * normal_count + 2U)
      : nullptr;
  auto *l_values = input->work_plane(
      2U * exponent_count + 2U * normal_count + second_stage_planes);
  auto *a_values = l_values + count;
  auto *k_values = need_cdf ? a_values + count : nullptr;
  auto *sqrt_x = need_cdf ? k_values + count : nullptr;
  auto *inv_sqrt_x = need_pdf
      ? a_values + (need_cdf ? 3U : 1U) * count
      : nullptr;
  const auto *x = input->elapsed();
  const auto *v = input->parameter(0U);
  const auto *B = input->parameter(1U);
  const auto *A = input->parameter(2U);
  const auto *s = input->parameter(3U);
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    const double inv_s = 1.0 / s[i];
    const double l = rdm_clamp_drift(v[i] * inv_s);
    const double a = std::max(kRdmAEpsilon, 0.5 * A[i] * inv_s);
    const double k = B[i] * inv_s + a;
    const double sqt = std::sqrt(x[i]);
    const double inv_sqt = 1.0 / sqt;
    const double inverse_x = 1.0 / x[i];
    l_values[lane] = l;
    a_values[lane] = a;
    if constexpr (need_cdf) {
      k_values[lane] = k;
      sqrt_x[lane] = sqt;
    }
    if constexpr (need_pdf) {
      inv_sqrt_x[lane] = inv_sqt;
    }
    if constexpr (need_pdf) {
      exponent_arguments[lane] =
          (-k + a) * inv_sqt + sqt * l;
      exponent_arguments[count + lane] =
          (k + a) * inv_sqt - sqt * l;
    }
    if constexpr (need_cdf) {
      constexpr std::size_t z_offset = combined ? 2U : 0U;
      exponent_arguments[z_offset * count + lane] =
          -(k - a + x[i] * l) * inv_sqt;
      exponent_arguments[(z_offset + 1U) * count + lane] =
          -(k + a + x[i] * l) * inv_sqt;
      if constexpr (combined) {
        exponent_arguments[4U * count + lane] =
            (k - a) * inv_sqt - sqt * l;
      } else {
        exponent_arguments[2U * count + lane] =
            (k + a) * inv_sqt - sqt * l;
        exponent_arguments[3U * count + lane] =
            (k - a) * inv_sqt - sqt * l;
      }
    }
    if constexpr (need_pdf) {
      const double delta_a = a - k + x[i] * l;
      const double delta_b = a + k - x[i] * l;
      exponent_arguments[normal_count * count + lane] =
          -0.5 * delta_a * delta_a * inverse_x;
      exponent_arguments[(normal_count + 1U) * count + lane] =
          -0.5 * delta_b * delta_b * inverse_x;
    } else {
      const double delta_a = k - a - x[i] * l;
      const double delta_b = a + k - x[i] * l;
      exponent_arguments[normal_count * count + lane] =
          -0.5 * delta_a * delta_a * inverse_x;
      exponent_arguments[(normal_count + 1U) * count + lane] =
          -0.5 * delta_b * delta_b * inverse_x;
    }
  }
  prepare_normal_lanes<true>(
      exponent_arguments, normal_count * count,
      factors, exponent_arguments);
  exp_lanes(exponent_arguments, exponentials, exponent_count * count);
  finish_normal_cdf_lanes(
      normal_count * count, factors, exponentials, normal_cdf);

  if constexpr (need_cdf) {
    auto *log_cdf = outer_arguments;
    constexpr std::size_t log_cdf_offset = combined ? 2U : 0U;
    const auto *log_cdf_input = normal_cdf + log_cdf_offset * count;
    log_lanes(log_cdf_input, log_cdf, 2U * count);
    for (std::size_t lane = 0U; lane < count; ++lane) {
      const double log_a = log_cdf_input[lane] > 0.0
                               ? log_cdf[lane]
                               : -1e30;
      const double log_b = log_cdf_input[count + lane] > 0.0
                               ? log_cdf[count + lane]
                               : -1e30;
      outer_arguments[lane] =
          2.0 * l_values[lane] * (k_values[lane] - a_values[lane]) +
          log_a;
      outer_arguments[count + lane] =
          2.0 * l_values[lane] * (k_values[lane] + a_values[lane]) +
          log_b;
    }
    exp_lanes(
        outer_arguments, outer_exponentials, 2U * count);
  }

  const auto *q = input->q();
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    const double l = l_values[lane];
    const double a = a_values[lane];
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (need_pdf) {
      const double t1 = kInverseSqrtTwoPi *
                        (exponentials[normal_count * count + lane] -
                         exponentials[(normal_count + 1U) * count + lane]) *
                        inv_sqrt_x[lane];
      const double t2a = 2.0 * normal_cdf[lane] - 1.0;
      const double t2b = 2.0 * normal_cdf[count + lane] - 1.0;
      const double t2 = 0.5 * l * (t2a + t2b);
      base_pdf = safe_density((t1 + t2) / (2.0 * a));
    }
    if constexpr (need_cdf) {
      const double k = k_values[lane];
      constexpr std::size_t t4a_index = combined ? 1U : 2U;
      constexpr std::size_t t4b_index = combined ? 4U : 3U;
      const double t1 = sqrt_x[lane] * kInverseSqrtTwoPi *
                        (exponentials[normal_count * count + lane] -
                         exponentials[(normal_count + 1U) * count + lane]);
      const double t2 = a +
          (outer_exponentials[count + lane] - outer_exponentials[lane]) /
              (2.0 * l);
      const double t4a =
          2.0 * normal_cdf[t4a_index * count + lane] - 1.0;
      const double t4b =
          2.0 * normal_cdf[t4b_index * count + lane] - 1.0;
      const double t4 =
          0.5 * (x[i] * l - a - k + 0.5 / l) * t4a +
          0.5 * (k - a - x[i] * l - 0.5 / l) * t4b;
      base_cdf = clamp_probability(0.5 * (t4 + t2 + t1) / a);
    }
    source_lane_store_fill_masked<Mask>(
        out, i,
        exact_source_finish_base_fill<Mask>(base_pdf, base_cdf, q[i]));
  }
}

template <std::uint8_t Mask>
inline void evaluate_rdm_leaf_batch(
    PreparedSourceLeafBatch *input,
    SourceLaneFill *out) {
  const auto impossible = exact_source_impossible_fill(Mask);
  const auto *x = input->elapsed();
  const auto *v = input->parameter(0U);
  const auto *B = input->parameter(1U);
  const auto *A = input->parameter(2U);
  const auto *s = input->parameter(3U);
  std::size_t regular_count = 0U;
  std::size_t degenerate_begin = input->count;
  for (std::size_t i = 0U; i < input->count; ++i) {
    source_lane_store_fill_masked<Mask>(out, i, impossible);
    if (!std::isfinite(x[i]) || x[i] <= 0.0 ||
        !std::isfinite(s[i]) || s[i] <= 0.0) {
      continue;
    }
    const double inv_s = 1.0 / s[i];
    const double l = rdm_clamp_drift(v[i] * inv_s);
    if (!std::isfinite(l) || l < 0.0) {
      continue;
    }
    if (A[i] < kRdmAEpsilon) {
      const double k = B[i] * inv_s;
      if (!std::isfinite(k) || k <= 0.0 || k > kRdmKMaximum) {
        continue;
      }
      input->positions[--degenerate_begin] = i;
    } else {
      if (!std::isfinite(B[i]) || !std::isfinite(A[i])) {
        continue;
      }
      const double a = std::max(kRdmAEpsilon, 0.5 * A[i] * inv_s);
      const double k = B[i] * inv_s + a;
      if (!std::isfinite(a) || !std::isfinite(k)) {
        continue;
      }
      input->positions[regular_count++] = i;
    }
  }
  if (regular_count != 0U) {
    evaluate_rdm_regular_group<Mask>(
        input, input->positions.data(), regular_count, out);
  }
  const std::size_t degenerate_count = input->count - degenerate_begin;
  if (degenerate_count != 0U) {
    evaluate_rdm_degenerate_group<Mask>(
        input, input->positions.data() + degenerate_begin,
        degenerate_count, out);
  }
}

template <leaf::DistKind Kind, std::uint8_t Mask>
// Keep one numerical kernel per distribution/mask rather than cloning it into
// each lane-mapping and elapsed-time caller.
[[gnu::noinline]] inline void evaluate_prepared_source_leaf_batch(
    PreparedSourceLeafBatch &input,
    SourceLaneFill *out) {
  if constexpr (Kind == leaf::DistKind::Lognormal) {
    evaluate_lognormal_leaf_batch<Mask>(&input, out);
    return;
  } else if constexpr (Kind == leaf::DistKind::Exgauss) {
    evaluate_exgauss_leaf_batch<Mask>(&input, out);
    return;
  } else if constexpr (Kind == leaf::DistKind::LBA) {
    evaluate_lba_leaf_batch<Mask>(&input, out);
    return;
  } else if constexpr (Kind == leaf::DistKind::RDM) {
    evaluate_rdm_leaf_batch<Mask>(&input, out);
    return;
  }
  const auto *elapsed = input.elapsed();
  const auto *q = input.q();
  const auto *p0 = input.parameter(0U);
  const auto *p1 = input.parameter(1U);
  if constexpr (Kind == leaf::DistKind::Gamma) {
    for (std::size_t i = 0U; i < input.count; ++i) {
      auto fill = exact_source_impossible_fill(Mask);
      const double x = elapsed[i];
      if (x > 0.0) {
        fill = exact_source_gamma_leaf_fill<Mask>(
            p0[i], p1[i], q[i], x);
      }
      source_lane_store_fill_masked<Mask>(out, i, fill);
    }
  }
}

template <leaf::DistKind Kind>
inline void evaluate_prepared_source_leaf_batch(
    PreparedSourceLeafBatch &input,
    SourceLaneFill *out) {
  switch (out->mask) {
  case kLeafChannelPdf:
    evaluate_prepared_source_leaf_batch<Kind, kLeafChannelPdf>(input, out);
    break;
  case kLeafChannelCdf:
    evaluate_prepared_source_leaf_batch<Kind, kLeafChannelCdf>(input, out);
    break;
  case kLeafChannelSurvival:
    evaluate_prepared_source_leaf_batch<Kind, kLeafChannelSurvival>(
        input, out);
    break;
  case kLeafChannelPdf | kLeafChannelCdf:
    evaluate_prepared_source_leaf_batch<
        Kind, kLeafChannelPdf | kLeafChannelCdf>(input, out);
    break;
  case kLeafChannelPdf | kLeafChannelSurvival:
    evaluate_prepared_source_leaf_batch<
        Kind, kLeafChannelPdf | kLeafChannelSurvival>(input, out);
    break;
  case kLeafChannelCdf | kLeafChannelSurvival:
    evaluate_prepared_source_leaf_batch<
        Kind, kLeafChannelCdf | kLeafChannelSurvival>(input, out);
    break;
  case kLeafChannelAll:
    evaluate_prepared_source_leaf_batch<Kind, kLeafChannelAll>(input, out);
    break;
  default:
    break;
  }
}

} // namespace detail
} // namespace accumulatr::eval
