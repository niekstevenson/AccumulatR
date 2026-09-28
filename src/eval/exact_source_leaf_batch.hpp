#pragma once

#include <Rcpp.h>

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

  double *q() noexcept {
    return values.data() + count;
  }

  double *parameter(const std::size_t slot) noexcept {
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
      factors[i] = prepare_normal_cdf(value);
    }
    exponent_arguments[i] = -0.5 * abs_z * abs_z;
  }
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
    z[i] = input->elapsed()[i] > 0.0 ? (z[i] - meanlog[i]) / sdlog[i] : 0.0;
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
    auto fill = exact_source_impossible_fill();
    const double x = input->elapsed()[i];
    if (x > 0.0) {
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
      const double lower = B[i] - x[i] * v[i];
      z[count + lane] = (lower + A[i]) / zs;
      z[2U * count + lane] = lower / zs;
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
      const double xx = B[i] - x[i] * v[i];
      const double cmz = xx + A[i];
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
  const auto impossible = exact_source_impossible_fill();
  const auto *x = input->elapsed();
  const auto *A = input->parameter(2U);
  std::size_t wide_count = 0U;
  std::size_t point_begin = input->count;
  for (std::size_t i = 0U; i < input->count; ++i) {
    if (!std::isfinite(x[i]) || x[i] <= 0.0) {
      source_lane_store_fill_masked<Mask>(out, i, impossible);
      continue;
    }
    if (A[i] > 1e-10) {
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
  const auto impossible = exact_source_impossible_fill();
  const auto *x = input->elapsed();
  const auto *mu = input->parameter(0U);
  const auto *sigma = input->parameter(1U);
  const auto *tau = input->parameter(2U);
  std::size_t count = 0U;
  for (std::size_t i = 0U; i < input->count; ++i) {
    if (std::isfinite(x[i]) && x[i] > 0.0) {
      input->positions[count++] = i;
    } else {
      source_lane_store_fill_masked<Mask>(out, i, impossible);
    }
  }
  if (count == 0U) {
    return;
  }

  input->ensure_work(count, 12U);
  auto *exponent_arguments = input->work_plane(0U);
  auto *factors = input->work_plane(4U);
  auto *tail_factors = input->work_plane(6U);
  auto *exponentials = input->work_plane(8U);
  std::size_t exponent_count = 2U * count;
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = input->positions[lane];
    const double ratio = sigma[i] / tau[i];
    const auto prepare_endpoint = [&](const std::size_t slot, const double z) {
      const auto position = slot * count + lane;
      const double w = z - ratio;
      exponent_arguments[position] = -0.5 * z * z;
      if (slot == 0U || need_cdf) {
        factors[position] = prepare_normal_cdf(-z);
      }
      const double factor = normal_tail_factor(std::fabs(w));
      tail_factors[position] = std::copysign(factor, w > 0.0 ? -1.0 : 1.0);
      // exp(r*(r/2-z)) * Phi(w) shares exp(-z*z/2) with Phi(-z).
      // For w > 0, normal symmetry also needs the unweighted exponential.
      if (w > 0.0) {
        exponent_arguments[exponent_count++] = ratio * (0.5 * ratio - z);
      }
    };
    prepare_endpoint(0U, -mu[i] / sigma[i]);
    prepare_endpoint(1U, (x[i] - mu[i]) / sigma[i]);
  }
  exp_lanes(exponent_arguments, exponentials, exponent_count);

  const auto *q = input->q();
  std::size_t extra_exponential = 2U * count;
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = input->positions[lane];
    const auto weighted_tail = [&](const std::size_t slot) {
      const auto position = slot * count + lane;
      double value = tail_factors[position] * exponentials[position];
      if (std::signbit(tail_factors[position])) {
        value += exponentials[extra_exponential++];
      }
      return value;
    };
    const double lower_tail = weighted_tail(0U);
    const double x_tail = weighted_tail(1U);
    const double lower_survival =
        finish_normal_cdf(factors[lane], exponentials[lane]) + lower_tail;
    if (!(lower_survival > 0.0) || !std::isfinite(lower_survival)) {
      source_lane_store_fill_masked<Mask>(out, i, impossible);
      continue;
    }
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (need_pdf) {
      base_pdf = (1.0 / tau[i]) * x_tail / lower_survival;
    }
    if constexpr (need_cdf) {
      const double survival =
          finish_normal_cdf(factors[count + lane], exponentials[count + lane]) +
          x_tail;
      base_cdf = clamp_probability(
          1.0 - survival / lower_survival);
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
  input->ensure_work(count, need_cdf ? 12U : 9U);
  auto *exponent_arguments = input->work_plane(0U);
  auto *factors = input->work_plane(2U);
  auto *exponentials = input->work_plane(4U);
  auto *l_values = input->work_plane(6U);
  auto *a_values = input->work_plane(7U);
  auto *sqrt_x = input->work_plane(8U);
  auto *k_values = need_cdf ? input->work_plane(9U) : nullptr;
  auto *tail_factors = need_cdf ? input->work_plane(10U) : nullptr;
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
    const double drift_time = sqt * l;
    const double z_a = (k - a) * inv_sqt - drift_time;
    const double z_b = (k + a) * inv_sqt - drift_time;
    exponent_arguments[lane] = -0.5 * z_a * z_a;
    exponent_arguments[count + lane] = -0.5 * z_b * z_b;
    factors[lane] = prepare_normal_cdf(z_a);
    factors[count + lane] = prepare_normal_cdf(z_b);
    l_values[lane] = l;
    a_values[lane] = a;
    sqrt_x[lane] = sqt;
    if constexpr (need_cdf) {
      k_values[lane] = k;
      // For w = (d + x*l)/sqrt(x) >= 0, with d = k +/- a:
      // exp(2*l*d) * Phi(-w) = tail_factor(w) * exp(-z*z/2),
      // where z = (d - x*l)/sqrt(x). Reuse the two Gaussians;
      // the weighted tails need neither logarithms nor extra exponentials.
      const double w_a = (k - a) * inv_sqt + drift_time;
      const double w_b = (k + a) * inv_sqt + drift_time;
      tail_factors[lane] =
          std::copysign(normal_tail_factor(std::fabs(w_a)), w_a);
      tail_factors[count + lane] =
          std::copysign(normal_tail_factor(std::fabs(w_b)), w_b);
    }
  }
  exp_lanes(exponent_arguments, exponentials, 2U * count);

  const auto *q = input->q();
  for (std::size_t lane = 0U; lane < count; ++lane) {
    const auto i = positions[lane];
    const double l = l_values[lane];
    const double a = a_values[lane];
    const double normal_a = finish_normal_cdf(factors[lane], exponentials[lane]);
    const double normal_b = finish_normal_cdf(
        factors[count + lane], exponentials[count + lane]);
    const double normal_difference = normal_b - normal_a;
    const double gaussian_difference = kInverseSqrtTwoPi *
        (exponentials[lane] - exponentials[count + lane]);
    double base_pdf = 0.0;
    double base_cdf = 0.0;
    if constexpr (need_pdf) {
      base_pdf = safe_density(
          (gaussian_difference / sqrt_x[lane] + l * normal_difference) /
          (2.0 * a));
    }
    if constexpr (need_cdf) {
      const double k = k_values[lane];
      double tail_a = tail_factors[lane] * exponentials[lane];
      double tail_b = tail_factors[count + lane] * exponentials[count + lane];
      // Normal symmetry covers negative w without taking a logarithm.
      if (std::signbit(tail_factors[lane])) {
        tail_a += std::exp(2.0 * l * (k - a));
      }
      if (std::signbit(tail_factors[count + lane])) {
        tail_b += std::exp(2.0 * l * (k + a));
      }
      const double t1 = sqrt_x[lane] * gaussian_difference;
      const double t2 = (tail_b - tail_a) / (2.0 * l);
      // Combine the normal terms before multiplying by time, so the far
      // tail does not subtract two large, nearly equal time coefficients.
      const double t4 = (x[i] * l - k + 0.5 / l) * normal_difference;
      base_cdf = clamp_probability(
          1.0 - 0.5 * (normal_a + normal_b) + (t1 + t2 + t4) / (2.0 * a));
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
  const auto impossible = exact_source_impossible_fill();
  const auto *x = input->elapsed();
  const auto *A = input->parameter(2U);
  std::size_t regular_count = 0U;
  std::size_t degenerate_begin = input->count;
  for (std::size_t i = 0U; i < input->count; ++i) {
    if (!std::isfinite(x[i]) || x[i] <= 0.0) {
      source_lane_store_fill_masked<Mask>(out, i, impossible);
      continue;
    }
    if (A[i] < kRdmAEpsilon) {
      input->positions[--degenerate_begin] = i;
    } else {
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
      auto fill = exact_source_impossible_fill();
      const double x = elapsed[i];
      if (x > 0.0) {
        const double scale = 1.0 / p1[i];
        fill = exact_source_finish_base_fill<Mask>(
            (Mask & kLeafChannelPdf) != 0U
                ? R::dgamma(x, p0[i], scale, 0) : 0.0,
            (Mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U
                ? R::pgamma(x, p0[i], scale, 1, 0) : 0.0,
            q[i]);
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
