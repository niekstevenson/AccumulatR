#pragma once

#include "accumulator.h"
#include "dist_vector.h"
#include "native_utils.h"
#include "quadrature_batch.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

double eval_pdf_single(const uuber::AccDistParams &cfg, double x);
double eval_cdf_single(const uuber::AccDistParams &cfg, double x);
struct LowerBoundTransform {
  double lower{0.0};
  double cdf_at_lower{0.0};
  double inv_tail_mass{1.0};
  bool active{false};
  bool valid{true};
};
LowerBoundTransform make_lower_bound_transform(const uuber::AccDistParams &cfg,
                                               double lower);
LowerBoundTransform
default_lower_bound_transform(const uuber::AccDistParams &cfg);
double eval_pdf_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform);
double eval_cdf_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform);
double eval_survival_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform);
void eval_pdf_vec_with_lower_bound(const uuber::AccDistParams &cfg,
                                   const LowerBoundTransform &transform,
                                   const double *x, std::size_t n,
                                   double *out);
void eval_cdf_vec_with_lower_bound(const uuber::AccDistParams &cfg,
                                   const LowerBoundTransform &transform,
                                   const double *x, std::size_t n,
                                   double *out);
double total_onset_with_t0(double onset, const uuber::AccDistParams &cfg);
double acc_density_from_cfg(double t, double onset, double q,
                            const uuber::AccDistParams &cfg);
double acc_survival_from_cfg(double t, double onset, double q,
                             const uuber::AccDistParams &cfg);
double acc_cdf_success_from_cfg(double t, double onset, double q,
                                const uuber::AccDistParams &cfg);

enum class EvalNeed : std::uint8_t {
  kDensity = 1 << 0,
  kSurvival = 1 << 1,
  kCDF = 1 << 2
};

constexpr EvalNeed kEvalAll =
    static_cast<EvalNeed>(static_cast<std::uint8_t>(EvalNeed::kDensity) |
                          static_cast<std::uint8_t>(EvalNeed::kSurvival) |
                          static_cast<std::uint8_t>(EvalNeed::kCDF));

inline EvalNeed operator|(EvalNeed lhs, EvalNeed rhs) {
  return static_cast<EvalNeed>(static_cast<std::uint8_t>(lhs) |
                               static_cast<std::uint8_t>(rhs));
}

inline EvalNeed operator&(EvalNeed lhs, EvalNeed rhs) {
  return static_cast<EvalNeed>(static_cast<std::uint8_t>(lhs) &
                               static_cast<std::uint8_t>(rhs));
}

inline bool needs_density(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kDensity) != 0;
}

inline bool needs_survival(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kSurvival) != 0;
}

inline bool needs_cdf(EvalNeed need) {
  return static_cast<std::uint8_t>(need & EvalNeed::kCDF) != 0;
}

struct NodeEvalResult {
  double density{std::numeric_limits<double>::quiet_NaN()};
  double survival{std::numeric_limits<double>::quiet_NaN()};
  double cdf{std::numeric_limits<double>::quiet_NaN()};
};

inline NodeEvalResult make_node_result(EvalNeed need, double density,
                                       double survival, double cdf) {
  NodeEvalResult out;
  out.density = needs_density(need) ? safe_density(density)
                                    : std::numeric_limits<double>::quiet_NaN();
  out.survival = needs_survival(need)
                     ? clamp_probability(survival)
                     : std::numeric_limits<double>::quiet_NaN();
  out.cdf = needs_cdf(need) ? clamp_probability(cdf)
                            : std::numeric_limits<double>::quiet_NaN();
  return out;
}

struct FusedIntegralResult {
  double density{0.0};
  double cdf{0.0};
  bool ok{false};
};

template <typename SourceDensityFn>
inline FusedIntegralResult integrate_fused_onset_terms(
    const uuber::AccDistParams &cfg, double t, double x_shift, double lower,
    double upper, bool need_density, bool need_cdf,
    const LowerBoundTransform *lower_bound_transform,
    SourceDensityFn &&source_density_fn) {
  FusedIntegralResult out;
  if (!(upper > lower) || (!need_density && !need_cdf)) {
    out.ok = true;
    return out;
  }
  constexpr int kSegments = 4;
  const double width = upper - lower;
  if (!std::isfinite(width) || width <= 0.0) {
    return out;
  }

  std::vector<double> nodes;
  std::vector<double> weights;
  nodes.reserve(static_cast<std::size_t>(kSegments * 15));
  weights.reserve(static_cast<std::size_t>(kSegments * 15));
  for (int seg = 0; seg < kSegments; ++seg) {
    const double seg_lo = lower + (width * static_cast<double>(seg) /
                                   static_cast<double>(kSegments));
    const double seg_hi = lower + (width * static_cast<double>(seg + 1) /
                                   static_cast<double>(kSegments));
    uuber::TimeBatch batch = uuber::build_time_batch(seg_lo, seg_hi);
    if (batch.nodes.empty() || batch.weights.size() != batch.nodes.size()) {
      return out;
    }
    nodes.insert(nodes.end(), batch.nodes.begin(), batch.nodes.end());
    weights.insert(weights.end(), batch.weights.begin(), batch.weights.end());
  }
  if (nodes.empty() || weights.size() != nodes.size()) {
    return out;
  }
  const std::size_t n = nodes.size();

  std::vector<double> source_density(n, 0.0);
  std::vector<double> shifted_times(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double u = nodes[i];
    const double fs = source_density_fn(u);
    source_density[i] = safe_density(fs);
    shifted_times[i] = t - u - x_shift;
  }

  std::vector<double> pdf_values;
  std::vector<double> cdf_values;
  if (need_density) {
    pdf_values.resize(n, 0.0);
    if (lower_bound_transform) {
      eval_pdf_vec_with_lower_bound(cfg, *lower_bound_transform,
                                    shifted_times.data(),
                                    shifted_times.size(), pdf_values.data());
    } else {
      uuber::eval_pdf_vec(cfg.code, cfg.p1, cfg.p2, cfg.p3, cfg.p4, cfg.p5,
                          cfg.p6, cfg.p7, cfg.p8, shifted_times.data(),
                          shifted_times.size(), pdf_values.data());
    }
  }
  if (need_cdf) {
    cdf_values.resize(n, 0.0);
    if (lower_bound_transform) {
      eval_cdf_vec_with_lower_bound(cfg, *lower_bound_transform,
                                    shifted_times.data(),
                                    shifted_times.size(), cdf_values.data());
    } else {
      uuber::eval_cdf_vec(cfg.code, cfg.p1, cfg.p2, cfg.p3, cfg.p4, cfg.p5,
                          cfg.p6, cfg.p7, cfg.p8, shifted_times.data(),
                          shifted_times.size(), cdf_values.data());
    }
  }

  double density_sum = 0.0;
  double cdf_sum = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double fs = source_density[i];
    if (!std::isfinite(fs) || fs <= 0.0) {
      continue;
    }
    const double w = weights[i];
    if (!std::isfinite(w) || w <= 0.0) {
      continue;
    }
    if (need_density) {
      const double fx = safe_density(pdf_values[i]);
      if (fx > 0.0) {
        density_sum += w * fs * fx;
      }
    }
    if (need_cdf) {
      const double Fx = clamp_probability(cdf_values[i]);
      if (Fx > 0.0) {
        cdf_sum += w * fs * Fx;
      }
    }
  }

  out.density = safe_density(density_sum);
  out.cdf = clamp_probability(cdf_sum);
  out.ok = true;
  return out;
}
