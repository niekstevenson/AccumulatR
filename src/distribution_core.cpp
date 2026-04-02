#include "distribution_core.h"
#include "distribution_kernels.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>

double eval_pdf_single(const uuber::AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return uuber::distkernels::lognormal_pdf(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    if (is_invalid_positive(cfg.p1) || is_invalid_positive(cfg.p2) ||
        Rcpp::NumericVector::is_na(x)) {
      return NA_REAL;
    }
    return uuber::distkernels::gamma_pdf(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    if (!std::isfinite(cfg.p2) || cfg.p2 <= 0.0 || !std::isfinite(cfg.p3) ||
        cfg.p3 <= 0.0 || Rcpp::NumericVector::is_na(x)) {
      return NA_REAL;
    }
    return uuber::distkernels::exgauss_pdf(x, cfg.p1, cfg.p2, cfg.p3);
  case uuber::ACC_DIST_LBA:
    return uuber::distkernels::lba_pdf(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  case uuber::ACC_DIST_RDM:
    return uuber::distkernels::rdm_pdf(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  default:
    return 0.0;
  }
}

double eval_cdf_single(const uuber::AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return uuber::distkernels::lognormal_cdf(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    if (is_invalid_positive(cfg.p1) || is_invalid_positive(cfg.p2) ||
        Rcpp::NumericVector::is_na(x)) {
      return NA_REAL;
    }
    return uuber::distkernels::gamma_cdf(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    if (!std::isfinite(cfg.p2) || cfg.p2 <= 0.0 || !std::isfinite(cfg.p3) ||
        cfg.p3 <= 0.0 || Rcpp::NumericVector::is_na(x)) {
      return NA_REAL;
    }
    return uuber::distkernels::exgauss_cdf(x, cfg.p1, cfg.p2, cfg.p3);
  case uuber::ACC_DIST_LBA:
    return uuber::distkernels::lba_cdf(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  case uuber::ACC_DIST_RDM:
    return uuber::distkernels::rdm_cdf(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  default:
    return 0.0;
  }
}

LowerBoundTransform make_lower_bound_transform(const uuber::AccDistParams &cfg,
                                               double lower) {
  LowerBoundTransform out;
  out.lower = lower;
  out.cdf_at_lower = 0.0;
  out.inv_tail_mass = 1.0;
  out.active = false;
  out.valid = true;

  if (!std::isfinite(lower)) {
    return out;
  }

  const double cdf_at_lower = eval_cdf_single(cfg, lower);
  if (Rcpp::NumericVector::is_na(cdf_at_lower) || !std::isfinite(cdf_at_lower)) {
    out.valid = false;
    return out;
  }
  const double clamped_cdf = clamp_probability(cdf_at_lower);
  const double tail_mass = 1.0 - clamped_cdf;
  if (!(tail_mass > 0.0)) {
    out.valid = false;
    return out;
  }

  out.cdf_at_lower = clamped_cdf;
  out.inv_tail_mass = 1.0 / tail_mass;
  out.active = (clamped_cdf > 0.0);
  return out;
}

LowerBoundTransform
default_lower_bound_transform(const uuber::AccDistParams &cfg) {
  if (cfg.code == uuber::ACC_DIST_EXGAUSS) {
    return make_lower_bound_transform(cfg, 0.0);
  }
  return LowerBoundTransform{};
}

double eval_pdf_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform) {
  if (!transform.valid) {
    return NA_REAL;
  }
  if (!transform.active) {
    return eval_pdf_single(cfg, x);
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x) || x < transform.lower) {
    return 0.0;
  }
  const double density = eval_pdf_single(cfg, x);
  if (Rcpp::NumericVector::is_na(density) || !std::isfinite(density)) {
    return NA_REAL;
  }
  return safe_density(density * transform.inv_tail_mass);
}

double eval_cdf_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform) {
  if (!transform.valid) {
    return NA_REAL;
  }
  if (!transform.active) {
    return eval_cdf_single(cfg, x);
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x)) {
    return x < 0.0 ? 0.0 : 1.0;
  }
  if (x < transform.lower) {
    return 0.0;
  }
  const double cdf = eval_cdf_single(cfg, x);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  return clamp_probability((clamp_probability(cdf) - transform.cdf_at_lower) *
                           transform.inv_tail_mass);
}

double eval_survival_single_with_lower_bound(
    const uuber::AccDistParams &cfg, double x,
    const LowerBoundTransform &transform) {
  const double cdf = eval_cdf_single_with_lower_bound(cfg, x, transform);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  return clamp_probability(1.0 - cdf);
}

void eval_pdf_vec_with_lower_bound(const uuber::AccDistParams &cfg,
                                   const LowerBoundTransform &transform,
                                   const double *x, std::size_t n,
                                   double *out) {
  uuber::eval_pdf_vec(cfg.code, cfg.p1, cfg.p2, cfg.p3, cfg.p4, cfg.p5,
                      cfg.p6, cfg.p7, cfg.p8, x, n, out);
  if (!transform.valid) {
    std::fill(out, out + n, NA_REAL);
    return;
  }
  if (!transform.active) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!std::isfinite(xi) || xi < transform.lower) {
      out[i] = 0.0;
      continue;
    }
    out[i] = safe_density(out[i] * transform.inv_tail_mass);
  }
}

void eval_cdf_vec_with_lower_bound(const uuber::AccDistParams &cfg,
                                   const LowerBoundTransform &transform,
                                   const double *x, std::size_t n,
                                   double *out) {
  uuber::eval_cdf_vec(cfg.code, cfg.p1, cfg.p2, cfg.p3, cfg.p4, cfg.p5,
                      cfg.p6, cfg.p7, cfg.p8, x, n, out);
  if (!transform.valid) {
    std::fill(out, out + n, NA_REAL);
    return;
  }
  if (!transform.active) {
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!std::isfinite(xi)) {
      out[i] = xi < 0.0 ? 0.0 : 1.0;
      continue;
    }
    if (xi < transform.lower) {
      out[i] = 0.0;
      continue;
    }
    out[i] = clamp_probability((clamp_probability(out[i]) -
                                transform.cdf_at_lower) *
                               transform.inv_tail_mass);
  }
}

double total_onset_with_t0(double onset, const uuber::AccDistParams &cfg) {
  return onset + cfg.t0;
}

double acc_density_from_cfg(double t, double onset, double q,
                            const uuber::AccDistParams &cfg) {
  const LowerBoundTransform lower_bound = default_lower_bound_transform(cfg);
  const double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t) || t < 0.0 || t < effective_onset) {
    return 0.0;
  }
  const double success_prob = 1.0 - q;
  if (success_prob <= 0.0) {
    return 0.0;
  }
  const double density =
      eval_pdf_single_with_lower_bound(cfg, t - effective_onset, lower_bound);
  if (Rcpp::NumericVector::is_na(density) || !std::isfinite(density)) {
    return NA_REAL;
  }
  return success_prob * density;
}

double acc_survival_from_cfg(double t, double onset, double q,
                             const uuber::AccDistParams &cfg) {
  const LowerBoundTransform lower_bound = default_lower_bound_transform(cfg);
  const double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t)) {
    return 0.0;
  }
  if (t < 0.0 || t < effective_onset) {
    return 1.0;
  }
  const double survival_underlying = eval_survival_single_with_lower_bound(
      cfg, t - effective_onset, lower_bound);
  if (Rcpp::NumericVector::is_na(survival_underlying) ||
      !std::isfinite(survival_underlying)) {
    return NA_REAL;
  }
  const double success_prob = 1.0 - q;
  const double result = q + success_prob * survival_underlying;
  return clamp(result, 0.0, 1.0);
}

double acc_cdf_success_from_cfg(double t, double onset, double q,
                                const uuber::AccDistParams &cfg) {
  const LowerBoundTransform lower_bound = default_lower_bound_transform(cfg);
  const double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t)) {
    return clamp_probability(1.0 - q);
  }
  if (t < 0.0 || t < effective_onset) {
    return 0.0;
  }
  const double cdf =
      eval_cdf_single_with_lower_bound(cfg, t - effective_onset, lower_bound);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  const double success_prob = 1.0 - q;
  return clamp(success_prob * cdf, 0.0, 1.0);
}
