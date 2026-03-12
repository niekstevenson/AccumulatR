#include "distribution_core.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>

namespace {

constexpr double kLogPi = 1.1447298858494001741434;

inline double lognormal_pdf_fast(double x, double meanlog, double sdlog) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) ||
      !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  const double norm = 0.3989422804014327 * inv_sigma;
  const double val = (norm / x) * std::exp(-0.5 * z * z);
  return std::isfinite(val) && val > 0.0 ? val : 0.0;
}

inline double lognormal_cdf_fast(double x, double meanlog, double sdlog) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) ||
      !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  const double arg = -z * 0.7071067811865475;
  const double val = 0.5 * std::erfc(arg);
  return clamp_probability(val);
}

inline double normal_cdf_fast(double z) {
  const double arg = -z * 0.7071067811865475;
  const double val = 0.5 * std::erfc(arg);
  return clamp(val, 0.0, 1.0);
}

inline double gamma_pdf_fast(double x, double shape, double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    return NA_REAL;
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x)) {
    return 0.0;
  }
  const double scale = 1.0 / rate;
  const double val = R::dgamma(x, shape, scale, 0);
  if (!std::isfinite(val) || val < 0.0) {
    return 0.0;
  }
  return val;
}

inline double gamma_cdf_fast(double x, double shape, double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    return NA_REAL;
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x)) {
    return x < 0.0 ? 0.0 : 1.0;
  }
  const double scale = 1.0 / rate;
  const double val = R::pgamma(x, shape, scale, 1, 0);
  if (!std::isfinite(val)) {
    return NA_REAL;
  }
  return clamp(val, 0.0, 1.0);
}

inline double exgauss_pdf_fast(double x, double mu, double sigma, double tau) {
  if (!std::isfinite(sigma) || sigma <= 0.0 || !std::isfinite(tau) ||
      tau <= 0.0) {
    return NA_REAL;
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x)) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = half_sigma_sq_over_tau_sq - (x - mu) * inv_tau;
  const double gaussian_tail = normal_cdf_fast(z - sigma_over_tau);
  const double val = inv_tau * std::exp(exponent) * gaussian_tail;
  if (!std::isfinite(val) || val < 0.0) {
    return 0.0;
  }
  return val;
}

inline double exgauss_cdf_fast(double x, double mu, double sigma, double tau) {
  if (!std::isfinite(sigma) || sigma <= 0.0 || !std::isfinite(tau) ||
      tau <= 0.0) {
    return NA_REAL;
  }
  if (Rcpp::NumericVector::is_na(x)) {
    return NA_REAL;
  }
  if (!std::isfinite(x)) {
    return x < 0.0 ? 0.0 : 1.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double base_cdf = normal_cdf_fast(z);
  const double tail = normal_cdf_fast(z - sigma_over_tau);
  const double exp_term =
      std::exp(half_sigma_sq_over_tau_sq - (x - mu) * inv_tau);
  const double val = base_cdf - exp_term * tail;
  if (!std::isfinite(val)) {
    return NA_REAL;
  }
  return clamp(val, 0.0, 1.0);
}

inline double lba_denom(double v, double sv) {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = R::pnorm(v / sv, 0.0, 1.0, 1, 0);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

inline double lba_pdf_fast(double x, double v, double sv, double B, double A) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) ||
      !std::isfinite(sv) || sv <= 0.0 || !std::isfinite(B) ||
      !std::isfinite(A)) {
    return 0.0;
  }
  const double denom = lba_denom(v, sv);
  if (!(denom > 0.0)) {
    return 0.0;
  }
  double pdf = 0.0;
  if (A > 1e-10) {
    const double zs = x * sv;
    if (!(zs > 0.0) || !std::isfinite(zs)) {
      return 0.0;
    }
    const double cmz = B - x * v;
    const double cz = cmz / zs;
    const double cz_max = (cmz - A) / zs;
    pdf = (v * (R::pnorm(cz, 0.0, 1.0, 1, 0) -
                R::pnorm(cz_max, 0.0, 1.0, 1, 0)) +
           sv * (R::dnorm(cz_max, 0.0, 1.0, 0) -
                 R::dnorm(cz, 0.0, 1.0, 0))) /
          (A * denom);
  } else {
    pdf = R::dnorm(B / x, v, sv, 0) * B / (x * x * denom);
  }
  return safe_density(pdf);
}

inline double lba_cdf_fast(double x, double v, double sv, double B, double A) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) ||
      !std::isfinite(sv) || sv <= 0.0 || !std::isfinite(B) ||
      !std::isfinite(A)) {
    return 0.0;
  }
  const double denom = lba_denom(v, sv);
  if (!(denom > 0.0)) {
    return 0.0;
  }
  double cdf = 0.0;
  if (A > 1e-10) {
    const double zs = x * sv;
    if (!(zs > 0.0) || !std::isfinite(zs)) {
      return 0.0;
    }
    const double cmz = B - x * v;
    const double xx = cmz - A;
    const double cz = cmz / zs;
    const double cz_max = xx / zs;
    cdf = (1.0 +
           (zs * (R::dnorm(cz_max, 0.0, 1.0, 0) -
                  R::dnorm(cz, 0.0, 1.0, 0)) +
            xx * R::pnorm(cz_max, 0.0, 1.0, 1, 0) -
            cmz * R::pnorm(cz, 0.0, 1.0, 1, 0)) /
               A) /
          denom;
  } else {
    cdf = R::pnorm(B / x, v, sv, 0, 0) / denom;
  }
  return clamp_probability(cdf);
}

inline double rdm_pigt0(double x, double k, double l) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  if (std::fabs(l) < 1e-12) {
    const double z = k / std::sqrt(x);
    return clamp_probability(2.0 * (1.0 - R::pnorm(z, 0.0, 1.0, 1, 0)));
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double p1 =
      1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 + x / mu), 0.0, 1.0, 1, 0);
  const double p2 =
      1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 - x / mu), 0.0, 1.0, 1, 0);
  const double part = std::exp(std::exp(std::log(2.0 * lambda) - std::log(mu)) +
                               std::log(std::max(1e-300, p1)));
  return clamp_probability(part + p2);
}

inline double rdm_digt0(double x, double k, double l) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l)) {
    return 0.0;
  }
  const double lambda = k * k;
  double e = 0.0;
  if (l == 0.0) {
    e = -0.5 * lambda / x;
  } else {
    const double mu = k / l;
    e = -(lambda / (2.0 * x)) * ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
  }
  return std::exp(e + 0.5 * std::log(lambda) -
                  0.5 * std::log(2.0 * x * x * x * M_PI));
}

inline double rdm_pigt(double x, double k, double l, double a,
                       double threshold = 1e-10) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return rdm_pigt0(x, k, l);
  }
  const double sqt = std::sqrt(x);
  const double lgt = std::log(x);
  double cdf = 0.0;
  if (l < threshold) {
    const double t5a = 2.0 * R::pnorm((k + a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
    const double t5b = 2.0 * R::pnorm((-k - a) / sqt, 0.0, 1.0, 1, 0) - 1.0;
    const double t6a =
        -0.5 * ((k + a) * (k + a) / x - M_LN2 - kLogPi + lgt) - std::log(a);
    const double t6b =
        -0.5 * ((k - a) * (k - a) / x - M_LN2 - kLogPi + lgt) - std::log(a);
    cdf = 1.0 + std::exp(t6a) - std::exp(t6b) +
          ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
  } else {
    const double t1a = std::exp(-0.5 * std::pow(k - a - x * l, 2.0) / x);
    const double t1b = std::exp(-0.5 * std::pow(a + k - x * l, 2.0) / x);
    const double t1 = std::exp(0.5 * (lgt - M_LN2 - kLogPi)) * (t1a - t1b);
    const double t2a =
        std::exp(2.0 * l * (k - a) +
                 R::pnorm(-(k - a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2b =
        std::exp(2.0 * l * (k + a) +
                 R::pnorm(-(k + a + x * l) / sqt, 0.0, 1.0, 1, 1));
    const double t2 = a + (t2b - t2a) / (2.0 * l);
    const double t4a =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t4b =
        2.0 * R::pnorm((k - a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t4 = 0.5 * (x * l - a - k + 0.5 / l) * t4a +
                      0.5 * (k - a - x * l - 0.5 / l) * t4b;
    cdf = 0.5 * (t4 + t2 + t1) / a;
  }
  if (!std::isfinite(cdf) || cdf < 0.0) {
    return 0.0;
  }
  return clamp_probability(cdf);
}

inline double rdm_digt(double x, double k, double l, double a,
                       double threshold = 1e-10) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) ||
      !std::isfinite(l) || !std::isfinite(a)) {
    return 0.0;
  }
  if (a < threshold) {
    return safe_density(rdm_digt0(x, k, l));
  }
  double pdf = 0.0;
  if (l < threshold) {
    const double term = std::exp(-(k - a) * (k - a) / (2.0 * x)) -
                        std::exp(-(k + a) * (k + a) / (2.0 * x));
    pdf = std::exp(-0.5 * (M_LN2 + kLogPi + std::log(x)) +
                   std::log(std::max(1e-300, term)) - M_LN2 - std::log(a));
  } else {
    const double sqt = std::sqrt(x);
    const double t1a = -std::pow(a - k + x * l, 2.0) / (2.0 * x);
    const double t1b = -std::pow(a + k - x * l, 2.0) / (2.0 * x);
    const double t1 =
        M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
    const double t2a =
        2.0 * R::pnorm((-k + a) / sqt + sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2b =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2 = std::exp(std::log(0.5) + std::log(l)) * (t2a + t2b);
    pdf = std::exp(std::log(std::max(1e-300, t1 + t2)) - M_LN2 - std::log(a));
  }
  return safe_density(pdf);
}

inline double rdm_pdf_fast(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_digt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double rdm_cdf_fast(double x, double v, double B, double A, double s) {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

} // namespace

double eval_pdf_single(const uuber::AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_pdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return gamma_pdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    return exgauss_pdf_fast(x, cfg.p1, cfg.p2, cfg.p3);
  case uuber::ACC_DIST_LBA:
    return lba_pdf_fast(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  case uuber::ACC_DIST_RDM:
    return rdm_pdf_fast(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  default:
    return 0.0;
  }
}

double eval_cdf_single(const uuber::AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_cdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return gamma_cdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    return exgauss_cdf_fast(x, cfg.p1, cfg.p2, cfg.p3);
  case uuber::ACC_DIST_LBA:
    return lba_cdf_fast(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
  case uuber::ACC_DIST_RDM:
    return rdm_cdf_fast(x, cfg.p1, cfg.p2, cfg.p3, cfg.p4);
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
