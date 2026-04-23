#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../leaf/channels.hpp"
#include "../runtime/direct_program.hpp"
#include "quadrature.hpp"
#include "trial_data.hpp"

namespace accumulatr::eval {

namespace detail {

inline double clamp_probability(double x) noexcept {
  if (!std::isfinite(x)) {
    return 0.0;
  }
  if (x <= 0.0) {
    return 0.0;
  }
  if (x >= 1.0) {
    return 1.0;
  }
  return x;
}

inline leaf::EventChannels impossible_channels() noexcept {
  return leaf::EventChannels::impossible();
}

inline double safe_density(double value) noexcept {
  return std::isfinite(value) && value > 0.0 ? value : 0.0;
}

inline double normal_cdf_fast(double z) noexcept {
  const double arg = -z * 0.7071067811865475244;
  return clamp_probability(0.5 * std::erfc(arg));
}

inline double exgauss_raw_pdf(double x, double mu, double sigma,
                              double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  const double tail = normal_cdf_fast(z - sigma_over_tau);
  return safe_density(inv_tau * std::exp(exponent) * tail);
}

inline double exgauss_raw_cdf(double x, double mu, double sigma,
                              double tau) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mu) || !std::isfinite(sigma) ||
      sigma <= 0.0 || !std::isfinite(tau) || tau <= 0.0) {
    return 0.0;
  }
  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double sigma_over_tau = sigma * inv_tau;
  const double z = (x - mu) / sigma;
  const double exponent = sigma_sq / (2.0 * tau_sq) - (x - mu) * inv_tau;
  const double tail = normal_cdf_fast(z - sigma_over_tau);
  const double base = normal_cdf_fast(z);
  const double exp_term = std::exp(exponent);
  return clamp_probability(base - exp_term * tail);
}

constexpr double kLogPi = 1.1447298858494001741434;

inline double lba_denom(double v, double sv) noexcept {
  if (!std::isfinite(sv) || sv <= 0.0) {
    return 0.0;
  }
  double denom = R::pnorm(v / sv, 0.0, 1.0, 1, 0);
  if (!std::isfinite(denom) || denom < 1e-10) {
    denom = 1e-10;
  }
  return denom;
}

inline double lba_pdf_fast(double x, double v, double B, double A, double sv) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) || !std::isfinite(sv) ||
      sv <= 0.0 || !std::isfinite(B) || !std::isfinite(A)) {
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

inline double lba_cdf_fast(double x, double v, double B, double A, double sv) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(v) || !std::isfinite(sv) ||
      sv <= 0.0 || !std::isfinite(B) || !std::isfinite(A)) {
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
    cdf = (1.0 + (zs * (R::dnorm(cz_max, 0.0, 1.0, 0) -
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

inline double rdm_pigt0(double x, double k, double l) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l)) {
    return 0.0;
  }
  if (std::fabs(l) < 1e-12) {
    const double z = k / std::sqrt(x);
    return clamp_probability(2.0 * (1.0 - R::pnorm(z, 0.0, 1.0, 1, 0)));
  }
  const double mu = k / l;
  const double lambda = k * k;
  const double p1 = 1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 + x / mu),
                                   0.0, 1.0, 1, 0);
  const double p2 = 1.0 - R::pnorm(std::sqrt(lambda / x) * (1.0 - x / mu),
                                   0.0, 1.0, 1, 0);
  const double part =
      std::exp(std::exp(std::log(2.0 * lambda) - std::log(mu)) +
               std::log(std::max(1e-300, p1)));
  return clamp_probability(part + p2);
}

inline double rdm_digt0(double x, double k, double l) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l)) {
    return 0.0;
  }
  const double lambda = k * k;
  double exponent = 0.0;
  if (l == 0.0) {
    exponent = -0.5 * lambda / x;
  } else {
    const double mu = k / l;
    exponent = -(lambda / (2.0 * x)) *
               ((x * x) / (mu * mu) - 2.0 * x / mu + 1.0);
  }
  return std::exp(exponent + 0.5 * std::log(lambda) -
                  0.5 * std::log(2.0 * x * x * x * M_PI));
}

inline double rdm_pigt(double x, double k, double l, double a,
                       double threshold = 1e-10) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l) ||
      !std::isfinite(a)) {
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
    const double t1 =
        std::exp(0.5 * (lgt - M_LN2 - kLogPi)) * (t1a - t1b);
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
                       double threshold = 1e-10) noexcept {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(k) || !std::isfinite(l) ||
      !std::isfinite(a)) {
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
        0.7071067811865475244 *
        (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
    const double t2a =
        2.0 * R::pnorm((-k + a) / sqt + sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2b =
        2.0 * R::pnorm((k + a) / sqt - sqt * l, 0.0, 1.0, 1, 0) - 1.0;
    const double t2 = std::exp(std::log(0.5) + std::log(l)) * (t2a + t2b);
    pdf = std::exp(std::log(std::max(1e-300, t1 + t2)) - M_LN2 - std::log(a));
  }
  return safe_density(pdf);
}

inline double rdm_pdf_fast(double x, double v, double B, double A, double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_digt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline double rdm_cdf_fast(double x, double v, double B, double A, double s) noexcept {
  if (!std::isfinite(s) || s <= 0.0) {
    return 0.0;
  }
  const double v_sc = v / s;
  const double B_sc = B / s;
  const double A_sc = A / s;
  return rdm_pigt(x, B_sc + 0.5 * A_sc, v_sc, 0.5 * A_sc);
}

inline leaf::EventChannels standard_leaf_channels(const std::uint8_t dist_kind,
                                                  const double *params,
                                                  const int n_params,
                                                  const double q,
                                                  const double t0,
                                                  const double t) {
  (void)n_params;
  if (!std::isfinite(t) || !std::isfinite(q) || q < 0.0 || q > 1.0 ||
      !std::isfinite(t0)) {
    return impossible_channels();
  }

  const double x = t - t0;
  if (!(x > 0.0)) {
    return impossible_channels();
  }

  double base_pdf = 0.0;
  double base_cdf = 0.0;
  double lower_cdf = 0.0;

  switch (static_cast<leaf::DistKind>(dist_kind)) {
  case leaf::DistKind::Lognormal: {
    const double m = params[0];
    const double s = params[1];
    if (!std::isfinite(m) || !std::isfinite(s) || s <= 0.0) {
      return impossible_channels();
    }
    base_pdf = R::dlnorm(x, m, s, 0);
    base_cdf = R::plnorm(x, m, s, 1, 0);
    break;
  }
  case leaf::DistKind::Gamma: {
    const double shape = params[0];
    const double rate = params[1];
    if (!std::isfinite(shape) || shape <= 0.0 || !std::isfinite(rate) ||
        rate <= 0.0) {
      return impossible_channels();
    }
    const double scale = 1.0 / rate;
    base_pdf = R::dgamma(x, shape, scale, 0);
    base_cdf = R::pgamma(x, shape, scale, 1, 0);
    break;
  }
  case leaf::DistKind::Exgauss: {
    const double mu = params[0];
    const double sigma = params[1];
    const double tau = params[2];
    if (!std::isfinite(mu) || !std::isfinite(sigma) || sigma <= 0.0 ||
        !std::isfinite(tau) || tau <= 0.0) {
      return impossible_channels();
    }
    base_pdf = exgauss_raw_pdf(x, mu, sigma, tau);
    base_cdf = exgauss_raw_cdf(x, mu, sigma, tau);
    lower_cdf = exgauss_raw_cdf(0.0, mu, sigma, tau);
    const double lower_survival = 1.0 - lower_cdf;
    if (!(lower_survival > 0.0)) {
      return impossible_channels();
    }
    base_pdf = safe_density(base_pdf / lower_survival);
    base_cdf = clamp_probability((base_cdf - lower_cdf) / lower_survival);
    break;
  }
  case leaf::DistKind::LBA: {
    base_pdf = lba_pdf_fast(x, params[0], params[1], params[2], params[3]);
    base_cdf = lba_cdf_fast(x, params[0], params[1], params[2], params[3]);
    break;
  }
  case leaf::DistKind::RDM: {
    base_pdf = rdm_pdf_fast(x, params[0], params[1], params[2], params[3]);
    base_cdf = rdm_cdf_fast(x, params[0], params[1], params[2], params[3]);
    break;
  }
  }

  base_pdf = std::isfinite(base_pdf) && base_pdf > 0.0 ? base_pdf : 0.0;
  base_cdf = clamp_probability(base_cdf);

  leaf::EventChannels out;
  const double start_prob = 1.0 - q;
  out.pdf = start_prob * base_pdf;
  out.cdf = clamp_probability(start_prob * base_cdf);
  out.survival = clamp_probability(1.0 - out.cdf);
  return out;
}

template <typename Fn>
double integrate_to_infinity(Fn &&density_fn) {
  return quadrature::integrate_tail_default(
      [&](const double t) {
        const double density = density_fn(t);
        return std::isfinite(density) && density > 0.0 ? density : 0.0;
      });
}

struct ParamView {
  const double *base;
  int nrow;
  const int *row_map{nullptr};

  explicit ParamView(SEXP paramsSEXP)
      : base(REAL(paramsSEXP)), nrow(Rf_nrows(paramsSEXP)) {}

  ParamView(SEXP paramsSEXP, const int *row_map_)
      : base(REAL(paramsSEXP)),
        nrow(Rf_nrows(paramsSEXP)),
        row_map(row_map_) {}

  inline int physical_row(const int row) const {
    return row_map == nullptr ? row : row_map[row];
  }

  inline double q(const int row) const {
    return base[physical_row(row)];
  }

  inline double t0(const int row) const {
    return base[nrow + physical_row(row)];
  }

  inline double p(const int row, const int slot) const {
    return base[(slot + 2) * nrow + physical_row(row)];
  }
};

struct ProbabilityQuery {
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
};

struct DirectLoglikQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
  double rt{NA_REAL};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

struct DirectProbabilityQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

struct DirectTriggerState {
  double weight{1.0};
  std::vector<std::uint8_t> shared_started;
};

struct WeightedOutcomeTerm {
  semantic::Index outcome_code{semantic::kInvalidIndex};
  double weight{0.0};
};

struct DirectMassQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  std::vector<WeightedOutcomeTerm> terms;
  bool include_total_finite{false};
  semantic::Index row_map_index{semantic::kInvalidIndex};
};

struct DirectMassResult {
  double weighted_mass{0.0};
  double total_finite_mass{0.0};
};

inline void append_weighted_term(std::vector<WeightedOutcomeTerm> *terms,
                                 const semantic::Index outcome_code,
                                 const double weight) {
  if (!(std::isfinite(weight) && weight > 0.0)) {
    return;
  }
  for (auto &term : *terms) {
    if (term.outcome_code == outcome_code) {
      term.weight += weight;
      return;
    }
  }
  terms->push_back(WeightedOutcomeTerm{outcome_code, weight});
}

inline void append_key_bytes(std::string *key,
                             const void *bytes,
                             const std::size_t n_bytes) {
  key->append(static_cast<const char *>(bytes), n_bytes);
}

template <typename T>
inline void append_key_bytes(std::string *key, const T &value) {
  append_key_bytes(key, &value, sizeof(T));
}

struct OutcomePlan {
  semantic::SourceKind source_kind{semantic::SourceKind::Leaf};
  semantic::Index source_index{semantic::kInvalidIndex};
  std::vector<semantic::Index> support;
};

struct VariantPlan {
  runtime::LoweredDirectVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<OutcomePlan> outcomes;
  std::vector<semantic::Index> shared_trigger_indices;
  bool leaf_outcome_partition{false};
};

inline std::string make_direct_mass_cache_key(const VariantPlan &plan,
                                              const ParamView &params,
                                              const int first_param_row,
                                              const DirectMassQuery &query,
                                              const semantic::Index *codes,
                                              const double *weights,
                                              const std::size_t n_terms) {
  std::string key;
  key.reserve(64U +
              static_cast<std::size_t>(plan.lowered.program.layout.n_leaves) * 96U +
              n_terms * 24U);
  append_key_bytes(&key, query.variant_index);
  append_key_bytes(&key, query.include_total_finite);
  append_key_bytes(&key, n_terms);
  for (std::size_t i = 0; i < n_terms; ++i) {
    append_key_bytes(&key, codes[i]);
    append_key_bytes(&key, weights[i]);
  }
  for (int leaf_index = 0; leaf_index < plan.lowered.program.layout.n_leaves;
       ++leaf_index) {
    const int row = first_param_row + leaf_index;
    append_key_bytes(&key, params.q(row));
    append_key_bytes(&key, params.t0(row));
    const auto &desc =
        plan.lowered.program.leaf_descriptors[static_cast<std::size_t>(leaf_index)];
    for (int j = 0; j < desc.param_count; ++j) {
      append_key_bytes(&key, params.p(row, j));
    }
  }
  return key;
}

inline bool variant_is_semantic_direct(const compile::CompiledVariant &variant) {
  return variant.semantic_backend == compile::BackendKind::Direct;
}

inline std::vector<semantic::Index> merge_support(
    std::vector<semantic::Index> merged,
    const std::vector<semantic::Index> &rhs,
    const std::string &what) {
  for (const auto idx : rhs) {
    if (std::find(merged.begin(), merged.end(), idx) != merged.end()) {
      throw std::runtime_error("direct evaluator requires disjoint " + what);
    }
    merged.push_back(idx);
  }
  std::sort(merged.begin(), merged.end());
  return merged;
}

inline std::vector<semantic::Index> collect_source_support(
    const runtime::DirectProgram &program,
    const semantic::SourceKind kind,
    const semantic::Index index,
    std::vector<std::vector<semantic::Index>> *pool_cache,
    std::vector<std::uint8_t> *pool_state) {
  if (kind == semantic::SourceKind::Leaf) {
    return {index};
  }
  if (kind != semantic::SourceKind::Pool) {
    throw std::runtime_error("direct evaluator encountered invalid source kind");
  }

  const auto pool_idx = static_cast<std::size_t>(index);
  if ((*pool_state)[pool_idx] == 2U) {
    return (*pool_cache)[pool_idx];
  }
  if ((*pool_state)[pool_idx] == 1U) {
    throw std::runtime_error("cyclic pool dependency in direct evaluator");
  }

  (*pool_state)[pool_idx] = 1U;
  std::vector<semantic::Index> merged;
  const auto begin = program.pool_member_offsets[pool_idx];
  const auto end = program.pool_member_offsets[pool_idx + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    const auto member_kind =
        static_cast<semantic::SourceKind>(program.pool_member_kind[static_cast<std::size_t>(i)]);
    const auto member_index = program.pool_member_indices[static_cast<std::size_t>(i)];
    merged = merge_support(
        std::move(merged),
        collect_source_support(program, member_kind, member_index, pool_cache, pool_state),
        "pool member supports");
  }
  (*pool_cache)[pool_idx] = merged;
  (*pool_state)[pool_idx] = 2U;
  return (*pool_cache)[pool_idx];
}

inline VariantPlan make_variant_plan(
    const runtime::LoweredDirectVariant &lowered,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_outcome_codes) {
  VariantPlan plan;
  plan.lowered = lowered;
  plan.outcome_index_by_code.assign(
      n_outcome_codes + 1U,
      semantic::kInvalidIndex);

  const auto &program = plan.lowered.program;
  std::vector<std::vector<semantic::Index>> pool_cache(
      static_cast<std::size_t>(program.layout.n_pools));
  std::vector<std::uint8_t> pool_state(
      static_cast<std::size_t>(program.layout.n_pools), 0U);

  plan.outcomes.reserve(static_cast<std::size_t>(program.layout.n_outcomes));
  for (int i = 0; i < program.layout.n_outcomes; ++i) {
    OutcomePlan outcome;
    outcome.source_kind =
        static_cast<semantic::SourceKind>(program.outcome_source_kind[static_cast<std::size_t>(i)]);
    outcome.source_index = program.outcome_source_index[static_cast<std::size_t>(i)];
    outcome.support = collect_source_support(
        program,
        outcome.source_kind,
        outcome.source_index,
        &pool_cache,
        &pool_state);
    plan.outcomes.push_back(std::move(outcome));
    const auto label =
        plan.lowered.outcome_labels[static_cast<std::size_t>(i)];
    const auto code_it = outcome_code_by_label.find(label);
    if (code_it == outcome_code_by_label.end()) {
      throw std::runtime_error(
          "direct evaluator found no prepared outcome code for '" + label + "'");
    }
    plan.outcome_index_by_code[static_cast<std::size_t>(code_it->second)] = i;
  }

  for (std::size_t i = 0; i < plan.outcomes.size(); ++i) {
    for (std::size_t j = i + 1U; j < plan.outcomes.size(); ++j) {
      std::vector<semantic::Index> overlap;
      std::set_intersection(
          plan.outcomes[i].support.begin(),
          plan.outcomes[i].support.end(),
          plan.outcomes[j].support.begin(),
          plan.outcomes[j].support.end(),
          std::back_inserter(overlap));
      if (!overlap.empty()) {
        throw std::runtime_error(
            "direct evaluator requires disjoint direct outcome supports");
      }
    }
  }

  for (int i = 0; i < program.layout.n_triggers; ++i) {
    const auto pos = static_cast<std::size_t>(i);
    if (static_cast<semantic::TriggerKind>(program.trigger_kind[pos]) ==
            semantic::TriggerKind::Shared &&
        program.trigger_member_offsets[pos + 1U] -
                program.trigger_member_offsets[pos] >
            1) {
      plan.shared_trigger_indices.push_back(i);
    }
  }

  if (static_cast<int>(plan.outcomes.size()) == program.layout.n_leaves) {
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(program.layout.n_leaves), 0U);
    bool all_leaf = true;
    for (const auto &outcome : plan.outcomes) {
      if (outcome.source_kind != semantic::SourceKind::Leaf ||
          outcome.source_index < 0 ||
          outcome.source_index >= program.layout.n_leaves) {
        all_leaf = false;
        break;
      }
      const auto pos = static_cast<std::size_t>(outcome.source_index);
      if (seen[pos] != 0U) {
        all_leaf = false;
        break;
      }
      seen[pos] = 1U;
    }
    if (all_leaf) {
      plan.leaf_outcome_partition = std::all_of(
          seen.begin(), seen.end(), [](const std::uint8_t value) {
            return value != 0U;
          });
    }
  }

  return plan;
}

inline std::vector<ProbabilityQuery> collapse_probability_queries(
    const Rcpp::DataFrame &data,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const PreparedTrialLayout &layout) {
  if (layout.spans.empty()) {
    return {};
  }

  const auto table = read_prepared_data_view(data);
  Rcpp::IntegerVector outcome = Rcpp::as<Rcpp::IntegerVector>(data["R"]);

  std::vector<ProbabilityQuery> out;
  out.reserve(layout.spans.size());
  for (const auto &span : layout.spans) {
    const auto row = static_cast<R_xlen_t>(span.start_row);
    const int component_code = table.component[row];
    out.push_back(ProbabilityQuery{
        variant_index_by_component_code[static_cast<std::size_t>(component_code)],
        outcome[row]});
  }
  return out;
}

inline void build_direct_plan_cache(
    const compile::CompiledModel &compiled,
    const std::unordered_map<std::string, semantic::Index> &component_code_by_id,
    const std::unordered_map<std::string, semantic::Index> &outcome_code_by_label,
    const std::size_t n_component_codes,
    const std::size_t n_outcome_codes,
    std::vector<semantic::Index> *variant_index_by_component_code,
    std::vector<VariantPlan> *plans) {
  variant_index_by_component_code->assign(
      n_component_codes + 1U,
      semantic::kInvalidIndex);
  plans->clear();
  plans->reserve(compiled.variants.size());
  for (const auto &variant : compiled.variants) {
    if (!variant_is_semantic_direct(variant)) {
      continue;
    }
    const auto plan_index = static_cast<semantic::Index>(plans->size());
    const auto component_it = component_code_by_id.find(variant.component_id);
    if (component_it == component_code_by_id.end()) {
      throw std::runtime_error(
          "direct evaluator found no prepared component code for '" +
          variant.component_id + "'");
    }
    (*variant_index_by_component_code)[static_cast<std::size_t>(component_it->second)] =
        plan_index;
    plans->push_back(
      make_variant_plan(
          runtime::lower_direct_variant(variant),
          outcome_code_by_label,
          n_outcome_codes));
  }
}

inline std::vector<DirectTriggerState> enumerate_direct_trigger_states(
    const VariantPlan &plan,
    const ParamView &params,
    const int first_param_row) {
  std::vector<DirectTriggerState> states(1);
  states.front().weight = 1.0;
  states.front().shared_started.assign(
      static_cast<std::size_t>(plan.lowered.program.layout.n_triggers), 2U);

  for (const auto trigger_index : plan.shared_trigger_indices) {
    const auto pos = static_cast<std::size_t>(trigger_index);
    double q = 0.0;
    if (plan.lowered.program.trigger_has_fixed_q[pos] != 0U) {
      q = plan.lowered.program.trigger_fixed_q[pos];
    } else {
      const auto member_begin = plan.lowered.program.trigger_member_offsets[pos];
      if (member_begin != plan.lowered.program.trigger_member_offsets[pos + 1U]) {
        q = params.q(first_param_row +
                     plan.lowered.program.trigger_member_indices[static_cast<std::size_t>(
                         member_begin)]);
      }
    }
    q = clamp_probability(q);

    std::vector<DirectTriggerState> next;
    next.reserve(states.size() * 2U);
    for (const auto &state : states) {
      if (q > 0.0) {
        auto fail = state;
        fail.weight *= q;
        fail.shared_started[pos] = 0U;
        next.push_back(std::move(fail));
      }
      if (q < 1.0) {
        auto start = state;
        start.weight *= (1.0 - q);
        start.shared_started[pos] = 1U;
        next.push_back(std::move(start));
      }
    }
    states.swap(next);
  }

  return states;
}

inline double direct_effective_leaf_q(const VariantPlan &plan,
                                      const ParamView &params,
                                      const int first_param_row,
                                      const DirectTriggerState &trigger_state,
                                      const semantic::Index leaf_index) {
  const auto pos = static_cast<std::size_t>(leaf_index);
  const auto trigger_index = plan.lowered.program.leaf_trigger_index[pos];
  if (trigger_index != semantic::kInvalidIndex &&
      static_cast<semantic::TriggerKind>(
          plan.lowered.program.trigger_kind[static_cast<std::size_t>(trigger_index)]) ==
          semantic::TriggerKind::Shared &&
      trigger_state.shared_started[static_cast<std::size_t>(trigger_index)] <= 1U) {
    return trigger_state.shared_started[static_cast<std::size_t>(trigger_index)] == 1U
               ? 0.0
               : 1.0;
  }
  return params.q(first_param_row + leaf_index);
}

inline bool direct_source_possible_under_state(
    const runtime::DirectProgram &program,
    const VariantPlan &plan,
    const ParamView &params,
    const int first_param_row,
    const DirectTriggerState &trigger_state,
    const semantic::SourceKind kind,
    const semantic::Index index,
    std::vector<std::int8_t> *pool_cache) {
  if (kind == semantic::SourceKind::Leaf) {
    return direct_effective_leaf_q(
               plan, params, first_param_row, trigger_state, index) < 1.0 - 1e-12;
  }
  const auto pos = static_cast<std::size_t>(index);
  if ((*pool_cache)[pos] >= 0) {
    return (*pool_cache)[pos] != 0;
  }
  int live_members = 0;
  const auto begin = program.pool_member_offsets[pos];
  const auto end = program.pool_member_offsets[pos + 1U];
  for (semantic::Index i = begin; i < end; ++i) {
    if (direct_source_possible_under_state(
            program,
            plan,
            params,
            first_param_row,
            trigger_state,
            static_cast<semantic::SourceKind>(
                program.pool_member_kind[static_cast<std::size_t>(i)]),
            program.pool_member_indices[static_cast<std::size_t>(i)],
            pool_cache)) {
      ++live_members;
      if (live_members >= program.pool_k[pos]) {
        (*pool_cache)[pos] = 1;
        return true;
      }
    }
  }
  (*pool_cache)[pos] = 0;
  return false;
}

inline bool direct_outcome_possible_under_state(const VariantPlan &plan,
                                                const ParamView &params,
                                                const int first_param_row,
                                                const DirectTriggerState &trigger_state,
                                                const semantic::Index outcome_code) {
  const auto outcome_index =
      plan.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
  const auto &outcome = plan.outcomes[static_cast<std::size_t>(outcome_index)];
  std::vector<std::int8_t> pool_cache(
      static_cast<std::size_t>(plan.lowered.program.layout.n_pools), -1);
  return direct_source_possible_under_state(
      plan.lowered.program,
      plan,
      params,
      first_param_row,
      trigger_state,
      outcome.source_kind,
      outcome.source_index,
      &pool_cache);
}

inline bool direct_any_outcome_possible_under_state(const VariantPlan &plan,
                                                    const ParamView &params,
                                                    const int first_param_row,
                                                    const DirectTriggerState &trigger_state) {
  std::vector<std::int8_t> pool_cache(
      static_cast<std::size_t>(plan.lowered.program.layout.n_pools), -1);
  for (const auto &outcome : plan.outcomes) {
    if (direct_source_possible_under_state(
            plan.lowered.program,
            plan,
            params,
            first_param_row,
            trigger_state,
            outcome.source_kind,
            outcome.source_index,
            &pool_cache)) {
      return true;
    }
  }
  return false;
}

inline bool direct_try_exact_total_finite_mass(const VariantPlan &plan,
                                               const ParamView &params,
                                               const int first_param_row,
                                               double *total_mass) {
  if (!plan.leaf_outcome_partition) {
    return false;
  }
  const auto states = enumerate_direct_trigger_states(plan, params, first_param_row);
  double total = 0.0;
  for (const auto &state : states) {
    if (!(state.weight > 0.0)) {
      continue;
    }
    double all_failed_prob = 1.0;
    for (semantic::Index leaf_index = 0;
         leaf_index < plan.lowered.program.layout.n_leaves;
         ++leaf_index) {
      const double q = clamp_probability(direct_effective_leaf_q(
          plan, params, first_param_row, state, leaf_index));
      all_failed_prob *= q;
    }
    total += state.weight * (1.0 - all_failed_prob);
  }
  *total_mass = clamp_probability(total);
  return true;
}

class TrialKernel {
public:
  struct LoadedLeafInput {
    std::array<double, 8> params{};
    double q{0.0};
    double t0{0.0};
  };

  TrialKernel(const VariantPlan &plan,
              const ParamView &params,
              const int first_param_row,
              const double rt,
              const std::vector<std::uint8_t> *shared_started = nullptr)
      : plan_(plan),
        program_(plan.lowered.program),
        params_(params),
        first_param_row_(first_param_row),
        rt_(rt),
        shared_started_(shared_started),
        leaf_inputs_(static_cast<std::size_t>(program_.layout.n_leaves)),
        leaf_channels_(static_cast<std::size_t>(program_.layout.n_leaves)),
        pool_channels_(static_cast<std::size_t>(program_.layout.n_pools)),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U) {
    for (int i = 0; i < program_.layout.n_leaves; ++i) {
      const auto pos = static_cast<std::size_t>(i);
      const auto &desc = program_.leaf_descriptors[pos];
      const int row = first_param_row_ + i;
      auto &loaded = leaf_inputs_[pos];
      const int n_local = std::min<int>(desc.param_count, 8);
      for (int j = 0; j < n_local; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = params_.p(row, j);
      }
      for (int j = n_local; j < 8; ++j) {
        loaded.params[static_cast<std::size_t>(j)] = 0.0;
      }
      loaded.q = params_.q(row);
      loaded.t0 = params_.t0(row);
    }
  }

  void set_rt(const double rt) {
    rt_ = rt;
    std::fill(pool_ready_.begin(), pool_ready_.end(), 0U);
  }

  double loglik_for_code(const semantic::Index outcome_code,
                         const double min_ll) {
    const auto outcome_index =
        plan_.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
    load_leaf_channels();
    const double prob =
        density_for_outcome_index(static_cast<std::size_t>(outcome_index));
    if (!std::isfinite(prob) || !(prob > 0.0)) {
      return min_ll;
    }
    return std::log(prob);
  }

  double density_for_code(const semantic::Index outcome_code) {
    const auto outcome_index =
        plan_.outcome_index_by_code[static_cast<std::size_t>(outcome_code)];
    load_leaf_channels();
    return density_for_outcome_index(static_cast<std::size_t>(outcome_index));
  }

  double total_density() {
    load_leaf_channels();
    double total = 0.0;
    for (std::size_t target_idx = 0; target_idx < plan_.outcomes.size();
         ++target_idx) {
      const double prob = density_for_outcome_index(target_idx);
      if (std::isfinite(prob) && prob > 0.0) {
        total += prob;
      }
    }
    return std::isfinite(total) && total > 0.0 ? total : 0.0;
  }

  double weighted_density_sum(const semantic::Index *outcome_codes,
                              const double *weights,
                              const std::size_t n_terms) {
    load_leaf_channels();
    double total = 0.0;
    for (std::size_t i = 0; i < n_terms; ++i) {
      const double weight = weights[i];
      if (!(weight > 0.0)) {
        continue;
      }
      const auto outcome_index = static_cast<std::size_t>(
          plan_.outcome_index_by_code[static_cast<std::size_t>(outcome_codes[i])]);
      total += weight * density_for_outcome_index(outcome_index);
    }
    return std::isfinite(total) && total > 0.0 ? total : 0.0;
  }

private:
  const VariantPlan &plan_;
  const runtime::DirectProgram &program_;
  const ParamView &params_;
  int first_param_row_;
  double rt_;
  const std::vector<std::uint8_t> *shared_started_;
  std::vector<LoadedLeafInput> leaf_inputs_;
  std::vector<leaf::EventChannels> leaf_channels_;
  std::vector<leaf::EventChannels> pool_channels_;
  std::vector<std::uint8_t> pool_ready_;

  double leaf_q(const semantic::Index leaf_index) const {
    if (shared_started_ == nullptr) {
      return leaf_inputs_[static_cast<std::size_t>(leaf_index)].q;
    }
    const auto trigger_index =
        program_.leaf_trigger_index[static_cast<std::size_t>(leaf_index)];
    if (trigger_index != semantic::kInvalidIndex &&
        static_cast<semantic::TriggerKind>(
            program_.trigger_kind[static_cast<std::size_t>(trigger_index)]) ==
            semantic::TriggerKind::Shared &&
        (*shared_started_)[static_cast<std::size_t>(trigger_index)] <= 1U) {
      return (*shared_started_)[static_cast<std::size_t>(trigger_index)] == 1U ? 0.0
                                                                               : 1.0;
    }
    return leaf_inputs_[static_cast<std::size_t>(leaf_index)].q;
  }

  double density_for_outcome_index(const std::size_t target_idx) {
    const auto target = source_channels(plan_.outcomes[target_idx].source_kind,
                                        plan_.outcomes[target_idx].source_index);
    double prob = target.pdf;
    for (std::size_t j = 0; j < plan_.outcomes.size(); ++j) {
      if (j == target_idx) {
        continue;
      }
      prob *= source_channels(plan_.outcomes[j].source_kind,
                              plan_.outcomes[j].source_index)
                  .survival;
      if (!(prob > 0.0)) {
        return 0.0;
      }
    }
    return std::isfinite(prob) && prob > 0.0 ? prob : 0.0;
  }

  void load_leaf_channels() {
    for (int i = 0; i < program_.layout.n_leaves; ++i) {
      const auto pos = static_cast<std::size_t>(i);
      const auto &desc = program_.leaf_descriptors[pos];
      const auto &loaded = leaf_inputs_[pos];
      const double onset_t0 =
          static_cast<semantic::OnsetKind>(desc.onset_kind) ==
                  semantic::OnsetKind::Absolute
              ? loaded.t0 + desc.onset_abs_value
              : loaded.t0;
      leaf_channels_[static_cast<std::size_t>(i)] = standard_leaf_channels(
          desc.dist_kind,
          loaded.params.data(),
          std::min(desc.param_count, 8),
          leaf_q(i),
          onset_t0,
          rt_);
    }
  }

  const leaf::EventChannels &source_channels(const semantic::SourceKind kind,
                                             const semantic::Index index) {
    if (kind == semantic::SourceKind::Leaf) {
      return leaf_channels_[static_cast<std::size_t>(index)];
    }
    return pool_channel(index);
  }

  const leaf::EventChannels &pool_channel(const semantic::Index pool_index) {
    const auto idx = static_cast<std::size_t>(pool_index);
    if (pool_ready_[idx] != 0U) {
      return pool_channels_[idx];
    }

    const auto begin = program_.pool_member_offsets[idx];
    const auto end = program_.pool_member_offsets[idx + 1U];
    const auto k = program_.pool_k[idx];
    const auto n_members = static_cast<int>(end - begin);

    std::vector<const leaf::EventChannels *> members;
    members.reserve(static_cast<std::size_t>(n_members));
    for (semantic::Index i = begin; i < end; ++i) {
      members.push_back(&source_channels(
          static_cast<semantic::SourceKind>(
              program_.pool_member_kind[static_cast<std::size_t>(i)]),
          program_.pool_member_indices[static_cast<std::size_t>(i)]));
    }

    std::vector<std::vector<double>> prefix(static_cast<std::size_t>(n_members + 1),
                                            std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    std::vector<std::vector<double>> suffix(static_cast<std::size_t>(n_members + 1),
                                            std::vector<double>(static_cast<std::size_t>(n_members + 1), 0.0));
    prefix[0][0] = 1.0;
    for (int i = 0; i < n_members; ++i) {
      for (int m = 0; m <= i; ++m) {
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)]->survival;
        prefix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m + 1)] +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)]->cdf;
      }
    }
    suffix[static_cast<std::size_t>(n_members)][0] = 1.0;
    for (int i = n_members - 1; i >= 0; --i) {
      const int count = n_members - i - 1;
      for (int m = 0; m <= count; ++m) {
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)]->survival;
        suffix[static_cast<std::size_t>(i)][static_cast<std::size_t>(m + 1)] +=
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(m)] *
            members[static_cast<std::size_t>(i)]->cdf;
      }
    }

    double surv = 0.0;
    for (int m = 0; m < k; ++m) {
      surv += prefix[static_cast<std::size_t>(n_members)][static_cast<std::size_t>(m)];
    }

    double pdf = 0.0;
    for (int i = 0; i < n_members; ++i) {
      double others_exact = 0.0;
      for (int left = 0; left < k; ++left) {
        const int right = (k - 1) - left;
        if (right < 0 || right > (n_members - i - 1)) {
          continue;
        }
        others_exact +=
            prefix[static_cast<std::size_t>(i)][static_cast<std::size_t>(left)] *
            suffix[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(right)];
      }
      pdf += members[static_cast<std::size_t>(i)]->pdf * others_exact;
    }

    pool_channels_[idx].pdf = std::isfinite(pdf) && pdf > 0.0 ? pdf : 0.0;
    pool_channels_[idx].survival = clamp_probability(surv);
    pool_channels_[idx].cdf =
        clamp_probability(1.0 - pool_channels_[idx].survival);
    pool_ready_[idx] = 1U;
    return pool_channels_[idx];
  }
};

class DirectMixtureKernel {
public:
  template <typename Predicate>
  DirectMixtureKernel(const VariantPlan &plan,
                      const ParamView &params,
                      const int first_param_row,
                      const double rt,
                      Predicate &&include_state)
      : states_(enumerate_direct_trigger_states(plan, params, first_param_row)) {
    kernels_.reserve(states_.size());
    weights_.reserve(states_.size());
    for (const auto &state : states_) {
      if (!(state.weight > 0.0) || !include_state(state)) {
        continue;
      }
      weights_.push_back(state.weight);
      kernels_.emplace_back(plan, params, first_param_row, rt, &state.shared_started);
    }
  }

  void set_rt(const double rt) {
    for (auto &kernel : kernels_) {
      kernel.set_rt(rt);
    }
  }

  double density_for_code(const semantic::Index outcome_code) {
    double total = 0.0;
    for (std::size_t i = 0; i < kernels_.size(); ++i) {
      total += weights_[i] * kernels_[i].density_for_code(outcome_code);
    }
    return std::isfinite(total) && total > 0.0 ? total : 0.0;
  }

  double total_density() {
    double total = 0.0;
    for (std::size_t i = 0; i < kernels_.size(); ++i) {
      total += weights_[i] * kernels_[i].total_density();
    }
    return std::isfinite(total) && total > 0.0 ? total : 0.0;
  }

  double weighted_density_sum(const semantic::Index *outcome_codes,
                              const double *weights,
                              const std::size_t n_terms) {
    double total = 0.0;
    for (std::size_t i = 0; i < kernels_.size(); ++i) {
      total += weights_[i] *
               kernels_[i].weighted_density_sum(
                   outcome_codes, weights, n_terms);
    }
    return std::isfinite(total) && total > 0.0 ? total : 0.0;
  }

  double loglik_for_code(const semantic::Index outcome_code,
                         const double min_ll) {
    const double density = density_for_code(outcome_code);
    return density > 0.0 ? std::log(density) : min_ll;
  }

private:
  std::vector<DirectTriggerState> states_;
  std::vector<double> weights_;
  std::vector<TrialKernel> kernels_;
};

inline Rcpp::NumericVector evaluate_direct_loglik_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectLoglikQuery> &queries,
    const double min_ll,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericVector out(queries.size(), min_ll);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    const int *row_map =
        row_maps != nullptr && query.row_map_index != semantic::kInvalidIndex
            ? row_maps->at(static_cast<std::size_t>(query.row_map_index)).data()
            : nullptr;
    ParamView params(paramsSEXP, row_map);
    const int first_param_row =
        row_map == nullptr
            ? static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row)
            : 0;
    if (plan.shared_trigger_indices.empty()) {
      TrialKernel kernel(plan, params, first_param_row, query.rt);
      out[static_cast<R_xlen_t>(i)] =
          kernel.loglik_for_code(query.outcome_code, min_ll);
    } else {
      DirectMixtureKernel kernel(
          plan,
          params,
          first_param_row,
          query.rt,
          [&](const DirectTriggerState &state) {
            return direct_outcome_possible_under_state(
                plan, params, first_param_row, state, query.outcome_code);
          });
      out[static_cast<R_xlen_t>(i)] =
          kernel.loglik_for_code(query.outcome_code, min_ll);
    }
  }
  return out;
}

inline Rcpp::NumericMatrix evaluate_direct_mass_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectMassQuery> &queries,
    const std::vector<std::vector<int>> *row_maps = nullptr) {
  Rcpp::NumericMatrix out(queries.size(), 2);
  std::unordered_map<std::string, DirectMassResult> memo;
  std::vector<semantic::Index> codes;
  std::vector<double> weights;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    const int *row_map =
        row_maps != nullptr && query.row_map_index != semantic::kInvalidIndex
            ? row_maps->at(static_cast<std::size_t>(query.row_map_index)).data()
            : nullptr;
    ParamView params(paramsSEXP, row_map);
    const bool need_weighted = !query.terms.empty();
    bool need_total = query.include_total_finite;
    double exact_total_mass = 0.0;
    const int first_param_row =
        row_map == nullptr
            ? static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row)
            : 0;
    if (need_total &&
        direct_try_exact_total_finite_mass(
            plan, params, first_param_row, &exact_total_mass)) {
      need_total = false;
      out(static_cast<R_xlen_t>(i), 1) = exact_total_mass;
    }
    if (!need_weighted && !need_total) {
      continue;
    }
    codes.resize(query.terms.size());
    weights.resize(query.terms.size());
    for (std::size_t j = 0; j < query.terms.size(); ++j) {
      codes[j] = query.terms[j].outcome_code;
      weights[j] = query.terms[j].weight;
    }
    if (codes.size() > 1U) {
      std::vector<std::size_t> order(codes.size());
      for (std::size_t j = 0; j < order.size(); ++j) {
        order[j] = j;
      }
      std::sort(order.begin(), order.end(), [&](const std::size_t lhs,
                                                const std::size_t rhs) {
        return codes[lhs] < codes[rhs];
      });
      std::vector<semantic::Index> codes_sorted(codes.size());
      std::vector<double> weights_sorted(weights.size());
      for (std::size_t j = 0; j < order.size(); ++j) {
        codes_sorted[j] = codes[order[j]];
        weights_sorted[j] = weights[order[j]];
      }
      codes.swap(codes_sorted);
      weights.swap(weights_sorted);
    }
    const auto cache_key = make_direct_mass_cache_key(
        plan,
        params,
        first_param_row,
        query,
        codes.data(),
        weights.data(),
        codes.size());
    const auto memo_it = memo.find(cache_key);
    if (memo_it != memo.end()) {
      out(static_cast<R_xlen_t>(i), 0) = memo_it->second.weighted_mass;
      if (query.include_total_finite) {
        out(static_cast<R_xlen_t>(i), 1) = memo_it->second.total_finite_mass;
      }
      continue;
    }
    DirectMassResult result;
    if (plan.shared_trigger_indices.empty()) {
      TrialKernel kernel(plan, params, first_param_row, 0.0);
      const auto &tail_rule = quadrature::canonical_tail_batch().nodes;
      for (std::size_t node = 0; node < tail_rule.nodes.size(); ++node) {
        kernel.set_rt(tail_rule.nodes[node]);
        if (need_weighted) {
          result.weighted_mass +=
              tail_rule.weights[node] * kernel.weighted_density_sum(
                                            codes.data(),
                                            weights.data(),
                                            codes.size());
        }
        if (need_total) {
          result.total_finite_mass +=
              tail_rule.weights[node] * kernel.total_density();
        }
      }
    } else {
      DirectMixtureKernel kernel(
          plan,
          params,
          first_param_row,
          0.0,
          [&](const DirectTriggerState &state) {
            if (need_total &&
                direct_any_outcome_possible_under_state(
                    plan, params, first_param_row, state)) {
              return true;
            }
            if (need_weighted) {
              for (const auto &term : query.terms) {
                if (term.weight > 0.0 &&
                    direct_outcome_possible_under_state(
                        plan,
                        params,
                        first_param_row,
                        state,
                        term.outcome_code)) {
                  return true;
                }
              }
            }
            return false;
          });
      const auto &tail_rule = quadrature::canonical_tail_batch().nodes;
      for (std::size_t node = 0; node < tail_rule.nodes.size(); ++node) {
        kernel.set_rt(tail_rule.nodes[node]);
        if (need_weighted) {
          result.weighted_mass +=
              tail_rule.weights[node] * kernel.weighted_density_sum(
                                            codes.data(),
                                            weights.data(),
                                            codes.size());
        }
        if (need_total) {
          result.total_finite_mass +=
              tail_rule.weights[node] * kernel.total_density();
        }
      }
    }
    result.weighted_mass =
        std::isfinite(result.weighted_mass) && result.weighted_mass > 0.0
            ? result.weighted_mass
            : 0.0;
    out(static_cast<R_xlen_t>(i), 0) = result.weighted_mass;
    if (query.include_total_finite) {
      const double total_finite = need_total ? result.total_finite_mass
                                             : exact_total_mass;
      result.total_finite_mass = clamp_probability(total_finite);
      out(static_cast<R_xlen_t>(i), 1) = result.total_finite_mass;
    }
    memo.emplace(std::move(cache_key), result);
  }
  return out;
}

inline Rcpp::List evaluate_direct_trials_cached(
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    const double min_ll,
    const std::vector<unsigned char> *selected = nullptr) {
  Rcpp::DataFrame data(dataSEXP);
  ParamView params(paramsSEXP);

  const auto table = read_prepared_data_view(data);
  const auto outcome = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);

  Rcpp::NumericVector loglik(layout.spans.size(), min_ll);
  std::vector<runtime::TrialBlock> blocks;
  std::size_t param_row = 0;
  runtime::TrialBlock current_block;
  bool have_block = false;
  for (std::size_t trial_index = 0; trial_index < layout.spans.size(); ++trial_index) {
    const auto &span = layout.spans[trial_index];
    const auto row = static_cast<R_xlen_t>(span.start_row);
    const int component_code = table.component[row];
    const auto variant_index =
        variant_index_by_component_code[static_cast<std::size_t>(component_code)];
    if (!have_block || current_block.variant_index != variant_index) {
      if (have_block) {
        blocks.push_back(current_block);
      }
      current_block.variant_index = variant_index;
      current_block.start_row = static_cast<int>(trial_index);
      current_block.row_count = 1;
      have_block = true;
    } else {
      ++current_block.row_count;
    }

    const auto &plan = plans.at(static_cast<std::size_t>(variant_index));
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    if (!trial_is_selected(selected, trial_index)) {
      param_row += leaf_count;
      continue;
    }
    if (plan.shared_trigger_indices.empty()) {
      TrialKernel kernel(plan, params, static_cast<int>(param_row), rt[row]);
      loglik[static_cast<R_xlen_t>(trial_index)] =
          kernel.loglik_for_code(outcome[row], min_ll);
    } else {
      DirectMixtureKernel kernel(
          plan,
          params,
          static_cast<int>(param_row),
          rt[row],
          [&](const DirectTriggerState &state) {
            return direct_outcome_possible_under_state(
                plan,
                params,
                static_cast<int>(param_row),
                state,
                outcome[row]);
          });
      loglik[static_cast<R_xlen_t>(trial_index)] =
          kernel.loglik_for_code(outcome[row], min_ll);
    }
    param_row += leaf_count;
  }
  if (have_block) {
    blocks.push_back(current_block);
  }
  Rcpp::List block_list(blocks.size());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    const auto &block = blocks[i];
    const auto &plan = plans.at(static_cast<std::size_t>(block.variant_index));
    block_list[i] = Rcpp::List::create(
        Rcpp::Named("variant_index") = block.variant_index,
        Rcpp::Named("component_id") = plan.lowered.component_id,
        Rcpp::Named("start_row") = block.start_row,
        Rcpp::Named("row_count") = block.row_count);
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < loglik.size(); ++i) {
    total_loglik += loglik[i];
  }

  return Rcpp::List::create(
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("total_loglik") = total_loglik,
      Rcpp::Named("blocks") = block_list,
      Rcpp::Named("n_trials") = static_cast<int>(layout.spans.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()),
      Rcpp::Named("n_model_leaves") = static_cast<int>(model.leaves.size()));
}

inline Rcpp::List evaluate_direct_outcome_probabilities_cached(
    const semantic::SemanticModel &model,
    const std::vector<semantic::Index> &variant_index_by_component_code,
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    SEXP dataSEXP) {
  Rcpp::DataFrame data(dataSEXP);

  const auto queries =
      collapse_probability_queries(data, variant_index_by_component_code, layout);

  Rcpp::NumericVector prob(queries.size());
  std::vector<DirectMassQuery> mass_queries;
  mass_queries.reserve(queries.size());
  for (const auto &query : queries) {
    DirectMassQuery mass_query;
    mass_query.trial_index = semantic::kInvalidIndex;
    mass_query.variant_index = query.variant_index;
    append_weighted_term(&mass_query.terms, query.outcome_code, 1.0);
    mass_queries.push_back(std::move(mass_query));
  }
  for (std::size_t i = 0; i < queries.size(); ++i) {
    mass_queries[i].trial_index = static_cast<semantic::Index>(i);
  }
  const auto mass = evaluate_direct_mass_queries_cached(
      plans, layout, paramsSEXP, mass_queries);
  for (R_xlen_t i = 0; i < prob.size(); ++i) {
    prob[i] = mass(i, 0);
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = prob,
      Rcpp::Named("n_trials") = static_cast<int>(queries.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()),
      Rcpp::Named("n_model_leaves") = static_cast<int>(model.leaves.size()));
}

} // namespace detail

} // namespace accumulatr::eval
