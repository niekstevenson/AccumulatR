#pragma once

#include <Rcpp.h>

#include <algorithm>
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

  explicit ParamView(SEXP paramsSEXP)
      : base(REAL(paramsSEXP)), nrow(Rf_nrows(paramsSEXP)) {}

  inline double q(const int row) const {
    return base[row];
  }

  inline double t0(const int row) const {
    return base[nrow + row];
  }

  inline double p(const int row, const int slot) const {
    return base[(slot + 2) * nrow + row];
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
};

struct DirectProbabilityQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  semantic::Index outcome_code{semantic::kInvalidIndex};
};

struct WeightedOutcomeTerm {
  semantic::Index outcome_code{semantic::kInvalidIndex};
  double weight{0.0};
};

struct DirectWeightedProbabilityQuery {
  semantic::Index trial_index{semantic::kInvalidIndex};
  semantic::Index variant_index{semantic::kInvalidIndex};
  std::vector<WeightedOutcomeTerm> terms;
};

struct OutcomePlan {
  semantic::SourceKind source_kind{semantic::SourceKind::Leaf};
  semantic::Index source_index{semantic::kInvalidIndex};
  std::vector<semantic::Index> support;
};

struct VariantPlan {
  runtime::LoweredDirectVariant lowered;
  std::vector<semantic::Index> outcome_index_by_code;
  std::vector<OutcomePlan> outcomes;
};

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

class TrialKernel {
public:
  TrialKernel(const VariantPlan &plan,
              const ParamView &params,
              const int first_param_row,
              const double rt)
      : plan_(plan),
        program_(plan.lowered.program),
        params_(params),
        first_param_row_(first_param_row),
        rt_(rt),
        leaf_channels_(static_cast<std::size_t>(program_.layout.n_leaves)),
        pool_channels_(static_cast<std::size_t>(program_.layout.n_pools)),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U) {}

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
  std::vector<leaf::EventChannels> leaf_channels_;
  std::vector<leaf::EventChannels> pool_channels_;
  std::vector<std::uint8_t> pool_ready_;

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
      const int row = first_param_row_ + i;
      const auto begin = program_.parameter_layout.leaf_param_offsets[static_cast<std::size_t>(i)];
      const auto end =
          program_.parameter_layout.leaf_param_offsets[static_cast<std::size_t>(i + 1)];
      double local_params[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      const int n_local = std::min<int>(end - begin, 8);
      for (int j = 0; j < n_local; ++j) {
        local_params[j] = params_.p(row, j);
      }
      double onset_t0 = params_.t0(row);
      if (static_cast<semantic::OnsetKind>(
              program_.onset_kind[static_cast<std::size_t>(i)]) ==
          semantic::OnsetKind::Absolute) {
        onset_t0 += program_.onset_abs_value[static_cast<std::size_t>(i)];
      }
      leaf_channels_[static_cast<std::size_t>(i)] = standard_leaf_channels(
          program_.leaf_dist_kind[static_cast<std::size_t>(i)],
          local_params,
          end - begin,
          params_.q(row),
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

inline Rcpp::NumericVector evaluate_direct_loglik_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectLoglikQuery> &queries,
    const double min_ll) {
  ParamView params(paramsSEXP);
  Rcpp::NumericVector out(queries.size(), min_ll);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    TrialKernel kernel(
        plan,
        params,
        static_cast<int>(layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row),
        query.rt);
    out[static_cast<R_xlen_t>(i)] =
        kernel.loglik_for_code(query.outcome_code, min_ll);
  }
  return out;
}

inline Rcpp::NumericVector evaluate_direct_probability_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectProbabilityQuery> &queries) {
  ParamView params(paramsSEXP);
  Rcpp::NumericVector out(queries.size(), 0.0);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    out[static_cast<R_xlen_t>(i)] = integrate_to_infinity(
        [&](const double rt) {
          TrialKernel kernel(
              plan,
              params,
              static_cast<int>(layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row),
              rt);
          return kernel.density_for_code(query.outcome_code);
        });
  }
  return out;
}

inline Rcpp::NumericVector evaluate_direct_weighted_probability_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectWeightedProbabilityQuery> &queries) {
  ParamView params(paramsSEXP);
  Rcpp::NumericVector out(queries.size(), 0.0);
  std::vector<semantic::Index> codes;
  std::vector<double> weights;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    codes.resize(query.terms.size());
    weights.resize(query.terms.size());
    for (std::size_t j = 0; j < query.terms.size(); ++j) {
      codes[j] = query.terms[j].outcome_code;
      weights[j] = query.terms[j].weight;
    }
    out[static_cast<R_xlen_t>(i)] = integrate_to_infinity(
        [&](const double rt) {
          TrialKernel kernel(
              plan,
              params,
              static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row),
              rt);
          return kernel.weighted_density_sum(
              codes.data(),
              weights.data(),
              codes.size());
        });
  }
  return out;
}

inline Rcpp::NumericVector evaluate_direct_total_probability_queries_cached(
    const std::vector<VariantPlan> &plans,
    const PreparedTrialLayout &layout,
    SEXP paramsSEXP,
    const std::vector<DirectProbabilityQuery> &queries) {
  ParamView params(paramsSEXP);
  Rcpp::NumericVector out(queries.size(), 0.0);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    const auto &plan = plans.at(static_cast<std::size_t>(query.variant_index));
    out[static_cast<R_xlen_t>(i)] = integrate_to_infinity(
        [&](const double rt) {
          TrialKernel kernel(
              plan,
              params,
              static_cast<int>(
                  layout.spans.at(static_cast<std::size_t>(query.trial_index)).start_row),
              rt);
          return kernel.total_density();
        });
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
    const double min_ll) {
  Rcpp::DataFrame data(dataSEXP);
  ParamView params(paramsSEXP);

  const auto table = read_prepared_data_view(data);
  const auto outcome = Rcpp::as<Rcpp::IntegerVector>(data["R"]);
  const auto rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);

  Rcpp::NumericVector loglik(layout.spans.size());
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
    TrialKernel kernel(plan, params, static_cast<int>(param_row), rt[row]);
    loglik[static_cast<R_xlen_t>(trial_index)] =
        kernel.loglik_for_code(outcome[row], min_ll);
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
  ParamView params(paramsSEXP);

  const auto queries =
      collapse_probability_queries(data, variant_index_by_component_code, layout);
  const auto blocks = build_variant_blocks(queries);

  Rcpp::NumericVector prob(queries.size());
  std::size_t param_row = 0;
  for (const auto &block : blocks) {
    const auto &plan = plans.at(static_cast<std::size_t>(block.variant_index));
    const auto leaf_count =
        static_cast<std::size_t>(plan.lowered.program.layout.n_leaves);
    for (int local = 0; local < block.row_count; ++local) {
      const auto query_idx = static_cast<std::size_t>(block.start_row + local);
      prob[static_cast<R_xlen_t>(query_idx)] = integrate_to_infinity(
          [&](const double rt) {
            TrialKernel kernel(
                plan,
                params,
                static_cast<int>(param_row),
                rt);
            return kernel.density_for_code(queries[query_idx].outcome_code);
          });
      param_row += leaf_count;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("probability") = prob,
      Rcpp::Named("n_trials") = static_cast<int>(queries.size()),
      Rcpp::Named("n_variants") = static_cast<int>(plans.size()),
      Rcpp::Named("n_model_leaves") = static_cast<int>(model.leaves.size()));
}

} // namespace detail

} // namespace accumulatr::eval
