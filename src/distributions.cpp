// [[Rcpp::depends(Rcpp, BH)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#define UUBER_HAVE_BOOST_GK 1
#include <boost/math/quadrature/gauss_kronrod.hpp>
#else
#define UUBER_HAVE_BOOST_GK 0
#endif

#include "accumulator.h"
#include "context.h"
#include "integrate.h"
#include "prep_builder.h"
#include "proto.h"

using uuber::AccDistParams;
using uuber::ComponentMap;
using uuber::LabelRef;
using uuber::resolve_acc_params_entries;

Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog);
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog);
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector &x, double shape,
                                   double rate);
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector &x, double shape,
                                   double rate);
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau);
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau);

// [[Rcpp::export]]
SEXP native_context_build(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  return uuber::build_native_context(prep);
}

// [[Rcpp::export]]
Rcpp::RawVector native_prep_serialize_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  uuber::NativePrepProto proto = uuber::build_prep_proto(prep);
  std::vector<std::uint8_t> bytes = uuber::serialize_native_prep(proto);
  Rcpp::RawVector out(bytes.size());
  if (!bytes.empty()) {
    std::copy(bytes.begin(), bytes.end(), out.begin());
  }
  return out;
}

// [[Rcpp::export]]
SEXP native_context_from_proto_cpp(Rcpp::RawVector blob) {
  uuber::NativePrepProto proto = uuber::deserialize_native_prep(
      reinterpret_cast<const std::uint8_t *>(blob.begin()),
      static_cast<std::size_t>(blob.size()));
  std::unique_ptr<uuber::NativeContext> ctx =
      uuber::build_context_from_proto(proto);
  return Rcpp::XPtr<uuber::NativeContext>(ctx.release(), true);
}

// [[Rcpp::export]]
bool native_ctx_invalid(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP)
    return true;
  return R_ExternalPtrAddr(ptr) == nullptr;
}

// [[Rcpp::export]]
double boost_integrate_cpp(Rcpp::Function integrand, double lower, double upper,
                           double rel_tol, double abs_tol, int max_depth) {
  return uuber::integrate_boost(integrand, lower, upper, rel_tol, abs_tol,
                                max_depth);
}

double acc_density_cpp(double, double, double, const std::string &,
                       const Rcpp::List &);
double acc_survival_cpp(double, double, double, const std::string &,
                        const Rcpp::List &);
double acc_cdf_success_cpp(double, double, double, const std::string &,
                           const Rcpp::List &);
double pool_density_fast_cpp(const Rcpp::NumericVector &,
                             const Rcpp::NumericVector &, int);
double pool_survival_fast_cpp(const Rcpp::NumericVector &, int);
double pool_density_fast(const std::vector<double> &density,
                         const std::vector<double> &survival, int k);
double pool_survival_fast(const std::vector<double> &survival, int k);
Rcpp::List pool_density_combine_native(
    const Rcpp::NumericVector &dens_vec, const Rcpp::NumericVector &cdf_vec,
    const Rcpp::NumericVector &surv_vec,
    const Rcpp::NumericVector &cdf_success_vec,
    const Rcpp::NumericVector &surv_success_vec,
    const std::vector<int> &shared_index,
    const std::vector<uuber::PoolTemplateEntry> &templates,
    const std::vector<int> &base_forced_complete,
    const std::vector<int> &base_forced_survive);
std::vector<uuber::ProtoParamEntry> params_from_rcpp(const Rcpp::List &params,
                                                     const std::string &dist) {
  Rcpp::CharacterVector names = params.names();
  std::vector<uuber::ProtoParamEntry> out;
  out.reserve(params.size());
  for (R_xlen_t i = 0; i < params.size(); ++i) {
    if (params[i] == R_NilValue)
      continue;
    uuber::ProtoParamEntry entry;
    entry.name = Rcpp::as<std::string>(names[i]);
    SEXP val = params[i];
    if (Rf_isLogical(val)) {
      Rcpp::LogicalVector lv(val);
      if (lv.size() == 1) {
        entry.tag = uuber::ParamValueTag::LogicalScalar;
        entry.logical_scalar = lv[0];
      } else {
        entry.tag = uuber::ParamValueTag::LogicalVector;
        entry.logical_values.assign(lv.begin(), lv.end());
      }
    } else {
      Rcpp::NumericVector nv(val);
      if (nv.size() == 1) {
        entry.tag = uuber::ParamValueTag::NumericScalar;
        entry.numeric_scalar = nv[0];
      } else {
        entry.tag = uuber::ParamValueTag::NumericVector;
        entry.numeric_values.assign(nv.begin(), nv.end());
      }
    }
    out.push_back(std::move(entry));
  }
  return out;
}
Rcpp::List pool_build_templates_cpp(int n,
                                    const Rcpp::IntegerVector &member_ids,
                                    int pool_idx, int k);
Rcpp::List native_component_plan_impl(const Rcpp::List &structure,
                                      const Rcpp::DataFrame *trial_rows,
                                      double trial_id,
                                      const std::string *forced_component);
Rcpp::List pool_density_combine_cpp(const Rcpp::NumericVector &dens_vec,
                                    const Rcpp::NumericVector &cdf_vec,
                                    const Rcpp::NumericVector &surv_vec,
                                    const Rcpp::NumericVector &cdf_success_vec,
                                    const Rcpp::NumericVector &surv_success_vec,
                                    const Rcpp::IntegerVector &shared_index,
                                    const Rcpp::List &templates,
                                    const Rcpp::IntegerVector &forced_complete,
                                    const Rcpp::IntegerVector &forced_survive);

namespace {

std::vector<std::vector<int>>
generate_combinations(const std::vector<int> &elements, int choose);
std::vector<int> survivors_from_combo(const std::vector<int> &others,
                                      const std::vector<int> &combo);
uuber::PoolTemplateCacheEntry
build_pool_template_cache(int n, const std::vector<int> &member_ids,
                          int pool_idx, int k);

void sort_unique(std::vector<int> &vec);
std::vector<int> integer_vector_to_std(const Rcpp::IntegerVector &src,
                                       bool subtract_one = false);
inline bool is_invalid_positive(double value);

constexpr double kDefaultRelTol = 1e-5;
constexpr double kDefaultAbsTol = 1e-6;
constexpr int kDefaultMaxDepth = 12;
constexpr int kOutcomeIdxNA = -2;

template <typename T> inline T clamp(T val, T lo, T hi) {
  if (!std::isfinite(val)) {
    return val;
  }
  if (val < lo)
    return lo;
  if (val > hi)
    return hi;
  return val;
}

inline double clamp_unit(double val) {
  if (!std::isfinite(val))
    return 0.0;
  if (val < 0.0)
    return 0.0;
  if (val > 1.0)
    return 1.0;
  return val;
}

inline double clamp_probability(double value) {
  if (!std::isfinite(value))
    return 0.0;
  if (value < 0.0)
    return 0.0;
  if (value > 1.0)
    return 1.0;
  return value;
}

inline double safe_density(double value) {
  if (!std::isfinite(value) || value <= 0.0)
    return 0.0;
  return value;
}

inline double lognormal_pdf_fast(double x, double meanlog, double sdlog) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) ||
      !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  const double norm = 0.3989422804014327 * inv_sigma; // 1/sqrt(2*pi) * 1/sdlog
  double val = (norm / x) * std::exp(-0.5 * z * z);
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
  // 0.5 * erfc(-z / sqrt(2))
  const double arg = -z * 0.7071067811865475;
  double val = 0.5 * std::erfc(arg);
  if (!std::isfinite(val))
    return 0.0;
  if (val < 0.0)
    return 0.0;
  if (val > 1.0)
    return 1.0;
  return val;
}

inline double normal_cdf_fast(double z) {
  const double arg = -z * 0.7071067811865475;
  double val = 0.5 * std::erfc(arg);
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
  double val = R::dgamma(x, shape, scale, /*give_log =*/0);
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
  double val = R::pgamma(x, shape, scale, /*lower_tail =*/1, /*log_p =*/0);
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
  double val = inv_tau * std::exp(exponent) * gaussian_tail;
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
  double val = base_cdf - exp_term * tail;
  if (!std::isfinite(val)) {
    return NA_REAL;
  }
  return clamp(val, 0.0, 1.0);
}

inline double eval_pdf_single(const AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_pdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return gamma_pdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    return exgauss_pdf_fast(x, cfg.p1, cfg.p2, cfg.p3);
  default:
    return 0.0;
  }
}

inline double eval_cdf_single(const AccDistParams &cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_cdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return gamma_cdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_EXGAUSS:
    return exgauss_cdf_fast(x, cfg.p1, cfg.p2, cfg.p3);
  default:
    return 0.0;
  }
}

inline double total_onset_with_t0(double onset, const AccDistParams &cfg) {
  return onset + cfg.t0;
}

inline double acc_density_from_cfg(double t, double onset, double q,
                                   const AccDistParams &cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t) || t < 0.0)
    return 0.0;
  if (t < effective_onset)
    return 0.0;
  double success_prob = 1.0 - q;
  if (success_prob <= 0.0)
    return 0.0;
  double dens = eval_pdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(dens) || !std::isfinite(dens)) {
    return NA_REAL;
  }
  return success_prob * dens;
}

inline double acc_survival_from_cfg(double t, double onset, double q,
                                    const AccDistParams &cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t))
    return 0.0;
  if (t < 0.0)
    return 1.0;
  if (t < effective_onset)
    return 1.0;
  double cdf = eval_cdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  double surv_underlying = clamp(1.0 - cdf, 0.0, 1.0);
  double success_prob = 1.0 - q;
  double result = q + success_prob * surv_underlying;
  return clamp(result, 0.0, 1.0);
}

inline double acc_cdf_success_from_cfg(double t, double onset,
                                       const AccDistParams &cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t))
    return 1.0;
  if (t < 0.0)
    return 0.0;
  if (t < effective_onset)
    return 0.0;
  double cdf = eval_cdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  return clamp(cdf, 0.0, 1.0);
}

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

inline std::vector<int>
set_to_sorted_vector(const std::unordered_set<int> &items) {
  if (items.empty())
    return {};
  std::vector<int> out(items.begin(), items.end());
  sort_unique(out);
  return out;
}

inline std::vector<int> union_vectors(const std::vector<int> &a,
                                      const std::vector<int> &b) {
  if (a.empty()) {
    std::vector<int> out = b;
    sort_unique(out);
    return out;
  }
  std::vector<int> out = a;
  out.insert(out.end(), b.begin(), b.end());
  sort_unique(out);
  return out;
}

inline bool vector_contains(const std::vector<int> &values, int target) {
  if (values.empty())
    return false;
  return std::find(values.begin(), values.end(), target) != values.end();
}

struct ScenarioRecord {
  double weight{0.0};
  std::vector<int> forced_complete;
  std::vector<int> forced_survive;
};

inline ScenarioRecord
make_scenario_record(double weight, const std::vector<int> &forced_complete,
                     const std::vector<int> &forced_survive) {
  ScenarioRecord rec;
  rec.weight = weight;
  rec.forced_complete = forced_complete;
  rec.forced_survive = forced_survive;
  return rec;
}

std::vector<ScenarioRecord>
rcpp_scenarios_to_records(const Rcpp::List &scenarios) {
  std::vector<ScenarioRecord> out;
  out.reserve(scenarios.size());
  for (R_xlen_t i = 0; i < scenarios.size(); ++i) {
    if (scenarios[i] == R_NilValue)
      continue;
    Rcpp::List sc(scenarios[i]);
    double weight = Rcpp::as<double>(sc["weight"]);
    if (!std::isfinite(weight) || weight <= 0.0)
      continue;
    std::vector<int> fc = integer_vector_to_std(sc["forced_complete"], false);
    std::vector<int> fs = integer_vector_to_std(sc["forced_survive"], false);
    out.push_back(make_scenario_record(weight, fc, fs));
  }
  return out;
}

std::vector<int> forced_vec_from_sexp(SEXP vec) {
  std::vector<int> out;
  if (Rf_isNull(vec))
    return out;
  Rcpp::IntegerVector iv(vec);
  out.reserve(iv.size());
  for (int val : iv) {
    if (val == NA_INTEGER)
      continue;
    out.push_back(val);
  }
  return out;
}

std::unordered_set<int> make_forced_set(const std::vector<int> &ids) {
  return std::unordered_set<int>(ids.begin(), ids.end());
}

struct ForcedScopeFilter {
  const std::vector<int> *source_ids{nullptr};
  const ForcedScopeFilter *parent{nullptr};
};

inline bool scope_filter_allows_id(const ForcedScopeFilter *scope_filter,
                                   int id) {
  if (id == NA_INTEGER)
    return false;
  for (const ForcedScopeFilter *it = scope_filter; it != nullptr;
       it = it->parent) {
    if (!it->source_ids)
      continue;
    if (!std::binary_search(it->source_ids->begin(), it->source_ids->end(),
                            id)) {
      return false;
    }
  }
  return true;
}

inline bool forced_contains_scoped(const std::unordered_set<int> &forced, int id,
                                   const ForcedScopeFilter *scope_filter) {
  if (id == NA_INTEGER || forced.empty())
    return false;
  if (!scope_filter_allows_id(scope_filter, id))
    return false;
  return forced.count(id) > 0;
}

inline bool
forced_set_intersects_scope(const std::unordered_set<int> &forced,
                            const ForcedScopeFilter *scope_filter) {
  if (forced.empty())
    return false;
  if (!scope_filter)
    return true;
  for (int id : forced) {
    if (id != NA_INTEGER && scope_filter_allows_id(scope_filter, id)) {
      return true;
    }
  }
  return false;
}

struct TrialParamSet;

inline int label_id(const uuber::NativeContext &ctx, const std::string &label) {
  auto it = ctx.label_to_id.find(label);
  if (it == ctx.label_to_id.end())
    return NA_INTEGER;
  return it->second;
}

inline int component_index_of(const uuber::NativeContext &ctx,
                              const std::string &component) {
  if (component.empty() || component == "__default__")
    return -1;
  auto it = ctx.component_index.find(component);
  if (it == ctx.component_index.end())
    return -1;
  return it->second;
}

inline int outcome_index_of(const uuber::NativeContext &ctx,
                            const std::string &label) {
  if (label.empty())
    return -1;
  auto it = ctx.outcome_index.find(label);
  if (it == ctx.outcome_index.end() || it->second.empty())
    return -1;
  return it->second[0];
}

inline bool outcome_allows_component(const uuber::OutcomeContextInfo &info,
                                     const std::string &component) {
  if (info.allowed_components.empty())
    return true;
  if (component.empty())
    return false;
  return std::find(info.allowed_components.begin(),
                   info.allowed_components.end(),
                   component) != info.allowed_components.end();
}

inline bool competitor_node_allowed(const uuber::NativeContext &ctx,
                                    int node_id,
                                    const std::string &component) {
  if (component.empty())
    return true;
  bool matched = false;
  for (const auto &info : ctx.outcome_info) {
    if (info.node_id != node_id)
      continue;
    matched = true;
    if (outcome_allows_component(info, component))
      return true;
  }
  return !matched;
}

inline const std::vector<int> &
filter_competitor_ids(const uuber::NativeContext &ctx,
                      const std::vector<int> &competitor_ids,
                      const std::string &component,
                      std::vector<int> &scratch) {
  if (competitor_ids.empty() || component.empty())
    return competitor_ids;
  bool all_allowed = true;
  for (int node_id : competitor_ids) {
    if (!competitor_node_allowed(ctx, node_id, component)) {
      all_allowed = false;
      break;
    }
  }
  if (all_allowed)
    return competitor_ids;
  scratch.clear();
  scratch.reserve(competitor_ids.size());
  for (int node_id : competitor_ids) {
    if (competitor_node_allowed(ctx, node_id, component)) {
      scratch.push_back(node_id);
    }
  }
  return scratch;
}

inline int resolve_outcome_index(const uuber::NativeContext &ctx,
                                 const std::string &label,
                                 const std::string &component) {
  if (label.empty())
    return -1;
  auto it = ctx.outcome_index.find(label);
  if (it == ctx.outcome_index.end() || it->second.empty())
    return -1;

  const std::vector<int> &indices = it->second;
  if (indices.size() == 1)
    return indices[0];

  // Try to find exact match
  if (!component.empty()) {
    int fallback_idx = -1;
    for (int idx : indices) {
      const auto &info = ctx.outcome_info[static_cast<std::size_t>(idx)];
      if (outcome_allows_component(info, component)) {
        if (!info.allowed_components.empty())
          return idx;
        if (fallback_idx < 0)
          fallback_idx = idx;
      }
    }
    return fallback_idx;
  }

  // Fallback: return first
  return indices[0];
}

inline LabelRef resolve_label_ref_runtime(const uuber::NativeContext &ctx,
                                          const std::string &label) {
  LabelRef ref;
  if (label.empty())
    return ref;
  ref.label_id = label_id(ctx, label);
  auto acc_it = ctx.accumulator_index.find(label);
  if (acc_it != ctx.accumulator_index.end()) {
    ref.acc_idx = acc_it->second;
  }
  auto pool_it = ctx.pool_index.find(label);
  if (pool_it != ctx.pool_index.end()) {
    ref.pool_idx = pool_it->second;
  }
  ref.outcome_idx = outcome_index_of(ctx, label);
  return ref;
}

inline std::vector<int> ensure_source_ids(const uuber::NativeContext &ctx,
                                          const uuber::NativeNode &node) {
  std::vector<int> ids = node.source_ids;
  if (ids.empty() && !node.source.empty()) {
    int label_idx = node.source_ref.label_id;
    if (label_idx != NA_INTEGER && label_idx >= 0) {
      ids.push_back(label_idx);
    }
  }
  sort_unique(ids);
  return ids;
}

struct TrialAccumulatorParams {
  double onset{0.0};
  int onset_kind{uuber::ONSET_ABSOLUTE};
  int onset_source_acc_idx{-1};
  int onset_source_pool_idx{-1};
  double onset_lag{0.0};
  double q{0.0};
  double shared_q{std::numeric_limits<double>::quiet_NaN()};
  AccDistParams dist_cfg{};
  // Flags/components unused in flat fast path but kept for compatibility with
  // component checks.
  bool has_components{false};
  std::vector<std::string> components;
  std::vector<int> component_indices;
  std::string shared_trigger_id;
  bool has_override{false};
};

struct TrialParamSet {
  std::vector<TrialAccumulatorParams> acc_params;
  bool shared_trigger_layout_matches_context{true};
};

constexpr std::uint64_t kFNV64Offset = 1469598103934665603ULL;
constexpr std::uint64_t kFNV64Prime = 1099511628211ULL;
constexpr std::uint64_t kCanonicalQuietNaNBits = 0x7ff8000000000000ULL;

inline std::uint64_t canonical_double_bits(double value) {
  if (value == 0.0) {
    // Normalize -0.0 and +0.0.
    return 0ULL;
  }
  if (std::isnan(value)) {
    return kCanonicalQuietNaNBits;
  }
  std::uint64_t bits = 0ULL;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline void hash_append_bytes(std::uint64_t &hash, const void *data,
                              std::size_t n_bytes) {
  const auto *bytes = static_cast<const unsigned char *>(data);
  for (std::size_t i = 0; i < n_bytes; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= kFNV64Prime;
  }
}

inline void append_i32(std::string &payload, std::uint64_t &hash, int value) {
  const std::int32_t v = static_cast<std::int32_t>(value);
  const std::size_t start = payload.size();
  payload.resize(start + sizeof(v));
  std::memcpy(&payload[start], &v, sizeof(v));
  hash_append_bytes(hash, &payload[start], sizeof(v));
}

inline void append_u32(std::string &payload, std::uint64_t &hash,
                       std::uint32_t value) {
  const std::size_t start = payload.size();
  payload.resize(start + sizeof(value));
  std::memcpy(&payload[start], &value, sizeof(value));
  hash_append_bytes(hash, &payload[start], sizeof(value));
}

inline void append_u64(std::string &payload, std::uint64_t &hash,
                       std::uint64_t value) {
  const std::size_t start = payload.size();
  payload.resize(start + sizeof(value));
  std::memcpy(&payload[start], &value, sizeof(value));
  hash_append_bytes(hash, &payload[start], sizeof(value));
}

inline void append_string_payload(std::string &payload, std::uint64_t &hash,
                                  const std::string &value) {
  append_u32(payload, hash, static_cast<std::uint32_t>(value.size()));
  if (value.empty()) {
    return;
  }
  const std::size_t start = payload.size();
  payload.resize(start + value.size());
  std::memcpy(&payload[start], value.data(), value.size());
  hash_append_bytes(hash, value.data(), value.size());
}

inline void append_double_payload(std::string &payload, std::uint64_t &hash,
                                  double value) {
  append_u64(payload, hash, canonical_double_bits(value));
}

uuber::NAMapCacheKey
build_na_map_cache_key(const uuber::OutcomeContextInfo &info,
                       const TrialParamSet *params_ptr,
                       const std::vector<std::string> &comp_labels,
                       const std::vector<double> &comp_weights,
                       const std::string &outcome_label_context) {
  uuber::NAMapCacheKey key;
  if (!params_ptr) {
    return key;
  }

  std::size_t payload_bytes =
      64 + outcome_label_context.size() + (params_ptr->acc_params.size() * 96) +
      (comp_labels.size() * 24);
  for (const auto &label : comp_labels) {
    payload_bytes += label.size();
  }
  key.payload.reserve(payload_bytes);

  std::uint64_t hash = kFNV64Offset;
  append_i32(key.payload, hash, info.node_id);
  append_string_payload(key.payload, hash, outcome_label_context);

  append_u32(key.payload, hash,
             static_cast<std::uint32_t>(params_ptr->acc_params.size()));
  for (const auto &acc : params_ptr->acc_params) {
    append_double_payload(key.payload, hash, acc.q);
    append_double_payload(key.payload, hash, acc.shared_q);
    append_double_payload(key.payload, hash, acc.onset);
    append_i32(key.payload, hash, acc.onset_kind);
    append_i32(key.payload, hash, acc.onset_source_acc_idx);
    append_i32(key.payload, hash, acc.onset_source_pool_idx);
    append_double_payload(key.payload, hash, acc.onset_lag);
    append_double_payload(key.payload, hash, acc.dist_cfg.t0);
    append_double_payload(key.payload, hash, acc.dist_cfg.p1);
    append_double_payload(key.payload, hash, acc.dist_cfg.p2);
    append_double_payload(key.payload, hash, acc.dist_cfg.p3);
  }

  append_u32(key.payload, hash, static_cast<std::uint32_t>(comp_labels.size()));
  for (std::size_t i = 0; i < comp_labels.size(); ++i) {
    append_string_payload(key.payload, hash, comp_labels[i]);
    const double w = (i < comp_weights.size()) ? comp_weights[i] : 0.0;
    append_double_payload(key.payload, hash, w);
  }

  key.hash = static_cast<std::size_t>(hash);
  return key;
}

inline TrialAccumulatorParams
base_params(const uuber::NativeAccumulator &base) {
  TrialAccumulatorParams p;
  p.onset = base.onset;
  p.onset_kind = base.onset_kind;
  p.onset_source_acc_idx = base.onset_source_acc_idx;
  p.onset_source_pool_idx = base.onset_source_pool_idx;
  p.onset_lag = base.onset_lag;
  p.q = base.q;
  p.dist_cfg = base.dist_cfg;
  p.shared_q = base.q;
  p.shared_trigger_id = base.shared_trigger_id;
  p.has_components = !base.components.empty();
  p.components = base.components;
  p.component_indices = base.component_indices;
  p.has_override = false;
  return p;
}

// Build a full parameter set from base context when no trial-level overrides
// are supplied.
inline TrialParamSet build_base_paramset(const uuber::NativeContext &ctx) {
  TrialParamSet ps;
  ps.acc_params.reserve(ctx.accumulators.size());
  for (const auto &acc : ctx.accumulators) {
    ps.acc_params.push_back(base_params(acc));
  }
  return ps;
}

struct SharedTriggerInfo {
  const std::vector<int> *acc_indices{nullptr};
  std::vector<int> owned_acc_indices;
  double q{0.0};

  const std::vector<int> &indices() const {
    return acc_indices ? *acc_indices : owned_acc_indices;
  }
};

inline bool shared_trigger_layout_uses_context_groups(
    const uuber::NativeContext &ctx, const TrialParamSet &params) {
  if (ctx.shared_trigger_groups.empty())
    return false;
  if (!params.shared_trigger_layout_matches_context)
    return false;
  return params.acc_params.size() == ctx.accumulators.size();
}

// Collect shared trigger definitions for the current trial parameters.
inline std::vector<SharedTriggerInfo>
collect_shared_triggers(const uuber::NativeContext &ctx,
                        const TrialParamSet &params) {
  std::vector<SharedTriggerInfo> out;
  if (shared_trigger_layout_uses_context_groups(ctx, params)) {
    out.reserve(ctx.shared_trigger_groups.size());
    for (const auto &group : ctx.shared_trigger_groups) {
      int q_idx = group.q_acc_idx;
      double q_val = 0.0;
      bool q_set = false;
      if (q_idx >= 0 && q_idx < static_cast<int>(params.acc_params.size())) {
        q_val = params.acc_params[static_cast<std::size_t>(q_idx)].shared_q;
        q_set = true;
      } else if (!group.acc_indices.empty()) {
        q_idx = group.acc_indices.front();
        if (q_idx >= 0 && q_idx < static_cast<int>(params.acc_params.size())) {
          q_val = params.acc_params[static_cast<std::size_t>(q_idx)].shared_q;
          q_set = true;
        }
      }
      if (!q_set && q_idx >= 0 &&
          q_idx < static_cast<int>(ctx.accumulators.size())) {
        q_val = ctx.accumulators[static_cast<std::size_t>(q_idx)].q;
        q_set = true;
      }
      if (!q_set)
        continue;
      SharedTriggerInfo info;
      info.acc_indices = &group.acc_indices;
      info.q = clamp_probability(q_val);
      out.push_back(std::move(info));
    }
    return out;
  }

  // Fallback: derive shared triggers directly from trial parameters when
  // shared_trigger_id overrides alter trigger membership.
  std::unordered_map<std::string, SharedTriggerInfo> derived;
  for (std::size_t i = 0; i < params.acc_params.size(); ++i) {
    const auto &acc = params.acc_params[i];
    if (acc.shared_trigger_id.empty())
      continue;
    SharedTriggerInfo &info = derived[acc.shared_trigger_id];
    info.owned_acc_indices.push_back(static_cast<int>(i));
    info.q = acc.shared_q;
  }
  out.reserve(derived.size());
  for (auto &kv : derived) {
    kv.second.acc_indices = &kv.second.owned_acc_indices;
    kv.second.q = clamp_probability(kv.second.q);
    out.push_back(std::move(kv.second));
  }
  return out;
}

inline double
shared_trigger_mask_weight(const std::vector<SharedTriggerInfo> &triggers,
                           std::uint64_t mask);

struct SharedTriggerPlan {
  std::vector<SharedTriggerInfo> triggers;
  std::vector<double> mask_weights;
};

inline SharedTriggerPlan
build_shared_trigger_plan(const uuber::NativeContext &ctx,
                          const TrialParamSet *params_ptr) {
  SharedTriggerPlan plan;
  if (!params_ptr)
    return plan;
  plan.triggers = collect_shared_triggers(ctx, *params_ptr);
  if (plan.triggers.empty())
    return plan;
  if (plan.triggers.size() >= 63) {
    Rcpp::stop("Too many shared triggers for density evaluation");
  }
  const std::uint64_t n_states = 1ULL << plan.triggers.size();
  constexpr std::uint64_t kMaxPrecomputedWeights = 1ULL << 20;
  if (n_states <= kMaxPrecomputedWeights) {
    plan.mask_weights.resize(n_states);
    for (std::uint64_t mask = 0; mask < n_states; ++mask) {
      plan.mask_weights[mask] = shared_trigger_mask_weight(plan.triggers, mask);
    }
  }
  return plan;
}

inline double
shared_trigger_mask_weight(const std::vector<SharedTriggerInfo> &triggers,
                           std::uint64_t mask) {
  double w = 1.0;
  for (std::size_t i = 0; i < triggers.size(); ++i) {
    double q = clamp_probability(triggers[i].q);
    bool fail = ((mask >> i) & 1ULL) != 0ULL;
    w *= fail ? q : (1.0 - q);
    if (!std::isfinite(w) || w <= 0.0)
      return 0.0;
  }
  return w;
}

inline void apply_trigger_state_inplace(TrialParamSet &params,
                                        const SharedTriggerInfo &trigger,
                                        bool fail) {
  double q_val = fail ? 1.0 : 0.0;
  const std::vector<int> &indices = trigger.indices();
  for (int acc_idx : indices) {
    if (acc_idx < 0 || acc_idx >= static_cast<int>(params.acc_params.size()))
      continue;
    TrialAccumulatorParams &p =
        params.acc_params[static_cast<std::size_t>(acc_idx)];
    p.q = q_val;
    p.has_override = true;
  }
}

struct ComponentCacheEntry;
const std::unordered_set<int> kEmptyForcedSet{};

double component_keep_weight(const uuber::NativeContext &ctx, int component_idx,
                             int outcome_idx);

double native_outcome_probability_impl(
    SEXP ctxSEXP, int node_id, double upper,
    Rcpp::Nullable<Rcpp::String> component, SEXP forced_complete,
    SEXP forced_survive, const Rcpp::IntegerVector &competitor_ids,
    double rel_tol, double abs_tol, int max_depth,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, const std::string &outcome_label_context);

std::string normalize_label_string(const std::string &label);

double node_density_with_competitors_internal(
    const uuber::NativeContext &ctx, int node_id, double t,
    const std::string &component_label, int component_idx,
    const std::unordered_set<int> &forced_complete,
    const std::unordered_set<int> &forced_survive,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context,
    const std::unordered_map<int, double> *exact_source_times = nullptr,
    const std::unordered_map<int, std::pair<double, double>> *source_time_bounds =
        nullptr);

double evaluate_outcome_density_idx(
    const uuber::NativeContext &ctx, int outcome_idx, double t,
    const std::string &component, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    bool include_na_donors, int outcome_label_context_idx = -1) {
  if (outcome_idx < 0 ||
      outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return 0.0;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
  if (info.node_id < 0)
    return 0.0;
  const std::string &outcome_label_context =
      (outcome_label_context_idx >= 0 &&
       outcome_label_context_idx < static_cast<int>(ctx.outcome_labels.size()))
          ? ctx.outcome_labels[static_cast<std::size_t>(
                outcome_label_context_idx)]
          : std::string();
  std::vector<int> comp_filtered;
  const std::vector<int> &comp_use =
      filter_competitor_ids(ctx, info.competitor_ids, component,
                            comp_filtered);
  return node_density_with_competitors_internal(
      ctx, info.node_id, t, component, component_idx, kEmptyForcedSet,
      kEmptyForcedSet, comp_use, trial_params, trial_type_key,
      include_na_donors, outcome_label_context, outcome_label_context_idx);
}

double accumulate_component_guess_density(const uuber::NativeContext &ctx,
                                          const std::string &target_label,
                                          int target_outcome_idx, double t,
                                          const std::string &component,
                                          int component_idx,
                                          const TrialParamSet *trial_params,
                                          const std::string &trial_type_key) {
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return 0.0;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid)
    return 0.0;
  if (!guess.target.empty()) {
    if (guess.target == "__GUESS__") {
      if (target_label != "__GUESS__")
        return 0.0;
    } else if (guess.target_outcome_idx >= 0) {
      if (target_outcome_idx != guess.target_outcome_idx)
        return 0.0;
    } else {
      std::string guess_target = normalize_label_string(guess.target);
      std::string target_key = normalize_label_string(target_label);
      if (guess_target != target_key)
        return 0.0;
    }
  }
  double total = 0.0;
  for (const auto &entry : guess.keep_weights) {
    int donor_outcome_idx = entry.first;
    if (donor_outcome_idx < 0)
      continue;
    double keep_prob = entry.second;
    double release = 1.0 - keep_prob;
    if (release <= 0.0)
      continue;
    total +=
        release * evaluate_outcome_density_idx(
                      ctx, donor_outcome_idx, t, component, component_idx,
                      trial_params, trial_type_key, false, donor_outcome_idx);
  }
  return total;
}

double accumulate_outcome_alias_density(const uuber::NativeContext &ctx,
                                        int target_outcome_idx, double t,
                                        const std::string &component,
                                        int component_idx,
                                        const TrialParamSet *trial_params,
                                        const std::string &trial_type_key,
                                        bool include_na_donors) {
  if (target_outcome_idx < 0 ||
      target_outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
    return 0.0;
  }
  const uuber::OutcomeContextInfo &info =
      ctx.outcome_info[static_cast<std::size_t>(target_outcome_idx)];
  double total = 0.0;
  for (int source_idx : info.alias_sources) {
    total += evaluate_outcome_density_idx(
        ctx, source_idx, t, component, component_idx, trial_params,
        trial_type_key, include_na_donors, source_idx);
  }
  for (const auto &donor : info.guess_donors) {
    if (donor.rt_policy == "na" && !include_na_donors)
      continue;
    if (donor.outcome_idx < 0)
      continue;
    total += donor.weight * evaluate_outcome_density_idx(
                                ctx, donor.outcome_idx, t, component,
                                component_idx, trial_params, trial_type_key,
                                include_na_donors, donor.outcome_idx);
  }
  return total;
}

inline const TrialAccumulatorParams *
get_trial_param_entry(const TrialParamSet *trial_params, int acc_index) {
  if (!trial_params)
    return nullptr;
  if (acc_index < 0 ||
      acc_index >= static_cast<int>(trial_params->acc_params.size())) {
    return nullptr;
  }
  const TrialAccumulatorParams &entry = trial_params->acc_params[acc_index];
  if (!entry.has_override)
    return nullptr;
  return &entry;
}

bool component_active(const uuber::NativeAccumulator &acc,
                      const std::string &component, int component_idx,
                      const TrialAccumulatorParams *override = nullptr) {
  if (component_idx < 0 && (component.empty() || component == "__default__")) {
    return true;
  }
  const std::vector<int> *comps_idx = nullptr;
  if (override && override->has_components &&
      !override->component_indices.empty()) {
    comps_idx = &override->component_indices;
  } else if (!acc.component_indices.empty()) {
    comps_idx = &acc.component_indices;
  }
  if (comps_idx && component_idx >= 0) {
    return std::find(comps_idx->begin(), comps_idx->end(), component_idx) !=
           comps_idx->end();
  }
  const std::vector<std::string> *comps = nullptr;
  if (override && override->has_components) {
    comps = &override->components;
  } else if (!acc.components.empty()) {
    comps = &acc.components;
  }
  if (!comps || comps->empty())
    return true;
  return std::find(comps->begin(), comps->end(), component) != comps->end();
}

using ExactSourceTimeMap = std::unordered_map<int, double>;
using SourceTimeBoundsMap = std::unordered_map<int, std::pair<double, double>>;

struct PoolEvalScratch {
  std::vector<std::size_t> active_indices;
  std::vector<double> density;
  std::vector<double> survival;
};

struct PoolEvalScratchStack {
  std::vector<PoolEvalScratch> stack;
  std::size_t depth{0};
};

inline PoolEvalScratchStack &pool_eval_scratch_stack() {
  thread_local PoolEvalScratchStack scratch;
  return scratch;
}

class PoolEvalScratchGuard {
public:
  PoolEvalScratchGuard()
      : stack(pool_eval_scratch_stack()), index(stack.depth++) {
    if (index >= stack.stack.size()) {
      stack.stack.emplace_back();
    }
  }

  ~PoolEvalScratchGuard() {
    if (stack.depth > 0) {
      --stack.depth;
    }
  }

  PoolEvalScratch &scratch() { return stack.stack[index]; }

private:
  PoolEvalScratchStack &stack;
  std::size_t index;
};

NodeEvalResult
eval_event_label(const uuber::NativeContext &ctx, const std::string &label,
                 double t, const std::string &component,
                 const std::unordered_set<int> &forced_complete,
                 const std::unordered_set<int> &forced_survive, EvalNeed need,
                 const TrialParamSet *trial_params = nullptr,
                 const std::string &trial_type_key = std::string(),
                 bool include_na_donors = false,
                 const std::string &outcome_label_context = std::string(),
                 const LabelRef *label_ref = nullptr, int component_idx = -1,
                 int outcome_idx_context = -1,
                 const ExactSourceTimeMap *exact_source_times = nullptr,
                 const SourceTimeBoundsMap *source_time_bounds = nullptr,
                 const ForcedScopeFilter *forced_scope_filter = nullptr) {
  LabelRef resolved_ref;
  if (!label_ref) {
    resolved_ref = resolve_label_ref_runtime(ctx, label);
    label_ref = &resolved_ref;
  }
  if (component_idx < 0) {
    component_idx = component_index_of(ctx, component);
  }
  if (outcome_idx_context < 0 && !outcome_label_context.empty()) {
    outcome_idx_context = outcome_index_of(ctx, outcome_label_context);
  }
  if (label.empty()) {
    return make_node_result(need, 0.0, 1.0, 0.0);
  }
  if (label == "__DEADLINE__") {
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (needs_survival(need) || needs_cdf(need)) {
      survival = 0.0;
      cdf = 1.0;
      if (std::isfinite(t)) {
        if (t < 0.0) {
          survival = 1.0;
          cdf = 0.0;
        } else {
          survival = 0.0;
          cdf = 1.0;
        }
      }
    }
    return make_node_result(need, 0.0, survival, cdf);
  }
  double donor_density = 0.0;
  if (needs_density(need)) {
    int target_outcome_idx = -1;
    if (outcome_idx_context >= 0) {
      target_outcome_idx = outcome_idx_context;
    } else if (outcome_label_context.empty()) {
      target_outcome_idx = label_ref->outcome_idx;
    }
    const std::string &target_label =
        outcome_label_context.empty() ? label : outcome_label_context;
    donor_density += accumulate_component_guess_density(
        ctx, target_label, target_outcome_idx, t, component, component_idx,
        trial_params, trial_type_key);
    donor_density += accumulate_outcome_alias_density(
        ctx, target_outcome_idx, t, component, component_idx, trial_params,
        trial_type_key, include_na_donors);
  }
  if (label == "__GUESS__") {
    return make_node_result(need, donor_density, 0.0, 1.0);
  }
  int label_idx = label_ref->label_id;
  if (label_idx >= 0 && label_idx != NA_INTEGER) {
    if (forced_contains_scoped(forced_complete, label_idx,
                               forced_scope_filter)) {
      return make_node_result(need, 0.0, 0.0, 1.0);
    }
    if (forced_contains_scoped(forced_survive, label_idx,
                               forced_scope_filter)) {
      return make_node_result(need, 0.0, 1.0, 0.0);
    }
  }

  int acc_idx = label_ref->acc_idx;
  if (acc_idx >= 0 && acc_idx < static_cast<int>(ctx.accumulators.size())) {
    const uuber::NativeAccumulator &acc = ctx.accumulators[acc_idx];
    const TrialAccumulatorParams *override =
        get_trial_param_entry(trial_params, acc_idx);
    if (!component_active(acc, component, component_idx, override)) {
      return make_node_result(need, 0.0, 1.0, 0.0);
    }
    int onset_kind = override ? override->onset_kind : acc.onset_kind;
    int onset_source_acc_idx =
        override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
    int onset_source_pool_idx =
        override ? override->onset_source_pool_idx : acc.onset_source_pool_idx;
    double onset_lag = override ? override->onset_lag : acc.onset_lag;
    double onset = override ? override->onset : acc.onset;
    double q = override ? override->q : acc.q;
    const AccDistParams &cfg = override ? override->dist_cfg : acc.dist_cfg;
    double density = std::numeric_limits<double>::quiet_NaN();
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (!ctx.has_chained_onsets || onset_kind == uuber::ONSET_ABSOLUTE) {
      if (needs_density(need)) {
        density = acc_density_from_cfg(t, onset, q, cfg);
      }
      if (needs_survival(need)) {
        survival = acc_survival_from_cfg(t, onset, q, cfg);
      }
      if (needs_cdf(need)) {
        cdf = acc_cdf_success_from_cfg(t, onset, cfg);
      }
    } else {
      LabelRef source_ref;
      std::string source_label;
      if (onset_kind == uuber::ONSET_AFTER_ACCUMULATOR) {
        if (onset_source_acc_idx >= 0 &&
            onset_source_acc_idx < static_cast<int>(ctx.accumulators.size())) {
          const auto &src_acc =
              ctx.accumulators[static_cast<std::size_t>(onset_source_acc_idx)];
          source_label = src_acc.id;
          source_ref.acc_idx = onset_source_acc_idx;
          source_ref.label_id = label_id(ctx, source_label);
        }
      } else if (onset_kind == uuber::ONSET_AFTER_POOL) {
        if (onset_source_pool_idx >= 0 &&
            onset_source_pool_idx < static_cast<int>(ctx.pools.size())) {
          const auto &src_pool =
              ctx.pools[static_cast<std::size_t>(onset_source_pool_idx)];
          source_label = src_pool.id;
          source_ref.pool_idx = onset_source_pool_idx;
          source_ref.label_id = label_id(ctx, source_label);
        }
      }
      if (source_label.empty()) {
        return make_node_result(need, 0.0, 1.0, 0.0);
      }

      const double lag_total = onset_lag + onset;
      const double x_shift = lag_total + cfg.t0;
      const int source_label_id = source_ref.label_id;

      // If source is forced to survive at time t, dependent accumulator cannot
      // have completed by t.
      if (source_label_id >= 0 && source_label_id != NA_INTEGER &&
          forced_contains_scoped(forced_survive, source_label_id,
                                 forced_scope_filter)) {
        return make_node_result(need, 0.0, 1.0, 0.0);
      }

      double bound_lower = 0.0;
      double bound_upper = std::numeric_limits<double>::infinity();
      bool has_conditioning_bounds = false;
      bool has_exact = false;
      double exact_time = std::numeric_limits<double>::quiet_NaN();

      if (source_label_id >= 0 && source_label_id != NA_INTEGER) {
        if (forced_contains_scoped(forced_complete, source_label_id,
                                   forced_scope_filter)) {
          bound_upper = std::min(bound_upper, t);
          has_conditioning_bounds = true;
        }
        if (exact_source_times) {
          auto it = exact_source_times->find(source_label_id);
          if (it != exact_source_times->end() && std::isfinite(it->second)) {
            has_exact = true;
            exact_time = it->second;
          }
        }
        if (source_time_bounds) {
          auto it = source_time_bounds->find(source_label_id);
          if (it != source_time_bounds->end()) {
            if (std::isfinite(it->second.first)) {
              bound_lower = std::max(bound_lower, it->second.first);
            }
            if (std::isfinite(it->second.second)) {
              bound_upper = std::min(bound_upper, it->second.second);
            }
            has_conditioning_bounds = true;
          }
        }
      }

      if (has_exact) {
        double x = t - exact_time - x_shift;
        double cdf_success = clamp_probability(eval_cdf_single(cfg, x));
        if (needs_density(need)) {
          density = (1.0 - q) * eval_pdf_single(cfg, x);
          density = safe_density(density);
        }
        if (needs_cdf(need)) {
          cdf = clamp_probability((1.0 - q) * cdf_success);
        }
        if (needs_survival(need)) {
          survival = clamp_probability(q + (1.0 - q) * (1.0 - cdf_success));
        }
      } else {
        if (bound_upper <= bound_lower) {
          return make_node_result(need, 0.0, 1.0, 0.0);
        }
        if (std::isfinite(t) && (t - x_shift <= bound_lower)) {
          return make_node_result(need, 0.0, 1.0, 0.0);
        }

        auto source_eval = [&](double u, EvalNeed source_need) -> NodeEvalResult {
          return eval_event_label(
              ctx, "__CHAIN_SOURCE__", u, component, forced_complete,
              forced_survive, source_need, trial_params, trial_type_key, false,
              std::string(), &source_ref, component_idx, -1,
              exact_source_times, source_time_bounds, forced_scope_filter);
        };
        auto source_cdf_at = [&](double u) -> double {
          if (u <= 0.0) {
            return 0.0;
          }
          NodeEvalResult source =
              source_eval(u, EvalNeed::kCDF);
          return clamp_probability(source.cdf);
        };

        double source_mass = 1.0;
        if (has_conditioning_bounds) {
          double lower_cdf =
              std::isfinite(bound_lower) ? source_cdf_at(bound_lower) : 0.0;
          double upper_cdf = std::isfinite(bound_upper)
                                 ? source_cdf_at(bound_upper)
                                 : source_cdf_at(
                                       std::numeric_limits<double>::infinity());
          source_mass = upper_cdf - lower_cdf;
          if (!std::isfinite(source_mass) || source_mass <= 0.0) {
            return make_node_result(need, 0.0, 1.0, 0.0);
          }
        }

        double upper_int =
            std::min(bound_upper, std::max(bound_lower, t - x_shift));
        if (upper_int <= bound_lower) {
          return make_node_result(need, 0.0, 1.0, 0.0);
        }

        auto source_pdf = [&](double u) -> double {
          NodeEvalResult source =
              source_eval(u, EvalNeed::kDensity);
          return safe_density(source.density);
        };

        if (!std::isfinite(t)) {
          double cdf_total = 0.0;
          if (has_conditioning_bounds) {
            cdf_total = clamp_probability(1.0 - q);
          } else {
            NodeEvalResult source_inf =
                source_eval(std::numeric_limits<double>::infinity(),
                            EvalNeed::kCDF);
            double source_finish = clamp_probability(source_inf.cdf);
            cdf_total = clamp_probability((1.0 - q) * source_finish);
          }
          if (needs_density(need)) {
            density = 0.0;
          }
          if (needs_cdf(need)) {
            cdf = cdf_total;
          }
          if (needs_survival(need)) {
            survival = clamp_probability(1.0 - cdf_total);
          }
          NodeEvalResult res = make_node_result(need, density, survival, cdf);
          if (needs_density(need) && donor_density > 0.0) {
            res.density = safe_density(res.density + donor_density);
          }
          return res;
        }

        auto cond_source_pdf = [&](double u) -> double {
          if (!std::isfinite(u) || u <= bound_lower || u > bound_upper) {
            return 0.0;
          }
          double dens = source_pdf(u);
          if (dens <= 0.0) {
            return 0.0;
          }
          if (!has_conditioning_bounds) {
            return dens;
          }
          return dens / source_mass;
        };

        if (needs_density(need)) {
          auto integrand_density = [&](double u) -> double {
            double fs = cond_source_pdf(u);
            if (fs <= 0.0) {
              return 0.0;
            }
            double x = t - u - x_shift;
            double fx = eval_pdf_single(cfg, x);
            if (!std::isfinite(fx) || fx <= 0.0) {
              return 0.0;
            }
            return fs * fx;
          };
          density = (1.0 - q) * uuber::integrate_boost_fn(
                                    integrand_density, bound_lower, upper_int,
                                    kDefaultRelTol, kDefaultAbsTol,
                                    kDefaultMaxDepth);
          density = safe_density(density);
        }

        if (needs_cdf(need) || needs_survival(need)) {
          auto integrand_cdf = [&](double u) -> double {
            double fs = cond_source_pdf(u);
            if (fs <= 0.0) {
              return 0.0;
            }
            double x = t - u - x_shift;
            double Fx = eval_cdf_single(cfg, x);
            if (!std::isfinite(Fx) || Fx <= 0.0) {
              return 0.0;
            }
            return fs * clamp_probability(Fx);
          };
          double cdf_cond = uuber::integrate_boost_fn(
              integrand_cdf, bound_lower, upper_int, kDefaultRelTol,
              kDefaultAbsTol, kDefaultMaxDepth);
          cdf_cond = clamp_probability(cdf_cond);
          double cdf_total = clamp_probability((1.0 - q) * cdf_cond);
          if (needs_cdf(need)) {
            cdf = cdf_total;
          }
          if (needs_survival(need)) {
            survival = clamp_probability(1.0 - cdf_total);
          }
        }
      }
    }
    NodeEvalResult res = make_node_result(need, density, survival, cdf);
    if (needs_density(need) && donor_density > 0.0) {
      res.density = safe_density(res.density + donor_density);
    }
    return res;
  }

  int pool_idx = label_ref->pool_idx;
  if (pool_idx >= 0 && pool_idx < static_cast<int>(ctx.pools.size())) {
    const uuber::NativePool &pool = ctx.pools[pool_idx];
    if (pool.members.empty()) {
      NodeEvalResult res = make_node_result(need, 0.0, 1.0, 0.0);
      if (needs_density(need) && donor_density > 0.0) {
        res.density = safe_density(res.density + donor_density);
      }
      return res;
    }
    PoolEvalScratchGuard scratch_guard;
    PoolEvalScratch &scratch = scratch_guard.scratch();
    scratch.active_indices.clear();
    scratch.active_indices.reserve(pool.members.size());
    for (std::size_t member_idx = 0; member_idx < pool.members.size();
         ++member_idx) {
      const LabelRef &member_ref = pool.member_refs[member_idx];
      int member_acc_idx = member_ref.acc_idx;
      if (member_acc_idx >= 0 &&
          member_acc_idx < static_cast<int>(ctx.accumulators.size())) {
        const uuber::NativeAccumulator &acc = ctx.accumulators[member_acc_idx];
        const TrialAccumulatorParams *override =
            get_trial_param_entry(trial_params, member_acc_idx);
        if (!component_active(acc, component, component_idx, override))
          continue;
      }
      scratch.active_indices.push_back(member_idx);
    }
    if (scratch.active_indices.empty()) {
      NodeEvalResult res = make_node_result(need, 0.0, 1.0, 0.0);
      if (needs_density(need) && donor_density > 0.0) {
        res.density = safe_density(res.density + donor_density);
      }
      return res;
    }
    const bool need_density = needs_density(need);
    const bool need_survival =
        needs_survival(need) || needs_cdf(need) || need_density;
    EvalNeed child_need = need_density
                              ? (EvalNeed::kDensity | EvalNeed::kSurvival)
                              : EvalNeed::kSurvival;
    if (need_density) {
      scratch.density.resize(scratch.active_indices.size());
    }
    if (need_survival) {
      scratch.survival.resize(scratch.active_indices.size());
    }
    for (std::size_t i = 0; i < scratch.active_indices.size(); ++i) {
      std::size_t member_idx = scratch.active_indices[i];
      const std::string &member = pool.members[member_idx];
      const LabelRef &member_ref = pool.member_refs[member_idx];
      NodeEvalResult child = eval_event_label(
          ctx, member, t, component, forced_complete, forced_survive,
          child_need, trial_params, std::string(), false, std::string(),
          &member_ref, component_idx, -1, exact_source_times,
          source_time_bounds, forced_scope_filter);
      if (need_density) {
        scratch.density[i] = child.density;
      }
      if (need_survival) {
        scratch.survival[i] = child.survival;
      }
    }
    double density = std::numeric_limits<double>::quiet_NaN();
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (need_density) {
      density = pool_density_fast(scratch.density, scratch.survival, pool.k);
      if (!std::isfinite(density) || density < 0.0)
        density = 0.0;
    }
    if (needs_survival(need) || needs_cdf(need)) {
      survival = pool_survival_fast(scratch.survival, pool.k);
      if (!std::isfinite(survival))
        survival = 0.0;
    }
    if (needs_cdf(need)) {
      cdf = clamp_probability(1.0 - survival);
    }
    NodeEvalResult res = make_node_result(need, density, survival, cdf);
    if (needs_density(need) && donor_density > 0.0) {
      res.density = safe_density(res.density + donor_density);
    }
    return res;
  }

  return make_node_result(need, 0.0, 1.0, 0.0);
}

struct ForcedKey {
  int complete_id{0};
  int survive_id{0};
};

struct TimeKey {
  double value{0.0};
};

std::string component_cache_key(const std::string &component,
                                const std::string &trial_key);

struct ComponentCacheEntry {
  std::string trial_type_key;
};

ComponentCacheEntry make_component_cache_entry(const std::string &component,
                                               const std::string &trial_key) {
  ComponentCacheEntry entry;
  entry.trial_type_key = component_cache_key(component, trial_key);
  return entry;
}

ComponentCacheEntry
default_component_cache_entry(const std::string &component) {
  return make_component_cache_entry(component, std::string());
}

double accumulate_plan_guess_density(const uuber::NativeContext &ctx,
                                     const Rcpp::List &donors, double t,
                                     const std::string &component,
                                     const std::string &trial_type_key,
                                     const TrialParamSet *trial_params,
                                     bool include_na_donors) {
  int component_idx = component_index_of(ctx, component);
  double total = 0.0;
  for (R_xlen_t i = 0; i < donors.size(); ++i) {
    Rcpp::RObject donor_obj = donors[i];
    if (donor_obj.isNULL())
      continue;
    Rcpp::List donor(donor_obj);
    std::string rt_policy = donor.containsElementNamed("rt_policy")
                                ? Rcpp::as<std::string>(donor["rt_policy"])
                                : std::string("keep");
    bool donor_na = rt_policy == "na";
    if (donor_na != include_na_donors)
      continue;
    double release = donor.containsElementNamed("release")
                         ? Rcpp::as<double>(donor["release"])
                         : 0.0;
    if (!std::isfinite(release))
      release = 0.0;
    if (release <= 0.0)
      continue;
    int donor_node = donor.containsElementNamed("node_id")
                         ? Rcpp::as<int>(donor["node_id"])
                         : NA_INTEGER;
    if (donor_node == NA_INTEGER)
      continue;
    Rcpp::IntegerVector comp_ids =
        donor.containsElementNamed("competitor_ids")
            ? Rcpp::IntegerVector(donor["competitor_ids"])
            : Rcpp::IntegerVector();
    std::vector<int> comp_vec = integer_vector_to_std(comp_ids, false);
    std::vector<int> comp_filtered;
    const std::vector<int> &comp_use =
        filter_competitor_ids(ctx, comp_vec, component, comp_filtered);
    std::string donor_label;
    if (donor.containsElementNamed("source_label") &&
        !Rf_isNull(donor["source_label"])) {
      donor_label = Rcpp::as<std::string>(donor["source_label"]);
    }
    int outcome_idx_context =
        donor_label.empty() ? -1 : outcome_index_of(ctx, donor_label);
    double density = node_density_with_competitors_internal(
        ctx, donor_node, t, component, component_idx, kEmptyForcedSet,
        kEmptyForcedSet, comp_use, trial_params, trial_type_key,
        include_na_donors, donor_label, outcome_idx_context);
    if (!std::isfinite(density) || density <= 0.0)
      continue;
    total += release * density;
  }
  return total;
}

std::vector<ComponentCacheEntry>
build_component_cache_entries(const std::vector<std::string> &components) {
  std::vector<ComponentCacheEntry> entries;
  entries.reserve(components.size());
  for (const auto &comp : components) {
    entries.push_back(default_component_cache_entry(comp));
  }
  return entries;
}

std::string component_cache_key(const std::string &component,
                                const std::string &trial_key = std::string()) {
  if (!trial_key.empty())
    return trial_key;
  if (component.empty())
    return "__default__";
  return component;
}

std::string normalize_label_string(const std::string &label) {
  if (label.empty())
    return std::string();
  return label;
}

double component_keep_weight(const uuber::NativeContext &ctx, int component_idx,
                             int outcome_idx) {
  if (component_idx < 0 ||
      component_idx >= static_cast<int>(ctx.component_info.size())) {
    return 1.0;
  }
  const uuber::ComponentGuessPolicy &guess =
      ctx.component_info[static_cast<std::size_t>(component_idx)].guess;
  if (!guess.valid)
    return 1.0;
  if (outcome_idx == kOutcomeIdxNA) {
    return guess.has_keep_weight_na ? guess.keep_weight_na : 1.0;
  }
  for (const auto &entry : guess.keep_weights) {
    if (entry.first == outcome_idx)
      return entry.second;
  }
  return 1.0;
}

inline double extract_trial_id(const Rcpp::DataFrame *trial_rows) {
  if (!trial_rows)
    return std::numeric_limits<double>::quiet_NaN();
  if (!trial_rows->containsElementNamed("trial"))
    return std::numeric_limits<double>::quiet_NaN();
  Rcpp::NumericVector trial_col((*trial_rows)["trial"]);
  if (trial_col.size() == 0)
    return std::numeric_limits<double>::quiet_NaN();
  double val = trial_col[0];
  if (Rcpp::NumericVector::is_na(val))
    return std::numeric_limits<double>::quiet_NaN();
  return static_cast<double>(val);
}

struct NodeEvalState {
  NodeEvalState(const uuber::NativeContext &ctx_, double time,
                const std::string &component_label,
                const std::unordered_set<int> &forced_complete_ref,
                const std::unordered_set<int> &forced_survive_ref,
                const TrialParamSet *params_ptr = nullptr,
                const std::string &trial_key = std::string(),
                bool include_na = false,
                const std::string &outcome_label_ctx = std::string(),
                int component_idx_val = -1, int outcome_idx_val = -1,
                const ExactSourceTimeMap *exact_source_times_ptr = nullptr,
                const SourceTimeBoundsMap *source_time_bounds_ptr = nullptr,
                const ForcedScopeFilter *forced_scope_filter_ptr = nullptr)
      : ctx(ctx_), t(time), component(component_label),
        forced_complete(forced_complete_ref),
        forced_survive(forced_survive_ref), trial_params(params_ptr),
        trial_type_key(component_cache_key(component_label, trial_key)),
        include_na_donors(include_na), outcome_label(outcome_label_ctx),
        component_idx(component_idx_val), outcome_idx(outcome_idx_val),
        exact_source_times(exact_source_times_ptr),
        source_time_bounds(source_time_bounds_ptr),
        forced_scope_filter(forced_scope_filter_ptr) {}

  const uuber::NativeContext &ctx;
  double t;
  std::string component;
  const std::unordered_set<int> &forced_complete;
  const std::unordered_set<int> &forced_survive;
  const TrialParamSet *trial_params;
  std::string trial_type_key;
  bool include_na_donors;
  std::string outcome_label;
  int component_idx;
  int outcome_idx;
  const ExactSourceTimeMap *exact_source_times;
  const SourceTimeBoundsMap *source_time_bounds;
  const ForcedScopeFilter *forced_scope_filter;
};

inline const TrialAccumulatorParams *
get_state_params(const NodeEvalState &state, int acc_index) {
  return get_trial_param_entry(state.trial_params, acc_index);
}

std::string member_ids_signature(const std::vector<int> &ids) {
  if (ids.empty())
    return ".";
  std::ostringstream oss;
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (i > 0)
      oss << ",";
    oss << ids[i];
  }
  return oss.str();
}

std::string shared_ids_signature(const std::vector<std::string> &shared_ids) {
  if (shared_ids.empty())
    return ".";
  std::ostringstream oss;
  for (std::size_t i = 0; i < shared_ids.size(); ++i) {
    if (i > 0)
      oss << ",";
    if (shared_ids[i].empty()) {
      oss << ".";
    } else {
      oss << shared_ids[i].size() << ":" << shared_ids[i];
    }
  }
  return oss.str();
}

std::vector<int>
build_shared_index(const std::vector<std::string> &shared_ids) {
  std::vector<int> shared_index(shared_ids.size(), -1);
  std::unordered_map<std::string, int> group_map;
  int next_group = 1;
  for (std::size_t i = 0; i < shared_ids.size(); ++i) {
    const std::string &sid = shared_ids[i];
    if (sid.empty())
      continue;
    auto it = group_map.find(sid);
    if (it == group_map.end()) {
      group_map[sid] = next_group;
      shared_index[i] = next_group;
      ++next_group;
    } else {
      shared_index[i] = it->second;
    }
  }
  return shared_index;
}

std::string pool_template_cache_key(const std::string &pool_label,
                                    const std::string &component, int k,
                                    const std::vector<int> &member_ids,
                                    const std::string &shared_signature) {
  std::ostringstream oss;
  oss << pool_label << "|" << component_cache_key(component) << "|" << k << "|"
      << member_ids_signature(member_ids) << "|" << shared_signature;
  return oss.str();
}

struct IntegrationSettings {
  double rel_tol{kDefaultRelTol};
  double abs_tol{kDefaultAbsTol};
  int max_depth{kDefaultMaxDepth};
};

struct GuardEvalInput {
  const uuber::NativeContext &ctx;
  const uuber::NativeNode &node;
  std::string component;
  int component_idx{-1};
  std::string trial_type_key;
  const std::unordered_set<int> *forced_complete{nullptr};
  const std::unordered_set<int> *forced_survive{nullptr};
  const ForcedScopeFilter *forced_scope_filter{nullptr};
  ForcedScopeFilter local_scope_filter{};
  bool has_scoped_forced{false};
  const TrialParamSet *trial_params;
  const ExactSourceTimeMap *exact_source_times{nullptr};
  const SourceTimeBoundsMap *source_time_bounds{nullptr};
};

const uuber::NativeNode &fetch_node(const uuber::NativeContext &ctx,
                                    int node_id);
NodeEvalResult
eval_node_with_forced(const uuber::NativeContext &ctx, int node_id, double time,
                      const std::string &component, int component_idx,
                      const std::unordered_set<int> &forced_complete,
                      const std::unordered_set<int> &forced_survive,
                      EvalNeed need, const TrialParamSet *trial_params,
                      const std::string &trial_key,
                      const ExactSourceTimeMap *exact_source_times = nullptr,
                      const SourceTimeBoundsMap *source_time_bounds = nullptr,
                      const ForcedScopeFilter *forced_scope_filter = nullptr);
GuardEvalInput make_guard_input(const uuber::NativeContext &ctx,
                                const uuber::NativeNode &node,
                                const std::string &component, int component_idx,
                                const std::unordered_set<int> &forced_complete,
                                const std::unordered_set<int> &forced_survive,
                                const std::string &trial_key,
                                const TrialParamSet *trial_params,
                                const ExactSourceTimeMap *exact_source_times,
                                const SourceTimeBoundsMap *source_time_bounds,
                                const ForcedScopeFilter *forced_scope_filter);
std::vector<int> gather_blocker_sources(const uuber::NativeContext &ctx,
                                        const uuber::NativeNode &guard_node);

std::vector<ScenarioRecord> compute_node_scenarios(int node_id,
                                                   NodeEvalState &state);

int first_missing_id(const std::vector<int> &forced_complete,
                     const std::vector<int> &required) {
  if (required.empty())
    return NA_INTEGER;
  for (int id : required) {
    if (!vector_contains(forced_complete, id)) {
      return id;
    }
  }
  return NA_INTEGER;
}

std::vector<ScenarioRecord>
compute_pool_event_scenarios(const uuber::NativeNode &node,
                             NodeEvalState &state) {
  const std::string &pool_label = node.source;
  if (pool_label.empty())
    return {};
  int pool_idx = node.source_ref.pool_idx;
  if (pool_idx < 0 || pool_idx >= static_cast<int>(state.ctx.pools.size()))
    return {};
  const uuber::NativePool &pool = state.ctx.pools[pool_idx];
  if (pool.members.empty())
    return {};
  int k = pool.k;
  if (k < 1)
    return {};

  std::vector<std::size_t> active_indices;
  active_indices.reserve(pool.members.size());
  for (std::size_t member_idx = 0; member_idx < pool.members.size();
       ++member_idx) {
    const LabelRef &member_ref = pool.member_refs[member_idx];
    int acc_idx = member_ref.acc_idx;
    if (acc_idx >= 0 &&
        acc_idx < static_cast<int>(state.ctx.accumulators.size())) {
      const uuber::NativeAccumulator &acc = state.ctx.accumulators[acc_idx];
      const TrialAccumulatorParams *override = get_state_params(state, acc_idx);
      if (!component_active(acc, state.component, state.component_idx,
                            override))
        continue;
    }
    active_indices.push_back(member_idx);
  }
  std::size_t n = active_indices.size();
  if (n == 0 || k > static_cast<int>(n))
    return {};

  Rcpp::NumericVector dens_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector cdf_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector surv_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector cdf_success_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector surv_success_vec(static_cast<R_xlen_t>(n));
  std::vector<int> member_ids_std(n, NA_INTEGER);
  std::vector<std::string> shared_ids(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t member_idx = active_indices[i];
    const LabelRef &member_ref = pool.member_refs[member_idx];
    if (member_ref.label_id >= 0 && member_ref.label_id != NA_INTEGER) {
      member_ids_std[i] = member_ref.label_id;
    }
    int acc_idx = member_ref.acc_idx;
    if (acc_idx >= 0 &&
        acc_idx < static_cast<int>(state.ctx.accumulators.size())) {
      const uuber::NativeAccumulator &acc = state.ctx.accumulators[acc_idx];
      const TrialAccumulatorParams *override = get_state_params(state, acc_idx);
      const std::string &shared_id =
          override ? override->shared_trigger_id : acc.shared_trigger_id;
      if (!shared_id.empty()) {
        shared_ids[i] = shared_id;
      }
    }
  }

  for (std::size_t i = 0; i < n; ++i) {
    std::size_t member_idx = active_indices[i];
    const std::string &member_label = pool.members[member_idx];
    const LabelRef &member_ref = pool.member_refs[member_idx];
    NodeEvalResult eval =
        eval_event_label(state.ctx, member_label, state.t, state.component,
                         state.forced_complete, state.forced_survive, kEvalAll,
                         state.trial_params, state.trial_type_key, false,
                         std::string(), &member_ref, state.component_idx, -1,
                         state.exact_source_times, state.source_time_bounds,
                         state.forced_scope_filter);
    double density = eval.density;
    double cdf = clamp_probability(eval.cdf);
    double survival = clamp_probability(eval.survival);
    double cdf_success = cdf;
    double surv_success = survival;
    int acc_idx = member_ref.acc_idx;
    if (acc_idx >= 0 &&
        acc_idx < static_cast<int>(state.ctx.accumulators.size())) {
      const uuber::NativeAccumulator &acc = state.ctx.accumulators[acc_idx];
      const TrialAccumulatorParams *override = get_state_params(state, acc_idx);
      const std::string &shared_id =
          override ? override->shared_trigger_id : acc.shared_trigger_id;
      if (!shared_id.empty()) {
        double onset = override ? override->onset : acc.onset;
        const AccDistParams &cfg = override ? override->dist_cfg : acc.dist_cfg;
        cdf_success =
            clamp_probability(acc_cdf_success_from_cfg(state.t, onset, cfg));
        surv_success = clamp_probability(1.0 - cdf_success);
        cdf = cdf_success;
        survival = surv_success;
      }
    }
    dens_vec[static_cast<R_xlen_t>(i)] = density;
    cdf_vec[static_cast<R_xlen_t>(i)] = cdf;
    surv_vec[static_cast<R_xlen_t>(i)] = survival;
    cdf_success_vec[static_cast<R_xlen_t>(i)] = cdf_success;
    surv_success_vec[static_cast<R_xlen_t>(i)] = surv_success;
  }

  int pool_label_id = node.source_ref.label_id;
  if (pool_label_id < 0) {
    pool_label_id = NA_INTEGER;
  }
  std::string shared_signature = shared_ids_signature(shared_ids);
  std::string template_key =
      pool_template_cache_key(pool_label, component_cache_key(state.component),
                              k, member_ids_std, shared_signature);
  uuber::PoolTemplateCacheEntry *template_entry = nullptr;
  auto tmpl_it = state.ctx.pool_template_cache.find(template_key);
  if (tmpl_it != state.ctx.pool_template_cache.end()) {
    template_entry = &tmpl_it->second;
  } else {
    uuber::PoolTemplateCacheEntry cache_entry = build_pool_template_cache(
        static_cast<int>(n), member_ids_std, pool_label_id, k);
    cache_entry.shared_index = build_shared_index(shared_ids);
    template_entry = &state.ctx.pool_template_cache
                          .emplace(template_key, std::move(cache_entry))
                          .first->second;
  }
  if (!template_entry) {
    return {};
  }
  if (template_entry->shared_index.size() != n) {
    template_entry->shared_index = build_shared_index(shared_ids);
  }
  if (template_entry->templates.empty()) {
    return {};
  }

  std::vector<int> base_forced_complete =
      set_to_sorted_vector(state.forced_complete);
  std::vector<int> base_forced_survive =
      set_to_sorted_vector(state.forced_survive);

  Rcpp::List combine_res = pool_density_combine_native(
      dens_vec, cdf_vec, surv_vec, cdf_success_vec, surv_success_vec,
      template_entry->shared_index, template_entry->templates,
      base_forced_complete, base_forced_survive);
  Rcpp::List scenario_list = combine_res["scenarios"];
  return rcpp_scenarios_to_records(scenario_list);
}

std::vector<ScenarioRecord>
compute_event_scenarios(const uuber::NativeNode &node, NodeEvalState &state) {
  std::vector<ScenarioRecord> out;
  if (node.source.empty()) {
    return out;
  }
  if (node.source_ref.pool_idx >= 0 &&
      node.source_ref.pool_idx < static_cast<int>(state.ctx.pools.size())) {
    std::vector<ScenarioRecord> pool_records =
        compute_pool_event_scenarios(node, state);
    if (!pool_records.empty()) {
      return pool_records;
    }
  }
  NodeEvalResult eval = eval_event_label(
      state.ctx, node.source, state.t, state.component, state.forced_complete,
      state.forced_survive, EvalNeed::kDensity, state.trial_params,
      state.trial_type_key, false, std::string(), &node.source_ref,
      state.component_idx, -1, state.exact_source_times,
      state.source_time_bounds, state.forced_scope_filter);
  double weight = eval.density;
  if (!std::isfinite(weight) || weight <= 0.0) {
    return out;
  }
  std::vector<int> base_fc = set_to_sorted_vector(state.forced_complete);
  std::vector<int> base_fs = set_to_sorted_vector(state.forced_survive);
  std::vector<int> source_ids = ensure_source_ids(state.ctx, node);
  std::vector<int> merged_fc = union_vectors(base_fc, source_ids);
  ScenarioRecord rec = make_scenario_record(weight, merged_fc, base_fs);
  if (rec.weight > 0.0) {
    out.push_back(std::move(rec));
  }
  return out;
}

std::vector<ScenarioRecord> compute_and_scenarios(const uuber::NativeNode &node,
                                                  NodeEvalState &state) {
  std::vector<ScenarioRecord> out;
  const std::vector<int> &child_ids = node.args;
  if (child_ids.empty())
    return out;
  std::vector<const uuber::NativeNode *> child_nodes;
  std::vector<std::vector<int>> child_sources;
  child_nodes.reserve(child_ids.size());
  child_sources.reserve(child_ids.size());
  for (int child_id : child_ids) {
    const uuber::NativeNode &child = fetch_node(state.ctx, child_id);
    child_nodes.push_back(&child);
    child_sources.push_back(ensure_source_ids(state.ctx, child));
  }
  for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
    std::vector<ScenarioRecord> child_scen =
        compute_node_scenarios(child_ids[idx], state);
    if (child_scen.empty())
      continue;
    for (const auto &sc : child_scen) {
      if (!std::isfinite(sc.weight) || sc.weight <= 0.0)
        continue;
      double weight = sc.weight;
      std::vector<int> forced_complete = sc.forced_complete;
      std::vector<int> forced_survive = sc.forced_survive;
      bool ok = true;
      for (std::size_t other_idx = 0; other_idx < child_ids.size();
           ++other_idx) {
        if (other_idx == idx)
          continue;
        std::unordered_set<int> fc_set = make_forced_set(forced_complete);
        std::unordered_set<int> fs_set = make_forced_set(forced_survive);
        NodeEvalResult sibling = eval_node_with_forced(
            state.ctx, child_ids[other_idx], state.t, state.component,
            state.component_idx, fc_set, fs_set, EvalNeed::kCDF,
            state.trial_params, state.trial_type_key,
            state.exact_source_times, state.source_time_bounds,
            state.forced_scope_filter);
        double Fj = clamp_probability(sibling.cdf);
        if (!std::isfinite(Fj) || Fj <= 0.0) {
          ok = false;
          break;
        }
        weight *= Fj;
        if (!std::isfinite(weight) || weight <= 0.0) {
          ok = false;
          break;
        }
        forced_complete =
            union_vectors(forced_complete, child_sources[other_idx]);
      }
      if (!ok)
        continue;
      ScenarioRecord merged =
          make_scenario_record(weight, forced_complete, forced_survive);
      if (merged.weight > 0.0) {
        out.push_back(std::move(merged));
      }
    }
  }
  return out;
}

std::vector<ScenarioRecord> compute_or_scenarios(const uuber::NativeNode &node,
                                                 NodeEvalState &state) {
  std::vector<ScenarioRecord> out;
  const std::vector<int> &child_ids = node.args;
  if (child_ids.empty())
    return out;
  std::vector<const uuber::NativeNode *> child_nodes;
  std::vector<std::vector<int>> child_sources;
  child_nodes.reserve(child_ids.size());
  child_sources.reserve(child_ids.size());
  for (int child_id : child_ids) {
    const uuber::NativeNode &child = fetch_node(state.ctx, child_id);
    child_nodes.push_back(&child);
    child_sources.push_back(ensure_source_ids(state.ctx, child));
  }
  for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
    std::vector<ScenarioRecord> child_scen =
        compute_node_scenarios(child_ids[idx], state);
    if (child_scen.empty())
      continue;
    for (const auto &sc : child_scen) {
      if (!std::isfinite(sc.weight) || sc.weight <= 0.0)
        continue;
      double weight = sc.weight;
      std::vector<int> forced_complete = sc.forced_complete;
      std::vector<int> forced_survive = sc.forced_survive;
      bool valid = true;
      for (std::size_t other_idx = 0; other_idx < child_ids.size();
           ++other_idx) {
        if (other_idx == idx)
          continue;
        const std::vector<int> &required = child_sources[other_idx];
        if (required.empty())
          continue;
        bool all_forced = true;
        for (int req : required) {
          if (!vector_contains(forced_complete, req)) {
            all_forced = false;
            break;
          }
        }
        if (all_forced) {
          valid = false;
          break;
        }
        int witness = first_missing_id(forced_complete, required);
        if (witness == NA_INTEGER) {
          valid = false;
          break;
        }
        std::vector<int> witness_vec{witness};
        forced_survive = union_vectors(forced_survive, witness_vec);
      }
      if (!valid)
        continue;
      ScenarioRecord merged =
          make_scenario_record(weight, forced_complete, forced_survive);
      if (merged.weight > 0.0) {
        out.push_back(std::move(merged));
      }
    }
  }
  return out;
}

double guard_effective_survival_internal(const GuardEvalInput &input, double t,
                                         const IntegrationSettings &settings);

double guard_density_internal(const GuardEvalInput &input, double t,
                              const IntegrationSettings &settings);

double guard_cdf_internal(const GuardEvalInput &input, double t,
                          const IntegrationSettings &settings);

std::vector<ScenarioRecord>
compute_guard_scenarios(const uuber::NativeNode &node, NodeEvalState &state) {
  std::vector<ScenarioRecord> out;
  if (node.reference_id < 0) {
    return out;
  }
  std::vector<ScenarioRecord> ref_scen =
      compute_node_scenarios(node.reference_id, state);
  if (ref_scen.empty())
    return out;
  std::vector<int> blocker_sources = gather_blocker_sources(state.ctx, node);
  IntegrationSettings settings;
  out.reserve(ref_scen.size());
  for (const auto &sc : ref_scen) {
    if (!std::isfinite(sc.weight) || sc.weight <= 0.0)
      continue;
    std::unordered_set<int> fc_set = make_forced_set(sc.forced_complete);
    std::unordered_set<int> fs_set = make_forced_set(sc.forced_survive);
    GuardEvalInput guard_input = make_guard_input(
        state.ctx, node, state.component, state.component_idx, fc_set, fs_set,
        state.trial_type_key, state.trial_params, state.exact_source_times,
        state.source_time_bounds, state.forced_scope_filter);
    double eff =
        guard_effective_survival_internal(guard_input, state.t, settings);
    if (!std::isfinite(eff) || eff <= 0.0)
      continue;
    double weight = sc.weight * eff;
    if (!std::isfinite(weight) || weight <= 0.0)
      continue;
    ScenarioRecord rec;
    rec.weight = weight;
    rec.forced_complete = sc.forced_complete;
    rec.forced_survive = union_vectors(sc.forced_survive, blocker_sources);
    out.push_back(std::move(rec));
  }
  return out;
}

std::vector<ScenarioRecord> compute_node_scenarios(int node_id,
                                                   NodeEvalState &state) {
  const uuber::NativeNode &node = fetch_node(state.ctx, node_id);
  std::vector<ScenarioRecord> result;
  if (node.kind == "event") {
    result = compute_event_scenarios(node, state);
  } else if (node.kind == "and") {
    result = compute_and_scenarios(node, state);
  } else if (node.kind == "or") {
    result = compute_or_scenarios(node, state);
  } else if (node.kind == "guard") {
    result = compute_guard_scenarios(node, state);
  }
  return result;
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState &state,
                                   EvalNeed need);

const uuber::NativeNode &fetch_node(const uuber::NativeContext &ctx,
                                    int node_id) {
  return ctx.nodes.at(ctx.node_index.at(node_id));
}

bool share_sources(const std::vector<int> &a, const std::vector<int> &b);

struct SharedGatePair {
  std::string x;
  std::string y;
  std::string c;
  LabelRef x_ref;
  LabelRef y_ref;
  LabelRef c_ref;
};

struct SharedGateNWay {
  std::string gate;
  std::string target_x;
  std::vector<std::string> competitor_xs;
  LabelRef gate_ref;
  LabelRef target_ref;
  std::vector<LabelRef> competitor_refs;
};

struct SharedGatePairCacheEntry {
  bool computed{false};
  bool is_shared{false};
  SharedGatePair spec;
};

struct SharedGateNWayCacheEntry {
  bool computed{false};
  bool is_shared{false};
  SharedGateNWay spec;
};

bool detect_shared_gate_pair(const uuber::NativeContext &ctx, int node_id,
                             int competitor_node_id, SharedGatePair &out) {
  if (node_id < 0 || competitor_node_id < 0)
    return false;
  const uuber::NativeNode &target = fetch_node(ctx, node_id);
  const uuber::NativeNode &competitor = fetch_node(ctx, competitor_node_id);
  if (target.kind != "and" || competitor.kind != "and")
    return false;
  if (target.args.size() != 2 || competitor.args.size() != 2)
    return false;

  const uuber::NativeNode &t0 = fetch_node(ctx, target.args[0]);
  const uuber::NativeNode &t1 = fetch_node(ctx, target.args[1]);
  const uuber::NativeNode &c0 = fetch_node(ctx, competitor.args[0]);
  const uuber::NativeNode &c1 = fetch_node(ctx, competitor.args[1]);
  if (t0.kind != "event" || t1.kind != "event" || c0.kind != "event" ||
      c1.kind != "event") {
    return false;
  }
  if (t0.source.empty() || t1.source.empty() || c0.source.empty() ||
      c1.source.empty()) {
    return false;
  }

  const std::string &target_a = t0.source;
  const std::string &target_b = t1.source;
  const std::string &comp_a = c0.source;
  const std::string &comp_b = c1.source;

  const uuber::NativeNode *shared_node = nullptr;
  const uuber::NativeNode *x_node = nullptr;
  const uuber::NativeNode *y_node = nullptr;
  if (target_a == comp_a) {
    shared_node = &t0;
    x_node = &t1;
    y_node = &c1;
  } else if (target_a == comp_b) {
    shared_node = &t0;
    x_node = &t1;
    y_node = &c0;
  } else if (target_b == comp_a) {
    shared_node = &t1;
    x_node = &t0;
    y_node = &c1;
  } else if (target_b == comp_b) {
    shared_node = &t1;
    x_node = &t0;
    y_node = &c0;
  } else {
    return false;
  }

  if (!shared_node || !x_node || !y_node)
    return false;
  const std::string &shared = shared_node->source;
  const std::string &x = x_node->source;
  const std::string &y = y_node->source;
  if (shared.empty() || x.empty() || y.empty())
    return false;
  if (x == y)
    return false;
  if (x == shared || y == shared)
    return false;

  out.x = x;
  out.y = y;
  out.c = shared;
  out.x_ref = x_node->source_ref;
  out.y_ref = y_node->source_ref;
  out.c_ref = shared_node->source_ref;
  return true;
}

bool detect_shared_gate_nway(const uuber::NativeContext &ctx, int node_id,
                             const std::vector<int> &competitor_node_ids,
                             SharedGateNWay &out) {
  if (node_id < 0)
    return false;
  if (competitor_node_ids.empty())
    return false;

  std::vector<int> node_ids;
  node_ids.reserve(1 + competitor_node_ids.size());
  node_ids.push_back(node_id);
  for (int cid : competitor_node_ids) {
    if (cid == NA_INTEGER || cid < 0)
      return false;
    node_ids.push_back(cid);
  }
  if (node_ids.size() < 3)
    return false; // N-way means >= 3 outcomes

  struct AndPair {
    std::string a;
    std::string b;
    int child0_id;
    int child1_id;
  };
  std::vector<AndPair> pairs;
  pairs.reserve(node_ids.size());
  std::unordered_map<std::string, LabelRef> label_ref_map;
  label_ref_map.reserve(node_ids.size() * 2);
  for (int id : node_ids) {
    const uuber::NativeNode &node = fetch_node(ctx, id);
    if (node.kind != "and")
      return false;
    if (node.args.size() != 2)
      return false;
    int c0_id = node.args[0];
    int c1_id = node.args[1];
    const uuber::NativeNode &c0 = fetch_node(ctx, c0_id);
    const uuber::NativeNode &c1 = fetch_node(ctx, c1_id);
    if (c0.kind != "event" || c1.kind != "event")
      return false;
    if (c0.source.empty() || c1.source.empty())
      return false;
    if (!c0.source.empty()) {
      label_ref_map.emplace(c0.source, c0.source_ref);
    }
    if (!c1.source.empty()) {
      label_ref_map.emplace(c1.source, c1.source_ref);
    }
    pairs.push_back(AndPair{c0.source, c1.source, c0_id, c1_id});
  }

  const std::string &cand1 = pairs[0].a;
  const std::string &cand2 = pairs[0].b;
  if (cand1.empty() || cand2.empty())
    return false;
  if (cand1 == cand2)
    return false;

  bool cand1_all = true;
  bool cand2_all = true;
  for (std::size_t i = 1; i < pairs.size(); ++i) {
    const AndPair &p = pairs[i];
    if (!(p.a == cand1 || p.b == cand1))
      cand1_all = false;
    if (!(p.a == cand2 || p.b == cand2))
      cand2_all = false;
  }
  if (cand1_all == cand2_all)
    return false; // either none shared or both shared
  const std::string shared_gate = cand1_all ? cand1 : cand2;

  std::unordered_set<std::string> seen_x;
  seen_x.reserve(pairs.size());
  std::string target_x;
  std::vector<std::string> competitor_xs;
  competitor_xs.reserve(pairs.size() - 1);

  std::vector<int> gate_source_ids;
  std::unordered_map<std::string, std::vector<int>> x_source_ids;
  x_source_ids.reserve(pairs.size());

  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const AndPair &p = pairs[i];
    const bool a_is_gate = (p.a == shared_gate);
    const bool b_is_gate = (p.b == shared_gate);
    if (a_is_gate == b_is_gate)
      return false;
    const std::string &x_lbl = a_is_gate ? p.b : p.a;
    if (x_lbl.empty() || x_lbl == shared_gate)
      return false;
    if (!seen_x.insert(x_lbl).second)
      return false;

    int gate_child_id = a_is_gate ? p.child0_id : p.child1_id;
    int x_child_id = a_is_gate ? p.child1_id : p.child0_id;
    const uuber::NativeNode &gate_child = fetch_node(ctx, gate_child_id);
    const uuber::NativeNode &x_child = fetch_node(ctx, x_child_id);
    if (gate_source_ids.empty()) {
      gate_source_ids = ensure_source_ids(ctx, gate_child);
    }
    x_source_ids.emplace(x_lbl, ensure_source_ids(ctx, x_child));

    if (i == 0) {
      target_x = x_lbl;
    } else {
      competitor_xs.push_back(x_lbl);
    }
  }

  if (shared_gate.empty() || target_x.empty() || competitor_xs.empty())
    return false;

  // Ensure the shared gate and all X leaves are independent (disjoint
  // underlying sources).
  std::vector<std::vector<int>> sets;
  sets.reserve(1 + x_source_ids.size());
  sets.push_back(gate_source_ids);
  for (const auto &kv : x_source_ids) {
    sets.push_back(kv.second);
  }
  for (std::size_t i = 0; i < sets.size(); ++i) {
    for (std::size_t j = i + 1; j < sets.size(); ++j) {
      if (share_sources(sets[i], sets[j])) {
        return false;
      }
    }
  }

  out.gate = shared_gate;
  out.target_x = std::move(target_x);
  out.competitor_xs = std::move(competitor_xs);
  out.competitor_refs.clear();
  out.competitor_refs.reserve(out.competitor_xs.size());
  auto ref_for_label = [&](const std::string &label) -> LabelRef {
    auto it = label_ref_map.find(label);
    if (it != label_ref_map.end()) {
      return it->second;
    }
    return resolve_label_ref_runtime(ctx, label);
  };
  out.gate_ref = ref_for_label(out.gate);
  out.target_ref = ref_for_label(out.target_x);
  for (const auto &lbl : out.competitor_xs) {
    out.competitor_refs.push_back(ref_for_label(lbl));
  }
  return true;
}

bool shared_gate_pair_cached(const uuber::NativeContext &ctx, int node_id,
                             int competitor_node_id, SharedGatePair &out) {
  if (node_id < 0 || competitor_node_id < 0)
    return false;
  static std::unordered_map<
      const uuber::NativeContext *,
      std::unordered_map<std::uint64_t, SharedGatePairCacheEntry>>
      cache_by_ctx;
  std::uint64_t key =
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(node_id)) << 32) ^
      static_cast<std::uint32_t>(competitor_node_id);
  SharedGatePairCacheEntry &entry = cache_by_ctx[&ctx][key];
  if (!entry.computed) {
    SharedGatePair spec;
    if (detect_shared_gate_pair(ctx, node_id, competitor_node_id, spec)) {
      entry.is_shared = true;
      entry.spec = std::move(spec);
    }
    entry.computed = true;
  }
  if (!entry.is_shared)
    return false;
  out = entry.spec;
  return true;
}

std::string competitor_cache_key(const std::vector<int> &ids);

bool shared_gate_nway_cached(const uuber::NativeContext &ctx, int node_id,
                             const std::vector<int> &competitor_node_ids,
                             SharedGateNWay &out) {
  if (node_id < 0)
    return false;
  if (competitor_node_ids.empty())
    return false;
  static std::unordered_map<
      const uuber::NativeContext *,
      std::unordered_map<std::string, SharedGateNWayCacheEntry>>
      cache_by_ctx;
  std::string key = std::to_string(node_id);
  key.push_back('|');
  key.append(competitor_cache_key(competitor_node_ids));
  SharedGateNWayCacheEntry &entry = cache_by_ctx[&ctx][key];
  if (!entry.computed) {
    SharedGateNWay spec;
    if (detect_shared_gate_nway(ctx, node_id, competitor_node_ids, spec)) {
      entry.is_shared = true;
      entry.spec = std::move(spec);
    }
    entry.computed = true;
  }
  if (!entry.is_shared)
    return false;
  out = entry.spec;
  return true;
}

NodeEvalResult eval_event_unforced(
    const uuber::NativeContext &ctx, const std::string &label,
    const LabelRef &label_ref, double t, const std::string &component,
    int component_idx, EvalNeed need, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context) {
  static const std::unordered_set<int> kEmptyForced;
  return eval_event_label(ctx, label, t, component, kEmptyForced, kEmptyForced,
                          need, trial_params, trial_type_key, include_na_donors,
                          outcome_label_context, &label_ref, component_idx,
                          outcome_idx_context);
}

double shared_gate_pair_density(const uuber::NativeContext &ctx,
                                const SharedGatePair &pair, double t,
                                const std::string &component, int component_idx,
                                const TrialParamSet *trial_params,
                                const std::string &trial_type_key,
                                bool include_na_donors,
                                const std::string &outcome_label_context,
                                int outcome_idx_context) {
  if (!std::isfinite(t) || t < 0.0)
    return 0.0;

  const EvalNeed need = EvalNeed::kDensity | EvalNeed::kSurvival;
  auto eval_unforced = [&](const std::string &lbl, const LabelRef &ref,
                           double u) -> NodeEvalResult {
    return eval_event_unforced(ctx, lbl, ref, u, component, component_idx, need,
                               trial_params, trial_type_key, include_na_donors,
                               outcome_label_context, outcome_idx_context);
  };
  NodeEvalResult ex = eval_unforced(pair.x, pair.x_ref, t);
  NodeEvalResult ey = eval_unforced(pair.y, pair.y_ref, t);
  NodeEvalResult ec = eval_unforced(pair.c, pair.c_ref, t);

  double fX = ex.density;
  double fC = ec.density;
  double SY = clamp_probability(ey.survival);
  double FY = clamp_probability(1.0 - SY);
  double SC = clamp_probability(ec.survival);
  double FC = clamp_probability(1.0 - SC);
  double SX = clamp_probability(ex.survival);
  double FX = clamp_probability(1.0 - SX);

  double term1 = (std::isfinite(fX) && fX > 0.0) ? (fX * FC * SY) : 0.0;
  double termCother = (std::isfinite(fC) && fC > 0.0) ? (fC * FX * SY) : 0.0;
  double term2 = 0.0;

  double denom = FX * FY;
  if (std::isfinite(fC) && fC > 0.0 && std::isfinite(denom) && denom > 0.0 &&
      t > 0.0) {
    auto integrand = [&](double u) -> double {
      if (!std::isfinite(u) || u < 0.0)
        return 0.0;
      NodeEvalResult ex_u = eval_unforced(pair.x, pair.x_ref, u);
      NodeEvalResult ey_u = eval_unforced(pair.y, pair.y_ref, u);
      double fx_u = ex_u.density;
      if (!std::isfinite(fx_u) || fx_u <= 0.0)
        return 0.0;
      double fy_u = clamp_probability(1.0 - clamp_probability(ey_u.survival));
      double val = fx_u * fy_u;
      if (!std::isfinite(val) || val <= 0.0)
        return 0.0;
      return val;
    };
    double integral = uuber::integrate_boost_fn(
        integrand, 0.0, t, kDefaultRelTol, kDefaultAbsTol, kDefaultMaxDepth);
    if (!std::isfinite(integral) || integral < 0.0)
      integral = 0.0;
    double term = denom - integral;
    if (term < 0.0)
      term = 0.0;
    term2 = fC * term;
  }

  double out = term1 + termCother + term2;
  if (!std::isfinite(out) || out < 0.0)
    return 0.0;
  return out;
}

double shared_gate_nway_density(const uuber::NativeContext &ctx,
                                const SharedGateNWay &gate, double t,
                                const std::string &component, int component_idx,
                                const TrialParamSet *trial_params,
                                const std::string &trial_type_key,
                                bool include_na_donors,
                                const std::string &outcome_label_context,
                                int outcome_idx_context) {
  if (!std::isfinite(t) || t < 0.0)
    return 0.0;
  if (gate.target_x.empty() || gate.gate.empty() || gate.competitor_xs.empty())
    return 0.0;
  const std::vector<LabelRef> *competitor_refs = &gate.competitor_refs;
  std::vector<LabelRef> resolved_competitor_refs;
  if (competitor_refs->size() != gate.competitor_xs.size()) {
    resolved_competitor_refs.resize(gate.competitor_xs.size());
    for (std::size_t i = 0; i < gate.competitor_xs.size(); ++i) {
      resolved_competitor_refs[i] =
          resolve_label_ref_runtime(ctx, gate.competitor_xs[i]);
    }
    competitor_refs = &resolved_competitor_refs;
  }

  const EvalNeed need = EvalNeed::kDensity | EvalNeed::kSurvival;
  auto eval_unforced = [&](const std::string &lbl, const LabelRef &ref,
                           double u) -> NodeEvalResult {
    return eval_event_unforced(ctx, lbl, ref, u, component, component_idx, need,
                               trial_params, trial_type_key, include_na_donors,
                               outcome_label_context, outcome_idx_context);
  };
  NodeEvalResult ex = eval_unforced(gate.target_x, gate.target_ref, t);
  NodeEvalResult ec = eval_unforced(gate.gate, gate.gate_ref, t);

  double prod_surv = 1.0;
  for (std::size_t i = 0; i < gate.competitor_xs.size(); ++i) {
    const std::string &lbl = gate.competitor_xs[i];
    NodeEvalResult ej = eval_unforced(lbl, (*competitor_refs)[i], t);
    double sj = clamp_probability(ej.survival);
    prod_surv *= sj;
    if (!std::isfinite(prod_surv) || prod_surv <= 0.0) {
      prod_surv = 0.0;
      break;
    }
  }

  double fX = ex.density;
  double fC = ec.density;
  double FC = clamp_probability(1.0 - clamp_probability(ec.survival));

  double term1 = (std::isfinite(fX) && fX > 0.0 && prod_surv > 0.0 && FC > 0.0)
                     ? (fX * FC * prod_surv)
                     : 0.0;

  double term2 = 0.0;
  if (std::isfinite(fC) && fC > 0.0 && t > 0.0) {
    auto integrand = [&](double u) -> double {
      if (!std::isfinite(u) || u < 0.0)
        return 0.0;
      NodeEvalResult ex_u = eval_unforced(gate.target_x, gate.target_ref, u);
      double fx_u = ex_u.density;
      if (!std::isfinite(fx_u) || fx_u <= 0.0)
        return 0.0;
      double prod = 1.0;
      for (std::size_t i = 0; i < gate.competitor_xs.size(); ++i) {
        const std::string &lbl = gate.competitor_xs[i];
        NodeEvalResult ej_u = eval_unforced(lbl, (*competitor_refs)[i], u);
        prod *= clamp_probability(ej_u.survival);
        if (!std::isfinite(prod) || prod <= 0.0)
          return 0.0;
      }
      double val = fx_u * prod;
      return (std::isfinite(val) && val > 0.0) ? val : 0.0;
    };
    double integral = uuber::integrate_boost_fn(
        integrand, 0.0, t, kDefaultRelTol, kDefaultAbsTol, kDefaultMaxDepth);
    if (!std::isfinite(integral) || integral < 0.0)
      integral = 0.0;
    term2 = fC * integral;
  }

  double out = term1 + term2;
  if (!std::isfinite(out) || out < 0.0)
    return 0.0;
  return out;
}

double shared_gate_nway_probability(
    const uuber::NativeContext &ctx, const SharedGateNWay &gate, double upper,
    const std::string &component, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    double rel_tol, double abs_tol, int max_depth, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context) {
  if (gate.target_x.empty() || gate.gate.empty() || gate.competitor_xs.empty())
    return 0.0;
  const std::vector<LabelRef> *competitor_refs = &gate.competitor_refs;
  std::vector<LabelRef> resolved_competitor_refs;
  if (competitor_refs->size() != gate.competitor_xs.size()) {
    resolved_competitor_refs.resize(gate.competitor_xs.size());
    for (std::size_t i = 0; i < gate.competitor_xs.size(); ++i) {
      resolved_competitor_refs[i] =
          resolve_label_ref_runtime(ctx, gate.competitor_xs[i]);
    }
    competitor_refs = &resolved_competitor_refs;
  }
  if (!std::isfinite(upper)) {
    upper = std::numeric_limits<double>::infinity();
  }
  if (upper <= 0.0)
    return 0.0;

  const EvalNeed need = EvalNeed::kDensity | EvalNeed::kSurvival;
  auto eval_unforced = [&](const std::string &lbl, const LabelRef &ref,
                           double t) -> NodeEvalResult {
    return eval_event_unforced(ctx, lbl, ref, t, component, component_idx, need,
                               trial_params, trial_type_key, include_na_donors,
                               outcome_label_context, outcome_idx_context);
  };

  // Gate must finish by `upper` for any outcome to be observed by `upper`.
  double gate_cdf = 0.0;
  if (std::isfinite(upper)) {
    NodeEvalResult ec_upper = eval_unforced(gate.gate, gate.gate_ref, upper);
    gate_cdf = clamp_probability(1.0 - clamp_probability(ec_upper.survival));
    if (!std::isfinite(gate_cdf) || gate_cdf <= 0.0)
      return 0.0;
  } else {
    // For upper = Inf, probability of finite gate time is 1 - S_gate(Inf).
    NodeEvalResult ec_inf = eval_unforced(
        gate.gate, gate.gate_ref, std::numeric_limits<double>::infinity());
    gate_cdf = clamp_probability(1.0 - clamp_probability(ec_inf.survival));
    if (gate_cdf <= 0.0)
      return 0.0;
  }

  auto order_integrand = [&](double t) -> double {
    NodeEvalResult ex = eval_unforced(gate.target_x, gate.target_ref, t);
    double fx = ex.density;
    if (!std::isfinite(fx) || fx <= 0.0)
      return 0.0;
    double prod = 1.0;
    for (std::size_t i = 0; i < gate.competitor_xs.size(); ++i) {
      const std::string &lbl = gate.competitor_xs[i];
      NodeEvalResult ej = eval_unforced(lbl, (*competitor_refs)[i], t);
      prod *= clamp_probability(ej.survival);
      if (!std::isfinite(prod) || prod <= 0.0)
        return 0.0;
    }
    double val = fx * prod;
    return (std::isfinite(val) && val > 0.0) ? val : 0.0;
  };

  double order_mass = 0.0;
  if (std::isfinite(upper)) {
    order_mass = uuber::integrate_boost_fn(order_integrand, 0.0, upper, rel_tol,
                                           abs_tol, max_depth);
  } else {
    auto transformed = [&](double x) -> double {
      if (x <= 0.0)
        return 0.0;
      if (x >= 1.0)
        x = std::nextafter(1.0, 0.0);
      double t = x / (1.0 - x);
      double jac = 1.0 / ((1.0 - x) * (1.0 - x));
      double val = order_integrand(t);
      if (!std::isfinite(val) || val <= 0.0)
        return 0.0;
      double out = val * jac;
      return std::isfinite(out) ? out : 0.0;
    };
    order_mass = uuber::integrate_boost_fn(transformed, 0.0, 1.0, rel_tol,
                                           abs_tol, max_depth);
  }
  if (!std::isfinite(order_mass) || order_mass < 0.0)
    order_mass = 0.0;

  double total = gate_cdf * order_mass;
  if (!std::isfinite(total) || total < 0.0)
    total = 0.0;
  return clamp_probability(total);
}

double shared_gate_pair_probability(
    const uuber::NativeContext &ctx, const SharedGatePair &pair, double upper,
    const std::string &component, int component_idx,
    const TrialParamSet *trial_params, const std::string &trial_type_key,
    double rel_tol, double abs_tol, int max_depth, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context) {
  if (!std::isfinite(upper)) {
    upper = std::numeric_limits<double>::infinity();
  }
  if (upper <= 0.0)
    return 0.0;

  auto integrate_signed = [&](auto &&integrand) -> double {
    if (std::isfinite(upper)) {
      return uuber::integrate_boost_fn(integrand, 0.0, upper, rel_tol, abs_tol,
                                       max_depth);
    }
    auto transformed = [&](double x) -> double {
      if (x <= 0.0)
        return 0.0;
      if (x >= 1.0)
        x = std::nextafter(1.0, 0.0);
      double t = x / (1.0 - x);
      double jac = 1.0 / ((1.0 - x) * (1.0 - x));
      double val = integrand(t);
      if (!std::isfinite(val))
        return 0.0;
      double out = val * jac;
      return std::isfinite(out) ? out : 0.0;
    };
    return uuber::integrate_boost_fn(transformed, 0.0, 1.0, rel_tol, abs_tol,
                                     max_depth);
  };

  const EvalNeed need = EvalNeed::kDensity | EvalNeed::kSurvival;
  auto eval_unforced = [&](const std::string &lbl, const LabelRef &ref,
                           double t) -> NodeEvalResult {
    return eval_event_unforced(ctx, lbl, ref, t, component, component_idx, need,
                               trial_params, trial_type_key, include_na_donors,
                               outcome_label_context, outcome_idx_context);
  };

  // Algebraic simplification of the previous 4-5 integral formulation:
  // total = ∫ [ fC(t)*FX(t) + fX(t)*(FC(t) - coeff*FY(t)) ] dt,
  // where coeff = FC(upper) for finite bounds and 1 for upper = Inf.
  double coeff = 1.0;
  if (std::isfinite(upper)) {
    NodeEvalResult ec_upper = eval_unforced(pair.c, pair.c_ref, upper);
    coeff = clamp_probability(1.0 - clamp_probability(ec_upper.survival));
  }

  auto total_integrand = [&](double t) -> double {
    NodeEvalResult ex = eval_unforced(pair.x, pair.x_ref, t);
    NodeEvalResult ec = eval_unforced(pair.c, pair.c_ref, t);

    double sx = clamp_probability(ex.survival);
    double Fx = clamp_probability(1.0 - sx);
    double fc_dens = ec.density;

    double val = 0.0;
    if (std::isfinite(fc_dens) && fc_dens > 0.0 && Fx > 0.0) {
      val += fc_dens * Fx;
    }

    double fx_dens = ex.density;
    if (std::isfinite(fx_dens) && fx_dens > 0.0) {
      NodeEvalResult ey = eval_unforced(pair.y, pair.y_ref, t);
      double Fy = clamp_probability(1.0 - clamp_probability(ey.survival));
      double Fc = clamp_probability(1.0 - clamp_probability(ec.survival));
      val += fx_dens * (Fc - coeff * Fy);
    }
    return std::isfinite(val) ? val : 0.0;
  };

  double total = integrate_signed(total_integrand);

  if (!std::isfinite(total) || total < 0.0)
    total = 0.0;
  return clamp_probability(total);
}

NodeEvalResult
eval_node_with_forced(const uuber::NativeContext &ctx, int node_id, double time,
                      const std::string &component, int component_idx,
                      const std::unordered_set<int> &forced_complete,
                      const std::unordered_set<int> &forced_survive,
                      EvalNeed need, const TrialParamSet *trial_params,
                      const std::string &trial_key,
                      const ExactSourceTimeMap *exact_source_times,
                      const SourceTimeBoundsMap *source_time_bounds,
                      const ForcedScopeFilter *forced_scope_filter) {
  NodeEvalState local(ctx, time, component, forced_complete, forced_survive,
                      trial_params, trial_key, false, std::string(),
                      component_idx, -1, exact_source_times,
                      source_time_bounds, forced_scope_filter);
  return eval_node_recursive(node_id, local, need);
}

NodeEvalResult eval_and_node(const uuber::NativeNode &node,
                             NodeEvalState &state, EvalNeed need) {
  const std::vector<int> &child_ids = node.args;
  if (child_ids.empty()) {
    return make_node_result(need, 0.0, 1.0, 0.0);
  }
  std::size_t n = child_ids.size();
  const bool want_density = needs_density(need);
  EvalNeed child_need = EvalNeed::kCDF;
  if (want_density) {
    child_need = child_need | EvalNeed::kDensity;
  }
  std::vector<NodeEvalResult> children;
  children.reserve(n);
  for (int child_id : child_ids) {
    children.push_back(eval_node_recursive(child_id, state, child_need));
  }
  std::vector<double> child_cdf(n);
  for (std::size_t i = 0; i < n; ++i) {
    child_cdf[i] = clamp_probability(children[i].cdf);
  }
  std::vector<double> prefix(n + 1, 1.0);
  std::vector<double> suffix(n + 1, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    prefix[i + 1] = prefix[i] * child_cdf[i];
    if (!std::isfinite(prefix[i + 1]))
      prefix[i + 1] = 0.0;
  }
  for (std::size_t i = n; i-- > 0;) {
    suffix[i] = suffix[i + 1] * child_cdf[i];
    if (!std::isfinite(suffix[i]))
      suffix[i] = 0.0;
  }
  double total_cdf = clamp_probability(prefix[n]);
  double density = std::numeric_limits<double>::quiet_NaN();
  if (want_density) {
    double accum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double child_density = children[i].density;
      if (!std::isfinite(child_density) || child_density <= 0.0)
        continue;
      double others = prefix[i] * suffix[i + 1];
      if (!std::isfinite(others) || others <= 0.0)
        continue;
      accum += child_density * clamp_probability(others);
    }
    density = (std::isfinite(accum) && accum >= 0.0) ? accum : 0.0;
  }
  double survival = needs_survival(need)
                        ? clamp_probability(1.0 - total_cdf)
                        : std::numeric_limits<double>::quiet_NaN();
  double cdf =
      needs_cdf(need) ? total_cdf : std::numeric_limits<double>::quiet_NaN();
  return make_node_result(need, density, survival, cdf);
}

NodeEvalResult eval_or_node(const uuber::NativeNode &node, NodeEvalState &state,
                            EvalNeed need) {
  const std::vector<int> &child_ids = node.args;
  if (child_ids.empty()) {
    return make_node_result(need, 0.0, 1.0, 0.0);
  }
  std::size_t n = child_ids.size();
  const bool want_density = needs_density(need);
  EvalNeed child_need = EvalNeed::kSurvival;
  if (want_density) {
    child_need = child_need | EvalNeed::kDensity;
  }
  std::vector<NodeEvalResult> children;
  children.reserve(n);
  for (int child_id : child_ids) {
    children.push_back(eval_node_recursive(child_id, state, child_need));
  }
  std::vector<double> child_surv(n);
  for (std::size_t i = 0; i < n; ++i) {
    child_surv[i] = clamp_probability(children[i].survival);
  }
  std::vector<double> prefix(n + 1, 1.0);
  std::vector<double> suffix(n + 1, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    prefix[i + 1] = prefix[i] * child_surv[i];
    if (!std::isfinite(prefix[i + 1]))
      prefix[i + 1] = 0.0;
  }
  for (std::size_t i = n; i-- > 0;) {
    suffix[i] = suffix[i + 1] * child_surv[i];
    if (!std::isfinite(suffix[i]))
      suffix[i] = 0.0;
  }
  double total_surv = clamp_probability(prefix[n]);
  double density = std::numeric_limits<double>::quiet_NaN();
  if (want_density) {
    double accum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double child_density = children[i].density;
      if (!std::isfinite(child_density) || child_density <= 0.0)
        continue;
      double others = prefix[i] * suffix[i + 1];
      if (!std::isfinite(others) || others <= 0.0)
        continue;
      accum += child_density * clamp_probability(others);
    }
    density = (std::isfinite(accum) && accum >= 0.0) ? accum : 0.0;
  }
  double survival = needs_survival(need)
                        ? total_surv
                        : std::numeric_limits<double>::quiet_NaN();
  double cdf = needs_cdf(need) ? clamp_probability(1.0 - total_surv)
                               : std::numeric_limits<double>::quiet_NaN();
  return make_node_result(need, density, survival, cdf);
}

NodeEvalResult eval_not_node(const uuber::NativeNode &node,
                             NodeEvalState &state, EvalNeed need) {
  int child_id = node.arg_id;
  if (child_id < 0) {
    return make_node_result(need, 0.0, 1.0, 0.0);
  }
  NodeEvalResult child = eval_node_recursive(child_id, state, EvalNeed::kCDF);
  double child_cdf = clamp_probability(child.cdf);
  double self_cdf = clamp_probability(1.0 - child_cdf);
  double survival = needs_survival(need)
                        ? child_cdf
                        : std::numeric_limits<double>::quiet_NaN();
  double cdf =
      needs_cdf(need) ? self_cdf : std::numeric_limits<double>::quiet_NaN();
  return make_node_result(need, 0.0, survival, cdf);
}

GuardEvalInput make_guard_input(const uuber::NativeContext &ctx,
                                const uuber::NativeNode &node,
                                const std::string &component, int component_idx,
                                const std::unordered_set<int> &forced_complete,
                                const std::unordered_set<int> &forced_survive,
                                const std::string &trial_key,
                                const TrialParamSet *trial_params,
                                const ExactSourceTimeMap *exact_source_times,
                                const SourceTimeBoundsMap *source_time_bounds,
                                const ForcedScopeFilter *forced_scope_filter) {
  GuardEvalInput input{ctx,
                       node,
                       component,
                       component_idx,
                       component_cache_key(component, trial_key),
                       &forced_complete,
                       &forced_survive,
                       nullptr,
                       {},
                       false,
                       trial_params,
                       exact_source_times,
                       source_time_bounds};
  input.local_scope_filter.source_ids = &node.source_ids;
  input.local_scope_filter.parent = forced_scope_filter;
  input.forced_scope_filter = &input.local_scope_filter;
  input.has_scoped_forced =
      forced_set_intersects_scope(forced_complete, input.forced_scope_filter) ||
      forced_set_intersects_scope(forced_survive, input.forced_scope_filter);
  return input;
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState &state,
                                   EvalNeed need) {
  const uuber::NativeNode &node = fetch_node(state.ctx, node_id);
  NodeEvalResult result;
  if (node.kind == "event") {
    result = eval_event_label(
        state.ctx, node.source, state.t, state.component, state.forced_complete,
        state.forced_survive, need, state.trial_params, state.trial_type_key,
        state.include_na_donors, state.outcome_label, &node.source_ref,
        state.component_idx, state.outcome_idx, state.exact_source_times,
        state.source_time_bounds, state.forced_scope_filter);
  } else if (node.kind == "and") {
    result = eval_and_node(node, state, need);
  } else if (node.kind == "or") {
    result = eval_or_node(node, state, need);
  } else if (node.kind == "not") {
    result = eval_not_node(node, state, need);
  } else if (node.kind == "guard") {
    GuardEvalInput guard_input =
        make_guard_input(state.ctx, node, state.component, state.component_idx,
                         state.forced_complete, state.forced_survive,
                         state.trial_type_key, state.trial_params,
                         state.exact_source_times, state.source_time_bounds,
                         state.forced_scope_filter);
    IntegrationSettings settings;
    double density = std::numeric_limits<double>::quiet_NaN();
    double survival = std::numeric_limits<double>::quiet_NaN();
    double cdf = std::numeric_limits<double>::quiet_NaN();
    if (needs_density(need)) {
      density = guard_density_internal(guard_input, state.t, settings);
    }
    if (needs_survival(need) || needs_cdf(need)) {
      cdf = guard_cdf_internal(guard_input, state.t, settings);
      if (needs_survival(need)) {
        survival = clamp_probability(1.0 - cdf);
      }
    }
    result = make_node_result(need, density, survival, cdf);
  } else {
    result = make_node_result(need, 0.0, 1.0, 0.0);
  }
  return result;
}

double guard_effective_survival_internal(const GuardEvalInput &input, double t,
                                         const IntegrationSettings &settings) {
  if (!std::isfinite(t)) {
    return 1.0;
  }
  if (t <= 0.0) {
    return 1.0;
  }
  int blocker_id = input.node.blocker_id;
  if (blocker_id < 0) {
    return 1.0;
  }
  NodeEvalResult block = eval_node_with_forced(
      input.ctx, blocker_id, t, input.component, input.component_idx,
      *input.forced_complete, *input.forced_survive, EvalNeed::kSurvival,
      input.trial_params, input.trial_type_key, input.exact_source_times,
      input.source_time_bounds, input.forced_scope_filter);
  return clamp_probability(block.survival);
}

double guard_reference_density(const GuardEvalInput &input, double t) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  int reference_id = input.node.reference_id;
  if (reference_id < 0) {
    return 0.0;
  }
  NodeEvalResult ref = eval_node_with_forced(
      input.ctx, reference_id, t, input.component, input.component_idx,
      *input.forced_complete, *input.forced_survive, EvalNeed::kDensity,
      input.trial_params, input.trial_type_key, input.exact_source_times,
      input.source_time_bounds, input.forced_scope_filter);
  double dens_ref = ref.density;
  if (!std::isfinite(dens_ref) || dens_ref <= 0.0)
    return 0.0;
  return dens_ref;
}

double guard_density_internal(const GuardEvalInput &input, double t,
                              const IntegrationSettings &settings) {
  double dens_ref = guard_reference_density(input, t);
  if (dens_ref <= 0.0)
    return 0.0;
  double s_eff = guard_effective_survival_internal(input, t, settings);
  double val = dens_ref * s_eff;
  if (!std::isfinite(val) || val <= 0.0)
    return 0.0;
  return val;
}

// --- Nested Guard Optimization ---

struct LinearGuardChain {
  std::vector<int> reference_ids; // [Ref0, Ref1, ...]
  int leaf_blocker_id{-1};
  bool valid{false};
};

struct FastEventInfo {
  const uuber::NativeAccumulator *acc{nullptr};
  const TrialAccumulatorParams *override{nullptr};
  int outcome_idx{-1};
  bool component_ok{false};
};

bool fast_event_info(const GuardEvalInput &input, int node_id,
                     FastEventInfo &info) {
  if (input.has_scoped_forced) {
    return false;
  }
  const uuber::NativeNode &node = fetch_node(input.ctx, node_id);
  if (node.kind != "event") {
    return false;
  }
  if (node.source.empty() || node.source == "__DEADLINE__" ||
      node.source == "__GUESS__") {
    return false;
  }
  int acc_idx = node.source_ref.acc_idx;
  if (acc_idx < 0 ||
      acc_idx >= static_cast<int>(input.ctx.accumulators.size())) {
    return false;
  }
  info.acc = &input.ctx.accumulators[acc_idx];
  info.override = get_trial_param_entry(input.trial_params, acc_idx);
  info.outcome_idx = node.source_ref.outcome_idx;
  info.component_ok =
      component_active(*info.acc, input.component, input.component_idx,
                       info.override);
  return true;
}

bool fast_event_density(const GuardEvalInput &input, int node_id, double t,
                        double &out) {
  FastEventInfo info;
  if (!fast_event_info(input, node_id, info)) {
    return false;
  }
  if (!info.component_ok) {
    out = 0.0;
    return true;
  }
  if (info.outcome_idx >= 0 &&
      info.outcome_idx < static_cast<int>(input.ctx.outcome_info.size())) {
    const uuber::OutcomeContextInfo &outcome =
        input.ctx.outcome_info[static_cast<std::size_t>(info.outcome_idx)];
    if (!outcome.alias_sources.empty() || !outcome.guess_donors.empty()) {
      return false;
    }
  }
  double onset = info.override ? info.override->onset : info.acc->onset;
  double q = info.override ? info.override->q : info.acc->q;
  const AccDistParams &cfg =
      info.override ? info.override->dist_cfg : info.acc->dist_cfg;
  out = safe_density(acc_density_from_cfg(t, onset, q, cfg));
  return true;
}

bool fast_event_cdf(const GuardEvalInput &input, int node_id, double t,
                    double &out) {
  FastEventInfo info;
  if (!fast_event_info(input, node_id, info)) {
    return false;
  }
  if (!info.component_ok) {
    out = 0.0;
    return true;
  }
  double onset = info.override ? info.override->onset : info.acc->onset;
  const AccDistParams &cfg =
      info.override ? info.override->dist_cfg : info.acc->dist_cfg;
  out = clamp_probability(acc_cdf_success_from_cfg(t, onset, cfg));
  return true;
}

bool fast_event_surv(const GuardEvalInput &input, int node_id, double t,
                     double &out) {
  FastEventInfo info;
  if (!fast_event_info(input, node_id, info)) {
    return false;
  }
  if (!info.component_ok) {
    out = 1.0;
    return true;
  }
  double onset = info.override ? info.override->onset : info.acc->onset;
  double q = info.override ? info.override->q : info.acc->q;
  const AccDistParams &cfg =
      info.override ? info.override->dist_cfg : info.acc->dist_cfg;
  out = clamp_probability(acc_survival_from_cfg(t, onset, q, cfg));
  return true;
}

// Detects A -> B -> C ... chain.
// Depth is number of GUARD nodes visited.
// A chain of depth 2 means: Guard0(A, Guard1(B, C)).
LinearGuardChain detect_linear_guard_chain(const uuber::NativeContext &ctx,
                                           const uuber::NativeNode &start_node,
                                           int max_depth_check) {
  LinearGuardChain chain;
  const uuber::NativeNode *curr = &start_node;
  int depth = 0;

  while (curr->kind == "guard" && depth < max_depth_check) {
    if (curr->reference_id < 0 || curr->blocker_id < 0) {
      return {}; // Invalid/Partial
    }
    chain.reference_ids.push_back(curr->reference_id);
    int next_id = curr->blocker_id;
    const uuber::NativeNode &next_node = fetch_node(ctx, next_id);

    // Move to next
    if (next_node.kind == "guard") {
      curr = &next_node;
      depth++;
    } else {
      // Leaf found
      chain.leaf_blocker_id = next_id;
      chain.valid = true;
      return chain;
    }
  }
  // If we exited because depth limit hit or kind mismatch (shouldn't happen
  // with leaf check) If curr is still guard, we exceeded depth or stopped? We
  // return invalid if we didn't find a non-guard leaf.
  return {};
}

// Helpers to evaluate components without creating full result structs
double quick_eval_density(const GuardEvalInput &input, int node_id, double t) {
  if (t < 0)
    return 0.0;
  double fast_val = 0.0;
  if (fast_event_density(input, node_id, t, fast_val)) {
    return fast_val;
  }
  NodeEvalResult res = eval_node_with_forced(
      input.ctx, node_id, t, input.component, input.component_idx,
      *input.forced_complete, *input.forced_survive, EvalNeed::kDensity,
      input.trial_params, input.trial_type_key, input.exact_source_times,
      input.source_time_bounds, input.forced_scope_filter);
  return res.density > 0 ? res.density : 0.0;
}

double quick_eval_cdf(const GuardEvalInput &input, int node_id, double t) {
  if (t <= 0)
    return 0.0;
  double fast_val = 0.0;
  if (fast_event_cdf(input, node_id, t, fast_val)) {
    return fast_val;
  }
  NodeEvalResult res = eval_node_with_forced(
      input.ctx, node_id, t, input.component, input.component_idx,
      *input.forced_complete, *input.forced_survive, EvalNeed::kCDF,
      input.trial_params, input.trial_type_key, input.exact_source_times,
      input.source_time_bounds, input.forced_scope_filter);
  return res.cdf; // clamp done in eval? usually yes
}

double quick_eval_surv(const GuardEvalInput &input, int node_id, double t) {
  if (t <= 0)
    return 1.0;
  double fast_val = 1.0;
  if (fast_event_surv(input, node_id, t, fast_val)) {
    return fast_val;
  }
  NodeEvalResult res = eval_node_with_forced(
      input.ctx, node_id, t, input.component, input.component_idx,
      *input.forced_complete, *input.forced_survive, EvalNeed::kSurvival,
      input.trial_params, input.trial_type_key, input.exact_source_times,
      input.source_time_bounds, input.forced_scope_filter);
  return res.survival;
}

// Depth 2 Optimization: Guard0(A, Guard1(B, C))
// P(t) = F_A(t) - Integral_0^t [ f_B(v) * S_C(v) * (F_A(t) - F_A(v)) ] dv
double eval_optimized_depth2(const GuardEvalInput &input, int id_A, int id_B,
                             int id_C, double t,
                             const IntegrationSettings &settings) {
  double FA_t = quick_eval_cdf(input, id_A, t);
  if (FA_t <= 0.0)
    return 0.0;

  auto integrand = [&](double v) -> double {
    double fB = quick_eval_density(input, id_B, v);
    if (fB <= 0)
      return 0.0;
    double SC = quick_eval_surv(input, id_C, v);
    if (SC <= 0)
      return 0.0;
    double FA_v = quick_eval_cdf(input, id_A, v);
    double term = FA_t - FA_v;
    if (term <= 0)
      return 0.0;
    return fB * SC * term;
  };

  double subtrahend =
      uuber::integrate_boost_fn(integrand, 0.0, t, settings.rel_tol,
                                settings.abs_tol, settings.max_depth);
  if (!std::isfinite(subtrahend))
    subtrahend = 0.0;

  return clamp_probability(FA_t - subtrahend);
}

// Depth 3 Optimization: Guard0(A, Guard1(B, Guard2(C, D)))
// P(t) = F_A(t) - T1 + T2
// T1 = Integral_0^t f_B(v) * (F_A(t) - F_A(v)) dv
// T2 = Integral_0^t f_C(w) * S_D(w) * K(w, t) dw
// K(w, t) = Integral_w^t f_B(v) * (F_A(t) - F_A(v)) dv
double eval_optimized_depth3(const GuardEvalInput &input, int id_A, int id_B,
                             int id_C, int id_D, double t,
                             const IntegrationSettings &settings) {
  double FA_t = quick_eval_cdf(input, id_A, t);
  if (FA_t <= 0.0)
    return 0.0;

  // T1: Effect of B if C were absent
  // T1 = Int_0^t f_B(v) * (FA(t) - FA(v))
  auto integrand_T1 = [&](double v) -> double {
    double fB = quick_eval_density(input, id_B, v);
    if (fB <= 0)
      return 0.0;
    double FA_v = quick_eval_cdf(input, id_A, v);
    double term = FA_t - FA_v;
    return (term > 0) ? fB * term : 0.0;
  };
  double T1 = uuber::integrate_boost_fn(integrand_T1, 0.0, t, settings.rel_tol,
                                        settings.abs_tol, settings.max_depth);

  // T2: Restoration by C (blocked by D)
  // T2 = Int_0^t f_C(w) * S_D(w) * K(w, t) dw
  auto integrand_T2 = [&](double w) -> double {
    double fC = quick_eval_density(input, id_C, w);
    if (fC <= 0)
      return 0.0;
    double SD = quick_eval_surv(input, id_D, w);
    if (SD <= 0)
      return 0.0;

    // Inner Integral K(w, t): same integrand as T1, but from w to t
    double K_val = uuber::integrate_boost_fn(integrand_T1, w, t,
                                             settings.rel_tol,
                                             settings.abs_tol,
                                             settings.max_depth);
    return fC * SD * K_val;
  };

  double T2 = uuber::integrate_boost_fn(integrand_T2, 0.0, t, settings.rel_tol,
                                        settings.abs_tol, settings.max_depth);

  if (!std::isfinite(T1))
    T1 = 0.0;
  if (!std::isfinite(T2))
    T2 = 0.0;

  return clamp_probability(FA_t - T1 + T2);
}

double guard_cdf_internal(const GuardEvalInput &input, double t,
                          const IntegrationSettings &settings) {
  if (!std::isfinite(t)) {
    return 1.0;
  }
  if (t <= 0.0) {
    return 0.0;
  }

  // Attempt optimization for Depth 2 and 3
  // Check depth up to 4 to warn
  LinearGuardChain chain = detect_linear_guard_chain(input.ctx, input.node, 10);

  // Chain.reference_ids has [RefA, RefB, RefC...]
  // Chain.leaf_blocker_id is Leaf
  // Depth is size of reference_ids.
  // Depth 2: size 2. (A, B). Leaf C.
  // Depth 3: size 3. (A, B, C). Leaf D.

  size_t depth = chain.reference_ids.size();

  if (chain.valid) {
    if (depth == 2) {
      // A blocked by B blocked by C
      int id_A = chain.reference_ids[0];
      int id_B = chain.reference_ids[1];
      int id_C = chain.leaf_blocker_id;
      return eval_optimized_depth2(input, id_A, id_B, id_C, t, settings);
    } else if (depth == 3) {
      // A->B->C->D
      int id_A = chain.reference_ids[0];
      int id_B = chain.reference_ids[1];
      int id_C = chain.reference_ids[2];
      int id_D = chain.leaf_blocker_id;
      return eval_optimized_depth3(input, id_A, id_B, id_C, id_D, t, settings);
    } else if (depth > 3) {
      // Fallback with warning
      // TODO: Implement ODE solver for O(N) scaling at depth > 3.
      // Current naive recursion is O(N^k).
      Rcpp::warning("Nested inhibition depth %d detected. Optimization "
                    "available only up to depth 3. Using slow recursive "
                    "integration. Implement ODE solver for O(N) scalability.",
                    depth);
    }
  }

  // Fallback to naive recursive integration

  auto integrand = [&](double u) -> double {
    return guard_density_internal(input, u, settings);
  };
  double val = uuber::integrate_boost_fn(integrand, 0.0, t, settings.rel_tol,
                                         settings.abs_tol, settings.max_depth);
  if (!std::isfinite(val) || val < 0.0)
    val = 0.0;
  return clamp_probability(val);
}

std::vector<int> gather_blocker_sources(const uuber::NativeContext &ctx,
                                        const uuber::NativeNode &guard_node) {
  std::vector<int> sources;
  const uuber::NativeNode &blocker = fetch_node(ctx, guard_node.blocker_id);
  sources = blocker.source_ids;
  if (blocker.kind == "event" && !blocker.source.empty()) {
    int label_idx = blocker.source_ref.label_id;
    if (label_idx < 0) {
      label_idx = NA_INTEGER;
    }
    sources.push_back(label_idx);
  }
  sort_unique(sources);
  return sources;
}

struct CompetitorMeta {
  int node_id;
  const uuber::NativeNode *node;
  std::vector<int> sources;
  bool has_guard{false};
  bool scenario_sensitive{false};
};

struct CompetitorClusterCacheEntry {
  std::vector<CompetitorMeta> metas;
  std::vector<std::vector<int>> clusters;
  std::vector<uint8_t> cluster_has_guard;
  std::vector<std::vector<int>> guard_orders;
};

using CompetitorCacheMap =
    std::unordered_map<std::string, CompetitorClusterCacheEntry>;

static std::unordered_map<const uuber::NativeContext *, CompetitorCacheMap>
    g_competitor_cache;

std::string competitor_cache_key(const std::vector<int> &ids) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (i > 0)
      oss << ",";
    oss << ids[i];
  }
  return oss.str();
}

// Forward declaration
std::vector<std::vector<int>>
build_competitor_clusters(const std::vector<CompetitorMeta> &metas);

bool node_contains_guard(int node_id, const uuber::NativeContext &ctx,
                         std::unordered_map<int, bool> &memo) {
  auto it = memo.find(node_id);
  if (it != memo.end())
    return it->second;
  const uuber::NativeNode &node = fetch_node(ctx, node_id);
  bool result = (node.kind == "guard");
  if (!result) {
    if (node.kind == "and" || node.kind == "or") {
      for (int child_id : node.args) {
        if (node_contains_guard(child_id, ctx, memo)) {
          result = true;
          break;
        }
      }
    } else if (node.kind == "not" && node.arg_id >= 0) {
      result = node_contains_guard(node.arg_id, ctx, memo);
    }
  }
  memo[node_id] = result;
  return result;
}

const CompetitorClusterCacheEntry &
fetch_competitor_cluster_cache(const uuber::NativeContext &ctx,
                               const std::vector<int> &competitor_ids) {
  CompetitorCacheMap &cache = g_competitor_cache[&ctx];
  const std::string key = competitor_cache_key(competitor_ids);
  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }

  CompetitorClusterCacheEntry entry;
  entry.metas.reserve(competitor_ids.size());
  std::unordered_map<int, bool> guard_memo;
  for (int node_id : competitor_ids) {
    const uuber::NativeNode &node = fetch_node(ctx, node_id);
    CompetitorMeta meta;
    meta.node_id = node_id;
    meta.node = &node;
    meta.sources = node.source_ids;
    sort_unique(meta.sources);
    meta.has_guard = node_contains_guard(node_id, ctx, guard_memo);
    meta.scenario_sensitive = node.scenario_sensitive;
    entry.metas.push_back(std::move(meta));
  }

  entry.clusters = build_competitor_clusters(entry.metas);
  entry.cluster_has_guard.reserve(entry.clusters.size());
  entry.guard_orders.resize(entry.clusters.size());
  for (std::size_t i = 0; i < entry.clusters.size(); ++i) {
    const auto &cluster = entry.clusters[i];
    bool has_guard = false;
    for (int idx : cluster) {
      if (entry.metas[static_cast<std::size_t>(idx)].has_guard) {
        has_guard = true;
        break;
      }
    }
    entry.cluster_has_guard.push_back(has_guard ? 1 : 0);
    if (has_guard) {
      struct OrderingMeta {
        int index;
        std::size_t source_count;
        bool scenario_sensitive;
      };
      std::vector<OrderingMeta> order;
      order.reserve(cluster.size());
      for (int idx : cluster) {
        const CompetitorMeta &meta = entry.metas[static_cast<std::size_t>(idx)];
        order.push_back({idx, meta.sources.size(), meta.scenario_sensitive});
      }
      std::sort(order.begin(), order.end(),
                [&](const OrderingMeta &a, const OrderingMeta &b) {
                  if (a.scenario_sensitive != b.scenario_sensitive) {
                    return !a.scenario_sensitive && b.scenario_sensitive;
                  }
                  if (a.source_count != b.source_count) {
                    return a.source_count < b.source_count;
                  }
                  return a.index < b.index;
                });
      std::vector<int> guard_order;
      guard_order.reserve(order.size());
      for (const auto &rec : order) {
        guard_order.push_back(rec.index);
      }
      entry.guard_orders[i] = std::move(guard_order);
    }
  }

  auto inserted = cache.emplace(key, std::move(entry));
  return inserted.first->second;
}

bool share_sources(const std::vector<int> &a, const std::vector<int> &b) {
  if (a.empty() || b.empty())
    return false;
  std::size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j])
      return true;
    if (a[i] < b[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return false;
}

std::vector<std::vector<int>>
build_competitor_clusters(const std::vector<CompetitorMeta> &metas) {
  std::vector<std::vector<int>> clusters;
  const std::size_t n = metas.size();
  if (n == 0)
    return clusters;
  std::vector<bool> visited(n, false);
  for (std::size_t i = 0; i < n; ++i) {
    if (visited[i])
      continue;
    std::vector<int> cluster;
    std::queue<std::size_t> q;
    q.push(i);
    visited[i] = true;
    while (!q.empty()) {
      std::size_t idx = q.front();
      q.pop();
      cluster.push_back(static_cast<int>(idx));
      for (std::size_t j = 0; j < n; ++j) {
        if (visited[j])
          continue;
        if (share_sources(metas[idx].sources, metas[j].sources)) {
          visited[j] = true;
          q.push(j);
        }
      }
    }
    clusters.push_back(cluster);
  }
  return clusters;
}

uuber::PoolTemplateCacheEntry
build_pool_template_cache(int n, const std::vector<int> &member_ids,
                          int pool_idx, int k) {
  uuber::PoolTemplateCacheEntry cache;
  if (n <= 0 || k < 1 || k > n) {
    return cache;
  }
  cache.finisher_map.assign(static_cast<std::size_t>(n), std::vector<int>());
  const int need = k - 1;
  int template_counter = 0;
  const bool pool_idx_valid = pool_idx != NA_INTEGER;

  for (int idx = 0; idx < n; ++idx) {
    std::vector<int> others;
    others.reserve(n - 1);
    for (int j = 0; j < n; ++j) {
      if (j == idx)
        continue;
      others.push_back(j + 1);
    }
    if (need > static_cast<int>(others.size())) {
      continue;
    }
    std::vector<std::vector<int>> combos = generate_combinations(others, need);
    for (const auto &combo : combos) {
      std::vector<int> survivors = survivors_from_combo(others, combo);

      std::vector<int> forced_complete_ids;
      int finisher_member_id = member_ids[idx];
      if (finisher_member_id != NA_INTEGER) {
        forced_complete_ids.push_back(finisher_member_id);
      }
      for (int c : combo) {
        int comp_id = member_ids[c - 1];
        if (comp_id != NA_INTEGER) {
          forced_complete_ids.push_back(comp_id);
        }
      }
      if (pool_idx_valid) {
        forced_complete_ids.push_back(pool_idx);
      }
      sort_unique(forced_complete_ids);

      std::vector<int> forced_survive_ids;
      for (int s : survivors) {
        int surv_id = member_ids[s - 1];
        if (surv_id != NA_INTEGER) {
          forced_survive_ids.push_back(surv_id);
        }
      }
      sort_unique(forced_survive_ids);

      std::vector<int> complete_zero;
      complete_zero.reserve(combo.size());
      for (int c : combo) {
        complete_zero.push_back(c - 1);
      }
      std::vector<int> survivor_zero;
      survivor_zero.reserve(survivors.size());
      for (int s : survivors) {
        survivor_zero.push_back(s - 1);
      }

      uuber::PoolTemplateEntry entry;
      entry.finisher_idx = idx;
      entry.complete_idx = std::move(complete_zero);
      entry.survivor_idx = std::move(survivor_zero);
      entry.forced_complete_ids = forced_complete_ids;
      entry.forced_survive_ids = forced_survive_ids;

      cache.templates.push_back(std::move(entry));
      cache.finisher_map[static_cast<std::size_t>(idx)].push_back(
          template_counter);
      ++template_counter;
    }
  }
  return cache;
}

double
evaluate_survival_with_forced(int node_id,
                              const std::unordered_set<int> &forced_complete,
                              const std::unordered_set<int> &forced_survive,
                              const std::string &component, int component_idx,
                              double t, const uuber::NativeContext &ctx,
                              const std::string &trial_key = std::string(),
                              const TrialParamSet *trial_params = nullptr,
                              const ExactSourceTimeMap *exact_source_times =
                                  nullptr,
                              const SourceTimeBoundsMap *source_time_bounds =
                                  nullptr) {
  NodeEvalState state(ctx, t, component, forced_complete, forced_survive,
                      trial_params, trial_key, false, std::string(),
                      component_idx, -1, exact_source_times,
                      source_time_bounds);
  NodeEvalResult res = eval_node_recursive(node_id, state, EvalNeed::kSurvival);
  return clamp_probability(res.survival);
}

double compute_guard_free_cluster_value(
    const std::vector<int> &cluster_indices,
    const std::vector<CompetitorMeta> &metas, const uuber::NativeContext &ctx,
    double t, const std::string &component, int component_idx,
    const std::unordered_set<int> &forced_complete,
    const std::unordered_set<int> &forced_survive,
    const std::string &trial_key = std::string(),
    const TrialParamSet *trial_params = nullptr,
    const ExactSourceTimeMap *exact_source_times = nullptr,
    const SourceTimeBoundsMap *source_time_bounds = nullptr) {
  double prod = 1.0;
  for (int idx : cluster_indices) {
    double surv = evaluate_survival_with_forced(
        metas[static_cast<std::size_t>(idx)].node_id, forced_complete,
        forced_survive, component, component_idx, t, ctx, trial_key,
        trial_params, exact_source_times, source_time_bounds);
    if (!std::isfinite(surv) || surv <= 0.0)
      return 0.0;
    prod *= surv;
    if (!std::isfinite(prod) || prod <= 0.0)
      return 0.0;
  }
  return clamp_probability(prod);
}

double compute_guard_cluster_value(
    const std::vector<int> &cluster_indices,
    const std::vector<CompetitorMeta> &metas, const uuber::NativeContext &ctx,
    double t, const std::string &component, int component_idx,
    const std::unordered_set<int> &base_forced_complete,
    const std::unordered_set<int> &base_forced_survive,
    const std::string &trial_key = std::string(),
    const TrialParamSet *trial_params = nullptr,
    const std::vector<int> *guard_order = nullptr,
    std::unordered_set<int> *forced_survive_scratch = nullptr,
    const ExactSourceTimeMap *exact_source_times = nullptr,
    const SourceTimeBoundsMap *source_time_bounds = nullptr) {
  const std::vector<int> &order = guard_order ? *guard_order : cluster_indices;
  std::unordered_set<int> local_forced;
  std::unordered_set<int> &forced_survive =
      forced_survive_scratch ? *forced_survive_scratch : local_forced;
  forced_survive = base_forced_survive;
  double prod = 1.0;
  for (int idx : order) {
    const CompetitorMeta &node_meta = metas[static_cast<std::size_t>(idx)];
    double surv = evaluate_survival_with_forced(
        node_meta.node_id, base_forced_complete, forced_survive, component,
        component_idx, t, ctx, trial_key, trial_params, exact_source_times,
        source_time_bounds);
    if (!std::isfinite(surv) || surv <= 0.0)
      return 0.0;
    prod *= surv;
    if (!std::isfinite(prod) || prod <= 0.0)
      return 0.0;
    for (int src : node_meta.sources) {
      forced_survive.insert(src);
    }
  }
  return clamp_probability(prod);
}

double competitor_survival_internal(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    double t, const std::string &component_label, int component_idx,
    const std::unordered_set<int> &forced_complete,
    const std::unordered_set<int> &forced_survive,
    const std::string &trial_type_key = std::string(),
    const TrialParamSet *trial_params = nullptr,
    const ExactSourceTimeMap *exact_source_times = nullptr,
    const SourceTimeBoundsMap *source_time_bounds = nullptr) {
  if (competitor_ids.empty())
    return 1.0;
  const CompetitorClusterCacheEntry &cache =
      fetch_competitor_cluster_cache(ctx, competitor_ids);
  if (cache.clusters.empty())
    return 1.0;

  double product = 1.0;
  std::unordered_set<int> forced_survive_scratch;
  for (std::size_t i = 0; i < cache.clusters.size(); ++i) {
    const auto &cluster = cache.clusters[i];
    const bool cluster_has_guard =
        (i < cache.cluster_has_guard.size()) && cache.cluster_has_guard[i];
    double cluster_val =
        cluster_has_guard
            ? compute_guard_cluster_value(
                  cluster, cache.metas, ctx, t, component_label, component_idx,
                  forced_complete, forced_survive, trial_type_key, trial_params,
                  (i < cache.guard_orders.size() ? &cache.guard_orders[i]
                                                 : nullptr),
                  &forced_survive_scratch, exact_source_times,
                  source_time_bounds)
            : compute_guard_free_cluster_value(cluster, cache.metas, ctx, t,
                                               component_label, component_idx,
                                               forced_complete, forced_survive,
                                               trial_type_key, trial_params,
                                               exact_source_times,
                                               source_time_bounds);
    if (!std::isfinite(cluster_val) || cluster_val <= 0.0) {
      return 0.0;
    }
    product *= cluster_val;
    if (!std::isfinite(product) || product <= 0.0) {
      return 0.0;
    }
  }
  return clamp_probability(product);
}

double node_density_with_competitors_internal(
    const uuber::NativeContext &ctx, int node_id, double t,
    const std::string &component_label, int component_idx,
    const std::unordered_set<int> &forced_complete,
    const std::unordered_set<int> &forced_survive,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context,
    const ExactSourceTimeMap *exact_source_times,
    const SourceTimeBoundsMap *source_time_bounds) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  if (competitor_ids.size() >= 2 && forced_complete.empty() &&
      forced_survive.empty()) {
    SharedGateNWay shared_gate;
    if (shared_gate_nway_cached(ctx, node_id, competitor_ids, shared_gate)) {
      return shared_gate_nway_density(
          ctx, shared_gate, t, component_label, component_idx, trial_params,
          trial_type_key, include_na_donors, outcome_label_context,
          outcome_idx_context);
    }
  }
  if (competitor_ids.size() == 1 && forced_complete.empty() &&
      forced_survive.empty() && competitor_ids[0] != NA_INTEGER) {
    SharedGatePair shared_gate;
    if (shared_gate_pair_cached(ctx, node_id, competitor_ids[0], shared_gate)) {
      return shared_gate_pair_density(
          ctx, shared_gate, t, component_label, component_idx, trial_params,
          trial_type_key, include_na_donors, outcome_label_context,
          outcome_idx_context);
    }
  }
  NodeEvalState state(ctx, t, component_label, forced_complete, forced_survive,
                      trial_params, trial_type_key, include_na_donors,
                      outcome_label_context, component_idx,
                      outcome_idx_context, exact_source_times,
                      source_time_bounds);
  NodeEvalResult base = eval_node_recursive(node_id, state, EvalNeed::kDensity);
  double density = base.density;
  if (!std::isfinite(density) || density <= 0.0) {
    return 0.0;
  }
  if (!competitor_ids.empty()) {
    double surv = competitor_survival_internal(
        ctx, competitor_ids, t, component_label, component_idx,
        forced_complete, forced_survive,
        state.trial_type_key, trial_params, exact_source_times,
        source_time_bounds);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    density *= surv;
  }
  return density;
}

inline double node_density_with_shared_triggers_plan(
    const uuber::NativeContext &ctx, int node_id, double t,
    const std::string &component_label, int component_idx,
    const std::unordered_set<int> &forced_complete,
    const std::unordered_set<int> &forced_survive,
    const std::vector<int> &competitor_ids, const TrialParamSet *trial_params,
    const std::string &trial_type_key, bool include_na_donors,
    const std::string &outcome_label_context, int outcome_idx_context,
    const SharedTriggerPlan *trigger_plan, TrialParamSet *scratch_params) {
  if (!trial_params) {
    return node_density_with_competitors_internal(
        ctx, node_id, t, component_label, component_idx, forced_complete,
        forced_survive, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_label_context, outcome_idx_context);
  }
  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_plan = build_shared_trigger_plan(ctx, trial_params);
    plan_ptr = &local_plan;
  }
  if (plan_ptr->triggers.empty()) {
    return node_density_with_competitors_internal(
        ctx, node_id, t, component_label, component_idx, forced_complete,
        forced_survive, competitor_ids, trial_params, trial_type_key,
        include_na_donors, outcome_label_context, outcome_idx_context);
  }
  const std::vector<SharedTriggerInfo> &triggers = plan_ptr->triggers;
  if (triggers.size() >= 63) {
    Rcpp::stop("Too many shared triggers for density evaluation");
  }
  TrialParamSet local_scratch;
  TrialParamSet &scratch = scratch_params ? *scratch_params : local_scratch;
  scratch = *trial_params;
  for (const auto &trigger : triggers) {
    apply_trigger_state_inplace(scratch, trigger, false);
  }
  const std::uint64_t n_states = 1ULL << triggers.size();
  double total = 0.0;
  std::uint64_t mask = 0;
  for (std::uint64_t idx = 0; idx < n_states; ++idx) {
    double w = plan_ptr->mask_weights.empty()
                   ? shared_trigger_mask_weight(triggers, mask)
                   : plan_ptr->mask_weights[mask];
    if (w > 0.0) {
      double d = node_density_with_competitors_internal(
          ctx, node_id, t, component_label, component_idx, forced_complete,
          forced_survive, competitor_ids, &scratch, trial_type_key,
          include_na_donors, outcome_label_context, outcome_idx_context);
      if (std::isfinite(d) && d > 0.0) {
        total += w * d;
      }
    }
    if (idx + 1 == n_states)
      break;
    std::uint64_t next = idx + 1;
    std::uint64_t next_gray = next ^ (next >> 1);
    std::uint64_t diff = next_gray ^ mask;
    if (diff != 0) {
      int bit = static_cast<int>(__builtin_ctzll(diff));
      bool fail = ((next_gray >> bit) & 1ULL) != 0ULL;
      apply_trigger_state_inplace(
          scratch, triggers[static_cast<std::size_t>(bit)], fail);
    }
    mask = next_gray;
  }
  if (!std::isfinite(total) || total <= 0.0)
    return 0.0;
  return total;
}

double native_outcome_probability_impl(
    SEXP ctxSEXP, int node_id, double upper,
    Rcpp::Nullable<Rcpp::String> component, SEXP forced_complete,
    SEXP forced_survive, const Rcpp::IntegerVector &competitor_ids,
    double rel_tol, double abs_tol, int max_depth,
    const TrialParamSet *trial_params,
    const std::string &trial_type_key = std::string(),
    bool include_na_donors = false,
    const std::string &outcome_label_context = std::string()) {
  if (upper <= 0.0) {
    return 0.0;
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label =
      component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  int component_idx = component_index_of(*ctx, component_label);
  int outcome_idx_context =
      outcome_label_context.empty()
          ? -1
          : resolve_outcome_index(*ctx, outcome_label_context, component_label);
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  std::vector<int> comp_vec_raw = integer_vector_to_std(competitor_ids, false);
  std::vector<int> comp_vec_filtered;
  const std::vector<int> &comp_vec =
      filter_competitor_ids(*ctx, comp_vec_raw, component_label,
                            comp_vec_filtered);
  SharedGatePair shared_gate;
  SharedGateNWay shared_gate_nway;
  bool use_shared_gate = false;
  bool use_shared_gate_nway = false;
  if (comp_vec.size() == 1 && comp_vec[0] != NA_INTEGER &&
      forced_complete_set.empty() && forced_survive_set.empty()) {
    use_shared_gate =
        detect_shared_gate_pair(*ctx, node_id, comp_vec[0], shared_gate);
  } else if (comp_vec.size() >= 2 && forced_complete_set.empty() &&
             forced_survive_set.empty()) {
    use_shared_gate_nway =
        detect_shared_gate_nway(*ctx, node_id, comp_vec, shared_gate_nway);
  }
  auto integrand = [&](double u, const TrialParamSet *params_ptr) -> double {
    if (!std::isfinite(u) || u < 0.0)
      return 0.0;
    double val = node_density_with_competitors_internal(
        *ctx, node_id, u, component_label, component_idx, forced_complete_set,
        forced_survive_set, comp_vec, params_ptr, trial_type_key,
        include_na_donors, outcome_label_context, outcome_idx_context);
    if (!std::isfinite(val) || val <= 0.0)
      return 0.0;
    return val;
  };
  auto integrate_with_params = [&](const TrialParamSet *params_ptr) -> double {
    if (use_shared_gate) {
      return shared_gate_pair_probability(
          *ctx, shared_gate, upper, component_label, component_idx, params_ptr,
          trial_type_key, rel_tol, abs_tol, max_depth, include_na_donors,
          outcome_label_context, outcome_idx_context);
    }
    if (use_shared_gate_nway) {
      return shared_gate_nway_probability(
          *ctx, shared_gate_nway, upper, component_label, component_idx,
          params_ptr, trial_type_key, rel_tol, abs_tol, max_depth,
          include_na_donors, outcome_label_context, outcome_idx_context);
    }
    double integral = 0.0;
    if (std::isfinite(upper)) {
      integral = uuber::integrate_boost_fn(
          [&](double u) { return integrand(u, params_ptr); }, 0.0, upper,
          rel_tol, abs_tol, max_depth);
    } else {
      auto transformed = [&](double x) -> double {
        if (x <= 0.0)
          return 0.0;
        if (x >= 1.0)
          x = std::nextafter(1.0, 0.0);
        double t = x / (1.0 - x);
        double jac = 1.0 / ((1.0 - x) * (1.0 - x));
        double val = integrand(t, params_ptr);
        if (!std::isfinite(val) || val <= 0.0)
          return 0.0;
        double out = val * jac;
        return std::isfinite(out) ? out : 0.0;
      };
      integral = uuber::integrate_boost_fn(transformed, 0.0, 1.0, rel_tol,
                                           abs_tol, max_depth);
    }
    if (!std::isfinite(integral))
      integral = 0.0;
    return clamp_probability(integral);
  };

  // Ensure we have a concrete parameter set to modify for shared triggers.
  TrialParamSet base_params_holder;
  const TrialParamSet *params_ptr = trial_params;
  if (!params_ptr) {
    base_params_holder = build_base_paramset(*ctx);
    params_ptr = &base_params_holder;
  }

  // Enumerate shared-trigger states to preserve correlation across linked
  // accumulators.
  std::vector<SharedTriggerInfo> triggers =
      collect_shared_triggers(*ctx, *params_ptr);
  if (triggers.empty()) {
    return integrate_with_params(params_ptr);
  }

  if (triggers.size() >= 63) {
    Rcpp::stop("Too many shared triggers for outcome probability evaluation");
  }
  const std::uint64_t n_states = 1ULL << triggers.size();
  std::vector<double> mask_weights;
  constexpr std::uint64_t kMaxPrecomputedWeights = 1ULL << 20;
  if (n_states <= kMaxPrecomputedWeights) {
    mask_weights.resize(static_cast<std::size_t>(n_states));
    for (std::uint64_t m = 0; m < n_states; ++m) {
      mask_weights[static_cast<std::size_t>(m)] =
          shared_trigger_mask_weight(triggers, m);
    }
  }
  TrialParamSet scratch = *params_ptr;
  for (const auto &trigger : triggers) {
    apply_trigger_state_inplace(scratch, trigger, false);
  }
  double total = 0.0;
  std::uint64_t mask = 0;
  for (std::uint64_t idx = 0; idx < n_states; ++idx) {
    double w = mask_weights.empty()
                   ? shared_trigger_mask_weight(triggers, mask)
                   : mask_weights[static_cast<std::size_t>(mask)];
    if (w > 0.0) {
      double p = integrate_with_params(&scratch);
      if (std::isfinite(p) && p > 0.0) {
        total += w * p;
      }
    }
    if (idx + 1 == n_states) {
      break;
    }
    std::uint64_t next = idx + 1;
    std::uint64_t next_gray = next ^ (next >> 1);
    std::uint64_t diff = next_gray ^ mask;
    if (diff != 0) {
      int bit = static_cast<int>(__builtin_ctzll(diff));
      bool fail = ((next_gray >> bit) & 1ULL) != 0ULL;
      apply_trigger_state_inplace(
          scratch, triggers[static_cast<std::size_t>(bit)], fail);
    }
    mask = next_gray;
  }
  return clamp_probability(total);
}

std::vector<std::string> string_vector_from_entry(SEXP entry) {
  if (Rf_isNull(entry))
    return {};
  Rcpp::CharacterVector vec(entry);
  std::vector<std::string> out;
  out.reserve(vec.size());
  for (R_xlen_t i = 0; i < vec.size(); ++i) {
    if (vec[i] == NA_STRING)
      continue;
    out.push_back(Rcpp::as<std::string>(vec[i]));
  }
  return out;
}

std::vector<int>
component_indices_from_labels(const uuber::NativeContext &ctx,
                              const std::vector<std::string> &labels) {
  std::vector<int> out;
  out.reserve(labels.size());
  for (const auto &label : labels) {
    auto it = ctx.component_index.find(label);
    if (it != ctx.component_index.end()) {
      out.push_back(it->second);
    }
  }
  return out;
}

std::unique_ptr<TrialParamSet>
build_trial_params_from_df(const uuber::NativeContext &ctx,
                           const Rcpp::Nullable<Rcpp::DataFrame> &rows_opt) {
  if (rows_opt.isNull())
    return nullptr;
  Rcpp::DataFrame rows(rows_opt.get());
  R_xlen_t n = rows.nrows();
  if (n == 0)
    return nullptr;

  bool has_acc = rows.containsElementNamed("accumulator");
  Rcpp::RObject acc_col_obj;
  if (has_acc) {
    acc_col_obj = rows["accumulator"];
  }

  Rcpp::CharacterVector dist_col = rows.containsElementNamed("dist")
                                       ? Rcpp::CharacterVector(rows["dist"])
                                       : Rcpp::CharacterVector();
  Rcpp::NumericVector onset_col = rows.containsElementNamed("onset")
                                      ? Rcpp::NumericVector(rows["onset"])
                                      : Rcpp::NumericVector();
  Rcpp::NumericVector q_col = rows.containsElementNamed("q")
                                  ? Rcpp::NumericVector(rows["q"])
                                  : Rcpp::NumericVector();
  Rcpp::NumericVector shared_q_col =
      rows.containsElementNamed("shared_trigger_q")
          ? Rcpp::NumericVector(rows["shared_trigger_q"])
          : Rcpp::NumericVector();
  Rcpp::CharacterVector shared_col =
      rows.containsElementNamed("shared_trigger_id")
          ? Rcpp::CharacterVector(rows["shared_trigger_id"])
          : Rcpp::CharacterVector();
  Rcpp::List comps_col = rows.containsElementNamed("components")
                             ? Rcpp::List(rows["components"])
                             : Rcpp::List();
  Rcpp::List params_list_col = rows.containsElementNamed("params")
                                   ? Rcpp::List(rows["params"])
                                   : Rcpp::List();

  std::unordered_set<std::string> base_cols = {
      "trial",   "component", "accumulator",      "type",
      "outcome", "rt",        "params",           "onset",
      "q",       "condition", "component_weight", "shared_trigger_id",
      "dist",    "components"};
  std::vector<std::string> param_cols;
  std::vector<SEXP> param_col_data;
  std::vector<int> param_col_types;
  Rcpp::CharacterVector df_names = rows.names();
  if (!df_names.isNULL()) {
    for (R_xlen_t i = 0; i < df_names.size(); ++i) {
      if (df_names[i] == NA_STRING)
        continue;
      std::string col_name = Rcpp::as<std::string>(df_names[i]);
      if (base_cols.find(col_name) != base_cols.end())
        continue;
      SEXP column = rows[col_name];
      param_cols.push_back(col_name);
      param_col_data.push_back(column);
      param_col_types.push_back(TYPEOF(column));
    }
  }

  auto params_set = std::make_unique<TrialParamSet>();
  params_set->acc_params.assign(ctx.accumulators.size(),
                                TrialAccumulatorParams{});

  auto upsert_param_entry = [](std::vector<uuber::ProtoParamEntry> &entries,
                               uuber::ProtoParamEntry entry) {
    for (auto &existing : entries) {
      if (existing.name == entry.name) {
        existing = std::move(entry);
        return;
      }
    }
    entries.push_back(std::move(entry));
  };

  auto append_scalar_entry = [&](std::vector<uuber::ProtoParamEntry> &entries,
                                 const std::string &name, double value,
                                 bool logical_scalar = false) {
    if (!std::isfinite(value))
      return;
    uuber::ProtoParamEntry entry;
    entry.name = name;
    if (logical_scalar) {
      entry.tag = uuber::ParamValueTag::LogicalScalar;
      entry.logical_scalar = static_cast<int>(value != 0.0);
    } else {
      entry.tag = uuber::ParamValueTag::NumericScalar;
      entry.numeric_scalar = value;
    }
    upsert_param_entry(entries, std::move(entry));
  };

  auto append_numeric_vector_entry =
      [&](std::vector<uuber::ProtoParamEntry> &entries, const std::string &name,
          const Rcpp::NumericVector &vec) {
        if (vec.size() == 0)
          return;
        uuber::ProtoParamEntry entry;
        entry.name = name;
        entry.tag = uuber::ParamValueTag::NumericVector;
        entry.numeric_values.reserve(vec.size());
        for (double val : vec) {
          entry.numeric_values.push_back(val);
        }
        upsert_param_entry(entries, std::move(entry));
      };

  auto append_logical_vector_entry =
      [&](std::vector<uuber::ProtoParamEntry> &entries, const std::string &name,
          const Rcpp::LogicalVector &vec) {
        if (vec.size() == 0)
          return;
        uuber::ProtoParamEntry entry;
        entry.name = name;
        entry.tag = uuber::ParamValueTag::LogicalVector;
        entry.logical_values.reserve(vec.size());
        for (int val : vec) {
          entry.logical_values.push_back(val);
        }
        upsert_param_entry(entries, std::move(entry));
      };

  for (R_xlen_t i = 0; i < n; ++i) {
    int acc_idx = -1;
    if (has_acc) {
      if (!Rf_isNumeric(acc_col_obj)) {
        Rcpp::stop("Column 'accumulator' must be numeric (1-based index)");
      }
      Rcpp::NumericVector acc_num(acc_col_obj);
      if (i >= acc_num.size())
        continue;
      double raw_val = acc_num[i];
      if (Rcpp::NumericVector::is_na(raw_val))
        continue;
      int idx_val = static_cast<int>(std::llround(raw_val)) - 1;
      if (idx_val < 0 || idx_val >= static_cast<int>(ctx.accumulators.size()))
        continue;
      acc_idx = idx_val;
    }
    if (acc_idx < 0 || acc_idx >= static_cast<int>(ctx.accumulators.size()))
      continue;
    const uuber::NativeAccumulator &base = ctx.accumulators[acc_idx];

    TrialAccumulatorParams override;
    override.onset_kind = base.onset_kind;
    override.onset_source_acc_idx = base.onset_source_acc_idx;
    override.onset_source_pool_idx = base.onset_source_pool_idx;
    override.onset_lag = base.onset_lag;
    override.onset =
        (i < onset_col.size() && !Rcpp::NumericVector::is_na(onset_col[i]))
            ? static_cast<double>(onset_col[i])
            : base.onset;
    // Resolve shared trigger first
    bool has_shared = false;
    if (i < shared_col.size() && shared_col[i] != NA_STRING) {
      override.shared_trigger_id = Rcpp::as<std::string>(shared_col[i]);
      has_shared = true;
    } else {
      override.shared_trigger_id = base.shared_trigger_id;
      has_shared = !override.shared_trigger_id.empty();
    }
    if (override.shared_trigger_id != base.shared_trigger_id) {
      params_set->shared_trigger_layout_matches_context = false;
    }
    // Shared gate probability from column or base; per-acc q=0 on success path.
    if (has_shared) {
      double gate_q = base.q;
      if (i < shared_q_col.size() &&
          !Rcpp::NumericVector::is_na(shared_q_col[i])) {
        gate_q = static_cast<double>(shared_q_col[i]);
      }
      gate_q = clamp_probability(gate_q);
      override.shared_q = gate_q;
      override.q = 0.0;
    } else {
      if (i < q_col.size() && !Rcpp::NumericVector::is_na(q_col[i])) {
        override.q = clamp_probability(static_cast<double>(q_col[i]));
      } else {
        override.q = base.q;
      }
      override.shared_q = override.q;
    }
    SEXP comp_entry = (i < comps_col.size()) ? comps_col[i] : R_NilValue;
    if (comp_entry != R_NilValue) {
      override.components = string_vector_from_entry(comp_entry);
      override.has_components = !override.components.empty();
      if (override.has_components) {
        override.component_indices =
            component_indices_from_labels(ctx, override.components);
      }
    } else {
      override.has_components = false;
    }

    std::string dist_name = base.dist;
    if (i < dist_col.size() && dist_col[i] != NA_STRING) {
      dist_name = Rcpp::as<std::string>(dist_col[i]);
    }
    override.dist_cfg = base.dist_cfg;
    std::vector<uuber::ProtoParamEntry> param_entries;
    if (params_list_col.size() > 0) {
      SEXP param_entry =
          (i < params_list_col.size()) ? params_list_col[i] : R_NilValue;
      if (param_entry != R_NilValue) {
        Rcpp::List param_list(param_entry);
        std::vector<uuber::ProtoParamEntry> entries =
            params_from_rcpp(param_list, dist_name);
        for (auto &entry : entries) {
          upsert_param_entry(param_entries, entry);
        }
      }
    }
    for (std::size_t pc = 0; pc < param_cols.size(); ++pc) {
      SEXP column = param_col_data[pc];
      int type = param_col_types[pc];
      const std::string &col_name = param_cols[pc];
      switch (type) {
      case REALSXP: {
        Rcpp::NumericVector vec(column);
        if (i < vec.size()) {
          double val = vec[i];
          if (!Rcpp::NumericVector::is_na(val)) {
            append_scalar_entry(param_entries, col_name, val, false);
          }
        }
        break;
      }
      case INTSXP: {
        Rcpp::IntegerVector vec(column);
        if (i < vec.size()) {
          int val = vec[i];
          if (val != NA_INTEGER) {
            append_scalar_entry(param_entries, col_name,
                                static_cast<double>(val), false);
          }
        }
        break;
      }
      case LGLSXP: {
        Rcpp::LogicalVector vec(column);
        if (i < vec.size()) {
          int val = vec[i];
          if (val != NA_LOGICAL) {
            append_scalar_entry(param_entries, col_name,
                                static_cast<double>(val), true);
          }
        }
        break;
      }
      case STRSXP:
        // Strings are not supported as distribution params; skip.
        break;
      case VECSXP: {
        Rcpp::List vec(column);
        if (i >= vec.size())
          break;
        SEXP cell = vec[i];
        if (cell == R_NilValue)
          break;
        if (Rf_isReal(cell)) {
          Rcpp::NumericVector nv(cell);
          if (nv.size() == 1) {
            double val = nv[0];
            if (!Rcpp::NumericVector::is_na(val)) {
              append_scalar_entry(param_entries, col_name, val, false);
            }
          } else {
            append_numeric_vector_entry(param_entries, col_name, nv);
          }
        } else if (Rf_isInteger(cell)) {
          Rcpp::IntegerVector iv(cell);
          if (iv.size() == 1) {
            int val = iv[0];
            if (val != NA_INTEGER) {
              append_scalar_entry(param_entries, col_name,
                                  static_cast<double>(val), false);
            }
          } else {
            Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(iv);
            append_numeric_vector_entry(param_entries, col_name, nv);
          }
        } else if (Rf_isLogical(cell)) {
          Rcpp::LogicalVector lv(cell);
          if (lv.size() == 1) {
            int val = lv[0];
            if (val != NA_LOGICAL) {
              append_scalar_entry(param_entries, col_name,
                                  static_cast<double>(val), true);
            }
          } else {
            append_logical_vector_entry(param_entries, col_name, lv);
          }
        }
        break;
      }
      default:
        break;
      }
    }
    if (!param_entries.empty() || dist_name != base.dist) {
      if (param_entries.empty()) {
        continue;
      }
      override.dist_cfg = resolve_acc_params_entries(dist_name, param_entries);
    }

    override.has_override = true;
    params_set->acc_params[acc_idx] = std::move(override);
  }

  return params_set;
}

} // namespace

Rcpp::List native_component_plan_impl(const Rcpp::List &structure,
                                      const Rcpp::DataFrame *trial_rows,
                                      double trial_id,
                                      const std::string *forced_component) {
  Rcpp::DataFrame comp_tbl(structure["components"]);
  Rcpp::CharacterVector component_ids = comp_tbl["component_id"];
  Rcpp::NumericVector base_weights = comp_tbl["weight"];
  std::string mode = "fixed";
  if (comp_tbl.containsElementNamed("mode")) {
    Rcpp::CharacterVector mode_vec(comp_tbl["mode"]);
    if (mode_vec.size() > 0 && mode_vec[0] != NA_STRING) {
      mode = Rcpp::as<std::string>(mode_vec[0]);
    }
  }
  std::string reference;
  if (comp_tbl.containsElementNamed("reference")) {
    Rcpp::CharacterVector ref_vec(comp_tbl["reference"]);
    if (ref_vec.size() > 0 && ref_vec[0] != NA_STRING) {
      reference = Rcpp::as<std::string>(ref_vec[0]);
    }
  }
  Rcpp::List attrs_list;
  if (comp_tbl.containsElementNamed("attrs")) {
    attrs_list = Rcpp::List(comp_tbl["attrs"]);
  }

  std::vector<std::string> all_components;
  all_components.reserve(component_ids.size());
  std::unordered_map<std::string, double> base_weight_map;
  base_weight_map.reserve(component_ids.size());
  for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
    std::string id = Rcpp::as<std::string>(component_ids[i]);
    double weight =
        (i < base_weights.size()) ? static_cast<double>(base_weights[i]) : 0.0;
    all_components.push_back(id);
    base_weight_map[id] = weight;
  }

  std::vector<std::string> selected_components;
  if (forced_component && !forced_component->empty()) {
    selected_components.push_back(*forced_component);
  } else if (trial_rows) {
    Rcpp::CharacterVector trial_comp = (*trial_rows)["component"];
    std::unordered_set<std::string> trial_filter;
    trial_filter.reserve(static_cast<std::size_t>(trial_comp.size()));
    for (R_xlen_t i = 0; i < trial_comp.size(); ++i) {
      trial_filter.insert(Rcpp::as<std::string>(trial_comp[i]));
    }
    if (!trial_filter.empty()) {
      for (const auto &comp : all_components) {
        if (trial_filter.count(comp)) {
          selected_components.push_back(comp);
        }
      }
    }
  }
  if (selected_components.empty()) {
    selected_components = all_components;
  }

  std::vector<double> weights(selected_components.size(), 0.0);

  if (mode == "sample") {
    // Weight parameters come from trial_rows for non-reference components.
    // Reference weight = 1 - sum(nonref), then normalize to guard against
    // rounding.
    std::unordered_map<std::string, int> comp_index;
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      comp_index[selected_components[i]] = static_cast<int>(i);
    }
    std::string ref_id =
        reference.empty() ? selected_components.front() : reference;
    double sum_nonref = 0.0;
    for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
      std::string cid = Rcpp::as<std::string>(component_ids[i]);
      if (cid == ref_id)
        continue;
      auto it_idx = comp_index.find(cid);
      if (it_idx == comp_index.end())
        continue;
      std::string wp;
      if (attrs_list.size() > 0 && i < attrs_list.size()) {
        Rcpp::List attr(attrs_list[i]);
        if (attr.containsElementNamed("weight_param")) {
          wp = Rcpp::as<std::string>(attr["weight_param"]);
        }
      }
      if (wp.empty()) {
        Rcpp::stop("Component '%s' missing weight_param for sampled mixture",
                   cid.c_str());
      }
      double val = std::numeric_limits<double>::quiet_NaN();
      if (trial_rows && trial_rows->containsElementNamed(wp.c_str())) {
        Rcpp::NumericVector col((*trial_rows)[wp]);
        if (col.size() > 0) {
          double cand = col[0];
          if (std::isfinite(cand))
            val = cand;
        }
      }
      if (!std::isfinite(val)) {
        auto it = base_weight_map.find(cid);
        if (it != base_weight_map.end())
          val = it->second;
      }
      if (!std::isfinite(val) || val < 0.0 || val > 1.0) {
        Rcpp::stop("Mixture weight '%s' must be a probability in [0,1]",
                   wp.c_str());
      }
      weights[static_cast<std::size_t>(it_idx->second)] = val;
      sum_nonref += val;
    }
    auto ref_it = comp_index.find(ref_id);
    if (ref_it != comp_index.end()) {
      double ref_weight = 1.0 - sum_nonref;
      if (ref_weight < -1e-8) {
        Rcpp::stop(
            "Mixture weights sum to >1; check non-reference weight params");
      }
      if (ref_weight < 0.0)
        ref_weight = 0.0;
      weights[static_cast<std::size_t>(ref_it->second)] = ref_weight;
    }
  } else {
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      auto it = base_weight_map.find(selected_components[i]);
      double w = (it != base_weight_map.end()) ? it->second : 0.0;
      if (!std::isfinite(w) || w < 0.0)
        w = 0.0;
      weights[i] = w;
    }
  }

  double total = 0.0;
  for (double w : weights)
    total += w;
  if (total <= 0.0) {
    double uniform = 1.0 / static_cast<double>(weights.size());
    std::fill(weights.begin(), weights.end(), uniform);
  } else {
    for (double &w : weights) {
      w /= total;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("components") = Rcpp::CharacterVector(
          selected_components.begin(), selected_components.end()),
      Rcpp::Named("weights") =
          Rcpp::NumericVector(weights.begin(), weights.end()));
}

double normalize_weights(std::vector<double> &weights) {
  if (weights.empty())
    return 0.0;
  bool need_uniform = false;
  double total = 0.0;
  for (double w : weights) {
    if (!std::isfinite(w) || w < 0.0) {
      need_uniform = true;
      break;
    }
    total += w;
  }
  if (!need_uniform && total > 0.0) {
    for (double &w : weights) {
      w /= total;
    }
    return total;
  }
  double uniform = 1.0 / static_cast<double>(weights.size());
  std::fill(weights.begin(), weights.end(), uniform);
  return 1.0;
}

bool build_component_mix(const Rcpp::CharacterVector &component_ids,
                         const Rcpp::NumericVector &weights_in,
                         Rcpp::Nullable<Rcpp::String> forced_component,
                         std::vector<std::string> &components_out,
                         std::vector<double> &weights_out) {
  if (component_ids.size() == 0)
    return false;
  components_out.clear();
  components_out.reserve(component_ids.size());
  for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
    Rcpp::String comp(component_ids[i]);
    if (comp == NA_STRING) {
      components_out.emplace_back("__default__");
    } else {
      components_out.emplace_back(comp.get_cstring());
    }
  }
  weights_out.assign(components_out.size(),
                     std::numeric_limits<double>::quiet_NaN());
  if (weights_in.size() == component_ids.size()) {
    for (R_xlen_t i = 0; i < weights_in.size(); ++i) {
      weights_out[static_cast<std::size_t>(i)] = weights_in[i];
    }
  }
  normalize_weights(weights_out);
  if (forced_component.isNotNull()) {
    std::string forced = Rcpp::as<std::string>(forced_component);
    auto it = std::find(components_out.begin(), components_out.end(), forced);
    if (it == components_out.end()) {
      return false;
    }
    std::string selected = *it;
    components_out.clear();
    components_out.push_back(selected);
    weights_out.clear();
    weights_out.push_back(1.0);
  }
  return !components_out.empty();
}

double native_trial_mixture_internal(
    uuber::NativeContext &ctx, int node_id, double t,
    const std::vector<std::string> &component_ids,
    const std::vector<double> &weights,
    Rcpp::Nullable<Rcpp::String> forced_component,
    const std::vector<int> &competitor_ids, const TrialParamSet *params_ptr,
    const std::vector<ComponentCacheEntry> &cache_entries,
    const std::string &keep_label = std::string(),
    const Rcpp::List &guess_donors = Rcpp::List()) {
  if (!std::isfinite(t) || t < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  static const std::unordered_set<int> kEmptyForced;
  double total = 0.0;
  const bool has_keep = !keep_label.empty();
  const std::string keep_norm =
      has_keep ? normalize_label_string(keep_label) : std::string();
  int keep_outcome_idx = -1;
  int keep_outcome_idx_context = -1;
  if (has_keep) {
    if (keep_label == "NA") {
      keep_outcome_idx = kOutcomeIdxNA;
    } else {
      keep_outcome_idx = outcome_index_of(ctx, keep_label);
      keep_outcome_idx_context = keep_outcome_idx;
    }
  }
  SharedTriggerPlan trigger_plan = build_shared_trigger_plan(ctx, params_ptr);
  TrialParamSet trigger_scratch;
  std::vector<int> component_indices;
  component_indices.reserve(component_ids.size());
  for (const auto &comp_id : component_ids) {
    component_indices.push_back(component_index_of(ctx, comp_id));
  }
  for (std::size_t i = 0; i < component_ids.size(); ++i) {
    int comp_idx = component_indices[i];
    const ComponentCacheEntry *cache_entry_ptr = nullptr;
    ComponentCacheEntry fallback;
    if (i < cache_entries.size()) {
      cache_entry_ptr = &cache_entries[i];
    } else {
      fallback = default_component_cache_entry(component_ids[i]);
      cache_entry_ptr = &fallback;
    }
    double contrib = 0.0;
    bool used_guess_shortcut = false;
    if (has_keep) {
      if (comp_idx >= 0 &&
          comp_idx < static_cast<int>(ctx.component_info.size())) {
        const uuber::ComponentGuessPolicy &guess =
            ctx.component_info[static_cast<std::size_t>(comp_idx)].guess;
        bool guess_applies = false;
        if (guess.valid) {
          if (guess.target == "__GUESS__") {
            guess_applies = (keep_label == "__GUESS__");
          } else if (guess.target_outcome_idx >= 0) {
            guess_applies = (keep_outcome_idx == guess.target_outcome_idx);
          } else {
            std::string target_norm = normalize_label_string(guess.target);
            guess_applies = (!keep_norm.empty() && keep_norm == target_norm);
          }
        }
        if (guess_applies) {
          contrib = accumulate_component_guess_density(
              ctx, keep_label, keep_outcome_idx, t, component_ids[i], comp_idx,
              params_ptr, cache_entry_ptr->trial_type_key);
          used_guess_shortcut = true;
        }
      }
    }
    if (!used_guess_shortcut) {
      std::vector<int> comp_competitors;
      const std::vector<int> &comp_use = filter_competitor_ids(
          ctx, competitor_ids, component_ids[i], comp_competitors);
      contrib = node_density_with_shared_triggers_plan(
          ctx, node_id, t, component_ids[i], comp_idx, kEmptyForced,
          kEmptyForced, comp_use, params_ptr,
          cache_entry_ptr->trial_type_key, false, keep_label,
          keep_outcome_idx_context, &trigger_plan, &trigger_scratch);
      if (!std::isfinite(contrib) || contrib < 0.0)
        contrib = 0.0;
      if (has_keep) {
        double keep_w = component_keep_weight(ctx, comp_idx, keep_outcome_idx);
        contrib *= keep_w;
      }
      if (guess_donors.size() > 0) {
        double donor_contrib = accumulate_plan_guess_density(
            ctx, guess_donors, t, component_ids[i],
            cache_entry_ptr->trial_type_key, params_ptr, false);
        contrib += donor_contrib;
      }
    }
    total += weights[i] * contrib;
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
}

// [[Rcpp::export(name = "native_trial_mixture_cpp")]]
double native_trial_mixture_driver(
    SEXP ctxSEXP, int node_id, double t, Rcpp::CharacterVector component_ids,
    Rcpp::NumericVector weights, Rcpp::Nullable<Rcpp::String> forced_component,
    Rcpp::IntegerVector competitor_ids,
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows,
    Rcpp::Nullable<Rcpp::List> guess_donors) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::vector<std::string> components;
  std::vector<double> mix_weights;
  if (!build_component_mix(component_ids, weights, forced_component, components,
                           mix_weights)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::vector<int> comp_ids = integer_vector_to_std(competitor_ids, false);
  std::unique_ptr<TrialParamSet> param_holder;
  const TrialParamSet *params_ptr = nullptr;
  if (!trial_rows.isNull()) {
    param_holder = build_trial_params_from_df(*ctx, trial_rows);
    params_ptr = param_holder.get();
  }
  std::vector<ComponentCacheEntry> cache_entries =
      build_component_cache_entries(components);
  return native_trial_mixture_internal(
      *ctx, node_id, t, components, mix_weights, forced_component, comp_ids,
      params_ptr, cache_entries, std::string(),
      guess_donors.isNotNull() ? Rcpp::List(guess_donors.get()) : Rcpp::List());
}

double mix_outcome_mass(SEXP ctxSEXP, const uuber::OutcomeContextInfo &info,
                        const std::vector<std::string> &comp_labels,
                        const std::vector<double> &comp_weights,
                        const TrialParamSet *params_ptr, double rel_tol,
                        double abs_tol, int max_depth,
                        const std::string &outcome_label_context) {
  if (info.node_id < 0)
    return 0.0;
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  // Cache NA-map integrals when no RT is provided (maps_to_na cases).
  uuber::NAMapCacheKey cache_key;
  bool use_cache = params_ptr != nullptr;
  if (use_cache) {
    cache_key = build_na_map_cache_key(info, params_ptr, comp_labels,
                                       comp_weights, outcome_label_context);
    auto hit = ctx->na_map_cache.find(cache_key);
    if (hit != ctx->na_map_cache.end()) {
      return hit->second;
    }
  }
  Rcpp::IntegerVector comp_ids(info.competitor_ids.begin(),
                               info.competitor_ids.end());
  int outcome_idx_context = -1;
  if (!outcome_label_context.empty()) {
    outcome_idx_context = (outcome_label_context == "NA")
                              ? kOutcomeIdxNA
                              : outcome_index_of(*ctx, outcome_label_context);
  }
  std::vector<int> comp_indices;
  comp_indices.reserve(comp_labels.size());
  for (const auto &comp : comp_labels) {
    comp_indices.push_back(component_index_of(*ctx, comp));
  }
  double total = 0.0;
  for (std::size_t i = 0; i < comp_labels.size(); ++i) {
    double w = (i < comp_weights.size()) ? comp_weights[i] : 0.0;
    if (w <= 0.0)
      continue;
    const std::string &comp = comp_labels[i];
    if (!outcome_allows_component(info, comp))
      continue;
    int comp_idx = comp_indices[i];
    Rcpp::Nullable<Rcpp::String> comp_val(comp.empty() ? R_NilValue
                                                       : Rcpp::wrap(comp));
    double p = native_outcome_probability_impl(
        ctxSEXP, info.node_id, std::numeric_limits<double>::infinity(),
        comp_val, R_NilValue, R_NilValue, comp_ids, rel_tol, abs_tol, max_depth,
        params_ptr, std::string(), false, std::string());
    if (std::isfinite(p) && p > 0.0) {
      double keep_w =
          component_keep_weight(*ctx, comp_idx, outcome_idx_context);
      total += w * p * keep_w;
    }
  }
  double out = clamp_probability(total);
  if (use_cache) {
    if (static_cast<int>(ctx->na_map_cache.size()) >= ctx->na_cache_limit) {
      ctx->na_map_cache.clear();
      ctx->na_cache_order.clear();
    }
    ctx->na_map_cache.emplace(std::move(cache_key), out);
  }
  return out;
}

struct RankedState {
  double weight{0.0};
  std::vector<int> forced_complete;
  std::vector<int> forced_survive;
  ExactSourceTimeMap exact_source_times;
  SourceTimeBoundsMap source_time_bounds;
};

inline std::string ranked_state_key(const std::vector<int> &forced_complete,
                                    const std::vector<int> &forced_survive,
                                    const ExactSourceTimeMap &exact_source_times,
                                    const SourceTimeBoundsMap &source_time_bounds) {
  std::ostringstream oss;
  oss << "c:";
  for (std::size_t i = 0; i < forced_complete.size(); ++i) {
    if (i > 0)
      oss << ",";
    oss << forced_complete[i];
  }
  oss << "|s:";
  for (std::size_t i = 0; i < forced_survive.size(); ++i) {
    if (i > 0)
      oss << ",";
    oss << forced_survive[i];
  }
  if (!exact_source_times.empty()) {
    std::vector<std::pair<int, double>> exact_pairs;
    exact_pairs.reserve(exact_source_times.size());
    for (const auto &kv : exact_source_times) {
      exact_pairs.push_back(kv);
    }
    std::sort(exact_pairs.begin(), exact_pairs.end(),
              [](const auto &lhs, const auto &rhs) {
                return lhs.first < rhs.first;
              });
    oss << "|x:";
    for (std::size_t i = 0; i < exact_pairs.size(); ++i) {
      if (i > 0)
        oss << ",";
      oss << exact_pairs[i].first << "@" << std::setprecision(17)
          << exact_pairs[i].second;
    }
  }
  if (!source_time_bounds.empty()) {
    std::vector<std::pair<int, std::pair<double, double>>> bound_pairs;
    bound_pairs.reserve(source_time_bounds.size());
    for (const auto &kv : source_time_bounds) {
      bound_pairs.push_back(kv);
    }
    std::sort(bound_pairs.begin(), bound_pairs.end(),
              [](const auto &lhs, const auto &rhs) {
                return lhs.first < rhs.first;
              });
    oss << "|b:";
    for (std::size_t i = 0; i < bound_pairs.size(); ++i) {
      if (i > 0)
        oss << ",";
      oss << bound_pairs[i].first << "@" << std::setprecision(17)
          << bound_pairs[i].second.first << ":" << bound_pairs[i].second.second;
    }
  }
  return oss.str();
}

std::vector<RankedState>
collapse_ranked_states(const std::vector<RankedState> &states, double eps) {
  if (states.empty()) {
    return {};
  }
  std::unordered_map<std::string, std::size_t> key_index;
  std::vector<RankedState> collapsed;
  collapsed.reserve(states.size());
  for (const RankedState &state : states) {
    if (!std::isfinite(state.weight) || state.weight <= eps) {
      continue;
    }
    std::string key =
        ranked_state_key(state.forced_complete, state.forced_survive,
                         state.exact_source_times, state.source_time_bounds);
    auto it = key_index.find(key);
    if (it == key_index.end()) {
      key_index.emplace(key, collapsed.size());
      collapsed.push_back(state);
    } else {
      collapsed[it->second].weight += state.weight;
    }
  }
  std::vector<RankedState> out;
  out.reserve(collapsed.size());
  for (RankedState &state : collapsed) {
    if (std::isfinite(state.weight) && state.weight > eps) {
      out.push_back(std::move(state));
    }
  }
  return out;
}

std::vector<int>
collect_competitor_sources(const uuber::NativeContext &ctx,
                           const std::vector<int> &competitor_ids) {
  if (competitor_ids.empty()) {
    return {};
  }
  std::vector<int> out;
  for (int node_id : competitor_ids) {
    if (node_id == NA_INTEGER) {
      continue;
    }
    const uuber::NativeNode &node = fetch_node(ctx, node_id);
    std::vector<int> source_ids = ensure_source_ids(ctx, node);
    out.insert(out.end(), source_ids.begin(), source_ids.end());
  }
  sort_unique(out);
  return out;
}

double ranked_prefix_density_single_params(
    const uuber::NativeContext &ctx, const std::vector<std::string> &labels,
    const std::vector<double> &times, const std::string &component_label,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key) {
  if (labels.empty() || labels.size() != times.size()) {
    return 0.0;
  }

  constexpr double kBranchEps = 1e-18;
  std::vector<int> outcome_indices;
  std::vector<int> node_sequence;
  outcome_indices.reserve(labels.size());
  node_sequence.reserve(labels.size());
  std::unordered_set<int> seen_node_ids;
  seen_node_ids.reserve(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i) {
    const std::string &label = labels[i];
    double t = times[i];
    if (label.empty() || !std::isfinite(t) || t < 0.0) {
      return 0.0;
    }
    int outcome_idx = resolve_outcome_index(ctx, label, component_label);
    if (outcome_idx < 0 ||
        outcome_idx >= static_cast<int>(ctx.outcome_info.size())) {
      return 0.0;
    }
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    if (info.node_id < 0) {
      return 0.0;
    }
    if (!seen_node_ids.insert(info.node_id).second) {
      return 0.0;
    }
    outcome_indices.push_back(outcome_idx);
    node_sequence.push_back(info.node_id);
  }

  std::unordered_set<int> onset_source_ids;
  onset_source_ids.reserve(ctx.accumulators.size());
  const std::vector<TrialAccumulatorParams> *acc_params =
      trial_params ? &trial_params->acc_params : nullptr;
  for (std::size_t acc_i = 0; acc_i < ctx.accumulators.size(); ++acc_i) {
    const uuber::NativeAccumulator &acc = ctx.accumulators[acc_i];
    const TrialAccumulatorParams *override =
        (acc_params && acc_i < acc_params->size())
            ? &((*acc_params)[acc_i])
            : nullptr;
    int onset_kind = override ? override->onset_kind : acc.onset_kind;
    if (onset_kind == uuber::ONSET_AFTER_ACCUMULATOR) {
      int src_idx =
          override ? override->onset_source_acc_idx : acc.onset_source_acc_idx;
      if (src_idx >= 0 && src_idx < static_cast<int>(ctx.accumulators.size())) {
        int src_id = label_id(ctx, ctx.accumulators[src_idx].id);
        if (src_id >= 0 && src_id != NA_INTEGER) {
          onset_source_ids.insert(src_id);
        }
      }
    } else if (onset_kind == uuber::ONSET_AFTER_POOL) {
      int src_idx =
          override ? override->onset_source_pool_idx : acc.onset_source_pool_idx;
      if (src_idx >= 0 && src_idx < static_cast<int>(ctx.pools.size())) {
        int src_id = label_id(ctx, ctx.pools[src_idx].id);
        if (src_id >= 0 && src_id != NA_INTEGER) {
          onset_source_ids.insert(src_id);
        }
      }
    }
  }

  std::vector<ExactSourceTimeMap> observed_prefix_exact(labels.size() + 1);
  for (std::size_t rank_idx = 0; rank_idx < labels.size(); ++rank_idx) {
    observed_prefix_exact[rank_idx + 1] = observed_prefix_exact[rank_idx];
    const uuber::NativeNode &node =
        fetch_node(ctx, node_sequence[rank_idx]);
    std::vector<int> src_ids = ensure_source_ids(ctx, node);
    for (int src_id : src_ids) {
      if (onset_source_ids.count(src_id) == 0) {
        continue;
      }
      observed_prefix_exact[rank_idx + 1][src_id] = times[rank_idx];
    }
  }

  std::vector<RankedState> states;
  RankedState init_state;
  init_state.weight = 1.0;
  states.push_back(std::move(init_state));
  std::unordered_set<int> observed_nodes;
  std::unordered_set<int> future_nodes(node_sequence.begin(),
                                       node_sequence.end());
  observed_nodes.reserve(labels.size());

  for (std::size_t rank_idx = 0; rank_idx < labels.size(); ++rank_idx) {
    double t = times[rank_idx];
    if (!std::isfinite(t) || t < 0.0) {
      return 0.0;
    }
    const ExactSourceTimeMap &prefix_exact = observed_prefix_exact[rank_idx];

    int outcome_idx = outcome_indices[rank_idx];
    const uuber::OutcomeContextInfo &info =
        ctx.outcome_info[static_cast<std::size_t>(outcome_idx)];
    observed_nodes.insert(info.node_id);
    future_nodes.erase(info.node_id);

    std::vector<int> filtered_competitors;
    const std::vector<int> &base_competitors =
        filter_competitor_ids(ctx, info.competitor_ids, component_label,
                              filtered_competitors);
    std::vector<int> competitors;
    competitors.reserve(base_competitors.size());
    for (int node_id : base_competitors) {
      if (node_id == NA_INTEGER || node_id == info.node_id) {
        continue;
      }
      if (observed_nodes.count(node_id) > 0) {
        continue;
      }
      competitors.push_back(node_id);
    }
    sort_unique(competitors);
    std::vector<int> persistent_competitors;
    persistent_competitors.reserve(competitors.size());
    for (int node_id : competitors) {
      if (future_nodes.count(node_id) == 0) {
        persistent_competitors.push_back(node_id);
      }
    }
    std::vector<int> persistent_sources =
        collect_competitor_sources(ctx, persistent_competitors);

    std::vector<RankedState> next_states;
    next_states.reserve(states.size() * 2);

    for (const RankedState &state : states) {
      if (!std::isfinite(state.weight) || state.weight <= kBranchEps) {
        continue;
      }
      ExactSourceTimeMap merged_exact = state.exact_source_times;
      if (!prefix_exact.empty()) {
        for (const auto &kv : prefix_exact) {
          if (std::isfinite(kv.second)) {
            merged_exact[kv.first] = kv.second;
          }
        }
      }
      SourceTimeBoundsMap merged_bounds = state.source_time_bounds;
      if (!merged_exact.empty() && !merged_bounds.empty()) {
        for (const auto &kv : merged_exact) {
          merged_bounds.erase(kv.first);
        }
      }
      std::unordered_set<int> forced_complete =
          make_forced_set(state.forced_complete);
      std::unordered_set<int> forced_survive =
          make_forced_set(state.forced_survive);
      double denom = 1.0;
      if (rank_idx > 0) {
        double lower_t = times[rank_idx - 1];
        denom = evaluate_survival_with_forced(
            info.node_id, forced_complete, forced_survive, component_label,
            component_idx, lower_t, ctx, trial_type_key, trial_params,
            &merged_exact, &merged_bounds);
        if (!std::isfinite(denom) || denom <= kBranchEps) {
          continue;
        }
      }
      NodeEvalState eval_state(
          ctx, t, component_label, forced_complete, forced_survive,
          trial_params, trial_type_key, false, std::string(), component_idx,
          outcome_idx, &merged_exact, &merged_bounds);
      std::vector<ScenarioRecord> scenarios =
          compute_node_scenarios(info.node_id, eval_state);
      if (scenarios.empty()) {
        continue;
      }
      for (const ScenarioRecord &scenario : scenarios) {
        if (!std::isfinite(scenario.weight) ||
            scenario.weight <= kBranchEps) {
          continue;
        }
        double weight = state.weight * (scenario.weight / denom);
        if (!std::isfinite(weight) || weight <= kBranchEps) {
          continue;
        }
        std::vector<int> next_complete = scenario.forced_complete;
        std::vector<int> next_survive = scenario.forced_survive;
        sort_unique(next_complete);
        sort_unique(next_survive);
        if (!competitors.empty()) {
          std::unordered_set<int> next_complete_set =
              make_forced_set(next_complete);
          std::unordered_set<int> next_survive_set =
              make_forced_set(next_survive);
          double surv = competitor_survival_internal(
              ctx, competitors, t, component_label, component_idx,
              next_complete_set, next_survive_set, trial_type_key,
              trial_params, &merged_exact, &merged_bounds);
          if (!std::isfinite(surv) || surv <= kBranchEps) {
            continue;
          }
          weight *= surv;
          if (!std::isfinite(weight) || weight <= kBranchEps) {
            continue;
          }
          if (!persistent_sources.empty()) {
            next_survive = union_vectors(next_survive, persistent_sources);
          }
        }
        RankedState next_state;
        next_state.weight = weight;
        next_state.forced_complete = std::move(next_complete);
        next_state.forced_survive = std::move(next_survive);
        next_state.exact_source_times = state.exact_source_times;
        next_state.source_time_bounds = state.source_time_bounds;
        if (!onset_source_ids.empty()) {
          bool bounds_ok = true;
          for (int src_id : next_state.forced_survive) {
            if (onset_source_ids.count(src_id) == 0) {
              continue;
            }
            if (next_state.exact_source_times.find(src_id) !=
                next_state.exact_source_times.end()) {
              continue;
            }
            if (prefix_exact.find(src_id) != prefix_exact.end()) {
              continue;
            }
            auto bound_it = next_state.source_time_bounds.find(src_id);
            if (bound_it == next_state.source_time_bounds.end()) {
              bound_it = next_state.source_time_bounds
                             .emplace(src_id, std::make_pair(
                                                  0.0,
                                                  std::numeric_limits<double>::infinity()))
                             .first;
            }
            auto &bound = bound_it->second;
            bound.first = std::max(bound.first, t);
            if (!(bound.second > bound.first)) {
              bounds_ok = false;
              break;
            }
          }
          if (!bounds_ok) {
            continue;
          }
          for (int src_id : next_state.forced_complete) {
            if (onset_source_ids.count(src_id) == 0) {
              continue;
            }
            if (next_state.exact_source_times.find(src_id) !=
                next_state.exact_source_times.end()) {
              continue;
            }
            if (prefix_exact.find(src_id) != prefix_exact.end()) {
              continue;
            }
            auto bound_it = next_state.source_time_bounds.find(src_id);
            if (bound_it == next_state.source_time_bounds.end()) {
              bound_it = next_state.source_time_bounds
                             .emplace(src_id, std::make_pair(
                                                  0.0,
                                                  std::numeric_limits<double>::infinity()))
                             .first;
            }
            auto &bound = bound_it->second;
            bound.second = std::min(bound.second, t);
            if (!(bound.second > bound.first)) {
              bounds_ok = false;
              break;
            }
          }
          if (!bounds_ok) {
            continue;
          }
        }
        next_states.push_back(std::move(next_state));
      }
    }

    states = collapse_ranked_states(next_states, kBranchEps);
    if (states.empty()) {
      return 0.0;
    }
  }

  double total = 0.0;
  for (const RankedState &state : states) {
    if (std::isfinite(state.weight) && state.weight > 0.0) {
      total += state.weight;
    }
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
}

double ranked_prefix_density(
    const uuber::NativeContext &ctx, const std::vector<std::string> &labels,
    const std::vector<double> &times, const std::string &component_label,
    int component_idx, const TrialParamSet *trial_params,
    const std::string &trial_type_key, const SharedTriggerPlan *trigger_plan,
    TrialParamSet *scratch_params) {
  if (!trial_params) {
    return ranked_prefix_density_single_params(
        ctx, labels, times, component_label, component_idx, trial_params,
        trial_type_key);
  }

  SharedTriggerPlan local_plan;
  const SharedTriggerPlan *plan_ptr = trigger_plan;
  if (!plan_ptr) {
    local_plan = build_shared_trigger_plan(ctx, trial_params);
    plan_ptr = &local_plan;
  }
  if (plan_ptr->triggers.empty()) {
    return ranked_prefix_density_single_params(
        ctx, labels, times, component_label, component_idx, trial_params,
        trial_type_key);
  }

  const std::vector<SharedTriggerInfo> &triggers = plan_ptr->triggers;
  if (triggers.size() >= 63) {
    Rcpp::stop("Too many shared triggers for ranked likelihood evaluation");
  }

  TrialParamSet local_scratch;
  TrialParamSet &scratch = scratch_params ? *scratch_params : local_scratch;
  scratch = *trial_params;
  for (const auto &trigger : triggers) {
    apply_trigger_state_inplace(scratch, trigger, false);
  }

  const std::uint64_t n_states = 1ULL << triggers.size();
  double total = 0.0;
  std::uint64_t mask = 0;
  for (std::uint64_t idx = 0; idx < n_states; ++idx) {
    double w = plan_ptr->mask_weights.empty()
                   ? shared_trigger_mask_weight(triggers, mask)
                   : plan_ptr->mask_weights[mask];
    if (w > 0.0) {
      double contrib = ranked_prefix_density_single_params(
          ctx, labels, times, component_label, component_idx, &scratch,
          trial_type_key);
      if (std::isfinite(contrib) && contrib > 0.0) {
        total += w * contrib;
      }
    }
    if (idx + 1 == n_states) {
      break;
    }
    std::uint64_t next = idx + 1;
    std::uint64_t next_gray = next ^ (next >> 1);
    std::uint64_t diff = next_gray ^ mask;
    if (diff != 0) {
      int bit = static_cast<int>(__builtin_ctzll(diff));
      bool fail = ((next_gray >> bit) & 1ULL) != 0ULL;
      apply_trigger_state_inplace(
          scratch, triggers[static_cast<std::size_t>(bit)], fail);
    }
    mask = next_gray;
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
}

// [[Rcpp::export]]
double cpp_loglik(SEXP ctxSEXP, Rcpp::NumericMatrix params_mat,
                  Rcpp::DataFrame data_df, Rcpp::LogicalVector ok,
                  Rcpp::IntegerVector expand, double min_ll, double rel_tol,
                  double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  ctx->na_map_cache.clear();
  ctx->na_cache_order.clear();
  const int q_col = 0;
  const int w_col = 1;
  const int t0_col = 2;
  const int p1_col = 3;
  const int p2_col = params_mat.ncol() > 4 ? 4 : -1;
  const int p3_col = params_mat.ncol() > 5 ? 5 : -1;

  Rcpp::NumericVector trial_col = data_df["trial"];
  Rcpp::CharacterVector outcome_col = data_df["R"];
  Rcpp::NumericVector rt_col = data_df["rt"];
  std::vector<Rcpp::CharacterVector> rank_outcome_cols;
  std::vector<Rcpp::NumericVector> rank_rt_cols;
  rank_outcome_cols.push_back(outcome_col);
  rank_rt_cols.push_back(rt_col);
  for (int rank = 2;; ++rank) {
    std::string r_name = "R" + std::to_string(rank);
    std::string rt_name = "rt" + std::to_string(rank);
    bool has_r = data_df.containsElementNamed(r_name.c_str());
    bool has_rt = data_df.containsElementNamed(rt_name.c_str());
    if (!has_r && !has_rt) {
      break;
    }
    if (has_r != has_rt) {
      Rcpp::stop("Ranked observations require paired columns '%s' and '%s'",
                 r_name.c_str(), rt_name.c_str());
    }
    rank_outcome_cols.push_back(Rcpp::CharacterVector(data_df[r_name]));
    rank_rt_cols.push_back(Rcpp::NumericVector(data_df[rt_name]));
  }
  const int rank_width = static_cast<int>(rank_outcome_cols.size());
  const bool ranked_mode = rank_width > 1;
  bool has_component = data_df.containsElementNamed("component");
  Rcpp::CharacterVector comp_col =
      has_component ? Rcpp::CharacterVector(data_df["component"])
                    : Rcpp::CharacterVector();
  bool has_onset = data_df.containsElementNamed("onset");
  Rcpp::NumericVector onset_col =
      has_onset ? Rcpp::NumericVector(data_df["onset"]) : Rcpp::NumericVector();

  const std::vector<TrialAccumulatorParams> base_acc_params =
      build_base_paramset(*ctx).acc_params;
  const ComponentMap &comp_map = ctx->components;
  const std::vector<std::string> &comp_ids = comp_map.ids;
  std::unordered_map<std::string, int> comp_index;
  for (std::size_t i = 0; i < comp_ids.size(); ++i)
    comp_index[comp_ids[i]] = static_cast<int>(i);
  std::vector<int> all_acc_indices(ctx->accumulators.size());
  for (std::size_t i = 0; i < all_acc_indices.size(); ++i)
    all_acc_indices[i] = static_cast<int>(i);
  std::vector<std::vector<int>> comp_acc_indices(comp_ids.size());
  for (std::size_t c = 0; c < comp_ids.size(); ++c) {
    const std::string &cid = comp_ids[c];
    for (std::size_t a = 0; a < ctx->accumulators.size(); ++a) {
      const auto &comps = ctx->accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        comp_acc_indices[c].push_back(static_cast<int>(a));
      }
    }
    if (comp_acc_indices[c].empty()) {
      comp_acc_indices[c] = all_acc_indices;
    }
  }

  std::vector<TrialParamSet> param_sets;
  std::vector<Rcpp::String> outcome_by_trial;
  std::vector<double> rt_by_trial;
  std::vector<std::string> component_by_trial;
  std::vector<int> comp_idx_by_trial;
  std::vector<std::vector<double>> weight_override_by_trial;
  std::vector<std::vector<Rcpp::String>> ranked_outcome_by_trial;
  std::vector<std::vector<double>> ranked_rt_by_trial;

  const std::size_t comp_count = comp_map.ids.size();

  const R_xlen_t n_rows = params_mat.nrow();
  if (trial_col.size() != n_rows) {
    Rcpp::stop(
        "Parameter matrix rows (%lld) must match nested likelihood rows (%lld)",
        static_cast<long long>(n_rows),
        static_cast<long long>(trial_col.size()));
  }
  int current_trial_label = std::numeric_limits<int>::min();
  int trial_idx = -1;
  std::size_t acc_cursor = 0;
  const std::vector<int> *current_acc_order = &all_acc_indices;
  for (R_xlen_t r = 0; r < n_rows; ++r) {
    int trial_label = static_cast<int>(trial_col[r]);
    if (r == 0 || trial_label != current_trial_label) {
      current_trial_label = trial_label;
      ++trial_idx;
      TrialParamSet ps;
      ps.acc_params = base_acc_params;
      param_sets.push_back(std::move(ps));
      outcome_by_trial.push_back(Rcpp::String(NA_STRING));
      rt_by_trial.push_back(std::numeric_limits<double>::quiet_NaN());
      component_by_trial.emplace_back();
      comp_idx_by_trial.push_back(-1);
      weight_override_by_trial.push_back(std::vector<double>(
          comp_count, std::numeric_limits<double>::quiet_NaN()));
      if (ranked_mode) {
        ranked_outcome_by_trial.push_back(
            std::vector<Rcpp::String>(static_cast<std::size_t>(rank_width),
                                      Rcpp::String(NA_STRING)));
        ranked_rt_by_trial.push_back(std::vector<double>(
            static_cast<std::size_t>(rank_width),
            std::numeric_limits<double>::quiet_NaN()));
      }
      acc_cursor = 0;
      current_acc_order = &all_acc_indices;
    }
    if (has_component && comp_idx_by_trial[trial_idx] < 0 &&
        r < comp_col.size()) {
      Rcpp::String comp_val(comp_col[r]);
      if (comp_val != NA_STRING) {
        std::string comp_str(comp_val.get_cstring());
        auto it = comp_index.find(comp_str);
        if (it != comp_index.end()) {
          comp_idx_by_trial[trial_idx] = it->second;
          current_acc_order =
              &comp_acc_indices[static_cast<std::size_t>(it->second)];
          component_by_trial[trial_idx] = comp_str;
        }
      }
    }
    int acc_idx = (acc_cursor < current_acc_order->size())
                      ? (*current_acc_order)[acc_cursor]
                      : all_acc_indices[acc_cursor % all_acc_indices.size()];
    ++acc_cursor;
    TrialAccumulatorParams &tap =
        param_sets[static_cast<std::size_t>(trial_idx)]
            .acc_params[static_cast<std::size_t>(acc_idx)];
    double q_val = clamp_probability(params_mat(r, q_col));
    if (!tap.shared_trigger_id.empty()) {
      tap.shared_q = q_val;
      tap.q = 0.0;
    } else {
      tap.q = q_val;
      tap.shared_q = q_val;
    }
    tap.dist_cfg.t0 = params_mat(r, t0_col);
    if (has_onset && r < onset_col.size()) {
      double onset_val = static_cast<double>(onset_col[r]);
      if (std::isfinite(onset_val)) {
        tap.onset = onset_val;
      }
    }
    tap.has_override = true;
    if (tap.dist_cfg.code == uuber::ACC_DIST_LOGNORMAL ||
        tap.dist_cfg.code == uuber::ACC_DIST_GAMMA ||
        tap.dist_cfg.code == uuber::ACC_DIST_EXGAUSS) {
      tap.dist_cfg.p1 = params_mat(r, p1_col);
      if (p2_col >= 0)
        tap.dist_cfg.p2 = params_mat(r, p2_col);
      if (p3_col >= 0)
        tap.dist_cfg.p3 = params_mat(r, p3_col);
    }

    for (std::size_t c = 0; c < comp_map.leader_idx.size(); ++c) {
      if (comp_map.leader_idx[c] == acc_idx) {
        weight_override_by_trial[static_cast<std::size_t>(trial_idx)][c] =
            params_mat(r, w_col);
      }
    }

    if (outcome_by_trial[trial_idx] == NA_STRING && r < outcome_col.size()) {
      outcome_by_trial[trial_idx] = Rcpp::String(outcome_col[r]);
    }
    if (!std::isfinite(rt_by_trial[trial_idx]) && r < rt_col.size()) {
      rt_by_trial[trial_idx] = static_cast<double>(rt_col[r]);
    }
    if (has_component && component_by_trial[trial_idx].empty() &&
        r < comp_col.size()) {
      Rcpp::String comp_val(comp_col[r]);
      if (comp_val != NA_STRING) {
        component_by_trial[trial_idx] = std::string(comp_val.get_cstring());
      }
    }
    if (ranked_mode) {
      std::vector<Rcpp::String> &trial_ranked_out =
          ranked_outcome_by_trial[trial_idx];
      std::vector<double> &trial_ranked_rt = ranked_rt_by_trial[trial_idx];
      for (int rank_idx = 0; rank_idx < rank_width; ++rank_idx) {
        if (trial_ranked_out[static_cast<std::size_t>(rank_idx)] == NA_STRING &&
            r < rank_outcome_cols[static_cast<std::size_t>(rank_idx)].size()) {
          Rcpp::String cand(
              rank_outcome_cols[static_cast<std::size_t>(rank_idx)][r]);
          if (cand != NA_STRING) {
            trial_ranked_out[static_cast<std::size_t>(rank_idx)] = cand;
          }
        }
        if (!std::isfinite(
                trial_ranked_rt[static_cast<std::size_t>(rank_idx)]) &&
            r < rank_rt_cols[static_cast<std::size_t>(rank_idx)].size()) {
          double cand_rt =
              rank_rt_cols[static_cast<std::size_t>(rank_idx)][r];
          if (std::isfinite(cand_rt)) {
            trial_ranked_rt[static_cast<std::size_t>(rank_idx)] = cand_rt;
          }
        }
      }
    }
  }

  int n_trials = static_cast<int>(param_sets.size());
  std::vector<SharedTriggerPlan> trigger_plans;
  trigger_plans.reserve(static_cast<std::size_t>(n_trials));
  for (int t = 0; t < n_trials; ++t) {
    trigger_plans.push_back(build_shared_trigger_plan(
        *ctx, &param_sets[static_cast<std::size_t>(t)]));
  }

  const std::vector<std::string> &default_components = comp_map.ids;
  const std::size_t n_components = default_components.size();
  const std::vector<ComponentCacheEntry> default_cache_entries =
      build_component_cache_entries(default_components);

  std::vector<std::vector<double>> weights_by_trial(
      static_cast<std::size_t>(n_trials),
      std::vector<double>(n_components, 0.0));
  for (int t = 0; t < n_trials; ++t) {
    std::vector<double> w = comp_map.base_weights;
    const auto &overrides =
        weight_override_by_trial[static_cast<std::size_t>(t)];
    for (std::size_t c = 0; c < w.size(); ++c) {
      double cand = overrides[c];
      if (std::isfinite(cand))
        w[c] = cand;
    }
    double total = 0.0;
    for (double v : w)
      total += v;
    if (!std::isfinite(total) || total <= 0.0) {
      double uni = 1.0 / static_cast<double>(w.size());
      std::fill(w.begin(), w.end(), uni);
    } else {
      for (double &v : w)
        v /= total;
    }
    weights_by_trial[static_cast<std::size_t>(t)] = std::move(w);
  }

  if (ok.size() == 0) {
    ok = Rcpp::LogicalVector(n_trials, true);
  } else if (ok.size() != n_trials) {
    Rcpp::stop("Length of ok must equal number of trials");
  }

  R_xlen_t n_expand = expand.size();
  if (n_expand == 0) {
    expand = Rcpp::seq(1, n_trials);
    n_expand = expand.size();
  }
  for (R_xlen_t i = 0; i < n_expand; ++i) {
    int comp_idx = expand[i];
    if (comp_idx == NA_INTEGER || comp_idx < 1 || comp_idx > n_trials) {
      Rcpp::stop("expand indices must reference trials");
    }
  }

  std::vector<double> trial_loglik(static_cast<std::size_t>(n_trials), min_ll);

  for (int t = 0; t < n_trials; ++t) {
    if (t >= ok.size() || !ok[t]) {
      trial_loglik[static_cast<std::size_t>(t)] = min_ll;
      continue;
    }

    TrialParamSet *params_ptr = &param_sets[static_cast<std::size_t>(t)];
    double rt = (t < static_cast<int>(rt_by_trial.size()))
                    ? rt_by_trial[static_cast<std::size_t>(t)]
                    : std::numeric_limits<double>::quiet_NaN();
    Rcpp::String out_str = outcome_by_trial[static_cast<std::size_t>(t)];
    bool outcome_is_na = (out_str == NA_STRING);
    std::string outcome_label =
        outcome_is_na ? std::string() : std::string(out_str.get_cstring());
    const std::vector<std::string> *comp_labels_ptr = &default_components;
    const std::vector<ComponentCacheEntry> *cache_entries_ptr =
        &default_cache_entries;
    const std::vector<double> *comp_weights_ptr =
        &weights_by_trial[static_cast<std::size_t>(t)];

    std::string forced_component =
        (t < static_cast<int>(component_by_trial.size()))
            ? component_by_trial[static_cast<std::size_t>(t)]
            : std::string();
    std::vector<std::string> forced_components;
    std::vector<double> forced_weights;
    std::vector<ComponentCacheEntry> forced_cache_entries;
    if (!forced_component.empty()) {
      forced_components = {forced_component};
      forced_weights = {1.0};
      forced_cache_entries = build_component_cache_entries(forced_components);
      comp_labels_ptr = &forced_components;
      comp_weights_ptr = &forced_weights;
      cache_entries_ptr = &forced_cache_entries;
    }
    std::vector<int> comp_indices;
    comp_indices.reserve(comp_labels_ptr->size());
    for (const auto &comp_id : *comp_labels_ptr) {
      comp_indices.push_back(component_index_of(*ctx, comp_id));
    }

    TrialParamSet trigger_scratch;
    const SharedTriggerPlan *trigger_plan_ptr =
        (t < static_cast<int>(trigger_plans.size()))
            ? &trigger_plans[static_cast<std::size_t>(t)]
            : nullptr;

    auto compute_nonresponse_prob = [&]() -> double {
      double non_na_total = 0.0;
      for (std::size_t oi = 0; oi < ctx->outcome_info.size(); ++oi) {
        const std::string &lbl = (oi < ctx->outcome_labels.size())
                                     ? ctx->outcome_labels[oi]
                                     : std::string();
        if (lbl == "__DEADLINE__")
          continue;
        const uuber::OutcomeContextInfo &info = ctx->outcome_info[oi];
        double p =
            mix_outcome_mass(ctxSEXP, info, *comp_labels_ptr, *comp_weights_ptr,
                             params_ptr, rel_tol, abs_tol, max_depth, lbl);
        if (!std::isfinite(p) || p <= 0.0)
          continue;
        if (!info.maps_to_na)
          non_na_total += p;
      }
      return clamp_probability(1.0 - non_na_total);
    };

    double prob = 0.0;
    if (ranked_mode) {
      std::vector<std::string> ranked_labels;
      std::vector<double> ranked_times;
      ranked_labels.reserve(static_cast<std::size_t>(rank_width));
      ranked_times.reserve(static_cast<std::size_t>(rank_width));
      bool ranked_valid = true;
      bool seen_truncation = false;
      std::unordered_set<std::string> seen_labels;
      if (t >= static_cast<int>(ranked_outcome_by_trial.size()) ||
          t >= static_cast<int>(ranked_rt_by_trial.size())) {
        ranked_valid = false;
      } else {
        const std::vector<Rcpp::String> &trial_ranked_out =
            ranked_outcome_by_trial[static_cast<std::size_t>(t)];
        const std::vector<double> &trial_ranked_rt =
            ranked_rt_by_trial[static_cast<std::size_t>(t)];
        for (int rank_i = 0; rank_i < rank_width; ++rank_i) {
          Rcpp::String label = trial_ranked_out[static_cast<std::size_t>(rank_i)];
          double rank_rt = trial_ranked_rt[static_cast<std::size_t>(rank_i)];
          bool has_label = (label != NA_STRING);
          bool has_rt = std::isfinite(rank_rt);
          if (!has_label && !has_rt) {
            seen_truncation = true;
            continue;
          }
          if (seen_truncation || (has_label != has_rt)) {
            ranked_valid = false;
            break;
          }
          std::string label_str(label.get_cstring());
          if (label_str.empty() || rank_rt < 0.0) {
            ranked_valid = false;
            break;
          }
          if (!ranked_times.empty() && !(rank_rt > ranked_times.back())) {
            ranked_valid = false;
            break;
          }
          if (!seen_labels.insert(label_str).second) {
            ranked_valid = false;
            break;
          }
          ranked_labels.push_back(label_str);
          ranked_times.push_back(rank_rt);
        }
      }

      if (!ranked_valid) {
        prob = 0.0;
      } else if (ranked_labels.empty()) {
        prob = compute_nonresponse_prob();
      } else {
        double total = 0.0;
        for (std::size_t c = 0; c < comp_labels_ptr->size(); ++c) {
          double w =
              (c < comp_weights_ptr->size()) ? (*comp_weights_ptr)[c] : 0.0;
          if (!std::isfinite(w) || w <= 0.0) {
            continue;
          }
          const std::string &comp_label = (*comp_labels_ptr)[c];
          int comp_idx = comp_indices[c];
          std::string trial_type_key =
              (c < cache_entries_ptr->size())
                  ? (*cache_entries_ptr)[c].trial_type_key
                  : std::string();

          double keep_mult = 1.0;
          bool component_ok = true;
          for (const std::string &label : ranked_labels) {
            int rank_outcome_idx =
                resolve_outcome_index(*ctx, label, comp_label);
            if (rank_outcome_idx < 0 ||
                rank_outcome_idx >=
                    static_cast<int>(ctx->outcome_info.size())) {
              component_ok = false;
              break;
            }
            double keep_w =
                component_keep_weight(*ctx, comp_idx, rank_outcome_idx);
            if (!std::isfinite(keep_w) || keep_w <= 0.0) {
              component_ok = false;
              break;
            }
            keep_mult *= keep_w;
            if (!std::isfinite(keep_mult) || keep_mult <= 0.0) {
              component_ok = false;
              break;
            }
          }
          if (!component_ok) {
            continue;
          }

          TrialParamSet ranked_trigger_scratch;
          double contrib = ranked_prefix_density(
              *ctx, ranked_labels, ranked_times, comp_label, comp_idx,
              params_ptr, trial_type_key, trigger_plan_ptr,
              &ranked_trigger_scratch);
          if (std::isfinite(contrib) && contrib > 0.0) {
            total += w * keep_mult * contrib;
          }
        }
        prob = total;
      }
    } else if (outcome_is_na) {
      prob = compute_nonresponse_prob();
    } else {
      bool has_forced_component = !forced_component.empty();
      auto out_it = ctx->outcome_index.find(outcome_label);
      bool has_multi_defs =
          (out_it != ctx->outcome_index.end() && out_it->second.size() > 1);
      if (!std::isfinite(rt) || rt < 0.0) {
        prob = 0.0;
      } else if (has_multi_defs && !has_forced_component) {
        double total = 0.0;
        for (std::size_t c = 0; c < comp_labels_ptr->size(); ++c) {
          double w =
              (c < comp_weights_ptr->size()) ? (*comp_weights_ptr)[c] : 0.0;
          if (!std::isfinite(w) || w <= 0.0)
            continue;
          const std::string &comp_label = (*comp_labels_ptr)[c];
          int comp_outcome_idx =
              resolve_outcome_index(*ctx, outcome_label, comp_label);
          if (comp_outcome_idx < 0 ||
              comp_outcome_idx >=
                  static_cast<int>(ctx->outcome_info.size())) {
            continue;
          }
          const uuber::OutcomeContextInfo &info =
              ctx->outcome_info[static_cast<std::size_t>(comp_outcome_idx)];
          if (info.node_id < 0)
            continue;
          std::vector<std::string> comp_only = {comp_label};
          std::vector<double> weight_only = {1.0};
          std::vector<ComponentCacheEntry> cache_only;
          cache_only.reserve(1);
          if (c < cache_entries_ptr->size()) {
            cache_only.push_back((*cache_entries_ptr)[c]);
          } else {
            cache_only.push_back(default_component_cache_entry(comp_label));
          }
          double contrib = native_trial_mixture_internal(
              *ctx, info.node_id, rt, comp_only, weight_only, R_NilValue,
              info.competitor_ids, params_ptr, cache_only, outcome_label,
              Rcpp::List());
          if (std::isfinite(contrib) && contrib > 0.0) {
            total += w * contrib;
          }
        }
        prob = total;
      } else {
        int outcome_idx = resolve_outcome_index(
            *ctx, outcome_label, component_by_trial[static_cast<std::size_t>(t)]);
        if (outcome_idx < 0 ||
            outcome_idx >= static_cast<int>(ctx->outcome_info.size()) ||
            ctx->outcome_info[static_cast<std::size_t>(outcome_idx)].node_id <
                0) {
          prob = 0.0;
        } else {
          const uuber::OutcomeContextInfo &info =
              ctx->outcome_info[static_cast<std::size_t>(outcome_idx)];
          bool can_fast = true;
          if (!info.alias_sources.empty()) {
            can_fast = false;
          } else if (!info.guess_donors.empty()) {
            can_fast = false;
          } else {
            for (std::size_t i = 0; i < comp_labels_ptr->size(); ++i) {
              const std::string &comp_label = (*comp_labels_ptr)[i];
              if (!outcome_allows_component(info, comp_label))
                continue;
              int comp_idx = comp_indices[i];
              if (comp_idx >= 0 &&
                  comp_idx < static_cast<int>(ctx->component_info.size()) &&
                  ctx->component_info[static_cast<std::size_t>(comp_idx)]
                      .guess.valid) {
                can_fast = false;
                break;
              }
              double kw = component_keep_weight(*ctx, comp_idx, outcome_idx);
              if (std::fabs(kw - 1.0) > 1e-12) {
                can_fast = false;
                break;
              }
            }
          }
          if (can_fast) {
            static const std::unordered_set<int> kEmptyForced;
            double total = 0.0;
            for (std::size_t c = 0; c < comp_labels_ptr->size(); ++c) {
              double w =
                  (c < comp_weights_ptr->size()) ? (*comp_weights_ptr)[c] : 0.0;
              if (!std::isfinite(w) || w <= 0.0)
                continue;
              const std::string &comp_label = (*comp_labels_ptr)[c];
              if (!outcome_allows_component(info, comp_label))
                continue;
              int comp_idx = comp_indices[c];
              std::vector<int> comp_competitors;
              const std::vector<int> &comp_use = filter_competitor_ids(
                  *ctx, info.competitor_ids, comp_label, comp_competitors);
              double contrib = node_density_with_shared_triggers_plan(
                  *ctx, info.node_id, rt, comp_label, comp_idx,
                  kEmptyForced, kEmptyForced, comp_use, params_ptr,
                  std::string(), false, std::string(), -1, trigger_plan_ptr,
                  &trigger_scratch);
              if (std::isfinite(contrib) && contrib > 0.0) {
                total += w * contrib;
              }
            }
            prob = total;
          } else {
            std::vector<std::string> allowed_components;
            std::vector<double> allowed_weights;
            std::vector<ComponentCacheEntry> allowed_cache_entries;
            allowed_components.reserve(comp_labels_ptr->size());
            allowed_weights.reserve(comp_labels_ptr->size());
            allowed_cache_entries.reserve(comp_labels_ptr->size());
            for (std::size_t c = 0; c < comp_labels_ptr->size(); ++c) {
              const std::string &comp_label = (*comp_labels_ptr)[c];
              if (!outcome_allows_component(info, comp_label))
                continue;
              double w = (c < comp_weights_ptr->size())
                             ? (*comp_weights_ptr)[c]
                             : 0.0;
              if (!std::isfinite(w) || w <= 0.0)
                continue;
              allowed_components.push_back(comp_label);
              allowed_weights.push_back(w);
              if (c < cache_entries_ptr->size()) {
                allowed_cache_entries.push_back((*cache_entries_ptr)[c]);
              } else {
                allowed_cache_entries.push_back(
                    default_component_cache_entry(comp_label));
              }
            }
            if (allowed_components.empty()) {
              prob = 0.0;
            } else {
              prob = native_trial_mixture_internal(
                  *ctx, info.node_id, rt, allowed_components, allowed_weights,
                  R_NilValue, info.competitor_ids, params_ptr,
                  allowed_cache_entries, outcome_label, Rcpp::List());
            }
          }
        }
      }
    }

    double ll_val = min_ll;
    if (!std::isfinite(prob) || prob <= 0.0) {
      ll_val = min_ll;
    } else {
      double lp = std::log(prob);
      ll_val = (std::isfinite(lp) && lp > min_ll) ? lp : min_ll;
    }
    trial_loglik[static_cast<std::size_t>(t)] = ll_val;
  }

  double total_loglik = 0.0;
  for (R_xlen_t i = 0; i < n_expand; ++i) {
    int comp_idx = expand[i];
    total_loglik += trial_loglik[static_cast<std::size_t>(comp_idx - 1)];
  }

  return total_loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_loglik_multiple(SEXP ctxSEXP, Rcpp::List params_list,
                                        Rcpp::DataFrame data_df,
                                        Rcpp::LogicalVector ok,
                                        Rcpp::IntegerVector expand,
                                        double min_ll, double rel_tol,
                                        double abs_tol, int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  Rcpp::NumericVector out(params_list.size());
  for (R_xlen_t i = 0; i < params_list.size(); ++i) {
    if (params_list[i] == R_NilValue) {
      out[i] = NA_REAL;
      continue;
    }
    Rcpp::NumericMatrix pm(params_list[i]);
    out[i] = cpp_loglik(ctxSEXP, pm, data_df, ok, expand, min_ll, rel_tol,
                        abs_tol, max_depth);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List native_component_plan_exported(SEXP structureSEXP,
                                          SEXP trial_rowsSEXP,
                                          SEXP forced_componentSEXP) {
  Rcpp::List structure = Rcpp::as<Rcpp::List>(structureSEXP);
  Rcpp::DataFrame trial_rows_df;
  const Rcpp::DataFrame *trial_rows_ptr = nullptr;
  if (!Rf_isNull(trial_rowsSEXP)) {
    trial_rows_df = Rcpp::DataFrame(trial_rowsSEXP);
    trial_rows_ptr = &trial_rows_df;
  }
  double trial_id = extract_trial_id(trial_rows_ptr);
  std::string forced_component_value;
  const std::string *forced_component_ptr = nullptr;
  if (!Rf_isNull(forced_componentSEXP)) {
    forced_component_value = Rcpp::as<std::string>(forced_componentSEXP);
    forced_component_ptr = &forced_component_value;
  }
  return native_component_plan_impl(structure, trial_rows_ptr, trial_id,
                                    forced_component_ptr);
}

// [[Rcpp::export]]
double native_outcome_probability_params_cpp(
    SEXP ctxSEXP, int node_id, double upper,
    Rcpp::Nullable<Rcpp::String> component, SEXP forced_complete,
    SEXP forced_survive, Rcpp::IntegerVector competitor_ids, double rel_tol,
    double abs_tol, int max_depth, Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder =
      build_trial_params_from_df(*ctx, trial_rows);
  return native_outcome_probability_impl(
      ctxSEXP, node_id, upper, component, forced_complete, forced_survive,
      competitor_ids, rel_tol, abs_tol, max_depth,
      params_holder ? params_holder.get() : nullptr);
}

// [[Rcpp::export]]
Rcpp::DataFrame native_outcome_labels_cpp(SEXP ctxSEXP) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  R_xlen_t n = static_cast<R_xlen_t>(ctx->outcome_info.size());
  Rcpp::CharacterVector labels(n);
  Rcpp::IntegerVector node_ids(n);
  Rcpp::IntegerVector comp_counts(n);
  Rcpp::LogicalVector maps_na(n);
  for (R_xlen_t idx = 0; idx < n; ++idx) {
    const uuber::OutcomeContextInfo &info =
        ctx->outcome_info[static_cast<std::size_t>(idx)];
    if (idx < static_cast<R_xlen_t>(ctx->outcome_labels.size())) {
      labels[idx] = ctx->outcome_labels[static_cast<std::size_t>(idx)];
    } else {
      labels[idx] = "";
    }
    node_ids[idx] = info.node_id;
    comp_counts[idx] = static_cast<int>(info.competitor_ids.size());
    maps_na[idx] = info.maps_to_na;
  }
  return Rcpp::DataFrame::create(Rcpp::Named("label") = labels,
                                 Rcpp::Named("node_id") = node_ids,
                                 Rcpp::Named("competitors") = comp_counts,
                                 Rcpp::Named("maps_to_na") = maps_na,
                                 Rcpp::Named("stringsAsFactors") = false);
}

// Forward declarations for native distribution helpers
Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog);
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog);
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector &x, double shape,
                                   double rate);
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector &x, double shape,
                                   double rate);
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau);
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau);

namespace {

inline bool is_invalid_positive(double value) {
  return !std::isfinite(value) || value <= 0.0;
}

inline std::vector<double> expand_poly(const std::vector<double> &coeff,
                                       double surv, double fail) {
  std::size_t len = coeff.size();
  std::vector<double> out(len + 1, 0.0);
  for (std::size_t i = 0; i < len; ++i) {
    double base = coeff[i];
    if (base == 0.0)
      continue;
    out[i] += base * surv;
    out[i + 1] += base * fail;
  }
  return out;
}

struct PoolPolyScratch {
  std::vector<std::vector<double>> prefix;
  std::vector<std::vector<double>> suffix;
  std::vector<double> surv;
  std::vector<double> fail;
};

inline void ensure_prefix_shape(std::vector<std::vector<double>> &prefix,
                                std::size_t n) {
  if (prefix.size() != n + 1) {
    prefix.resize(n + 1);
  }
  for (std::size_t i = 0; i <= n; ++i) {
    if (prefix[i].size() != i + 1) {
      prefix[i].assign(i + 1, 0.0);
    }
  }
}

inline void ensure_suffix_shape(std::vector<std::vector<double>> &suffix,
                                std::size_t n) {
  if (suffix.size() != n + 1) {
    suffix.resize(n + 1);
  }
  for (std::size_t i = 0; i <= n; ++i) {
    std::size_t len = (n - i) + 1;
    if (suffix[i].size() != len) {
      suffix[i].assign(len, 0.0);
    }
  }
}

inline PoolPolyScratch &fetch_pool_poly_scratch(std::size_t n,
                                                bool need_suffix) {
  thread_local std::unordered_map<std::size_t, PoolPolyScratch> scratch_map;
  PoolPolyScratch &scratch = scratch_map[n];
  ensure_prefix_shape(scratch.prefix, n);
  if (need_suffix) {
    ensure_suffix_shape(scratch.suffix, n);
  }
  return scratch;
}

inline void ensure_prob_buffers(PoolPolyScratch &scratch, std::size_t n) {
  if (scratch.surv.size() != n) {
    scratch.surv.resize(n);
  }
  if (scratch.fail.size() != n) {
    scratch.fail.resize(n);
  }
}

inline void fill_prefix_buffers(std::vector<std::vector<double>> &prefix,
                                const std::vector<double> &surv,
                                const std::vector<double> &fail) {
  const std::size_t n = surv.size();
  std::fill(prefix[0].begin(), prefix[0].end(), 0.0);
  prefix[0][0] = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::vector<double> &prev = prefix[i];
    std::vector<double> &out = prefix[i + 1];
    std::fill(out.begin(), out.end(), 0.0);
    double surv_i = surv[i];
    double fail_i = fail[i];
    for (std::size_t j = 0; j < prev.size(); ++j) {
      double base = prev[j];
      if (base == 0.0)
        continue;
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline void fill_suffix_buffers(std::vector<std::vector<double>> &suffix,
                                const std::vector<double> &surv,
                                const std::vector<double> &fail) {
  const std::size_t n = surv.size();
  std::fill(suffix[n].begin(), suffix[n].end(), 0.0);
  suffix[n][0] = 1.0;
  for (std::size_t idx = n; idx-- > 0;) {
    const std::vector<double> &next = suffix[idx + 1];
    std::vector<double> &out = suffix[idx];
    std::fill(out.begin(), out.end(), 0.0);
    double surv_i = surv[idx];
    double fail_i = fail[idx];
    for (std::size_t j = 0; j < next.size(); ++j) {
      double base = next[j];
      if (base == 0.0)
        continue;
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline double coefficient_for_order(const std::vector<double> &pref,
                                    const std::vector<double> &suff,
                                    int order) {
  if (order < 0)
    return 0.0;
  double total = 0.0;
  int max_pref = static_cast<int>(pref.size()) - 1;
  int max_index = std::min(order, max_pref);
  for (int i = 0; i <= max_index; ++i) {
    int s_idx = order - i;
    if (s_idx < 0 || s_idx >= static_cast<int>(suff.size()))
      continue;
    total += pref[static_cast<std::size_t>(i)] *
             suff[static_cast<std::size_t>(s_idx)];
  }
  return total;
}

struct PoolTemplate {
  int finisher_idx;
  std::vector<int> complete_idx;
  std::vector<int> survivor_idx;
  std::vector<int> forced_complete_ids;
  std::vector<int> forced_survive_ids;
};

inline void sort_unique(std::vector<int> &vec) {
  if (vec.empty())
    return;
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

inline std::vector<int> integer_vector_to_std(const Rcpp::IntegerVector &src,
                                              bool subtract_one) {
  std::vector<int> out;
  out.reserve(src.size());
  for (int val : src) {
    if (val == NA_INTEGER)
      continue;
    out.push_back(subtract_one ? val - 1 : val);
  }
  sort_unique(out);
  return out;
}

inline std::vector<int> merge_forced_vectors(const std::vector<int> &base,
                                             const std::vector<int> &addition) {
  std::vector<int> result = base;
  for (int val : addition) {
    if (val == NA_INTEGER)
      continue;
    result.push_back(val);
  }
  sort_unique(result);
  return result;
}

void combinations_recursive(const std::vector<int> &elements, int choose,
                            std::size_t start, std::vector<int> &current,
                            std::vector<std::vector<int>> &output) {
  if (static_cast<int>(current.size()) == choose) {
    output.push_back(current);
    return;
  }
  for (std::size_t i = start; i < elements.size(); ++i) {
    current.push_back(elements[i]);
    combinations_recursive(elements, choose, i + 1, current, output);
    current.pop_back();
  }
}

inline std::vector<std::vector<int>>
generate_combinations(const std::vector<int> &elements, int choose) {
  std::vector<std::vector<int>> combos;
  if (choose <= 0) {
    combos.emplace_back();
    return combos;
  }
  if (choose >= static_cast<int>(elements.size())) {
    combos.push_back(elements);
    return combos;
  }
  if (choose == 1) {
    combos.reserve(elements.size());
    for (int el : elements) {
      combos.push_back(std::vector<int>{el});
    }
    return combos;
  }
  std::vector<int> current;
  combinations_recursive(elements, choose, 0, current, combos);
  return combos;
}

inline std::vector<int> survivors_from_combo(const std::vector<int> &others,
                                             const std::vector<int> &combo) {
  if (combo.empty())
    return others;
  std::vector<int> survivors;
  std::vector<int> combo_sorted = combo;
  std::sort(combo_sorted.begin(), combo_sorted.end());
  survivors.reserve(others.size());
  for (int candidate : others) {
    if (!std::binary_search(combo_sorted.begin(), combo_sorted.end(),
                            candidate)) {
      survivors.push_back(candidate);
    }
  }
  return survivors;
}

} // namespace

// ------------------------------------------------------------------
// Lognormal
// ------------------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sdlog)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!R_FINITE(xi)) {
      out[i] = 0.0;
      continue;
    }
    out[i] = R::dlnorm(xi, meanlog, sdlog, /*give_log =*/0);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sdlog)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i];
    if (Rcpp::NumericVector::is_na(xi)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!R_FINITE(xi)) {
      out[i] = xi < 0.0 ? 0.0 : 1.0;
      continue;
    }
    out[i] = R::plnorm(xi, meanlog, sdlog, /*lower_tail =*/1, /*log_p =*/0);
    out[i] = clamp(out[i], 0.0, 1.0);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lognormal_rng(int n, double meanlog, double sdlog) {
  if (n <= 0 || is_invalid_positive(sdlog)) {
    return Rcpp::NumericVector();
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector draws = Rcpp::rnorm(n, meanlog, sdlog);
  for (R_xlen_t i = 0; i < draws.size(); ++i) {
    draws[i] = std::exp(draws[i]);
  }
  return draws;
}

// ------------------------------------------------------------------
// Gamma
// ------------------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector &x, double shape,
                                   double rate) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = gamma_pdf_fast(x[i], shape, rate);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector &x, double shape,
                                   double rate) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = gamma_cdf_fast(x[i], shape, rate);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_rng(int n, double shape, double rate) {
  if (n <= 0 || is_invalid_positive(shape) || is_invalid_positive(rate)) {
    return Rcpp::NumericVector();
  }
  Rcpp::RNGScope scope;
  return Rcpp::rgamma(n, shape, 1.0 / rate);
}

// ------------------------------------------------------------------
// Ex-Gaussian (Normal + Exponential)
// ------------------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = exgauss_pdf_fast(x[i], mu, sigma, tau);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = exgauss_cdf_fast(x[i], mu, sigma, tau);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_rng(int n, double mu, double sigma,
                                     double tau) {
  if (n <= 0 || is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    return Rcpp::NumericVector();
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector normals = Rcpp::rnorm(n, mu, sigma);
  Rcpp::NumericVector expo = Rcpp::rexp(n, 1.0 / tau);
  for (R_xlen_t i = 0; i < normals.size(); ++i) {
    normals[i] += expo[i];
  }
  return normals;
}

// ------------------------------------------------------------------
// Pool helpers
// ------------------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector pool_coeffs_cpp(const Rcpp::NumericVector &Svec,
                                    const Rcpp::NumericVector &Fvec) {
  std::size_t n = Svec.size();
  std::vector<double> coeff{1.0};
  coeff.reserve(n + 1);
  for (std::size_t i = 0; i < n; ++i) {
    double surv = clamp_unit(Svec[i]);
    double fail = clamp_unit(Fvec[i]);
    coeff = expand_poly(coeff, surv, fail);
  }
  Rcpp::NumericVector out(coeff.size());
  for (std::size_t i = 0; i < coeff.size(); ++i) {
    out[i] = coeff[i];
  }
  return out;
}

static inline void fill_surv_fail_buffers(PoolPolyScratch &scratch,
                                          const double *survival,
                                          std::size_t n) {
  ensure_prob_buffers(scratch, n);
  std::vector<double> &surv = scratch.surv;
  std::vector<double> &fail = scratch.fail;
  for (std::size_t i = 0; i < n; ++i) {
    double s = clamp_unit(survival[i]);
    surv[i] = s;
    fail[i] = clamp_unit(1.0 - s);
  }
}

static double pool_density_fast_impl(const double *density,
                                     const double *survival, std::size_t n,
                                     int k) {
  if (n == 0)
    return 0.0;
  if (k < 1)
    return 0.0;
  if (k > static_cast<int>(n))
    return 0.0;

  PoolPolyScratch &scratch = fetch_pool_poly_scratch(n, true);
  fill_surv_fail_buffers(scratch, survival, n);
  fill_prefix_buffers(scratch.prefix, scratch.surv, scratch.fail);
  fill_suffix_buffers(scratch.suffix, scratch.surv, scratch.fail);

  double total = 0.0;
  for (std::size_t idx = 0; idx < n; ++idx) {
    double dens_val = density[idx];
    if (!std::isfinite(dens_val) || dens_val <= 0.0)
      continue;
    double coeff =
        coefficient_for_order(scratch.prefix[idx], scratch.suffix[idx + 1],
                               k - 1);
    if (!std::isfinite(coeff) || coeff <= 0.0)
      continue;
    total += dens_val * coeff;
  }
  if (!std::isfinite(total) || total < 0.0)
    return 0.0;
  return total;
}

static double pool_survival_fast_impl(const double *survival, std::size_t n,
                                      int k) {
  if (n == 0)
    return 1.0;
  if (k < 1)
    return 0.0;

  PoolPolyScratch &scratch = fetch_pool_poly_scratch(n, false);
  fill_surv_fail_buffers(scratch, survival, n);
  fill_prefix_buffers(scratch.prefix, scratch.surv, scratch.fail);
  const std::vector<double> &poly = scratch.prefix.back();
  int upto = std::min<int>(k, static_cast<int>(poly.size()));
  if (upto <= 0)
    return 0.0;
  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += poly[static_cast<std::size_t>(i)];
  }
  return clamp(total, 0.0, 1.0);
}

double pool_density_fast(const std::vector<double> &density,
                         const std::vector<double> &survival, int k) {
  if (density.size() != survival.size())
    return 0.0;
  return pool_density_fast_impl(density.data(), survival.data(), density.size(),
                                k);
}

double pool_survival_fast(const std::vector<double> &survival, int k) {
  return pool_survival_fast_impl(survival.data(), survival.size(), k);
}

double pool_density_fast_cpp(const Rcpp::NumericVector &density,
                             const Rcpp::NumericVector &survival, int k) {
  if (density.size() != survival.size())
    return 0.0;
  return pool_density_fast_impl(density.begin(), survival.begin(),
                                density.size(), k);
}

double pool_survival_fast_cpp(const Rcpp::NumericVector &survival, int k) {
  return pool_survival_fast_impl(survival.begin(), survival.size(), k);
}

Rcpp::List pool_build_templates_cpp(int n,
                                    const Rcpp::IntegerVector &member_ids,
                                    int pool_idx, int k) {
  std::vector<int> member_std(member_ids.begin(), member_ids.end());
  uuber::PoolTemplateCacheEntry cache =
      build_pool_template_cache(n, member_std, pool_idx, k);
  Rcpp::List finisher_map(n);
  Rcpp::List templates_out(cache.templates.size());
  for (std::size_t i = 0; i < cache.templates.size(); ++i) {
    const uuber::PoolTemplateEntry &tpl = cache.templates[i];
    Rcpp::IntegerVector complete_idx;
    if (!tpl.complete_idx.empty()) {
      complete_idx =
          Rcpp::IntegerVector(tpl.complete_idx.begin(), tpl.complete_idx.end());
      for (auto &val : complete_idx) {
        val += 1;
      }
    } else {
      complete_idx = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector survivor_idx;
    if (!tpl.survivor_idx.empty()) {
      survivor_idx =
          Rcpp::IntegerVector(tpl.survivor_idx.begin(), tpl.survivor_idx.end());
      for (auto &val : survivor_idx) {
        val += 1;
      }
    } else {
      survivor_idx = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector forced_complete_vec;
    if (!tpl.forced_complete_ids.empty()) {
      forced_complete_vec = Rcpp::IntegerVector(tpl.forced_complete_ids.begin(),
                                                tpl.forced_complete_ids.end());
    } else {
      forced_complete_vec = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector forced_survive_vec;
    if (!tpl.forced_survive_ids.empty()) {
      forced_survive_vec = Rcpp::IntegerVector(tpl.forced_survive_ids.begin(),
                                               tpl.forced_survive_ids.end());
    } else {
      forced_survive_vec = Rcpp::IntegerVector(0);
    }
    templates_out[i] = Rcpp::List::create(
        Rcpp::Named("finisher_idx") = tpl.finisher_idx + 1,
        Rcpp::Named("complete_idx") = complete_idx,
        Rcpp::Named("survivor_idx") = survivor_idx,
        Rcpp::Named("forced_complete_ids") = forced_complete_vec,
        Rcpp::Named("forced_survive_ids") = forced_survive_vec);
  }
  for (int idx = 0; idx < n; ++idx) {
    const std::vector<int> &entries =
        cache.finisher_map[static_cast<std::size_t>(idx)];
    if (entries.empty()) {
      finisher_map[idx] = Rcpp::IntegerVector(0);
      continue;
    }
    Rcpp::IntegerVector fm(entries.size());
    for (std::size_t j = 0; j < entries.size(); ++j) {
      fm[j] = entries[j] + 1;
    }
    finisher_map[idx] = fm;
  }
  return Rcpp::List::create(Rcpp::Named("templates") = templates_out,
                            Rcpp::Named("finisher_map") = finisher_map);
}

Rcpp::List pool_density_combine_native(
    const Rcpp::NumericVector &dens_vec, const Rcpp::NumericVector &cdf_vec,
    const Rcpp::NumericVector &surv_vec,
    const Rcpp::NumericVector &cdf_success_vec,
    const Rcpp::NumericVector &surv_success_vec,
    const std::vector<int> &shared_index,
    const std::vector<uuber::PoolTemplateEntry> &templates,
    const std::vector<int> &base_forced_complete,
    const std::vector<int> &base_forced_survive) {
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<Rcpp::List> scenario_vec;
  scenario_vec.reserve(templates.size());
  double total = 0.0;

  auto same_shared = [&](int i, int j) {
    if (i < 0 || j < 0 || i >= static_cast<int>(shared_index.size()) ||
        j >= static_cast<int>(shared_index.size())) {
      return false;
    }
    int si = shared_index[static_cast<std::size_t>(i)];
    int sj = shared_index[static_cast<std::size_t>(j)];
    return si > 0 && sj > 0 && si == sj;
  };

  for (const auto &tpl : templates) {
    int finisher_idx = tpl.finisher_idx;
    if (finisher_idx < 0 || finisher_idx >= dens_vec.size())
      continue;
    double dens_mid = dens_vec[finisher_idx];
    if (!std::isfinite(dens_mid) || dens_mid <= 0.0)
      continue;
    double weight = dens_mid;

    for (int j : tpl.complete_idx) {
      if (j < 0 || j >= cdf_vec.size())
        continue;
      weight *= cdf_vec[j];
    }
    for (int j : tpl.survivor_idx) {
      if (j < 0 || j >= surv_vec.size())
        continue;
      weight *= surv_vec[j];
    }
    if (!std::isfinite(weight) || weight <= 0.0)
      continue;

    for (int j : tpl.complete_idx) {
      if (j < 0 || j >= cdf_vec.size())
        continue;
      if (!same_shared(finisher_idx, j))
        continue;
      double denom = std::max(cdf_vec[j], eps);
      double ratio = cdf_success_vec[j] / denom;
      weight *= ratio;
    }
    for (int j : tpl.survivor_idx) {
      if (j < 0 || j >= surv_vec.size())
        continue;
      if (!same_shared(finisher_idx, j))
        continue;
      double denom = std::max(surv_vec[j], eps);
      double ratio = surv_success_vec[j] / denom;
      weight *= ratio;
    }
    if (!std::isfinite(weight) || weight <= 0.0)
      continue;

    std::vector<int> merged_fc =
        merge_forced_vectors(base_forced_complete, tpl.forced_complete_ids);
    std::vector<int> merged_fs =
        merge_forced_vectors(base_forced_survive, tpl.forced_survive_ids);

    total += weight;

    Rcpp::IntegerVector fc_vec;
    if (!merged_fc.empty()) {
      fc_vec = Rcpp::IntegerVector(merged_fc.begin(), merged_fc.end());
    } else {
      fc_vec = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector fs_vec;
    if (!merged_fs.empty()) {
      fs_vec = Rcpp::IntegerVector(merged_fs.begin(), merged_fs.end());
    } else {
      fs_vec = Rcpp::IntegerVector(0);
    }

    scenario_vec.emplace_back(Rcpp::List::create(
        Rcpp::Named("weight") = weight, Rcpp::Named("forced_complete") = fc_vec,
        Rcpp::Named("forced_survive") = fs_vec));
  }

  Rcpp::List scenarios_out(scenario_vec.size());
  for (std::size_t i = 0; i < scenario_vec.size(); ++i) {
    scenarios_out[i] = scenario_vec[i];
  }

  return Rcpp::List::create(Rcpp::Named("value") = total,
                            Rcpp::Named("scenarios") = scenarios_out);
}

Rcpp::List pool_density_combine_cpp(const Rcpp::NumericVector &dens_vec,
                                    const Rcpp::NumericVector &cdf_vec,
                                    const Rcpp::NumericVector &surv_vec,
                                    const Rcpp::NumericVector &cdf_success_vec,
                                    const Rcpp::NumericVector &surv_success_vec,
                                    const Rcpp::IntegerVector &shared_index,
                                    const Rcpp::List &templates,
                                    const Rcpp::IntegerVector &forced_complete,
                                    const Rcpp::IntegerVector &forced_survive) {
  std::vector<int> shared_group(shared_index.begin(), shared_index.end());
  std::vector<int> base_fc = integer_vector_to_std(forced_complete, false);
  std::vector<int> base_fs = integer_vector_to_std(forced_survive, false);
  std::vector<uuber::PoolTemplateEntry> native_templates;
  native_templates.reserve(templates.size());
  for (R_xlen_t i = 0; i < templates.size(); ++i) {
    Rcpp::List tpl(templates[i]);
    uuber::PoolTemplateEntry entry;
    entry.finisher_idx = tpl["finisher_idx"];
    entry.finisher_idx -= 1;
    entry.complete_idx = integer_vector_to_std(tpl["complete_idx"], true);
    entry.survivor_idx = integer_vector_to_std(tpl["survivor_idx"], true);
    entry.forced_complete_ids =
        integer_vector_to_std(tpl["forced_complete_ids"], false);
    entry.forced_survive_ids =
        integer_vector_to_std(tpl["forced_survive_ids"], false);
    native_templates.push_back(std::move(entry));
  }
  return pool_density_combine_native(
      dens_vec, cdf_vec, surv_vec, cdf_success_vec, surv_success_vec,
      shared_group, native_templates, base_fc, base_fs);
}

double pool_survival_general_cpp(const Rcpp::NumericVector &Fvec, int k) {
  const std::size_t n = Fvec.size();
  if (n == 0)
    return 1.0;
  if (k < 1)
    return 0.0;

  Rcpp::NumericVector Fclamp(n);
  Rcpp::NumericVector Svec(n);
  for (std::size_t i = 0; i < n; ++i) {
    double val = Fvec[i];
    if (!std::isfinite(val))
      val = 0.0;
    val = clamp(val, 0.0, 1.0);
    Fclamp[i] = val;
    Svec[i] = 1.0 - val;
  }

  Rcpp::NumericVector coeffs = pool_coeffs_cpp(Svec, Fclamp);
  int upto = std::min<int>(coeffs.size(), k);
  if (upto <= 0)
    return 0.0;
  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += coeffs[i];
  }
  return clamp(total, 0.0, 1.0);
}

double guard_effective_survival_cpp(Rcpp::Function integrand, double upper,
                                    double rel_tol, double abs_tol,
                                    int max_depth) {
  if (!std::isfinite(upper))
    return 0.0;
  if (upper <= 0.0)
    return 1.0;
  if (rel_tol <= 0.0)
    rel_tol = 1e-5;
  if (abs_tol <= 0.0)
    abs_tol = 1e-6;
  double integral = 0.0;
  try {
    integral = uuber::integrate_boost(integrand, 0.0, upper, rel_tol, abs_tol,
                                      max_depth);
  } catch (...) {
    integral = 0.0;
  }
  if (!std::isfinite(integral) || integral < 0.0)
    integral = 0.0;
  double surv = 1.0 - integral;
  if (!std::isfinite(surv))
    surv = 0.0;
  if (surv < 0.0)
    surv = 0.0;
  if (surv > 1.0)
    surv = 1.0;
  return surv;
}
