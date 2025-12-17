// [[Rcpp::depends(Rcpp, RcppParallel, BH)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>
#include <vector>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <memory>
#include <queue>
#include <deque>
#include <sstream>
#include <numeric>
#include <iomanip>
#include <cstring>
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#define UUBER_HAVE_BOOST_GK 1
#include <boost/math/quadrature/gauss_kronrod.hpp>
#else
#define UUBER_HAVE_BOOST_GK 0
#endif

#include "context.h"
#include "accumulator.h"
#include "prep_builder.h"
#include "proto.h"
#include "integrate.h"

using uuber::AccDistParams;
using uuber::resolve_acc_params_entries;

Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog);
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog);
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate);
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate);
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau);
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau);

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
    reinterpret_cast<const std::uint8_t*>(blob.begin()),
    static_cast<std::size_t>(blob.size())
  );
  std::unique_ptr<uuber::NativeContext> ctx = uuber::build_context_from_proto(proto);
  return Rcpp::XPtr<uuber::NativeContext>(ctx.release(), true);
}

// [[Rcpp::export]]
double boost_integrate_cpp(Rcpp::Function integrand,
                           double lower,
                           double upper,
                           double rel_tol,
                           double abs_tol,
                           int max_depth) {
  return uuber::integrate_boost(integrand, lower, upper, rel_tol, abs_tol, max_depth);
}

double acc_density_cpp(double, double, double, const std::string&, const Rcpp::List&);
double acc_survival_cpp(double, double, double, const std::string&, const Rcpp::List&);
double acc_cdf_success_cpp(double, double, double, const std::string&, const Rcpp::List&);
double pool_density_fast_cpp(const Rcpp::NumericVector&, const Rcpp::NumericVector&, int);
double pool_survival_fast_cpp(const Rcpp::NumericVector&, int);
Rcpp::List pool_density_combine_native(const Rcpp::NumericVector& dens_vec,
                                       const Rcpp::NumericVector& cdf_vec,
                                       const Rcpp::NumericVector& surv_vec,
                                       const Rcpp::NumericVector& cdf_success_vec,
                                       const Rcpp::NumericVector& surv_success_vec,
                                       const std::vector<int>& shared_index,
                                       const std::vector<uuber::PoolTemplateEntry>& templates,
                                       const std::vector<int>& base_forced_complete,
                                       const std::vector<int>& base_forced_survive);
std::vector<uuber::ProtoParamEntry> params_from_rcpp(const Rcpp::List& params,
                                                     const std::string& dist) {
  Rcpp::CharacterVector names = params.names();
  std::vector<uuber::ProtoParamEntry> out;
  out.reserve(params.size());
  for (R_xlen_t i = 0; i < params.size(); ++i) {
    if (params[i] == R_NilValue) continue;
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
                                    const Rcpp::IntegerVector& member_ids,
                                    int pool_idx,
                                    int k);
Rcpp::List native_component_plan_impl(const Rcpp::List& structure,
                                      const Rcpp::DataFrame* trial_rows,
                                      double trial_id,
                                      const std::string* forced_component);
Rcpp::List pool_density_combine_cpp(const Rcpp::NumericVector& dens_vec,
                                    const Rcpp::NumericVector& cdf_vec,
                                    const Rcpp::NumericVector& surv_vec,
                                    const Rcpp::NumericVector& cdf_success_vec,
                                    const Rcpp::NumericVector& surv_success_vec,
                                    const Rcpp::IntegerVector& shared_index,
                                    const Rcpp::List& templates,
                                    const Rcpp::IntegerVector& forced_complete,
                                    const Rcpp::IntegerVector& forced_survive);

namespace {

std::vector<std::vector<int>> generate_combinations(const std::vector<int>& elements,
                                                    int choose);
std::vector<int> survivors_from_combo(const std::vector<int>& others,
                                      const std::vector<int>& combo);
uuber::PoolTemplateCacheEntry build_pool_template_cache(int n,
                                                       const std::vector<int>& member_ids,
                                                       int pool_idx,
                                                       int k);

void sort_unique(std::vector<int>& vec);
std::vector<int> integer_vector_to_std(const Rcpp::IntegerVector& src,
                                       bool subtract_one = false);

constexpr double kDefaultRelTol = 1e-5;
constexpr double kDefaultAbsTol = 1e-6;
constexpr int kDefaultMaxDepth = 12;
template <typename T>
inline T clamp(T val, T lo, T hi) {
  if (!std::isfinite(val)) {
    return val;
  }
  if (val < lo) return lo;
  if (val > hi) return hi;
  return val;
}

inline double clamp_unit(double val) {
  if (!std::isfinite(val)) return 0.0;
  if (val < 0.0) return 0.0;
  if (val > 1.0) return 1.0;
  return val;
}

inline double clamp_probability(double value) {
  if (!std::isfinite(value)) return 0.0;
  if (value < 0.0) return 0.0;
  if (value > 1.0) return 1.0;
  return value;
}

inline double safe_density(double value) {
  if (!std::isfinite(value) || value <= 0.0) return 0.0;
  return value;
}

inline double lognormal_pdf_fast(double x, double meanlog, double sdlog) {
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) || !std::isfinite(sdlog) || sdlog <= 0.0) {
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
  if (!std::isfinite(x) || x <= 0.0 || !std::isfinite(meanlog) || !std::isfinite(sdlog) || sdlog <= 0.0) {
    return 0.0;
  }
  const double inv_sigma = 1.0 / sdlog;
  const double logx = std::log(x);
  const double z = (logx - meanlog) * inv_sigma;
  // 0.5 * erfc(-z / sqrt(2))
  const double arg = -z * 0.7071067811865475;
  double val = 0.5 * std::erfc(arg);
  if (!std::isfinite(val)) return 0.0;
  if (val < 0.0) return 0.0;
  if (val > 1.0) return 1.0;
  return val;
}

inline double eval_pdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_pdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return dist_gamma_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_EXGAUSS:
    return dist_exgauss_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    return 0.0;
  }
}

inline double eval_cdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return lognormal_cdf_fast(x, cfg.p1, cfg.p2);
  case uuber::ACC_DIST_GAMMA:
    return dist_gamma_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_EXGAUSS:
    return dist_exgauss_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    return 0.0;
  }
}

inline double total_onset_with_t0(double onset,
                                  const AccDistParams& cfg) {
  return onset + cfg.t0;
}

inline double acc_density_from_cfg(double t,
                                   double onset,
                                   double q,
                                   const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t) || t < 0.0) return 0.0;
  if (t < effective_onset) return 0.0;
  double success_prob = 1.0 - q;
  if (success_prob <= 0.0) return 0.0;
  double dens = eval_pdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(dens) || !std::isfinite(dens)) {
    return NA_REAL;
  }
  return success_prob * dens;
}

inline double acc_survival_from_cfg(double t,
                                    double onset,
                                    double q,
                                    const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t)) return 0.0;
  if (t < 0.0) return 1.0;
  if (t < effective_onset) return 1.0;
  double cdf = eval_cdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  double surv_underlying = clamp(1.0 - cdf, 0.0, 1.0);
  double success_prob = 1.0 - q;
  double result = q + success_prob * surv_underlying;
  return clamp(result, 0.0, 1.0);
}

inline double acc_cdf_success_from_cfg(double t,
                                       double onset,
                                       const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t)) return 1.0;
  if (t < 0.0) return 0.0;
  if (t < effective_onset) return 0.0;
  double cdf = eval_cdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  return clamp(cdf, 0.0, 1.0);
}

struct NodeEvalResult {
  double density{0.0};
  double survival{1.0};
  double cdf{0.0};
};

inline NodeEvalResult make_node_result(double density,
                                       double survival,
                                       double cdf) {
  NodeEvalResult out;
  out.density = safe_density(density);
  out.survival = clamp_probability(survival);
  out.cdf = clamp_probability(cdf);
  return out;
}

inline std::unordered_set<int> make_scope_set(const std::vector<int>& ids) {
  std::unordered_set<int> scope;
  scope.reserve(ids.size());
  for (int id : ids) {
    if (id == NA_INTEGER) continue;
    scope.insert(id);
  }
  return scope;
}

inline std::vector<int> set_to_sorted_vector(const std::unordered_set<int>& items) {
  if (items.empty()) return {};
  std::vector<int> out(items.begin(), items.end());
  sort_unique(out);
  return out;
}

inline std::vector<int> union_vectors(const std::vector<int>& a,
                                      const std::vector<int>& b) {
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

inline bool vector_contains(const std::vector<int>& values, int target) {
  if (values.empty()) return false;
  return std::find(values.begin(), values.end(), target) != values.end();
}

struct ScenarioRecord {
  double weight{0.0};
  std::vector<int> forced_complete;
  std::vector<int> forced_survive;
};

inline ScenarioRecord make_scenario_record(double weight,
                                           const std::vector<int>& forced_complete,
                                           const std::vector<int>& forced_survive) {
  ScenarioRecord rec;
  rec.weight = weight;
  rec.forced_complete = forced_complete;
  rec.forced_survive = forced_survive;
  return rec;
}

std::vector<ScenarioRecord> rcpp_scenarios_to_records(const Rcpp::List& scenarios) {
  std::vector<ScenarioRecord> out;
  out.reserve(scenarios.size());
  for (R_xlen_t i = 0; i < scenarios.size(); ++i) {
    if (scenarios[i] == R_NilValue) continue;
    Rcpp::List sc(scenarios[i]);
    double weight = Rcpp::as<double>(sc["weight"]);
    if (!std::isfinite(weight) || weight <= 0.0) continue;
    std::vector<int> fc = integer_vector_to_std(sc["forced_complete"], false);
    std::vector<int> fs = integer_vector_to_std(sc["forced_survive"], false);
    out.push_back(make_scenario_record(weight, fc, fs));
  }
  return out;
}

std::unordered_set<int> filter_forced_scope(const std::unordered_set<int>& forced,
                                            const std::vector<int>& scope_ids) {
  if (forced.empty() || scope_ids.empty()) {
    return std::unordered_set<int>();
  }
  std::unordered_set<int> scope = make_scope_set(scope_ids);
  if (scope.empty()) {
    return std::unordered_set<int>();
  }
  std::unordered_set<int> out;
  out.reserve(scope.size());
  for (int id : forced) {
    if (scope.count(id)) {
      out.insert(id);
    }
  }
  return out;
}

std::vector<int> forced_vec_from_sexp(SEXP vec) {
  std::vector<int> out;
  if (Rf_isNull(vec)) return out;
  Rcpp::IntegerVector iv(vec);
  out.reserve(iv.size());
  for (int val : iv) {
    if (val == NA_INTEGER) continue;
    out.push_back(val);
  }
  return out;
}

std::unordered_set<int> make_forced_set(const std::vector<int>& ids) {
  return std::unordered_set<int>(ids.begin(), ids.end());
}

struct TrialParamSet;

inline int label_id(const uuber::NativeContext& ctx, const std::string& label) {
  auto it = ctx.label_to_id.find(label);
  if (it == ctx.label_to_id.end()) return NA_INTEGER;
  return it->second;
}

inline std::vector<int> ensure_source_ids(const uuber::NativeContext& ctx,
                                          const uuber::NativeNode& node) {
  std::vector<int> ids = node.source_ids;
  if (ids.empty() && !node.source.empty()) {
    int label_idx = label_id(ctx, node.source);
    if (label_idx != NA_INTEGER) {
      ids.push_back(label_idx);
    }
  }
  sort_unique(ids);
  return ids;
}

struct TrialAccumulatorParams {
  double onset{0.0};
  double q{0.0};
  double shared_q{std::numeric_limits<double>::quiet_NaN()};
  AccDistParams dist_cfg{};
  // Flags/components unused in flat fast path but kept for compatibility with component checks.
  bool has_components{false};
  std::vector<std::string> components;
  std::string shared_trigger_id;
  bool has_override{false};
};

struct TrialParamSet {
  std::vector<TrialAccumulatorParams> acc_params;
};

inline TrialAccumulatorParams base_params(const uuber::NativeAccumulator& base) {
  TrialAccumulatorParams p;
  p.onset = base.onset;
  p.q = base.q;
  p.dist_cfg = base.dist_cfg;
  p.shared_q = base.q;
  p.shared_trigger_id = base.shared_trigger_id;
  p.has_override = false;
  return p;
}

// Build a full parameter set from base context when no trial-level overrides are supplied.
inline TrialParamSet build_base_paramset(const uuber::NativeContext& ctx) {
  TrialParamSet ps;
  ps.acc_params.reserve(ctx.accumulators.size());
  for (const auto& acc : ctx.accumulators) {
    ps.acc_params.push_back(base_params(acc));
  }
  return ps;
}

struct SharedTriggerInfo {
  std::string id;
  std::vector<int> acc_indices;
  double q{0.0};
};

// Collect shared trigger definitions for the current trial parameters.
inline std::vector<SharedTriggerInfo> collect_shared_triggers(
  const uuber::NativeContext& ctx,
  const TrialParamSet& params) {
  std::vector<SharedTriggerInfo> out;
  if (!ctx.shared_trigger_map.empty()) {
    for (const auto& kv : ctx.shared_trigger_map) {
      SharedTriggerInfo info;
      info.id = kv.first;
      info.acc_indices = kv.second;
      double q_val = 0.0;
      bool q_set = false;
      for (int idx : info.acc_indices) {
        if (idx < 0 || idx >= static_cast<int>(params.acc_params.size())) continue;
        q_val = params.acc_params[static_cast<std::size_t>(idx)].shared_q;
        q_set = true;
        break;
      }
      if (!q_set && !info.acc_indices.empty()) {
        int idx = info.acc_indices.front();
        if (idx >= 0 && idx < static_cast<int>(ctx.accumulators.size())) {
          q_val = ctx.accumulators[static_cast<std::size_t>(idx)].q;
          q_set = true;
        }
      }
      if (!q_set) continue;
      q_val = clamp_probability(q_val);
      info.q = q_val;
      out.push_back(info);
    }
  } else {
    // Fallback: derive shared triggers directly from the trial parameters.
    std::unordered_map<std::string, SharedTriggerInfo> derived;
    for (std::size_t i = 0; i < params.acc_params.size(); ++i) {
      const auto& acc = params.acc_params[i];
      if (acc.shared_trigger_id.empty()) continue;
      SharedTriggerInfo& info = derived[acc.shared_trigger_id];
      info.id = acc.shared_trigger_id;
      info.acc_indices.push_back(static_cast<int>(i));
      info.q = acc.q;
    }
    for (auto& kv : derived) {
      out.push_back(std::move(kv.second));
    }
  }
  return out;
}

// Apply a trigger success/fail mask to a base parameter set.
inline TrialParamSet apply_trigger_state(const TrialParamSet& base,
                                         const std::vector<SharedTriggerInfo>& triggers,
                                         std::uint64_t mask) {
  TrialParamSet modified = base;
  for (std::size_t i = 0; i < triggers.size(); ++i) {
    bool fail = (mask >> i) & 1ULL;
    double succeed_q = 0.0;  // when gate succeeds, drop the shared q to avoid double-counting
    double fail_q = 1.0;     // when gate fails, accumulator never finishes
    for (int acc_idx : triggers[i].acc_indices) {
      if (acc_idx < 0 || acc_idx >= static_cast<int>(modified.acc_params.size())) continue;
      TrialAccumulatorParams& p = modified.acc_params[static_cast<std::size_t>(acc_idx)];
      p.q = fail ? fail_q : succeed_q;
      p.has_override = true;
    }
  }
  return modified;
}

struct EvalCache;
struct ComponentCacheEntry;
const std::unordered_set<int> kEmptyForcedSet{};

double component_keep_weight(const uuber::NativeContext& ctx,
                             const std::string& component,
                             const std::string& outcome_label);

double native_outcome_probability_impl(SEXP ctxSEXP,
                                       int node_id,
                                       double upper,
                                       Rcpp::Nullable<Rcpp::String> component,
                                       SEXP forced_complete,
                                       SEXP forced_survive,
                                       const Rcpp::IntegerVector& competitor_ids,
                                       double rel_tol,
                                       double abs_tol,
                                       int max_depth,
                                       const TrialParamSet* trial_params,
                                       const std::string& trial_type_key,
                                       bool include_na_donors,
                                       const std::string& outcome_label_context);

std::string normalize_label_string(const std::string& label);

double node_density_with_competitors_internal(
  const uuber::NativeContext& ctx,
  int node_id,
  double t,
  const std::string& component_label,
  const std::unordered_set<int>& forced_complete,
  const std::unordered_set<int>& forced_survive,
  const std::vector<int>& competitor_ids,
  const TrialParamSet* trial_params,
  const std::string& trial_type_key,
  EvalCache* eval_cache_override,
  bool include_na_donors = false,
  const std::string& outcome_label_context = std::string());

double evaluate_outcome_density_label(const uuber::NativeContext& ctx,
                                      const std::string& donor_label,
                                      double t,
                                      const std::string& component,
                                      const TrialParamSet* trial_params,
                                      const std::string& trial_type_key,
                                      EvalCache* eval_cache_override,
                                      bool include_na_donors,
                                      const std::string& outcome_label_context = std::string()) {
  auto meta_it = ctx.outcome_info.find(donor_label);
  if (meta_it == ctx.outcome_info.end()) return 0.0;
  const uuber::OutcomeContextInfo& info = meta_it->second;
  if (info.node_id < 0) return 0.0;
  return node_density_with_competitors_internal(
    ctx,
    info.node_id,
    t,
    component,
    kEmptyForcedSet,
    kEmptyForcedSet,
    info.competitor_ids,
    trial_params,
    trial_type_key,
    eval_cache_override,
    include_na_donors,
    outcome_label_context);
}

double accumulate_component_guess_density(const uuber::NativeContext& ctx,
                                          const std::string& target_label,
                                          double t,
                                          const std::string& component,
                                          const TrialParamSet* trial_params,
                                          const std::string& trial_type_key,
                                          EvalCache* eval_cache_override) {
  auto comp_it = ctx.component_info.find(component);
  if (comp_it == ctx.component_info.end()) return 0.0;
  const uuber::ComponentGuessPolicy& guess = comp_it->second.guess;
  if (!guess.valid) return 0.0;
  std::string guess_target = normalize_label_string(guess.target);
  std::string target_key = normalize_label_string(target_label);
  if (!guess_target.empty() && guess_target != target_key) return 0.0;
  double total = 0.0;
  for (const auto& entry : guess.keep_weights) {
    double keep_prob = entry.second;
    double release = 1.0 - keep_prob;
    if (release <= 0.0) continue;
      total += release * evaluate_outcome_density_label(
        ctx,
        entry.first,
        t,
        component,
        trial_params,
        trial_type_key,
        eval_cache_override,
        false,
        entry.first);
  }
  return total;
}

double accumulate_outcome_alias_density(const uuber::NativeContext& ctx,
                                        const std::string& target_label,
                                        double t,
                                        const std::string& component,
                                        const TrialParamSet* trial_params,
                                        const std::string& trial_type_key,
                                        EvalCache* eval_cache_override,
                                        bool include_na_donors) {
  double total = 0.0;
  auto alias_it = ctx.alias_sources.find(target_label);
  if (alias_it != ctx.alias_sources.end()) {
    for (const std::string& source_label : alias_it->second) {
      total += evaluate_outcome_density_label(
        ctx,
        source_label,
        t,
        component,
        trial_params,
        trial_type_key,
        eval_cache_override,
        include_na_donors,
        source_label);
    }
  }
  auto guess_it = ctx.guess_target_map.find(target_label);
  if (guess_it != ctx.guess_target_map.end()) {
    for (const auto& donor : guess_it->second) {
      if (donor.rt_policy == "na" && !include_na_donors) continue;
      total += donor.weight * evaluate_outcome_density_label(
        ctx,
        donor.label,
        t,
        component,
        trial_params,
        trial_type_key,
        eval_cache_override,
        include_na_donors,
        donor.label);
    }
  }
  return total;
}


inline const TrialAccumulatorParams* get_trial_param_entry(const TrialParamSet* trial_params,
                                                          int acc_index) {
  if (!trial_params) return nullptr;
  if (acc_index < 0 || acc_index >= static_cast<int>(trial_params->acc_params.size())) {
    return nullptr;
  }
  const TrialAccumulatorParams& entry = trial_params->acc_params[acc_index];
  if (!entry.has_override) return nullptr;
  return &entry;
}

bool component_active(const uuber::NativeAccumulator& acc,
                      const std::string& component,
                      const TrialAccumulatorParams* override = nullptr) {
  if (component.empty() || component == "__default__") return true;
  const std::vector<std::string>* comps = nullptr;
  if (override && override->has_components) {
    comps = &override->components;
  } else if (!acc.components.empty()) {
    comps = &acc.components;
  }
  if (!comps || comps->empty()) return true;
  return std::find(comps->begin(), comps->end(), component) != comps->end();
}

NodeEvalResult eval_event_label(const uuber::NativeContext& ctx,
                                const std::string& label,
                                double t,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive,
                                const TrialParamSet* trial_params = nullptr,
                                const std::string& trial_type_key = std::string(),
                                EvalCache* eval_cache_override = nullptr,
                                bool include_na_donors = false,
                                const std::string& outcome_label_context = std::string()) {
  if (label.empty()) {
    return make_node_result(0.0, 1.0, 0.0);
  }
  if (label == "__DEADLINE__") {
    double survival = 0.0;
    double cdf = 1.0;
    if (std::isfinite(t)) {
      if (t < 0.0) {
        survival = 1.0;
        cdf = 0.0;
      } else {
        survival = 0.0;
        cdf = 1.0;
      }
    }
    return make_node_result(0.0, survival, cdf);
  }
  double donor_density = 0.0;
  donor_density += accumulate_component_guess_density(ctx,
                                                      label,
                                                      t,
                                                      component,
                                                      trial_params,
                                                      trial_type_key,
                                                      eval_cache_override);
  const std::string& target_label = outcome_label_context.empty() ? label : outcome_label_context;
  donor_density += accumulate_outcome_alias_density(ctx,
                                                    target_label,
                                                    t,
                                                    component,
                                                    trial_params,
                                                    trial_type_key,
                                                    eval_cache_override,
                                                    include_na_donors);
  if (label == "__GUESS__") {
    NodeEvalResult guess_res = make_node_result(donor_density, 0.0, 1.0);
    return guess_res;
  }
  int label_idx = label_id(ctx, label);
  if (label_idx != NA_INTEGER) {
    if (forced_complete.count(label_idx)) {
      return make_node_result(0.0, 0.0, 1.0);
    }
    if (forced_survive.count(label_idx)) {
      return make_node_result(0.0, 1.0, 0.0);
    }
  }

  auto acc_it = ctx.accumulator_index.find(label);
  if (acc_it != ctx.accumulator_index.end()) {
    int acc_idx = acc_it->second;
    const uuber::NativeAccumulator& acc = ctx.accumulators[acc_idx];
    const TrialAccumulatorParams* override = get_trial_param_entry(trial_params, acc_idx);
    if (!component_active(acc, component, override)) {
      return make_node_result(0.0, 1.0, 0.0);
    }
    double onset = override ? override->onset : acc.onset;
    double q = override ? override->q : acc.q;
    const AccDistParams& cfg = override ? override->dist_cfg : acc.dist_cfg;
    double density = acc_density_from_cfg(t, onset, q, cfg);
    double survival = acc_survival_from_cfg(t, onset, q, cfg);
    double cdf = acc_cdf_success_from_cfg(t, onset, cfg);
    NodeEvalResult res = make_node_result(density, survival, cdf);
    if (donor_density > 0.0) res.density += donor_density;
    return res;
  }

  auto pool_it = ctx.pool_index.find(label);
  if (pool_it != ctx.pool_index.end()) {
    const uuber::NativePool& pool = ctx.pools[pool_it->second];
    if (pool.members.empty()) {
      NodeEvalResult res = make_node_result(0.0, 1.0, 0.0);
      if (donor_density > 0.0) {
        res.density += donor_density;
      }
      return res;
    }
    std::vector<std::string> active_members;
    active_members.reserve(pool.members.size());
    for (const auto& member : pool.members) {
      auto member_acc = ctx.accumulator_index.find(member);
      if (member_acc != ctx.accumulator_index.end()) {
        const uuber::NativeAccumulator& acc = ctx.accumulators[member_acc->second];
        const TrialAccumulatorParams* override = get_trial_param_entry(trial_params, member_acc->second);
        if (!component_active(acc, component, override)) continue;
      }
      active_members.push_back(member);
    }
    if (active_members.empty()) {
      NodeEvalResult res = make_node_result(0.0, 1.0, 0.0);
      if (donor_density > 0.0) res.density += donor_density;
      return res;
    }
    Rcpp::NumericVector density_vec(active_members.size());
    Rcpp::NumericVector survival_vec(active_members.size());
    for (std::size_t i = 0; i < active_members.size(); ++i) {
      NodeEvalResult child = eval_event_label(
        ctx,
        active_members[i],
        t,
        component,
        forced_complete,
        forced_survive,
        trial_params,
        std::string(),
        nullptr
      );
      density_vec[i] = child.density;
      survival_vec[i] = child.survival;
    }
    double density = pool_density_fast_cpp(density_vec, survival_vec, pool.k);
    double survival = pool_survival_fast_cpp(survival_vec, pool.k);
    double cdf = clamp_probability(1.0 - survival);
    if (!std::isfinite(density) || density < 0.0) density = 0.0;
    if (!std::isfinite(survival)) survival = 0.0;
    NodeEvalResult res = make_node_result(density, survival, cdf);
    if (donor_density > 0.0) res.density += donor_density;
    return res;
  }

  return make_node_result(0.0, 1.0, 0.0);
}

struct ForcedKey {
  int complete_id{0};
  int survive_id{0};
};

struct TimeKey {
  double value{0.0};
};

struct EvalCache {};

std::string component_cache_key(const std::string& component,
                                const std::string& trial_key);

struct ComponentCacheEntry {
  std::string trial_type_key;
};

ComponentCacheEntry make_component_cache_entry(const std::string& component,
                                               const std::string& trial_key) {
  ComponentCacheEntry entry;
  entry.trial_type_key = component_cache_key(component, trial_key);
  return entry;
}

ComponentCacheEntry default_component_cache_entry(const std::string& component) {
  return make_component_cache_entry(component, std::string());
}

double accumulate_plan_guess_density(const uuber::NativeContext& ctx,
                                     const Rcpp::List& donors,
                                     double t,
                                     const std::string& component,
                                     const std::string& trial_type_key,
                                     EvalCache& eval_cache,
                                     const TrialParamSet* trial_params,
                                     bool include_na_donors) {
  double total = 0.0;
  for (R_xlen_t i = 0; i < donors.size(); ++i) {
    Rcpp::RObject donor_obj = donors[i];
    if (donor_obj.isNULL()) continue;
    Rcpp::List donor(donor_obj);
    std::string rt_policy = donor.containsElementNamed("rt_policy")
      ? Rcpp::as<std::string>(donor["rt_policy"])
      : std::string("keep");
    bool donor_na = rt_policy == "na";
    if (donor_na != include_na_donors) continue;
    double release = donor.containsElementNamed("release")
      ? Rcpp::as<double>(donor["release"])
      : 0.0;
    if (!std::isfinite(release)) release = 0.0;
    if (release <= 0.0) continue;
    int donor_node = donor.containsElementNamed("node_id")
      ? Rcpp::as<int>(donor["node_id"])
      : NA_INTEGER;
    if (donor_node == NA_INTEGER) continue;
    Rcpp::IntegerVector comp_ids = donor.containsElementNamed("competitor_ids")
      ? Rcpp::IntegerVector(donor["competitor_ids"])
      : Rcpp::IntegerVector();
    std::vector<int> comp_vec = integer_vector_to_std(comp_ids, false);
    std::string donor_label;
    if (donor.containsElementNamed("source_label") && !Rf_isNull(donor["source_label"])) {
      donor_label = Rcpp::as<std::string>(donor["source_label"]);
    }
    double density = node_density_with_competitors_internal(
      ctx,
      donor_node,
      t,
      component,
      kEmptyForcedSet,
      kEmptyForcedSet,
      comp_vec,
      trial_params,
      trial_type_key,
      &eval_cache,
      include_na_donors,
      donor_label
    );
    if (!std::isfinite(density) || density <= 0.0) continue;
    total += release * density;
  }
  return total;
}

std::vector<ComponentCacheEntry> build_component_cache_entries(
  const std::vector<std::string>& components) {
  std::vector<ComponentCacheEntry> entries;
  entries.reserve(components.size());
  for (const auto& comp : components) {
    entries.push_back(default_component_cache_entry(comp));
  }
  return entries;
}

std::string component_cache_key(const std::string& component,
                                const std::string& trial_key = std::string()) {
  if (!trial_key.empty()) return trial_key;
  if (component.empty()) return "__default__";
  return component;
}

std::string normalize_label_string(const std::string& label) {
  if (label.empty()) return std::string();
  return label;
}

double component_keep_weight(const uuber::NativeContext& ctx,
                             const std::string& component_label,
                             const std::string& outcome_label) {
  auto comp_it = ctx.component_info.find(component_label);
  if (comp_it == ctx.component_info.end()) return 1.0;
  const uuber::ComponentGuessPolicy& guess = comp_it->second.guess;
  if (!guess.valid) return 1.0;
  auto keep_it = guess.keep_weights.find(outcome_label);
  if (keep_it != guess.keep_weights.end()) return keep_it->second;
  if (outcome_label == "NA") {
    auto na_it = guess.keep_weights.find("<NA>");
    if (na_it != guess.keep_weights.end()) return na_it->second;
  }
  return 1.0;
}

inline double extract_trial_id(const Rcpp::DataFrame* trial_rows) {
  if (!trial_rows) return std::numeric_limits<double>::quiet_NaN();
  if (!trial_rows->containsElementNamed("trial")) return std::numeric_limits<double>::quiet_NaN();
  Rcpp::NumericVector trial_col((*trial_rows)["trial"]);
  if (trial_col.size() == 0) return std::numeric_limits<double>::quiet_NaN();
  double val = trial_col[0];
  if (Rcpp::NumericVector::is_na(val)) return std::numeric_limits<double>::quiet_NaN();
  return static_cast<double>(val);
}


inline std::uint64_t double_to_bits(double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline std::vector<std::uint64_t> build_forced_bits(const uuber::NativeContext& ctx,
                                                    const std::unordered_set<int>& forced) {
  std::size_t node_count = ctx.nodes.size();
  std::size_t word_count = (node_count + 63) / 64;
  std::vector<std::uint64_t> bits(word_count, 0);
  if (forced.empty() || word_count == 0) return bits;
  for (int id : forced) {
    auto it = ctx.node_index.find(id);
    if (it == ctx.node_index.end()) continue;
    std::size_t idx = static_cast<std::size_t>(it->second);
    std::size_t word = idx >> 6;
    std::size_t offset = idx & 63;
    if (word < bits.size()) {
      bits[word] |= (std::uint64_t{1} << offset);
    }
  }
  return bits;
}

inline std::uint64_t hash_combine_u64(std::uint64_t seed, std::uint64_t value) {
  seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
  return seed;
}

inline std::uint64_t hash_forced_bits(const std::vector<std::uint64_t>& bits) {
  std::uint64_t h = 0x9e3779b185ebca87ULL;
  for (std::uint64_t v : bits) {
    h = hash_combine_u64(h, v);
  }
  return h;
}

struct NodeEvalState {
  NodeEvalState(const uuber::NativeContext& ctx_,
                EvalCache& cache_,
                double time,
                const std::string& component_label,
                const std::unordered_set<int>& forced_complete_ref,
                const std::unordered_set<int>& forced_survive_ref,
                const TrialParamSet* params_ptr = nullptr,
                const std::string& trial_key = std::string(),
                bool include_na = false,
                const std::string& outcome_label_ctx = std::string())
    : ctx(ctx_),
      cache(cache_),
      t(time),
      component(component_label),
      forced_complete(forced_complete_ref),
      forced_survive(forced_survive_ref),
      trial_params(params_ptr),
      trial_type_key(component_cache_key(component_label, trial_key)),
      include_na_donors(include_na),
      outcome_label(outcome_label_ctx),
      forced_complete_bits(build_forced_bits(ctx_, forced_complete_ref)),
      forced_survive_bits(build_forced_bits(ctx_, forced_survive_ref)),
      forced_complete_hash(hash_forced_bits(forced_complete_bits)),
      forced_survive_hash(hash_forced_bits(forced_survive_bits)),
      time_bits(double_to_bits(time)) {}

  const uuber::NativeContext& ctx;
  EvalCache& cache;
  double t;
  std::string component;
  const std::unordered_set<int>& forced_complete;
  const std::unordered_set<int>& forced_survive;
  const TrialParamSet* trial_params;
  std::string trial_type_key;
  bool include_na_donors;
  std::string outcome_label;
  std::vector<std::uint64_t> forced_complete_bits;
  std::vector<std::uint64_t> forced_survive_bits;
  std::uint64_t forced_complete_hash;
  std::uint64_t forced_survive_hash;
  std::uint64_t time_bits;
};

inline const TrialAccumulatorParams* get_state_params(const NodeEvalState& state,
                                                          int acc_index) {
  return get_trial_param_entry(state.trial_params, acc_index);
}

std::string member_ids_signature(const std::vector<int>& ids) {
  if (ids.empty()) return ".";
  std::ostringstream oss;
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (i > 0) oss << ",";
    oss << ids[i];
  }
  return oss.str();
}

std::string shared_ids_signature(const std::vector<std::string>& shared_ids) {
  if (shared_ids.empty()) return ".";
  std::ostringstream oss;
  for (std::size_t i = 0; i < shared_ids.size(); ++i) {
    if (i > 0) oss << ",";
    if (shared_ids[i].empty()) {
      oss << ".";
    } else {
      oss << shared_ids[i].size() << ":" << shared_ids[i];
    }
  }
  return oss.str();
}

std::vector<int> build_shared_index(const std::vector<std::string>& shared_ids) {
  std::vector<int> shared_index(shared_ids.size(), -1);
  std::unordered_map<std::string, int> group_map;
  int next_group = 1;
  for (std::size_t i = 0; i < shared_ids.size(); ++i) {
    const std::string& sid = shared_ids[i];
    if (sid.empty()) continue;
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

std::string pool_template_cache_key(const std::string& pool_label,
                                    const std::string& component,
                                    int k,
                                    const std::vector<int>& member_ids,
                                    const std::string& shared_signature) {
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
  const uuber::NativeContext& ctx;
  const uuber::NativeNode& node;
  EvalCache& cache;
  std::string component;
  std::string trial_type_key;
  std::unordered_set<int> forced_complete;
  std::unordered_set<int> forced_survive;
  const TrialParamSet* trial_params;
};

const uuber::NativeNode& fetch_node(const uuber::NativeContext& ctx, int node_id);
NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                     EvalCache& cache,
                                     int node_id,
                                     double time,
                                     const std::string& component,
                                     const std::unordered_set<int>& forced_complete,
                                     const std::unordered_set<int>& forced_survive,
                                     const TrialParamSet* trial_params,
                                     const std::string& trial_key);
inline NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                            EvalCache& cache,
                                            int node_id,
                                            double time,
                                            const std::string& component,
                                            const std::unordered_set<int>& forced_complete,
                                            const std::unordered_set<int>& forced_survive) {
  return eval_node_with_forced(ctx,
                               cache,
                               node_id,
                               time,
                               component,
                               forced_complete,
                               forced_survive,
                               nullptr,
                               std::string());
}
GuardEvalInput make_guard_input(const uuber::NativeContext& ctx,
                                const uuber::NativeNode& node,
                                EvalCache& cache,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive,
                                const std::string& trial_key,
                                const TrialParamSet* trial_params);
std::vector<int> gather_blocker_sources(const uuber::NativeContext& ctx,
                                        const uuber::NativeNode& guard_node);

std::vector<ScenarioRecord> compute_node_scenarios(int node_id, NodeEvalState& state);

int first_missing_id(const std::vector<int>& forced_complete,
                     const std::vector<int>& required) {
  if (required.empty()) return NA_INTEGER;
  for (int id : required) {
    if (!vector_contains(forced_complete, id)) {
      return id;
    }
  }
  return NA_INTEGER;
}

std::vector<ScenarioRecord> compute_pool_event_scenarios(const uuber::NativeNode& node,
                                                         NodeEvalState& state) {
  const std::string& pool_label = node.source;
  if (pool_label.empty()) return {};
  auto pool_it = state.ctx.pool_index.find(pool_label);
  if (pool_it == state.ctx.pool_index.end()) return {};
  const uuber::NativePool& pool = state.ctx.pools[pool_it->second];
  if (pool.members.empty()) return {};
  int k = pool.k;
  if (k < 1) return {};

  std::vector<std::string> active_members;
  active_members.reserve(pool.members.size());
  for (const auto& member : pool.members) {
    auto acc_it = state.ctx.accumulator_index.find(member);
    if (acc_it != state.ctx.accumulator_index.end()) {
      const uuber::NativeAccumulator& acc = state.ctx.accumulators[acc_it->second];
      const TrialAccumulatorParams* override = get_state_params(state, acc_it->second);
      if (!component_active(acc, state.component, override)) continue;
    }
    active_members.push_back(member);
  }
  std::size_t n = active_members.size();
  if (n == 0 || k > static_cast<int>(n)) return {};

  Rcpp::NumericVector dens_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector cdf_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector surv_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector cdf_success_vec(static_cast<R_xlen_t>(n));
  Rcpp::NumericVector surv_success_vec(static_cast<R_xlen_t>(n));
  std::vector<int> member_ids_std(n, NA_INTEGER);
  std::vector<std::string> shared_ids(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::string& member_label = active_members[i];
    member_ids_std[i] = label_id(state.ctx, member_label);
    auto acc_it = state.ctx.accumulator_index.find(member_label);
    if (acc_it != state.ctx.accumulator_index.end()) {
      const uuber::NativeAccumulator& acc = state.ctx.accumulators[acc_it->second];
      const TrialAccumulatorParams* override = get_state_params(state, acc_it->second);
      const std::string& shared_id = override ? override->shared_trigger_id : acc.shared_trigger_id;
      if (!shared_id.empty()) {
        shared_ids[i] = shared_id;
      }
    }
  }


  for (std::size_t i = 0; i < n; ++i) {
    const std::string& member_label = active_members[i];
    NodeEvalResult eval = eval_event_label(
      state.ctx,
      member_label,
      state.t,
      state.component,
      state.forced_complete,
      state.forced_survive,
      state.trial_params,
      state.trial_type_key,
      &state.cache
    );
    double density = eval.density;
    double cdf = clamp_probability(eval.cdf);
    double survival = clamp_probability(eval.survival);
    double cdf_success = cdf;
    double surv_success = survival;
    auto acc_it = state.ctx.accumulator_index.find(member_label);
    if (acc_it != state.ctx.accumulator_index.end()) {
      const uuber::NativeAccumulator& acc = state.ctx.accumulators[acc_it->second];
      const TrialAccumulatorParams* override = get_state_params(state, acc_it->second);
      const std::string& shared_id = override ? override->shared_trigger_id : acc.shared_trigger_id;
      if (!shared_id.empty()) {
        double onset = override ? override->onset : acc.onset;
        const AccDistParams& cfg = override ? override->dist_cfg : acc.dist_cfg;
        cdf_success = clamp_probability(acc_cdf_success_from_cfg(
          state.t,
          onset,
          cfg
        ));
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

  int pool_label_id = label_id(state.ctx, pool.id);
  std::string shared_signature = shared_ids_signature(shared_ids);
  std::string template_key = pool_template_cache_key(
    pool_label,
    component_cache_key(state.component),
    k,
    member_ids_std,
    shared_signature
  );
  uuber::PoolTemplateCacheEntry* template_entry = nullptr;
  auto tmpl_it = state.ctx.pool_template_cache.find(template_key);
  if (tmpl_it != state.ctx.pool_template_cache.end()) {
    template_entry = &tmpl_it->second;
  } else {
    uuber::PoolTemplateCacheEntry cache_entry = build_pool_template_cache(
      static_cast<int>(n),
      member_ids_std,
      pool_label_id,
      k
    );
    cache_entry.shared_index = build_shared_index(shared_ids);
    template_entry = &state.ctx.pool_template_cache.emplace(
      template_key,
      std::move(cache_entry)
    ).first->second;
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

  std::vector<int> base_forced_complete = set_to_sorted_vector(state.forced_complete);
  std::vector<int> base_forced_survive = set_to_sorted_vector(state.forced_survive);

  Rcpp::List combine_res = pool_density_combine_native(
    dens_vec,
    cdf_vec,
    surv_vec,
    cdf_success_vec,
    surv_success_vec,
    template_entry->shared_index,
    template_entry->templates,
    base_forced_complete,
    base_forced_survive
  );
  Rcpp::List scenario_list = combine_res["scenarios"];
  return rcpp_scenarios_to_records(scenario_list);
}

std::vector<ScenarioRecord> compute_event_scenarios(const uuber::NativeNode& node,
                                                    NodeEvalState& state) {
  std::vector<ScenarioRecord> out;
  if (node.source.empty()) {
    return out;
  }
  auto pool_lookup = state.ctx.pool_index.find(node.source);
  if (pool_lookup != state.ctx.pool_index.end()) {
    std::vector<ScenarioRecord> pool_records = compute_pool_event_scenarios(node, state);
    if (!pool_records.empty()) {
      return pool_records;
    }
  }
  NodeEvalResult eval = eval_event_label(
    state.ctx,
    node.source,
    state.t,
    state.component,
    state.forced_complete,
    state.forced_survive,
    state.trial_params,
    state.trial_type_key,
    &state.cache
  );
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

std::vector<ScenarioRecord> compute_and_scenarios(const uuber::NativeNode& node,
                                                  NodeEvalState& state) {
  std::vector<ScenarioRecord> out;
  const std::vector<int>& child_ids = node.args;
  if (child_ids.empty()) return out;
  std::vector<const uuber::NativeNode*> child_nodes;
  std::vector<std::vector<int>> child_sources;
  child_nodes.reserve(child_ids.size());
  child_sources.reserve(child_ids.size());
  for (int child_id : child_ids) {
    const uuber::NativeNode& child = fetch_node(state.ctx, child_id);
    child_nodes.push_back(&child);
    child_sources.push_back(ensure_source_ids(state.ctx, child));
  }
  for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
    std::vector<ScenarioRecord> child_scen = compute_node_scenarios(child_ids[idx], state);
    if (child_scen.empty()) continue;
    for (const auto& sc : child_scen) {
      if (!std::isfinite(sc.weight) || sc.weight <= 0.0) continue;
      double weight = sc.weight;
      std::vector<int> forced_complete = sc.forced_complete;
      std::vector<int> forced_survive = sc.forced_survive;
      bool ok = true;
      for (std::size_t other_idx = 0; other_idx < child_ids.size(); ++other_idx) {
        if (other_idx == idx) continue;
        std::unordered_set<int> fc_set = make_forced_set(forced_complete);
        std::unordered_set<int> fs_set = make_forced_set(forced_survive);
        NodeEvalResult sibling = eval_node_with_forced(state.ctx,
                                                       state.cache,
                                                       child_ids[other_idx],
                                                       state.t,
                                                       state.component,
                                                       fc_set,
                                                       fs_set);
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
        forced_complete = union_vectors(forced_complete, child_sources[other_idx]);
      }
      if (!ok) continue;
      ScenarioRecord merged = make_scenario_record(weight, forced_complete, forced_survive);
      if (merged.weight > 0.0) {
        out.push_back(std::move(merged));
      }
    }
  }
  return out;
}

std::vector<ScenarioRecord> compute_or_scenarios(const uuber::NativeNode& node,
                                                 NodeEvalState& state) {
  std::vector<ScenarioRecord> out;
  const std::vector<int>& child_ids = node.args;
  if (child_ids.empty()) return out;
  std::vector<const uuber::NativeNode*> child_nodes;
  std::vector<std::vector<int>> child_sources;
  child_nodes.reserve(child_ids.size());
  child_sources.reserve(child_ids.size());
  for (int child_id : child_ids) {
    const uuber::NativeNode& child = fetch_node(state.ctx, child_id);
    child_nodes.push_back(&child);
    child_sources.push_back(ensure_source_ids(state.ctx, child));
  }
  for (std::size_t idx = 0; idx < child_ids.size(); ++idx) {
    std::vector<ScenarioRecord> child_scen = compute_node_scenarios(child_ids[idx], state);
    if (child_scen.empty()) continue;
    for (const auto& sc : child_scen) {
      if (!std::isfinite(sc.weight) || sc.weight <= 0.0) continue;
      double weight = sc.weight;
      std::vector<int> forced_complete = sc.forced_complete;
      std::vector<int> forced_survive = sc.forced_survive;
      bool valid = true;
      for (std::size_t other_idx = 0; other_idx < child_ids.size(); ++other_idx) {
        if (other_idx == idx) continue;
        const std::vector<int>& required = child_sources[other_idx];
        if (required.empty()) continue;
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
      if (!valid) continue;
      ScenarioRecord merged = make_scenario_record(weight, forced_complete, forced_survive);
      if (merged.weight > 0.0) {
        out.push_back(std::move(merged));
      }
    }
  }
  return out;
}

struct GuardMemoKey {
  int node_id;
  std::uint64_t t_bits;
  bool operator==(const GuardMemoKey& other) const noexcept {
    return node_id == other.node_id && t_bits == other.t_bits;
  }
};

struct GuardMemoKeyHash {
  std::size_t operator()(const GuardMemoKey& key) const noexcept {
    std::size_t h = static_cast<std::size_t>(key.node_id);
    h ^= static_cast<std::size_t>(key.t_bits + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
    return h;
  }
};

using GuardMemo = std::unordered_map<GuardMemoKey, NodeEvalResult, GuardMemoKeyHash>;

inline NodeEvalResult eval_node_guard_cached(const GuardEvalInput& input,
                                             int node_id,
                                             double t,
                                             GuardMemo* memo) {
  if (!memo) {
    return eval_node_with_forced(
      input.ctx,
      input.cache,
      node_id,
      t,
      input.component,
      input.forced_complete,
      input.forced_survive,
      input.trial_params,
      input.trial_type_key
    );
  }
  GuardMemoKey key{node_id, double_to_bits(t)};
  auto it = memo->find(key);
  if (it != memo->end()) {
    return it->second;
  }
  NodeEvalResult res = eval_node_with_forced(
    input.ctx,
    input.cache,
    node_id,
    t,
    input.component,
    input.forced_complete,
    input.forced_survive,
    input.trial_params,
    input.trial_type_key
  );
  memo->emplace(key, res);
  return res;
}

double guard_effective_survival_internal(const GuardEvalInput& input,
                                         double t,
                                         const IntegrationSettings& settings,
                                         GuardMemo* memo = nullptr);

double guard_density_internal(const GuardEvalInput& input,
                              double t,
                              const IntegrationSettings& settings,
                              GuardMemo* memo = nullptr);

double guard_cdf_internal(const GuardEvalInput& input,
                          double t,
                          const IntegrationSettings& settings);

std::vector<ScenarioRecord> compute_guard_scenarios(const uuber::NativeNode& node,
                                                    NodeEvalState& state) {
  std::vector<ScenarioRecord> out;
  if (node.reference_id < 0) {
    return out;
  }
  std::vector<ScenarioRecord> ref_scen = compute_node_scenarios(node.reference_id, state);
  if (ref_scen.empty()) return out;
  std::vector<int> blocker_sources = gather_blocker_sources(state.ctx, node);
  IntegrationSettings settings;
  out.reserve(ref_scen.size());
  for (const auto& sc : ref_scen) {
    if (!std::isfinite(sc.weight) || sc.weight <= 0.0) continue;
    std::unordered_set<int> fc_set = make_forced_set(sc.forced_complete);
    std::unordered_set<int> fs_set = make_forced_set(sc.forced_survive);
    GuardEvalInput guard_input = make_guard_input(state.ctx,
                                                  node,
                                                  state.cache,
                                                  state.component,
                                                  fc_set,
                                                  fs_set,
                                                  state.trial_type_key,
                                                  state.trial_params);
    double eff = guard_effective_survival_internal(guard_input, state.t, settings);
    if (!std::isfinite(eff) || eff <= 0.0) continue;
    double weight = sc.weight * eff;
    if (!std::isfinite(weight) || weight <= 0.0) continue;
    ScenarioRecord rec;
    rec.weight = weight;
    rec.forced_complete = sc.forced_complete;
    rec.forced_survive = union_vectors(sc.forced_survive, blocker_sources);
    out.push_back(std::move(rec));
  }
  return out;
}

std::vector<ScenarioRecord> compute_node_scenarios(int node_id, NodeEvalState& state) {
  const uuber::NativeNode& node = fetch_node(state.ctx, node_id);
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

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state);

const uuber::NativeNode& fetch_node(const uuber::NativeContext& ctx, int node_id) {
  return ctx.nodes.at(ctx.node_index.at(node_id));
}

NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                     EvalCache& cache,
                                     int node_id,
                                     double time,
                                     const std::string& component,
                                     const std::unordered_set<int>& forced_complete,
                                     const std::unordered_set<int>& forced_survive,
                                     const TrialParamSet* trial_params,
                                     const std::string& trial_key) {
  NodeEvalState local(ctx,
                      cache,
                      time,
                      component,
                      forced_complete,
                      forced_survive,
                      trial_params,
                      trial_key);
  return eval_node_recursive(node_id, local);
}

NodeEvalResult eval_and_node(const uuber::NativeNode& node, NodeEvalState& state) {
  const std::vector<int>& child_ids = node.args;
  if (child_ids.empty()) {
    return make_node_result(0.0, 1.0, 0.0);
  }
  std::size_t n = child_ids.size();
  std::vector<NodeEvalResult> children;
  children.reserve(n);
  for (int child_id : child_ids) {
    children.push_back(eval_node_recursive(child_id, state));
  }
  std::vector<double> child_cdf(n);
  for (std::size_t i = 0; i < n; ++i) {
    child_cdf[i] = clamp_probability(children[i].cdf);
  }
  std::vector<double> prefix(n + 1, 1.0);
  std::vector<double> suffix(n + 1, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    prefix[i + 1] = prefix[i] * child_cdf[i];
    if (!std::isfinite(prefix[i + 1])) prefix[i + 1] = 0.0;
  }
  for (std::size_t i = n; i-- > 0;) {
    suffix[i] = suffix[i + 1] * child_cdf[i];
    if (!std::isfinite(suffix[i])) suffix[i] = 0.0;
  }
  double total_cdf = clamp_probability(prefix[n]);
  double survival = clamp_probability(1.0 - total_cdf);
  double density = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double child_density = children[i].density;
    if (!std::isfinite(child_density) || child_density <= 0.0) continue;
    double others = prefix[i] * suffix[i + 1];
    if (!std::isfinite(others) || others <= 0.0) continue;
    density += child_density * clamp_probability(others);
  }
  if (!std::isfinite(density) || density < 0.0) density = 0.0;
  return make_node_result(density, survival, total_cdf);
}

NodeEvalResult eval_or_node(const uuber::NativeNode& node, NodeEvalState& state) {
  const std::vector<int>& child_ids = node.args;
  if (child_ids.empty()) {
    return make_node_result(0.0, 1.0, 0.0);
  }
  std::size_t n = child_ids.size();
  std::vector<NodeEvalResult> children;
  children.reserve(n);
  for (int child_id : child_ids) {
    children.push_back(eval_node_recursive(child_id, state));
  }
  std::vector<double> child_surv(n);
  for (std::size_t i = 0; i < n; ++i) {
    child_surv[i] = clamp_probability(children[i].survival);
  }
  std::vector<double> prefix(n + 1, 1.0);
  std::vector<double> suffix(n + 1, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    prefix[i + 1] = prefix[i] * child_surv[i];
    if (!std::isfinite(prefix[i + 1])) prefix[i + 1] = 0.0;
  }
  for (std::size_t i = n; i-- > 0;) {
    suffix[i] = suffix[i + 1] * child_surv[i];
    if (!std::isfinite(suffix[i])) suffix[i] = 0.0;
  }
  double total_surv = clamp_probability(prefix[n]);
  double cdf = clamp_probability(1.0 - total_surv);
  double density = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double child_density = children[i].density;
    if (!std::isfinite(child_density) || child_density <= 0.0) continue;
    double others = prefix[i] * suffix[i + 1];
    if (!std::isfinite(others) || others <= 0.0) continue;
    density += child_density * clamp_probability(others);
  }
  if (!std::isfinite(density) || density < 0.0) density = 0.0;
  return make_node_result(density, total_surv, cdf);
}

NodeEvalResult eval_not_node(const uuber::NativeNode& node, NodeEvalState& state) {
  int child_id = node.arg_id;
  if (child_id < 0) {
    return make_node_result(0.0, 1.0, 0.0);
  }
  NodeEvalResult child = eval_node_recursive(child_id, state);
  double child_cdf = clamp_probability(child.cdf);
  double self_cdf = clamp_probability(1.0 - child_cdf);
  double survival = clamp_probability(1.0 - self_cdf);
  return make_node_result(0.0, survival, self_cdf);
}

GuardEvalInput make_guard_input(const uuber::NativeContext& ctx,
                                const uuber::NativeNode& node,
                                EvalCache& cache,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive,
                                const std::string& trial_key = std::string(),
                                const TrialParamSet* trial_params = nullptr) {
  GuardEvalInput input{
    ctx,
    node,
    cache,
    component,
    component_cache_key(component, trial_key),
    {},
    {},
    trial_params
  };
  input.forced_complete = filter_forced_scope(forced_complete, node.source_ids);
  input.forced_survive = filter_forced_scope(forced_survive, node.source_ids);
  return input;
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state) {
  const uuber::NativeNode& node = fetch_node(state.ctx, node_id);
  NodeEvalResult result;
  if (node.kind == "event") {
    result = eval_event_label(
      state.ctx,
      node.source,
      state.t,
      state.component,
      state.forced_complete,
      state.forced_survive,
      state.trial_params,
      state.trial_type_key,
      &state.cache,
      state.include_na_donors,
      state.outcome_label
    );
  } else if (node.kind == "and") {
    result = eval_and_node(node, state);
  } else if (node.kind == "or") {
    result = eval_or_node(node, state);
  } else if (node.kind == "not") {
    result = eval_not_node(node, state);
  } else if (node.kind == "guard") {
    GuardEvalInput guard_input = make_guard_input(state.ctx,
                                                  node,
                                                  state.cache,
                                                  state.component,
                                                  state.forced_complete,
                                                  state.forced_survive,
                                                  state.trial_type_key,
                                                  state.trial_params);
    IntegrationSettings settings;
    double density = guard_density_internal(guard_input, state.t, settings);
    double cdf = guard_cdf_internal(guard_input, state.t, settings);
    double survival = clamp_probability(1.0 - cdf);
    result = make_node_result(density, survival, cdf);
  } else {
    result = make_node_result(0.0, 1.0, 0.0);
  }
  return result;
}

double guard_effective_survival_internal(const GuardEvalInput& input,
                                         double t,
                                         const IntegrationSettings& settings,
                                         GuardMemo* memo) {
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
  if (input.node.unless_ids.empty()) {
    NodeEvalResult block = eval_node_guard_cached(
      input,
      blocker_id,
      t,
      memo
    );
    return clamp_probability(block.survival);
  }
  auto integrand = [&](double u) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    NodeEvalResult block = eval_node_guard_cached(
      input,
      blocker_id,
      u,
      memo
    );
    double density = block.density;
    if (!std::isfinite(density) || density <= 0.0) return 0.0;
    double prod = 1.0;
    for (int unl_id : input.node.unless_ids) {
      if (unl_id < 0) continue;
      NodeEvalResult protector = eval_node_guard_cached(
        input,
        unl_id,
        u,
        memo
      );
      prod *= clamp_probability(protector.survival);
      if (!std::isfinite(prod) || prod <= 0.0) {
        prod = 0.0;
        break;
      }
    }
    if (prod <= 0.0) return 0.0;
    double val = density * prod;
    if (!std::isfinite(val) || val <= 0.0) return 0.0;
    return val;
  };
  double integral = uuber::integrate_fixed_gauss15(integrand, 0.0, t);
  double surv = 1.0 - integral;
  return clamp_probability(surv);
}

double guard_reference_density(const GuardEvalInput& input,
                               double t,
                               GuardMemo* memo = nullptr) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  int reference_id = input.node.reference_id;
  if (reference_id < 0) {
    return 0.0;
  }
  NodeEvalResult ref = eval_node_guard_cached(
    input,
    reference_id,
    t,
    memo
  );
  double dens_ref = ref.density;
  if (!std::isfinite(dens_ref) || dens_ref <= 0.0) return 0.0;
  return dens_ref;
}

double guard_density_internal(const GuardEvalInput& input,
                              double t,
                              const IntegrationSettings& settings,
                              GuardMemo* memo) {
  double dens_ref = guard_reference_density(input, t, memo);
  if (dens_ref <= 0.0) return 0.0;
  double s_eff = guard_effective_survival_internal(input, t, settings, memo);
  double val = dens_ref * s_eff;
  if (!std::isfinite(val) || val <= 0.0) return 0.0;
  return val;
}

double guard_cdf_internal(const GuardEvalInput& input,
                          double t,
                          const IntegrationSettings& settings) {
  if (!std::isfinite(t)) {
    return 1.0;
  }
  if (t <= 0.0) {
    return 0.0;
  }
  GuardMemo memo;
  auto integrand = [&](double u) -> double {
    return guard_density_internal(input, u, settings, &memo);
  };
  double val = uuber::integrate_fixed_gauss15(integrand, 0.0, t);
  if (!std::isfinite(val) || val < 0.0) val = 0.0;
  return clamp_probability(val);
}

std::vector<int> gather_blocker_sources(const uuber::NativeContext& ctx,
                                        const uuber::NativeNode& guard_node) {
  std::vector<int> sources;
  const uuber::NativeNode& blocker = fetch_node(ctx, guard_node.blocker_id);
  sources = blocker.source_ids;
  if (blocker.kind == "event" && !blocker.source.empty()) {
    sources.push_back(label_id(ctx, blocker.source));
  }
  sort_unique(sources);
  return sources;
}

struct CompetitorMeta {
  int node_id;
  const uuber::NativeNode* node;
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

using CompetitorCacheMap = std::unordered_map<std::string, CompetitorClusterCacheEntry>;

static std::unordered_map<const uuber::NativeContext*, CompetitorCacheMap> g_competitor_cache;

std::string competitor_cache_key(const std::vector<int>& ids) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (i > 0) oss << ",";
    oss << ids[i];
  }
  return oss.str();
}

// Forward declaration
std::vector<std::vector<int>> build_competitor_clusters(const std::vector<CompetitorMeta>& metas);

bool node_contains_guard(int node_id,
                         const uuber::NativeContext& ctx,
                         std::unordered_map<int, bool>& memo) {
  auto it = memo.find(node_id);
  if (it != memo.end()) return it->second;
  const uuber::NativeNode& node = fetch_node(ctx, node_id);
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

const CompetitorClusterCacheEntry& fetch_competitor_cluster_cache(
  const uuber::NativeContext& ctx,
  const std::vector<int>& competitor_ids) {
  CompetitorCacheMap& cache = g_competitor_cache[&ctx];
  const std::string key = competitor_cache_key(competitor_ids);
  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }

  CompetitorClusterCacheEntry entry;
  entry.metas.reserve(competitor_ids.size());
  std::unordered_map<int, bool> guard_memo;
  for (int node_id : competitor_ids) {
    const uuber::NativeNode& node = fetch_node(ctx, node_id);
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
    const auto& cluster = entry.clusters[i];
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
        const CompetitorMeta& meta = entry.metas[static_cast<std::size_t>(idx)];
        order.push_back({idx, meta.sources.size(), meta.scenario_sensitive});
      }
      std::sort(order.begin(), order.end(), [&](const OrderingMeta& a, const OrderingMeta& b) {
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
      for (const auto& rec : order) {
        guard_order.push_back(rec.index);
      }
      entry.guard_orders[i] = std::move(guard_order);
    }
  }

  auto inserted = cache.emplace(key, std::move(entry));
  return inserted.first->second;
}

bool share_sources(const std::vector<int>& a, const std::vector<int>& b) {
  if (a.empty() || b.empty()) return false;
  std::size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j]) return true;
    if (a[i] < b[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return false;
}

std::vector<std::vector<int>> build_competitor_clusters(const std::vector<CompetitorMeta>& metas) {
  std::vector<std::vector<int>> clusters;
  const std::size_t n = metas.size();
  if (n == 0) return clusters;
  std::vector<bool> visited(n, false);
  for (std::size_t i = 0; i < n; ++i) {
    if (visited[i]) continue;
    std::vector<int> cluster;
    std::queue<std::size_t> q;
    q.push(i);
    visited[i] = true;
    while (!q.empty()) {
      std::size_t idx = q.front();
      q.pop();
      cluster.push_back(static_cast<int>(idx));
      for (std::size_t j = 0; j < n; ++j) {
        if (visited[j]) continue;
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

uuber::PoolTemplateCacheEntry build_pool_template_cache(int n,
                                                       const std::vector<int>& member_ids,
                                                       int pool_idx,
                                                       int k) {
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
      if (j == idx) continue;
      others.push_back(j + 1);
    }
    if (need > static_cast<int>(others.size())) {
      continue;
    }
    std::vector<std::vector<int>> combos = generate_combinations(others, need);
    for (const auto& combo : combos) {
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
      cache.finisher_map[static_cast<std::size_t>(idx)].push_back(template_counter);
      ++template_counter;
    }
  }
  return cache;
}

double evaluate_survival_with_forced(int node_id,
                                     const std::unordered_set<int>& forced_survive,
                                     const std::string& component,
                                     double t,
                                     const uuber::NativeContext& ctx,
                                     const std::string& trial_key = std::string(),
                                     const TrialParamSet* trial_params = nullptr,
                                     EvalCache* eval_cache_override = nullptr) {
  static const std::unordered_set<int> empty_forced_complete;
  EvalCache local_cache;
  EvalCache& cache = eval_cache_override ? *eval_cache_override : local_cache;
  NodeEvalState state(ctx,
                      cache,
                      t,
                      component,
                      empty_forced_complete,
                      forced_survive,
                      trial_params,
                      trial_key);
  NodeEvalResult res = eval_node_recursive(node_id, state);
  return clamp_probability(res.survival);
}

double compute_guard_free_cluster_value(const std::vector<int>& cluster_indices,
                                        const std::vector<CompetitorMeta>& metas,
                                        const uuber::NativeContext& ctx,
                                        double t,
                                        const std::string& component,
                                        const std::string& trial_key = std::string(),
                                        const TrialParamSet* trial_params = nullptr,
                                        EvalCache* eval_cache_override = nullptr) {
  std::unordered_set<int> empty_forced;
  double prod = 1.0;
  for (int idx : cluster_indices) {
    double surv = evaluate_survival_with_forced(metas[static_cast<std::size_t>(idx)].node_id,
                                                empty_forced,
                                                component,
                                                t,
                                                ctx,
                                                trial_key,
                                                trial_params,
                                                eval_cache_override);
    if (!std::isfinite(surv) || surv <= 0.0) return 0.0;
    prod *= surv;
    if (!std::isfinite(prod) || prod <= 0.0) return 0.0;
  }
  return clamp_probability(prod);
}

double compute_guard_cluster_value(const std::vector<int>& cluster_indices,
                                   const std::vector<CompetitorMeta>& metas,
                                   const uuber::NativeContext& ctx,
                                   double t,
                                   const std::string& component,
                                   const std::string& trial_key = std::string(),
                                   const TrialParamSet* trial_params = nullptr,
                                   EvalCache* eval_cache_override = nullptr,
                                   const std::vector<int>* guard_order = nullptr,
                                   std::unordered_set<int>* forced_survive_scratch = nullptr) {
  const std::vector<int>& order = guard_order ? *guard_order : cluster_indices;
  std::unordered_set<int> local_forced;
  std::unordered_set<int>& forced_survive = forced_survive_scratch ? *forced_survive_scratch : local_forced;
  forced_survive.clear();
  double prod = 1.0;
  for (int idx : order) {
    const CompetitorMeta& node_meta = metas[static_cast<std::size_t>(idx)];
    double surv = evaluate_survival_with_forced(node_meta.node_id,
                                                forced_survive,
                                                component,
                                                t,
                                                ctx,
                                                trial_key,
                                                trial_params,
                                                eval_cache_override);
    if (!std::isfinite(surv) || surv <= 0.0) return 0.0;
    prod *= surv;
    if (!std::isfinite(prod) || prod <= 0.0) return 0.0;
    for (int src : node_meta.sources) {
      forced_survive.insert(src);
    }
  }
  return clamp_probability(prod);
}

double competitor_survival_internal(const uuber::NativeContext& ctx,
                                    const std::vector<int>& competitor_ids,
                                    double t,
                                    const std::string& component_label,
                                    const std::string& trial_type_key = std::string(),
                                    const TrialParamSet* trial_params = nullptr,
                                    EvalCache* eval_cache_override = nullptr) {
  if (competitor_ids.empty()) return 1.0;
  const CompetitorClusterCacheEntry& cache = fetch_competitor_cluster_cache(ctx, competitor_ids);
  if (cache.clusters.empty()) return 1.0;

  double product = 1.0;
  std::unordered_set<int> forced_survive_scratch;
  for (std::size_t i = 0; i < cache.clusters.size(); ++i) {
    const auto& cluster = cache.clusters[i];
    const bool cluster_has_guard = (i < cache.cluster_has_guard.size()) && cache.cluster_has_guard[i];
    double cluster_val = cluster_has_guard
      ? compute_guard_cluster_value(cluster,
                                    cache.metas,
                                    ctx,
                                    t,
                                    component_label,
                                    trial_type_key,
                                    trial_params,
                                    eval_cache_override,
                                    (i < cache.guard_orders.size() ? &cache.guard_orders[i] : nullptr),
                                    &forced_survive_scratch)
      : compute_guard_free_cluster_value(cluster,
                                         cache.metas,
                                         ctx,
                                         t,
                                         component_label,
                                         trial_type_key,
                                         trial_params,
                                         eval_cache_override);
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
  const uuber::NativeContext& ctx,
  int node_id,
  double t,
  const std::string& component_label,
  const std::unordered_set<int>& forced_complete,
  const std::unordered_set<int>& forced_survive,
  const std::vector<int>& competitor_ids,
  const TrialParamSet* trial_params,
  const std::string& trial_type_key = std::string(),
  EvalCache* eval_cache_override = nullptr,
  bool include_na_donors,
  const std::string& outcome_label_context) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  EvalCache local_cache;
  EvalCache& eval_cache = eval_cache_override ? *eval_cache_override : local_cache;
  EvalCache* cache_ptr = eval_cache_override ? eval_cache_override : &eval_cache;
  NodeEvalState state(ctx,
                      eval_cache,
                      t,
                      component_label,
                      forced_complete,
                      forced_survive,
                      trial_params,
                      trial_type_key,
                      include_na_donors,
                      outcome_label_context);
  NodeEvalResult base = eval_node_recursive(node_id, state);
  double density = base.density;
  if (!std::isfinite(density) || density <= 0.0) {
    return 0.0;
  }
  if (!competitor_ids.empty()) {
    double surv = competitor_survival_internal(ctx,
                                               competitor_ids,
                                               t,
                                               component_label,
                                               state.trial_type_key,
                                               trial_params,
                                               cache_ptr);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    density *= surv;
  }
  return density;
}

double native_outcome_probability_impl(SEXP ctxSEXP,
                                       int node_id,
                                       double upper,
                                       Rcpp::Nullable<Rcpp::String> component,
                                       SEXP forced_complete,
                                       SEXP forced_survive,
                                       const Rcpp::IntegerVector& competitor_ids,
                                       double rel_tol,
                                       double abs_tol,
                                       int max_depth,
                                       const TrialParamSet* trial_params,
                                       const std::string& trial_type_key = std::string(),
                                       bool include_na_donors = false,
                                       const std::string& outcome_label_context = std::string()) {
  if (upper <= 0.0) {
    return 0.0;
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull()
    ? Rcpp::as<std::string>(component)
    : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  std::vector<int> comp_vec = integer_vector_to_std(competitor_ids, false);
  auto integrand = [&](double u, const TrialParamSet* params_ptr) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    double val = node_density_with_competitors_internal(
      *ctx,
      node_id,
      u,
      component_label,
      forced_complete_set,
      forced_survive_set,
      comp_vec,
      params_ptr,
      trial_type_key,
      nullptr,
      include_na_donors,
      outcome_label_context
    );
    if (!std::isfinite(val) || val <= 0.0) return 0.0;
    return val;
  };
  auto integrate_with_params = [&](const TrialParamSet* params_ptr) -> double {
    double integral = 0.0;
    if (std::isfinite(upper)) {
      integral = uuber::integrate_boost_fn(
        [&](double u) { return integrand(u, params_ptr); },
        0.0,
        upper,
        rel_tol,
        abs_tol,
        max_depth
      );
    } else {
      auto transformed = [&](double x) -> double {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) x = std::nextafter(1.0, 0.0);
        double t = x / (1.0 - x);
        double jac = 1.0 / ((1.0 - x) * (1.0 - x));
        double val = integrand(t, params_ptr);
        if (!std::isfinite(val) || val <= 0.0) return 0.0;
        double out = val * jac;
        return std::isfinite(out) ? out : 0.0;
      };
      integral = uuber::integrate_boost_fn(
        transformed,
        0.0,
        1.0,
        rel_tol,
        abs_tol,
        max_depth
      );
    }
    if (!std::isfinite(integral)) integral = 0.0;
    return clamp_probability(integral);
  };

  // Ensure we have a concrete parameter set to modify for shared triggers.
  TrialParamSet base_params_holder;
  const TrialParamSet* params_ptr = trial_params;
  if (!params_ptr) {
    base_params_holder = build_base_paramset(*ctx);
    params_ptr = &base_params_holder;
  }

  // Enumerate shared-trigger states to preserve correlation across linked accumulators.
  std::vector<SharedTriggerInfo> triggers = collect_shared_triggers(*ctx, *params_ptr);
  if (triggers.empty()) {
    return integrate_with_params(params_ptr);
  }

  // Joint gate: trigger must succeed once per shared id. When it fails, all linked
  // accumulators are disabled; when it succeeds, they race with q set to 0 to avoid
  // double-counting the gate probability.
  double gate_keep = 1.0;
  for (const auto& trig : triggers) {
    double q = trig.q;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    gate_keep *= (1.0 - q);
  }
  TrialParamSet success_params = apply_trigger_state(*params_ptr, triggers, 0ULL);
  double success_prob = integrate_with_params(&success_params);
  if (!std::isfinite(success_prob) || success_prob < 0.0) success_prob = 0.0;
  double total = gate_keep * success_prob;
  return clamp_probability(total);
}

std::vector<std::string> string_vector_from_entry(SEXP entry) {
  if (Rf_isNull(entry)) return {};
  Rcpp::CharacterVector vec(entry);
  std::vector<std::string> out;
  out.reserve(vec.size());
  for (R_xlen_t i = 0; i < vec.size(); ++i) {
    if (vec[i] == NA_STRING) continue;
    out.push_back(Rcpp::as<std::string>(vec[i]));
  }
  return out;
}

std::unique_ptr<TrialParamSet> build_trial_params_from_df(
  const uuber::NativeContext& ctx,
  const Rcpp::Nullable<Rcpp::DataFrame>& rows_opt) {
  if (rows_opt.isNull()) return nullptr;
  Rcpp::DataFrame rows(rows_opt.get());
  R_xlen_t n = rows.nrows();
  if (n == 0) return nullptr;

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
  Rcpp::NumericVector shared_q_col = rows.containsElementNamed("shared_trigger_q")
    ? Rcpp::NumericVector(rows["shared_trigger_q"])
    : Rcpp::NumericVector();
  Rcpp::CharacterVector shared_col = rows.containsElementNamed("shared_trigger_id")
    ? Rcpp::CharacterVector(rows["shared_trigger_id"])
    : Rcpp::CharacterVector();
  Rcpp::List comps_col = rows.containsElementNamed("components")
    ? Rcpp::List(rows["components"])
    : Rcpp::List();
  Rcpp::List params_list_col = rows.containsElementNamed("params")
    ? Rcpp::List(rows["params"])
    : Rcpp::List();

  std::unordered_set<std::string> base_cols = {
    "trial", "component", "accumulator",
    "type", "outcome", "rt", "params", "onset", "q",
    "condition", "component_weight", "shared_trigger_id",
    "dist", "components"
  };
  std::vector<std::string> param_cols;
  std::vector<SEXP> param_col_data;
  std::vector<int> param_col_types;
  Rcpp::CharacterVector df_names = rows.names();
  if (!df_names.isNULL()) {
    for (R_xlen_t i = 0; i < df_names.size(); ++i) {
      if (df_names[i] == NA_STRING) continue;
      std::string col_name = Rcpp::as<std::string>(df_names[i]);
      if (base_cols.find(col_name) != base_cols.end()) continue;
      SEXP column = rows[col_name];
      param_cols.push_back(col_name);
      param_col_data.push_back(column);
      param_col_types.push_back(TYPEOF(column));
    }
  }

  auto params_set = std::make_unique<TrialParamSet>();
    params_set->acc_params.assign(ctx.accumulators.size(), TrialAccumulatorParams{});

  auto upsert_param_entry = [](std::vector<uuber::ProtoParamEntry>& entries,
                               uuber::ProtoParamEntry entry) {
    for (auto& existing : entries) {
      if (existing.name == entry.name) {
        existing = std::move(entry);
        return;
      }
    }
    entries.push_back(std::move(entry));
  };

  auto append_scalar_entry = [&](std::vector<uuber::ProtoParamEntry>& entries,
                                 const std::string& name,
                                 double value,
                                 bool logical_scalar = false) {
    if (!std::isfinite(value)) return;
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

  auto append_numeric_vector_entry = [&](std::vector<uuber::ProtoParamEntry>& entries,
                                         const std::string& name,
                                         const Rcpp::NumericVector& vec) {
    if (vec.size() == 0) return;
    uuber::ProtoParamEntry entry;
    entry.name = name;
    entry.tag = uuber::ParamValueTag::NumericVector;
    entry.numeric_values.reserve(vec.size());
    for (double val : vec) {
      entry.numeric_values.push_back(val);
    }
    upsert_param_entry(entries, std::move(entry));
  };

  auto append_logical_vector_entry = [&](std::vector<uuber::ProtoParamEntry>& entries,
                                         const std::string& name,
                                         const Rcpp::LogicalVector& vec) {
    if (vec.size() == 0) return;
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
      if (i >= acc_num.size()) continue;
      double raw_val = acc_num[i];
      if (Rcpp::NumericVector::is_na(raw_val)) continue;
      int idx_val = static_cast<int>(std::llround(raw_val)) - 1;
      if (idx_val < 0 || idx_val >= static_cast<int>(ctx.accumulators.size())) continue;
      acc_idx = idx_val;
    }
    if (acc_idx < 0 || acc_idx >= static_cast<int>(ctx.accumulators.size())) continue;
    const uuber::NativeAccumulator& base = ctx.accumulators[acc_idx];

    TrialAccumulatorParams override;
    override.onset = (i < onset_col.size() && !Rcpp::NumericVector::is_na(onset_col[i]))
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
    // Shared gate probability from column or base; per-acc q=0 on success path.
    if (has_shared) {
      double gate_q = base.q;
      if (i < shared_q_col.size() && !Rcpp::NumericVector::is_na(shared_q_col[i])) {
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
      SEXP param_entry = (i < params_list_col.size()) ? params_list_col[i] : R_NilValue;
      if (param_entry != R_NilValue) {
        Rcpp::List param_list(param_entry);
        std::vector<uuber::ProtoParamEntry> entries = params_from_rcpp(param_list, dist_name);
        for (auto& entry : entries) {
          upsert_param_entry(param_entries, entry);
        }
      }
    }
    for (std::size_t pc = 0; pc < param_cols.size(); ++pc) {
      SEXP column = param_col_data[pc];
      int type = param_col_types[pc];
      const std::string& col_name = param_cols[pc];
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
            append_scalar_entry(param_entries, col_name, static_cast<double>(val), false);
          }
        }
        break;
      }
      case LGLSXP: {
        Rcpp::LogicalVector vec(column);
        if (i < vec.size()) {
          int val = vec[i];
          if (val != NA_LOGICAL) {
            append_scalar_entry(param_entries, col_name, static_cast<double>(val), true);
          }
        }
        break;
      }
      case STRSXP:
        // Strings are not supported as distribution params; skip.
        break;
      case VECSXP: {
        Rcpp::List vec(column);
        if (i >= vec.size()) break;
        SEXP cell = vec[i];
        if (cell == R_NilValue) break;
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
              append_scalar_entry(param_entries, col_name, static_cast<double>(val), false);
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
              append_scalar_entry(param_entries, col_name, static_cast<double>(val), true);
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

Rcpp::List native_component_plan_impl(const Rcpp::List& structure,
                                      const Rcpp::DataFrame* trial_rows,
                                      double trial_id,
                                      const std::string* forced_component) {
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
    double weight = (i < base_weights.size()) ? static_cast<double>(base_weights[i]) : 0.0;
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
      for (const auto& comp : all_components) {
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
    // Reference weight = 1 - sum(nonref), then normalize to guard against rounding.
    std::unordered_map<std::string, int> comp_index;
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      comp_index[selected_components[i]] = static_cast<int>(i);
    }
    std::string ref_id = reference.empty() ? selected_components.front() : reference;
    double sum_nonref = 0.0;
    for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
      std::string cid = Rcpp::as<std::string>(component_ids[i]);
      if (cid == ref_id) continue;
      auto it_idx = comp_index.find(cid);
      if (it_idx == comp_index.end()) continue;
      std::string wp;
      if (attrs_list.size() > 0 && i < attrs_list.size()) {
        Rcpp::List attr(attrs_list[i]);
        if (attr.containsElementNamed("weight_param")) {
          wp = Rcpp::as<std::string>(attr["weight_param"]);
        }
      }
      if (wp.empty()) {
        Rcpp::stop("Component '%s' missing weight_param for sampled mixture", cid.c_str());
      }
      double val = std::numeric_limits<double>::quiet_NaN();
      if (trial_rows && trial_rows->containsElementNamed(wp.c_str())) {
        Rcpp::NumericVector col((*trial_rows)[wp]);
        if (col.size() > 0) {
          double cand = col[0];
          if (std::isfinite(cand)) val = cand;
        }
      }
      if (!std::isfinite(val)) {
        auto it = base_weight_map.find(cid);
        if (it != base_weight_map.end()) val = it->second;
      }
      if (!std::isfinite(val) || val < 0.0 || val > 1.0) {
        Rcpp::stop("Mixture weight '%s' must be a probability in [0,1]", wp.c_str());
      }
      weights[static_cast<std::size_t>(it_idx->second)] = val;
      sum_nonref += val;
    }
    auto ref_it = comp_index.find(ref_id);
    if (ref_it != comp_index.end()) {
      double ref_weight = 1.0 - sum_nonref;
      if (ref_weight < -1e-8) {
        Rcpp::stop("Mixture weights sum to >1; check non-reference weight params");
      }
      if (ref_weight < 0.0) ref_weight = 0.0;
      weights[static_cast<std::size_t>(ref_it->second)] = ref_weight;
    }
  } else {
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      auto it = base_weight_map.find(selected_components[i]);
      double w = (it != base_weight_map.end()) ? it->second : 0.0;
      if (!std::isfinite(w) || w < 0.0) w = 0.0;
      weights[i] = w;
    }
  }

  double total = 0.0;
  for (double w : weights) total += w;
  if (total <= 0.0) {
    double uniform = 1.0 / static_cast<double>(weights.size());
    std::fill(weights.begin(), weights.end(), uniform);
  } else {
    for (double& w : weights) {
      w /= total;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("components") = Rcpp::CharacterVector(selected_components.begin(), selected_components.end()),
    Rcpp::Named("weights") = Rcpp::NumericVector(weights.begin(), weights.end())
  );
}

double normalize_weights(std::vector<double>& weights) {
  if (weights.empty()) return 0.0;
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
    for (double& w : weights) {
      w /= total;
    }
    return total;
  }
  double uniform = 1.0 / static_cast<double>(weights.size());
  std::fill(weights.begin(), weights.end(), uniform);
  return 1.0;
}

bool build_component_mix(const Rcpp::CharacterVector& component_ids,
                         const Rcpp::NumericVector& weights_in,
                         Rcpp::Nullable<Rcpp::String> forced_component,
                         std::vector<std::string>& components_out,
                         std::vector<double>& weights_out) {
  if (component_ids.size() == 0) return false;
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
  weights_out.assign(components_out.size(), std::numeric_limits<double>::quiet_NaN());
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

double native_trial_mixture_internal(uuber::NativeContext& ctx,
                                     int node_id,
                                     double t,
                                     const std::vector<std::string>& component_ids,
                                     const std::vector<double>& weights,
                                     Rcpp::Nullable<Rcpp::String> forced_component,
                                     const std::vector<int>& competitor_ids,
                                     const TrialParamSet* params_ptr,
                                     const std::vector<ComponentCacheEntry>& cache_entries,
                                     EvalCache* eval_cache_override = nullptr,
                                     const std::string& keep_label = std::string(),
                                     const Rcpp::List& guess_donors = Rcpp::List()) {
  if (!std::isfinite(t) || t < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  static const std::unordered_set<int> kEmptyForced;
  EvalCache local_eval_cache;
  EvalCache& shared_cache = eval_cache_override ? *eval_cache_override : local_eval_cache;
  double total = 0.0;
  const bool has_keep = !keep_label.empty();
  const std::string keep_norm = has_keep ? normalize_label_string(keep_label) : std::string();
  for (std::size_t i = 0; i < component_ids.size(); ++i) {
    const ComponentCacheEntry* cache_entry_ptr = nullptr;
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
      auto comp_meta = ctx.component_info.find(component_ids[i]);
      if (comp_meta != ctx.component_info.end()) {
        std::string target_norm = normalize_label_string(comp_meta->second.guess.target);
        if (comp_meta->second.guess.valid && !keep_norm.empty() && keep_norm == target_norm) {
          contrib = accumulate_component_guess_density(
            ctx,
            keep_label,
            t,
            component_ids[i],
            params_ptr,
            cache_entry_ptr->trial_type_key,
            &shared_cache
          );
          used_guess_shortcut = true;
        }
      }
    }
    if (!used_guess_shortcut) {
      contrib = node_density_with_competitors_internal(
        ctx,
        node_id,
        t,
        component_ids[i],
        kEmptyForced,
        kEmptyForced,
        competitor_ids,
        params_ptr,
        cache_entry_ptr->trial_type_key,
        &shared_cache,
        false,
        keep_label
      );
      if (!std::isfinite(contrib) || contrib < 0.0) contrib = 0.0;
      if (has_keep) {
        double keep_w = component_keep_weight(ctx, component_ids[i], keep_label);
        contrib *= keep_w;
      }
      if (guess_donors.size() > 0) {
        double donor_contrib = accumulate_plan_guess_density(
          ctx,
          guess_donors,
          t,
          component_ids[i],
          cache_entry_ptr->trial_type_key,
          shared_cache,
          params_ptr,
          false
        );
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
double native_trial_mixture_driver(SEXP ctxSEXP,
                                   int node_id,
                                   double t,
                                   Rcpp::CharacterVector component_ids,
                                   Rcpp::NumericVector weights,
                                   Rcpp::Nullable<Rcpp::String> forced_component,
                                   Rcpp::IntegerVector competitor_ids,
                                   Rcpp::Nullable<Rcpp::DataFrame> trial_rows,
                                   Rcpp::Nullable<Rcpp::List> guess_donors) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::vector<std::string> components;
  std::vector<double> mix_weights;
  if (!build_component_mix(component_ids, weights, forced_component, components, mix_weights)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::vector<int> comp_ids = integer_vector_to_std(competitor_ids, false);
  std::unique_ptr<TrialParamSet> param_holder;
  const TrialParamSet* params_ptr = nullptr;
  if (!trial_rows.isNull()) {
    param_holder = build_trial_params_from_df(*ctx, trial_rows);
    params_ptr = param_holder.get();
  }
  std::vector<ComponentCacheEntry> cache_entries = build_component_cache_entries(components);
  EvalCache trial_eval_cache;
  return native_trial_mixture_internal(*ctx,
                                       node_id,
                                       t,
                                       components,
                                       mix_weights,
                                       forced_component,
                                       comp_ids,
                                       params_ptr,
                                       cache_entries,
                                       &trial_eval_cache,
                                       std::string(),
                                       guess_donors.isNotNull() ? Rcpp::List(guess_donors.get()) : Rcpp::List());
}

struct ComponentMap {
  std::vector<std::string> ids;
  std::vector<int> leader_idx;
  std::vector<double> base_weights;
};

int col_index(const Rcpp::CharacterVector& cols, const std::string& name) {
  for (R_xlen_t i = 0; i < cols.size(); ++i) {
    Rcpp::String s(cols[i]);
    if (s == NA_STRING) continue;
    if (name == std::string(s.get_cstring())) return static_cast<int>(i);
  }
  return -1;
}

ComponentMap build_component_map(const uuber::NativeContext& ctx, const Rcpp::List& structure) {
  ComponentMap map;
  if (structure.containsElementNamed("components")) {
    Rcpp::DataFrame comp_df(structure["components"]);
    Rcpp::CharacterVector ids = comp_df["component_id"];
    Rcpp::NumericVector weights = comp_df["weight"];
    map.ids.reserve(ids.size());
    map.base_weights.reserve(ids.size());
    for (R_xlen_t i = 0; i < ids.size(); ++i) {
      if (ids[i] == NA_STRING) continue;
      map.ids.push_back(Rcpp::as<std::string>(ids[i]));
      double w = (i < weights.size()) ? static_cast<double>(weights[i]) : 1.0;
      map.base_weights.push_back(std::isfinite(w) && w >= 0.0 ? w : 1.0);
    }
  }
  if (map.ids.empty()) {
    map.ids.push_back("__default__");
    map.base_weights.push_back(1.0);
  }
  map.leader_idx.assign(map.ids.size(), -1);
  for (std::size_t c = 0; c < map.ids.size(); ++c) {
    const std::string& cid = map.ids[c];
    for (std::size_t a = 0; a < ctx.accumulators.size(); ++a) {
      const auto& comps = ctx.accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        map.leader_idx[c] = static_cast<int>(a);
        break;
      }
    }
    if (map.leader_idx[c] < 0 && !ctx.accumulators.empty()) {
      map.leader_idx[c] = 0;
    }
  }
  return map;
}

std::vector<double> component_weights_for_trial(const ComponentMap& map,
                                                const Rcpp::NumericMatrix& params,
                                                int trial_idx,
                                                int n_acc,
                                                int w_col) {
  std::vector<double> weights(map.ids.size(), 0.0);
  for (std::size_t i = 0; i < map.ids.size(); ++i) {
    double w = map.base_weights[i];
    int leader = map.leader_idx[i];
    if (leader >= 0 && w_col >= 0) {
      int row_idx = trial_idx * n_acc + leader;
      if (row_idx >= 0 && row_idx < params.nrow()) {
        double cand = params(row_idx, w_col);
        if (std::isfinite(cand) && cand >= 0.0) w = cand;
      }
    }
    if (!std::isfinite(w) || w < 0.0) w = 0.0;
    weights[i] = w;
  }
  double total = 0.0;
  for (double v : weights) total += v;
  if (!std::isfinite(total) || total <= 0.0) {
    double uni = 1.0 / static_cast<double>(weights.size());
    std::fill(weights.begin(), weights.end(), uni);
  } else {
    for (double& v : weights) v /= total;
  }
  return weights;
}

std::vector<TrialParamSet> params_from_matrix(const uuber::NativeContext& ctx,
                                              const Rcpp::NumericMatrix& params_mat,
                                              int n_trials) {
  int n_acc = static_cast<int>(ctx.accumulators.size());
  Rcpp::RObject dimnames = params_mat.attr("dimnames");
  Rcpp::CharacterVector cols;
  if (!dimnames.isNULL()) {
    Rcpp::List dn(dimnames);
    if (dn.size() > 1) {
      cols = Rcpp::as<Rcpp::CharacterVector>(dn[1]);
    }
  }
  int q_col = col_index(cols, "q");
  int t0_col = col_index(cols, "t0");
  int p1_col = col_index(cols, "p1");
  int p2_col = col_index(cols, "p2");
  int p3_col = col_index(cols, "p3");
  if (q_col < 0 || t0_col < 0 || p1_col < 0) {
    Rcpp::stop("Parameter matrix must include columns q, t0, and p1");
  }
  std::vector<TrialParamSet> out(static_cast<std::size_t>(n_trials));
  for (int t = 0; t < n_trials; ++t) {
    TrialParamSet ps;
    ps.acc_params.reserve(ctx.accumulators.size());
    for (int a = 0; a < n_acc; ++a) {
      int row_idx = t * n_acc + a;
    TrialAccumulatorParams tap = base_params(ctx.accumulators[a]);
    tap.q = clamp_probability(params_mat(row_idx, q_col));
    tap.shared_q = tap.q;
    tap.shared_trigger_id = ctx.accumulators[a].shared_trigger_id;
    tap.dist_cfg.t0 = params_mat(row_idx, t0_col);
    tap.has_components = !ctx.accumulators[a].components.empty();
    tap.components = ctx.accumulators[a].components;
    tap.has_override = true;
      if (tap.dist_cfg.code == uuber::ACC_DIST_LOGNORMAL ||
          tap.dist_cfg.code == uuber::ACC_DIST_GAMMA ||
          tap.dist_cfg.code == uuber::ACC_DIST_EXGAUSS) {
        tap.dist_cfg.p1 = params_mat(row_idx, p1_col);
        if (p2_col >= 0) tap.dist_cfg.p2 = params_mat(row_idx, p2_col);
        if (p3_col >= 0) tap.dist_cfg.p3 = params_mat(row_idx, p3_col);
      }
      ps.acc_params.push_back(tap);
    }
    out[static_cast<std::size_t>(t)] = std::move(ps);
  }
  return out;
}

double mix_outcome_prob(SEXP ctxSEXP,
                        const uuber::OutcomeContextInfo& info,
                        double rt,
                        const std::vector<std::string>& comp_labels,
                        const std::vector<double>& comp_weights,
                        const TrialParamSet* params_ptr,
                        double rel_tol,
                        double abs_tol,
                        int max_depth,
                        const std::string& outcome_label_context) {
  if (info.node_id < 0) return 0.0;
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::vector<ComponentCacheEntry> cache_entries = build_component_cache_entries(comp_labels);
  EvalCache eval_cache;
  double density = native_trial_mixture_internal(
    *ctx,
    info.node_id,
    rt,
    comp_labels,
    comp_weights,
    R_NilValue,
    info.competitor_ids,
    params_ptr,
    cache_entries,
    &eval_cache,
    outcome_label_context,
    Rcpp::List()
  );
  if (!std::isfinite(density) || density < 0.0) return 0.0;
  return density;
}

// [[Rcpp::export]]
Rcpp::List cpp_loglik(SEXP ctxSEXP,
                      Rcpp::List structure,
                      Rcpp::NumericMatrix params_mat,
                      Rcpp::DataFrame data_df,
                      double rel_tol,
                      double abs_tol,
                      int max_depth) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  int n_acc = static_cast<int>(ctx->accumulators.size());
  if (n_acc == 0) Rcpp::stop("No accumulators in native context");
  if (params_mat.nrow() % n_acc != 0) {
    Rcpp::stop("Parameter rows must be a multiple of the number of accumulators");
  }
  int n_trials = params_mat.nrow() / n_acc;
  ComponentMap comp_map = build_component_map(*ctx, structure);
  Rcpp::RObject dimnames = params_mat.attr("dimnames");
  Rcpp::CharacterVector col_names;
  if (!dimnames.isNULL()) {
    Rcpp::List dn(dimnames);
    if (dn.size() > 1) {
      col_names = Rcpp::as<Rcpp::CharacterVector>(dn[1]);
    }
  }
  int w_col = col_index(col_names, "w");
  if (w_col < 0) Rcpp::stop("Parameter matrix missing column 'w'");
  std::vector<TrialParamSet> param_sets = params_from_matrix(*ctx, params_mat, n_trials);

  Rcpp::NumericVector trial_col = data_df["trial"];
  Rcpp::CharacterVector outcome_col = data_df["outcome"];
  Rcpp::NumericVector rt_col = data_df["rt"];
  bool has_component = data_df.containsElementNamed("component");
  Rcpp::CharacterVector comp_col = has_component ? Rcpp::CharacterVector(data_df["component"]) : Rcpp::CharacterVector();

  R_xlen_t n_obs = outcome_col.size();
  Rcpp::NumericVector per_trial(n_obs);
  double total_loglik = 0.0;
  bool hit_neg_inf = false;

  for (R_xlen_t i = 0; i < n_obs; ++i) {
    int trial_idx = static_cast<int>(trial_col[i]);
    if (trial_idx == NA_INTEGER) {
      Rcpp::stop("Data frame contains NA trial index");
    }
    trial_idx -= 1;
    if (trial_idx < 0 || trial_idx >= n_trials) {
      Rcpp::stop("Trial index out of range for parameter matrix");
    }
    TrialParamSet* params_ptr = &param_sets[static_cast<std::size_t>(trial_idx)];
    double rt = rt_col[i];
    Rcpp::String out_str(outcome_col[i]);
    bool outcome_is_na = (out_str == NA_STRING);
    std::string outcome_label = outcome_is_na ? std::string() : std::string(out_str.get_cstring());
    std::string forced_component;
    if (has_component && comp_col.size() > i) {
      Rcpp::String comp_val(comp_col[i]);
      if (comp_val != NA_STRING) {
        forced_component = std::string(comp_val.get_cstring());
      }
    }

    std::vector<std::string> comp_labels;
    std::vector<double> comp_weights;
    if (!forced_component.empty()) {
      comp_labels.push_back(forced_component);
      comp_weights.push_back(1.0);
    } else {
      comp_labels = comp_map.ids;
      comp_weights = component_weights_for_trial(comp_map, params_mat, trial_idx, n_acc, w_col);
    }

    double prob = 0.0;
    if (outcome_is_na) {
      double na_explicit = 0.0;
      double non_na = 0.0;
      for (const auto& kv : ctx->outcome_info) {
        const std::string& lbl = kv.first;
        if (lbl == "__DEADLINE__") continue;
        const uuber::OutcomeContextInfo& info = kv.second;
        double p = mix_outcome_prob(ctxSEXP,
                                    info,
                                    rt,
                                    comp_labels,
                                    comp_weights,
                                    params_ptr,
                                    rel_tol,
                                    abs_tol,
                                    max_depth,
                                    lbl);
        if (!std::isfinite(p) || p <= 0.0) continue;
        if (info.maps_to_na) {
          na_explicit += p;
        } else {
          non_na += p;
        }
      }
      double residual = 1.0 - na_explicit - non_na;
      if (!std::isfinite(residual)) residual = 0.0;
      if (residual < 0.0) residual = 0.0;
      prob = clamp_probability(na_explicit + residual);
    } else {
      auto it = ctx->outcome_info.find(outcome_label);
      if (it == ctx->outcome_info.end() || it->second.node_id < 0) {
        prob = 0.0;
      } else if (!std::isfinite(rt) || rt < 0.0) {
        prob = 0.0;
      } else {
        prob = mix_outcome_prob(ctxSEXP,
                                it->second,
                                rt,
                                comp_labels,
                                comp_weights,
                                params_ptr,
                                rel_tol,
                                abs_tol,
                                max_depth,
                                outcome_label);
      }
    }

    double log_contrib = R_NegInf;
    if (std::isfinite(prob) && prob > 0.0) {
      log_contrib = std::log(prob);
    } else {
      hit_neg_inf = true;
    }
    per_trial[i] = log_contrib;
    if (!hit_neg_inf && std::isfinite(log_contrib)) {
      total_loglik += log_contrib;
    }
  }

  if (hit_neg_inf) {
    total_loglik = R_NegInf;
  }

  return Rcpp::List::create(
    Rcpp::Named("loglik") = total_loglik,
    Rcpp::Named("per_trial") = per_trial
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_loglik_multiple(SEXP ctxSEXP,
                                        Rcpp::List structure,
                                        Rcpp::List params_list,
                                        Rcpp::DataFrame data_df,
                                        double rel_tol,
                                        double abs_tol,
                                        int max_depth) {
  Rcpp::NumericVector out(params_list.size());
  for (R_xlen_t i = 0; i < params_list.size(); ++i) {
    if (params_list[i] == R_NilValue) {
      out[i] = NA_REAL;
      continue;
    }
    Rcpp::NumericMatrix pm(params_list[i]);
    Rcpp::List res = cpp_loglik(ctxSEXP, structure, pm, data_df, rel_tol, abs_tol, max_depth);
    out[i] = Rcpp::as<double>(res["loglik"]);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List native_component_plan_exported(SEXP structureSEXP,
                                          SEXP trial_rowsSEXP,
                                          SEXP forced_componentSEXP) {
  Rcpp::List structure = Rcpp::as<Rcpp::List>(structureSEXP);
  Rcpp::DataFrame trial_rows_df;
  const Rcpp::DataFrame* trial_rows_ptr = nullptr;
  if (!Rf_isNull(trial_rowsSEXP)) {
    trial_rows_df = Rcpp::DataFrame(trial_rowsSEXP);
    trial_rows_ptr = &trial_rows_df;
  }
  double trial_id = extract_trial_id(trial_rows_ptr);
  std::string forced_component_value;
  const std::string* forced_component_ptr = nullptr;
  if (!Rf_isNull(forced_componentSEXP)) {
    forced_component_value = Rcpp::as<std::string>(forced_componentSEXP);
    forced_component_ptr = &forced_component_value;
  }
  return native_component_plan_impl(structure,
                                    trial_rows_ptr,
                                    trial_id,
                                    forced_component_ptr);
}

// [[Rcpp::export]]
double native_outcome_probability_params_cpp(SEXP ctxSEXP,
                                             int node_id,
                                             double upper,
                                             Rcpp::Nullable<Rcpp::String> component,
                                             SEXP forced_complete,
                                             SEXP forced_survive,
                                             Rcpp::IntegerVector competitor_ids,
                                             double rel_tol,
                                             double abs_tol,
                                             int max_depth,
                                             Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder = build_trial_params_from_df(*ctx, trial_rows);
  return native_outcome_probability_impl(ctxSEXP,
                                         node_id,
                                         upper,
                                         component,
                                         forced_complete,
                                         forced_survive,
                                         competitor_ids,
                                         rel_tol,
                                         abs_tol,
                                         max_depth,
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
  R_xlen_t idx = 0;
  for (const auto& kv : ctx->outcome_info) {
    if (idx >= n) break;
    labels[idx] = kv.first;
    node_ids[idx] = kv.second.node_id;
    comp_counts[idx] = static_cast<int>(kv.second.competitor_ids.size());
    maps_na[idx] = kv.second.maps_to_na;
    ++idx;
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("label") = labels,
    Rcpp::Named("node_id") = node_ids,
    Rcpp::Named("competitors") = comp_counts,
    Rcpp::Named("maps_to_na") = maps_na,
    Rcpp::Named("stringsAsFactors") = false
  );
}

// Forward declarations for native distribution helpers
Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog);
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog);
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate);
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate);
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau);
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau);

namespace {

inline bool is_invalid_positive(double value) {
  return !std::isfinite(value) || value <= 0.0;
}

inline std::vector<double> expand_poly(const std::vector<double>& coeff,
                                       double surv,
                                       double fail) {
  std::size_t len = coeff.size();
  std::vector<double> out(len + 1, 0.0);
  for (std::size_t i = 0; i < len; ++i) {
    double base = coeff[i];
    if (base == 0.0) continue;
    out[i] += base * surv;
    out[i + 1] += base * fail;
  }
  return out;
}

struct PoolPolyScratch {
  std::vector<std::vector<double>> prefix;
  std::vector<std::vector<double>> suffix;
};

inline void ensure_prefix_shape(std::vector<std::vector<double>>& prefix,
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

inline void ensure_suffix_shape(std::vector<std::vector<double>>& suffix,
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

inline PoolPolyScratch& fetch_pool_poly_scratch(std::size_t n,
                                                bool need_suffix) {
  thread_local std::unordered_map<std::size_t, PoolPolyScratch> scratch_map;
  PoolPolyScratch& scratch = scratch_map[n];
  ensure_prefix_shape(scratch.prefix, n);
  if (need_suffix) {
    ensure_suffix_shape(scratch.suffix, n);
  }
  return scratch;
}

inline void fill_prefix_buffers(std::vector<std::vector<double>>& prefix,
                                const std::vector<double>& surv,
                                const std::vector<double>& fail) {
  const std::size_t n = surv.size();
  std::fill(prefix[0].begin(), prefix[0].end(), 0.0);
  prefix[0][0] = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::vector<double>& prev = prefix[i];
    std::vector<double>& out = prefix[i + 1];
    std::fill(out.begin(), out.end(), 0.0);
    double surv_i = surv[i];
    double fail_i = fail[i];
    for (std::size_t j = 0; j < prev.size(); ++j) {
      double base = prev[j];
      if (base == 0.0) continue;
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline void fill_suffix_buffers(std::vector<std::vector<double>>& suffix,
                                const std::vector<double>& surv,
                                const std::vector<double>& fail) {
  const std::size_t n = surv.size();
  std::fill(suffix[n].begin(), suffix[n].end(), 0.0);
  suffix[n][0] = 1.0;
  for (std::size_t idx = n; idx-- > 0;) {
    const std::vector<double>& next = suffix[idx + 1];
    std::vector<double>& out = suffix[idx];
    std::fill(out.begin(), out.end(), 0.0);
    double surv_i = surv[idx];
    double fail_i = fail[idx];
    for (std::size_t j = 0; j < next.size(); ++j) {
      double base = next[j];
      if (base == 0.0) continue;
      out[j] += base * surv_i;
      out[j + 1] += base * fail_i;
    }
  }
}

inline double coefficient_for_order(const std::vector<double>& pref,
                                    const std::vector<double>& suff,
                                    int order) {
  if (order < 0) return 0.0;
  double total = 0.0;
  int max_pref = static_cast<int>(pref.size()) - 1;
  int max_index = std::min(order, max_pref);
  for (int i = 0; i <= max_index; ++i) {
    int s_idx = order - i;
    if (s_idx < 0 || s_idx >= static_cast<int>(suff.size())) continue;
    total += pref[static_cast<std::size_t>(i)] * suff[static_cast<std::size_t>(s_idx)];
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

inline void sort_unique(std::vector<int>& vec) {
  if (vec.empty()) return;
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

inline std::vector<int> integer_vector_to_std(const Rcpp::IntegerVector& src,
                                              bool subtract_one) {
  std::vector<int> out;
  out.reserve(src.size());
  for (int val : src) {
    if (val == NA_INTEGER) continue;
    out.push_back(subtract_one ? val - 1 : val);
  }
  sort_unique(out);
  return out;
}

inline std::vector<int> merge_forced_vectors(const std::vector<int>& base,
                                             const std::vector<int>& addition) {
  std::vector<int> result = base;
  for (int val : addition) {
    if (val == NA_INTEGER) continue;
    result.push_back(val);
  }
  sort_unique(result);
  return result;
}

void combinations_recursive(const std::vector<int>& elements,
                            int choose,
                            std::size_t start,
                            std::vector<int>& current,
                            std::vector<std::vector<int>>& output) {
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

inline std::vector<std::vector<int>> generate_combinations(const std::vector<int>& elements,
                                                           int choose) {
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

inline std::vector<int> survivors_from_combo(const std::vector<int>& others,
                                             const std::vector<int>& combo) {
  if (combo.empty()) return others;
  std::vector<int> survivors;
  std::vector<int> combo_sorted = combo;
  std::sort(combo_sorted.begin(), combo_sorted.end());
  survivors.reserve(others.size());
  for (int candidate : others) {
    if (!std::binary_search(combo_sorted.begin(), combo_sorted.end(), candidate)) {
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
Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog) {
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
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector& x,
                                       double meanlog,
                                       double sdlog) {
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
Rcpp::NumericVector dist_lognormal_rng(int n,
                                       double meanlog,
                                       double sdlog) {
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
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  const double scale = 1.0 / rate;
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
    out[i] = R::dgamma(xi, shape, scale, /*give_log =*/0);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector& x,
                                   double shape,
                                   double rate) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  const double scale = 1.0 / rate;
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
    out[i] = R::pgamma(xi, shape, scale, /*lower_tail =*/1, /*log_p =*/0);
    out[i] = clamp(out[i], 0.0, 1.0);
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_rng(int n,
                                   double shape,
                                   double rate) {
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
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }

  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;

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
    const double z = (xi - mu) / sigma;
    const double exponent = half_sigma_sq_over_tau_sq - (xi - mu) * inv_tau;
    const double gaussian_tail = R::pnorm(z - sigma_over_tau, 0.0, 1.0, 1, 0);
    double val = inv_tau * std::exp(exponent) * gaussian_tail;
    if (!std::isfinite(val) || val < 0.0) {
      val = 0.0;
    }
    out[i] = val;
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector& x,
                                     double mu,
                                     double sigma,
                                     double tau) {
  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }

  const double inv_tau = 1.0 / tau;
  const double sigma_sq = sigma * sigma;
  const double tau_sq = tau * tau;
  const double half_sigma_sq_over_tau_sq = sigma_sq / (2.0 * tau_sq);
  const double sigma_over_tau = sigma * inv_tau;

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
    const double z = (xi - mu) / sigma;
    const double base_cdf = R::pnorm(z, 0.0, 1.0, 1, 0);
    const double tail = R::pnorm(z - sigma_over_tau, 0.0, 1.0, 1, 0);
    const double exp_term = std::exp(half_sigma_sq_over_tau_sq - (xi - mu) * inv_tau);
    double val = base_cdf - exp_term * tail;
    if (!std::isfinite(val)) {
      val = NA_REAL;
    } else {
      val = clamp(val, 0.0, 1.0);
    }
    out[i] = val;
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_rng(int n,
                                     double mu,
                                     double sigma,
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

struct PoolDensityWorker : public RcppParallel::Worker {
  const RcppParallel::RVector<double> density;
  const std::vector<std::vector<double>>* prefix;
  const std::vector<std::vector<double>>* suffix;
  const int order;
  double total;

  PoolDensityWorker(const Rcpp::NumericVector& densityVec,
                    const std::vector<std::vector<double>>& prefixRef,
                    const std::vector<std::vector<double>>& suffixRef,
                    int order_)
    : density(densityVec),
      prefix(&prefixRef),
      suffix(&suffixRef),
      order(order_),
      total(0.0) {}

  PoolDensityWorker(const PoolDensityWorker& other, RcppParallel::Split)
    : density(other.density),
      prefix(other.prefix),
      suffix(other.suffix),
      order(other.order),
      total(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local = 0.0;
    const auto& pref = *prefix;
    const auto& suff = *suffix;
    for (std::size_t idx = begin; idx < end; ++idx) {
      double dens_val = density[idx];
      if (!std::isfinite(dens_val) || dens_val <= 0.0) continue;
      double coeff = coefficient_for_order(pref[idx], suff[idx + 1], order);
      if (!std::isfinite(coeff) || coeff <= 0.0) continue;
      local += dens_val * coeff;
    }
    total += local;
  }

  void join(const PoolDensityWorker& rhs) {
    total += rhs.total;
  }
};

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector pool_coeffs_cpp(const Rcpp::NumericVector& Svec,
                                    const Rcpp::NumericVector& Fvec) {
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

double pool_density_fast_cpp(const Rcpp::NumericVector& density,
                             const Rcpp::NumericVector& survival,
                             int k) {
  const std::size_t n = density.size();
  if (n == 0) return 0.0;
  if (k < 1) return 0.0;
  if (k > static_cast<int>(n)) return 0.0;

  std::vector<double> surv(n);
  std::vector<double> fail(n);
  for (std::size_t i = 0; i < n; ++i) {
    double s = clamp_unit(survival[i]);
    surv[i] = s;
    fail[i] = clamp_unit(1.0 - s);
  }

  PoolPolyScratch& scratch = fetch_pool_poly_scratch(n, true);
  fill_prefix_buffers(scratch.prefix, surv, fail);
  fill_suffix_buffers(scratch.suffix, surv, fail);

  PoolDensityWorker worker(density, scratch.prefix, scratch.suffix, k - 1);
#if defined(RCPP_PARALLEL_USE_TBB) && RCPP_PARALLEL_USE_TBB
  RcppParallel::parallelReduce(0, n, worker);
#else
  worker(0, n);
#endif

  double total = worker.total;
  if (!std::isfinite(total) || total < 0.0) return 0.0;
  return total;
}

double pool_survival_fast_cpp(const Rcpp::NumericVector& survival,
                              int k) {
  const std::size_t n = survival.size();
  if (n == 0) return 1.0;
  if (k < 1) return 0.0;

  std::vector<double> surv(n);
  std::vector<double> fail(n);
  for (std::size_t i = 0; i < n; ++i) {
    double s = clamp_unit(survival[i]);
    surv[i] = s;
    fail[i] = clamp_unit(1.0 - s);
  }

  PoolPolyScratch& scratch = fetch_pool_poly_scratch(n, false);
  fill_prefix_buffers(scratch.prefix, surv, fail);
  const std::vector<double>& poly = scratch.prefix.back();
  int upto = std::min<int>(k, static_cast<int>(poly.size()));
  if (upto <= 0) return 0.0;
  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += poly[static_cast<std::size_t>(i)];
  }
  return clamp(total, 0.0, 1.0);
}

Rcpp::List pool_build_templates_cpp(int n,
                                    const Rcpp::IntegerVector& member_ids,
                                    int pool_idx,
                                    int k) {
  std::vector<int> member_std(member_ids.begin(), member_ids.end());
  uuber::PoolTemplateCacheEntry cache = build_pool_template_cache(n, member_std, pool_idx, k);
  Rcpp::List finisher_map(n);
  Rcpp::List templates_out(cache.templates.size());
  for (std::size_t i = 0; i < cache.templates.size(); ++i) {
    const uuber::PoolTemplateEntry& tpl = cache.templates[i];
    Rcpp::IntegerVector complete_idx;
    if (!tpl.complete_idx.empty()) {
      complete_idx = Rcpp::IntegerVector(tpl.complete_idx.begin(), tpl.complete_idx.end());
      for (auto& val : complete_idx) {
        val += 1;
      }
    } else {
      complete_idx = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector survivor_idx;
    if (!tpl.survivor_idx.empty()) {
      survivor_idx = Rcpp::IntegerVector(tpl.survivor_idx.begin(), tpl.survivor_idx.end());
      for (auto& val : survivor_idx) {
        val += 1;
      }
    } else {
      survivor_idx = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector forced_complete_vec;
    if (!tpl.forced_complete_ids.empty()) {
      forced_complete_vec = Rcpp::IntegerVector(tpl.forced_complete_ids.begin(), tpl.forced_complete_ids.end());
    } else {
      forced_complete_vec = Rcpp::IntegerVector(0);
    }
    Rcpp::IntegerVector forced_survive_vec;
    if (!tpl.forced_survive_ids.empty()) {
      forced_survive_vec = Rcpp::IntegerVector(tpl.forced_survive_ids.begin(), tpl.forced_survive_ids.end());
    } else {
      forced_survive_vec = Rcpp::IntegerVector(0);
    }
    templates_out[i] = Rcpp::List::create(
      Rcpp::Named("finisher_idx") = tpl.finisher_idx + 1,
      Rcpp::Named("complete_idx") = complete_idx,
      Rcpp::Named("survivor_idx") = survivor_idx,
      Rcpp::Named("forced_complete_ids") = forced_complete_vec,
      Rcpp::Named("forced_survive_ids") = forced_survive_vec
    );
  }
  for (int idx = 0; idx < n; ++idx) {
    const std::vector<int>& entries = cache.finisher_map[static_cast<std::size_t>(idx)];
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
  return Rcpp::List::create(
    Rcpp::Named("templates") = templates_out,
    Rcpp::Named("finisher_map") = finisher_map
  );
}

Rcpp::List pool_density_combine_native(const Rcpp::NumericVector& dens_vec,
                                       const Rcpp::NumericVector& cdf_vec,
                                       const Rcpp::NumericVector& surv_vec,
                                       const Rcpp::NumericVector& cdf_success_vec,
                                       const Rcpp::NumericVector& surv_success_vec,
                                       const std::vector<int>& shared_index,
                                       const std::vector<uuber::PoolTemplateEntry>& templates,
                                       const std::vector<int>& base_forced_complete,
                                       const std::vector<int>& base_forced_survive) {
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

  for (const auto& tpl : templates) {
    int finisher_idx = tpl.finisher_idx;
    if (finisher_idx < 0 || finisher_idx >= dens_vec.size()) continue;
    double dens_mid = dens_vec[finisher_idx];
    if (!std::isfinite(dens_mid) || dens_mid <= 0.0) continue;
    double weight = dens_mid;

    for (int j : tpl.complete_idx) {
      if (j < 0 || j >= cdf_vec.size()) continue;
      weight *= cdf_vec[j];
    }
    for (int j : tpl.survivor_idx) {
      if (j < 0 || j >= surv_vec.size()) continue;
      weight *= surv_vec[j];
    }
    if (!std::isfinite(weight) || weight <= 0.0) continue;

    for (int j : tpl.complete_idx) {
      if (j < 0 || j >= cdf_vec.size()) continue;
      if (!same_shared(finisher_idx, j)) continue;
      double denom = std::max(cdf_vec[j], eps);
      double ratio = cdf_success_vec[j] / denom;
      weight *= ratio;
    }
    for (int j : tpl.survivor_idx) {
      if (j < 0 || j >= surv_vec.size()) continue;
      if (!same_shared(finisher_idx, j)) continue;
      double denom = std::max(surv_vec[j], eps);
      double ratio = surv_success_vec[j] / denom;
      weight *= ratio;
    }
    if (!std::isfinite(weight) || weight <= 0.0) continue;

    std::vector<int> merged_fc = merge_forced_vectors(base_forced_complete, tpl.forced_complete_ids);
    std::vector<int> merged_fs = merge_forced_vectors(base_forced_survive, tpl.forced_survive_ids);

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

    scenario_vec.emplace_back(
      Rcpp::List::create(
        Rcpp::Named("weight") = weight,
        Rcpp::Named("forced_complete") = fc_vec,
        Rcpp::Named("forced_survive") = fs_vec
      )
    );
  }

  Rcpp::List scenarios_out(scenario_vec.size());
  for (std::size_t i = 0; i < scenario_vec.size(); ++i) {
    scenarios_out[i] = scenario_vec[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = total,
    Rcpp::Named("scenarios") = scenarios_out
  );
}

Rcpp::List pool_density_combine_cpp(const Rcpp::NumericVector& dens_vec,
                                    const Rcpp::NumericVector& cdf_vec,
                                    const Rcpp::NumericVector& surv_vec,
                                    const Rcpp::NumericVector& cdf_success_vec,
                                    const Rcpp::NumericVector& surv_success_vec,
                                    const Rcpp::IntegerVector& shared_index,
                                    const Rcpp::List& templates,
                                    const Rcpp::IntegerVector& forced_complete,
                                    const Rcpp::IntegerVector& forced_survive) {
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
    entry.forced_complete_ids = integer_vector_to_std(tpl["forced_complete_ids"], false);
    entry.forced_survive_ids = integer_vector_to_std(tpl["forced_survive_ids"], false);
    native_templates.push_back(std::move(entry));
  }
  return pool_density_combine_native(dens_vec,
                                     cdf_vec,
                                     surv_vec,
                                     cdf_success_vec,
                                     surv_success_vec,
                                     shared_group,
                                     native_templates,
                                     base_fc,
                                     base_fs);
}

double pool_survival_general_cpp(const Rcpp::NumericVector& Fvec,
                                 int k) {
  const std::size_t n = Fvec.size();
  if (n == 0) return 1.0;
  if (k < 1) return 0.0;

  Rcpp::NumericVector Fclamp(n);
  Rcpp::NumericVector Svec(n);
  for (std::size_t i = 0; i < n; ++i) {
    double val = Fvec[i];
    if (!std::isfinite(val)) val = 0.0;
    val = clamp(val, 0.0, 1.0);
    Fclamp[i] = val;
    Svec[i] = 1.0 - val;
  }

  Rcpp::NumericVector coeffs = pool_coeffs_cpp(Svec, Fclamp);
  int upto = std::min<int>(coeffs.size(), k);
  if (upto <= 0) return 0.0;
  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += coeffs[i];
  }
  return clamp(total, 0.0, 1.0);
}

double guard_effective_survival_cpp(Rcpp::Function integrand,
                                    double upper,
                                    double rel_tol,
                                    double abs_tol,
                                    int max_depth) {
  if (!std::isfinite(upper)) return 0.0;
  if (upper <= 0.0) return 1.0;
  if (rel_tol <= 0.0) rel_tol = 1e-5;
  if (abs_tol <= 0.0) abs_tol = 1e-6;
  double integral = 0.0;
  try {
    integral = uuber::integrate_boost(integrand, 0.0, upper, rel_tol, abs_tol, max_depth);
  } catch (...) {
    integral = 0.0;
  }
  if (!std::isfinite(integral) || integral < 0.0) integral = 0.0;
  double surv = 1.0 - integral;
  if (!std::isfinite(surv)) surv = 0.0;
  if (surv < 0.0) surv = 0.0;
  if (surv > 1.0) surv = 1.0;
  return surv;
}
