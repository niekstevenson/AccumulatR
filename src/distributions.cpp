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
#include <sstream>
#include <numeric>
#include <iomanip>
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#define UUBER_HAVE_BOOST_GK 1
#include <boost/math/quadrature/gauss_kronrod.hpp>
#else
#define UUBER_HAVE_BOOST_GK 0
#endif

#include "native_context.h"
#include "native_accumulator.hpp"
#include "native_prep_builder.hpp"
#include "native_proto.hpp"
#include "native_integrate.hpp"

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
  if (blob.size() == 0) {
    Rcpp::stop("native_context_from_proto: buffer is empty");
  }
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
std::vector<uuber::ProtoParamEntry> params_from_rcpp(const Rcpp::List& params,
                                                     const std::string& dist) {
  if (params.size() == 0) {
    Rcpp::stop("Distribution '%s' received empty parameter list", dist);
  }
  Rcpp::CharacterVector names = params.names();
  if (names.isNULL()) {
    Rcpp::stop("Distribution '%s' parameter list must be named", dist);
  }
  std::vector<uuber::ProtoParamEntry> out;
  out.reserve(params.size());
  for (R_xlen_t i = 0; i < params.size(); ++i) {
    if (params[i] == R_NilValue) continue;
    if (names[i] == NA_STRING) {
      Rcpp::stop("Distribution '%s' missing parameter name at index %d",
                 dist, static_cast<int>(i + 1));
    }
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
        entry.logical_values.reserve(lv.size());
        for (int v : lv) entry.logical_values.push_back(v);
      }
    } else if (Rf_isInteger(val) || Rf_isReal(val)) {
      Rcpp::NumericVector nv(val);
      if (nv.size() == 1) {
        entry.tag = uuber::ParamValueTag::NumericScalar;
        entry.numeric_scalar = nv[0];
      } else {
        entry.tag = uuber::ParamValueTag::NumericVector;
        entry.numeric_values.reserve(nv.size());
        for (double v : nv) entry.numeric_values.push_back(v);
      }
    } else {
      Rcpp::stop("Distribution '%s' parameter '%s' must be numeric/logical",
                 dist, entry.name);
    }
    out.push_back(std::move(entry));
  }
  return out;
}
Rcpp::List pool_build_templates_cpp(int n,
                                    const Rcpp::IntegerVector& member_ids,
                                    int pool_idx,
                                    int k);
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

inline double eval_pdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return dist_lognormal_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_GAMMA:
    return dist_gamma_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_EXGAUSS:
    return dist_exgauss_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    Rcpp::stop("Invalid accumulator distribution code '%d'", cfg.code);
  }
}

inline double eval_cdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case uuber::ACC_DIST_LOGNORMAL:
    return dist_lognormal_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_GAMMA:
    return dist_gamma_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case uuber::ACC_DIST_EXGAUSS:
    return dist_exgauss_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    Rcpp::stop("Invalid accumulator distribution code '%d'", cfg.code);
  }
}

inline double total_onset_with_t0(double onset,
                                  const AccDistParams& cfg) {
  if (!std::isfinite(onset)) {
    Rcpp::stop("Accumulator parameter 'onset' must be finite");
  }
  if (!std::isfinite(cfg.t0)) {
    Rcpp::stop("Accumulator parameter 't0' must be finite");
  }
  return onset + cfg.t0;
}

inline double acc_density_from_cfg(double t,
                                   double onset,
                                   double q,
                                   const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(q)) {
    Rcpp::stop("Accumulator parameter 'q' must be finite");
  }
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

inline double acc_density_success_from_cfg(double t,
                                           double onset,
                                           const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(t) || t < 0.0) return 0.0;
  if (t < effective_onset) return 0.0;
  double dens = eval_pdf_single(cfg, t - effective_onset);
  if (Rcpp::NumericVector::is_na(dens) || !std::isfinite(dens)) {
    return NA_REAL;
  }
  return dens;
}

inline double acc_survival_from_cfg(double t,
                                    double onset,
                                    double q,
                                    const AccDistParams& cfg) {
  double effective_onset = total_onset_with_t0(onset, cfg);
  if (!std::isfinite(q)) {
    Rcpp::stop("Accumulator parameter 'q' must be finite");
  }
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

inline Rcpp::List node_result_to_list(const NodeEvalResult& res) {
  return Rcpp::List::create(
    Rcpp::Named("density") = res.density,
    Rcpp::Named("survival") = res.survival,
    Rcpp::Named("cdf") = res.cdf
  );
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

inline Rcpp::List scenario_records_to_r(const std::vector<ScenarioRecord>& records) {
  Rcpp::List out(records.size());
  for (std::size_t i = 0; i < records.size(); ++i) {
    const ScenarioRecord& rec = records[i];
    out[i] = Rcpp::List::create(
      Rcpp::Named("weight") = rec.weight,
      Rcpp::Named("forced_complete") = Rcpp::IntegerVector(rec.forced_complete.begin(), rec.forced_complete.end()),
      Rcpp::Named("forced_survive") = Rcpp::IntegerVector(rec.forced_survive.begin(), rec.forced_survive.end())
    );
  }
  return out;
}

inline Rcpp::IntegerVector forced_set_to_vector(const std::unordered_set<int>& ids) {
  std::vector<int> sorted = set_to_sorted_vector(ids);
  return Rcpp::IntegerVector(sorted.begin(), sorted.end());
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

struct TrialRowColumns {
  R_xlen_t nrows{0};
  bool has_acc_id{false};
  Rcpp::CharacterVector acc_id_col;
  bool has_acc_column{false};
  bool acc_is_numeric{false};
  Rcpp::NumericVector acc_num_col;
  Rcpp::CharacterVector acc_chr_col;
  Rcpp::CharacterVector dist_col;
  Rcpp::NumericVector onset_col;
  Rcpp::NumericVector q_col;
  Rcpp::CharacterVector shared_col;
  Rcpp::List comps_col;
  Rcpp::List params_list_col;
  std::vector<std::string> param_cols;
  std::vector<SEXP> param_col_data;
  std::vector<int> param_col_types;
};

TrialRowColumns capture_trial_row_columns(const Rcpp::DataFrame& df) {
  TrialRowColumns cols;
  cols.nrows = df.nrows();
  cols.has_acc_id = df.containsElementNamed("accumulator_id");
  if (cols.has_acc_id) {
    cols.acc_id_col = Rcpp::CharacterVector(df["accumulator_id"]);
  }
  cols.has_acc_column = df.containsElementNamed("accumulator");
  if (cols.has_acc_column) {
    SEXP acc_col = df["accumulator"];
    if (Rf_isNumeric(acc_col)) {
      cols.acc_is_numeric = true;
      cols.acc_num_col = Rcpp::NumericVector(acc_col);
    } else {
      cols.acc_is_numeric = false;
      cols.acc_chr_col = Rcpp::CharacterVector(acc_col);
    }
  }
  cols.dist_col = df.containsElementNamed("dist")
    ? Rcpp::CharacterVector(df["dist"])
    : Rcpp::CharacterVector();
  cols.onset_col = df.containsElementNamed("onset")
    ? Rcpp::NumericVector(df["onset"])
    : Rcpp::NumericVector();
  cols.q_col = df.containsElementNamed("q")
    ? Rcpp::NumericVector(df["q"])
    : Rcpp::NumericVector();
  cols.shared_col = df.containsElementNamed("shared_trigger_id")
    ? Rcpp::CharacterVector(df["shared_trigger_id"])
    : Rcpp::CharacterVector();
  cols.comps_col = df.containsElementNamed("components")
    ? Rcpp::List(df["components"])
    : Rcpp::List();
  cols.params_list_col = df.containsElementNamed("params")
    ? Rcpp::List(df["params"])
    : Rcpp::List();

  static const std::unordered_set<std::string> kBaseCols = {
    "trial", "component", "accumulator", "accumulator_id",
    "accumulator_index", "acc_idx", "type", "role",
    "outcome", "rt", "params", "onset", "q",
    "condition", "component_weight", "shared_trigger_id",
    "dist", "components"
  };
  Rcpp::CharacterVector df_names = df.names();
  if (!df_names.isNULL()) {
    for (R_xlen_t i = 0; i < df_names.size(); ++i) {
      if (df_names[i] == NA_STRING) continue;
      std::string col_name = Rcpp::as<std::string>(df_names[i]);
      if (kBaseCols.find(col_name) != kBaseCols.end()) continue;
      SEXP column = df[col_name];
      cols.param_cols.push_back(col_name);
      cols.param_col_data.push_back(column);
      cols.param_col_types.push_back(TYPEOF(column));
    }
  }
  return cols;
}

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
  AccDistParams dist_cfg{};
  bool has_components{false};
  std::vector<std::string> components;
  std::string shared_trigger_id;
};

struct TrialParamSet {
  std::unordered_map<int, TrialAccumulatorParams> acc_params;
};

inline const TrialAccumulatorParams* get_trial_param_entry(const TrialParamSet* trial_params,
                                                          int acc_index) {
  if (!trial_params) return nullptr;
  auto it = trial_params->acc_params.find(acc_index);
  if (it == trial_params->acc_params.end()) return nullptr;
  return &it->second;
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
                                const TrialParamSet* trial_params = nullptr) {
  if (label.empty()) {
    Rcpp::stop("native_event_eval: empty label");
  }
  if (label == "__DEADLINE__" || label == "__GUESS__") {
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
    return make_node_result(density, survival, cdf);
  }

  auto pool_it = ctx.pool_index.find(label);
  if (pool_it != ctx.pool_index.end()) {
    const uuber::NativePool& pool = ctx.pools[pool_it->second];
    if (pool.members.empty()) {
      return make_node_result(0.0, 1.0, 0.0);
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
      return make_node_result(0.0, 1.0, 0.0);
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
        trial_params
      );
      density_vec[i] = child.density;
      survival_vec[i] = child.survival;
    }
    double density = pool_density_fast_cpp(density_vec, survival_vec, pool.k);
    double survival = pool_survival_fast_cpp(survival_vec, pool.k);
    double cdf = clamp_probability(1.0 - survival);
    if (!std::isfinite(density) || density < 0.0) density = 0.0;
    if (!std::isfinite(survival)) survival = 0.0;
    return make_node_result(density, survival, cdf);
  }

  Rcpp::stop("native_event_eval: unknown label '%s'", label);
}

struct ForcedKey {
  int complete_id{0};
  int survive_id{0};
};

struct TimeKey {
  double value{0.0};
};

struct CacheEntry {
  bool has_density{false};
  bool has_survival{false};
  bool has_cdf{false};
  bool scenarios_cached{false};
  NodeEvalResult result;
  std::vector<ScenarioRecord> scenarios;
};

struct TimeCache {
  std::unordered_map<std::string, CacheEntry> entries;
};

struct ParamCache {
  std::unordered_map<std::string, TimeCache> times;
};

struct ComponentCache {
  std::unordered_map<std::string, ParamCache> params;
};

struct NodeCache {
  std::unordered_map<std::string, ComponentCache> components;
};

struct EvalCache {
  std::unordered_map<int, NodeCache> nodes;
};

std::string component_cache_key(const std::string& component,
                                const std::string& trial_key);

std::string params_cache_key(const std::string& params_hash);

struct ComponentCacheEntry {
  std::string trial_type_key;
  std::string params_hash;
};

ComponentCacheEntry make_component_cache_entry(const std::string& component,
                                               const std::string& trial_key,
                                               const std::string& hash) {
  ComponentCacheEntry entry;
  entry.trial_type_key = component_cache_key(component, trial_key);
  entry.params_hash = params_cache_key(hash);
  return entry;
}

ComponentCacheEntry default_component_cache_entry(const std::string& component) {
  return make_component_cache_entry(component, std::string(), std::string());
}

using ComponentCacheMap = std::unordered_map<std::string, ComponentCacheEntry>;

ComponentCacheMap build_component_cache_map(const Rcpp::Nullable<Rcpp::List>& cache_opt) {
  ComponentCacheMap out;
  if (cache_opt.isNull()) return out;
  Rcpp::List cache_list(cache_opt.get());
  Rcpp::CharacterVector cache_names = cache_list.names();
  for (R_xlen_t i = 0; i < cache_list.size(); ++i) {
    Rcpp::RObject obj = cache_list[i];
    if (obj.isNULL()) continue;
    Rcpp::List entry(obj);
    std::string comp_label;
    if (!cache_names.isNULL() && cache_names.size() > i && cache_names[i] != NA_STRING) {
      comp_label = Rcpp::as<std::string>(cache_names[i]);
    }
    if (entry.containsElementNamed("component")) {
      Rcpp::RObject comp_obj = entry["component"];
      if (!comp_obj.isNULL()) {
        comp_label = Rcpp::as<std::string>(comp_obj);
      }
    }
    std::string trial_key;
    std::string params_hash;
    if (entry.containsElementNamed("key")) {
      Rcpp::RObject key_obj = entry["key"];
      if (!key_obj.isNULL()) {
        trial_key = Rcpp::as<std::string>(key_obj);
      }
    }
    if (entry.containsElementNamed("hash")) {
      Rcpp::RObject hash_obj = entry["hash"];
      if (!hash_obj.isNULL()) {
        params_hash = Rcpp::as<std::string>(hash_obj);
      }
    }
    ComponentCacheEntry cache_entry = make_component_cache_entry(comp_label, trial_key, params_hash);
    out[comp_label] = cache_entry;
  }
  return out;
}

ComponentCacheEntry lookup_component_cache_entry(const ComponentCacheMap& map,
                                                 const std::string& component) {
  auto it = map.find(component);
  if (it != map.end()) return it->second;
  return default_component_cache_entry(component);
}

std::vector<ComponentCacheEntry> build_component_cache_entries(
  const std::vector<std::string>& components,
  const ComponentCacheMap& cache_map) {
  std::vector<ComponentCacheEntry> entries;
  entries.reserve(components.size());
  for (const auto& comp : components) {
    if (cache_map.empty()) {
      entries.push_back(default_component_cache_entry(comp));
    } else {
      entries.push_back(lookup_component_cache_entry(cache_map, comp));
    }
  }
  return entries;
}

std::vector<std::string> character_vector_to_std(const Rcpp::CharacterVector& vec) {
  std::vector<std::string> out;
  out.reserve(vec.size());
  for (R_xlen_t i = 0; i < vec.size(); ++i) {
    Rcpp::String s(vec[i]);
    if (s == NA_STRING) {
      out.emplace_back("__default__");
    } else {
      out.emplace_back(s.get_cstring());
    }
  }
  return out;
}

std::string component_cache_key(const std::string& component,
                                const std::string& trial_key = std::string()) {
  if (!trial_key.empty()) return trial_key;
  if (component.empty()) return "__default__";
  return component;
}

std::string params_cache_key(const std::string& params_hash) {
  if (params_hash.empty()) return "__base__";
  return params_hash;
}

std::string time_cache_key(double t) {
  if (std::isnan(t)) return "NA";
  if (!std::isfinite(t)) {
    if (t > 0) return "Inf";
    if (t < 0) return "-Inf";
    return "NA";
  }
  std::ostringstream oss;
  oss << std::setprecision(15) << t;
  return oss.str();
}

std::string ids_to_string(const std::unordered_set<int>& ids) {
  if (ids.empty()) return ".";
  std::vector<int> vec(ids.begin(), ids.end());
  sort_unique(vec);
  std::ostringstream oss;
  for (std::size_t i = 0; i < vec.size(); ++i) {
    if (i > 0) oss << ",";
    oss << vec[i];
  }
  return oss.str();
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
                const std::string& params_hash_key = std::string())
    : ctx(ctx_),
      cache(cache_),
      t(time),
      component(component_label),
      forced_complete(forced_complete_ref),
      forced_survive(forced_survive_ref),
      trial_params(params_ptr),
      trial_type_key(component_cache_key(component_label, trial_key)),
      params_hash(params_cache_key(params_hash_key)) {}

  const uuber::NativeContext& ctx;
  EvalCache& cache;
  double t;
  std::string component;
  const std::unordered_set<int>& forced_complete;
  const std::unordered_set<int>& forced_survive;
  const TrialParamSet* trial_params;
  std::string trial_type_key;
  std::string params_hash;
};

std::string forced_key_string(const std::unordered_set<int>& forced_complete,
                              const std::unordered_set<int>& forced_survive) {
  return ids_to_string(forced_complete) + "|" + ids_to_string(forced_survive);
}

inline const TrialAccumulatorParams* get_state_params(const NodeEvalState& state,
                                                          int acc_index) {
  return get_trial_param_entry(state.trial_params, acc_index);
}

std::string member_ids_signature(const Rcpp::IntegerVector& ids) {
  if (ids.size() == 0) return ".";
  std::ostringstream oss;
  for (R_xlen_t i = 0; i < ids.size(); ++i) {
    if (i > 0) oss << ",";
    oss << static_cast<int>(ids[i]);
  }
  return oss.str();
}

std::string pool_template_cache_key(const std::string& pool_label,
                                    const std::string& component,
                                    int k,
                                    const Rcpp::IntegerVector& member_ids) {
  std::ostringstream oss;
  oss << pool_label << "|" << component_cache_key(component) << "|" << k << "|"
      << member_ids_signature(member_ids);
  return oss.str();
}

CacheEntry& fetch_cache_entry(NodeEvalState& state, int node_id) {
  NodeCache& node_cache = state.cache.nodes[node_id];
  ComponentCache& comp_cache = node_cache.components[state.trial_type_key];
  ParamCache& param_cache = comp_cache.params[state.params_hash];
  TimeCache& time_cache = param_cache.times[time_cache_key(state.t)];
  std::string forced_key = forced_key_string(state.forced_complete, state.forced_survive);
  return time_cache.entries[forced_key];
}

bool fetch_cached_result(NodeEvalState& state,
                         int node_id,
                         NodeEvalResult& out) {
  CacheEntry& entry = fetch_cache_entry(state, node_id);
  if (!entry.has_density && !entry.has_survival && !entry.has_cdf) {
    return false;
  }
  out = entry.result;
  return true;
}

void store_cached_result(NodeEvalState& state,
                         int node_id,
                         const NodeEvalResult& res) {
  CacheEntry& entry = fetch_cache_entry(state, node_id);
  entry.result = res;
  entry.has_density = true;
  entry.has_survival = true;
  entry.has_cdf = true;
}

bool fetch_cached_scenarios(NodeEvalState& state,
                            int node_id,
                            std::vector<ScenarioRecord>& out) {
  CacheEntry& entry = fetch_cache_entry(state, node_id);
  if (!entry.scenarios_cached) return false;
  out = entry.scenarios;
  return true;
}

void store_cached_scenarios(NodeEvalState& state,
                            int node_id,
                            const std::vector<ScenarioRecord>& records) {
  CacheEntry& entry = fetch_cache_entry(state, node_id);
  entry.scenarios = records;
  entry.scenarios_cached = true;
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
  std::string params_hash;
  std::unordered_set<int> forced_complete;
  std::unordered_set<int> forced_survive;
  const TrialParamSet* trial_params;
};

std::string guard_survival_cache_key(const GuardEvalInput& input, double t) {
  std::ostringstream oss;
  oss << input.node.id << "|"
      << component_cache_key(input.component, input.trial_type_key) << "|"
      << params_cache_key(input.params_hash) << "|"
      << time_cache_key(t) << "|" << forced_key_string(input.forced_complete, input.forced_survive);
  return oss.str();
}

const uuber::NativeNode& fetch_node(const uuber::NativeContext& ctx, int node_id);
NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                     EvalCache& cache,
                                     int node_id,
                                     double time,
                                     const std::string& component,
                                     const std::unordered_set<int>& forced_complete,
                                     const std::unordered_set<int>& forced_survive);
GuardEvalInput make_guard_input(const uuber::NativeContext& ctx,
                                const uuber::NativeNode& node,
                                EvalCache& cache,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive,
                                const std::string& trial_key,
                                const std::string& params_hash,
                                const TrialParamSet* trial_params);
double guard_effective_survival_internal(const GuardEvalInput& input,
                                         double t,
                                         const IntegrationSettings& settings);
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
  Rcpp::IntegerVector shared_index_vec(static_cast<R_xlen_t>(n));
  Rcpp::IntegerVector member_ids_vec(static_cast<R_xlen_t>(n));

  std::vector<std::string> shared_ids(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::string& member_label = active_members[i];
    member_ids_vec[static_cast<R_xlen_t>(i)] = label_id(state.ctx, member_label);
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

  std::unordered_map<std::string, int> shared_group_map;
  int next_shared_group = 1;
  for (std::size_t i = 0; i < n; ++i) {
    if (shared_ids[i].empty()) {
      shared_index_vec[static_cast<R_xlen_t>(i)] = NA_INTEGER;
    } else {
      auto it = shared_group_map.find(shared_ids[i]);
      if (it == shared_group_map.end()) {
        shared_group_map[shared_ids[i]] = next_shared_group;
        shared_index_vec[static_cast<R_xlen_t>(i)] = next_shared_group;
        ++next_shared_group;
      } else {
        shared_index_vec[static_cast<R_xlen_t>(i)] = it->second;
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
      state.trial_params
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
  std::string template_key = pool_template_cache_key(
    pool_label,
    component_cache_key(state.component),
    k,
    member_ids_vec
  );
  Rcpp::List template_info;
  auto tmpl_it = state.ctx.pool_template_cache.find(template_key);
  if (tmpl_it != state.ctx.pool_template_cache.end()) {
    template_info = tmpl_it->second;
  } else {
    template_info = pool_build_templates_cpp(
      static_cast<int>(n),
      member_ids_vec,
      pool_label_id,
      k
    );
    state.ctx.pool_template_cache[template_key] = template_info;
  }
  Rcpp::List templates = template_info["templates"];
  if (templates.size() == 0) {
    return {};
  }

  Rcpp::IntegerVector forced_complete_vec = forced_set_to_vector(state.forced_complete);
  Rcpp::IntegerVector forced_survive_vec = forced_set_to_vector(state.forced_survive);

  Rcpp::List combine_res = pool_density_combine_cpp(
    dens_vec,
    cdf_vec,
    surv_vec,
    cdf_success_vec,
    surv_success_vec,
    shared_index_vec,
    templates,
    forced_complete_vec,
    forced_survive_vec
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
    state.forced_survive
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
                                                  state.params_hash,
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
  std::vector<ScenarioRecord> cached;
  if (fetch_cached_scenarios(state, node_id, cached)) {
    return cached;
  }
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
  store_cached_scenarios(state, node_id, result);
  return result;
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state);

const uuber::NativeNode& fetch_node(const uuber::NativeContext& ctx, int node_id) {
  auto node_it = ctx.node_index.find(node_id);
  if (node_it == ctx.node_index.end()) {
    Rcpp::stop("native_node_eval: unknown node id %d", node_id);
  }
  return ctx.nodes[node_it->second];
}

NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                     EvalCache& cache,
                                     int node_id,
                                     double time,
                                     const std::string& component,
                                     const std::unordered_set<int>& forced_complete,
                                     const std::unordered_set<int>& forced_survive) {
  NodeEvalState local(ctx, cache, time, component, forced_complete, forced_survive);
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
                                const std::string& params_hash = std::string(),
                                const TrialParamSet* trial_params = nullptr) {
  GuardEvalInput input{
    ctx,
    node,
    cache,
    component,
    component_cache_key(component, trial_key),
    params_cache_key(params_hash),
    {},
    {},
    trial_params
  };
  input.forced_complete = filter_forced_scope(forced_complete, node.source_ids);
  input.forced_survive = filter_forced_scope(forced_survive, node.source_ids);
  return input;
}

double guard_effective_survival_internal(const GuardEvalInput& input,
                                         double t,
                                         const IntegrationSettings& settings);

double guard_density_internal(const GuardEvalInput& input,
                              double t,
                              const IntegrationSettings& settings);

double guard_cdf_internal(const GuardEvalInput& input,
                          double t,
                          const IntegrationSettings& settings);

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state) {
  NodeEvalResult cached;
  if (fetch_cached_result(state, node_id, cached)) {
    return cached;
  }
  const uuber::NativeNode& node = fetch_node(state.ctx, node_id);
  NodeEvalResult result;
  if (node.kind == "event") {
    if (node.source.empty()) {
      Rcpp::stop("native_node_eval: event node missing source id");
    }
    result = eval_event_label(
      state.ctx,
      node.source,
      state.t,
      state.component,
      state.forced_complete,
      state.forced_survive,
      state.trial_params
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
                                                  state.params_hash,
                                                  state.trial_params);
    IntegrationSettings settings;
    double density = guard_density_internal(guard_input, state.t, settings);
    double cdf = guard_cdf_internal(guard_input, state.t, settings);
    double survival = clamp_probability(1.0 - cdf);
    result = make_node_result(density, survival, cdf);
  } else {
    Rcpp::stop("native_node_eval: unsupported node kind '%s'", node.kind);
  }
  store_cached_result(state, node_id, result);
  return result;
}

double guard_effective_survival_internal(const GuardEvalInput& input,
                                         double t,
                                         const IntegrationSettings& settings) {
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
  std::string cache_key = guard_survival_cache_key(input, t);
  auto cached_surv = input.ctx.guard_survival_cache.find(cache_key);
  if (cached_surv != input.ctx.guard_survival_cache.end()) {
    return cached_surv->second;
  }
  if (input.node.unless_ids.empty()) {
    NodeEvalResult block = eval_node_with_forced(
      input.ctx,
      input.cache,
      blocker_id,
      t,
      input.component,
      input.forced_complete,
      input.forced_survive
    );
    double surv = clamp_probability(block.survival);
    input.ctx.guard_survival_cache[cache_key] = surv;
    return surv;
  }
  auto integrand = [&](double u) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    NodeEvalResult block = eval_node_with_forced(
      input.ctx,
      input.cache,
      blocker_id,
      u,
      input.component,
      input.forced_complete,
      input.forced_survive
    );
    double density = block.density;
    if (!std::isfinite(density) || density <= 0.0) return 0.0;
    double prod = 1.0;
    for (int unl_id : input.node.unless_ids) {
      if (unl_id < 0) continue;
      NodeEvalResult protector = eval_node_with_forced(
        input.ctx,
        input.cache,
        unl_id,
        u,
        input.component,
        input.forced_complete,
        input.forced_survive
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
  double integral = uuber::integrate_boost_fn(
    integrand,
    0.0,
    t,
    settings.rel_tol,
    settings.abs_tol,
    settings.max_depth
  );
  double surv = 1.0 - integral;
  surv = clamp_probability(surv);
  input.ctx.guard_survival_cache[cache_key] = surv;
  return surv;
}

double guard_reference_density(const GuardEvalInput& input,
                               double t) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  int reference_id = input.node.reference_id;
  if (reference_id < 0) {
    return 0.0;
  }
  NodeEvalResult ref = eval_node_with_forced(
    input.ctx,
    input.cache,
    reference_id,
    t,
    input.component,
    input.forced_complete,
    input.forced_survive
  );
  double dens_ref = ref.density;
  if (!std::isfinite(dens_ref) || dens_ref <= 0.0) return 0.0;
  return dens_ref;
}

double guard_density_internal(const GuardEvalInput& input,
                              double t,
                              const IntegrationSettings& settings) {
  double dens_ref = guard_reference_density(input, t);
  if (dens_ref <= 0.0) return 0.0;
  double s_eff = guard_effective_survival_internal(input, t, settings);
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
  auto integrand = [&](double u) -> double {
    return guard_density_internal(input, u, settings);
  };
  double val = uuber::integrate_boost_fn(
    integrand,
    0.0,
    t,
    settings.rel_tol,
    settings.abs_tol,
    settings.max_depth
  );
  if (!std::isfinite(val) || val < 0.0) val = 0.0;
  return clamp_probability(val);
}

IntegrationSettings resolve_integration_settings(double rel_tol,
                                                 double abs_tol,
                                                 int max_depth) {
  IntegrationSettings settings;
  if (std::isfinite(rel_tol) && rel_tol > 0.0) {
    settings.rel_tol = rel_tol;
  }
  if (std::isfinite(abs_tol) && abs_tol > 0.0) {
    settings.abs_tol = abs_tol;
  }
  if (max_depth > 0) {
    settings.max_depth = max_depth;
  }
  return settings;
}

GuardEvalInput resolve_guard_input(uuber::NativeContext& ctx,
                                   int node_id,
                                   EvalCache& cache,
                                   const std::string& component,
                                   const std::unordered_set<int>& forced_complete,
                                   const std::unordered_set<int>& forced_survive) {
  auto node_it = ctx.node_index.find(node_id);
  if (node_it == ctx.node_index.end()) {
    Rcpp::stop("native_guard_eval: unknown node id %d", node_id);
  }
  const uuber::NativeNode& node = ctx.nodes[node_it->second];
  if (node.kind != "guard") {
    Rcpp::stop("native_guard_eval: node id %d is not a guard expression", node_id);
  }
  return make_guard_input(ctx, node, cache, component, forced_complete, forced_survive);
}

std::vector<int> gather_blocker_sources(const uuber::NativeContext& ctx,
                                        const uuber::NativeNode& guard_node) {
  std::vector<int> sources;
  if (guard_node.blocker_id < 0) return sources;
  const uuber::NativeNode& blocker = fetch_node(ctx, guard_node.blocker_id);
  sources = blocker.source_ids;
  if (blocker.kind == "event" && !blocker.source.empty()) {
    int lid = label_id(ctx, blocker.source);
    if (lid != NA_INTEGER) {
      sources.push_back(lid);
    }
  }
  sort_unique(sources);
  return sources;
}

Rcpp::List native_guard_eval_impl(SEXP ctxSEXP,
                                   int guard_node_id,
                                   double t,
                                   Rcpp::Nullable<Rcpp::String> component,
                                   SEXP forced_complete,
                                   SEXP forced_survive,
                                   double rel_tol,
                                   double abs_tol,
                                   int max_depth) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_guard_eval: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> fc_set = make_forced_set(fc_vec);
  std::unordered_set<int> fs_set = make_forced_set(fs_vec);
  EvalCache guard_cache;
  GuardEvalInput guard_input = resolve_guard_input(*ctx,
                                                   guard_node_id,
                                                   guard_cache,
                                                   component_label,
                                                   fc_set,
                                                   fs_set);
  IntegrationSettings settings = resolve_integration_settings(rel_tol, abs_tol, max_depth);
  double density = guard_density_internal(guard_input, t, settings);
  double cdf = guard_cdf_internal(guard_input, t, settings);
  double survival = clamp_probability(1.0 - cdf);
  double eff_surv = guard_effective_survival_internal(guard_input, t, settings);
  NodeEvalResult result = make_node_result(density, survival, cdf);
  Rcpp::List out = node_result_to_list(result);
  out["effective_survival"] = eff_surv;
  return out;
}

double native_guard_effective_survival_impl(SEXP ctxSEXP,
                                            int guard_node_id,
                                            double t,
                                            Rcpp::Nullable<Rcpp::String> component,
                                            SEXP forced_complete,
                                            SEXP forced_survive,
                                            double rel_tol,
                                            double abs_tol,
                                            int max_depth) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_guard_effective_survival: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> fc_set = make_forced_set(fc_vec);
  std::unordered_set<int> fs_set = make_forced_set(fs_vec);
  EvalCache guard_cache2;
  GuardEvalInput guard_input = resolve_guard_input(*ctx,
                                                   guard_node_id,
                                                   guard_cache2,
                                                   component_label,
                                                   fc_set,
                                                   fs_set);
  IntegrationSettings settings = resolve_integration_settings(rel_tol, abs_tol, max_depth);
  return guard_effective_survival_internal(guard_input, t, settings);
}

struct CompetitorMeta {
  int node_id;
  const uuber::NativeNode* node;
  std::vector<int> sources;
  bool has_guard{false};
  bool scenario_sensitive{false};
};

bool node_contains_guard(int node_id,
                         const uuber::NativeContext& ctx,
                         std::unordered_map<int, bool>& memo) {
  auto it = memo.find(node_id);
  if (it != memo.end()) return it->second;
  const uuber::NativeNode& node = fetch_node(ctx, node_id);
  bool result = false;
  if (node.kind == "guard") {
    result = true;
  } else if (node.kind == "and" || node.kind == "or") {
    for (int child_id : node.args) {
      if (child_id < 0) continue;
      if (node_contains_guard(child_id, ctx, memo)) {
        result = true;
        break;
      }
    }
  } else if (node.kind == "not") {
    if (node.arg_id >= 0) {
      result = node_contains_guard(node.arg_id, ctx, memo);
    }
  } else if (node.kind == "guard") {
    result = true;
  }
  memo[node_id] = result;
  return result;
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

double evaluate_survival_with_forced(int node_id,
                                     const std::unordered_set<int>& forced_survive,
                                     const std::string& component,
                                     double t,
                                     const uuber::NativeContext& ctx,
                                     const std::string& trial_key = std::string(),
                                     const std::string& params_hash = std::string(),
                                     const TrialParamSet* trial_params = nullptr) {
  static const std::unordered_set<int> empty_forced_complete;
  EvalCache cache;
  NodeEvalState state(ctx,
                      cache,
                      t,
                      component,
                      empty_forced_complete,
                      forced_survive,
                      trial_params,
                      trial_key,
                      params_hash);
  NodeEvalResult res = eval_node_recursive(node_id, state);
  return clamp_probability(res.survival);
}

double compute_guard_free_cluster_value(const std::vector<int>& cluster_indices,
                                        const std::vector<CompetitorMeta>& metas,
                                        const uuber::NativeContext& ctx,
                                        double t,
                                        const std::string& component,
                                        const std::string& trial_key = std::string(),
                                        const std::string& params_hash = std::string(),
                                        const TrialParamSet* trial_params = nullptr) {
  std::unordered_set<int> empty_forced;
  double prod = 1.0;
  for (int idx : cluster_indices) {
    double surv = evaluate_survival_with_forced(metas[static_cast<std::size_t>(idx)].node_id,
                                                empty_forced,
                                                component,
                                                t,
                                                ctx,
                                                trial_key,
                                                params_hash,
                                                trial_params);
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
                                   const std::string& params_hash = std::string(),
                                   const TrialParamSet* trial_params = nullptr) {
  struct OrderingMeta {
    int index;
    std::size_t source_count;
    bool scenario_sensitive;
  };
  std::vector<OrderingMeta> order;
  order.reserve(cluster_indices.size());
  for (int idx : cluster_indices) {
    const CompetitorMeta& meta = metas[static_cast<std::size_t>(idx)];
    order.push_back({idx,
                     meta.sources.size(),
                     meta.scenario_sensitive});
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

  std::unordered_set<int> forced_survive;
  double prod = 1.0;
  for (const OrderingMeta& meta : order) {
    const CompetitorMeta& node_meta = metas[static_cast<std::size_t>(meta.index)];
    double surv = evaluate_survival_with_forced(node_meta.node_id,
                                                forced_survive,
                                                component,
                                                t,
                                                ctx,
                                                trial_key,
                                                params_hash,
                                                trial_params);
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
                                    const std::string& params_hash = std::string(),
                                    const TrialParamSet* trial_params = nullptr) {
  if (competitor_ids.empty()) return 1.0;
  std::vector<CompetitorMeta> metas;
  metas.reserve(competitor_ids.size());
  std::unordered_map<int, bool> guard_memo;
  for (int node_id : competitor_ids) {
    if (node_id == NA_INTEGER) {
      Rcpp::stop("native_competitor_survival: invalid node id");
    }
    const uuber::NativeNode& node = fetch_node(ctx, node_id);
    CompetitorMeta meta;
    meta.node_id = node_id;
    meta.node = &node;
    meta.sources = node.source_ids;
    sort_unique(meta.sources);
    meta.has_guard = node_contains_guard(node_id, ctx, guard_memo);
    meta.scenario_sensitive = node.scenario_sensitive;
    metas.push_back(std::move(meta));
  }

  std::vector<std::vector<int>> clusters = build_competitor_clusters(metas);
  if (clusters.empty()) return 1.0;

  double product = 1.0;
  for (const auto& cluster : clusters) {
    bool cluster_has_guard = false;
    for (int idx : cluster) {
      if (metas[static_cast<std::size_t>(idx)].has_guard) {
        cluster_has_guard = true;
        break;
      }
    }
    double cluster_val = cluster_has_guard
      ? compute_guard_cluster_value(cluster, metas, ctx, t, component_label, trial_type_key, params_hash, trial_params)
      : compute_guard_free_cluster_value(cluster, metas, ctx, t, component_label, trial_type_key, params_hash, trial_params);
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

double native_competitor_survival_impl(SEXP ctxSEXP,
                                       const Rcpp::IntegerVector& competitor_ids,
                                       double t,
                                       Rcpp::Nullable<Rcpp::String> component) {
  if (competitor_ids.size() == 0) return 1.0;
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_competitor_survival: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> ids = integer_vector_to_std(competitor_ids, false);
  return competitor_survival_internal(*ctx, ids, t, component_label);
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
  const std::string& params_hash = std::string()) {
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  EvalCache eval_cache;
  NodeEvalState state(ctx,
                      eval_cache,
                      t,
                      component_label,
                      forced_complete,
                      forced_survive,
                      trial_params,
                      trial_type_key,
                      params_hash);
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
                                               state.params_hash,
                                               trial_params);
    if (!std::isfinite(surv) || surv <= 0.0) {
      return 0.0;
    }
    density *= surv;
  }
  return density;
}

Rcpp::List native_density_with_competitors_impl(SEXP ctxSEXP,
                                                int node_id,
                                                double t,
                                                Rcpp::Nullable<Rcpp::String> component,
                                                SEXP forced_complete,
                                                SEXP forced_survive,
                                                const Rcpp::IntegerVector& competitor_ids,
                                                const TrialParamSet* trial_params) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_density_with_competitors: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  std::vector<int> comp_vec = integer_vector_to_std(competitor_ids, false);
  double density = node_density_with_competitors_internal(
    *ctx,
    node_id,
    t,
    component_label,
    forced_complete_set,
    forced_survive_set,
    comp_vec,
    trial_params
  );
  if (!std::isfinite(density) || density <= 0.0) {
    density = 0.0;
  }
  return Rcpp::List::create(Rcpp::Named("density") = density);
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
                                       const TrialParamSet* trial_params) {
  if (upper <= 0.0) {
    return 0.0;
  }
  auto integrand = [&](double u) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    Rcpp::List dens = native_density_with_competitors_impl(ctxSEXP,
                                                          node_id,
                                                          u,
                                                          component,
                                                          forced_complete,
                                                          forced_survive,
                                                          competitor_ids,
                                                          trial_params);
    double val = Rcpp::as<double>(dens["density"]);
    if (!std::isfinite(val) || val <= 0.0) return 0.0;
    return val;
  };
  double integral = 0.0;
  if (std::isfinite(upper)) {
    integral = uuber::integrate_boost_fn(
      integrand,
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
      double val = integrand(t);
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
}

Rcpp::List native_component_plan_impl(const Rcpp::List& structure,
                                      const Rcpp::Nullable<Rcpp::DataFrame>& trial_rows_opt,
                                      const Rcpp::Nullable<Rcpp::String>& forced_component_opt,
                                      const Rcpp::Nullable<Rcpp::DataFrame>& component_weights_opt) {
  if (!structure.containsElementNamed("components")) {
    Rcpp::stop("native_component_plan: structure is missing 'components'");
  }
  Rcpp::DataFrame comp_tbl(structure["components"]);
  if (!comp_tbl.containsElementNamed("component_id") ||
      !comp_tbl.containsElementNamed("weight")) {
    Rcpp::stop("native_component_plan: components table must have component_id and weight");
  }
  Rcpp::CharacterVector component_ids = comp_tbl["component_id"];
  Rcpp::NumericVector base_weights = comp_tbl["weight"];
  if (component_ids.size() == 0) {
    Rcpp::stop("native_component_plan: no components supplied");
  }

  std::unordered_map<std::string, double> base_weight_map;
  base_weight_map.reserve(component_ids.size());
  for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
    if (component_ids[i] == NA_STRING) continue;
    std::string id = Rcpp::as<std::string>(component_ids[i]);
    double weight = (i < base_weights.size()) ? static_cast<double>(base_weights[i]) : NA_REAL;
    base_weight_map[id] = weight;
  }

  std::unordered_set<std::string> component_set;
  for (auto& kv : base_weight_map) {
    component_set.insert(kv.first);
  }

  std::string forced_component;
  if (forced_component_opt.isNotNull()) {
    forced_component = Rcpp::as<std::string>(forced_component_opt.get());
  }

  std::unordered_set<std::string> trial_component_filter;
  double trial_id = NA_REAL;
  if (trial_rows_opt.isNotNull()) {
    Rcpp::DataFrame trial_df(trial_rows_opt.get());
    if (trial_df.containsElementNamed("trial")) {
      Rcpp::NumericVector trial_col(trial_df["trial"]);
      if (trial_col.size() > 0) {
        trial_id = trial_col[0];
      }
    }
    if (trial_df.containsElementNamed("component")) {
      Rcpp::CharacterVector trial_comp = trial_df["component"];
      for (R_xlen_t i = 0; i < trial_comp.size(); ++i) {
        if (trial_comp[i] == NA_STRING) continue;
        std::string comp = Rcpp::as<std::string>(trial_comp[i]);
        if (component_set.find(comp) != component_set.end()) {
          trial_component_filter.insert(comp);
        }
      }
    }
  }

  std::vector<std::string> selected_components;
  if (!forced_component.empty()) {
    if (component_set.find(forced_component) != component_set.end()) {
      selected_components.push_back(forced_component);
    }
  } else if (!trial_component_filter.empty()) {
    for (const auto& comp : component_set) {
      if (trial_component_filter.find(comp) != trial_component_filter.end()) {
        selected_components.push_back(comp);
      }
    }
  }
  if (selected_components.empty()) {
    selected_components.assign(component_set.begin(), component_set.end());
  }
  std::sort(selected_components.begin(), selected_components.end(), [&](const std::string& a, const std::string& b) {
    return a < b;
  });

  if (selected_components.empty()) {
    Rcpp::stop("native_component_plan: unable to determine available components");
  }

  std::vector<double> weights(selected_components.size(), NA_REAL);
  for (std::size_t i = 0; i < selected_components.size(); ++i) {
    auto it = base_weight_map.find(selected_components[i]);
    if (it != base_weight_map.end()) {
      weights[i] = it->second;
    }
  }

  if (component_weights_opt.isNotNull() && !std::isnan(trial_id)) {
    Rcpp::DataFrame cw_df(component_weights_opt.get());
    if (cw_df.containsElementNamed("trial") &&
        cw_df.containsElementNamed("component") &&
        cw_df.containsElementNamed("weight")) {
      Rcpp::NumericVector cw_trial(cw_df["trial"]);
      Rcpp::CharacterVector cw_component(cw_df["component"]);
      Rcpp::NumericVector cw_weight(cw_df["weight"]);
      for (R_xlen_t i = 0; i < cw_trial.size(); ++i) {
        if (cw_component[i] == NA_STRING) continue;
        double trial_val = cw_trial[i];
        if (std::isnan(trial_val) || std::fabs(trial_val - trial_id) > 1e-9) continue;
        std::string comp = Rcpp::as<std::string>(cw_component[i]);
        auto it = std::find(selected_components.begin(), selected_components.end(), comp);
        if (it == selected_components.end()) continue;
        std::size_t idx = static_cast<std::size_t>(std::distance(selected_components.begin(), it));
        double w = cw_weight[i];
        if (std::isfinite(w) && w >= 0.0) {
          weights[idx] = w;
        }
      }
    }
  }

  bool need_uniform = weights.size() == 0;
  if (!need_uniform) {
    bool has_finite = false;
    bool invalid = false;
    for (double w : weights) {
      if (std::isfinite(w) && w >= 0.0) {
        has_finite = true;
      } else if (!std::isnan(w)) {
        invalid = true;
        break;
      }
    }
    if (!has_finite || invalid) {
      need_uniform = true;
    } else {
      double total = 0.0;
      for (double w : weights) {
        total += std::isfinite(w) ? w : 0.0;
      }
      if (!std::isfinite(total) || total <= 0.0) {
        need_uniform = true;
      } else {
        for (double& w : weights) {
          w = std::isfinite(w) ? w / total : 0.0;
        }
      }
    }
  }
  if (need_uniform) {
    double uniform = 1.0 / static_cast<double>(selected_components.size());
    std::fill(weights.begin(), weights.end(), uniform);
  }

  Rcpp::CharacterVector comp_out(selected_components.begin(), selected_components.end());
  Rcpp::NumericVector weight_out(weights.begin(), weights.end());
  return Rcpp::List::create(
    Rcpp::Named("components") = comp_out,
    Rcpp::Named("weights") = weight_out
  );
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

  bool has_acc_id = rows.containsElementNamed("accumulator_id");
  bool has_acc = rows.containsElementNamed("accumulator");
  if (!has_acc_id && !has_acc) {
    Rcpp::stop("trial parameter rows must include 'accumulator_id' or 'accumulator'");
  }
  Rcpp::CharacterVector acc_id_col;
  if (has_acc_id) {
    acc_id_col = Rcpp::CharacterVector(rows["accumulator_id"]);
  }
  Rcpp::RObject acc_col_obj;
  if (has_acc && !has_acc_id) {
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
    "trial", "component", "accumulator", "accumulator_id",
    "accumulator_index", "acc_idx", "type", "role",
    "outcome", "rt", "params", "onset", "q",
    "condition", "component_weight", "shared_trigger_id",
    "dist", "components", "params_hash"
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
    if (has_acc_id) {
      if (acc_id_col[i] == NA_STRING) continue;
      std::string acc_label = Rcpp::as<std::string>(acc_id_col[i]);
      auto it = ctx.accumulator_index.find(acc_label);
      if (it == ctx.accumulator_index.end()) {
        Rcpp::stop("trial parameter rows reference unknown accumulator '%s'", acc_label);
      }
      acc_idx = it->second;
    } else if (has_acc) {
      if (Rf_isNumeric(acc_col_obj)) {
        Rcpp::NumericVector acc_num(acc_col_obj);
        if (i >= acc_num.size()) continue;
        double raw_val = acc_num[i];
        if (Rcpp::NumericVector::is_na(raw_val)) continue;
        int idx_val = static_cast<int>(std::llround(raw_val)) - 1;
        if (idx_val < 0 || idx_val >= static_cast<int>(ctx.accumulators.size())) {
          Rcpp::stop("trial parameter rows reference invalid accumulator index '%d'",
                     idx_val + 1);
        }
        acc_idx = idx_val;
      } else {
        Rcpp::CharacterVector acc_chr(acc_col_obj);
        if (acc_chr[i] == NA_STRING) continue;
        std::string acc_label = Rcpp::as<std::string>(acc_chr[i]);
        auto it = ctx.accumulator_index.find(acc_label);
        if (it == ctx.accumulator_index.end()) {
          Rcpp::stop("trial parameter rows reference unknown accumulator '%s'", acc_label);
        }
        acc_idx = it->second;
      }
    }
    if (acc_idx < 0 || acc_idx >= static_cast<int>(ctx.accumulators.size())) continue;
    const uuber::NativeAccumulator& base = ctx.accumulators[acc_idx];

    TrialAccumulatorParams override;
    override.onset = (i < onset_col.size() && !Rcpp::NumericVector::is_na(onset_col[i]))
      ? static_cast<double>(onset_col[i])
      : base.onset;
    override.q = (i < q_col.size() && !Rcpp::NumericVector::is_na(q_col[i]))
      ? static_cast<double>(q_col[i])
      : base.q;
    if (i < shared_col.size() && shared_col[i] != NA_STRING) {
      override.shared_trigger_id = Rcpp::as<std::string>(shared_col[i]);
    } else {
      override.shared_trigger_id = base.shared_trigger_id;
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
        Rcpp::stop("trial parameter rows for '%s' missing distribution parameters", base.id);
      }
      override.dist_cfg = resolve_acc_params_entries(dist_name, param_entries);
    }

    params_set->acc_params[acc_idx] = std::move(override);
  }

  if (params_set->acc_params.empty()) {
    return nullptr;
  }
  return params_set;
}


Rcpp::List native_guard_scenarios_impl(SEXP ctxSEXP,
                                       int guard_node_id,
                                       double t,
                                       const Rcpp::List& reference_scenarios,
                                       Rcpp::Nullable<Rcpp::String> component,
                                       SEXP forced_complete,
                                       SEXP forced_survive,
                                       double rel_tol,
                                       double abs_tol,
                                       int max_depth) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_guard_scenarios: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  auto guard_it = ctx->node_index.find(guard_node_id);
  if (guard_it == ctx->node_index.end()) {
    Rcpp::stop("native_guard_scenarios: unknown guard node id %d", guard_node_id);
  }
  const uuber::NativeNode& guard_node = ctx->nodes[guard_it->second];
  if (guard_node.kind != "guard") {
    Rcpp::stop("native_guard_scenarios: node %d is not a guard", guard_node_id);
  }
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> base_fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> base_fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> base_fc_set = make_forced_set(base_fc_vec);
  std::unordered_set<int> base_fs_set = make_forced_set(base_fs_vec);
  IntegrationSettings settings = resolve_integration_settings(rel_tol, abs_tol, max_depth);
  std::vector<int> blocker_sources = gather_blocker_sources(*ctx, guard_node);
  std::vector<ScenarioRecord> records;
  records.reserve(static_cast<std::size_t>(reference_scenarios.size()));
  for (R_xlen_t i = 0; i < reference_scenarios.size(); ++i) {
    Rcpp::List sc(reference_scenarios[i]);
    if (sc.isNULL()) continue;
    double weight = Rcpp::as<double>(sc["weight"]);
    if (!std::isfinite(weight) || weight <= 0.0) continue;
    std::vector<int> sc_fc_vec = forced_vec_from_sexp(sc["forced_complete"]);
    std::vector<int> sc_fs_vec = forced_vec_from_sexp(sc["forced_survive"]);
    std::unordered_set<int> sc_fc_set = make_forced_set(sc_fc_vec);
    std::unordered_set<int> sc_fs_set = make_forced_set(sc_fs_vec);
    sc_fc_set.insert(base_fc_set.begin(), base_fc_set.end());
    sc_fs_set.insert(base_fs_set.begin(), base_fs_set.end());
    EvalCache scenario_cache;
    GuardEvalInput scenario_input = make_guard_input(*ctx,
                                                     guard_node,
                                                     scenario_cache,
                                                     component_label,
                                                     sc_fc_set,
                                                     sc_fs_set);
    double s_eff = guard_effective_survival_internal(scenario_input, t, settings);
    if (!std::isfinite(s_eff) || s_eff <= 0.0) continue;
    double final_weight = weight * s_eff;
    if (!std::isfinite(final_weight) || final_weight <= 0.0) continue;
    ScenarioRecord record;
    record.weight = final_weight;
    record.forced_complete = set_to_sorted_vector(scenario_input.forced_complete);
    record.forced_survive = union_vectors(set_to_sorted_vector(scenario_input.forced_survive), blocker_sources);
    records.push_back(std::move(record));
  }
  Rcpp::List out(records.size());
  for (std::size_t i = 0; i < records.size(); ++i) {
    const ScenarioRecord& rec = records[i];
    out[i] = Rcpp::List::create(
      Rcpp::Named("weight") = rec.weight,
      Rcpp::Named("forced_complete") = Rcpp::IntegerVector(rec.forced_complete.begin(), rec.forced_complete.end()),
      Rcpp::Named("forced_survive") = Rcpp::IntegerVector(rec.forced_survive.begin(), rec.forced_survive.end())
    );
  }
  return out;
}


Rcpp::List native_node_eval_impl(SEXP ctxSEXP,
                                 int node_id,
                                 double t,
                                 Rcpp::Nullable<Rcpp::String> component,
                                 SEXP forced_complete,
                                 SEXP forced_survive) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_node_eval: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  EvalCache eval_cache;
  NodeEvalState state(*ctx,
                      eval_cache,
                      t,
                      component_label,
                      forced_complete_set,
                      forced_survive_set);
  NodeEvalResult result = eval_node_recursive(node_id, state);
  return node_result_to_list(result);
}

Rcpp::List native_node_scenarios_impl(SEXP ctxSEXP,
                                      int node_id,
                                      double t,
                                      Rcpp::Nullable<Rcpp::String> component,
                                      SEXP forced_complete,
                                      SEXP forced_survive) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_node_scenarios: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  EvalCache eval_cache;
  NodeEvalState state(*ctx,
                      eval_cache,
                      t,
                      component_label,
                      forced_complete_set,
                      forced_survive_set);
  std::vector<ScenarioRecord> records = compute_node_scenarios(node_id, state);
  return scenario_records_to_r(records);
}

Rcpp::List native_likelihood_driver_impl(SEXP ctxSEXP,
                                         Rcpp::IntegerVector node_ids,
                                         Rcpp::NumericVector times,
                                         Rcpp::Nullable<Rcpp::String> component,
                                         SEXP forced_complete,
                                         SEXP forced_survive) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_likelihood_driver: context pointer is null");
  }
  if (node_ids.size() == 0 || times.size() == 0) {
    return Rcpp::List();
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  Rcpp::List node_results(node_ids.size());
  for (R_xlen_t i = 0; i < node_ids.size(); ++i) {
    int node_id = node_ids[i];
    if (node_id == NA_INTEGER) {
      node_results[i] = R_NilValue;
      continue;
    }
    Rcpp::List evals(times.size());
    for (R_xlen_t j = 0; j < times.size(); ++j) {
      double t = static_cast<double>(times[j]);
      EvalCache eval_cache;
      NodeEvalState state(*ctx,
                          eval_cache,
                          t,
                          component_label,
                          forced_complete_set,
                          forced_survive_set);
      NodeEvalResult result = eval_node_recursive(node_id, state);
      evals[j] = Rcpp::List::create(
        Rcpp::Named("time") = t,
        Rcpp::Named("density") = result.density,
        Rcpp::Named("survival") = result.survival,
        Rcpp::Named("cdf") = result.cdf
      );
    }
    node_results[i] = Rcpp::List::create(
      Rcpp::Named("node_id") = node_id,
      Rcpp::Named("evaluations") = evals
    );
  }
  return node_results;
}

Rcpp::List native_likelihood_eval_impl(SEXP ctxSEXP,
                                       const Rcpp::List& task_list) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_likelihood_eval: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  Rcpp::List out(task_list.size());
  for (R_xlen_t i = 0; i < task_list.size(); ++i) {
    Rcpp::List task(task_list[i]);
    if (task.isNULL() || task.size() == 0) {
      out[i] = R_NilValue;
      continue;
    }
    int node_id = task.containsElementNamed("node_id")
      ? Rcpp::as<int>(task["node_id"])
      : NA_INTEGER;
    if (node_id == NA_INTEGER) {
      out[i] = R_NilValue;
      continue;
    }
    Rcpp::NumericVector times = task.containsElementNamed("times")
      ? Rcpp::NumericVector(task["times"])
      : Rcpp::NumericVector();
    if (times.size() == 0) {
      out[i] = Rcpp::List::create(
        Rcpp::Named("node_id") = node_id,
        Rcpp::Named("evaluations") = Rcpp::List()
      );
      continue;
    }
    std::string component_label;
    if (task.containsElementNamed("component") && !Rf_isNull(task["component"])) {
      component_label = Rcpp::as<std::string>(task["component"]);
    }
    std::vector<int> fc_vec = forced_vec_from_sexp(task.containsElementNamed("forced_complete")
                                                   ? task["forced_complete"] : R_NilValue);
    std::vector<int> fs_vec = forced_vec_from_sexp(task.containsElementNamed("forced_survive")
                                                   ? task["forced_survive"] : R_NilValue);
    std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
    std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
    EvalCache eval_cache;
    Rcpp::List eval_list(times.size());
    for (R_xlen_t j = 0; j < times.size(); ++j) {
      double t = static_cast<double>(times[j]);
      NodeEvalState state(*ctx,
                          eval_cache,
                          t,
                          component_label,
                          forced_complete_set,
                          forced_survive_set);
      NodeEvalResult result = eval_node_recursive(node_id, state);
      eval_list[j] = Rcpp::List::create(
        Rcpp::Named("time") = t,
        Rcpp::Named("density") = result.density,
        Rcpp::Named("survival") = result.survival,
        Rcpp::Named("cdf") = result.cdf
      );
    }
    out[i] = Rcpp::List::create(
      Rcpp::Named("node_id") = node_id,
      Rcpp::Named("component") = component_label,
      Rcpp::Named("evaluations") = eval_list
    );
  }
  return out;
}

} // namespace

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

struct SharedGateSpec {
  std::string x_label;
  std::string y_label;
  std::string c_label;
};

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
                                     const std::vector<ComponentCacheEntry>& cache_entries) {
  if (!std::isfinite(t) || t < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::unordered_set<int> forced_empty;
  double total = 0.0;
  for (std::size_t i = 0; i < component_ids.size(); ++i) {
    ComponentCacheEntry cache_entry;
    if (i < cache_entries.size()) {
      cache_entry = cache_entries[i];
    } else {
      cache_entry = default_component_cache_entry(component_ids[i]);
    }
    double contrib = node_density_with_competitors_internal(
      ctx,
      node_id,
      t,
      component_ids[i],
      forced_empty,
      forced_empty,
      competitor_ids,
      params_ptr,
      cache_entry.trial_type_key,
      cache_entry.params_hash
    );
    if (!std::isfinite(contrib) || contrib < 0.0) contrib = 0.0;
    total += weights[i] * contrib;
  }
  if (!std::isfinite(total) || total <= 0.0) {
    return 0.0;
  }
  return total;
}

double native_trial_mixture_driver(SEXP ctxSEXP,
                                   int node_id,
                                   double t,
                                   Rcpp::CharacterVector component_ids,
                                   Rcpp::NumericVector weights,
                                   Rcpp::Nullable<Rcpp::String> forced_component,
                                   Rcpp::IntegerVector competitor_ids,
                                   SEXP trial_params);

std::unordered_map<std::string, double> build_component_deadlines(const Rcpp::List& structure,
                                                                  double default_deadline) {
  std::unordered_map<std::string, double> out;
  if (!structure.containsElementNamed("components")) {
    return out;
  }
  Rcpp::DataFrame comp_df(structure["components"]);
  if (!comp_df.containsElementNamed("component_id")) {
    return out;
  }
  Rcpp::CharacterVector comp_ids(comp_df["component_id"]);
  Rcpp::List attrs = comp_df.containsElementNamed("attrs")
    ? Rcpp::List(comp_df["attrs"])
    : Rcpp::List(comp_ids.size());
  for (R_xlen_t i = 0; i < comp_ids.size(); ++i) {
    std::string comp_label = "__default__";
    if (comp_ids[i] != NA_STRING) {
      comp_label = Rcpp::as<std::string>(comp_ids[i]);
    }
    double deadline = default_deadline;
    if (i < attrs.size()) {
      Rcpp::RObject attr_obj = attrs[i];
      if (!attr_obj.isNULL()) {
        Rcpp::List attr_list(attr_obj);
        if (attr_list.containsElementNamed("deadline")) {
          Rcpp::RObject val = attr_list["deadline"];
          if (!Rf_isNull(val)) {
            double attr_deadline = Rcpp::as<double>(val);
            if (std::isfinite(attr_deadline)) {
              deadline = attr_deadline;
            }
          }
        }
      }
    }
    out[comp_label] = deadline;
  }
  return out;
}

double lookup_component_deadline(const std::unordered_map<std::string, double>& map,
                                 const std::string& component,
                                 double default_deadline) {
  auto it = map.find(component);
  if (it == map.end()) return default_deadline;
  return it->second;
}

double evaluate_alias_mix(SEXP ctxSEXP,
                          const Rcpp::CharacterVector& components,
                          const Rcpp::NumericVector& weights,
                          const Rcpp::List& alias_sources,
                          double t,
                          Rcpp::Nullable<Rcpp::String> forced_component,
                          const TrialParamSet* params_ptr,
                          const ComponentCacheMap& cache_map) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::vector<std::string> comp_labels;
  std::vector<double> mix_weights;
  if (!build_component_mix(components, weights, forced_component, comp_labels, mix_weights)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::vector<ComponentCacheEntry> cache_entries = build_component_cache_entries(comp_labels, cache_map);
  double mix = 0.0;
  for (R_xlen_t i = 0; i < alias_sources.size(); ++i) {
    Rcpp::List src(alias_sources[i]);
    if (!src.containsElementNamed("node_id")) continue;
    int node_id = Rcpp::as<int>(src["node_id"]);
    Rcpp::IntegerVector comp_ids = src.containsElementNamed("competitor_ids")
      ? Rcpp::IntegerVector(src["competitor_ids"])
      : Rcpp::IntegerVector();
    double contrib = native_trial_mixture_internal(*ctx,
                                                   node_id,
                                                   t,
                                                   comp_labels,
                                                   mix_weights,
                                                   forced_component,
                                                   integer_vector_to_std(comp_ids, false),
                                                   params_ptr,
                                                   cache_entries);
    mix += contrib;
  }
  return mix;
}

double evaluate_na_map_mix(SEXP ctxSEXP,
                           const Rcpp::CharacterVector& components,
                           const Rcpp::NumericVector& weights,
                           const Rcpp::List& na_sources,
                           double rt,
                           const std::unordered_map<std::string, double>& component_deadlines,
                           double default_deadline,
                           double rel_tol,
                           double abs_tol,
                           int max_depth,
                           Rcpp::Nullable<Rcpp::String> forced_component,
                           const TrialParamSet* params_ptr,
                           const ComponentCacheMap& cache_map) {
  if (std::isfinite(rt) && rt >= 0.0) {
    return evaluate_alias_mix(ctxSEXP,
                              components,
                              weights,
                              na_sources,
                              rt,
                              forced_component,
                              params_ptr,
                              cache_map);
  }
  double mix = 0.0;
  for (R_xlen_t i = 0; i < components.size(); ++i) {
    SEXP comp_sexp = components[i];
    bool comp_is_na = comp_sexp == NA_STRING;
    Rcpp::String comp_val(comp_sexp);
    std::string comp_label = "__default__";
    if (!comp_is_na) {
      comp_label = static_cast<std::string>(comp_val);
    }
    double deadline = lookup_component_deadline(component_deadlines, comp_label, default_deadline);
    if (deadline <= 0.0) continue;
    Rcpp::Nullable<Rcpp::String> component_str;
    if (!comp_is_na) {
      component_str = Rcpp::Nullable<Rcpp::String>(comp_sexp);
    }
    double component_sum = 0.0;
    for (R_xlen_t s = 0; s < na_sources.size(); ++s) {
      Rcpp::List src(na_sources[s]);
      if (!src.containsElementNamed("node_id")) continue;
      int node_id = Rcpp::as<int>(src["node_id"]);
      Rcpp::IntegerVector comp_ids = src.containsElementNamed("competitor_ids")
        ? Rcpp::IntegerVector(src["competitor_ids"])
        : Rcpp::IntegerVector();
      double prob = native_outcome_probability_impl(ctxSEXP,
                                                    node_id,
                                                    deadline,
                                                    component_str,
                                                    R_NilValue,
                                                    R_NilValue,
                                                    comp_ids,
                                                    rel_tol,
                                                    abs_tol,
                                                    max_depth,
                                                    params_ptr);
      if (std::isfinite(prob) && prob > 0.0) {
        component_sum += prob;
      }
    }
    double weight = (i < weights.size()) ? static_cast<double>(weights[i]) : 0.0;
    mix += weight * component_sum;
  }
  return mix;
}

bool extract_shared_gate_spec(const Rcpp::List& entry, SharedGateSpec& spec) {
  if (!entry.containsElementNamed("shared_gate")) return false;
  Rcpp::RObject obj = entry["shared_gate"];
  if (obj.isNULL()) return false;
  Rcpp::List sg(obj);
  spec = SharedGateSpec();
  auto fetch_label = [&](const char* key, std::string& out) -> bool {
    if (!sg.containsElementNamed(key)) return false;
    Rcpp::RObject val = sg[key];
    if (Rf_isNull(val)) return false;
    Rcpp::CharacterVector cv(val);
    if (cv.size() == 0 || cv[0] == NA_STRING) return false;
    out = Rcpp::as<std::string>(cv[0]);
    return true;
  };
  if (!fetch_label("x_label", spec.x_label)) return false;
  if (!fetch_label("y_label", spec.y_label)) return false;
  if (!fetch_label("c_label", spec.c_label)) return false;
  return !(spec.x_label.empty() || spec.y_label.empty() || spec.c_label.empty());
}

const std::unordered_set<int> kEmptyForcedSet{};

double compute_shared_gate_palloc(const uuber::NativeContext& ctx,
                                  const SharedGateSpec& spec,
                                  double limit,
                                  double denom,
                                  const std::string& component,
                                  const TrialParamSet* trial_params,
                                  double rel_tol,
                                  double abs_tol,
                                  int max_depth) {
  if (!std::isfinite(limit) || limit <= 0.0) return 0.0;
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  auto integrand = [&](double u) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    NodeEvalResult x_res = eval_event_label(ctx,
                                            spec.x_label,
                                            u,
                                            component,
                                            kEmptyForcedSet,
                                            kEmptyForcedSet,
                                            trial_params);
    NodeEvalResult y_res = eval_event_label(ctx,
                                            spec.y_label,
                                            u,
                                            component,
                                            kEmptyForcedSet,
                                            kEmptyForcedSet,
                                            trial_params);
    double fX = safe_density(x_res.density);
    double FY = clamp_probability(1.0 - clamp_unit(y_res.survival));
    double val = fX * FY;
    if (!std::isfinite(val) || val <= 0.0) return 0.0;
    return val;
  };
  double integral = uuber::integrate_boost_fn(
    integrand,
    0.0,
    limit,
    rel_tol,
    abs_tol,
    max_depth
  );
  if (!std::isfinite(integral) || integral <= 0.0) {
    return 0.0;
  }
  double out = 1.0 - (integral / denom);
  if (!std::isfinite(out)) return 0.0;
  if (out < 0.0) out = 0.0;
  if (out > 1.0) out = 1.0;
  return out;
}

double shared_gate_density_component(const uuber::NativeContext& ctx,
                                     const SharedGateSpec& spec,
                                     double t,
                                     const std::string& component,
                                     const TrialParamSet* trial_params,
                                     double rel_tol,
                                     double abs_tol,
                                     int max_depth) {
  if (spec.x_label.empty() || spec.y_label.empty() || spec.c_label.empty()) {
    return 0.0;
  }
  if (!std::isfinite(t) || t < 0.0) {
    return 0.0;
  }
  NodeEvalResult x_res = eval_event_label(ctx,
                                          spec.x_label,
                                          t,
                                          component,
                                          kEmptyForcedSet,
                                          kEmptyForcedSet,
                                          trial_params);
  NodeEvalResult y_res = eval_event_label(ctx,
                                          spec.y_label,
                                          t,
                                          component,
                                          kEmptyForcedSet,
                                          kEmptyForcedSet,
                                          trial_params);
  NodeEvalResult c_res = eval_event_label(ctx,
                                          spec.c_label,
                                          t,
                                          component,
                                          kEmptyForcedSet,
                                          kEmptyForcedSet,
                                          trial_params);
  double fX = safe_density(x_res.density);
  double fC = safe_density(c_res.density);
  double FX = clamp_probability(1.0 - clamp_unit(x_res.survival));
  double FY = clamp_probability(1.0 - clamp_unit(y_res.survival));
  double FC = clamp_probability(1.0 - clamp_unit(c_res.survival));
  double SY = clamp_unit(y_res.survival);
  double term1 = fX * FC * SY;
  double termCother = fC * FX * SY;
  double term2 = 0.0;
  double denom = FX * FY;
  if (denom > 0.0 && fC > 0.0) {
    double palloc = compute_shared_gate_palloc(ctx,
                                               spec,
                                               t,
                                               denom,
                                               component,
                                               trial_params,
                                               rel_tol,
                                               abs_tol,
                                               max_depth);
    if (palloc > 0.0) {
      term2 = fC * FX * FY * palloc;
    }
  }
  double density = term1 + termCother + term2;
  if (!std::isfinite(density) || density <= 0.0) {
    return 0.0;
  }
  return density;
}

double evaluate_shared_gate_mix(const uuber::NativeContext& ctx,
                                const SharedGateSpec& spec,
                                double t,
                                const Rcpp::CharacterVector& component_ids,
                                const Rcpp::NumericVector& weights,
                                Rcpp::Nullable<Rcpp::String> forced_component,
                                const TrialParamSet* trial_params,
                                double rel_tol,
                                double abs_tol,
                                int max_depth,
                                const ComponentCacheMap& cache_map) {
  (void) cache_map;
  if (!std::isfinite(t) || t < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::vector<std::string> components;
  std::vector<double> mix_weights;
  if (!build_component_mix(component_ids, weights, forced_component, components, mix_weights)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double total = 0.0;
  for (std::size_t i = 0; i < components.size(); ++i) {
    double contrib = shared_gate_density_component(ctx,
                                                   spec,
                                                   t,
                                                   components[i],
                                                   trial_params,
                                                   rel_tol,
                                                   abs_tol,
                                                   max_depth);
    if (!std::isfinite(contrib) || contrib <= 0.0) continue;
    total += mix_weights[i] * contrib;
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
                                   Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_trial_mixture: context pointer is null");
  }
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
  ComponentCacheMap cache_map;
  std::vector<ComponentCacheEntry> cache_entries = build_component_cache_entries(components, cache_map);
  return native_trial_mixture_internal(*ctx,
                                       node_id,
                                       t,
                                       components,
                                       mix_weights,
                                       forced_component,
                                       comp_ids,
                                       params_ptr,
                                       cache_entries);
}

// [[Rcpp::export]]
Rcpp::List native_loglik_from_params_cpp(SEXP ctxSEXP,
                                         Rcpp::List structure,
                                         Rcpp::List trial_entries,
                                         Rcpp::Nullable<Rcpp::DataFrame> component_weights_opt,
                                         double default_deadline,
                                         double rel_tol,
                                         double abs_tol,
                                         int max_depth) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_loglik_from_params_cpp: context pointer is null");
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  R_xlen_t n_trials = trial_entries.size();
  Rcpp::NumericVector per_trial(n_trials);
  double total_loglik = 0.0;
  bool hit_neg_inf = false;
  std::unordered_map<std::string, double> component_deadlines =
    build_component_deadlines(structure, default_deadline);

  std::vector<std::unique_ptr<TrialParamSet>> params_cache;
  params_cache.reserve(static_cast<std::size_t>(n_trials));

  for (R_xlen_t i = 0; i < n_trials; ++i) {
    Rcpp::List entry(trial_entries[i]);
    std::string entry_type = entry.containsElementNamed("type")
      ? Rcpp::as<std::string>(entry["type"])
      : "direct";
    double t = entry.containsElementNamed("rt")
      ? Rcpp::as<double>(entry["rt"])
      : NA_REAL;
    SEXP trial_rows_sexp = entry.containsElementNamed("trial_rows")
      ? entry["trial_rows"]
      : R_NilValue;
    Rcpp::Nullable<Rcpp::DataFrame> trial_rows(trial_rows_sexp);
    SEXP forced_component_sexp = entry.containsElementNamed("forced_component")
      ? entry["forced_component"]
      : R_NilValue;
    Rcpp::Nullable<Rcpp::String> forced_component(forced_component_sexp);
    Rcpp::IntegerVector competitor_ids = entry.containsElementNamed("competitor_ids")
      ? Rcpp::IntegerVector(entry["competitor_ids"])
      : Rcpp::IntegerVector();
    const TrialParamSet* params_ptr = nullptr;
    if (!trial_rows.isNull()) {
      std::unique_ptr<TrialParamSet> holder = build_trial_params_from_df(*ctx, trial_rows);
      if (holder) {
        params_ptr = holder.get();
        params_cache.emplace_back(std::move(holder));
      }
    }
    ComponentCacheMap component_cache_map;
    if (entry.containsElementNamed("component_cache")) {
      Rcpp::RObject cache_obj = entry["component_cache"];
      if (!cache_obj.isNULL()) {
        component_cache_map = build_component_cache_map(Rcpp::Nullable<Rcpp::List>(cache_obj));
      }
    }

    Rcpp::List plan = native_component_plan_impl(structure,
                                                 trial_rows,
                                                 forced_component,
                                                 component_weights_opt);
    Rcpp::CharacterVector components(plan["components"]);
    Rcpp::NumericVector weights(plan["weights"]);
    std::vector<std::string> component_labels = character_vector_to_std(components);
    std::vector<double> weight_vec(component_labels.size());
    for (std::size_t j = 0; j < component_labels.size(); ++j) {
      weight_vec[j] = (j < static_cast<std::size_t>(weights.size()))
        ? static_cast<double>(weights[static_cast<R_xlen_t>(j)])
        : 0.0;
    }
    std::vector<ComponentCacheEntry> plan_cache_entries = build_component_cache_entries(component_labels,
                                                                                       component_cache_map);
    double mix = std::numeric_limits<double>::quiet_NaN();
    if (entry_type == "alias_sum") {
      if (!entry.containsElementNamed("alias_sources") || !std::isfinite(t) || t < 0.0) {
        per_trial[i] = R_NegInf;
        hit_neg_inf = true;
        continue;
      }
      Rcpp::List alias_sources(entry["alias_sources"]);
      mix = evaluate_alias_mix(ctxSEXP,
                               components,
                               weights,
                               alias_sources,
                               t,
                               forced_component,
                               params_ptr,
                               component_cache_map);
    } else if (entry_type == "na_map") {
      if (!entry.containsElementNamed("na_sources")) {
        per_trial[i] = R_NegInf;
        hit_neg_inf = true;
        continue;
      }
      Rcpp::List na_sources(entry["na_sources"]);
      mix = evaluate_na_map_mix(ctxSEXP,
                                components,
                                weights,
                                na_sources,
                                t,
                                component_deadlines,
                                default_deadline,
                                rel_tol,
                                abs_tol,
                                max_depth,
                                forced_component,
                                params_ptr,
                                component_cache_map);
    } else {
      SharedGateSpec shared_spec;
      bool has_shared_gate = extract_shared_gate_spec(entry, shared_spec);
      if (has_shared_gate) {
        mix = evaluate_shared_gate_mix(*ctx,
                                       shared_spec,
                                       t,
                                       components,
                                       weights,
                                       forced_component,
                                        params_ptr,
                                        rel_tol,
                                        abs_tol,
                                        max_depth,
                                        component_cache_map);
      } else {
        if (!entry.containsElementNamed("node_id") ||
            !entry.containsElementNamed("competitor_ids") ||
            !std::isfinite(t) || t < 0.0) {
          per_trial[i] = R_NegInf;
          hit_neg_inf = true;
          continue;
        }
        int node_id = Rcpp::as<int>(entry["node_id"]);
        Rcpp::IntegerVector comp_ids(entry["competitor_ids"]);
        mix = native_trial_mixture_internal(*ctx,
                                            node_id,
                                            t,
                                            component_labels,
                                            weight_vec,
                                            forced_component,
                                            integer_vector_to_std(comp_ids, false),
                                            params_ptr,
                                            plan_cache_entries);
      }
    }
    double log_contrib = R_NegInf;
    if (std::isfinite(mix) && mix > 0.0) {
      log_contrib = std::log(mix);
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
Rcpp::List native_loglik_from_buffer_cpp(SEXP ctxSEXP,
                                         Rcpp::List structure,
                                         Rcpp::List trial_entries,
                                         Rcpp::DataFrame params_df,
                                         Rcpp::Nullable<Rcpp::DataFrame> component_weights_opt,
                                         double default_deadline,
                                         double rel_tol,
                                         double abs_tol,
                                         int max_depth) {
  (void) params_df;
  return native_loglik_from_params_cpp(ctxSEXP,
                                       structure,
                                       trial_entries,
                                       component_weights_opt,
                                       default_deadline,
                                       rel_tol,
                                       abs_tol,
                                       max_depth);
}

// [[Rcpp::export]]
Rcpp::List native_component_plan_exported(SEXP structureSEXP,
                                          SEXP trial_rowsSEXP,
                                          SEXP forced_componentSEXP,
                                          SEXP component_weightsSEXP) {
  Rcpp::List structure = Rcpp::as<Rcpp::List>(structureSEXP);
  Rcpp::Nullable<Rcpp::DataFrame> trial_rows;
  if (!Rf_isNull(trial_rowsSEXP)) {
    trial_rows = Rcpp::Nullable<Rcpp::DataFrame>(trial_rowsSEXP);
  }
  Rcpp::Nullable<Rcpp::String> forced_component;
  if (!Rf_isNull(forced_componentSEXP)) {
    forced_component = Rcpp::Nullable<Rcpp::String>(forced_componentSEXP);
  }
  Rcpp::Nullable<Rcpp::DataFrame> component_weights;
  if (!Rf_isNull(component_weightsSEXP)) {
    component_weights = Rcpp::Nullable<Rcpp::DataFrame>(component_weightsSEXP);
  }
  return native_component_plan_impl(structure,
                                    trial_rows,
                                    forced_component,
                                    component_weights);
}

// [[Rcpp::export]]
Rcpp::List native_guard_eval_cpp(SEXP ctxSEXP,
                                 int guard_node_id,
                                 double t,
                                 Rcpp::Nullable<Rcpp::String> component,
                                 SEXP forced_complete,
                                 SEXP forced_survive,
                                 double rel_tol,
                                 double abs_tol,
                                 int max_depth) {
  return native_guard_eval_impl(ctxSEXP,
                                guard_node_id,
                                t,
                                component,
                                forced_complete,
                                forced_survive,
                                rel_tol,
                                abs_tol,
                                max_depth);
}

// [[Rcpp::export]]
double native_guard_effective_survival_cpp(SEXP ctxSEXP,
                                           int guard_node_id,
                                           double t,
                                           Rcpp::Nullable<Rcpp::String> component,
                                           SEXP forced_complete,
                                           SEXP forced_survive,
                                           double rel_tol,
                                           double abs_tol,
                                           int max_depth) {
  return native_guard_effective_survival_impl(ctxSEXP,
                                              guard_node_id,
                                              t,
                                              component,
                                              forced_complete,
                                              forced_survive,
                                              rel_tol,
                                              abs_tol,
                                              max_depth);
}

// [[Rcpp::export]]
Rcpp::List native_guard_scenarios_cpp(SEXP ctxSEXP,
                                      int guard_node_id,
                                      double t,
                                      Rcpp::List reference_scenarios,
                                      Rcpp::Nullable<Rcpp::String> component,
                                      SEXP forced_complete,
                                      SEXP forced_survive,
                                      double rel_tol,
                                      double abs_tol,
                                      int max_depth) {
  return native_guard_scenarios_impl(ctxSEXP,
                                     guard_node_id,
                                     t,
                                     reference_scenarios,
                                     component,
                                     forced_complete,
                                     forced_survive,
                                     rel_tol,
                                     abs_tol,
                                     max_depth);
}

// [[Rcpp::export]]
Rcpp::List native_node_eval_cpp(SEXP ctxSEXP,
                                int node_id,
                                double t,
                                Rcpp::Nullable<Rcpp::String> component,
                                SEXP forced_complete,
                                SEXP forced_survive) {
  return native_node_eval_impl(ctxSEXP,
                               node_id,
                               t,
                               component,
                               forced_complete,
                               forced_survive);
}

// [[Rcpp::export]]
Rcpp::List native_node_scenarios_cpp(SEXP ctxSEXP,
                                     int node_id,
                                     double t,
                                     Rcpp::Nullable<Rcpp::String> component,
                                     SEXP forced_complete,
                                     SEXP forced_survive) {
  return native_node_scenarios_impl(ctxSEXP,
                                    node_id,
                                    t,
                                    component,
                                    forced_complete,
                                    forced_survive);
}

// [[Rcpp::export]]
Rcpp::List native_likelihood_driver_cpp(SEXP ctxSEXP,
                                        Rcpp::IntegerVector node_ids,
                                        Rcpp::NumericVector times,
                                        Rcpp::Nullable<Rcpp::String> component,
                                        SEXP forced_complete,
                                        SEXP forced_survive) {
  return native_likelihood_driver_impl(ctxSEXP,
                                       node_ids,
                                       times,
                                       component,
                                       forced_complete,
                                       forced_survive);
}

// [[Rcpp::export]]
Rcpp::List native_likelihood_eval_cpp(SEXP ctxSEXP,
                                      Rcpp::List task_list) {
  return native_likelihood_eval_impl(ctxSEXP, task_list);
}

// [[Rcpp::export]]
double native_competitor_survival_cpp(SEXP ctxSEXP,
                                      Rcpp::IntegerVector competitor_ids,
                                      double t,
                                      Rcpp::Nullable<Rcpp::String> component) {
  return native_competitor_survival_impl(ctxSEXP,
                                         competitor_ids,
                                         t,
                                         component);
}

// [[Rcpp::export]]
Rcpp::List native_density_with_competitors_cpp(SEXP ctxSEXP,
                                               int node_id,
                                               double t,
                                               Rcpp::Nullable<Rcpp::String> component,
                                               SEXP forced_complete,
                                               SEXP forced_survive,
                                               Rcpp::IntegerVector competitor_ids,
                                               double rel_tol,
                                               double abs_tol,
                                               int max_depth) {
  (void)rel_tol;
  (void)abs_tol;
  (void)max_depth;
  return native_density_with_competitors_impl(ctxSEXP,
                                              node_id,
                                              t,
                                              component,
                                              forced_complete,
                                              forced_survive,
                                              competitor_ids,
                                              nullptr);
}

// [[Rcpp::export]]
double native_outcome_probability_cpp(SEXP ctxSEXP,
                                      int node_id,
                                      double upper,
                                      Rcpp::Nullable<Rcpp::String> component,
                                      SEXP forced_complete,
                                      SEXP forced_survive,
                                      Rcpp::IntegerVector competitor_ids,
                                      double rel_tol,
                                      double abs_tol,
                                      int max_depth) {
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
                                         nullptr);
}

// [[Rcpp::export]]
Rcpp::List native_density_with_competitors_params_cpp(SEXP ctxSEXP,
                                                      int node_id,
                                                      double t,
                                                      Rcpp::Nullable<Rcpp::String> component,
                                                      SEXP forced_complete,
                                                      SEXP forced_survive,
                                                      Rcpp::IntegerVector competitor_ids,
                                                      Rcpp::Nullable<Rcpp::DataFrame> trial_rows) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::unique_ptr<TrialParamSet> params_holder = build_trial_params_from_df(*ctx, trial_rows);
  return native_density_with_competitors_impl(ctxSEXP,
                                              node_id,
                                              t,
                                              component,
                                              forced_complete,
                                              forced_survive,
                                              competitor_ids,
                                              params_holder ? params_holder.get() : nullptr);
}

Rcpp::NumericVector native_density_with_competitors_vector_impl(SEXP ctxSEXP,
                                                                int node_id,
                                                                const Rcpp::NumericVector& times,
                                                                Rcpp::Nullable<Rcpp::String> component,
                                                                SEXP forced_complete,
                                                                SEXP forced_survive,
                                                                const Rcpp::IntegerVector& competitor_ids,
                                                                const TrialParamSet* trial_params) {
  if (Rf_isNull(ctxSEXP)) {
    Rcpp::stop("native_density_with_competitors_vector: context pointer is null");
  }
  if (times.size() == 0) {
    return Rcpp::NumericVector();
  }
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  std::string component_label = component.isNotNull() ? Rcpp::as<std::string>(component) : std::string();
  std::vector<int> fc_vec = forced_vec_from_sexp(forced_complete);
  std::vector<int> fs_vec = forced_vec_from_sexp(forced_survive);
  std::unordered_set<int> forced_complete_set = make_forced_set(fc_vec);
  std::unordered_set<int> forced_survive_set = make_forced_set(fs_vec);
  EvalCache eval_cache;
  Rcpp::NumericVector out(times.size());
  for (R_xlen_t i = 0; i < times.size(); ++i) {
    double t = times[i];
    if (!std::isfinite(t) || t < 0.0) {
      out[i] = 0.0;
      continue;
    }
    NodeEvalState state(*ctx,
                        eval_cache,
                        t,
                        component_label,
                        forced_complete_set,
                        forced_survive_set,
                        trial_params);
    NodeEvalResult base = eval_node_recursive(node_id, state);
    double density = base.density;
    if (!std::isfinite(density) || density <= 0.0) {
      out[i] = 0.0;
      continue;
    }
    if (competitor_ids.size() > 0) {
      double surv = native_competitor_survival_impl(ctxSEXP,
                                                    competitor_ids,
                                                    t,
                                                    component);
      if (!std::isfinite(surv) || surv <= 0.0) {
        density = 0.0;
      } else {
        density *= surv;
      }
    }
    out[i] = density;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector native_density_with_competitors_vec_cpp(SEXP ctxSEXP,
                                                            int node_id,
                                                            Rcpp::NumericVector times,
                                                            Rcpp::Nullable<Rcpp::String> component,
                                                            SEXP forced_complete,
                                                            SEXP forced_survive,
                                                            Rcpp::IntegerVector competitor_ids) {
  return native_density_with_competitors_vector_impl(ctxSEXP,
                                                     node_id,
                                                     times,
                                                     component,
                                                     forced_complete,
                                                     forced_survive,
                                                     competitor_ids,
                                                     nullptr);
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

inline AccDistParams resolve_acc_params(const std::string& dist,
                                        const Rcpp::List& params) {
  std::vector<uuber::ProtoParamEntry> entries = params_from_rcpp(params, dist);
  try {
    return resolve_acc_params_entries(dist, entries);
  } catch (const std::exception& e) {
    Rcpp::stop("%s", e.what());
  }
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

inline std::vector<std::vector<double>> build_prefix(const std::vector<double>& surv,
                                                     const std::vector<double>& fail) {
  std::size_t n = surv.size();
  std::vector<std::vector<double>> prefix(n + 1);
  prefix[0] = {1.0};
  for (std::size_t i = 0; i < n; ++i) {
    prefix[i + 1] = expand_poly(prefix[i], surv[i], fail[i]);
  }
  return prefix;
}

inline std::vector<std::vector<double>> build_suffix(const std::vector<double>& surv,
                                                     const std::vector<double>& fail) {
  std::size_t n = surv.size();
  std::vector<std::vector<double>> suffix(n + 1);
  suffix[n] = {1.0};
  for (std::size_t i = n; i-- > 0;) {
    suffix[i] = expand_poly(suffix[i + 1], surv[i], fail[i]);
  }
  return suffix;
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
                                             const Rcpp::IntegerVector& addition) {
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
  if (n < 0) {
    Rcpp::stop("n must be non-negative");
  }
  if (is_invalid_positive(sdlog)) {
    Rcpp::stop("sdlog must be positive and finite");
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
  if (n < 0) {
    Rcpp::stop("n must be non-negative");
  }
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    Rcpp::stop("shape and rate must be positive and finite");
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
  if (n < 0) {
    Rcpp::stop("n must be non-negative");
  }
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    Rcpp::stop("sigma and tau must be positive and finite");
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector normals = Rcpp::rnorm(n, mu, sigma);
  Rcpp::NumericVector expo = Rcpp::rexp(n, 1.0 / tau);
  for (R_xlen_t i = 0; i < normals.size(); ++i) {
    normals[i] += expo[i];
  }
  return normals;
}

//' @noRd
// [[Rcpp::export]]
double acc_density_cpp(double t,
                       double onset,
                       double q,
                       const std::string& dist,
                       const Rcpp::List& params) {
  AccDistParams cfg = resolve_acc_params(dist, params);
  return acc_density_from_cfg(t, onset, q, cfg);
}

//' @noRd
// [[Rcpp::export]]
double acc_density_success_cpp(double t,
                               double onset,
                               double q,
                               const std::string& dist,
                               const Rcpp::List& params) {
  (void)q;
  AccDistParams cfg = resolve_acc_params(dist, params);
  return acc_density_success_from_cfg(t, onset, cfg);
}

//' @noRd
// [[Rcpp::export]]
double acc_survival_cpp(double t,
                        double onset,
                        double q,
                        const std::string& dist,
                        const Rcpp::List& params) {
  AccDistParams cfg = resolve_acc_params(dist, params);
  return acc_survival_from_cfg(t, onset, q, cfg);
}

//' @noRd
// [[Rcpp::export]]
double acc_cdf_success_cpp(double t,
                           double onset,
                           double q,
                           const std::string& dist,
                           const Rcpp::List& params) {
  (void)q;
  AccDistParams cfg = resolve_acc_params(dist, params);
  return acc_cdf_success_from_cfg(t, onset, cfg);
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
  if (Svec.size() != Fvec.size()) {
    Rcpp::stop("pool_coeffs expects Svec and Fvec of equal length");
  }
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

//' @noRd
// [[Rcpp::export]]
double pool_density_fast_cpp(const Rcpp::NumericVector& density,
                             const Rcpp::NumericVector& survival,
                             int k) {
  const std::size_t n = density.size();
  if (survival.size() != density.size()) {
    Rcpp::stop("pool_density_fast_cpp expects matching vector lengths");
  }
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

  auto prefix = build_prefix(surv, fail);
  auto suffix = build_suffix(surv, fail);

  PoolDensityWorker worker(density, prefix, suffix, k - 1);
  RcppParallel::parallelReduce(0, n, worker);

  double total = worker.total;
  if (!std::isfinite(total) || total < 0.0) return 0.0;
  return total;
}

//' @noRd
// [[Rcpp::export]]
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

  auto prefix = build_prefix(surv, fail);
  const std::vector<double>& poly = prefix.back();
  int upto = std::min<int>(k, static_cast<int>(poly.size()));
  if (upto <= 0) return 0.0;
  double total = 0.0;
  for (int i = 0; i < upto; ++i) {
    total += poly[static_cast<std::size_t>(i)];
  }
  return clamp(total, 0.0, 1.0);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::List pool_build_templates_cpp(int n,
                                    const Rcpp::IntegerVector& member_ids,
                                    int pool_idx,
                                    int k) {
  Rcpp::List finisher_map(n);
  std::vector<Rcpp::List> templates_vec;
  if (n <= 0 || k < 1 || k > n) {
    return Rcpp::List::create(
      Rcpp::Named("templates") = Rcpp::List(),
      Rcpp::Named("finisher_map") = finisher_map
    );
  }
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
      finisher_map[idx] = Rcpp::IntegerVector(0);
      continue;
    }
    std::vector<std::vector<int>> combos = generate_combinations(others, need);
    std::vector<int> idx_entries;
    idx_entries.reserve(combos.size());

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

      Rcpp::IntegerVector complete_idx;
      if (!combo.empty()) {
        complete_idx = Rcpp::IntegerVector(combo.begin(), combo.end());
      } else {
        complete_idx = Rcpp::IntegerVector(0);
      }
      Rcpp::IntegerVector survivor_idx;
      if (!survivors.empty()) {
        survivor_idx = Rcpp::IntegerVector(survivors.begin(), survivors.end());
      } else {
        survivor_idx = Rcpp::IntegerVector(0);
      }
      Rcpp::IntegerVector forced_complete_vec;
      if (!forced_complete_ids.empty()) {
        forced_complete_vec = Rcpp::IntegerVector(forced_complete_ids.begin(), forced_complete_ids.end());
      } else {
        forced_complete_vec = Rcpp::IntegerVector(0);
      }
      Rcpp::IntegerVector forced_survive_vec;
      if (!forced_survive_ids.empty()) {
        forced_survive_vec = Rcpp::IntegerVector(forced_survive_ids.begin(), forced_survive_ids.end());
      } else {
        forced_survive_vec = Rcpp::IntegerVector(0);
      }

      ++template_counter;
      idx_entries.push_back(template_counter);

      templates_vec.emplace_back(
        Rcpp::List::create(
          Rcpp::Named("finisher_idx") = idx + 1,
          Rcpp::Named("complete_idx") = complete_idx,
          Rcpp::Named("survivor_idx") = survivor_idx,
          Rcpp::Named("forced_complete_ids") = forced_complete_vec,
          Rcpp::Named("forced_survive_ids") = forced_survive_vec
        )
      );
    }
    if (!idx_entries.empty()) {
      finisher_map[idx] = Rcpp::IntegerVector(idx_entries.begin(), idx_entries.end());
    } else {
      finisher_map[idx] = Rcpp::IntegerVector(0);
    }
  }

  Rcpp::List templates_out(templates_vec.size());
  for (std::size_t i = 0; i < templates_vec.size(); ++i) {
    templates_out[i] = templates_vec[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("templates") = templates_out,
    Rcpp::Named("finisher_map") = finisher_map
  );
}

//' @noRd
// [[Rcpp::export]]
Rcpp::List pool_density_combine_cpp(const Rcpp::NumericVector& dens_vec,
                                    const Rcpp::NumericVector& cdf_vec,
                                    const Rcpp::NumericVector& surv_vec,
                                    const Rcpp::NumericVector& cdf_success_vec,
                                    const Rcpp::NumericVector& surv_success_vec,
                                    const Rcpp::IntegerVector& shared_index,
                                    const Rcpp::List& templates,
                                    const Rcpp::IntegerVector& forced_complete,
                                    const Rcpp::IntegerVector& forced_survive) {
  const std::size_t n = dens_vec.size();
  const std::size_t template_count = templates.size();
  const double eps = std::numeric_limits<double>::epsilon();

  std::vector<int> base_fc = integer_vector_to_std(forced_complete, false);
  std::vector<int> base_fs = integer_vector_to_std(forced_survive, false);
  std::vector<int> shared_group(n, 0);
  for (std::size_t i = 0; i < n && i < static_cast<std::size_t>(shared_index.size()); ++i) {
    int val = shared_index[i];
    shared_group[i] = (val == NA_INTEGER) ? 0 : val;
  }

  auto same_shared = [&](int i, int j) {
    if (i < 0 || j < 0 || i >= static_cast<int>(shared_group.size()) ||
        j >= static_cast<int>(shared_group.size())) {
      return false;
    }
    int si = shared_group[static_cast<std::size_t>(i)];
    int sj = shared_group[static_cast<std::size_t>(j)];
    return si > 0 && sj > 0 && si == sj;
  };

  std::vector<Rcpp::List> scenario_vec;
  scenario_vec.reserve(template_count);
  double total = 0.0;

  for (std::size_t t_idx = 0; t_idx < template_count; ++t_idx) {
    Rcpp::List tpl = templates[t_idx];
    int finisher_raw = tpl["finisher_idx"];
    if (finisher_raw == NA_INTEGER) continue;
    int finisher_idx = finisher_raw - 1;
    if (finisher_idx < 0 || finisher_idx >= static_cast<int>(n)) continue;

    double dens_mid = dens_vec[finisher_idx];
    if (!std::isfinite(dens_mid) || dens_mid <= 0.0) continue;

    double weight = dens_mid;

    std::vector<int> complete_idx = integer_vector_to_std(tpl["complete_idx"], true);
    std::vector<int> survivor_idx = integer_vector_to_std(tpl["survivor_idx"], true);

    for (int j : complete_idx) {
      if (j < 0 || j >= static_cast<int>(cdf_vec.size())) continue;
      weight *= cdf_vec[j];
    }
    for (int j : survivor_idx) {
      if (j < 0 || j >= static_cast<int>(surv_vec.size())) continue;
      weight *= surv_vec[j];
    }
    if (!std::isfinite(weight) || weight <= 0.0) continue;

    for (int j : complete_idx) {
      if (j < 0 || j >= static_cast<int>(cdf_vec.size())) continue;
      if (!same_shared(finisher_idx, j)) continue;
      double denom = std::max(cdf_vec[j], eps);
      double ratio = cdf_success_vec[j] / denom;
      weight *= ratio;
    }
    for (int j : survivor_idx) {
      if (j < 0 || j >= static_cast<int>(surv_vec.size())) continue;
      if (!same_shared(finisher_idx, j)) continue;
      double denom = std::max(surv_vec[j], eps);
      double ratio = surv_success_vec[j] / denom;
      weight *= ratio;
    }
    if (!std::isfinite(weight) || weight <= 0.0) continue;

    std::vector<int> merged_fc = merge_forced_vectors(base_fc, tpl["forced_complete_ids"]);
    std::vector<int> merged_fs = merge_forced_vectors(base_fs, tpl["forced_survive_ids"]);

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

//' @noRd
// [[Rcpp::export]]
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

//' @noRd
// [[Rcpp::export]]
double guard_effective_survival_cpp(Rcpp::Function integrand,
                                    double upper,
                                    double rel_tol,
                                    double abs_tol,
                                    int max_depth) {
  if (!std::isfinite(upper)) return 0.0;
  if (upper <= 0.0) return 1.0;
  if (rel_tol <= 0.0) rel_tol = 1e-6;
  if (abs_tol <= 0.0) abs_tol = 1e-8;
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
