// [[Rcpp::depends(Rcpp, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <sstream>
#if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#define UUBER_HAVE_BOOST_GK 1
#include <boost/math/quadrature/gauss_kronrod.hpp>
#else
#define UUBER_HAVE_BOOST_GK 0
#endif

#include "native_context.hpp"
#include "native_integrate.hpp"

// [[Rcpp::export]]
SEXP native_context_build(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  return uuber::build_native_context(prep);
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

namespace {

void sort_unique(std::vector<int>& vec);

constexpr double kDefaultRelTol = 1e-5;
constexpr double kDefaultAbsTol = 1e-6;
constexpr int kDefaultMaxDepth = 12;

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

inline int label_id(const uuber::NativeContext& ctx, const std::string& label) {
  auto it = ctx.label_to_id.find(label);
  if (it == ctx.label_to_id.end()) return NA_INTEGER;
  return it->second;
}

bool component_active(const uuber::NativeAccumulator& acc, const std::string& component) {
  if (component.empty() || component == "__default__") return true;
  if (acc.components.empty()) return true;
  return std::find(acc.components.begin(), acc.components.end(), component) != acc.components.end();
}

NodeEvalResult eval_event_label(const uuber::NativeContext& ctx,
                                const std::string& label,
                                double t,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive) {
  if (label.empty()) {
    Rcpp::stop("native_event_eval: empty label");
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
    const uuber::NativeAccumulator& acc = ctx.accumulators[acc_it->second];
    if (!component_active(acc, component)) {
      return make_node_result(0.0, 1.0, 0.0);
    }
    double density = acc_density_cpp(t, acc.onset, acc.q, acc.dist, acc.params);
    double survival = acc_survival_cpp(t, acc.onset, acc.q, acc.dist, acc.params);
    double cdf = acc_cdf_success_cpp(t, acc.onset, acc.q, acc.dist, acc.params);
    if (!std::isfinite(cdf)) {
      cdf = clamp_probability(1.0 - survival);
    }
    if (!std::isfinite(survival)) {
      survival = 0.0;
    }
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
        if (!component_active(acc, component)) continue;
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
        forced_survive
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

struct NodeEvalState {
  NodeEvalState(const uuber::NativeContext& ctx_,
                double time,
                const std::string& component_label,
                const std::unordered_set<int>& forced_complete_ref,
                const std::unordered_set<int>& forced_survive_ref)
    : ctx(ctx_),
      t(time),
      component(component_label),
      forced_complete(forced_complete_ref),
      forced_survive(forced_survive_ref) {}

  const uuber::NativeContext& ctx;
  double t;
  std::string component;
  const std::unordered_set<int>& forced_complete;
  const std::unordered_set<int>& forced_survive;
  std::unordered_map<int, NodeEvalResult> cache;
};

struct IntegrationSettings {
  double rel_tol{kDefaultRelTol};
  double abs_tol{kDefaultAbsTol};
  int max_depth{kDefaultMaxDepth};
};

struct ScenarioRecord {
  double weight{0.0};
  std::vector<int> forced_complete;
  std::vector<int> forced_survive;
};

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state);
NodeEvalResult eval_guard_node(const uuber::NativeNode& node, NodeEvalState& parent_state);

const uuber::NativeNode& fetch_node(const uuber::NativeContext& ctx, int node_id) {
  auto node_it = ctx.node_index.find(node_id);
  if (node_it == ctx.node_index.end()) {
    Rcpp::stop("native_node_eval: unknown node id %d", node_id);
  }
  return ctx.nodes[node_it->second];
}

NodeEvalResult eval_node_with_forced(const uuber::NativeContext& ctx,
                                     int node_id,
                                     double time,
                                     const std::string& component,
                                     const std::unordered_set<int>& forced_complete,
                                     const std::unordered_set<int>& forced_survive) {
  NodeEvalState local(ctx, time, component, forced_complete, forced_survive);
  return eval_node_recursive(node_id, local);
}

NodeEvalResult eval_node_recursive(int node_id, NodeEvalState& state);
NodeEvalResult eval_guard_node(const uuber::NativeNode& node, NodeEvalState& parent_state);

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

struct GuardEvalInput {
  const uuber::NativeContext& ctx;
  const uuber::NativeNode& node;
  std::string component;
  std::unordered_set<int> forced_complete;
  std::unordered_set<int> forced_survive;
};

GuardEvalInput make_guard_input(const uuber::NativeContext& ctx,
                                const uuber::NativeNode& node,
                                const std::string& component,
                                const std::unordered_set<int>& forced_complete,
                                const std::unordered_set<int>& forced_survive) {
  GuardEvalInput input{ctx, node, component, {}, {}};
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
  auto cached = state.cache.find(node_id);
  if (cached != state.cache.end()) {
    return cached->second;
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
      state.forced_survive
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
                                                  state.component,
                                                  state.forced_complete,
                                                  state.forced_survive);
    IntegrationSettings settings;
    double density = guard_density_internal(guard_input, state.t, settings);
    double cdf = guard_cdf_internal(guard_input, state.t, settings);
    double survival = clamp_probability(1.0 - cdf);
    result = make_node_result(density, survival, cdf);
  } else {
    Rcpp::stop("native_node_eval: unsupported node kind '%s'", node.kind);
  }
  state.cache.emplace(node_id, result);
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
  if (input.node.unless_ids.empty()) {
    NodeEvalResult block = eval_node_with_forced(
      input.ctx,
      blocker_id,
      t,
      input.component,
      input.forced_complete,
      input.forced_survive
    );
    return clamp_probability(block.survival);
  }
  auto integrand = [&](double u) -> double {
    if (!std::isfinite(u) || u < 0.0) return 0.0;
    NodeEvalResult block = eval_node_with_forced(
      input.ctx,
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
  return clamp_probability(surv);
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
  return make_guard_input(ctx, node, component, forced_complete, forced_survive);
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
  GuardEvalInput guard_input = resolve_guard_input(*ctx,
                                                   guard_node_id,
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
  GuardEvalInput guard_input = resolve_guard_input(*ctx,
                                                   guard_node_id,
                                                   component_label,
                                                   fc_set,
                                                   fs_set);
  IntegrationSettings settings = resolve_integration_settings(rel_tol, abs_tol, max_depth);
  return guard_effective_survival_internal(guard_input, t, settings);
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
    GuardEvalInput scenario_input = make_guard_input(*ctx,
                                                     guard_node,
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
  NodeEvalState state(*ctx,
                      t,
                      component_label,
                      forced_complete_set,
                      forced_survive_set);
  NodeEvalResult result = eval_node_recursive(node_id, state);
  return node_result_to_list(result);
}

} // namespace

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

enum AccDistKind {
  ACC_DIST_LOGNORMAL = 1,
  ACC_DIST_GAMMA = 2,
  ACC_DIST_EXGAUSS = 3
};

struct AccDistParams {
  int code;
  double p1;
  double p2;
  double p3;
};

inline std::string normalize_dist_name(const std::string& dist) {
  std::string out(dist);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

inline double get_required_param(const Rcpp::List& params,
                                 const char* name,
                                 const std::string& dist) {
  if (!params.containsElementNamed(name)) {
    Rcpp::stop("Distribution '%s' missing required parameter '%s'",
               dist, name);
  }
  SEXP val = params[name];
  if (val == R_NilValue) {
    Rcpp::stop("Distribution '%s' missing required parameter '%s'",
               dist, name);
  }
  return Rcpp::as<double>(val);
}

inline AccDistParams resolve_acc_params(const std::string& dist,
                                        const Rcpp::List& params) {
  if (params.size() == 0) {
    Rcpp::stop("Distribution '%s' received empty parameter list", dist);
  }
  AccDistParams cfg{};
  std::string dist_name = normalize_dist_name(dist);
  if (dist_name == "lognormal") {
    cfg.code = ACC_DIST_LOGNORMAL;
    cfg.p1 = get_required_param(params, "meanlog", dist_name);
    cfg.p2 = get_required_param(params, "sdlog", dist_name);
    cfg.p3 = 0.0;
  } else if (dist_name == "gamma") {
    cfg.code = ACC_DIST_GAMMA;
    cfg.p1 = get_required_param(params, "shape", dist_name);
    cfg.p2 = get_required_param(params, "rate", dist_name);
    cfg.p3 = 0.0;
  } else if (dist_name == "exgauss") {
    cfg.code = ACC_DIST_EXGAUSS;
    cfg.p1 = get_required_param(params, "mu", dist_name);
    cfg.p2 = get_required_param(params, "sigma", dist_name);
    cfg.p3 = get_required_param(params, "tau", dist_name);
  } else {
    Rcpp::stop("Unsupported accumulator distribution '%s'", dist);
  }
  return cfg;
}

inline double eval_pdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case ACC_DIST_LOGNORMAL:
    return dist_lognormal_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case ACC_DIST_GAMMA:
    return dist_gamma_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case ACC_DIST_EXGAUSS:
    return dist_exgauss_pdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    Rcpp::stop("Invalid accumulator distribution code '%d'", cfg.code);
  }
}

inline double eval_cdf_single(const AccDistParams& cfg, double x) {
  switch (cfg.code) {
  case ACC_DIST_LOGNORMAL:
    return dist_lognormal_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case ACC_DIST_GAMMA:
    return dist_gamma_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2)[0];
  case ACC_DIST_EXGAUSS:
    return dist_exgauss_cdf(Rcpp::NumericVector::create(x), cfg.p1, cfg.p2, cfg.p3)[0];
  default:
    Rcpp::stop("Invalid accumulator distribution code '%d'", cfg.code);
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
                                              bool subtract_one = false) {
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

// ------------------------------------------------------------------
// Accumulator helpers built on the native distributions
// ------------------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
double acc_density_cpp(double t,
                       double onset,
                       double q,
                       const std::string& dist,
                       const Rcpp::List& params) {
  if (!std::isfinite(onset)) {
    Rcpp::stop("Accumulator parameter 'onset' must be finite");
  }
  if (!std::isfinite(q)) {
    Rcpp::stop("Accumulator parameter 'q' must be finite");
  }
  if (!std::isfinite(t) || t < 0.0) return 0.0;
  if (t < onset) return 0.0;

  double success_prob = 1.0 - q;
  if (success_prob <= 0.0) return 0.0;

  AccDistParams cfg = resolve_acc_params(dist, params);
  double dens = eval_pdf_single(cfg, t - onset);
  if (Rcpp::NumericVector::is_na(dens) || !std::isfinite(dens)) {
    return NA_REAL;
  }
  return success_prob * dens;
}

//' @noRd
// [[Rcpp::export]]
double acc_density_success_cpp(double t,
                               double onset,
                               double q,
                               const std::string& dist,
                               const Rcpp::List& params) {
  if (!std::isfinite(onset)) {
    Rcpp::stop("Accumulator parameter 'onset' must be finite");
  }
  (void)q;
  if (!std::isfinite(t) || t < 0.0) return 0.0;
  if (t < onset) return 0.0;

  AccDistParams cfg = resolve_acc_params(dist, params);
  double dens = eval_pdf_single(cfg, t - onset);
  if (Rcpp::NumericVector::is_na(dens) || !std::isfinite(dens)) {
    return NA_REAL;
  }
  return dens;
}

//' @noRd
// [[Rcpp::export]]
double acc_survival_cpp(double t,
                        double onset,
                        double q,
                        const std::string& dist,
                        const Rcpp::List& params) {
  if (!std::isfinite(onset)) {
    Rcpp::stop("Accumulator parameter 'onset' must be finite");
  }
  if (!std::isfinite(q)) {
    Rcpp::stop("Accumulator parameter 'q' must be finite");
  }
  if (!std::isfinite(t)) return 0.0;
  if (t < 0.0) return 1.0;
  if (t < onset) return 1.0;

  AccDistParams cfg = resolve_acc_params(dist, params);
  double cdf = eval_cdf_single(cfg, t - onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  double surv_underlying = 1.0 - cdf;
  surv_underlying = clamp(surv_underlying, 0.0, 1.0);
  double success_prob = 1.0 - q;
  double result = q + success_prob * surv_underlying;
  return clamp(result, 0.0, 1.0);
}

//' @noRd
// [[Rcpp::export]]
double acc_cdf_success_cpp(double t,
                           double onset,
                           double q,
                           const std::string& dist,
                           const Rcpp::List& params) {
  if (!std::isfinite(onset)) {
    Rcpp::stop("Accumulator parameter 'onset' must be finite");
  }
  (void)q;
  if (!std::isfinite(t)) return 1.0;
  if (t < 0.0) return 0.0;
  if (t < onset) return 0.0;

  AccDistParams cfg = resolve_acc_params(dist, params);
  double cdf = eval_cdf_single(cfg, t - onset);
  if (Rcpp::NumericVector::is_na(cdf) || !std::isfinite(cdf)) {
    return NA_REAL;
  }
  return clamp(cdf, 0.0, 1.0);
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
