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

double eval_integrand(Rcpp::Function& fn, double x) {
  Rcpp::NumericVector res = fn(Rcpp::NumericVector::create(x));
  if (res.size() == 0) return 0.0;
  double val = res[0];
  if (!std::isfinite(val) || val < 0.0) return 0.0;
  return val;
}

inline double simpson_estimate(double a, double b, double fa, double fb, double fm) {
  return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

double adaptive_simpson_guard(Rcpp::Function& integrand,
                              double a,
                              double b,
                              double fa,
                              double fb,
                              double fm,
                              double whole,
                              double abs_tol,
                              double rel_tol,
                              int depth) {
  double m = 0.5 * (a + b);
  double left_mid = 0.5 * (a + m);
  double right_mid = 0.5 * (m + b);
  double f_left_mid = eval_integrand(integrand, left_mid);
  double f_right_mid = eval_integrand(integrand, right_mid);
  double left = simpson_estimate(a, m, fa, fm, f_left_mid);
  double right = simpson_estimate(m, b, fm, fb, f_right_mid);
  double result = left + right;
  double error = std::fabs(result - whole);
  double tol = std::max(abs_tol, rel_tol * std::fabs(result));
  if (error <= 15.0 * tol || depth <= 0) {
    return result + (result - whole) / 15.0;
  }
  double half_abs = abs_tol * 0.5;
  return adaptive_simpson_guard(integrand, a, m, fa, fm, f_left_mid, left,
                                half_abs, rel_tol, depth - 1) +
         adaptive_simpson_guard(integrand, m, b, fm, fb, f_right_mid, right,
                                half_abs, rel_tol, depth - 1);
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
  if (max_depth <= 0) max_depth = 12;

  double a = 0.0;
  double b = upper;
  double fa = eval_integrand(integrand, a);
  double fb = eval_integrand(integrand, b);
  double fm = eval_integrand(integrand, 0.5 * (a + b));
  double whole = simpson_estimate(a, b, fa, fb, fm);
  double integral = adaptive_simpson_guard(integrand, a, b, fa, fb, fm, whole,
                                           abs_tol, rel_tol, max_depth);
  if (!std::isfinite(integral) || integral < 0.0) integral = 0.0;
  double surv = 1.0 - integral;
  if (!std::isfinite(surv)) surv = 0.0;
  if (surv < 0.0) surv = 0.0;
  if (surv > 1.0) surv = 1.0;
  return surv;
}
