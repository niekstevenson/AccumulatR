#include <Rcpp.h>
#include <cmath>

namespace {

inline double clamp(double val, double lo, double hi) {
  if (!std::isfinite(val)) {
    return val;
  }
  if (val < lo) {
    return lo;
  }
  if (val > hi) {
    return hi;
  }
  return val;
}

inline bool is_invalid_positive(double value) {
  return !std::isfinite(value) || value <= 0.0;
}

inline double normal_cdf_fast(double z) {
  const double arg = -z * 0.7071067811865475;
  double val = 0.5 * std::erfc(arg);
  return clamp(val, 0.0, 1.0);
}

inline double exgauss_pdf_fast(double x, double mu, double sigma, double tau) {
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
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
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
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

} // namespace

// [[Rcpp::export]]
Rcpp::NumericVector exgauss_pdf_cpp(const Rcpp::NumericVector &x, double mu,
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

// [[Rcpp::export]]
Rcpp::NumericVector exgauss_cdf_cpp(const Rcpp::NumericVector &x, double mu,
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

// [[Rcpp::export]]
Rcpp::NumericVector exgauss_rng_cpp(int n, double mu, double sigma,
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
