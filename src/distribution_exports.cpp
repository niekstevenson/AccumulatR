#include <Rcpp.h>

#include <cmath>
#include <vector>

#include "accumulator.h"
#include "distribution_core.h"
#include "dist_vector.h"

namespace {

inline uuber::AccDistParams make_cfg(int dist_code, double p1, double p2,
                                     double p3, double p4, double p5 = 0.0,
                                     double p6 = 0.0, double p7 = 0.0,
                                     double p8 = 0.0) {
  uuber::AccDistParams cfg{};
  cfg.code = dist_code;
  cfg.p1 = p1;
  cfg.p2 = p2;
  cfg.p3 = p3;
  cfg.p4 = p4;
  cfg.p5 = p5;
  cfg.p6 = p6;
  cfg.p7 = p7;
  cfg.p8 = p8;
  return cfg;
}

Rcpp::NumericVector evaluate_pdf_vector(const Rcpp::NumericVector &x,
                                        int dist_code, double p1, double p2,
                                        double p3, double p4, double p5 = 0.0,
                                        double p6 = 0.0, double p7 = 0.0,
                                        double p8 = 0.0) {
  Rcpp::NumericVector out(x.size());
  std::vector<double> finite_values;
  std::vector<R_xlen_t> finite_indices;
  finite_values.reserve(static_cast<std::size_t>(x.size()));
  finite_indices.reserve(static_cast<std::size_t>(x.size()));

  for (R_xlen_t i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (Rcpp::NumericVector::is_na(value)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!R_FINITE(value)) {
      out[i] = 0.0;
      continue;
    }
    finite_indices.push_back(i);
    finite_values.push_back(value);
  }

  if (!finite_values.empty()) {
    std::vector<double> finite_out(finite_values.size(), 0.0);
    uuber::eval_pdf_vec(dist_code, p1, p2, p3, p4, p5, p6, p7, p8,
                        finite_values.data(), finite_values.size(),
                        finite_out.data());
    for (std::size_t i = 0; i < finite_indices.size(); ++i) {
      out[finite_indices[i]] = finite_out[i];
    }
  }

  return out;
}

Rcpp::NumericVector evaluate_cdf_vector(const Rcpp::NumericVector &x,
                                        int dist_code, double p1, double p2,
                                        double p3, double p4, double p5 = 0.0,
                                        double p6 = 0.0, double p7 = 0.0,
                                        double p8 = 0.0) {
  Rcpp::NumericVector out(x.size());
  std::vector<double> finite_values;
  std::vector<R_xlen_t> finite_indices;
  finite_values.reserve(static_cast<std::size_t>(x.size()));
  finite_indices.reserve(static_cast<std::size_t>(x.size()));

  for (R_xlen_t i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (Rcpp::NumericVector::is_na(value)) {
      out[i] = NA_REAL;
      continue;
    }
    if (!R_FINITE(value)) {
      out[i] = value < 0.0 ? 0.0 : 1.0;
      continue;
    }
    finite_indices.push_back(i);
    finite_values.push_back(value);
  }

  if (!finite_values.empty()) {
    std::vector<double> finite_out(finite_values.size(), 0.0);
    uuber::eval_cdf_vec(dist_code, p1, p2, p3, p4, p5, p6, p7, p8,
                        finite_values.data(), finite_values.size(),
                        finite_out.data());
    for (std::size_t i = 0; i < finite_indices.size(); ++i) {
      out[finite_indices[i]] = finite_out[i];
    }
  }

  return out;
}

Rcpp::NumericVector evaluate_pdf_vector_with_lower_bound(
    const Rcpp::NumericVector &x, int dist_code, double lower_bound, double p1,
    double p2, double p3, double p4, double p5 = 0.0, double p6 = 0.0,
    double p7 = 0.0, double p8 = 0.0) {
  Rcpp::NumericVector out(x.size());
  std::vector<double> finite_values;
  std::vector<R_xlen_t> finite_indices;
  finite_values.reserve(static_cast<std::size_t>(x.size()));
  finite_indices.reserve(static_cast<std::size_t>(x.size()));

  for (R_xlen_t i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (Rcpp::NumericVector::is_na(value)) {
      out[i] = NA_REAL;
      continue;
    }
    finite_indices.push_back(i);
    finite_values.push_back(value);
  }

  if (!finite_values.empty()) {
    const uuber::AccDistParams cfg =
        make_cfg(dist_code, p1, p2, p3, p4, p5, p6, p7, p8);
    const LowerBoundTransform transform =
        make_lower_bound_transform(cfg, lower_bound);
    std::vector<double> finite_out(finite_values.size(), 0.0);
    eval_pdf_vec_with_lower_bound(cfg, transform, finite_values.data(),
                                  finite_values.size(), finite_out.data());
    for (std::size_t i = 0; i < finite_indices.size(); ++i) {
      out[finite_indices[i]] = finite_out[i];
    }
  }

  return out;
}

Rcpp::NumericVector evaluate_cdf_vector_with_lower_bound(
    const Rcpp::NumericVector &x, int dist_code, double lower_bound, double p1,
    double p2, double p3, double p4, double p5 = 0.0, double p6 = 0.0,
    double p7 = 0.0, double p8 = 0.0) {
  Rcpp::NumericVector out(x.size());
  std::vector<double> finite_values;
  std::vector<R_xlen_t> finite_indices;
  finite_values.reserve(static_cast<std::size_t>(x.size()));
  finite_indices.reserve(static_cast<std::size_t>(x.size()));

  for (R_xlen_t i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (Rcpp::NumericVector::is_na(value)) {
      out[i] = NA_REAL;
      continue;
    }
    finite_indices.push_back(i);
    finite_values.push_back(value);
  }

  if (!finite_values.empty()) {
    const uuber::AccDistParams cfg =
        make_cfg(dist_code, p1, p2, p3, p4, p5, p6, p7, p8);
    const LowerBoundTransform transform =
        make_lower_bound_transform(cfg, lower_bound);
    std::vector<double> finite_out(finite_values.size(), 0.0);
    eval_cdf_vec_with_lower_bound(cfg, transform, finite_values.data(),
                                  finite_values.size(), finite_out.data());
    for (std::size_t i = 0; i < finite_indices.size(); ++i) {
      out[finite_indices[i]] = finite_out[i];
    }
  }

  return out;
}

inline double rdm_rng_single(double k, double l, double tiny = 1e-6) {
  if (!std::isfinite(k) || !std::isfinite(l) || l < 0.0) {
    return R_PosInf;
  }
  if (l <= tiny) {
    const double q = std::max(1e-15, 1.0 - R::runif(0.0, 1.0) / 2.0);
    const double z = R::qnorm(q, 0.0, 1.0, 1, 0);
    if (!std::isfinite(z) || std::fabs(z) < 1e-15) {
      return R_PosInf;
    }
    return (k * k) / (z * z);
  }

  const double mu = k / l;
  const double lambda = k * k;
  if (!std::isfinite(mu) || !std::isfinite(lambda) || mu <= 0.0 ||
      lambda <= 0.0) {
    return R_PosInf;
  }

  const double y = std::pow(R::rnorm(0.0, 1.0), 2.0);
  const double root = std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  double x = mu + (mu * mu * y) / (2.0 * lambda) -
             (mu * root) / (2.0 * lambda);
  const double u = R::runif(0.0, 1.0);
  const double accept = mu / (mu + x);
  if (u > accept) {
    x = (mu * mu) / x;
  }
  if (!std::isfinite(x) || x <= 0.0) {
    return R_PosInf;
  }
  return x;
}

} // namespace

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lognormal_pdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog) {
  if (is_invalid_positive(sdlog)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_pdf_vector(x, uuber::ACC_DIST_LOGNORMAL, meanlog, sdlog, 0.0,
                             0.0);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lognormal_cdf(const Rcpp::NumericVector &x,
                                       double meanlog, double sdlog) {
  if (is_invalid_positive(sdlog)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_cdf_vector(x, uuber::ACC_DIST_LOGNORMAL, meanlog, sdlog, 0.0,
                             0.0);
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

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_pdf(const Rcpp::NumericVector &x, double shape,
                                   double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_pdf_vector(x, uuber::ACC_DIST_GAMMA, shape, rate, 0.0, 0.0);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_gamma_cdf(const Rcpp::NumericVector &x, double shape,
                                   double rate) {
  if (is_invalid_positive(shape) || is_invalid_positive(rate)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_cdf_vector(x, uuber::ACC_DIST_GAMMA, shape, rate, 0.0, 0.0);
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

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_pdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau) {
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_pdf_vector(x, uuber::ACC_DIST_EXGAUSS, mu, sigma, tau, 0.0);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_exgauss_cdf(const Rcpp::NumericVector &x, double mu,
                                     double sigma, double tau) {
  if (is_invalid_positive(sigma) || is_invalid_positive(tau)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_cdf_vector(x, uuber::ACC_DIST_EXGAUSS, mu, sigma, tau, 0.0);
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

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_pdf_lower_bound(const Rcpp::NumericVector &x,
                                         int dist_code, double lower_bound,
                                         double p1, double p2, double p3,
                                         double p4, double p5 = 0.0,
                                         double p6 = 0.0, double p7 = 0.0,
                                         double p8 = 0.0) {
  return evaluate_pdf_vector_with_lower_bound(x, dist_code, lower_bound, p1, p2,
                                              p3, p4, p5, p6, p7, p8);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_cdf_lower_bound(const Rcpp::NumericVector &x,
                                         int dist_code, double lower_bound,
                                         double p1, double p2, double p3,
                                         double p4, double p5 = 0.0,
                                         double p6 = 0.0, double p7 = 0.0,
                                         double p8 = 0.0) {
  return evaluate_cdf_vector_with_lower_bound(x, dist_code, lower_bound, p1, p2,
                                              p3, p4, p5, p6, p7, p8);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lba_pdf(const Rcpp::NumericVector &x, double v,
                                 double sv, double B, double A) {
  if (!std::isfinite(v) || is_invalid_positive(sv) || !std::isfinite(B) ||
      !std::isfinite(A)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_pdf_vector(x, uuber::ACC_DIST_LBA, v, sv, B, A);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lba_cdf(const Rcpp::NumericVector &x, double v,
                                 double sv, double B, double A) {
  if (!std::isfinite(v) || is_invalid_positive(sv) || !std::isfinite(B) ||
      !std::isfinite(A)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_cdf_vector(x, uuber::ACC_DIST_LBA, v, sv, B, A);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_lba_rng(int n, double v, double sv, double B,
                                 double A) {
  if (n <= 0 || !std::isfinite(v) || is_invalid_positive(sv) ||
      !std::isfinite(B) || !std::isfinite(A)) {
    return Rcpp::NumericVector();
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector out(n);
  const double lower_mass = R::pnorm(0.0, v, sv, 1, 0);
  const double upper_mass = 1.0 - lower_mass;
  if (!std::isfinite(upper_mass) || upper_mass <= 0.0) {
    std::fill(out.begin(), out.end(), R_PosInf);
    return out;
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    const double start = B - (A * R::runif(0.0, 1.0));
    const double u = R::runif(0.0, 1.0);
    const double q = std::min(1.0 - 1e-15, lower_mass + u * upper_mass);
    const double drift = R::qnorm(q, v, sv, 1, 0);
    const double dt = start / drift;
    out[i] = (std::isfinite(dt) && dt > 0.0) ? dt : R_PosInf;
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_rdm_pdf(const Rcpp::NumericVector &x, double v,
                                 double B, double A, double s) {
  if (!std::isfinite(v) || !std::isfinite(B) || !std::isfinite(A) ||
      is_invalid_positive(s)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_pdf_vector(x, uuber::ACC_DIST_RDM, v, B, A, s);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_rdm_cdf(const Rcpp::NumericVector &x, double v,
                                 double B, double A, double s) {
  if (!std::isfinite(v) || !std::isfinite(B) || !std::isfinite(A) ||
      is_invalid_positive(s)) {
    Rcpp::NumericVector out(x.size(), NA_REAL);
    return out;
  }
  return evaluate_cdf_vector(x, uuber::ACC_DIST_RDM, v, B, A, s);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector dist_rdm_rng(int n, double v, double B, double A,
                                 double s) {
  if (n <= 0 || !std::isfinite(v) || !std::isfinite(B) || !std::isfinite(A) ||
      is_invalid_positive(s)) {
    return Rcpp::NumericVector();
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector out(n);
  const double v_scaled = v / s;
  const double B_scaled = std::max(0.0, B / s);
  const double A_scaled = std::max(0.0, A / s);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double start = B_scaled + R::runif(0.0, A_scaled);
    out[i] = rdm_rng_single(start, v_scaled);
  }
  return out;
}
