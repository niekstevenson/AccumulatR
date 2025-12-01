#pragma once

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>
#include <utility>

// Detect availability of Boost Gauss-Kronrod header if not already set
#if !defined(UUBER_HAVE_BOOST_GK)
# if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#  define UUBER_HAVE_BOOST_GK 1
# else
#  define UUBER_HAVE_BOOST_GK 0
# endif
#endif

#if UUBER_HAVE_BOOST_GK
#include <boost/math/quadrature/gauss_kronrod.hpp>
#endif

namespace uuber {

// Fixed 15-point Gauss-Legendre rule on [lower, upper]; no adaptivity.
template <typename Fn>
inline double integrate_fixed_gauss15(Fn&& integrand,
                                      double lower,
                                      double upper) {
  if (!std::isfinite(lower) || !std::isfinite(upper)) return 0.0;
  if (lower == upper) return 0.0;
  static const double x[15] = {
    -0.9879925180204854, -0.9372733924007060, -0.8482065834104272,
    -0.7244177313601700, -0.5709721726085388, -0.3941513470775634,
    -0.2011940939974345,  0.0,
     0.2011940939974345,  0.3941513470775634,  0.5709721726085388,
     0.7244177313601700,  0.8482065834104272,  0.9372733924007060,
     0.9879925180204854
  };
  static const double w[15] = {
    0.0307532419961173, 0.0703660474881081, 0.1071592204671719,
    0.1395706779261543, 0.1662692058169939, 0.1861610000155622,
    0.1984314853271116, 0.2025782419255613, 0.1984314853271116,
    0.1861610000155622, 0.1662692058169939, 0.1395706779261543,
    0.1071592204671719, 0.0703660474881081, 0.0307532419961173
  };
  const double a = 0.5 * (upper - lower);
  const double b = 0.5 * (upper + lower);
  double sum = 0.0;
  for (int i = 0; i < 15; ++i) {
    double xi = b + a * x[i];
    double vi = integrand(xi);
    if (!std::isfinite(vi)) vi = 0.0;
    sum += w[i] * vi;
  }
  double res = a * sum;
  return std::isfinite(res) ? res : 0.0;
}

inline double integrate_boost(Rcpp::Function integrand,
                              double lower,
                              double upper,
                              double rel_tol,
                              double abs_tol,
                              int max_depth) {
  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    Rcpp::stop("boost_integrate_cpp requires finite bounds");
  }
  if (max_depth <= 0) max_depth = 10;
  if (rel_tol <= 0.0) rel_tol = 1e-4;
  if (abs_tol <= 0.0) abs_tol = 1e-5;
  double tol = std::max(abs_tol, rel_tol);
  if (!(tol > 0.0)) tol = 1e-6;

  double a = lower;
  double b = upper;
  double sign = 1.0;
  if (b < a) {
    std::swap(a, b);
    sign = -1.0;
  }

  auto wrapper = [&](double x) -> double {
    if (!std::isfinite(x)) return 0.0;
    Rcpp::NumericVector res = integrand(Rcpp::NumericVector::create(x));
    if (res.size() == 0) return 0.0;
    double val = res[0];
    if (!std::isfinite(val)) return 0.0;
    return val;
  };

  double integral = 0.0;
  #if UUBER_HAVE_BOOST_GK
  try {
    integral = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
      wrapper,
      a,
      b,
      tol,
      max_depth
    );
  } catch (...) {
    integral = 0.0;
  }
  #else
  Rcpp::stop("Boost Gauss-Kronrod header not available; install BH or set include path to Boost");
  #endif
  if (!std::isfinite(integral)) integral = 0.0;
  return sign * integral;
}

template <typename Fn>
inline double integrate_boost_fn(Fn&& integrand,
                                 double lower,
                                 double upper,
                                 double rel_tol,
                                 double abs_tol,
                                 int max_depth) {
  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    Rcpp::stop("integrate_boost_fn requires finite bounds");
  }
  if (max_depth <= 0) max_depth = 10;
  if (rel_tol <= 0.0) rel_tol = 1e-5;
  if (abs_tol <= 0.0) abs_tol = 1e-6;
  double tol = std::max(abs_tol, rel_tol);
  if (!(tol > 0.0)) tol = 1e-6;
  double a = lower;
  double b = upper;
  double sign = 1.0;
  if (b < a) {
    std::swap(a, b);
    sign = -1.0;
  }
  auto wrapper = [&](double x) -> double {
    if (!std::isfinite(x)) return 0.0;
    double val = integrand(x);
    if (!std::isfinite(val)) return 0.0;
    return val;
  };
  double integral = 0.0;
  #if UUBER_HAVE_BOOST_GK
  try {
    integral = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
      wrapper,
      a,
      b,
      tol,
      max_depth
    );
  } catch (...) {
    integral = 0.0;
  }
  #else
  Rcpp::stop("Boost Gauss-Kronrod header not available; install BH or set include path to Boost");
  #endif
  if (!std::isfinite(integral)) integral = 0.0;
  return sign * integral;
}

// Variant using a lower-order Gauss-Kronrod rule (21-point) to reduce per-interval evaluations.

} // namespace uuber
