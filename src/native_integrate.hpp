#pragma once

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <functional>

namespace uuber {

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
  if (rel_tol <= 0.0) rel_tol = 1e-6;
  if (abs_tol <= 0.0) abs_tol = 1e-8;
  double tol = std::max(abs_tol, rel_tol);
  if (!(tol > 0.0)) tol = 1e-8;

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
  if (rel_tol <= 0.0) rel_tol = 1e-6;
  if (abs_tol <= 0.0) abs_tol = 1e-8;
  double tol = std::max(abs_tol, rel_tol);
  if (!(tol > 0.0)) tol = 1e-8;
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
  if (!std::isfinite(integral)) integral = 0.0;
  return sign * integral;
}

} // namespace uuber
