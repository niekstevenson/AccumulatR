#pragma once

#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <string>
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

#if !defined(UUBER_HAVE_BOOST_TANH_SINH)
# if __has_include(<boost/math/quadrature/tanh_sinh.hpp>)
#  define UUBER_HAVE_BOOST_TANH_SINH 1
# else
#  define UUBER_HAVE_BOOST_TANH_SINH 0
# endif
#endif

#if UUBER_HAVE_BOOST_TANH_SINH
#include <boost/math/quadrature/tanh_sinh.hpp>
#endif

#if !defined(UUBER_HAVE_BOOST_TRAPEZOIDAL)
# if __has_include(<boost/math/quadrature/trapezoidal.hpp>)
#  define UUBER_HAVE_BOOST_TRAPEZOIDAL 1
# else
#  define UUBER_HAVE_BOOST_TRAPEZOIDAL 0
# endif
#endif

#if UUBER_HAVE_BOOST_TRAPEZOIDAL
#include <boost/math/quadrature/trapezoidal.hpp>
#endif

namespace uuber {

enum class BoostQuadratureStrategy {
  kAuto,
  kGk21,
  kGk61,
  kTanhSinh,
  kTrapezoidal,
  kRobust
};

inline std::string ascii_lower(std::string text) {
  for (char &ch : text) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }
  return text;
}

inline BoostQuadratureStrategy quadrature_strategy_from_env() {
  static BoostQuadratureStrategy cached = []() {
    const char *raw = std::getenv("ACCUMULATR_QUADRATURE_STRATEGY");
    if (!raw || *raw == '\0') {
      return BoostQuadratureStrategy::kAuto;
    }
    const std::string val = ascii_lower(std::string(raw));
    if (val == "auto") return BoostQuadratureStrategy::kAuto;
    if (val == "gk21") return BoostQuadratureStrategy::kGk21;
    if (val == "gk61") return BoostQuadratureStrategy::kGk61;
    if (val == "tanh_sinh" || val == "tanhsinh")
      return BoostQuadratureStrategy::kTanhSinh;
    if (val == "trapezoidal" || val == "trap")
      return BoostQuadratureStrategy::kTrapezoidal;
    if (val == "robust")
      return BoostQuadratureStrategy::kRobust;
    return BoostQuadratureStrategy::kAuto;
  }();
  return cached;
}

template <typename Fn>
inline bool try_integrate_gk21(Fn&& integrand, double a, double b, double tol,
                               int max_depth, double &out) {
#if UUBER_HAVE_BOOST_GK
  try {
    out = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
        integrand, a, b, tol, max_depth);
    return std::isfinite(out);
  } catch (...) {
    out = 0.0;
    return false;
  }
#else
  (void)integrand;
  (void)a;
  (void)b;
  (void)tol;
  (void)max_depth;
  out = 0.0;
  return false;
#endif
}

template <typename Fn>
inline bool try_integrate_gk61(Fn&& integrand, double a, double b, double tol,
                               int max_depth, double &out) {
#if UUBER_HAVE_BOOST_GK
  try {
    out = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
        integrand, a, b, tol, max_depth);
    return std::isfinite(out);
  } catch (...) {
    out = 0.0;
    return false;
  }
#else
  (void)integrand;
  (void)a;
  (void)b;
  (void)tol;
  (void)max_depth;
  out = 0.0;
  return false;
#endif
}

template <typename Fn>
inline bool try_integrate_tanh_sinh(Fn&& integrand, double a, double b,
                                    double tol, int max_depth, double &out) {
#if UUBER_HAVE_BOOST_TANH_SINH
  try {
    boost::math::quadrature::tanh_sinh<double> q(
        static_cast<std::size_t>(std::max(8, max_depth + 3)));
    double error_est = 0.0;
    double l1 = 0.0;
    std::size_t levels = 0;
    out = q.integrate(integrand, a, b, tol, &error_est, &l1, &levels);
    return std::isfinite(out);
  } catch (...) {
    out = 0.0;
    return false;
  }
#else
  (void)integrand;
  (void)a;
  (void)b;
  (void)tol;
  (void)max_depth;
  out = 0.0;
  return false;
#endif
}

template <typename Fn>
inline bool try_integrate_trapezoidal(Fn&& integrand, double a, double b,
                                      double tol, int max_depth, double &out) {
#if UUBER_HAVE_BOOST_TRAPEZOIDAL
  try {
    double error_est = 0.0;
    double l1 = 0.0;
    const std::size_t refinements =
        static_cast<std::size_t>(std::max(8, std::min(20, max_depth + 4)));
    out = boost::math::quadrature::trapezoidal(
        integrand, a, b, tol, refinements, &error_est, &l1);
    return std::isfinite(out);
  } catch (...) {
    out = 0.0;
    return false;
  }
#else
  (void)integrand;
  (void)a;
  (void)b;
  (void)tol;
  (void)max_depth;
  out = 0.0;
  return false;
#endif
}

template <typename Fn>
inline double integrate_boost_core(Fn&& integrand, double lower, double upper,
                                   double rel_tol, double abs_tol,
                                   int max_depth) {
  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    Rcpp::stop("integrate_boost_fn requires finite bounds");
  }
  if (max_depth <= 0) max_depth = 10;
  if (rel_tol <= 0.0) rel_tol = 1e-5;
  if (abs_tol <= 0.0) abs_tol = 1e-6;

  double a = lower;
  double b = upper;
  double sign = 1.0;
  if (b < a) {
    std::swap(a, b);
    sign = -1.0;
  }

  const double tol = std::max(std::max(abs_tol, rel_tol), 1e-10);
  const BoostQuadratureStrategy strategy = quadrature_strategy_from_env();
  double integral = 0.0;

  auto robust_auto = [&]() -> bool {
    if (tol >= 1e-4 &&
        try_integrate_gk21(integrand, a, b, tol, max_depth, integral)) {
      return true;
    }
    if (try_integrate_gk61(integrand, a, b, tol, max_depth, integral)) {
      return true;
    }
    if (try_integrate_tanh_sinh(integrand, a, b, tol, max_depth, integral)) {
      return true;
    }
    if (try_integrate_trapezoidal(integrand, a, b, tol, max_depth, integral)) {
      return true;
    }
    return false;
  };

  bool ok = false;
  switch (strategy) {
  case BoostQuadratureStrategy::kGk21:
    ok = try_integrate_gk21(integrand, a, b, tol, max_depth, integral);
    if (!ok) ok = robust_auto();
    break;
  case BoostQuadratureStrategy::kGk61:
    ok = try_integrate_gk61(integrand, a, b, tol, max_depth, integral);
    if (!ok) ok = robust_auto();
    break;
  case BoostQuadratureStrategy::kTanhSinh:
    ok = try_integrate_tanh_sinh(integrand, a, b, tol, max_depth, integral);
    if (!ok) ok = robust_auto();
    break;
  case BoostQuadratureStrategy::kTrapezoidal:
    ok = try_integrate_trapezoidal(integrand, a, b, tol, max_depth, integral);
    if (!ok) ok = robust_auto();
    break;
  case BoostQuadratureStrategy::kRobust:
    ok = robust_auto();
    break;
  case BoostQuadratureStrategy::kAuto:
  default:
    if (tol >= 1e-4) {
      ok = try_integrate_gk21(integrand, a, b, tol, max_depth, integral);
      if (!ok) {
        ok = try_integrate_gk61(integrand, a, b, tol, max_depth, integral);
      }
    } else {
      ok = try_integrate_gk61(integrand, a, b, tol, max_depth, integral);
    }
    if (!ok) {
      ok = robust_auto();
    }
    break;
  }

  if (!ok || !std::isfinite(integral)) {
    integral = 0.0;
  }
  return sign * integral;
}

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
  auto wrapper = [&](double x) -> double {
    if (!std::isfinite(x)) return 0.0;
    Rcpp::NumericVector res = integrand(Rcpp::NumericVector::create(x));
    if (res.size() == 0) return 0.0;
    double val = res[0];
    if (!std::isfinite(val)) return 0.0;
    return val;
  };
  return integrate_boost_core(wrapper, lower, upper, rel_tol, abs_tol,
                              max_depth);
}

template <typename Fn>
inline double integrate_boost_fn(Fn&& integrand,
                                 double lower,
                                 double upper,
                                 double rel_tol,
                                 double abs_tol,
                                 int max_depth) {
  auto wrapper = [&](double x) -> double {
    if (!std::isfinite(x)) return 0.0;
    double val = integrand(x);
    if (!std::isfinite(val)) return 0.0;
    return val;
  };
  return integrate_boost_core(wrapper, lower, upper, rel_tol, abs_tol,
                              max_depth);
}

// Variant using a lower-order Gauss-Kronrod rule (21-point) to reduce per-interval evaluations.

} // namespace uuber
