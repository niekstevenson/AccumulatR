#pragma once

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace accumulatr::eval::quadrature {

constexpr std::size_t kDefaultFiniteOrder = 31;
constexpr std::size_t kDefaultTailOrder = 63;
constexpr unsigned kRobustFiniteOrder = 61;
constexpr unsigned kRobustMaxDepth = 12;
constexpr double kRobustTolerance = 1e-4;

template <std::size_t N>
struct MappedRule {
  std::array<double, N> nodes{};
  std::array<double, N> weights{};
};

struct FiniteBatch {
  double lower{0.0};
  double upper{0.0};
  MappedRule<kDefaultFiniteOrder> nodes;
};

struct TailBatch {
  MappedRule<kDefaultTailOrder> nodes;
};

template <std::size_t N>
struct Rule {
  std::array<double, N> nodes{};
  std::array<double, N> weights{};
};

template <std::size_t N>
inline Rule<N> build_gauss_legendre_rule() {
  Rule<N> rule;
  constexpr double kEps = 1e-14;
  constexpr double kPi = 3.141592653589793238462643383279502884;
  constexpr std::size_t kHalf = (N + 1U) / 2U;
  for (std::size_t i = 0; i < kHalf; ++i) {
    double z = std::cos(kPi * (static_cast<double>(i) + 0.75) /
                        (static_cast<double>(N) + 0.5));
    double z_prev = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
    double pp = 0.0;
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (std::size_t j = 1; j <= N; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * static_cast<double>(j) - 1.0) * z * p2 -
              (static_cast<double>(j) - 1.0) * p3) /
             static_cast<double>(j);
      }
      pp = static_cast<double>(N) * (z * p1 - p2) / (z * z - 1.0);
      z_prev = z;
      z = z_prev - p1 / pp;
    } while (std::fabs(z - z_prev) > kEps);

    const double w = 2.0 / ((1.0 - z * z) * pp * pp);
    rule.nodes[i] = -z;
    rule.nodes[N - 1U - i] = z;
    rule.weights[i] = w;
    rule.weights[N - 1U - i] = w;
  }
  return rule;
}

template <std::size_t N>
inline const Rule<N> &gauss_legendre_rule() {
  static const Rule<N> rule = build_gauss_legendre_rule<N>();
  return rule;
}

template <std::size_t N>
inline MappedRule<N> map_rule_to_finite_interval(const double lower,
                                                 const double upper) {
  MappedRule<N> mapped;
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return mapped;
  }
  const auto &rule = gauss_legendre_rule<N>();
  const double scale = 0.5 * (upper - lower);
  const double shift = 0.5 * (upper + lower);
  for (std::size_t i = 0; i < N; ++i) {
    mapped.nodes[i] = shift + scale * rule.nodes[i];
    mapped.weights[i] = scale * rule.weights[i];
  }
  return mapped;
}

template <std::size_t N>
inline MappedRule<N> map_rule_to_positive_tail() {
  MappedRule<N> mapped;
  const auto &rule = gauss_legendre_rule<N>();
  for (std::size_t i = 0; i < N; ++i) {
    const double u = 0.5 * (rule.nodes[i] + 1.0);
    const double base_weight = 0.5 * rule.weights[i];
    const double one_minus_u = 1.0 - u;
    const double t = u / one_minus_u;
    const double jacobian = 1.0 / (one_minus_u * one_minus_u);
    mapped.nodes[i] = t;
    mapped.weights[i] = base_weight * jacobian;
  }
  return mapped;
}

inline FiniteBatch build_finite_batch(const double lower, const double upper) {
  FiniteBatch batch;
  batch.lower = lower;
  batch.upper = upper;
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return batch;
  }
  batch.nodes = map_rule_to_finite_interval<kDefaultFiniteOrder>(lower, upper);
  return batch;
}

inline const TailBatch &canonical_tail_batch() {
  static const TailBatch batch{map_rule_to_positive_tail<kDefaultTailOrder>()};
  return batch;
}

template <std::size_t N, typename Fn>
inline double integrate_rule(const MappedRule<N> &rule, Fn &&fn) {
  double sum = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    const double value = fn(rule.nodes[i]);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    sum += rule.weights[i] * value;
  }
  return std::isfinite(sum) ? sum : 0.0;
}

template <typename Fn>
inline double integrate_finite_default(const FiniteBatch &batch, Fn &&fn) {
  return integrate_rule(batch.nodes, std::forward<Fn>(fn));
}

template <typename Fn>
inline double integrate_finite_default(const double lower,
                                       const double upper,
                                       Fn &&fn) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return 0.0;
  }
  const auto batch = build_finite_batch(lower, upper);
  return integrate_finite_default(batch, std::forward<Fn>(fn));
}

template <typename Fn>
inline double integrate_finite_default(Fn &&fn,
                                       const double lower,
                                       const double upper) {
  return integrate_finite_default(lower, upper, std::forward<Fn>(fn));
}

template <typename Fn>
inline double integrate_finite_robust(const double lower,
                                      const double upper,
                                      Fn &&fn,
                                      const double tol = kRobustTolerance) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return 0.0;
  }
  double error = 0.0;
  double l1 = 0.0;
  const double value =
      boost::math::quadrature::gauss_kronrod<double, kRobustFiniteOrder>::integrate(
          std::forward<Fn>(fn),
          lower,
          upper,
          kRobustMaxDepth,
          tol,
          &error,
          &l1);
  return std::isfinite(value) ? value : 0.0;
}

template <typename Fn>
inline double integrate_finite_robust(Fn &&fn,
                                      const double lower,
                                      const double upper,
                                      const double tol = kRobustTolerance) {
  return integrate_finite_robust(
      lower, upper, std::forward<Fn>(fn), tol);
}

template <typename Fn>
inline double integrate_tail_default(Fn &&fn) {
  return integrate_rule(canonical_tail_batch().nodes, std::forward<Fn>(fn));
}

template <typename Fn>
inline double integrate_tail_robust(Fn &&fn,
                                    const double lower = 0.0,
                                    const double tol = kRobustTolerance) {
  if (!std::isfinite(lower)) {
    return 0.0;
  }
  boost::math::quadrature::exp_sinh<double> integrator(12);
  double error = 0.0;
  double l1 = 0.0;
  std::size_t levels = 0;
  const double value = integrator.integrate(
      std::forward<Fn>(fn),
      lower,
      std::numeric_limits<double>::infinity(),
      tol,
      &error,
      &l1,
      &levels);
  return std::isfinite(value) ? value : 0.0;
}

} // namespace accumulatr::eval::quadrature
