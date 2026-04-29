#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include "src/eval/quadrature.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace aq = accumulatr::eval::quadrature;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr long double kRefTol = 1e-18L;
constexpr int kRefDepth = 15;

struct Problem {
  std::string name;
  std::string domain;
  std::string description;
  double lower{0.0};
  double upper{0.0};
  std::function<double(double)> fn;
  long double reference{0.0L};
  long double alt_reference{0.0L};
};

struct Method {
  std::string name;
  std::string family;
  int nodes{0};
  std::function<double(const Problem&)> eval;
};

inline double dlnorm_cpp(const double x, const double meanlog, const double sdlog) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 0.0;
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return std::exp(-0.5 * z * z) /
         (x * sdlog * std::sqrt(2.0 * M_PI));
}

inline double plnorm_cpp(const double x, const double meanlog, const double sdlog) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 0.0;
  }
  const double z = (std::log(x) - meanlog) / (sdlog * std::sqrt(2.0));
  return 0.5 * (1.0 + std::erf(z));
}

inline double slnorm_cpp(const double x, const double meanlog, const double sdlog) {
  return 1.0 - plnorm_cpp(x, meanlog, sdlog);
}

template <std::size_t N, typename Fn>
double fixed_finite_integral(const double lower, const double upper, Fn&& fn) {
  const auto rule = aq::map_rule_to_finite_interval<N>(lower, upper);
  return aq::integrate_rule(rule, std::forward<Fn>(fn));
}

template <std::size_t N, typename Fn>
double fixed_tail_shifted_integral(const double lower, Fn&& fn) {
  const auto& rule = aq::gauss_legendre_rule<N>();
  double sum = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    const double u = 0.5 * (rule.nodes[i] + 1.0);
    const double base_weight = 0.5 * rule.weights[i];
    const double one_minus_u = 1.0 - u;
    const double t = lower + u / one_minus_u;
    const double jacobian = 1.0 / (one_minus_u * one_minus_u);
    const double value = fn(t);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    sum += base_weight * jacobian * value;
  }
  return std::isfinite(sum) ? sum : 0.0;
}

template <typename Fn>
double composite_simpson_finite_integral(const double lower,
                                         const double upper,
                                         const int nodes,
                                         Fn&& fn) {
  const int panels = nodes - 1;
  if (!(upper > lower) || nodes < 3 || (panels % 2) != 0) {
    return 0.0;
  }
  const double h = (upper - lower) / static_cast<double>(panels);
  double sum = 0.0;
  for (int i = 0; i < nodes; ++i) {
    const double x = lower + h * static_cast<double>(i);
    const double value = fn(x);
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    const double weight =
        (i == 0 || i == panels) ? 1.0 : ((i % 2) != 0 ? 4.0 : 2.0);
    sum += weight * value;
  }
  const double result = sum * h / 3.0;
  return std::isfinite(result) ? result : 0.0;
}

template <typename Fn>
double composite_simpson_tail_shifted_integral(const double lower,
                                               const int nodes,
                                               Fn&& fn) {
  const int panels = nodes - 1;
  if (nodes < 3 || (panels % 2) != 0) {
    return 0.0;
  }
  const double h = 1.0 / static_cast<double>(panels);
  double sum = 0.0;
  for (int i = 0; i < nodes; ++i) {
    const double u = h * static_cast<double>(i);
    double value = 0.0;
    if (u < 1.0) {
      const double one_minus_u = 1.0 - u;
      const double t = lower + u / one_minus_u;
      const double jacobian = 1.0 / (one_minus_u * one_minus_u);
      value = fn(t) * jacobian;
    }
    if (!std::isfinite(value) || value == 0.0) {
      continue;
    }
    const double weight =
        (i == 0 || i == panels) ? 1.0 : ((i % 2) != 0 ? 4.0 : 2.0);
    sum += weight * value;
  }
  const double result = sum * h / 3.0;
  return std::isfinite(result) ? result : 0.0;
}

long double reference_finite_gk(const Problem& problem) {
  return boost::math::quadrature::gauss_kronrod<long double, 61>::integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.lower),
      static_cast<long double>(problem.upper),
      kRefDepth,
      kRefTol);
}

long double reference_finite_tanh(const Problem& problem) {
  boost::math::quadrature::tanh_sinh<long double> integrator(15);
  return integrator.integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.lower),
      static_cast<long double>(problem.upper),
      kRefTol);
}

long double reference_tail_exp(const Problem& problem) {
  boost::math::quadrature::exp_sinh<long double> integrator(15);
  return integrator.integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.lower),
      std::numeric_limits<long double>::infinity(),
      kRefTol);
}

long double reference_tail_transformed_gk(const Problem& problem) {
  return boost::math::quadrature::gauss_kronrod<long double, 61>::integrate(
      [&](const long double u) {
        if (u <= 0.0L) {
          return static_cast<long double>(
              problem.fn(static_cast<double>(problem.lower)));
        }
        if (u >= 1.0L) {
          return 0.0L;
        }
        const long double one_minus_u = 1.0L - u;
        const long double t = static_cast<long double>(problem.lower) +
                              u / one_minus_u;
        const long double jacobian = 1.0L / (one_minus_u * one_minus_u);
        return static_cast<long double>(problem.fn(static_cast<double>(t))) *
               jacobian;
      },
      0.0L,
      1.0L,
      kRefDepth,
      kRefTol);
}

std::vector<Problem> build_problems() {
  std::vector<Problem> problems;

  {
    const double a_m = std::log(0.30), a_s = 0.16;
    const double b_m = std::log(0.22), b_s = 0.16;
    const double t = 0.60;
    problems.push_back(Problem{
        "finite_convolution_pdf",
        "finite",
        "after/onset convolution density",
        0.0,
        t,
        [=](const double u) {
          return dlnorm_cpp(u, a_m, a_s) * dlnorm_cpp(t - u, b_m, b_s);
        }});
  }

  {
    const double a_m = std::log(0.28), a_s = 0.18;
    const double b_m = std::log(0.32), b_s = 0.18;
    const double c_m = std::log(0.30), c_s = 0.18;
    const double t = 0.60;
    problems.push_back(Problem{
        "finite_guard_truth",
        "finite",
        "density times two survival gates",
        0.0,
        t,
        [=](const double u) {
          return dlnorm_cpp(u, a_m, a_s) *
                 slnorm_cpp(u, b_m, b_s) *
                 slnorm_cpp(u, c_m, c_s);
        }});
  }

  {
    const double a_m = std::log(0.28), a_s = 0.18;
    const double b_m = std::log(0.31), b_s = 0.19;
    const double c_m = std::log(0.34), c_s = 0.17;
    const double t = 0.58;
    problems.push_back(Problem{
        "finite_all_of_three_way_density",
        "finite",
        "three-way all_of derivative terms at fixed time",
        0.0,
        t,
        [=](const double u) {
          return dlnorm_cpp(u, a_m, a_s) *
                     plnorm_cpp(u, b_m, b_s) *
                     plnorm_cpp(u, c_m, c_s) +
                 dlnorm_cpp(u, b_m, b_s) *
                     plnorm_cpp(u, a_m, a_s) *
                     plnorm_cpp(u, c_m, c_s) +
                 dlnorm_cpp(u, c_m, c_s) *
                     plnorm_cpp(u, a_m, a_s) *
                     plnorm_cpp(u, b_m, b_s);
        }});
  }

  {
    const double go_left_m = std::log(0.30), go_left_s = 0.20;
    const double go_right_m = std::log(0.32), go_right_s = 0.20;
    const double stop_m = std::log(0.15), stop_s = 0.18;
    const double stop_onset = 0.15;
    problems.push_back(Problem{
        "tail_stop_like_na_mass",
        "tail",
        "shifted stop density gated by go survivals",
        stop_onset,
        kInf,
        [=](const double t) {
          return dlnorm_cpp(t - stop_onset, stop_m, stop_s) *
                 slnorm_cpp(t, go_left_m, go_left_s) *
                 slnorm_cpp(t, go_right_m, go_right_s);
        }});
  }

  {
    const double a_m = std::log(0.30), a_s = 0.15;
    const double b_m = std::log(0.30), b_s = 0.20;
    const double c_m = std::log(0.30), c_s = 0.18;
    problems.push_back(Problem{
        "tail_pooled_total_mass",
        "tail",
        "pooled two-way branch against competitor",
        0.0,
        kInf,
        [=](const double t) {
          const double pool_pdf =
              dlnorm_cpp(t, a_m, a_s) * slnorm_cpp(t, b_m, b_s) +
              dlnorm_cpp(t, b_m, b_s) * slnorm_cpp(t, a_m, a_s);
          const double pool_surv =
              slnorm_cpp(t, a_m, a_s) * slnorm_cpp(t, b_m, b_s);
          const double c_pdf = dlnorm_cpp(t, c_m, c_s);
          const double c_surv = slnorm_cpp(t, c_m, c_s);
          return pool_pdf * c_surv + c_pdf * pool_surv;
        }});
  }

  {
    const double a_m = std::log(0.28), a_s = 0.18;
    const double c_m = std::log(0.30), c_s = 0.18;
    problems.push_back(Problem{
        "tail_all_of_probability",
        "tail",
        "two-child all_of probability mass",
        0.0,
        kInf,
        [=](const double t) {
          return dlnorm_cpp(t, a_m, a_s) * plnorm_cpp(t, c_m, c_s) +
                 dlnorm_cpp(t, c_m, c_s) * plnorm_cpp(t, a_m, a_s);
        }});
  }

  for (auto& problem : problems) {
    if (problem.domain == "tail") {
      problem.reference = reference_tail_exp(problem);
      problem.alt_reference = reference_tail_transformed_gk(problem);
    } else {
      problem.reference = reference_finite_gk(problem);
      problem.alt_reference = reference_finite_tanh(problem);
    }
  }

  return problems;
}

template <std::size_t N>
void add_gauss_legendre_method(std::vector<Method>& methods) {
  methods.push_back(Method{
      "gl" + std::to_string(N),
      "fixed_gauss_legendre",
      static_cast<int>(N),
      [](const Problem& problem) {
        if (problem.domain == "finite") {
          return fixed_finite_integral<N>(problem.lower, problem.upper, problem.fn);
        }
        return fixed_tail_shifted_integral<N>(problem.lower, problem.fn);
      }});
}

std::vector<Method> build_methods() {
  std::vector<Method> methods;

  add_gauss_legendre_method<7>(methods);
  add_gauss_legendre_method<11>(methods);
  add_gauss_legendre_method<15>(methods);
  add_gauss_legendre_method<21>(methods);
  add_gauss_legendre_method<31>(methods);
  add_gauss_legendre_method<47>(methods);
  add_gauss_legendre_method<63>(methods);
  add_gauss_legendre_method<95>(methods);
  add_gauss_legendre_method<127>(methods);
  add_gauss_legendre_method<191>(methods);
  add_gauss_legendre_method<255>(methods);

  for (const int nodes : {7, 11, 15, 21, 31, 47, 63, 95, 127, 191, 255}) {
    methods.push_back(Method{
        "simpson" + std::to_string(nodes),
        "composite_simpson",
        nodes,
        [nodes](const Problem& problem) {
          if (problem.domain == "finite") {
            return composite_simpson_finite_integral(
                problem.lower, problem.upper, nodes, problem.fn);
          }
          return composite_simpson_tail_shifted_integral(
              problem.lower, nodes, problem.fn);
        }});
  }

  return methods;
}

int choose_repetitions(const std::function<double()>& fn) {
  int reps = 1;
  while (reps < 1 << 20) {
    volatile double sink = 0.0;
    const auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < reps; ++i) {
      sink += fn() * 1e-300;
    }
    const auto end = std::chrono::high_resolution_clock::now();
    const double elapsed =
        std::chrono::duration<double>(end - start).count();
    if (elapsed >= 0.04) {
      return reps;
    }
    reps *= 2;
  }
  return reps;
}

bool is_engine_default(const Problem& problem, const Method& method) {
  if (method.family != "fixed_gauss_legendre") {
    return false;
  }
  return (problem.domain == "finite" && method.nodes == 31) ||
         (problem.domain == "tail" && method.nodes == 47);
}

int function_evals_for_row(const Problem& problem, const Method& method) {
  if (method.family == "composite_simpson" && problem.domain == "tail") {
    return method.nodes - 1;
  }
  return method.nodes;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "Usage: %s <output.csv>\n", argv[0]);
    return 1;
  }

  const auto problems = build_problems();
  const auto methods = build_methods();

  std::ofstream out(argv[1]);
  if (!out) {
    std::fprintf(stderr, "Failed to open output file\n");
    return 1;
  }

  out << "problem,domain,description,method_family,method,nodes,function_evals,"
         "engine_default,reference,alt_reference,reference_disagreement,value,"
         "abs_prob_err,rel_prob_err,loglik_err,usec_per_eval,repetitions\n";
  out << std::setprecision(17);

  for (const auto& problem : problems) {
    for (const auto& method : methods) {
      const auto eval_once = [&]() {
        return method.eval(problem);
      };
      const double value = eval_once();
      const int reps = choose_repetitions(eval_once);

      volatile double sink = 0.0;
      const auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < reps; ++i) {
        sink += eval_once() * 1e-300;
      }
      const auto end = std::chrono::high_resolution_clock::now();
      (void)sink;

      const double usec_per_eval =
          std::chrono::duration<double, std::micro>(end - start).count() /
          static_cast<double>(reps);
      const long double abs_prob_err =
          std::fabs(static_cast<long double>(value) - problem.reference);
      const long double rel_prob_err =
          abs_prob_err / std::max(std::fabs(problem.reference), 1e-30L);
      const double loglik_err =
          (value > 0.0 && problem.reference > 0.0L)
              ? std::fabs(std::log(value) -
                          std::log(static_cast<double>(problem.reference)))
              : std::numeric_limits<double>::infinity();
      const long double ref_disagreement =
          std::fabs(problem.reference - problem.alt_reference);

      out << problem.name << ','
          << problem.domain << ','
          << '"' << problem.description << '"' << ','
          << method.family << ','
          << method.name << ','
          << method.nodes << ','
          << function_evals_for_row(problem, method) << ','
          << (is_engine_default(problem, method) ? "true" : "false") << ','
          << static_cast<double>(problem.reference) << ','
          << static_cast<double>(problem.alt_reference) << ','
          << static_cast<double>(ref_disagreement) << ','
          << value << ','
          << static_cast<double>(abs_prob_err) << ','
          << static_cast<double>(rel_prob_err) << ','
          << loglik_err << ','
          << usec_per_eval << ','
          << reps << '\n';
    }
  }

  return 0;
}
