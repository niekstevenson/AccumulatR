#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include "src/eval/quadrature.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace aq = accumulatr::eval::quadrature;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr long double kRefTol = 1e-18L;
constexpr int kRefDepth = 15;

struct Problem {
  std::string family;
  std::string model;
  std::string integral_type;
  std::string domain;
  double lower{0.0};
  double upper{0.0};
  std::function<double(double)> fn;
  long double reference{0.0L};
  long double alt_reference{0.0L};
};

struct Method {
  std::string name;
  std::string family;
  std::optional<int> order;
  std::optional<double> tol;
  bool supports_finite{true};
  bool supports_tail{true};
  std::function<double(const Problem&, double*, double*, std::size_t*)> eval;
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
double fixed_finite_integral(const double a, const double b, Fn&& fn) {
  const auto rule = aq::map_rule_to_finite_interval<N>(a, b);
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

template <unsigned N>
double gk_finite_integral(const Problem& problem,
                          const double tol,
                          double* error,
                          double* l1) {
  return boost::math::quadrature::gauss_kronrod<double, N>::integrate(
      problem.fn,
      problem.lower,
      problem.upper,
      12,
      tol,
      error,
      l1);
}

template <unsigned N>
double gk_tail_shifted_integral(const Problem& problem,
                                const double tol,
                                double* error,
                                double* l1) {
  auto transformed = [&](const double u) {
    if (u <= 0.0) {
      return problem.fn(problem.lower);
    }
    if (u >= 1.0) {
      return 0.0;
    }
    const double one_minus_u = 1.0 - u;
    const double t = problem.lower + u / one_minus_u;
    const double jacobian = 1.0 / (one_minus_u * one_minus_u);
    return problem.fn(t) * jacobian;
  };
  return boost::math::quadrature::gauss_kronrod<double, N>::integrate(
      transformed, 0.0, 1.0, 12, tol, error, l1);
}

double tanh_sinh_finite_integral(const Problem& problem,
                                 const double tol,
                                 double* error,
                                 double* l1,
                                 std::size_t* levels) {
  boost::math::quadrature::tanh_sinh<double> integrator(15);
  return integrator.integrate(
      problem.fn, problem.lower, problem.upper, tol, error, l1, levels);
}

double exp_sinh_tail_integral(const Problem& problem,
                              const double tol,
                              double* error,
                              double* l1,
                              std::size_t* levels) {
  boost::math::quadrature::exp_sinh<double> integrator(12);
  return integrator.integrate(
      problem.fn,
      problem.lower,
      std::numeric_limits<double>::infinity(),
      tol,
      error,
      l1,
      levels);
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

std::vector<Problem> build_problems() {
  std::vector<Problem> problems;

  {
    // example_3_stop_na: mapped-NA stop outcome mass
    const double go_left_m = std::log(0.30);
    const double go_left_s = 0.20;
    const double go_right_m = std::log(0.32);
    const double go_right_s = 0.20;
    const double stop_m = std::log(0.15);
    const double stop_s = 0.18;
    const double stop_onset = 0.15;
    problems.push_back(Problem{
        "direct_tail_outcome_mass",
        "example_3_stop_na",
        "mapped-NA outcome tail mass",
        "tail",
        stop_onset,
        kInf,
        [=](const double t) {
          return dlnorm_cpp(t - stop_onset, stop_m, stop_s) *
                 slnorm_cpp(t, go_left_m, go_left_s) *
                 slnorm_cpp(t, go_right_m, go_right_s);
        }});
  }

  {
    // example_4_two_on_one: direct total finite mass with overlapping pool outcome
    const double a_m = std::log(0.30), a_s = 0.15;
    const double b_m = std::log(0.30), b_s = 0.20;
    const double c_m = std::log(0.30), c_s = 0.18;
    problems.push_back(Problem{
        "direct_tail_total_mass",
        "example_4_two_on_one",
        "total finite mass with pooled top-1 outcome",
        "tail",
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
    // exact onset convolution: b starts after a
    const double a_m = std::log(0.30), a_s = 0.16;
    const double b_m = std::log(0.22), b_s = 0.16;
    const double lag = 0.0;
    const double t = 0.60;
    problems.push_back(Problem{
        "exact_finite_convolution",
        "example_23_ranked_chain",
        "onset convolution pdf",
        "finite",
        0.0,
        t - lag,
        [=](const double u) {
          return dlnorm_cpp(u, a_m, a_s) *
                 dlnorm_cpp(t - u - lag, b_m, b_s);
        }});
  }

  {
    // exact guard/truth: guarded A against blockers B and C
    const double a_m = std::log(0.28), a_s = 0.18;
    const double b_m = std::log(0.32), b_s = 0.18;
    const double c_m = std::log(0.30), c_s = 0.18;
    const double t = 0.60;
    problems.push_back(Problem{
        "exact_finite_truth",
        "guarded_exclusion_typical",
        "guard/truth cdf integral",
        "finite",
        0.0,
        t,
        [=](const double u) {
          return dlnorm_cpp(u, a_m, a_s) *
                 slnorm_cpp(u, b_m, b_s) *
                 slnorm_cpp(u, c_m, c_s);
        }});
  }

  {
    // exact tail probability: single conjunctive all_of outcome
    const double a_m = std::log(0.28), a_s = 0.18;
    const double c_m = std::log(0.30), c_s = 0.18;
    problems.push_back(Problem{
        "exact_tail_probability",
        "single_all_of_tail",
        "all_of outcome tail probability",
        "tail",
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
      problem.alt_reference = problem.reference;
    } else {
      problem.reference = reference_finite_gk(problem);
      problem.alt_reference = reference_finite_tanh(problem);
    }
  }

  return problems;
}

std::vector<Method> build_methods() {
  std::vector<Method> methods;

  methods.push_back(Method{
      "fixed_gl31",
      "fixed",
      31,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "finite") {
          return fixed_finite_integral<31>(problem.lower, problem.upper, problem.fn);
        }
        return fixed_tail_shifted_integral<31>(problem.lower, problem.fn);
      }});

  methods.push_back(Method{
      "fixed_gl47",
      "fixed",
      47,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "finite") {
          return fixed_finite_integral<47>(problem.lower, problem.upper, problem.fn);
        }
        return fixed_tail_shifted_integral<47>(problem.lower, problem.fn);
      }});

  methods.push_back(Method{
      "fixed_gl63",
      "fixed",
      63,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "finite") {
          return fixed_finite_integral<63>(problem.lower, problem.upper, problem.fn);
        }
        return fixed_tail_shifted_integral<63>(problem.lower, problem.fn);
      }});

  methods.push_back(Method{
      "gk21_tol1e-4",
      "gauss_kronrod",
      21,
      1e-4,
      true,
      true,
      [](const Problem& problem, double* error, double* l1, std::size_t*) {
        if (problem.domain == "finite") {
          return gk_finite_integral<21>(problem, 1e-4, error, l1);
        }
        return gk_tail_shifted_integral<21>(problem, 1e-4, error, l1);
      }});

  methods.push_back(Method{
      "gk21_tol1e-6",
      "gauss_kronrod",
      21,
      1e-6,
      true,
      true,
      [](const Problem& problem, double* error, double* l1, std::size_t*) {
        if (problem.domain == "finite") {
          return gk_finite_integral<21>(problem, 1e-6, error, l1);
        }
        return gk_tail_shifted_integral<21>(problem, 1e-6, error, l1);
      }});

  methods.push_back(Method{
      "gk61_tol1e-4",
      "gauss_kronrod",
      61,
      1e-4,
      true,
      true,
      [](const Problem& problem, double* error, double* l1, std::size_t*) {
        if (problem.domain == "finite") {
          return gk_finite_integral<61>(problem, 1e-4, error, l1);
        }
        return gk_tail_shifted_integral<61>(problem, 1e-4, error, l1);
      }});

  methods.push_back(Method{
      "tanh_sinh_tol1e-4",
      "double_exponential",
      std::nullopt,
      1e-4,
      true,
      false,
      [](const Problem& problem, double* error, double* l1, std::size_t* levels) {
        return tanh_sinh_finite_integral(problem, 1e-4, error, l1, levels);
      }});

  methods.push_back(Method{
      "tanh_sinh_tol1e-6",
      "double_exponential",
      std::nullopt,
      1e-6,
      true,
      false,
      [](const Problem& problem, double* error, double* l1, std::size_t* levels) {
        return tanh_sinh_finite_integral(problem, 1e-6, error, l1, levels);
      }});

  methods.push_back(Method{
      "exp_sinh_tol1e-4",
      "double_exponential",
      std::nullopt,
      1e-4,
      false,
      true,
      [](const Problem& problem, double* error, double* l1, std::size_t* levels) {
        return exp_sinh_tail_integral(problem, 1e-4, error, l1, levels);
      }});

  methods.push_back(Method{
      "exp_sinh_tol1e-6",
      "double_exponential",
      std::nullopt,
      1e-6,
      false,
      true,
      [](const Problem& problem, double* error, double* l1, std::size_t* levels) {
        return exp_sinh_tail_integral(problem, 1e-6, error, l1, levels);
      }});

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

std::string format_optional_int(const std::optional<int>& value) {
  if (!value.has_value()) {
    return "";
  }
  return std::to_string(*value);
}

std::string format_optional_double(const std::optional<double>& value) {
  if (!value.has_value()) {
    return "";
  }
  std::ostringstream oss;
  oss << std::setprecision(8) << *value;
  return oss.str();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "Usage: %s <output.csv>\\n", argv[0]);
    return 1;
  }

  const auto problems = build_problems();
  const auto methods = build_methods();

  std::ofstream out(argv[1]);
  if (!out) {
    std::fprintf(stderr, "Failed to open output file\\n");
    return 1;
  }

  out << "family,model,integral_type,domain,method,method_family,order,tol,"
         "reference,alt_reference,reference_disagreement,value,abs_prob_err,"
         "rel_prob_err,loglik_err,usec_per_eval,repetitions\\n";
  out << std::setprecision(17);

  for (const auto& problem : problems) {
    const long double reference_disagreement =
        std::fabs(problem.reference - problem.alt_reference);
    for (const auto& method : methods) {
      if (problem.domain == "finite" && !method.supports_finite) {
        continue;
      }
      if (problem.domain == "tail" && !method.supports_tail) {
        continue;
      }

      const auto eval_once = [&]() {
        double error = 0.0;
        double l1 = 0.0;
        std::size_t levels = 0;
        return method.eval(problem, &error, &l1, &levels);
      };
      const int reps = choose_repetitions(eval_once);
      volatile double sink = 0.0;
      double value = 0.0;
      const auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < reps; ++i) {
        value = eval_once();
        sink += value * 1e-300;
      }
      const auto end = std::chrono::high_resolution_clock::now();
      (void)sink;

      const double usec_per_eval =
          std::chrono::duration<double, std::micro>(end - start).count() /
          static_cast<double>(reps);
      const double abs_prob_err =
          std::fabs(static_cast<long double>(value) - problem.reference);
      const double rel_prob_err =
          abs_prob_err / std::fabs(static_cast<double>(problem.reference));
      const double loglik_err =
          std::fabs(std::log(value) - std::log(static_cast<double>(problem.reference)));

      out << problem.family << ','
          << problem.model << ','
          << problem.integral_type << ','
          << problem.domain << ','
          << method.name << ','
          << method.family << ','
          << format_optional_int(method.order) << ','
          << format_optional_double(method.tol) << ','
          << static_cast<double>(problem.reference) << ','
          << static_cast<double>(problem.alt_reference) << ','
          << static_cast<double>(reference_disagreement) << ','
          << value << ','
          << abs_prob_err << ','
          << rel_prob_err << ','
          << loglik_err << ','
          << usec_per_eval << ','
          << reps << '\n';
    }
  }

  return 0;
}
