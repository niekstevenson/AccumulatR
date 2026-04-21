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
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace aq = accumulatr::eval::quadrature;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr long double kRefTol = 1e-18L;
constexpr int kRefDepth = 15;

struct Problem {
  std::string name;
  std::string domain;
  double a;
  double b;
  std::function<double(double)> fn;
  long double reference{0.0L};
  long double ref_alt{0.0L};
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

inline double dlnorm(const double x, const double meanlog, const double sdlog) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 0.0;
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return std::exp(-0.5 * z * z) /
         (x * sdlog * std::sqrt(2.0 * M_PI));
}

inline double plnorm(const double x, const double meanlog, const double sdlog) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    return 0.0;
  }
  const double z = (std::log(x) - meanlog) / (sdlog * std::sqrt(2.0));
  return 0.5 * (1.0 + std::erf(z));
}

inline double slnorm(const double x, const double meanlog, const double sdlog) {
  return 1.0 - plnorm(x, meanlog, sdlog);
}

inline long double dlnorm_ld(const long double x,
                             const long double meanlog,
                             const long double sdlog) {
  if (!(x > 0.0L) || !std::isfinite(static_cast<double>(x))) {
    return 0.0L;
  }
  const long double z = (std::log(x) - meanlog) / sdlog;
  return std::exp(-0.5L * z * z) /
         (x * sdlog * std::sqrt(2.0L * static_cast<long double>(M_PI)));
}

inline long double plnorm_ld(const long double x,
                             const long double meanlog,
                             const long double sdlog) {
  if (!(x > 0.0L) || !std::isfinite(static_cast<double>(x))) {
    return 0.0L;
  }
  const long double z =
      (std::log(x) - meanlog) / (sdlog * std::sqrt(2.0L));
  return 0.5L * (1.0L + std::erf(z));
}

inline long double slnorm_ld(const long double x,
                             const long double meanlog,
                             const long double sdlog) {
  return 1.0L - plnorm_ld(x, meanlog, sdlog);
}

template <std::size_t N, typename Fn>
double fixed_finite_integral(const double a, const double b, Fn&& fn) {
  const auto rule = aq::map_rule_to_finite_interval<N>(a, b);
  return aq::integrate_rule(rule, std::forward<Fn>(fn));
}

template <std::size_t N, typename Fn>
double fixed_tail_integral(Fn&& fn) {
  const auto rule = aq::map_rule_to_positive_tail<N>();
  return aq::integrate_rule(rule, std::forward<Fn>(fn));
}

template <unsigned N>
double gk_integral(const Problem& problem,
                   const double tol,
                   double* error,
                   double* L1) {
  return boost::math::quadrature::gauss_kronrod<double, N>::integrate(
      problem.fn,
      problem.a,
      problem.b,
      12,
      tol,
      error,
      L1);
}

double tanh_sinh_integral(const Problem& problem,
                          const double tol,
                          double* error,
                          double* L1,
                          std::size_t* levels) {
  boost::math::quadrature::tanh_sinh<double> integrator(15);
  return integrator.integrate(problem.fn,
                              problem.a,
                              problem.b,
                              tol,
                              error,
                              L1,
                              levels);
}

double exp_sinh_integral(const Problem& problem,
                         const double tol,
                         double* error,
                         double* L1,
                         std::size_t* levels) {
  boost::math::quadrature::exp_sinh<double> integrator(12);
  return integrator.integrate(problem.fn,
                              problem.a,
                              problem.b,
                              tol,
                              error,
                              L1,
                              levels);
}

long double reference_tail(const Problem& problem) {
  boost::math::quadrature::exp_sinh<long double> integrator(15);
  return integrator.integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.a),
      static_cast<long double>(problem.b),
      kRefTol);
}

long double reference_finite_gk(const Problem& problem) {
  return boost::math::quadrature::gauss_kronrod<long double, 61>::integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.a),
      static_cast<long double>(problem.b),
      kRefDepth,
      kRefTol);
}

long double reference_finite_tanh(const Problem& problem) {
  boost::math::quadrature::tanh_sinh<long double> integrator(15);
  return integrator.integrate(
      [&](const long double x) {
        return static_cast<long double>(problem.fn(static_cast<double>(x)));
      },
      static_cast<long double>(problem.a),
      static_cast<long double>(problem.b),
      kRefTol);
}

std::vector<Problem> build_problems() {
  const double mu1 = std::log(0.30);
  const double mu2 = std::log(0.32);
  const double sigma = 0.18;
  const double q = 0.10;
  const double rt = 0.60;

  std::vector<Problem> problems;

  problems.push_back(Problem{
      "tail_single_lognorm",
      "tail",
      0.0,
      kInf,
      [=](const double t) {
        return (1.0 - q) * dlnorm(t, mu1, sigma);
      }});

  problems.push_back(Problem{
      "tail_shared_race",
      "tail",
      0.0,
      kInf,
      [=](const double t) {
        return (1.0 - q) *
               (dlnorm(t, mu1, sigma) * slnorm(t, mu2, sigma) +
                dlnorm(t, mu2, sigma) * slnorm(t, mu1, sigma));
      }});

  problems.push_back(Problem{
      "tail_indep_race",
      "tail",
      0.0,
      kInf,
      [=](const double t) {
        return (1.0 - q) * dlnorm(t, mu1, sigma) *
                   (q + (1.0 - q) * slnorm(t, mu2, sigma)) +
               (1.0 - q) * dlnorm(t, mu2, sigma) *
                   (q + (1.0 - q) * slnorm(t, mu1, sigma));
      }});

  problems.push_back(Problem{
      "finite_race_mass_to_rt",
      "finite",
      0.0,
      rt,
      [=](const double t) {
        return dlnorm(t, mu1, sigma) * slnorm(t, mu2, sigma) +
               dlnorm(t, mu2, sigma) * slnorm(t, mu1, sigma);
      }});

  problems.push_back(Problem{
      "finite_conv_pdf",
      "finite",
      0.0,
      rt,
      [=](const double s) {
        return dlnorm(s, mu1, sigma) * dlnorm(rt - s, mu2, sigma);
      }});

  problems.push_back(Problem{
      "finite_conv_cdf",
      "finite",
      0.0,
      rt,
      [=](const double s) {
        return dlnorm(s, mu1, sigma) * plnorm(rt - s, mu2, sigma);
      }});

  for (auto& problem : problems) {
    if (problem.domain == "tail") {
      problem.reference = reference_tail(problem);
      problem.ref_alt =
          boost::math::quadrature::gauss_kronrod<long double, 61>::integrate(
              [&](const long double x) {
                return static_cast<long double>(
                    problem.fn(static_cast<double>(x)));
              },
              0.0L,
              std::numeric_limits<long double>::infinity(),
              kRefDepth,
              kRefTol);
    } else {
      problem.reference = reference_finite_gk(problem);
      problem.ref_alt = reference_finite_tanh(problem);
    }
  }

  return problems;
}

std::vector<Method> build_methods() {
  std::vector<Method> methods;

  methods.push_back(Method{
      "fixed_gl15",
      "fixed_gauss_legendre",
      15,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "tail") {
          return fixed_tail_integral<15>(problem.fn);
        }
        return fixed_finite_integral<15>(problem.a, problem.b, problem.fn);
      }});

  methods.push_back(Method{
      "fixed_gl31",
      "fixed_gauss_legendre",
      31,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "tail") {
          return fixed_tail_integral<31>(problem.fn);
        }
        return fixed_finite_integral<31>(problem.a, problem.b, problem.fn);
      }});

  methods.push_back(Method{
      "fixed_gl63",
      "fixed_gauss_legendre",
      63,
      std::nullopt,
      true,
      true,
      [](const Problem& problem, double*, double*, std::size_t*) {
        if (problem.domain == "tail") {
          return fixed_tail_integral<63>(problem.fn);
        }
        return fixed_finite_integral<63>(problem.a, problem.b, problem.fn);
      }});

  for (double tol : {1e-6, 1e-8, 1e-10}) {
    {
      std::ostringstream name;
      name << "boost_gk21_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_gauss_kronrod",
          21,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t*) {
            return gk_integral<21>(problem, tol, error, L1);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_gk15_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_gauss_kronrod",
          15,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t*) {
            return gk_integral<15>(problem, tol, error, L1);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_gk61_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_gauss_kronrod",
          61,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t*) {
            return gk_integral<61>(problem, tol, error, L1);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_tanh_sinh_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_tanh_sinh",
          std::nullopt,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t* levels) {
            return tanh_sinh_integral(problem, tol, error, L1, levels);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_exp_sinh_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_exp_sinh",
          std::nullopt,
          tol,
          false,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t* levels) {
            return exp_sinh_integral(problem, tol, error, L1, levels);
          }});
    }
  }

  for (double tol : {1e-4}) {
    {
      std::ostringstream name;
      name << "boost_gk21_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_gauss_kronrod",
          21,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t*) {
            return gk_integral<21>(problem, tol, error, L1);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_gk61_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_gauss_kronrod",
          61,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t*) {
            return gk_integral<61>(problem, tol, error, L1);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_tanh_sinh_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_tanh_sinh",
          std::nullopt,
          tol,
          true,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t* levels) {
            return tanh_sinh_integral(problem, tol, error, L1, levels);
          }});
    }
    {
      std::ostringstream name;
      name << "boost_exp_sinh_tol" << std::scientific << tol;
      methods.push_back(Method{
          name.str(),
          "boost_exp_sinh",
          std::nullopt,
          tol,
          false,
          true,
          [tol](const Problem& problem, double* error, double* L1, std::size_t* levels) {
            return exp_sinh_integral(problem, tol, error, L1, levels);
          }});
    }
  }

  return methods;
}

std::size_t benchmark_iterations(const Method& method, const Problem& problem) {
  volatile double sink = 0.0;
  const auto pilot_start = std::chrono::steady_clock::now();
  for (int i = 0; i < 3; ++i) {
    sink += method.eval(problem, nullptr, nullptr, nullptr);
  }
  const auto pilot_end = std::chrono::steady_clock::now();
  const double pilot_sec =
      std::chrono::duration<double>(pilot_end - pilot_start).count() / 3.0;
  if (!(pilot_sec > 0.0) || !std::isfinite(pilot_sec)) {
    return 1;
  }
  const std::size_t reps = static_cast<std::size_t>(
      std::clamp(0.2 / pilot_sec, 1.0, 100000.0));
  (void)sink;
  return reps;
}

std::string fmt_opt_double(const std::optional<double>& x) {
  if (!x.has_value()) {
    return "";
  }
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(2) << *x;
  return oss.str();
}

std::string fmt_opt_int(const std::optional<int>& x) {
  if (!x.has_value()) {
    return "";
  }
  return std::to_string(*x);
}

}  // namespace

int main(int argc, char** argv) {
  const std::string out_csv =
      argc > 1 ? argv[1]
               : "dev/scripts/scratch_outputs/quadrature_method_benchmark.csv";

  auto problems = build_problems();
  auto methods = build_methods();

  std::ofstream out(out_csv);
  if (!out) {
    std::cerr << "failed to open output file: " << out_csv << "\n";
    return 1;
  }

  out << "problem,domain,method,family,order,tolerance,value,reference,ref_alt,"
         "abs_error,rel_error,abs_ref_disagreement,ms_per_eval,reps,error_estimate,"
         "L1,condition_number,levels\n";

  for (const auto& problem : problems) {
    for (const auto& method : methods) {
      if (problem.domain == "finite" && !method.supports_finite) {
        continue;
      }
      if (problem.domain == "tail" && !method.supports_tail) {
        continue;
      }

      double error_est = std::numeric_limits<double>::quiet_NaN();
      double L1 = std::numeric_limits<double>::quiet_NaN();
      std::size_t levels = 0;
      const double value = method.eval(problem, &error_est, &L1, &levels);
      const long double abs_error =
          std::fabs(static_cast<long double>(value) - problem.reference);
      const long double rel_error =
          abs_error / std::max(std::fabsl(problem.reference), 1e-30L);
      const long double abs_ref_disagreement =
          std::fabs(problem.reference - problem.ref_alt);

      const auto reps = benchmark_iterations(method, problem);
      volatile double sink = 0.0;
      const auto start = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < reps; ++i) {
        sink += method.eval(problem, nullptr, nullptr, nullptr);
      }
      const auto end = std::chrono::steady_clock::now();
      const double ms_per_eval =
          std::chrono::duration<double, std::milli>(end - start).count() /
          static_cast<double>(reps);
      (void)sink;

      const double condition_number =
          (std::isfinite(L1) && std::fabs(value) > 0.0)
              ? L1 / std::fabs(value)
              : std::numeric_limits<double>::quiet_NaN();

      out << problem.name << ','
          << problem.domain << ','
          << method.name << ','
          << method.family << ','
          << fmt_opt_int(method.order) << ','
          << fmt_opt_double(method.tol) << ','
          << std::setprecision(17) << value << ','
          << std::setprecision(21) << static_cast<double>(problem.reference) << ','
          << std::setprecision(21) << static_cast<double>(problem.ref_alt) << ','
          << std::setprecision(12) << static_cast<double>(abs_error) << ','
          << std::setprecision(12) << static_cast<double>(rel_error) << ','
          << std::setprecision(12) << static_cast<double>(abs_ref_disagreement) << ','
          << std::setprecision(9) << ms_per_eval << ','
          << reps << ','
          << std::setprecision(12) << error_est << ','
          << std::setprecision(12) << L1 << ','
          << std::setprecision(12) << condition_number << ','
          << levels << '\n';
    }
  }

  std::cerr << "Wrote quadrature benchmark to " << out_csv << "\n";
  return 0;
}
