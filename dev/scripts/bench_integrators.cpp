// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <chrono>
#include "src/semantic_bridge.hpp"

// All particle evaluation and timing happens in C++; R prepares parameters once.
// [[Rcpp::export]]
Rcpp::List benchmark_integrators(SEXP context, SEXP data, Rcpp::List parameters,
    int method, double absolute, double relative, int repetitions = 3) {
  using namespace accumulatr::eval::detail;
  benchmark_method = method;
  kAdaptiveAbsoluteTolerance = absolute;
  kAdaptiveRelativeTolerance = relative;
  const auto &ctx = likelihood_context_from_xptr(context);
  const auto n = XLENGTH(VECTOR_ELT(data, 0)) / ctx.global_leaf_count;
  Rcpp::NumericMatrix loglik(n, parameters.size());
  Rcpp::NumericVector seconds(repetitions);
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(120);
  loglik_trials_context(ctx, parameters[0], data, R_NilValue, std::log(1e-10), &loglik(0, 0));
  benchmark_evaluations = 0;
  for (int repeat = 0; repeat < repetitions; ++repeat) {
    const auto start = std::chrono::steady_clock::now();
    for (int particle = 0; particle < parameters.size(); ++particle) {
      loglik_trials_context(ctx, parameters[particle], data, R_NilValue,
                           std::log(1e-10), &loglik(0, particle));
      if (std::chrono::steady_clock::now() > deadline) Rcpp::stop("Exceeded two minutes");
    }
    seconds[repeat] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  }
  return Rcpp::List::create(Rcpp::Named("seconds") = seconds,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("nodes") = static_cast<double>(benchmark_evaluations) / repetitions);
}
