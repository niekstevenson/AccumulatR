// [[Rcpp::plugins(cpp17)]]
#include <chrono>
#include "src/semantic_bridge.hpp"

// Build with the same headers as the evaluator, not an installed-package pointer.
// [[Rcpp::export]]
SEXP benchmark_context(Rcpp::List prep) {
  return semantic_make_likelihood_context_prep_cpp(prep, Rcpp::wrap(false));
}

struct BenchmarkCall {
  const accumulatr::eval::detail::NativeLikelihoodContext *context;
  SEXP data, parameters;
  double *output;
};

// Parameters, contexts and output buffers are prepared before timing.
Rcpp::NumericVector benchmark_calls(const std::vector<BenchmarkCall> &calls, int repetitions) {
  Rcpp::NumericVector seconds(repetitions);
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(110);
  const auto evaluate = [&](const BenchmarkCall &call) {
    loglik_trials_context(*call.context, call.parameters, call.data, R_NilValue,
                         std::log(1e-10), call.output);
    if (std::chrono::steady_clock::now() > deadline) Rcpp::stop("110 second limit");
  };
  for (const auto &call : calls) evaluate(call);
  for (int repeat = 0; repeat < repetitions; ++repeat) {
    const auto start = std::chrono::steady_clock::now();
    for (const auto &call : calls) {
      evaluate(call);
    }
    seconds[repeat] = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  }
  return seconds;
}

// [[Rcpp::export]]
Rcpp::List benchmark_particles(SEXP context, SEXP data, Rcpp::List parameters,
                              int repetitions = 5) {
  const auto &ctx = accumulatr::eval::detail::likelihood_context_from_xptr(context);
  const int trials = XLENGTH(VECTOR_ELT(data, 0)) / ctx.global_leaf_count;
  Rcpp::NumericMatrix ll(trials, parameters.size());
  std::vector<BenchmarkCall> calls;
  for (int particle = 0; particle < parameters.size(); ++particle)
    calls.push_back({&ctx, data, parameters[particle], &ll(0, particle)});
  return Rcpp::List::create(Rcpp::_["seconds"]=benchmark_calls(calls,repetitions),
                            Rcpp::_["loglik"]=ll);
}

// Round-robin independent calls to several contexts, rather than one hot context.
// [[Rcpp::export]]
Rcpp::List benchmark_interleaved(Rcpp::List contexts, Rcpp::List data,
                                Rcpp::List parameters, int repetitions = 5) {
  Rcpp::List ll(contexts.size());
  std::vector<std::vector<BenchmarkCall>> groups(contexts.size());
  for (int model = 0; model < contexts.size(); ++model) {
    const auto &ctx = accumulatr::eval::detail::likelihood_context_from_xptr(contexts[model]);
    Rcpp::List particles = parameters[model];
    const int trials = XLENGTH(VECTOR_ELT(data[model], 0)) / ctx.global_leaf_count;
    Rcpp::NumericMatrix values(trials, particles.size());
    ll[model] = values;
    for (int p = 0; p < particles.size(); ++p)
      groups[model].push_back({&ctx, data[model], particles[p], &values(0,p)});
  }
  std::vector<BenchmarkCall> calls;
  std::size_t rounds = 0;
  for (const auto &group : groups) rounds = std::max(rounds, group.size());
  for (std::size_t p = 0; p < rounds; ++p)
    for (const auto &group : groups)
      if (p < group.size()) calls.push_back(group[p]);
  return Rcpp::List::create(Rcpp::_["seconds"]=benchmark_calls(calls,repetitions),
                            Rcpp::_["loglik"]=ll);
}
