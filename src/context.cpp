#include "context.h"

#include "prep_builder.h"

#include <atomic>
#include <cmath>

namespace {

std::atomic<std::uint64_t> g_native_context_runtime_cache_instance_counter{1};

int resolve_na_cache_limit() {
  int default_limit = 128;
  Rcpp::Function getOption("getOption");
  Rcpp::RObject opt =
      getOption("uuber.cache.na.max_per_trial", Rcpp::wrap(default_limit));
  if (opt.isNULL()) {
    return default_limit;
  }
  try {
    double val = Rcpp::as<double>(opt);
    if (!std::isfinite(val) || val < 0.0) {
      return 0;
    }
    return static_cast<int>(std::floor(val));
  } catch (...) {
    return default_limit;
  }
}

} // namespace

namespace uuber {

void clear_ranked_distribution_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept;
void clear_exact_outcome_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept;

std::uint64_t next_native_context_runtime_cache_instance_id() noexcept {
  return g_native_context_runtime_cache_instance_counter.fetch_add(
      1, std::memory_order_relaxed);
}

void clear_native_context_runtime_caches(
    std::uint64_t runtime_cache_instance_id) noexcept {
  if (runtime_cache_instance_id == 0) {
    return;
  }
  clear_ranked_distribution_runtime_caches(runtime_cache_instance_id);
  clear_exact_outcome_runtime_caches(runtime_cache_instance_id);
}

NativeContext::~NativeContext() {
  clear_native_context_runtime_caches(runtime_cache_instance_id);
}

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  NativePrepProto proto = build_prep_proto(prep);
  std::unique_ptr<NativeContext> ctx = build_context_from_proto(proto);
  ctx->na_cache_limit = resolve_na_cache_limit();
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber
