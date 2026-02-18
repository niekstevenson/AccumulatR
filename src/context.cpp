#include "context.h"

#include "prep_builder.h"

#include <cmath>

namespace {

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

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  NativePrepProto proto = build_prep_proto(prep);
  std::unique_ptr<NativeContext> ctx = build_context_from_proto(proto);
  ctx->na_cache_limit = resolve_na_cache_limit();
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber
