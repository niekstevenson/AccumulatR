#include "native_context.h"

#include "native_prep_builder.hpp"

namespace uuber {

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  NativePrepProto proto = build_prep_proto(prep);
  std::unique_ptr<NativeContext> ctx = build_context_from_proto(proto);
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber

