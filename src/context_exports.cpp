#include <Rcpp.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "context.h"
#include "prep_builder.h"
#include "proto.h"

// [[Rcpp::export]]
SEXP native_context_build(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  return uuber::build_native_context(prep);
}

// [[Rcpp::export]]
Rcpp::RawVector native_prep_serialize_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  uuber::NativePrepProto proto = uuber::build_prep_proto(prep);
  std::vector<std::uint8_t> bytes = uuber::serialize_native_prep(proto);
  Rcpp::RawVector out(bytes.size());
  if (!bytes.empty()) {
    std::copy(bytes.begin(), bytes.end(), out.begin());
  }
  return out;
}

// [[Rcpp::export]]
SEXP native_context_from_proto_cpp(Rcpp::RawVector blob) {
  uuber::NativePrepProto proto = uuber::deserialize_native_prep(
      reinterpret_cast<const std::uint8_t *>(blob.begin()),
      static_cast<std::size_t>(blob.size()));
  std::unique_ptr<uuber::NativeContext> ctx =
      uuber::build_context_from_proto(proto);
  return Rcpp::XPtr<uuber::NativeContext>(ctx.release(), true);
}

// [[Rcpp::export]]
bool native_ctx_invalid(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP) {
    return true;
  }
  return R_ExternalPtrAddr(ptr) == nullptr;
}
