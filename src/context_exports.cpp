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

// [[Rcpp::export]]
Rcpp::DataFrame native_outcome_labels_cpp(SEXP ctxSEXP) {
  Rcpp::XPtr<uuber::NativeContext> ctx(ctxSEXP);
  const R_xlen_t n = static_cast<R_xlen_t>(ctx->outcome_info.size());
  Rcpp::CharacterVector labels(n);
  Rcpp::IntegerVector label_ids(n);
  Rcpp::IntegerVector node_ids(n);
  Rcpp::IntegerVector competitor_counts(n);
  Rcpp::LogicalVector maps_to_na(n);

  for (R_xlen_t idx = 0; idx < n; ++idx) {
    const uuber::OutcomeContextInfo &info =
        ctx->outcome_info[static_cast<std::size_t>(idx)];
    if (idx < static_cast<R_xlen_t>(ctx->outcome_labels.size())) {
      labels[idx] = ctx->outcome_labels[static_cast<std::size_t>(idx)];
    } else {
      labels[idx] = "";
    }
    if (idx < static_cast<R_xlen_t>(ctx->outcome_label_ids.size())) {
      label_ids[idx] = ctx->outcome_label_ids[static_cast<std::size_t>(idx)];
    } else {
      label_ids[idx] = NA_INTEGER;
    }
    node_ids[idx] = info.node_id;
    competitor_counts[idx] = static_cast<int>(info.competitor_ids.size());
    maps_to_na[idx] = info.maps_to_na;
  }

  return Rcpp::DataFrame::create(Rcpp::Named("label") = labels,
                                 Rcpp::Named("label_id") = label_ids,
                                 Rcpp::Named("node_id") = node_ids,
                                 Rcpp::Named("competitors") =
                                     competitor_counts,
                                 Rcpp::Named("maps_to_na") = maps_to_na,
                                 Rcpp::Named("stringsAsFactors") = false);
}
