// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

#include "../compile/project_semantic.hpp"

// [[Rcpp::export]]
SEXP semantic_compile_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  return accumulatr::compile::detail::to_r_list(model);
}

// [[Rcpp::export]]
Rcpp::CharacterVector semantic_validate_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto issues = accumulatr::semantic::validate_basic(model);
  Rcpp::CharacterVector out(issues.size());
  for (std::size_t i = 0; i < issues.size(); ++i) {
    out[i] = issues[i].message;
  }
  return out;
}

// [[Rcpp::export]]
SEXP semantic_project_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  return accumulatr::compile::to_r_list(compiled);
}
