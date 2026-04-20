// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

#include "eval/direct_kernel.hpp"
#include "eval/exact_kernel.hpp"
#include "eval/observed_kernel.hpp"
#include "eval/observed_plan.hpp"
#include "runtime/direct_program.hpp"
#include "runtime/exact_program.hpp"

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

// [[Rcpp::export]]
SEXP semantic_lower_direct_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  const auto lowered = accumulatr::runtime::lower_direct_model(compiled);
  return accumulatr::runtime::to_r_list(lowered);
}

// [[Rcpp::export]]
SEXP semantic_lower_exact_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  const auto lowered = accumulatr::runtime::lower_exact_model(compiled);
  return accumulatr::runtime::to_r_list(lowered);
}

// [[Rcpp::export]]
SEXP semantic_compile_observed_plan_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  return accumulatr::eval::detail::to_r_observed_plan(
      accumulatr::eval::detail::build_component_observation_plans(prep, compiled));
}

// [[Rcpp::export]]
SEXP semantic_direct_loglik_prep_cpp(SEXP prepSEXP,
                                     SEXP paramsSEXP,
                                     SEXP dataSEXP,
                                     SEXP minLLSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  const double min_ll = Rcpp::as<double>(minLLSEXP);
  return accumulatr::eval::detail::evaluate_direct_trials(
      compiled,
      model,
      paramsSEXP,
      dataSEXP,
      min_ll);
}

// [[Rcpp::export]]
SEXP semantic_direct_prob_prep_cpp(SEXP prepSEXP,
                                   SEXP paramsSEXP,
                                   SEXP dataSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  return accumulatr::eval::detail::evaluate_direct_outcome_probabilities(
      compiled,
      model,
      paramsSEXP,
      dataSEXP);
}

// [[Rcpp::export]]
SEXP semantic_exact_loglik_prep_cpp(SEXP prepSEXP,
                                    SEXP paramsSEXP,
                                    SEXP dataSEXP,
                                    SEXP minLLSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  const double min_ll = Rcpp::as<double>(minLLSEXP);
  return accumulatr::eval::detail::evaluate_exact_trials(
      compiled,
      paramsSEXP,
      dataSEXP,
      min_ll);
}

// [[Rcpp::export]]
SEXP semantic_exact_prob_prep_cpp(SEXP prepSEXP,
                                  SEXP paramsSEXP,
                                  SEXP dataSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  return accumulatr::eval::detail::evaluate_exact_outcome_probabilities(
      compiled,
      paramsSEXP,
      dataSEXP);
}

// [[Rcpp::export]]
SEXP semantic_observed_loglik_prep_cpp(SEXP prepSEXP,
                                       SEXP observedPlanSEXP,
                                       SEXP paramsSEXP,
                                       SEXP dataSEXP,
                                       SEXP minLLSEXP) {
  Rcpp::List prep(prepSEXP);
  const auto model = accumulatr::compile::compile_prep(prep);
  const auto compiled = accumulatr::compile::project_semantic_model(model);
  const double min_ll = Rcpp::as<double>(minLLSEXP);
  return accumulatr::eval::detail::evaluate_observed_trials(
      observedPlanSEXP,
      model,
      compiled,
      paramsSEXP,
      dataSEXP,
      min_ll);
}
