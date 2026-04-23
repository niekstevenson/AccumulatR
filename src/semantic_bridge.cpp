// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

#include "eval/likelihood_context.hpp"
#include "eval/observed_kernel.hpp"

// [[Rcpp::export]]
SEXP semantic_make_likelihood_context_prep_cpp(SEXP prepSEXP) {
  Rcpp::List prep(prepSEXP);
  auto ctx = accumulatr::eval::detail::build_native_likelihood_context(prep);
  auto ptr = Rcpp::XPtr<accumulatr::eval::detail::NativeLikelihoodContext>(
      new accumulatr::eval::detail::NativeLikelihoodContext(std::move(ctx)),
      true);
  return Rcpp::List::create(
      Rcpp::Named("native") = ptr,
      Rcpp::Named("observed_identity") = ptr->observed_identity,
      Rcpp::Named("identity_backend") =
          ptr->identity_backend == accumulatr::compile::BackendKind::Direct
              ? "direct"
              : "exact",
      Rcpp::Named("ranked_supported") = ptr->ranked_supported);
}

// [[Rcpp::export]]
SEXP semantic_prepare_trial_layout_cpp(SEXP dataSEXP) {
  return Rcpp::XPtr<accumulatr::eval::detail::PreparedTrialLayout>(
      new accumulatr::eval::detail::PreparedTrialLayout(
          accumulatr::eval::detail::build_native_trial_layout(dataSEXP)),
      true);
}

// [[Rcpp::export]]
SEXP semantic_loglik_context_cpp(SEXP contextSEXP,
                                 SEXP layoutSEXP,
                                 SEXP paramsSEXP,
                                 SEXP dataSEXP,
                                 SEXP minLLSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  const auto &layout =
      accumulatr::eval::detail::trial_layout_from_xptr(layoutSEXP);
  const double min_ll = Rcpp::as<double>(minLLSEXP);
  return accumulatr::eval::detail::evaluate_observed_trials_cached(
      ctx.observed_plans_by_component_code,
      ctx.observed_identity,
      ctx.identity_backend,
      ctx.model,
      ctx.direct_variant_index_by_component_code,
      ctx.direct_plans,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll);
}

// [[Rcpp::export]]
SEXP semantic_probability_context_cpp(SEXP contextSEXP,
                                      SEXP layoutSEXP,
                                      SEXP paramsSEXP,
                                      SEXP dataSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  const auto &layout =
      accumulatr::eval::detail::trial_layout_from_xptr(layoutSEXP);
  return accumulatr::eval::detail::evaluate_outcome_queries_cached(
      ctx.observed_plans_by_component_code,
      ctx.direct_variant_index_by_component_code,
      ctx.direct_plans,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP);
}
