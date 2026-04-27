// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R_ext/Rdynload.h>

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
      Rcpp::Named("identity_backend") = "exact",
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
                                 SEXP okSEXP,
                                 SEXP expandSEXP,
                                 SEXP minLLSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  const auto &layout =
      accumulatr::eval::detail::trial_layout_from_xptr(layoutSEXP);
  const Rcpp::LogicalVector ok =
      Rf_isNull(okSEXP) ? Rcpp::LogicalVector() : Rcpp::LogicalVector(okSEXP);
  const Rcpp::IntegerVector expand =
      Rf_isNull(expandSEXP) ? Rcpp::IntegerVector() : Rcpp::IntegerVector(expandSEXP);
  const double min_ll = Rcpp::as<double>(minLLSEXP);
  const auto control = accumulatr::eval::detail::build_trial_eval_control(
      layout.spans.size(),
      ok,
      expand);
  Rcpp::List observed = accumulatr::eval::detail::evaluate_observed_trials_cached(
      ctx.observed_plans_by_component_code,
      ctx.observed_identity,
      ctx.model,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll,
      control.selected.empty() ? nullptr : &control.selected);
  if (!control.weights.empty()) {
    const Rcpp::NumericVector loglik = observed["loglik"];
    observed["total_loglik"] =
        accumulatr::eval::detail::aggregate_trial_loglik(
            loglik,
            control);
  }
  return observed;
}

extern "C" {

double accumulatr_cpp_loglik_ccallable(SEXP contextSEXP,
                                       SEXP paramsSEXP,
                                       SEXP dataSEXP,
                                       SEXP okSEXP,
                                       SEXP expandSEXP,
                                       double min_ll) {
  try {
    SEXP layoutSEXP = Rf_getAttrib(dataSEXP, Rf_install("cpp_layout"));
    if (layoutSEXP == R_NilValue) {
      Rcpp::stop("prepared data are missing native layout metadata");
    }
    Rcpp::List observed = semantic_loglik_context_cpp(
        contextSEXP,
        layoutSEXP,
        paramsSEXP,
        dataSEXP,
        okSEXP,
        expandSEXP,
        Rcpp::wrap(min_ll));
    return Rcpp::as<double>(observed["total_loglik"]);
  } catch (const std::exception &e) {
    ::Rf_error("%s", e.what());
  } catch (...) {
    ::Rf_error("Unknown C++ exception in AccumulatR::cpp_loglik");
  }
  return NA_REAL;
}

} // extern "C"

// [[Rcpp::init]]
void accumulatr_register_ccallables(DllInfo *dll) {
  (void)dll;
  R_RegisterCCallable(
      "AccumulatR",
      "cpp_loglik",
      reinterpret_cast<DL_FUNC>(accumulatr_cpp_loglik_ccallable));
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
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      layout,
      paramsSEXP,
      dataSEXP);
}
