// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include <algorithm>
#include <utility>

#include "../inst/include/AccumulatR/api.hpp"
#include "eval/likelihood_context.hpp"
#include "eval/observation_trial_loop.hpp"

namespace accumulatr {

struct PreparedLikelihood {
  PreparedLikelihood(SEXP context_, SEXP data_, int rows, int columns, double floor)
      : context_owner(context_), data_owner(data_),
        context(eval::detail::likelihood_context_from_xptr(context_)),
        parameter_rows(rows), parameter_columns(columns), min_ll(floor) {
    eval::detail::build_observation_likelihood_schedule(
        context.observation_plans_by_component_code,
        context.observation_is_identity, context.component_mixture,
        context.exact_variant_index_by_component_code, context.exact_plans,
        context.exact_leaf_row_offsets_by_variant, context.global_leaf_count,
        data_, &schedule);
  }

  Rcpp::RObject context_owner, data_owner;
  const eval::detail::NativeLikelihoodContext &context;
  eval::detail::ObservationLikelihoodSchedule schedule;
  int parameter_rows, parameter_columns;
  double min_ll;
};

struct LikelihoodWorkspace {
  explicit LikelihoodWorkspace(const PreparedLikelihood &likelihood)
      : lanes(likelihood.context.exact_plans.size()) {}
  eval::detail::ObservationLikelihoodLaneWorkspace lanes;
};

} // namespace accumulatr

namespace {

Rcpp::List complexity_metrics_list(
    const accumulatr::eval::detail::NativeLikelihoodContext &ctx) {
  const auto n = ctx.exact_complexity_metrics.size();
  if (n == 0U && !ctx.exact_plans.empty()) {
    Rcpp::stop(
        "complexity metrics were not collected; create the context with diagnostics = TRUE");
  }
  using Metrics = accumulatr::eval::detail::ExactComplexityMetrics;
  struct Column {
    const char *name;
    accumulatr::semantic::Index Metrics::*member;
    bool maximum{false};
  };
  const Column columns[] = {
      {"symbolic_regions", &Metrics::symbolic_region_count},
      {"symbolic_cells", &Metrics::symbolic_cell_count},
      {"max_symbolic_cells_per_region", &Metrics::max_symbolic_cells_per_region, true},
      {"negative_symbolic_cells", &Metrics::negative_symbolic_cell_count},
      {"overlapping_symbolic_cell_pairs", &Metrics::overlapping_symbolic_cell_pair_count},
      {"expr_relation_atoms", &Metrics::expr_relation_atom_count},
      {"compiled_roots", &Metrics::compiled_root_count},
      {"compiled_nodes", &Metrics::compiled_node_count},
      {"integral_nodes", &Metrics::integral_node_count},
      {"integral_kernels", &Metrics::integral_kernel_count},
      {"source_product_integral_kernels", &Metrics::source_product_integral_kernel_count},
      {"generic_integral_kernels", &Metrics::generic_integral_kernel_count},
      {"max_integral_depth", &Metrics::max_integral_depth, true}};
  Rcpp::IntegerVector indices(n);
  for (std::size_t i = 0; i < n; ++i) indices[i] = i;
  Rcpp::List variants = Rcpp::List::create(Rcpp::Named("variant_index") = indices);
  Rcpp::List totals;
  for (const auto &column : columns) {
    Rcpp::IntegerVector values(n);
    int total = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const int value = ctx.exact_complexity_metrics[i].*column.member;
      values[i] = value;
      total = column.maximum ? std::max(total, value) : total + value;
    }
    variants.push_back(values, column.name);
    totals.push_back(total, column.name);
  }
  return Rcpp::List::create(
      Rcpp::Named("variants") = Rcpp::DataFrame(variants),
      Rcpp::Named("total") = totals);
}

void loglik_trials_context(
    const accumulatr::eval::detail::NativeLikelihoodContext &ctx,
    SEXP paramsSEXP,
    SEXP dataSEXP,
    SEXP okSEXP,
    const double min_ll,
    double *out) {
  accumulatr::eval::detail::LikelihoodWorkspacePool::Lease lease(*ctx.workspace_pool);
  lease.retain_data(dataSEXP);
  auto &workspace = lease.get();
  auto &schedule = lease.slot->schedule;
  if (lease.slot->data != dataSEXP) {
    lease.slot->data = R_NilValue;
    accumulatr::eval::detail::ObservationLikelihoodSchedule replacement;
    accumulatr::eval::detail::build_observation_likelihood_schedule(
        ctx.observation_plans_by_component_code, ctx.observation_is_identity,
        ctx.component_mixture, ctx.exact_variant_index_by_component_code,
        ctx.exact_plans, ctx.exact_leaf_row_offsets_by_variant,
        ctx.global_leaf_count, dataSEXP, &replacement);
    schedule = std::move(replacement);
    lease.slot->data = dataSEXP;
  }
  const int *ok = Rf_isNull(okSEXP) ? nullptr : LOGICAL(okSEXP);
  const accumulatr::eval::detail::ParamMatrixView parameters(
      REAL(paramsSEXP), Rf_nrows(paramsSEXP), schedule.onset);
  const accumulatr::eval::detail::TrustedParamMatrix weights(
      parameters.base, parameters.nrow, Rf_ncols(paramsSEXP),
      ctx.component_mixture.weight_param_count);
  accumulatr::eval::detail::evaluate_observation_likelihood_trial_values_lanes(
      ctx.component_mixture,
      ctx.exact_plans,
      ctx.exact_leaf_row_offsets_by_variant,
      ctx.global_leaf_count,
      schedule,
      parameters,
      weights,
      min_ll,
      ok,
      &workspace,
      out);
}

const accumulatr::LikelihoodAPI *prepared_likelihood_api() {
  static const accumulatr::LikelihoodAPI api{
      [](SEXP context, SEXP data, int rows, int columns, double min_ll) {
        return new accumulatr::PreparedLikelihood(context, data, rows, columns, min_ll);
      },
      [](accumulatr::PreparedLikelihood *likelihood) { delete likelihood; },
      [](const accumulatr::PreparedLikelihood *likelihood) {
        return new accumulatr::LikelihoodWorkspace(*likelihood);
      },
      [](accumulatr::LikelihoodWorkspace *workspace) { delete workspace; },
      [](const accumulatr::PreparedLikelihood *likelihood,
         accumulatr::LikelihoodWorkspace *workspace, const double *values,
         const int *ok, double *out) {
        const auto &ctx = likelihood->context;
        const accumulatr::eval::detail::ParamMatrixView parameters(
            values, likelihood->parameter_rows, likelihood->schedule.onset);
        const accumulatr::eval::detail::TrustedParamMatrix weights(
            values, likelihood->parameter_rows, likelihood->parameter_columns,
            ctx.component_mixture.weight_param_count);
        accumulatr::eval::detail::evaluate_observation_likelihood_trial_values_lanes(
            ctx.component_mixture, ctx.exact_plans,
            ctx.exact_leaf_row_offsets_by_variant, ctx.global_leaf_count,
            likelihood->schedule, parameters, weights, likelihood->min_ll,
            ok, &workspace->lanes, out);
      }};
  return &api;
}

Rcpp::NumericVector loglik_context(SEXP contextSEXP,
                                   SEXP paramsSEXP,
                                   SEXP dataSEXP,
                                   SEXP okSEXP,
                                   const double min_ll) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  Rcpp::NumericVector compact(
      XLENGTH(VECTOR_ELT(dataSEXP, 0)) /
      static_cast<R_xlen_t>(ctx.global_leaf_count));
  loglik_trials_context(
      ctx,
      paramsSEXP,
      dataSEXP,
      okSEXP,
      min_ll,
      REAL(compact));

  const SEXP expandSEXP =
      accumulatr::eval::detail::trusted_data_attr(dataSEXP, "expand");
  if (expandSEXP == R_NilValue || XLENGTH(expandSEXP) == 0) {
    return compact;
  }

  const int *expand = INTEGER(expandSEXP);
  Rcpp::NumericVector out(XLENGTH(expandSEXP));
  for (R_xlen_t i = 0; i < out.size(); ++i) {
    out[i] = compact[expand[i] - 1];
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
SEXP semantic_make_likelihood_context_prep_cpp(SEXP prepSEXP,
                                               SEXP diagnosticsSEXP) {
  Rcpp::List prep(prepSEXP);
  const bool diagnostics = Rcpp::as<bool>(diagnosticsSEXP);
  auto ctx = accumulatr::eval::detail::build_native_likelihood_context(
      prep,
      diagnostics);
  auto ptr = Rcpp::XPtr<accumulatr::eval::detail::NativeLikelihoodContext>(
      new accumulatr::eval::detail::NativeLikelihoodContext(std::move(ctx)),
      true);
  ptr->workspace_pool->owner = ptr;
  return ptr;
}

// [[Rcpp::export]]
SEXP semantic_complexity_metrics_context_cpp(SEXP contextSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  return complexity_metrics_list(ctx);
}

// [[Rcpp::export]]
SEXP semantic_loglik_context_cpp(SEXP contextSEXP,
                                 SEXP paramsSEXP,
                                 SEXP dataSEXP,
                                 SEXP okSEXP,
                                 SEXP minLLSEXP) {
  return loglik_context(
      contextSEXP,
      paramsSEXP,
      dataSEXP,
      okSEXP,
      REAL(minLLSEXP)[0]);
}

// [[Rcpp::init]]
void accumulatr_register_ccallables(DllInfo *dll) {
  (void)dll;
  R_RegisterCCallable(
      "AccumulatR", "likelihood_api",
      reinterpret_cast<DL_FUNC>(prepared_likelihood_api));
}

// [[Rcpp::export]]
SEXP semantic_response_probabilities_context_cpp(SEXP contextSEXP,
                                                 SEXP paramsSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  accumulatr::eval::detail::LikelihoodWorkspacePool::Lease lease(*ctx.workspace_pool);
  auto &workspace = lease.get();
  return accumulatr::eval::detail::evaluate_response_probabilities_cached(
      ctx.component_mixture,
      ctx.observation_plans_by_component_code,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      ctx.exact_leaf_row_offsets_by_variant,
      ctx.outcome_count,
      ctx.global_leaf_count,
      paramsSEXP,
      &workspace.observation);
}
