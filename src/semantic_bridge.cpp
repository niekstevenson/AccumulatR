// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "eval/likelihood_context.hpp"
#include "eval/observation_likelihood.hpp"

namespace {

Rcpp::List complexity_metrics_list(
    const accumulatr::eval::detail::NativeLikelihoodContext &ctx) {
  const auto n = ctx.exact_complexity_metrics.size();
  if (n == 0U && !ctx.exact_plans.empty()) {
    Rcpp::stop(
        "complexity metrics were not collected; create the context with diagnostics = TRUE");
  }
  Rcpp::IntegerVector variant_index(n);
  Rcpp::IntegerVector symbolic_regions(n);
  Rcpp::IntegerVector symbolic_cells(n);
  Rcpp::IntegerVector max_symbolic_cells_per_region(n);
  Rcpp::IntegerVector negative_symbolic_cells(n);
  Rcpp::IntegerVector overlapping_symbolic_cell_pairs(n);
  Rcpp::IntegerVector expr_relation_atoms(n);
  Rcpp::IntegerVector compiled_roots(n);
  Rcpp::IntegerVector compiled_nodes(n);
  Rcpp::IntegerVector integral_nodes(n);
  Rcpp::IntegerVector integral_kernels(n);
  Rcpp::IntegerVector source_product_integral_kernels(n);
  Rcpp::IntegerVector generic_integral_kernels(n);
  Rcpp::IntegerVector max_integral_depth(n);

  int total_symbolic_regions = 0;
  int total_symbolic_cells = 0;
  int total_negative_symbolic_cells = 0;
  int total_overlapping_symbolic_cell_pairs = 0;
  int total_expr_relation_atoms = 0;
  int total_compiled_roots = 0;
  int total_compiled_nodes = 0;
  int total_integral_nodes = 0;
  int total_integral_kernels = 0;
  int total_source_product_integral_kernels = 0;
  int total_generic_integral_kernels = 0;
  int aggregate_max_symbolic_cells = 0;
  int aggregate_max_integral_depth = 0;

  for (std::size_t i = 0; i < n; ++i) {
    const auto &metrics = ctx.exact_complexity_metrics[i];
    variant_index[i] = static_cast<int>(i);
    symbolic_regions[i] = metrics.symbolic_region_count;
    symbolic_cells[i] = metrics.symbolic_cell_count;
    max_symbolic_cells_per_region[i] =
        metrics.max_symbolic_cells_per_region;
    negative_symbolic_cells[i] = metrics.negative_symbolic_cell_count;
    overlapping_symbolic_cell_pairs[i] =
        metrics.overlapping_symbolic_cell_pair_count;
    expr_relation_atoms[i] = metrics.expr_relation_atom_count;
    compiled_roots[i] = metrics.compiled_root_count;
    compiled_nodes[i] = metrics.compiled_node_count;
    integral_nodes[i] = metrics.integral_node_count;
    integral_kernels[i] = metrics.integral_kernel_count;
    source_product_integral_kernels[i] =
        metrics.source_product_integral_kernel_count;
    generic_integral_kernels[i] = metrics.generic_integral_kernel_count;
    max_integral_depth[i] = metrics.max_integral_depth;

    total_symbolic_regions += metrics.symbolic_region_count;
    total_symbolic_cells += metrics.symbolic_cell_count;
    total_negative_symbolic_cells += metrics.negative_symbolic_cell_count;
    total_overlapping_symbolic_cell_pairs +=
        metrics.overlapping_symbolic_cell_pair_count;
    total_expr_relation_atoms += metrics.expr_relation_atom_count;
    total_compiled_roots += metrics.compiled_root_count;
    total_compiled_nodes += metrics.compiled_node_count;
    total_integral_nodes += metrics.integral_node_count;
    total_integral_kernels += metrics.integral_kernel_count;
    total_source_product_integral_kernels +=
        metrics.source_product_integral_kernel_count;
    total_generic_integral_kernels += metrics.generic_integral_kernel_count;
    aggregate_max_symbolic_cells =
        std::max(
            aggregate_max_symbolic_cells,
            static_cast<int>(metrics.max_symbolic_cells_per_region));
    aggregate_max_integral_depth =
        std::max(
            aggregate_max_integral_depth,
            static_cast<int>(metrics.max_integral_depth));
  }

  return Rcpp::List::create(
      Rcpp::Named("variants") = Rcpp::DataFrame::create(
          Rcpp::Named("variant_index") = variant_index,
          Rcpp::Named("symbolic_regions") = symbolic_regions,
          Rcpp::Named("symbolic_cells") = symbolic_cells,
          Rcpp::Named("max_symbolic_cells_per_region") =
              max_symbolic_cells_per_region,
          Rcpp::Named("negative_symbolic_cells") =
              negative_symbolic_cells,
          Rcpp::Named("overlapping_symbolic_cell_pairs") =
              overlapping_symbolic_cell_pairs,
          Rcpp::Named("expr_relation_atoms") = expr_relation_atoms,
          Rcpp::Named("compiled_roots") = compiled_roots,
          Rcpp::Named("compiled_nodes") = compiled_nodes,
          Rcpp::Named("integral_nodes") = integral_nodes,
          Rcpp::Named("integral_kernels") = integral_kernels,
          Rcpp::Named("source_product_integral_kernels") =
              source_product_integral_kernels,
          Rcpp::Named("generic_integral_kernels") =
              generic_integral_kernels,
          Rcpp::Named("max_integral_depth") = max_integral_depth),
      Rcpp::Named("total") = Rcpp::List::create(
          Rcpp::Named("symbolic_regions") = total_symbolic_regions,
          Rcpp::Named("symbolic_cells") = total_symbolic_cells,
          Rcpp::Named("max_symbolic_cells_per_region") =
              aggregate_max_symbolic_cells,
          Rcpp::Named("negative_symbolic_cells") =
              total_negative_symbolic_cells,
          Rcpp::Named("overlapping_symbolic_cell_pairs") =
              total_overlapping_symbolic_cell_pairs,
          Rcpp::Named("expr_relation_atoms") =
              total_expr_relation_atoms,
          Rcpp::Named("compiled_roots") = total_compiled_roots,
          Rcpp::Named("compiled_nodes") = total_compiled_nodes,
          Rcpp::Named("integral_nodes") = total_integral_nodes,
          Rcpp::Named("integral_kernels") = total_integral_kernels,
          Rcpp::Named("source_product_integral_kernels") =
              total_source_product_integral_kernels,
          Rcpp::Named("generic_integral_kernels") =
              total_generic_integral_kernels,
          Rcpp::Named("max_integral_depth") =
              aggregate_max_integral_depth));
}

void loglik_trials_context(SEXP contextSEXP,
                           SEXP paramsSEXP,
                           SEXP dataSEXP,
                           SEXP okSEXP,
                           const double min_ll,
                           double *out) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  const auto layout =
      accumulatr::eval::detail::read_prepared_trial_layout(dataSEXP);
  const int *ok = Rf_isNull(okSEXP) ? nullptr : LOGICAL(okSEXP);
  accumulatr::eval::detail::evaluate_observation_likelihood_trial_values_cached(
      ctx.observation_plans_by_component_code,
      ctx.observation_is_identity,
      ctx.component_mixture,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      ctx.exact_leaf_row_offsets_by_variant,
      layout,
      paramsSEXP,
      dataSEXP,
      min_ll,
      ok,
      out);
}

Rcpp::NumericVector loglik_context(SEXP contextSEXP,
                                   SEXP paramsSEXP,
                                   SEXP dataSEXP,
                                   SEXP okSEXP,
                                   const double min_ll) {
  const auto layout =
      accumulatr::eval::detail::read_prepared_trial_layout(dataSEXP);
  Rcpp::NumericVector compact(static_cast<R_xlen_t>(layout.trials.size()));
  loglik_trials_context(
      contextSEXP,
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
  return Rcpp::List::create(
      Rcpp::Named("native") = ptr,
      Rcpp::Named("has_complexity_metrics") = diagnostics);
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

extern "C" {

void accumulatr_loglik_trials_ccallable(SEXP contextSEXP,
                                        SEXP paramsSEXP,
                                        SEXP dataSEXP,
                                        SEXP okSEXP,
                                        double min_ll,
                                        double *out) {
  try {
    loglik_trials_context(
        contextSEXP,
        paramsSEXP,
        dataSEXP,
        okSEXP,
        min_ll,
        out);
    return;
  } catch (const std::exception &e) {
    ::Rf_error("%s", e.what());
  } catch (...) {
    ::Rf_error("Unknown C++ exception in AccumulatR::loglik_trials");
  }
}

} // extern "C"

// [[Rcpp::init]]
void accumulatr_register_ccallables(DllInfo *dll) {
  (void)dll;
  R_RegisterCCallable(
      "AccumulatR",
      "loglik_trials",
      reinterpret_cast<DL_FUNC>(accumulatr_loglik_trials_ccallable));
}

// [[Rcpp::export]]
SEXP semantic_response_probabilities_context_cpp(SEXP contextSEXP,
                                                 SEXP paramsSEXP,
                                                 SEXP layoutSEXP) {
  const auto &ctx =
      accumulatr::eval::detail::likelihood_context_from_xptr(contextSEXP);
  return accumulatr::eval::detail::evaluate_response_probabilities_cached(
      ctx.component_mixture,
      ctx.observation_plans_by_component_code,
      ctx.exact_variant_index_by_component_code,
      ctx.exact_plans,
      ctx.exact_leaf_row_offsets_by_variant,
      ctx.outcome_count,
      paramsSEXP,
      layoutSEXP);
}
