#include <Rcpp.h>
#include <limits>

// [[Rcpp::depends(Rcpp)]]

namespace {

Rcpp::Environment fetch_native_env() {
  Rcpp::Environment global = Rcpp::Environment::global_env();
  if (!global.exists(".dist_native_env")) {
    Rcpp::stop("`.dist_native_env` not found; source('R/dist.R') before calling this demo.");
  }
  Rcpp::Environment native_env = global[".dist_native_env"];
  if (!native_env.exists("loaded") || !Rcpp::as<bool>(native_env["loaded"])) {
    Rcpp::stop("Native kernels are not loaded. Call `.load_dist_native()` first.");
  }
  return native_env;
}

Rcpp::Function fetch_native_fn(const std::string& name) {
  Rcpp::Environment native_env = fetch_native_env();
  if (!native_env.exists(name)) {
    Rcpp::stop("Native function '%s' is unavailable; recompile src/distributions.cpp.", name);
  }
  return native_env[name];
}

} // namespace

//' Batch-evaluate buffer likelihoods for multiple parameter tables.
//'
//' This demo shows how an Rcpp-based sampler could stay entirely in C++ once the
//' prep context, structure, and buffer-style trial plan are prepared in R. The
//' helper loops over each parameter data frame, calls the native buffer bridge,
//' and returns the summed log-likelihood per proposal.
//'
//' @param ctx_ptr External pointer returned by `likelihood_build_native_bundle()$context`.
//' @param structure Generator structure (output of `build_generator_structure()`).
//' @param trial_entries Per-trial evaluation records produced by
//'   `build_buffer_trial_entries()` (these are the exact lists R sends to
//'   `native_loglik_from_buffer_cpp`).
//' @param params_list List of parameter data frames (each matching the buffer plan schema).
//' @param component_weights Optional component-weight overrides (data frame).
//' @param default_deadline Default deadline (usually `prep$default_deadline`).
//' @param rel_tol Relative tolerance for guard/NA integrals.
//' @param abs_tol Absolute tolerance for guard/NA integrals.
//' @param max_depth Maximum recursion depth for integrals.
//' @return Numeric vector of total log-likelihoods, one per parameter table.
// [[Rcpp::export]]
Rcpp::NumericVector demo_buffer_batch_loglik_cpp(
    SEXP ctx_ptr,
    Rcpp::List structure,
    Rcpp::List trial_entries,
    Rcpp::List params_list,
    Rcpp::Nullable<Rcpp::DataFrame> component_weights = R_NilValue,
    double default_deadline = NA_REAL,
    double rel_tol = 1e-8,
    double abs_tol = 1e-10,
    int max_depth = 20) {

  if (Rf_isNull(ctx_ptr)) {
    Rcpp::stop("Context pointer is NULL; build it via likelihood_build_native_bundle().");
  }
  R_xlen_t n_tables = params_list.size();
  if (n_tables == 0) {
    Rcpp::stop("params_list must contain at least one parameter data frame.");
  }

  Rcpp::Function native_fn = fetch_native_fn("native_loglik_from_buffer_cpp");
  double effective_deadline = default_deadline;
  if (Rcpp::NumericVector::is_na(effective_deadline)) {
    effective_deadline = std::numeric_limits<double>::infinity();
  }
  Rcpp::NumericVector totals(n_tables);

  for (R_xlen_t i = 0; i < n_tables; ++i) {
    SEXP entry = params_list[i];
    if (Rf_isNull(entry)) {
      Rcpp::stop("params_list[%d] is NULL; supply a data frame.", static_cast<int>(i + 1));
    }
    Rcpp::DataFrame params_df(entry);
    Rcpp::List result = native_fn(ctx_ptr,
                                  structure,
                                  trial_entries,
                                  params_df,
                                  component_weights,
                                  effective_deadline,
                                  rel_tol,
                                  abs_tol,
                                  max_depth);
    totals[i] = Rcpp::as<double>(result["loglik"]);
  }

  return totals;
}
