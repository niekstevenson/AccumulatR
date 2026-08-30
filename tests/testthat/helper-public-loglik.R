run_public_loglik <- function(spec, trial_df, params, min_ll = log(1e-10)) {
  structure <- finalize_model(spec)
  prepared <- prepare_data(structure, trial_df)
  parameter_matrix <- build_param_matrix(structure, params, trial_df = prepared)
  as.numeric(log_likelihood(
    make_context(structure), prepared, parameter_matrix, min_ll = min_ll
  ))
}
