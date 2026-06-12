# Replication code for the basic workflow in the AccumulatR JSS manuscript.

library(AccumulatR)

basic_spec <- race_spec() |>
  add_accumulator("left", "lognormal") |>
  add_accumulator("right", "lognormal") |>
  add_outcome("left", "left") |>
  add_outcome("right", "right")

print(par_names(basic_spec))

basic_model <- basic_spec |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()

basic_params <- c(
  left.m = log(0.28),
  left.s = 0.16,
  right.m = log(0.35),
  right.s = 0.18
)

basic_param_df <- build_param_matrix(
  basic_model,
  basic_params,
  n_trials = 120
)

basic_sim <- simulate(
  basic_model,
  basic_param_df,
  seed = 20260612
)

print(head(basic_sim, 4))

basic_prepared <- prepare_data(
  basic_model,
  basic_sim[c("trials", "R", "rt")]
)
basic_ctx <- make_context(basic_model)

basic_loglik <- log_likelihood(
  basic_ctx,
  basic_prepared,
  build_param_matrix(basic_model, basic_params, trial_df = basic_prepared)
)
print(basic_loglik)

basic_one_trial <- build_param_matrix(
  basic_model,
  basic_params,
  n_trials = 1
)
print(round(response_probabilities(basic_model, basic_one_trial), 3))
