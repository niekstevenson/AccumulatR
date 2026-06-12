# Mixture and likelihood example for the AccumulatR JSS manuscript.

library(AccumulatR)

mix_model <- race_spec() |>
  add_accumulator("target_fast", "lognormal") |>
  add_accumulator("target_slow", "lognormal") |>
  add_accumulator("competitor", "lognormal") |>
  add_pool("target", c("target_fast", "target_slow")) |>
  add_outcome("R1", "target") |>
  add_outcome("R2", "competitor") |>
  add_component("fast", members = c("target_fast", "competitor")) |>
  add_component("slow", members = c("target_slow", "competitor")) |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  set_mixture(mode = "sample", reference = "slow") |>
  finalize_model()

mix_params <- c(
  target_fast.m = log(0.25),
  target_fast.s = 0.15,
  target_slow.m = log(0.45),
  target_slow.s = 0.20,
  competitor.m = log(0.35),
  competitor.s = 0.18,
  p.fast = 0.35
)

mix_one_trial <- build_param_matrix(mix_model, mix_params, n_trials = 1)
print(round(response_probabilities(mix_model, mix_one_trial), 3))

mix_param_df <- build_param_matrix(
  mix_model,
  mix_params,
  n_trials = 160
)
mix_sim <- simulate(
  mix_model,
  mix_param_df,
  seed = 20260612
)
mix_data <- mix_sim[c("trials", "R", "rt")]

mix_prepared <- prepare_data(mix_model, mix_data)
mix_ctx <- make_context(mix_model)

neg_loglik_p <- function(theta) {
  est <- mix_params
  est["p.fast"] <- plogis(theta[["logit_p_fast"]])
  params <- build_param_matrix(
    mix_model,
    est,
    trial_df = mix_prepared
  )
  -as.numeric(log_likelihood(mix_ctx, mix_prepared, params))
}

fit_p <- optim(
  c(logit_p_fast = qlogis(0.55)),
  neg_loglik_p,
  method = "BFGS"
)

print(c(
  true = mix_params[["p.fast"]],
  estimated = plogis(fit_p$par[["logit_p_fast"]])
))
