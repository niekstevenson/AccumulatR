# Structured model examples for the AccumulatR JSS manuscript.

library(AccumulatR)

pool_model <- race_spec() |>
  add_accumulator("A1", "lognormal") |>
  add_accumulator("A2", "lognormal") |>
  add_accumulator("A3", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_pool("A_pool", c("A1", "A2", "A3"), k = 2L) |>
  add_outcome("A", "A_pool") |>
  add_outcome("B", "B") |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()

pool_params <- c(
  A1.m = log(0.28), A1.s = 0.16,
  A2.m = log(0.28), A2.s = 0.16,
  A3.m = log(0.28), A3.s = 0.16,
  B.m = log(0.34), B.s = 0.18
)
pool_sim <- simulate(
  pool_model,
  build_param_matrix(pool_model, pool_params, n_trials = 80),
  seed = 20260612
)
print(table(pool_sim$R))

logical_model <- race_spec() |>
  add_accumulator("encoding", "lognormal") |>
  add_accumulator("left_route", "lognormal") |>
  add_accumulator("right_route", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome(
    "choice",
    all_of(
      "encoding",
      first_of("left_route", "right_route"),
      none_of("stop")
    )
  ) |>
  finalize_model()
print(names(logical_model$prep$outcomes))

stop_model <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("respond", inhibit("go", "stop")) |>
  finalize_model()
print(names(stop_model$prep$outcomes))

onset_model <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_accumulator("C", "lognormal", onset = after("B")) |>
  add_outcome("A", "A") |>
  add_outcome("C", "C") |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()
print(par_names(onset_model))

trigger_model <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("go2", "lognormal") |>
  add_outcome("R1", "go1") |>
  add_outcome("R2", "go2") |>
  add_trigger("shared_trigger", members = c("go1", "go2")) |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()

trigger_params <- c(
  go1.m = log(0.30),
  go1.s = 0.18,
  go2.m = log(0.35),
  go2.s = 0.18,
  shared_trigger = 0.15
)
trigger_sim <- simulate(
  trigger_model,
  build_param_matrix(trigger_model, trigger_params, n_trials = 80),
  seed = 20260612
)
print(table(trigger_sim$R, useNA = "ifany"))

rank_model <- race_spec(n_outcomes = 2L) |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_outcome("A", "A") |>
  add_outcome("B", "B") |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()

rank_params <- c(
  A.m = log(0.30),
  A.s = 0.18,
  B.m = log(0.38),
  B.s = 0.22
)
print(simulate(
  rank_model,
  build_param_matrix(rank_model, rank_params, n_trials = 3),
  seed = 7
))
