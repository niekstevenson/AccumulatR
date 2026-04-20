source("dev/scripts/shared_gate_pair_r_helpers.R")

shared_gate_model <- race_spec() |>
  add_accumulator("x1", "lognormal") |>
  add_accumulator("x2", "lognormal") |>
  add_accumulator("gate", "lognormal") |>
  add_outcome("RESP", all_of("x2", "gate")) |>
  add_outcome("NR_RAW", all_of("x1", "gate"),
              options = list(map_outcome_to = NA_character_)) |>
  finalize_model()

shared_gate_params <- c(
  x1.m = log(0.32), x1.s = 0.18, x1.q = 0.00, x1.t0 = 0.00,
  x2.m = log(0.36), x2.s = 0.18, x2.q = 0.00, x2.t0 = 0.00,
  gate.m = log(0.24), gate.s = 0.14, gate.q = 0.00, gate.t0 = 0.00
)

x1_par <- acc_parts("x1", shared_gate_params)
x2_par <- acc_parts("x2", shared_gate_params)
gate_par <- acc_parts("gate", shared_gate_params)

eval_x1 <- function(t) {
  list(
    density = acc_pdf_scalar(t, x1_par),
    survival = acc_survival_scalar(t, x1_par)
  )
}

eval_x2 <- function(t) {
  list(
    density = acc_pdf_scalar(t, x2_par),
    survival = acc_survival_scalar(t, x2_par)
  )
}

eval_gate <- function(t) {
  list(
    density = acc_pdf_scalar(t, gate_par),
    survival = acc_survival_scalar(t, gate_par)
  )
}

compare_density <- function(rt) {
  data_df <- data.frame(
    trial = 1L,
    R = "RESP",
    rt = rt,
    stringsAsFactors = FALSE
  )
  engine <- engine_trial_density_or_mass(shared_gate_model, shared_gate_params, data_df)
  manual <- shared_gate_pair_density_r(eval_x2, eval_x1, eval_gate, rt)
  data.frame(
    response = "RESP",
    rt = rt,
    engine = engine,
    manual = manual,
    diff = engine - manual,
    stringsAsFactors = FALSE
  )
}

density_comparison <- rbind(
  compare_density(0.30),
  compare_density(0.42),
  compare_density(0.47)
)

cat("Shared-gate pair density comparison\n")
print(density_comparison)
