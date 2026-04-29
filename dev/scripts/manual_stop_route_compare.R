source("dev/scripts/shared_gate_pair_r_helpers.R")

stop_model <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_accumulator("S1", "lognormal") |>
  add_accumulator("IS", "lognormal") |>
  add_accumulator("S2", "lognormal") |>
  add_outcome(
    "A",
    first_of(
      inhibit("A", by = "S1"),
      all_of("A", "S1", inhibit("IS", by = "S2"))
    )
  ) |>
  add_outcome(
    "B",
    first_of(
      inhibit("B", by = "S1"),
      all_of("B", "S1", inhibit("IS", by = "S2"))
    )
  ) |>
  add_outcome(
    "STOP",
    all_of("S1", inhibit("S2", by = "IS")),
    options = list(map_outcome_to = NA_character_)
  ) |>
  add_component("go", members = c("A", "B")) |>
  add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
  set_parameters(list(
    m_go = c("A.m", "B.m"),
    s_go = c("A.s", "B.s"),
    t0_go = c("A.t0", "B.t0")
  )) |>
  finalize_model()

stop_params <- c(
  m_go = log(0.30), s_go = 0.18, t0_go = 0.05,
  S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
  IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
  S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
)

go_a <- go_parts(stop_params)
go_b <- go_parts(stop_params)
s1_par <- acc_parts("S1", stop_params)
is_par <- acc_parts("IS", stop_params)
s2_par <- acc_parts("S2", stop_params)

integrate_scalar <- function(fn,
                             lower,
                             upper,
                             rel_tol = 1e-10,
                             abs_tol = 1e-12,
                             subdivisions = 400L) {
  if (is.finite(upper) && upper <= lower) {
    return(0.0)
  }
  integrate(
    function(x) {
      vapply(
        x,
        function(xx) {
          value <- fn(xx)
          if (!is.finite(value)) 0.0 else value
        },
        numeric(1)
      )
    },
    lower = lower,
    upper = upper,
    rel.tol = rel_tol,
    abs.tol = abs_tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )$value
}

ignore_density <- function(t) {
  acc_pdf_scalar(t, is_par) * acc_survival_scalar(t, s2_par)
}

ignore_cdf <- function(t) {
  integrate_scalar(ignore_density, 0.0, t)
}

late_target_integrand <- function(t, target, competitor) {
  acc_pdf_scalar(t, target) *
    acc_cdf_scalar(t, s1_par) *
    acc_survival_scalar(t, competitor)
}

early_target_integrand <- function(t, target, competitor) {
  acc_pdf_scalar(t, target) *
    acc_survival_scalar(t, s1_par) *
    acc_survival_scalar(t, competitor)
}

late_target_cdf_part <- function(t, target, competitor) {
  integrate_scalar(
    function(u) late_target_integrand(u, target, competitor),
    0.0,
    t
  )
}

# Correct stimulus-selective stop semantics:
# S1 gates the visible route; it is not an independent response clock.
# If the target accumulator already beat its competitor, later ignore release
# can reveal that target identity.
visible_density_terms <- function(target, competitor, t) {
  target_at_t <- acc_pdf_scalar(t, target) *
    acc_survival_scalar(t, competitor) *
    (
      acc_survival_scalar(t, s1_par) +
        acc_cdf_scalar(t, s1_par) * ignore_cdf(t)
    )

  ignore_release_at_t <- ignore_density(t) *
    late_target_cdf_part(t, target, competitor)

  list(
    target_at_t = target_at_t,
    ignore_release_at_t = ignore_release_at_t,
    total = target_at_t + ignore_release_at_t
  )
}

visible_density <- function(response, t) {
  if (identical(response, "A")) {
    visible_density_terms(go_a, go_b, t)$total
  } else {
    visible_density_terms(go_b, go_a, t)$total
  }
}

visible_probability_fubini <- function(target, competitor) {
  early_mass <- integrate_scalar(
    function(t) early_target_integrand(t, target, competitor),
    0.0,
    Inf
  )
  ignore_mass <- integrate_scalar(ignore_density, 0.0, Inf)
  late_target_mass <- integrate_scalar(
    function(t) late_target_integrand(t, target, competitor),
    0.0,
    Inf
  )
  early_mass + ignore_mass * late_target_mass
}

compare_observed <- function(response, rt) {
  data_df <- data.frame(
    trial = 1L,
    component = "stop",
    R = response,
    rt = rt,
    stringsAsFactors = FALSE
  )
  terms <- if (identical(response, "A")) {
    visible_density_terms(go_a, go_b, rt)
  } else {
    visible_density_terms(go_b, go_a, rt)
  }
  engine <- engine_trial_density_or_mass(stop_model, stop_params, data_df)
  data.frame(
    component = "stop",
    response = response,
    rt = rt,
    target_at_t = terms$target_at_t,
    ignore_release_at_t = terms$ignore_release_at_t,
    manual_total = terms$total,
    engine_total = engine,
    diff = engine - terms$total,
    stringsAsFactors = FALSE
  )
}

compare_stop_mass <- function() {
  data_df <- data.frame(
    trial = 1L,
    component = "stop",
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  engine <- engine_trial_density_or_mass(stop_model, stop_params, data_df)
  prob_a <- visible_probability_fubini(go_a, go_b)
  prob_b <- visible_probability_fubini(go_b, go_a)
  manual <- max(0.0, 1.0 - prob_a - prob_b)
  data.frame(
    component = "stop",
    response = NA_character_,
    rt = NA_real_,
    prob_a = prob_a,
    prob_b = prob_b,
    complement_na = manual,
    engine_total = engine,
    diff = engine - manual,
    stringsAsFactors = FALSE
  )
}

observed_comparison <- rbind(
  compare_observed("A", 0.42),
  compare_observed("B", 0.47)
)

stop_mass_comparison <- compare_stop_mass()

cat("Observed stop-component decomposition\n")
print(observed_comparison)
cat("\nSTOP / NA mass comparison\n")
print(stop_mass_comparison)
