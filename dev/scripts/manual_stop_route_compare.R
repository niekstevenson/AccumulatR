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
  m_go = log(0.30), s_go = 0.18, t0_go = 0.00,
  S1.m = log(0.26), S1.s = 0.18, S1.q = 0.00, S1.t0 = 0.00,
  IS.m = log(0.35), IS.s = 0.18, IS.q = 0.00, IS.t0 = 0.00,
  S2.m = log(0.32), S2.s = 0.18, S2.q = 0.00, S2.t0 = 0.00
)

# Primitive accumulator parameters used by the hand derivation.
go_a <- go_parts(stop_params)
go_b <- go_parts(stop_params)
s1_par <- acc_parts("S1", stop_params)
is_par <- acc_parts("IS", stop_params)
s2_par <- acc_parts("S2", stop_params)

# Leaf event: return density/cdf/survival for one accumulator at time t.
leaf_node <- function(p) {
  force(p)
  function(t) {
    list(
      density = acc_pdf_scalar(t, p),
      cdf = acc_cdf_scalar(t, p),
      survival = acc_survival_scalar(t, p)
    )
  }
}

# Guard node: reference succeeds only while the blocker has not yet finished.
# This mirrors the generic engine route:
#   density(t) = f_ref(t) * S_blocker(t)
#   cdf(t)     = integral_0^t f_ref(u) * S_blocker(u) du
guard_node <- function(reference, blocker) {
  function(t) {
    if (!is.finite(t) || t < 0.0) {
      return(list(density = 0.0, cdf = 0.0, survival = 1.0))
    }
    density <- reference(t)$density * blocker(t)$survival
    cdf <- if (t <= 0.0) {
      0.0
    } else {
      integrate(
        function(u) {
          vapply(
            u,
            function(uu) {
              reference(uu)$density * blocker(uu)$survival
            },
            numeric(1)
          )
        },
        lower = 0.0,
        upper = t,
        rel.tol = 1e-8,
        abs.tol = 1e-10,
        subdivisions = 200L
      )$value
    }
    cdf <- min(max(cdf, 0.0), 1.0)
    list(density = density, cdf = cdf, survival = 1.0 - cdf)
  }
}

# AND node: event occurs once all children have occurred.
# For independent child times:
#   F_and(t) = prod_i F_i(t)
#   f_and(t) = sum_i f_i(t) * prod_{j != i} F_j(t)
and_node <- function(...) {
  children <- list(...)
  function(t) {
    vals <- lapply(children, function(fn) fn(t))
    cdfs <- vapply(vals, `[[`, numeric(1), "cdf")
    densities <- vapply(vals, `[[`, numeric(1), "density")
    cdf <- prod(cdfs)
    density <- sum(vapply(
      seq_along(densities),
      function(i) {
        others <- if (length(densities) == 1L) 1.0 else prod(cdfs[-i])
        densities[i] * others
      },
      numeric(1)
    ))
    cdf <- min(max(cdf, 0.0), 1.0)
    list(density = density, cdf = cdf, survival = 1.0 - cdf)
  }
}

# OR / first_of node: event occurs when the first child occurs.
# For independent child times:
#   S_or(t) = prod_i S_i(t)
#   f_or(t) = sum_i f_i(t) * prod_{j != i} S_j(t)
or_node <- function(...) {
  children <- list(...)
  function(t) {
    vals <- lapply(children, function(fn) fn(t))
    survivals <- vapply(vals, `[[`, numeric(1), "survival")
    densities <- vapply(vals, `[[`, numeric(1), "density")
    survival <- prod(survivals)
    density <- sum(vapply(
      seq_along(densities),
      function(i) {
        others <- if (length(densities) == 1L) 1.0 else prod(survivals[-i])
        densities[i] * others
      },
      numeric(1)
    ))
    survival <- min(max(survival, 0.0), 1.0)
    list(density = density, cdf = 1.0 - survival, survival = survival)
  }
}

event_a <- leaf_node(go_a)
event_b <- leaf_node(go_b)
event_s1 <- leaf_node(s1_par)
event_is <- leaf_node(is_par)
event_s2 <- leaf_node(s2_par)

# Exact node graph for the stop component.
early_a <- guard_node(event_a, event_s1)
early_b <- guard_node(event_b, event_s1)
ignore_gate <- guard_node(event_is, event_s2)
late_a <- and_node(event_a, event_s1, ignore_gate)
late_b <- and_node(event_b, event_s1, ignore_gate)

a_out <- or_node(early_a, late_a)
b_out <- or_node(early_b, late_b)

# Split the observed response density into the two branches of first_of().
observed_density_terms <- function(response, rt) {
  if (identical(response, "A")) {
    early_target <- early_a(rt)
    late_target <- late_a(rt)
    competitor_early <- early_b(rt)
    competitor_late <- late_b(rt)
  } else {
    early_target <- early_b(rt)
    late_target <- late_b(rt)
    competitor_early <- early_a(rt)
    competitor_late <- late_a(rt)
  }

  competitor_survival <- competitor_early$survival * competitor_late$survival
  early_term <- early_target$density * late_target$survival * competitor_survival
  late_term <- late_target$density * early_target$survival * competitor_survival
  total <- early_term + late_term

  list(
    early_term = early_term,
    late_term = late_term,
    total = total
  )
}

# Exact observed density for one visible response:
#   p(A_obs at t) = f_A_out(t) * S_B_out(t)
#   p(B_obs at t) = f_B_out(t) * S_A_out(t)
observed_density_exact <- function(response, rt) {
  if (identical(response, "A")) {
    a_out(rt)$density * b_out(rt)$survival
  } else {
    b_out(rt)$density * a_out(rt)$survival
  }
}

probability_exact <- function(response) {
  integrate(
    function(t) {
      vapply(
        t,
        function(tt) observed_density_exact(response, tt),
        numeric(1)
      )
    },
    lower = 0.0,
    upper = Inf,
    rel.tol = 1e-7,
    abs.tol = 1e-9,
    subdivisions = 200L,
    stop.on.error = FALSE
  )$value
}

# Compare the hand-derived visible-response density against the engine.
compare_observed <- function(response, rt) {
  data_df <- data.frame(
    trial = 1L,
    component = "stop",
    R = response,
    rt = rt,
    stringsAsFactors = FALSE
  )
  terms <- observed_density_terms(response, rt)
  manual <- observed_density_exact(response, rt)
  engine <- engine_trial_density_or_mass(stop_model, stop_params, data_df)
  data.frame(
    component = "stop",
    response = response,
    rt = rt,
    early_term = terms$early_term,
    late_term = terms$late_term,
    manual_total = terms$total,
    exact_total = manual,
    engine_total = engine,
    diff = engine - manual,
    stringsAsFactors = FALSE
  )
}

# Under the current engine semantics, map_outcome_to = NA is an observation map.
# The NA likelihood is therefore the leftover observed mass after A and B.
compare_stop_mass <- function() {
  data_df <- data.frame(
    trial = 1L,
    component = "stop",
    R = NA_character_,
    rt = NA_real_,
    stringsAsFactors = FALSE
  )
  engine <- engine_trial_density_or_mass(stop_model, stop_params, data_df)
  prob_a <- probability_exact("A")
  prob_b <- probability_exact("B")
  manual <- 1.0 - prob_a - prob_b
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
