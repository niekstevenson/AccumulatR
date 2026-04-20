# Simplified declarative API for building race-model specifications
# -----------------------------------------------------------------
# Helper functions now live in R/model_definition.R so this file focuses on
# concrete model examples built with the exported DSL.

# --- Example library ----------------------------------------------------------

# Example 1 – simple two-response race built with the helper pipeline
example_1_simple <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("go2", "lognormal") |>
  add_outcome("R1", "go1") |>
  add_outcome("R2", "go2") |>
  finalize_model()
params_example_1_simple <- c(
  go1.m = log(0.30),
  go1.s = 0.18,
  go2.m = log(0.32),
  go2.s = 0.18
)

# Example 2 – stop/go mixture with gating and inhibition
example_2_stop_mixture <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("stop", "exgauss", onset = 0.20) |>
  add_accumulator("go2", "lognormal", onset = 0.20) |>
  add_outcome("R1", inhibit("go1", by = "stop"),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_outcome("R2", all_of("go2", "stop"), options = list(component = "go_stop")) |>
  add_component("go_only", members = c("go1"), attrs = list(component = "go_only"), weight = 0.5) |>
  add_component("go_stop", members = c("go1", "stop", "go2"), attrs = list(component = "go_stop"), weight = 0.5) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_2_stop_mixture <- c(
  go1.m = log(0.35),
  go1.s = 0.2,
  stop.mu = 0.1,
  stop.sigma = 0.04,
  stop.tau = 0.1,
  go2.m = log(0.60),
  go2.s = 0.18
)

# Example 3 – stop outcome mapped to NA
example_3_stop_na <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_accumulator("stop", "lognormal", onset = 0.15) |>
  add_outcome("Left", "go_left") |>
  add_outcome("Right", "go_right") |>
  add_outcome("STOP", "stop", options = list(map_outcome_to = NA_character_)) |>
  finalize_model()
params_example_3_stop_na <- c(
  go_left.m = log(0.30),
  go_left.s = 0.20,
  go_right.m = log(0.32),
  go_right.s = 0.20,
  stop.m = log(0.15),
  stop.s = 0.18
)

# Example 4 – two accumulators vs one
example_4_two_on_one <- race_spec() |>
  add_accumulator("R1_A", "lognormal") |>
  add_accumulator("R1_B", "lognormal") |>
  add_accumulator("R2", "lognormal") |>
  add_pool("R1", c("R1_A", "R1_B")) |>
  add_outcome("R1", "R1") |>
  add_outcome("R2", "R2") |>
  finalize_model()
params_example_4_two_on_one <- c(
  R1_A.m = log(0.30),
  R1_A.s = 0.15,
  R1_B.m = log(0.30),
  R1_B.s = 0.20,
  R2.m = log(0.30),
  R2.s = 0.18
)

# Example 5 – timeout with weighted guess
example_5_timeout_guess <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_accumulator("timeout", "lognormal", onset = 0.05) |>
  add_outcome("Left", "go_left") |>
  add_outcome("Right", "go_right") |>
  add_outcome("TIMEOUT", "timeout", options = list(
    guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
  )) |>
  finalize_model()
params_example_5_timeout_guess <- c(
  go_left.m = log(0.30),
  go_left.s = 0.18,
  go_right.m = log(0.325),
  go_right.s = 0.18,
  timeout.m = log(0.25),
  timeout.s = 0.10
)

# Example 6 – dual path (A & C) OR (B & C)
example_6_dual_path <- race_spec() |>
  add_accumulator("acc_taskA", "lognormal") |>
  add_accumulator("acc_taskB", "lognormal") |>
  add_accumulator("acc_gateC", "lognormal") |>
  add_outcome("Outcome_via_A", all_of("acc_taskA", "acc_gateC")) |>
  add_outcome("Outcome_via_B", all_of("acc_taskB", "acc_gateC")) |>
  finalize_model()
params_example_6_dual_path <- c(
  acc_taskA.m = log(0.28),
  acc_taskA.s = 0.18,
  acc_taskB.m = log(0.32),
  acc_taskB.s = 0.18,
  acc_gateC.m = log(0.30),
  acc_gateC.s = 0.18
)

# Example 7 – fast/slow mixture
example_7_mixture <- race_spec() |>
  add_accumulator("target_fast", "lognormal") |>
  add_accumulator("target_slow", "lognormal") |>
  add_accumulator("competitor", "lognormal") |>
  add_pool("TARGET", c("target_fast", "target_slow")) |>
  add_outcome("R1", "TARGET") |>
  add_outcome("R2", "competitor") |>
  add_component("fast", members = c("target_fast", "competitor"), weight_param = "p_fast") |>
  add_component("slow", members = c("target_slow", "competitor")) |>
  set_mixture_options(mode = "sample", reference = "slow") |>
  finalize_model()
params_example_7_mixture <- c(
  target_fast.m = log(0.25),
  target_fast.s = 0.15,
  target_slow.m = log(0.45),
  target_slow.s = 0.20,
  competitor.m = log(0.35),
  competitor.s = 0.18,
  p_fast = 0.2
)

# Example 8 – shared parameter group
example_8_shared_params <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_outcome("Left", "go_left") |>
  add_outcome("Right", "go_right") |>
  set_parameters(list(shared_m = c("go_left.m", "go_right.m"))) |>
  finalize_model()
params_example_8_shared_params <- c(
  go_left.s = 0.18,
  go_right.s = 0.24,
  shared_m = log(0.32)
)

# Example 9 – advanced k-of-n logic
example_9_advanced_k <- race_spec() |>
  add_accumulator("a1", "lognormal") |>
  add_accumulator("a2", "lognormal") |>
  add_accumulator("a3", "lognormal") |>
  add_accumulator("b1_1", "lognormal") |>
  add_accumulator("b1_2", "lognormal") |>
  add_accumulator("b2_1", "lognormal") |>
  add_accumulator("b2_2", "lognormal") |>
  add_pool("A", c("a1", "a2", "a3"), k = 2L) |>
  add_pool("B1", c("b1_1", "b1_2")) |>
  add_pool("B2", c("b2_1", "b2_2")) |>
  add_outcome("A", "A") |>
  add_outcome("B", all_of("B1", "B2")) |>
  set_parameters(list(
    A_shared_m = c("a2.m", "a3.m"),
    A_shared_s = c("a2.s", "a3.s"),
    Cross_shared_m = c("a1.m", "b1_1.m", "b2_2.m"),
    Cross_shared_s = c("a1.s", "b1_1.s", "b2_2.s")
  )) |>
  finalize_model()
params_example_9_advanced_k <- c(
  b1_2.m = log(0.50),
  b1_2.s = 0.25,
  b2_1.m = log(0.48),
  b2_1.s = 0.25,
  A_shared_m = log(0.575),
  A_shared_s = 0.25,
  Cross_shared_m = log(0.50),
  Cross_shared_s = 0.25
)

# Example 10 – exclusion via guard
example_10_exclusion <- race_spec() |>
  add_accumulator("R1_acc", "lognormal") |>
  add_accumulator("R2_acc", "lognormal") |>
  add_accumulator("X_acc", "lognormal") |>
  add_outcome("R1", inhibit("R1_acc", by = "X_acc")) |>
  add_outcome("R2", "R2_acc") |>
  finalize_model()
params_example_10_exclusion <- c(
  R1_acc.m = log(0.35),
  R1_acc.s = 0.18,
  R2_acc.m = log(0.45),
  R2_acc.s = 0.18,
  X_acc.m = log(0.35),
  X_acc.s = 0.18
)

# Example 11 – inhibitor with protector
example_11_inhibitor_with_protector <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("stop", "exgauss", onset = 0.20) |>
  add_accumulator("go2", "lognormal", onset = 0.20) |>
  add_accumulator("safety", "lognormal", onset = 0.125) |>
  add_outcome("R1", inhibit("go1", by = inhibit("stop", by = "safety")),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_outcome("R2", all_of("go2", "stop"), options = list(component = "go_stop")) |>
  add_component("go_only", members = c("go1"), attrs = list(component = "go_only"), weight = 0.5) |>
  add_component("go_stop", members = c("go1", "stop", "go2", "safety"), attrs = list(component = "go_stop"), weight = 0.5) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_11_inhibitor_with_protector <- c(
  go1.m = log(0.35),
  go1.s = 0.2,
  stop.mu = 0.1,
  stop.sigma = 0.04,
  stop.tau = 0.1,
  go2.m = log(0.60),
  go2.s = 0.18,
  safety.m = log(0.22),
  safety.s = 0.13
)

# Example 12 – nested pools
example_12_nested_pools <- race_spec() |>
  add_accumulator("a1", "lognormal") |>
  add_accumulator("a2", "lognormal") |>
  add_accumulator("a3", "lognormal") |>
  add_accumulator("b1", "lognormal") |>
  add_accumulator("b2", "lognormal") |>
  add_pool("A_line", c("a1", "a2")) |>
  add_pool("A_team", c("A_line", "a3"), k = 2L) |>
  add_pool("B_team", c("b1", "b2")) |>
  add_outcome("TeamA", "A_team") |>
  add_outcome("TeamB", "B_team") |>
  finalize_model()
params_example_12_nested_pools <- c(
  a1.m = log(0.27),
  a1.s = 0.15,
  a2.m = log(0.29),
  a2.s = 0.15,
  a3.m = log(0.33),
  a3.s = 0.16,
  b1.m = log(0.31),
  b1.s = 0.17,
  b2.m = log(0.36),
  b2.s = 0.19
)

# Example 13 – weighted pool
example_13_weighted_pool <- race_spec() |>
  add_accumulator("w1", "lognormal") |>
  add_accumulator("w2", "lognormal") |>
  add_accumulator("w3", "lognormal") |>
  # Compete two items in WEIGHTED against an external competitor w3
  # Compete two items in WEIGHTED against an external competitor w3
  add_pool("WEIGHTED", c("w1", "w2"), k = 1L) |>
  add_outcome("WeightedChoice", "WEIGHTED") |>
  add_outcome("Competitor", "w3") |>
  finalize_model()
params_example_13_weighted_pool <- c(
  w1.m = log(0.30),
  w1.s = 0.18,
  w2.m = log(0.30),
  w2.s = 0.18,
  w3.m = log(0.30),
  w3.s = 0.20
)

# Example 14 – mixture with component metadata
example_14_component_metadata <- race_spec() |>
  # Separate shared accumulators per component to enhance identifiability
  add_accumulator("acc_shared_fast", "lognormal") |>
  add_accumulator("acc_shared_slow", "lognormal") |>
  add_accumulator("acc_fast", "lognormal") |>
  add_accumulator("acc_slow", "lognormal") |>
  add_pool("FAST", c("acc_fast", "acc_shared_fast")) |>
  add_pool("SLOW", c("acc_slow", "acc_shared_slow")) |>
  add_outcome("Response", first_of("FAST", "SLOW")) |>
  add_component("fast", members = c("acc_fast", "acc_shared_fast"), weight = 0.6, attrs = list(
    component = "fast"
  )) |>
  add_component("slow", members = c("acc_slow", "acc_shared_slow"), weight = 0.4, attrs = list(component = "slow")) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_14_component_metadata <- c(
  acc_shared_fast.m = log(0.32),
  acc_shared_fast.s = 0.10,
  acc_shared_slow.m = log(0.60),
  acc_shared_slow.s = 0.12,
  acc_fast.m = log(0.24),
  acc_fast.s = 0.10,
  acc_slow.m = log(0.50),
  acc_slow.s = 0.15
)

# Example 15 – guard tie with stop control
example_15_guard_tie_simple <- race_spec() |>
  add_accumulator("go_fast", "lognormal") |>
  add_accumulator("go_slow", "lognormal") |>
  add_accumulator("gate_shared", "lognormal") |>
  add_accumulator("stop_control", "lognormal") |>
  add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
  add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
  finalize_model()
params_example_15_guard_tie_simple <- c(
  go_fast.m = log(0.28),
  go_fast.s = 0.18,
  go_slow.m = log(0.34),
  go_slow.s = 0.18,
  gate_shared.m = log(0.30),
  gate_shared.s = 0.16,
  stop_control.m = log(0.27),
  stop_control.s = 0.15
)

example_16_k_of_n_inhibitors <- race_spec() |>
  add_accumulator("a1", "lognormal") |>
  add_accumulator("a2", "lognormal") |>
  add_accumulator("a3", "lognormal") |>
  # Define A as a k-of-n pool with 2 members
  add_pool("A", c("a1", "a2", "a3"), k = 2L) |>
  add_accumulator("b1", "lognormal") |>
  add_accumulator("b2", "lognormal") |>
  add_accumulator("b3", "lognormal") |>
  # Define B as a k-of-n pool with 2 members
  add_pool("B", c("b1", "b2", "b3"), k = 2L) |>
  add_pool("B_other_k1", c("b2", "b3"), k = 1L) |>
  add_pool("B_no_b1", c("b2", "b3"), k = 2L) |>
  add_accumulator("c1", "lognormal") |>
  add_accumulator("c2", "lognormal") |>
  add_accumulator("c3", "lognormal") |>
  # Define C as a k-of-n pool with 2 members
  add_pool("C", c("c1", "c2", "c3"), k = 2L) |>
  add_pool("C23", c("c2", "c3"), k = 1L) |>
  add_pool("C13", c("c1", "c3"), k = 1L) |>
  add_accumulator("d1", "lognormal") |>
  add_accumulator("d2", "lognormal") |>
  add_accumulator("d3", "lognormal") |>
  # Define D as a k-of-n pool with 2 members
  add_pool("D", c("d1", "d2", "d3"), k = 2L) |>
  add_outcome("A", "A") |>
  add_outcome("B", first_of(
    inhibit(all_of("b1", "B_other_k1"), by = "a1"),
    "B_no_b1"
  )) |>
  add_outcome("C", first_of(
    inhibit(all_of("c1", "C23"), by = "b2"),
    inhibit(all_of("c2", "C13"), by = "b2")
  )) |>
  add_outcome("D", inhibit("D", by = "c3")) |>
  set_parameters(list(
    shared_A_m = c("a1.m", "a2.m", "a3.m"),
    shared_A_s = c("a1.s", "a2.s", "a3.s"),
    shared_B_m = c("b1.m", "b2.m", "b3.m"),
    shared_B_s = c("b1.s", "b2.s", "b3.s"),
    shared_C_m = c("c1.m", "c2.m", "c3.m"),
    shared_C_s = c("c1.s", "c2.s", "c3.s"),
    shared_D_m = c("d1.m", "d2.m", "d3.m"),
    shared_D_s = c("d1.s", "d2.s", "d3.s")
  )) |>
  finalize_model()
params_example_16_k_of_n_inhibitors <- c(
  shared_A_m = log(0.30),
  shared_A_s = 0.15,
  shared_B_m = log(0.30),
  shared_B_s = 0.15,
  shared_C_m = log(0.30),
  shared_C_s = 0.15,
  shared_D_m = log(0.30),
  shared_D_s = 0.15
)

# Example 18 univalent stop change - bimanual
example_18_univalent_stop_change <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_accumulator("stop", "lognormal", onset = 0.15) |>
  add_accumulator("change2left", "lognormal", onset = 0.20) |>
  add_accumulator("change2right", "lognormal", onset = 0.20) |>
  add_outcome(
    "Left",
    first_of(
      inhibit("go_left", by = "stop"),
      all_of("stop", "change2left")
    ),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_outcome(
    "Right",
    first_of(
      inhibit("go_right", by = "stop"),
      all_of("stop", "change2right")
    ),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_component("go_only", members = c("go_left", "go_right"), attrs = list(component = "go_only"), weight = 0.75) |>
  add_component("go_stop", members = c("go_left", "go_right", "stop", "change2left", "change2right"), attrs = list(component = "go_stop"), weight = 0.25) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_18_univalent_stop_change <- c(
  go_left.m = log(0.35),
  go_left.s = 0.2,
  go_right.m = log(0.35),
  go_right.s = 0.2,
  stop.m = log(.2),
  stop.s = .1,
  change2left.m = log(0.3),
  change2left.s = 0.18,
  change2right.m = log(0.3),
  change2right.s = 0.18
)

# Helper: shared accumulator definitions for stimulus-selective race
.stim_base_spec <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("ignore_high", "lognormal", onset = 0.20) |>
    add_accumulator("ignore_low", "lognormal", onset = 0.20) |>
    add_accumulator("go2", "lognormal", onset = 0.20) |>
    add_accumulator("stop1", "exgauss", onset = 0.20) |>
    add_accumulator("stop2", "exgauss", onset = 0.20)
}

# Stimulus selective stopping
example_19_stim_select_stop <- .stim_base_spec() |>
  add_pool("IGNORE", c("ignore_high", "ignore_low")) |>
  add_outcome(
    "GO1",
    inhibit("go1", by = "stop1"),
    options = list(component = c("go_only", "stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "GO2",
    all_of("go2", inhibit("IGNORE", by = "stop2")),
    options = list(component = c("stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "StopSuccess",
    all_of(
      inhibit("stop1", by = "go1"),
      inhibit("stop2", by = "IGNORE")
    ),
    options = list(component = c("stop_ignore", "stop_true"), map_outcome_to = NA_character_)
  ) |>
  add_component("go_only", members = c("go1"), attrs = list(component = "go_only"), weight = 0.2) |>
  add_component("stop_ignore", members = c("go1", "ignore_high", "go2", "stop1", "stop2"), attrs = list(component = "stop_ignore"), weight = 0.4) |>
  add_component("stop_true", members = c("go1", "ignore_low", "go2", "stop1", "stop2"), attrs = list(component = "stop_true"), weight = 0.4) |>
  set_parameters(list(
    shared_ignore_s = c("ignore_high.s", "ignore_low.s"),
    shared_stop_sigma = c("stop1.sigma", "stop2.sigma"),
    shared_stop_tau = c("stop1.tau", "stop2.tau")
  )) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()

params_example_19_stim_select_stop <- c(
  go1.m = log(0.50),
  go1.s = 0.18,
  ignore_high.m = log(0.20),
  ignore_low.m = log(0.30),
  go2.m = log(0.25),
  go2.s = 0.16,
  stop1.mu = 0.2,
  stop2.mu = 0.15,
  shared_ignore_s = 0.16,
  shared_stop_sigma = 0.05,
  shared_stop_tau = 0.08
)


# Example 20 simple stop change
example_20_simple_stop_change <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_accumulator("change2right", "lognormal", onset = 0.20) |>
  add_outcome("Left", "go_left", options = list(component = "go_only")) |>
  add_outcome(
    "Left",
    inhibit("go_left", by = "stop"),
    options = list(component = "go_stop")
  ) |>
  add_outcome("Right", "go_right", options = list(component = "go_only")) |>
  add_outcome(
    "Right",
    first_of(
      inhibit("go_right", by = "stop"),
      all_of("stop", "change2right")
    ),
    options = list(component = "go_stop")
  ) |>
  add_component("go_only", members = c("go_left", "go_right"), attrs = list(component = "go_only"), weight = 0.5) |>
  add_component("go_stop", members = c("go_left", "go_right", "stop", "change2right"), attrs = list(component = "go_stop"), weight = 0.5) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_20_simple_stop_change <- c(
  go_left.m = log(0.35),
  go_left.s = 0.2,
  go_left.t0 = .1,
  go_right.m = log(0.35),
  go_right.s = 0.2,
  stop.m = log(.35),
  stop.s = 0.2,
  change2right.m = log(0.5),
  change2right.s = 0.18
)


# Example 21 – simple two-response race with trigger failure
example_21_simple_q <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("go2", "lognormal") |>
  add_outcome("R1", "go1") |>
  add_outcome("R2", "go2") |>
  add_trigger(
    "shared_trigger",
    members = c("go1", "go2"),
    q = 0.10,
    draw = "shared"
  ) |>
  finalize_model()

params_example_21_simple_q <- c(
  go1.m = log(0.30),
  go1.s = 0.18,
  go2.m = log(0.32),
  go2.s = 0.18,
  shared_trigger = 0.10 # estimable shared gate parameter (probability)
)

# Example 22 – independent shared q (per-acc Bernoulli) for go accumulators
example_22_shared_q <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_outcome("Left", "go_left") |>
  add_outcome("Right", "go_right") |>
  add_trigger("q_shared",
    members = c("go_left", "go_right"),
    q = 0.10, draw = "independent"
  ) |>
  finalize_model()

params_example_22_shared_q <- c(
  go_left.m = log(0.30),
  go_left.s = 0.18,
  go_right.m = log(0.32),
  go_right.s = 0.18,
  q_shared = 0.10
)

# Example 23 univalent stop change - bimanual
example_23_univalent_stop_change <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("stop", "lognormal", onset = 0.15) |>
  add_accumulator("change2right", "lognormal", onset = 0.20) |>
  add_accumulator("change2right", "lognormal", onset = 0.20) |>
  add_outcome(
    "Left",
    inhibit("go_left", by = "stop"),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_outcome(
    "Right",
    all_of("stop", "change2right"),
    options = list(component = c("go_only", "go_stop"))
  ) |>
  add_component("go_only", members = c("go_left"), attrs = list(component = "go_only"), weight = 0.5) |>
  add_component("go_stop", members = c("go_left", "stop", "change2right"), attrs = list(component = "go_stop"), weight = 0.5) |>
  add_trigger("go_shared_q",
    members = c("go_left"),
    q = 0.10, draw = "independent"
  ) |>
  add_trigger("stop_shared_q",
    members = c("stop", "change2right"),
    q = 0.10, draw = "independent"
  ) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()
params_example_23_univalent_stop_change <- c(
  go_left.m = log(0.35),
  go_left.s = 0.2,
  stop.m = log(.2),
  stop.s = .1,
  change2right.m = log(0.3),
  change2right.s = 0.18,
  go_shared_q = 0.10,
  stop_shared_q = 0.10
)


new_api_examples <- list(
  example_1_simple = example_1_simple,
  example_2_stop_mixture = example_2_stop_mixture,
  example_3_stop_na = example_3_stop_na,
  example_4_two_on_one = example_4_two_on_one,
  example_5_timeout_guess = example_5_timeout_guess,
  example_6_dual_path = example_6_dual_path,
  example_7_mixture = example_7_mixture,
  example_8_shared_params = example_8_shared_params,
  example_9_advanced_k = example_9_advanced_k,
  example_10_exclusion = example_10_exclusion,
  example_11_inhibitor_with_protector = example_11_inhibitor_with_protector,
  example_12_nested_pools = example_12_nested_pools,
  example_13_weighted_pool = example_13_weighted_pool,
  example_14_component_metadata = example_14_component_metadata,
  example_15_guard_tie_simple = example_15_guard_tie_simple,
  example_16_k_of_n_inhibitors = example_16_k_of_n_inhibitors,
  example_18_univalent_stop_change = example_18_univalent_stop_change,
  example_19_stim_select_stop = example_19_stim_select_stop,
  example_20_simple_stop_change = example_20_simple_stop_change,
  example_21_simple_q = example_21_simple_q,
  example_22_shared_q = example_22_shared_q,
  example_23_univalent_stop_change = example_23_univalent_stop_change
)



new_api_example_params <- list(
  example_1_simple = params_example_1_simple,
  example_2_stop_mixture = params_example_2_stop_mixture,
  example_3_stop_na = params_example_3_stop_na,
  example_4_two_on_one = params_example_4_two_on_one,
  example_5_timeout_guess = params_example_5_timeout_guess,
  example_6_dual_path = params_example_6_dual_path,
  example_7_mixture = params_example_7_mixture,
  example_8_shared_params = params_example_8_shared_params,
  example_9_advanced_k = params_example_9_advanced_k,
  example_10_exclusion = params_example_10_exclusion,
  example_11_inhibitor_with_protector = params_example_11_inhibitor_with_protector,
  example_12_nested_pools = params_example_12_nested_pools,
  example_13_weighted_pool = params_example_13_weighted_pool,
  example_14_component_metadata = params_example_14_component_metadata,
  example_15_guard_tie_simple = params_example_15_guard_tie_simple,
  example_16_k_of_n_inhibitors = params_example_16_k_of_n_inhibitors,
  example_18_univalent_stop_change = params_example_18_univalent_stop_change,
  example_19_stim_select_stop = params_example_19_stim_select_stop,
  example_20_simple_stop_change = params_example_20_simple_stop_change,
  example_21_simple_q = params_example_21_simple_q,
  example_22_shared_q = params_example_22_shared_q,
  example_23_univalent_stop_change = params_example_23_univalent_stop_change
)
