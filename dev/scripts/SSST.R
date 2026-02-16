stim_selective <- race_spec() |>
  add_accumulator("A1", "lognormal") |>
  add_accumulator("B1", "lognormal") |>
  add_accumulator("S1", "lognormal") |>
  add_accumulator("IS", "lognormal") |>
  add_accumulator("S2", "lognormal") |>
  add_accumulator("A2", "lognormal") |>
  add_accumulator("B2", "lognormal") |>
  add_outcome("A", first_of(
      inhibit("A1", by = "S1"),"A2")
  ) |>
  add_outcome("B", first_of(
    inhibit("B1", by = "S1"),"B2")
  ) |>
  add_outcome(
    "STOP", inhibit("S2", by = "IS"), 
    options = list(map_outcome_to = NA_character_)
  ) |>
  finalize_model()

# A vs B, both subject to S1; if stop branch wins (S1 and S2 beats IS) return NA.
# All accumulators start at time 0 (no onset = after).
s1_gate_is_vs_s2 <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_accumulator("S1", "lognormal") |>
  add_accumulator("IS", "lognormal") |>
  add_accumulator("S2", "lognormal") |>
  add_outcome(
    "A",
    inhibit("A", by = all_of("S1", inhibit("S2", by = "IS")))
  ) |>
  add_outcome(
    "B",
    inhibit("B", by = all_of("S1", inhibit("S2", by = "IS")))
  ) |>
  add_outcome(
    "STOP",
    all_of("S1", inhibit("S2", by = "IS")),
    options = list(map_outcome_to = NA_character_)
  ) |>
  add_component("go", members = c("A", "B")) |>
  add_component("stop", members = c("A", "B", "S1", "IS", "S2")) |>
  finalize_model()
