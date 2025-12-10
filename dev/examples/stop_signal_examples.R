source("dev/examples/new_API.R")

# Stop-signal focused examples derived from the new API library

# Example 1 – stimulus-selective stopping with ignore gating
example_1_stimulus_selective <- race_spec() |>
  add_accumulator("go1", "lognormal", meanlog = log(0.28), sdlog = 0.18) |>
  add_accumulator("ignore_high", "lognormal", onset = 0.10,
                  meanlog = log(0.32), sdlog = 0.16) |>
  add_accumulator("ignore_low", "lognormal", onset = 0.10,
                  meanlog = log(0.24), sdlog = 0.16) |>
  add_accumulator("go2", "lognormal", onset = 0.23,
                  meanlog = log(0.38), sdlog = 0.16) |>
  add_accumulator("stop1", "exgauss", onset = 0.125,
                  params = list(mu = 0.06, sigma = 0.05, tau = 0.08)) |>
  add_accumulator("stop2", "exgauss", onset = 0.20,
                  params = list(mu = 0.12, sigma = 0.05, tau = 0.08)) |>
  add_pool("GO1", "go1") |>
  add_pool("IGNORE", c("ignore_high", "ignore_low")) |>
  add_pool("GO2", "go2") |>
  add_pool("STOP1", "stop1") |>
  add_pool("STOP2", "stop2") |>
  add_outcome(
    "GO1",
    inhibit("GO1", by = "STOP1"),
    options = list(component = c("go_only", "stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "GO2",
    all_of("GO2", inhibit("IGNORE", by = "STOP2")),
    options = list(component = c("stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "StopSuccess",
    all_of(
      inhibit("STOP1", by = "GO1"),
      inhibit("STOP2", by = "IGNORE")
    ),
    options = list(component = c("stop_ignore", "stop_true"), map_outcome_to = NA_character_)
  ) |>
  add_group("component:go_only",
            members = c("go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:stop_ignore",
            members = c("go1", "ignore_high", "go2", "stop1", "stop2"),
            attrs = list(component = "stop_ignore")) |>
  add_group("component:stop_true",
            members = c("go1", "ignore_low", "go2", "stop1", "stop2"),
            attrs = list(component = "stop_true")) |>
  add_group("shared_ignore_sdlog",
            members = c("ignore_high", "ignore_low"),
            attrs = list(shared_params = list(sdlog = 0.16))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2"),
            attrs = list(shared_params = list(sigma = 0.05, tau = 0.08))) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.2),
      component("stop_ignore", weight = 0.4),
      component("stop_true", weight = 0.4)
    )
  )) |>
  finalize_model()

# Example 2 – go/stop mixture with late go2 activation when STOP wins
example_2_stop_mixture <- race_spec() |>
  add_accumulator("go1", "lognormal", meanlog = log(0.35), sdlog = 0.2) |>
  add_accumulator("stop", "exgauss", onset = 0.20,
                  params = list(mu = 0.1, sigma = 0.04, tau = 0.1)) |>
  add_accumulator("go2", "lognormal", meanlog = log(0.60), sdlog = 0.18,
                  onset = 0.20) |>
  add_pool("GO1", "go1") |>
  add_pool("STOP", "stop") |>
  add_pool("GO2", "go2") |>
  add_outcome("R1", inhibit("GO1", by = "STOP"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("R2", all_of("GO2", "STOP"), options = list(component = "go_stop")) |>
  add_group("component:go_only", members = c("go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:go_stop", members = c("go1", "stop", "go2"),
            attrs = list(component = "go_stop")) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.5),
      component("go_stop", weight = 0.5)
    )
  )) |>
  finalize_model()

# Example 3 – stop outcome mapped to NA (censor)
example_3_stop_na <- race_spec() |>
  add_accumulator("go_left", "lognormal",
                  params = list(meanlog = log(0.30), sdlog = 0.20)) |>
  add_accumulator("go_right", "lognormal",
                  params = list(meanlog = log(0.32), sdlog = 0.20)) |>
  add_accumulator("stop", "lognormal", onset = 0.15,
                  params = list(meanlog = log(0.15), sdlog = 0.18)) |>
  add_pool("L", "go_left") |>
  add_pool("R", "go_right") |>
  add_pool("STOP", "stop") |>
  add_outcome("Left", "L") |>
  add_outcome("Right", "R") |>
  add_outcome("STOP", "STOP", options = list(map_outcome_to = NA_character_)) |>
  finalize_model()

# Example 4 – mixture with inhibitor vs censor components
# Half the trials behave like Example 2 (stop permits slow response), the other half like Example 3 (stop censors to NA).
example_4_stop_inhibit_vs_censor <- race_spec() |>
  add_accumulator("go_fast", "lognormal", meanlog = log(0.35), sdlog = 0.2) |>
  add_accumulator("go_slow", "lognormal", meanlog = log(0.50), sdlog = 0.18,
                  onset = 0.20) |>
  add_accumulator("stop_inhibit", "exgauss", onset = 0.175) |>
  add_accumulator("stop_censor", "exgauss", onset = 0.175) |>
  add_pool("GO_FAST", "go_fast") |>
  add_pool("GO_SLOW", "go_slow") |>
  add_pool("STOP_INHIBIT", "stop_inhibit") |>
  add_pool("STOP_CENSOR", "stop_censor") |>
  add_pool("STOP_ANY", c("stop_inhibit", "stop_censor")) |>
  add_outcome(
    "FastResponse",
    inhibit("GO_FAST", by = "STOP_ANY"),
    options = list(component = c("inhibit", "censor"))
  ) |>
  add_outcome(
    "SlowResponse",
    all_of("GO_SLOW", "STOP_INHIBIT"),
    options = list(component = "inhibit")
  ) |>
  add_outcome(
    "Stopped",
    "STOP_CENSOR",
    options = list(component = "censor", map_outcome_to = NA_character_)
  ) |>
  add_group("component:inhibit",
            members = c("go_fast", "stop_inhibit", "go_slow"),
            attrs = list(component = "inhibit")) |>
  add_group("component:censor",
            members = c("go_fast", "stop_censor"),
            attrs = list(component = "censor")) |>
  add_group("shared_stop_params",
            members = c("stop_inhibit", "stop_censor"),
            attrs = list(shared_params = list(mu = 0.1, sigma = 0.075, tau = 0.1))) |>
  set_metadata(mixture = list(
    components = list(
      component("inhibit", weight = 0.5),
      component("censor", weight = 0.5)
    )
  )) |>
  finalize_model()

stop_signal_examples <- list(
  example_1_stimulus_selective = example_1_stimulus_selective,
  example_2_stop_mixture = example_2_stop_mixture,
  example_3_stop_na = example_3_stop_na,
  example_4_stop_inhibit_vs_censor = example_4_stop_inhibit_vs_censor
)