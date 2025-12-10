source("dev/examples/new_API.R")

# Helper: shared accumulator definitions for stimulus-selective race
.stim_base_spec <- function() {
  race_spec() |>
    add_accumulator("go1", "lognormal") |>
    add_accumulator("ignore_high", "lognormal") |>
    add_accumulator("ignore_low", "lognormal") |>
    add_accumulator("go2", "lognormal") |>
    add_accumulator("stop1", "exgauss") |>
    add_accumulator("stop2", "exgauss")
}

# Version 1: Original guard coupling between GO2 and StopSuccess
stim_sel_guarded <- .stim_base_spec() |>
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
            attrs = list(shared_params = list("sdlog"))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2"),
            attrs = list(shared_params = list("sigma", "tau"))) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.2),
      component("stop_ignore", weight = 0.4),
      component("stop_true", weight = 0.4)
    )
  )) |>
  finalize_model()

# Version 2: GO2 requires IGNORE completion (no STOP2 coupling)
stim_sel_ignore_then_go2 <- .stim_base_spec() |>
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
    all_of("IGNORE", "GO2"),
    options = list(component = c("stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "StopSuccess",
    all_of("STOP2", inhibit("STOP1", by = "GO1")),
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
            attrs = list(shared_params = list("sdlog"))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2"),
            attrs = list(shared_params = list("sigma", "tau"))) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.2),
      component("stop_ignore", weight = 0.4),
      component("stop_true", weight = 0.4)
    )
  )) |>
  finalize_model()

# Version 3: Duplicate STOP2 accumulators so guards depend on disjoint primitives
stim_sel_split_stop2 <- race_spec() |>
  add_accumulator("go1", "lognormal") |>
  add_accumulator("ignore_high", "lognormal") |>
  add_accumulator("ignore_low", "lognormal") |>
  add_accumulator("go2", "lognormal") |>
  add_accumulator("stop1", "exgauss") |>
  add_accumulator("stop2_go2", "exgauss") |>
  add_accumulator("stop2_stop", "exgauss") |>
  add_pool("GO1", "go1") |>
  add_pool("IGNORE", c("ignore_high", "ignore_low")) |>
  add_pool("GO2", "go2") |>
  add_pool("STOP1", "stop1") |>
  add_pool("STOP2_GO2", "stop2_go2") |>
  add_pool("STOP2_STOP", "stop2_stop") |>
  add_outcome(
    "GO1",
    inhibit("GO1", by = "STOP1"),
    options = list(component = c("go_only", "stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "GO2",
    all_of("GO2", inhibit("IGNORE", by = "STOP2_GO2")),
    options = list(component = c("stop_ignore", "stop_true"))
  ) |>
  add_outcome(
    "StopSuccess",
    all_of("STOP2_STOP", inhibit("STOP1", by = "GO1")),
    options = list(component = c("stop_ignore", "stop_true"), map_outcome_to = NA_character_)
  ) |>
  add_group("component:go_only",
            members = c("go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:stop_ignore",
            members = c("go1", "ignore_high", "go2", "stop1",
                        "stop2_go2", "stop2_stop"),
            attrs = list(component = "stop_ignore")) |>
  add_group("component:stop_true",
            members = c("go1", "ignore_low", "go2", "stop1",
                        "stop2_go2", "stop2_stop"),
            attrs = list(component = "stop_true")) |>
  add_group("shared_ignore_sdlog",
            members = c("ignore_high", "ignore_low"),
            attrs = list(shared_params = list("sdlog"))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2_go2", "stop2_stop"),
            attrs = list(shared_params = list("sigma", "tau"))) |>
  add_group("shared_ignore_sdlog",
            members = c("ignore_high", "ignore_low"),
            attrs = list(shared_params = list("sdlog"))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2_go2", "stop2_stop"),
            attrs = list(shared_params = list("sigma", "tau"))) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.2),
      component("stop_ignore", weight = 0.4),
      component("stop_true", weight = 0.4)
    )
  )) |>
  finalize_model()

# Version 4: Guarded model with shared trigger failures across response groups
stim_sel_guarded_triggered <- .stim_base_spec() |>
  add_accumulator("trigger_go1", "lognormal") |>
  add_accumulator("trigger_stop", "lognormal") |>
  add_accumulator("trigger_ignore_go2", "lognormal") |>
  add_pool("GO1", c("trigger_go1", "go1"), k = 2) |>
  add_pool("IGNORE_HIGH_TRIGGERED", c("trigger_ignore_go2", "ignore_high"), k = 2) |>
  add_pool("IGNORE_LOW_TRIGGERED", c("trigger_ignore_go2", "ignore_low"), k = 2) |>
  add_pool("IGNORE", c("IGNORE_HIGH_TRIGGERED", "IGNORE_LOW_TRIGGERED")) |>
  add_pool("GO2", c("trigger_ignore_go2", "go2"), k = 2) |>
  add_pool("STOP1", c("trigger_stop", "stop1"), k = 2) |>
  add_pool("STOP2", c("trigger_stop", "stop2"), k = 2) |>
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
            members = c("go1", "trigger_go1"),
            attrs = list(component = "go_only")) |>
  add_group("component:stop_ignore",
            members = c("go1", "ignore_high", "go2", "stop1", "stop2",
                        "trigger_go1", "trigger_stop", "trigger_ignore_go2"),
            attrs = list(component = "stop_ignore")) |>
  add_group("component:stop_true",
            members = c("go1", "ignore_low", "go2", "stop1", "stop2",
                        "trigger_go1", "trigger_stop", "trigger_ignore_go2"),
            attrs = list(component = "stop_true")) |>
  add_group("shared_ignore_sdlog",
            members = c("ignore_high", "ignore_low"),
            attrs = list(shared_params = list("sdlog"))) |>
  add_group("shared_stop_params",
            members = c("stop1", "stop2"),
            attrs = list(shared_params = list("sigma", "tau"))) |>
  add_group("shared_trigger_timing",
            members = c("trigger_go1", "trigger_stop", "trigger_ignore_go2"),
            attrs = list(shared_params = list("meanlog", "sdlog"))) |>
  set_metadata(mixture = list(
    components = list(
      component("go_only", weight = 0.2),
      component("stop_ignore", weight = 0.4),
      component("stop_true", weight = 0.4)
    )
  )) |>
  finalize_model()

stim_selective_versions <- list(
  guarded = stim_sel_guarded,
  ignore_then_go2 = stim_sel_ignore_then_go2,
  split_stop2 = stim_sel_split_stop2,
  guarded_triggered = stim_sel_guarded_triggered
)

stim_selective_base_params <- c(
  go1.meanlog = log(0.50),
  go1.sdlog = 0.18,
  ignore_high.meanlog = log(0.20),
  ignore_high.onset = 0.20,
  ignore_low.meanlog = log(0.30),
  ignore_low.onset = 0.20,
  go2.meanlog = log(0.25),
  go2.sdlog = 0.16,
  go2.onset = 0.20,
  stop1.mu = 0.2,
  stop1.onset = 0.20,
  stop2.mu = 0.15,
  stop2.onset = 0.20,
  `shared_ignore_sdlog.sdlog` = 0.16,
  `shared_stop_params.sigma` = 0.05,
  `shared_stop_params.tau` = 0.08
)

params_stim_sel_guarded <- stim_selective_base_params
params_stim_sel_ignore_then_go2 <- stim_selective_base_params

params_stim_sel_split_stop2 <- c(
  go1.meanlog = log(0.28),
  go1.sdlog = 0.18,
  ignore_high.meanlog = log(0.30),
  ignore_high.onset = 0.10,
  ignore_low.meanlog = log(0.22),
  ignore_low.onset = 0.10,
  go2.meanlog = log(0.38),
  go2.sdlog = 0.16,
  go2.onset = 0.23,
  stop1.mu = 0.06,
  stop1.onset = 0.125,
  stop2_go2.mu = 0.12,
  stop2_go2.onset = 0.20,
  stop2_stop.mu = 0.12,
  stop2_stop.onset = 0.20,
  `shared_ignore_sdlog.sdlog` = 0.16,
  `shared_stop_params.sigma` = 0.05,
  `shared_stop_params.tau` = 0.08
)

params_stim_sel_guarded_triggered <- c(
  stim_selective_base_params,
  trigger_go1.q = 0.1,
  trigger_stop.q = 0.1,
  trigger_ignore_go2.q = 0.1,
  `shared_trigger_timing.meanlog` = log(0.02),
  `shared_trigger_timing.sdlog` = 0.05
)