benchmark_models <- list(
  example_1_simple = local({
    structure <- race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("go2", "lognormal") |>
      add_outcome("R1", "go1") |>
      add_outcome("R2", "go2") |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go1.m = log(0.30), go1.s = 0.18,
      go2.m = log(0.32), go2.s = 0.18
    )
    list(structure = structure, pars = pars)
  }),

  example_2_stop_mixture = local({
    structure <- race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("stop", "exgauss", onset = 0.20) |>
      add_accumulator("go2", "lognormal", onset = 0.20) |>
      add_outcome(
        "R1", inhibit("go1", by = "stop"),
        options = list(component = c("go_only", "go_stop"))
      ) |>
      add_outcome(
        "R2", all_of("go2", "stop"),
        options = list(component = "go_stop")
      ) |>
      add_component("go_only", members = "go1") |>
      add_component("go_stop", members = c("go1", "stop", "go2")) |>
      set_mixture(mode = "fixed", weights = c(go_only = 0.5, go_stop = 0.5)) |>
      set_parameters(
        separate = list(m = TRUE, s = TRUE, mu = TRUE, sigma = TRUE, tau = TRUE)
      ) |>
      finalize_model()
    pars <- c(
      go1.m = log(0.35), go1.s = 0.20,
      stop.mu = 0.10, stop.sigma = 0.04, stop.tau = 0.10,
      go2.m = log(0.60), go2.s = 0.18
    )
    list(structure = structure, pars = pars)
  }),

  example_3_stop_na = local({
    structure <- race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_accumulator("stop", "lognormal", onset = 0.15) |>
      add_outcome("Left", "go_left") |>
      add_outcome("Right", "go_right") |>
      add_outcome(
        "STOP", "stop",
        options = list(map_outcome_to = NA_character_)
      ) |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go_left.m = log(0.30), go_left.s = 0.20,
      go_right.m = log(0.32), go_right.s = 0.20,
      stop.m = log(0.15), stop.s = 0.18
    )
    list(structure = structure, pars = pars)
  }),

  example_5_timeout_guess = local({
    structure <- race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_accumulator("timeout", "lognormal", onset = 0.05) |>
      add_outcome("Left", "go_left") |>
      add_outcome("Right", "go_right") |>
      add_outcome(
        "TIMEOUT", "timeout",
        options = list(guess = list(
          labels = c("Left", "Right"),
          weights = c(0.2, 0.8),
          rt_policy = "keep"
        ))
      ) |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go_left.m = log(0.30), go_left.s = 0.18,
      go_right.m = log(0.325), go_right.s = 0.18,
      timeout.m = log(0.25), timeout.s = 0.10
    )
    list(structure = structure, pars = pars)
  }),

  example_6_dual_path = local({
    structure <- race_spec() |>
      add_accumulator("acc_taskA", "lognormal") |>
      add_accumulator("acc_taskB", "lognormal") |>
      add_accumulator("acc_gateC", "lognormal") |>
      add_outcome("Outcome_via_A", all_of("acc_taskA", "acc_gateC")) |>
      add_outcome("Outcome_via_B", all_of("acc_taskB", "acc_gateC")) |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      acc_taskA.m = log(0.28), acc_taskA.s = 0.18,
      acc_taskB.m = log(0.32), acc_taskB.s = 0.18,
      acc_gateC.m = log(0.30), acc_gateC.s = 0.18
    )
    list(structure = structure, pars = pars)
  }),

  example_7_mixture = local({
    structure <- race_spec() |>
      add_accumulator("target_fast", "lognormal") |>
      add_accumulator("target_slow", "lognormal") |>
      add_accumulator("competitor", "lognormal") |>
      add_pool("TARGET", c("target_fast", "target_slow")) |>
      add_outcome("R1", "TARGET") |>
      add_outcome("R2", "competitor") |>
      add_component("fast", members = c("target_fast", "competitor")) |>
      add_component("slow", members = c("target_slow", "competitor")) |>
      set_mixture(mode = "sample", reference = "slow") |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      target_fast.m = log(0.25), target_fast.s = 0.15,
      target_slow.m = log(0.45), target_slow.s = 0.20,
      competitor.m = log(0.35), competitor.s = 0.18,
      p.fast = 0.20
    )
    list(structure = structure, pars = pars)
  }),

  example_10_exclusion = local({
    structure <- race_spec() |>
      add_accumulator("R1_acc", "lognormal") |>
      add_accumulator("R2_acc", "lognormal") |>
      add_accumulator("X_acc", "lognormal") |>
      add_outcome("R1", inhibit("R1_acc", by = "X_acc")) |>
      add_outcome("R2", "R2_acc") |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      R1_acc.m = log(0.35), R1_acc.s = 0.18,
      R2_acc.m = log(0.45), R2_acc.s = 0.18,
      X_acc.m = log(0.35), X_acc.s = 0.18
    )
    list(structure = structure, pars = pars)
  }),

  example_16_guard_tie_simple = local({
    structure <- race_spec() |>
      add_accumulator("go_fast", "lognormal") |>
      add_accumulator("go_slow", "lognormal") |>
      add_accumulator("gate_shared", "lognormal") |>
      add_accumulator("stop_control", "lognormal") |>
      add_outcome(
        "Fast",
        inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")
      ) |>
      add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go_fast.m = log(0.28), go_fast.s = 0.18,
      go_slow.m = log(0.34), go_slow.s = 0.18,
      gate_shared.m = log(0.30), gate_shared.s = 0.16,
      stop_control.m = log(0.27), stop_control.s = 0.15
    )
    list(structure = structure, pars = pars)
  }),

  example_21_simple_q = local({
    structure <- race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("go2", "lognormal") |>
      add_outcome("R1", "go1") |>
      add_outcome("R2", "go2") |>
      add_trigger("shared_trigger", members = c("go1", "go2")) |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go1.m = log(0.30), go1.s = 0.18,
      go2.m = log(0.32), go2.s = 0.18,
      shared_trigger = 0.10
    )
    list(structure = structure, pars = pars)
  }),

  example_22_shared_q = local({
    structure <- race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_outcome("Left", "go_left") |>
      add_outcome("Right", "go_right") |>
      add_trigger("q_left", members = "go_left") |>
      add_trigger("q_right", members = "go_right") |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      go_left.m = log(0.30), go_left.s = 0.18,
      go_right.m = log(0.32), go_right.s = 0.18,
      q_left = 0.10, q_right = 0.10
    )
    list(structure = structure, pars = pars)
  }),

  example_23_ranked_chain = local({
    structure <- race_spec(n_outcomes = 2L) |>
      add_accumulator("a", "lognormal") |>
      add_accumulator("b", "lognormal", onset = after("a")) |>
      add_outcome("A", "a") |>
      add_outcome("B", "b") |>
      set_parameters(separate = list(m = TRUE, s = TRUE)) |>
      finalize_model()
    pars <- c(
      a.m = log(0.30), a.s = 0.16,
      b.m = log(0.22), b.s = 0.16
    )
    list(structure = structure, pars = pars)
  }),

  stop_change_shared_trigger = local({
    structure <- race_spec() |>
      add_accumulator("S", "lognormal") |>
      add_accumulator("stop", "lognormal") |>
      add_accumulator("change", "lognormal") |>
      add_outcome("S", inhibit("S", by = "stop")) |>
      add_outcome("X", all_of("change", "stop")) |>
      add_component("go_only", members = "S") |>
      add_component("go_stop", members = c("S", "stop", "change")) |>
      add_trigger("stop_trigger", members = c("stop", "change")) |>
      set_mixture(mode = "fixed", weights = c(go_only = 0.75, go_stop = 0.25)) |>
      set_parameters(
        separate = list(
          m = c("S", "stop", "change"),
          s = c("S", "stop", "change"),
          t0 = c("S", "change")
        ),
        rename = c(
          S.m = "m_go", stop.m = "m_stop", change.m = "m_change",
          S.s = "s_go", stop.s = "s_stop", change.s = "s_change",
          S.t0 = "t0_go", change.t0 = "t0_change",
          stop_trigger = "q"
        )
      ) |>
      finalize_model()
    pars <- c(
      m_go = log(0.30), s_go = 0.18, t0_go = 0,
      m_stop = log(0.22), s_stop = 0.18,
      m_change = log(0.40), s_change = 0.18, t0_change = 0,
      q = 0.05
    )
    list(structure = structure, pars = pars)
  }),

  stim_selective_stop = local({
    structure <- race_spec() |>
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
      set_parameters(
        separate = list(
          m = c("S1", "IS", "S2"),
          s = c("S1", "IS", "S2"),
          t0 = c("S1", "IS", "S2")
        ),
        share = list(
          m_go = c("A.m", "B.m"),
          s_go = c("A.s", "B.s"),
          t0_go = c("A.t0", "B.t0")
        )
      ) |>
      finalize_model()
    pars <- c(
      m_go = log(0.30), s_go = 0.18, t0_go = 0.05,
      S1.m = log(0.26), S1.s = 0.18, S1.t0 = 0,
      IS.m = log(0.35), IS.s = 0.18, IS.t0 = 0,
      S2.m = log(0.32), S2.s = 0.18, S2.t0 = 0
    )
    list(structure = structure, pars = pars)
  })
)
