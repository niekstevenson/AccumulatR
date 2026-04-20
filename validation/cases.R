validation_cases <- function() {
  list(
    independent_trigger_two_way = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_outcome("A", "a") |>
        add_outcome("B", "b") |>
        finalize_model()
      params <- c(
        a.m = log(0.32), a.s = 0.18, a.q = 0.15, a.t0 = 0.02,
        b.m = log(0.39), b.s = 0.16, b.q = 0.30, b.t0 = 0.01
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      rows <- list()
      for (check in list(list(label = "A", rt = 0.45), list(label = "B", rt = 0.52))) {
        target <- if (identical(check$label, "A")) a else b
        competitor <- if (identical(check$label, "A")) b else a
        manual <- acc_pdf_scalar(check$rt, target) * acc_survival_scalar(check$rt, competitor)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = check$label, rt = check$rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "independent_trigger_two_way",
          paste0(check$label, "_rt_", format(check$rt, nsmall = 2)),
          engine,
          manual,
          1e-9,
          "Two-way race with independent q semantics"
        )
      }
      do.call(rbind, rows)
    },

    pool_vs_competitor = function() {
      structure <- race_spec() |>
        add_accumulator("a1", "lognormal") |>
        add_accumulator("a2", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_pool("A_pool", c("a1", "a2"), k = 1L) |>
        add_outcome("A", "A_pool") |>
        add_outcome("B", "b") |>
        finalize_model()
      params <- c(
        a1.m = log(0.28), a1.s = 0.17, a1.q = 0.00, a1.t0 = 0.00,
        a2.m = log(0.34), a2.s = 0.15, a2.q = 0.00, a2.t0 = 0.00,
        b.m = log(0.42), b.s = 0.18, b.q = 0.00, b.t0 = 0.00
      )
      a1 <- acc_parts("a1", params)
      a2 <- acc_parts("a2", params)
      b <- acc_parts("b", params)
      rows <- list()
      rows[[1L]] <- {
        rt <- 0.40
        manual <- (
          acc_pdf_scalar(rt, a1) * acc_survival_scalar(rt, a2) +
            acc_pdf_scalar(rt, a2) * acc_survival_scalar(rt, a1)
        ) * acc_survival_scalar(rt, b)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "A", rt = rt, stringsAsFactors = FALSE)
        )
        check_row("pool_vs_competitor", "A_rt_0.40", engine, manual, 1e-9, "k = 1 pool winner against competitor")
      }
      rows[[2L]] <- {
        rt <- 0.48
        manual <- acc_pdf_scalar(rt, b) *
          acc_survival_scalar(rt, a1) *
          acc_survival_scalar(rt, a2)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "B", rt = rt, stringsAsFactors = FALSE)
        )
        check_row("pool_vs_competitor", "B_rt_0.48", engine, manual, 1e-9, "Competitor outruns pooled branch")
      }
      do.call(rbind, rows)
    },

    all_of_three_way = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("c", "lognormal") |>
        add_accumulator("d", "lognormal") |>
        add_outcome("ABC", all_of("a", "b", "c")) |>
        add_outcome("D", "d") |>
        finalize_model()
      params <- c(
        a.m = log(0.29), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.34), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
        c.m = log(0.38), c.s = 0.14, c.q = 0.00, c.t0 = 0.00,
        d.m = log(0.44), d.s = 0.17, d.q = 0.00, d.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      c <- acc_parts("c", params)
      d <- acc_parts("d", params)
      rt <- 0.46
      manual <- (
        acc_pdf_scalar(rt, a) * acc_cdf_scalar(rt, b) * acc_cdf_scalar(rt, c) +
          acc_pdf_scalar(rt, b) * acc_cdf_scalar(rt, a) * acc_cdf_scalar(rt, c) +
          acc_pdf_scalar(rt, c) * acc_cdf_scalar(rt, a) * acc_cdf_scalar(rt, b)
      ) * acc_survival_scalar(rt, d)
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "ABC", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "all_of_three_way",
        "ABC_rt_0.46",
        engine,
        manual,
        2e-3,
        "Three-child all_of observed against competitor"
      )
    },

    first_of_three_way = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("c", "lognormal") |>
        add_accumulator("d", "lognormal") |>
        add_outcome("ABC_FIRST", first_of("a", "b", "c")) |>
        add_outcome("D", "d") |>
        finalize_model()
      params <- c(
        a.m = log(0.27), a.s = 0.14, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.31), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
        c.m = log(0.36), c.s = 0.15, c.q = 0.00, c.t0 = 0.00,
        d.m = log(0.43), d.s = 0.17, d.q = 0.00, d.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      c <- acc_parts("c", params)
      d <- acc_parts("d", params)
      rt <- 0.33
      manual <- (
        acc_pdf_scalar(rt, a) * acc_survival_scalar(rt, b) * acc_survival_scalar(rt, c) +
          acc_pdf_scalar(rt, b) * acc_survival_scalar(rt, a) * acc_survival_scalar(rt, c) +
          acc_pdf_scalar(rt, c) * acc_survival_scalar(rt, a) * acc_survival_scalar(rt, b)
      ) * acc_survival_scalar(rt, d)
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "ABC_FIRST", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "first_of_three_way",
        "ABC_FIRST_rt_0.33",
        engine,
        manual,
        2e-3,
        "Three-child first_of observed against competitor"
      )
    },

    chained_onset_single_outcome = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal", onset = after("a")) |>
        add_outcome("B", "b") |>
        finalize_model()
      params <- c(
        a.m = log(0.25), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.35), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      rt <- 0.85
      manual <- integrate_scalar(
        function(u) acc_pdf_scalar(u, a) * acc_pdf_scalar(rt - u, b),
        0.0,
        rt
      )
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "B", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "chained_onset_single_outcome",
        "B_rt_0.85",
        engine,
        manual,
        2e-3,
        "Single-outcome after(a) convolution"
      )
    },

    ranked_independent = function() {
      structure <- race_spec(n_outcomes = 2L) |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_outcome("A", "a") |>
        add_outcome("B", "b") |>
        finalize_model()
      params <- c(
        a.m = log(0.30), a.s = 0.20, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.45), b.s = 0.25, b.q = 0.00, b.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      data_df <- data.frame(
        trial = 1L,
        R = "A",
        rt = 0.30,
        R2 = "B",
        rt2 = 0.55,
        stringsAsFactors = FALSE
      )
      manual <- acc_pdf_scalar(0.30, a) * acc_pdf_scalar(0.55, b)
      engine <- engine_density_or_mass(structure, params, data_df)
      check_row(
        "ranked_independent",
        "A_0.30_then_B_0.55",
        engine,
        manual,
        1e-9,
        "Ranked independent outcomes use joint exact times"
      )
    },

    ranked_chained_onset = function() {
      structure <- race_spec(n_outcomes = 2L) |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal", onset = after("a")) |>
        add_outcome("A", "a") |>
        add_outcome("B", "b") |>
        finalize_model()
      params <- c(
        a.m = log(0.25), a.s = 0.15, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.35), b.s = 0.20, b.q = 0.00, b.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      data_df <- data.frame(
        trial = 1L,
        R = "A",
        rt = 0.30,
        R2 = "B",
        rt2 = 0.60,
        stringsAsFactors = FALSE
      )
      manual <- acc_pdf_scalar(0.30, a) * acc_pdf_scalar(0.30, b)
      engine <- engine_density_or_mass(structure, params, data_df)
      check_row(
        "ranked_chained_onset",
        "A_0.30_then_B_0.60",
        engine,
        manual,
        2e-3,
        "Ranked after(a) uses observed prerequisite time"
      )
    },

    shared_trigger_conditioning = function() {
      structure <- race_spec() |>
        add_accumulator("s", "lognormal") |>
        add_accumulator("stop", "lognormal") |>
        add_accumulator("change", "lognormal") |>
        add_outcome("S", inhibit("s", by = "stop")) |>
        add_outcome("X", all_of("change", "stop")) |>
        add_trigger("tg", members = c("stop", "change"), q = 0.05, draw = "shared") |>
        finalize_model()
      params <- c(
        s.m = log(0.28), s.s = 0.12, s.q = 0.00, s.t0 = 0.00,
        stop.m = log(0.35), stop.s = 0.15, stop.t0 = 0.00,
        change.m = log(0.40), change.s = 0.18, change.t0 = 0.00,
        tg = 0.05
      )
      s <- acc_parts("s", params)
      stop <- acc_parts("stop", params)
      rt <- 0.30
      manual <- params[["tg"]] * acc_pdf_scalar(rt, s) +
        (1.0 - params[["tg"]]) * acc_pdf_scalar(rt, s) * acc_survival_scalar(rt, stop)
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "S", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "shared_trigger_conditioning",
        "S_rt_0.30",
        engine,
        manual,
        1e-9,
        "Shared trigger competitors conditioned on realized state"
      )
    },

    shared_gate_pair = function() {
      structure <- race_spec() |>
        add_accumulator("x1", "lognormal") |>
        add_accumulator("x2", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_outcome("RESP", all_of("x2", "gate")) |>
        add_outcome("NR_RAW", all_of("x1", "gate"), options = list(map_outcome_to = NA_character_)) |>
        finalize_model()
      params <- c(
        x1.m = log(0.32), x1.s = 0.18, x1.q = 0.00, x1.t0 = 0.00,
        x2.m = log(0.36), x2.s = 0.18, x2.q = 0.00, x2.t0 = 0.00,
        gate.m = log(0.24), gate.s = 0.14, gate.q = 0.00, gate.t0 = 0.00
      )
      x1 <- acc_parts("x1", params)
      x2 <- acc_parts("x2", params)
      gate <- acc_parts("gate", params)
      eval_x1 <- function(t) list(density = acc_pdf_scalar(t, x1), survival = acc_survival_scalar(t, x1))
      eval_x2 <- function(t) list(density = acc_pdf_scalar(t, x2), survival = acc_survival_scalar(t, x2))
      eval_gate <- function(t) list(density = acc_pdf_scalar(t, gate), survival = acc_survival_scalar(t, gate))

      resp_rt <- 0.42
      resp_manual <- shared_gate_pair_density_r(eval_x2, eval_x1, eval_gate, resp_rt)
      resp_engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "RESP", rt = resp_rt, stringsAsFactors = FALSE)
      )
      resp_row <- check_row(
        "shared_gate_pair",
        "RESP_rt_0.42",
        resp_engine,
        resp_manual,
        1e-8,
        "Composite max event with mapped-NA competitor"
      )

      resp_prob <- shared_gate_pair_probability_r(eval_x2, eval_x1, eval_gate, Inf)
      miss_manual <- 1.0 - resp_prob
      miss_engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = NA_character_, rt = NA_real_, stringsAsFactors = FALSE)
      )
      miss_row <- check_row(
        "shared_gate_pair",
        "NA_mass",
        miss_engine,
        miss_manual,
        2e-6,
        "Missing-response mass for mapped NA outcome"
      )

      rbind(resp_row, miss_row)
    },

    guarded_positive_mass_tie = function() {
      structure <- race_spec() |>
        add_accumulator("go_fast", "lognormal") |>
        add_accumulator("go_slow", "lognormal") |>
        add_accumulator("gate_shared", "lognormal") |>
        add_accumulator("stop_control", "lognormal") |>
        add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
        add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
        finalize_model()
      params <- c(
        go_fast.m = log(0.28), go_fast.s = 0.18, go_fast.q = 0.00, go_fast.t0 = 0.00,
        go_slow.m = log(0.34), go_slow.s = 0.18, go_slow.q = 0.00, go_slow.t0 = 0.00,
        gate_shared.m = log(0.30), gate_shared.s = 0.16, gate_shared.q = 0.00, gate_shared.t0 = 0.00,
        stop_control.m = log(0.27), stop_control.s = 0.15, stop_control.q = 0.00, stop_control.t0 = 0.00
      )
      rt <- 0.30
      f_fast <- function(x) dlnorm(x, params[["go_fast.m"]], params[["go_fast.s"]])
      F_fast <- function(x) plnorm(x, params[["go_fast.m"]], params[["go_fast.s"]])
      f_slow <- function(x) dlnorm(x, params[["go_slow.m"]], params[["go_slow.s"]])
      F_slow <- function(x) plnorm(x, params[["go_slow.m"]], params[["go_slow.s"]])
      f_gate <- function(x) dlnorm(x, params[["gate_shared.m"]], params[["gate_shared.s"]])
      F_gate <- function(x) plnorm(x, params[["gate_shared.m"]], params[["gate_shared.s"]])
      S_stop <- function(x) 1 - plnorm(x, params[["stop_control.m"]], params[["stop_control.s"]])

      manual <- S_stop(rt) * (
        f_fast(rt) * F_gate(rt) * (1 - F_slow(rt)) +
          f_gate(rt) * F_fast(rt) * (1 - F_slow(rt)) +
          f_gate(rt) * integrate_scalar(
            function(u) f_fast(u) * (F_slow(rt) - F_slow(u)),
            0.0,
            rt
          )
      )
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "Fast", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "guarded_positive_mass_tie",
        "Fast_rt_0.30",
        engine,
        manual,
        2e-8,
        "Guarded tie mass carried into observed density"
      )
    },

    shared_gate_three_way_tie = function() {
      structure <- race_spec() |>
        add_accumulator("x1", "lognormal") |>
        add_accumulator("x2", "lognormal") |>
        add_accumulator("x3", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_outcome("R1", all_of("x1", "gate")) |>
        add_outcome("R2", all_of("x2", "gate")) |>
        add_outcome("R3", all_of("x3", "gate")) |>
        finalize_model()
      params <- c(
        x1.m = log(0.31), x1.s = 0.16, x1.q = 0.00, x1.t0 = 0.00,
        x2.m = log(0.36), x2.s = 0.18, x2.q = 0.00, x2.t0 = 0.00,
        x3.m = log(0.41), x3.s = 0.17, x3.q = 0.00, x3.t0 = 0.00,
        gate.m = log(0.24), gate.s = 0.14, gate.q = 0.00, gate.t0 = 0.00
      )
      x1 <- acc_parts("x1", params)
      x2 <- acc_parts("x2", params)
      x3 <- acc_parts("x3", params)
      gate <- acc_parts("gate", params)

      rows <- list()
      for (rt in c(0.28, 0.42, 0.55)) {
        manual <- acc_pdf_scalar(rt, x1) *
          acc_cdf_scalar(rt, gate) *
          acc_survival_scalar(rt, x2) *
          acc_survival_scalar(rt, x3) +
          acc_pdf_scalar(rt, gate) * integrate_scalar(
            function(u) {
              acc_pdf_scalar(u, x1) *
                acc_survival_scalar(u, x2) *
                acc_survival_scalar(u, x3)
            },
            0.0,
            rt
          )
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "R1", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "shared_gate_three_way_tie",
          paste0("R1_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Three-way gated tie with unique observed branch"
        )
      }
      do.call(rbind, rows)
    },

    nested_guard_pair = function() {
      structure <- race_spec() |>
        add_accumulator("s1", "lognormal") |>
        add_accumulator("is", "lognormal") |>
        add_accumulator("s2", "lognormal") |>
        add_outcome("STOP", all_of("s1", inhibit("s2", by = "is"))) |>
        finalize_model()
      params <- c(
        s1.m = log(0.35), s1.s = 0.16, s1.q = 0.00, s1.t0 = 0.00,
        is.m = log(0.30), is.s = 0.14, is.q = 0.00, is.t0 = 0.00,
        s2.m = log(0.28), s2.s = 0.14, s2.q = 0.00, s2.t0 = 0.00
      )
      rt <- 0.42
      f_s1 <- function(x) dlnorm(x, params[["s1.m"]], params[["s1.s"]])
      F_s1 <- function(x) plnorm(x, params[["s1.m"]], params[["s1.s"]])
      f_s2 <- function(x) dlnorm(x, params[["s2.m"]], params[["s2.s"]])
      S_is <- function(x) 1 - plnorm(x, params[["is.m"]], params[["is.s"]])
      guard_cdf <- integrate_scalar(function(u) f_s2(u) * S_is(u), 0.0, rt)
      guard_density <- f_s2(rt) * S_is(rt)
      manual <- f_s1(rt) * guard_cdf + F_s1(rt) * guard_density
      engine <- engine_density_or_mass(
        structure,
        params,
        data.frame(trial = 1L, R = "STOP", rt = rt, stringsAsFactors = FALSE)
      )
      check_row(
        "nested_guard_pair",
        "STOP_rt_0.42",
        engine,
        manual,
        2e-3,
        "Nested all_of with guarded subexpression"
      )
    },

    deep_guard_chain = function() {
      structure <- race_spec() |>
        add_accumulator("plain", "lognormal") |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("c", "lognormal") |>
        add_accumulator("d", "lognormal") |>
        add_outcome("PLAIN", "plain") |>
        add_outcome("GUARD", inhibit("a", by = inhibit("b", by = inhibit("c", by = "d")))) |>
        finalize_model()
      params <- c(
        plain.m = log(0.42), plain.s = 0.18, plain.q = 0.00, plain.t0 = 0.00,
        a.m = log(0.34), a.s = 0.20, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.30), b.s = 0.20, b.q = 0.00, b.t0 = 0.00,
        c.m = log(0.28), c.s = 0.18, c.q = 0.00, c.t0 = 0.00,
        d.m = log(0.26), d.s = 0.16, d.q = 0.00, d.t0 = 0.00
      )
      plain <- acc_parts("plain", params)
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      c <- acc_parts("c", params)
      d <- acc_parts("d", params)

      g3_cdf <- function(t) {
        inhibit_cdf_scalar(
          function(u) acc_pdf_scalar(u, c),
          function(u) acc_cdf_scalar(u, d),
          t
        )
      }
      g2_cdf <- function(t) {
        inhibit_cdf_scalar(
          function(u) acc_pdf_scalar(u, b),
          g3_cdf,
          t
        )
      }
      g1_cdf <- function(t) {
        inhibit_cdf_scalar(
          function(u) acc_pdf_scalar(u, a),
          g2_cdf,
          t
        )
      }
      guard_density <- function(t) {
        inhibit_density_scalar(
          function(u) acc_pdf_scalar(u, a),
          g2_cdf,
          t
        )
      }
      rows <- list()
      for (rt in c(0.30, 0.42, 0.55)) {
        manual_plain <- acc_pdf_scalar(rt, plain) * max(0.0, 1.0 - g1_cdf(rt))
        engine_plain <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "PLAIN", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "deep_guard_chain",
          paste0("PLAIN_rt_", format(rt, nsmall = 2)),
          engine_plain,
          manual_plain,
          2e-3,
          "Deep nested guard chain against competitor"
        )

        manual_guard <- guard_density(rt) * acc_survival_scalar(rt, plain)
        engine_guard <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "GUARD", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "deep_guard_chain",
          paste0("GUARD_rt_", format(rt, nsmall = 2)),
          engine_guard,
          manual_guard,
          2e-3,
          "Deep nested guard chain observed branch"
        )
      }
      do.call(rbind, rows)
    },

    pooled_shared_gate_tie = function() {
      structure <- race_spec() |>
        add_accumulator("a1", "lognormal") |>
        add_accumulator("a2", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_pool("A_pool", c("a1", "a2"), k = 1L) |>
        add_outcome("A", all_of("A_pool", "gate")) |>
        add_outcome("B", all_of("b", "gate")) |>
        finalize_model()
      params <- c(
        a1.m = log(0.29), a1.s = 0.16, a1.q = 0.00, a1.t0 = 0.00,
        a2.m = log(0.35), a2.s = 0.17, a2.q = 0.00, a2.t0 = 0.00,
        b.m = log(0.38), b.s = 0.18, b.q = 0.00, b.t0 = 0.00,
        gate.m = log(0.25), gate.s = 0.14, gate.q = 0.00, gate.t0 = 0.00
      )
      a1 <- acc_parts("a1", params)
      a2 <- acc_parts("a2", params)
      b <- acc_parts("b", params)
      gate <- acc_parts("gate", params)

      eval_pool <- function(t) {
        list(
          density = acc_pdf_scalar(t, a1) * acc_survival_scalar(t, a2) +
            acc_pdf_scalar(t, a2) * acc_survival_scalar(t, a1),
          survival = acc_survival_scalar(t, a1) * acc_survival_scalar(t, a2)
        )
      }
      eval_b <- function(t) {
        list(
          density = acc_pdf_scalar(t, b),
          survival = acc_survival_scalar(t, b)
        )
      }
      eval_gate <- function(t) {
        list(
          density = acc_pdf_scalar(t, gate),
          survival = acc_survival_scalar(t, gate)
        )
      }

      rows <- list()
      for (rt in c(0.28, 0.44)) {
        manual <- shared_gate_pair_density_r(eval_pool, eval_b, eval_gate, rt)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "A", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "pooled_shared_gate_tie",
          paste0("A_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Shared-gate tie with pooled readiness branch"
        )
      }
      do.call(rbind, rows)
    },

    pooled_guarded_shared_gate_tie = function() {
      structure <- race_spec() |>
        add_accumulator("a1", "lognormal") |>
        add_accumulator("a2", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_accumulator("stop", "lognormal") |>
        add_pool("A_pool", c("a1", "a2"), k = 1L) |>
        add_outcome("A", inhibit(all_of("A_pool", "gate"), by = "stop")) |>
        add_outcome("B", all_of("b", "gate")) |>
        finalize_model()
      params <- c(
        a1.m = log(0.27), a1.s = 0.15, a1.q = 0.00, a1.t0 = 0.00,
        a2.m = log(0.33), a2.s = 0.16, a2.q = 0.00, a2.t0 = 0.00,
        b.m = log(0.36), b.s = 0.18, b.q = 0.00, b.t0 = 0.00,
        gate.m = log(0.24), gate.s = 0.13, gate.q = 0.00, gate.t0 = 0.00,
        stop.m = log(0.31), stop.s = 0.15, stop.q = 0.00, stop.t0 = 0.00
      )
      a1 <- acc_parts("a1", params)
      a2 <- acc_parts("a2", params)
      b <- acc_parts("b", params)
      gate <- acc_parts("gate", params)
      stop <- acc_parts("stop", params)

      eval_pool <- function(t) {
        list(
          density = acc_pdf_scalar(t, a1) * acc_survival_scalar(t, a2) +
            acc_pdf_scalar(t, a2) * acc_survival_scalar(t, a1),
          survival = acc_survival_scalar(t, a1) * acc_survival_scalar(t, a2)
        )
      }
      eval_b <- function(t) {
        list(
          density = acc_pdf_scalar(t, b),
          survival = acc_survival_scalar(t, b)
        )
      }
      eval_gate <- function(t) {
        list(
          density = acc_pdf_scalar(t, gate),
          survival = acc_survival_scalar(t, gate)
        )
      }

      rows <- list()
      for (rt in c(0.29, 0.41)) {
        manual <- acc_survival_scalar(rt, stop) *
          shared_gate_pair_density_r(eval_pool, eval_b, eval_gate, rt)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "A", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "pooled_guarded_shared_gate_tie",
          paste0("A_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Guarded shared-gate tie with pooled readiness branch"
        )
      }
      do.call(rbind, rows)
    },

    overlapping_composite_competitors = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_accumulator("d", "lognormal") |>
        add_outcome("C1", all_of("a", "gate")) |>
        add_outcome("C2", all_of("b", "gate")) |>
        add_outcome("D", "d") |>
        finalize_model()
      params <- c(
        a.m = log(0.28), a.s = 0.16, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.34), b.s = 0.17, b.q = 0.00, b.t0 = 0.00,
        gate.m = log(0.23), gate.s = 0.14, gate.q = 0.00, gate.t0 = 0.00,
        d.m = log(0.47), d.s = 0.19, d.q = 0.00, d.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      gate <- acc_parts("gate", params)
      d <- acc_parts("d", params)

      rows <- list()
      for (rt in c(0.26, 0.40, 0.58)) {
        manual <- acc_pdf_scalar(rt, d) * (
          acc_survival_scalar(rt, gate) +
            acc_cdf_scalar(rt, gate) *
            acc_survival_scalar(rt, a) *
            acc_survival_scalar(rt, b)
        )
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "D", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "overlapping_composite_competitors",
          paste0("D_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Simple target against two overlapping shared-gate competitors"
        )
      }
      do.call(rbind, rows)
    },

    guarded_overlapping_competitors = function() {
      structure <- race_spec() |>
        add_accumulator("a", "lognormal") |>
        add_accumulator("b", "lognormal") |>
        add_accumulator("stop", "lognormal") |>
        add_accumulator("d", "lognormal") |>
        add_outcome("C1", inhibit("a", by = "stop")) |>
        add_outcome("C2", inhibit("b", by = "stop")) |>
        add_outcome("D", "d") |>
        finalize_model()
      params <- c(
        a.m = log(0.29), a.s = 0.17, a.q = 0.00, a.t0 = 0.00,
        b.m = log(0.35), b.s = 0.16, b.q = 0.00, b.t0 = 0.00,
        stop.m = log(0.24), stop.s = 0.14, stop.q = 0.00, stop.t0 = 0.00,
        d.m = log(0.46), d.s = 0.19, d.q = 0.00, d.t0 = 0.00
      )
      a <- acc_parts("a", params)
      b <- acc_parts("b", params)
      stop <- acc_parts("stop", params)
      d <- acc_parts("d", params)

      rows <- list()
      for (rt in c(0.27, 0.43)) {
        manual <- acc_pdf_scalar(rt, d) * (
          integrate_scalar(
            function(u) {
              acc_pdf_scalar(u, stop) *
                acc_survival_scalar(u, a) *
                acc_survival_scalar(u, b)
            },
            0.0,
            rt
          ) +
            acc_survival_scalar(rt, stop) *
            acc_survival_scalar(rt, a) *
            acc_survival_scalar(rt, b)
        )
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "D", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "guarded_overlapping_competitors",
          paste0("D_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Simple target against two overlapping guarded competitors"
        )
      }
      do.call(rbind, rows)
    },

    shared_gate_four_way_tie = function() {
      structure <- race_spec() |>
        add_accumulator("x1", "lognormal") |>
        add_accumulator("x2", "lognormal") |>
        add_accumulator("x3", "lognormal") |>
        add_accumulator("x4", "lognormal") |>
        add_accumulator("gate", "lognormal") |>
        add_outcome("R1", all_of("x1", "gate")) |>
        add_outcome("R2", all_of("x2", "gate")) |>
        add_outcome("R3", all_of("x3", "gate")) |>
        add_outcome("R4", all_of("x4", "gate")) |>
        finalize_model()
      params <- c(
        x1.m = log(0.30), x1.s = 0.15, x1.q = 0.00, x1.t0 = 0.00,
        x2.m = log(0.34), x2.s = 0.16, x2.q = 0.00, x2.t0 = 0.00,
        x3.m = log(0.39), x3.s = 0.17, x3.q = 0.00, x3.t0 = 0.00,
        x4.m = log(0.44), x4.s = 0.18, x4.q = 0.00, x4.t0 = 0.00,
        gate.m = log(0.22), gate.s = 0.13, gate.q = 0.00, gate.t0 = 0.00
      )
      x1 <- acc_parts("x1", params)
      x2 <- acc_parts("x2", params)
      x3 <- acc_parts("x3", params)
      x4 <- acc_parts("x4", params)
      gate <- acc_parts("gate", params)

      rows <- list()
      for (rt in c(0.27, 0.41, 0.56)) {
        manual <- acc_pdf_scalar(rt, x1) *
          acc_cdf_scalar(rt, gate) *
          acc_survival_scalar(rt, x2) *
          acc_survival_scalar(rt, x3) *
          acc_survival_scalar(rt, x4) +
          acc_pdf_scalar(rt, gate) * integrate_scalar(
            function(u) {
              acc_pdf_scalar(u, x1) *
                acc_survival_scalar(u, x2) *
                acc_survival_scalar(u, x3) *
                acc_survival_scalar(u, x4)
            },
            0.0,
            rt
          )
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "R1", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "shared_gate_four_way_tie",
          paste0("R1_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Four-way gated tie with unique observed branch"
        )
      }
      do.call(rbind, rows)
    },

    none_of_conjunction = function() {
      structure <- race_spec() |>
        add_accumulator("go", "lognormal") |>
        add_accumulator("stop", "lognormal") |>
        add_accumulator("other", "lognormal") |>
        add_outcome("RESPOND", all_of("go", none_of("stop"))) |>
        add_outcome("OTHER", "other") |>
        finalize_model()
      params <- c(
        go.m = log(0.31), go.s = 0.15, go.q = 0.00, go.t0 = 0.00,
        stop.m = log(0.27), stop.s = 0.13, stop.q = 0.00, stop.t0 = 0.00,
        other.m = log(0.45), other.s = 0.17, other.q = 0.00, other.t0 = 0.00
      )
      go <- acc_parts("go", params)
      stop <- acc_parts("stop", params)
      other <- acc_parts("other", params)

      rows <- list()
      for (rt in c(0.29, 0.41)) {
        manual <- acc_pdf_scalar(rt, go) *
          acc_survival_scalar(rt, stop) *
          acc_survival_scalar(rt, other)
        engine <- engine_density_or_mass(
          structure,
          params,
          data.frame(trial = 1L, R = "RESPOND", rt = rt, stringsAsFactors = FALSE)
        )
        rows[[length(rows) + 1L]] <- check_row(
          "none_of_conjunction",
          paste0("RESPOND_rt_", format(rt, nsmall = 2)),
          engine,
          manual,
          2e-3,
          "Conjunction with none_of as an absence condition"
        )
      }
      do.call(rbind, rows)
    }
  )
}
