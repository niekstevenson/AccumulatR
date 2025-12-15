testthat::test_that("selected examples agree across simulate/probability/likelihood", {
  # Helper to compute manual log-likelihood from response probabilities
  manual_loglik <- function(probs, outcomes) {
    labels <- as.character(outcomes)
    labels[is.na(labels)] <- "NA"
    if (!all(labels %in% names(probs))) {
      missing <- setdiff(unique(labels), names(probs))
      stop("Missing probabilities for outcomes: ", paste(missing, collapse = ", "))
    }
    sum(log(probs[labels]))
  }

  # Model definitions (copied from dev/examples/new_API.R) ---------------------
  example_1_simple <- function() {
    race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("go2", "lognormal") |>
      add_pool("R1", "go1") |>
      add_pool("R2", "go2") |>
      add_outcome("R1", "R1") |>
      add_outcome("R2", "R2") |>
      finalize_model()
  }
  params_example_1_simple <- c(
    go1.meanlog = log(0.30),
    go1.sdlog = 0.18,
    go2.meanlog = log(0.32),
    go2.sdlog = 0.18
  )

  example_2_stop_mixture <- function() {
    race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("stop", "exgauss", onset = 0.20) |>
      add_accumulator("go2", "lognormal", onset = 0.20) |>
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
        components = list(component("go_only", weight = 0.5),
                          component("go_stop", weight = 0.5))
      )) |>
      finalize_model()
  }
  params_example_2_stop_mixture <- c(
    go1.meanlog = log(0.35),
    go1.sdlog = 0.2,
    stop.mu = 0.1,
    stop.sigma = 0.04,
    stop.tau = 0.1,
    go2.meanlog = log(0.60),
    go2.sdlog = 0.18
  )

  example_3_stop_na <- function() {
    race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_accumulator("stop", "lognormal", onset = 0.15) |>
      add_pool("L", "go_left") |>
      add_pool("R", "go_right") |>
      add_pool("STOP", "stop") |>
      add_outcome("Left", "L") |>
      add_outcome("Right", "R") |>
      add_outcome("STOP", "STOP", options = list(map_outcome_to = NA_character_)) |>
      finalize_model()
  }
  params_example_3_stop_na <- c(
    go_left.meanlog = log(0.30),
    go_left.sdlog = 0.20,
    go_right.meanlog = log(0.32),
    go_right.sdlog = 0.20,
    stop.meanlog = log(0.15),
    stop.sdlog = 0.18
  )

  example_5_timeout_guess <- function() {
    race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_accumulator("timeout", "lognormal", onset = 0.05) |>
      add_pool("L", "go_left") |>
      add_pool("R", "go_right") |>
      add_pool("TO", "timeout") |>
      add_outcome("Left", "L") |>
      add_outcome("Right", "R") |>
      add_outcome("TIMEOUT", "TO", options = list(
        guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
      )) |>
      finalize_model()
  }
  params_example_5_timeout_guess <- c(
    go_left.meanlog = log(0.30),
    go_left.sdlog = 0.18,
    go_right.meanlog = log(0.325),
    go_right.sdlog = 0.18,
    timeout.meanlog = log(0.25),
    timeout.sdlog = 0.10
  )

  example_6_dual_path <- function() {
    race_spec() |>
      add_accumulator("acc_taskA", "lognormal") |>
      add_accumulator("acc_taskB", "lognormal") |>
      add_accumulator("acc_gateC", "lognormal") |>
      add_pool("TaskA", "acc_taskA") |>
      add_pool("TaskB", "acc_taskB") |>
      add_pool("GateC", "acc_gateC") |>
      add_outcome("Outcome_via_A", all_of("TaskA", "GateC")) |>
      add_outcome("Outcome_via_B", all_of("TaskB", "GateC")) |>
      finalize_model()
  }
  params_example_6_dual_path <- c(
    acc_taskA.meanlog = log(0.28),
    acc_taskA.sdlog = 0.18,
    acc_taskB.meanlog = log(0.32),
    acc_taskB.sdlog = 0.18,
    acc_gateC.meanlog = log(0.30),
    acc_gateC.sdlog = 0.18
  )

  example_7_mixture <- function() {
    race_spec() |>
      add_accumulator("target_fast", "lognormal") |>
      add_accumulator("target_slow", "lognormal") |>
      add_accumulator("competitor", "lognormal") |>
      add_pool("TARGET", c("target_fast", "target_slow")) |>
      add_pool("COMP", "competitor") |>
      add_outcome("R1", "TARGET") |>
      add_outcome("R2", "COMP") |>
      add_group("component:fast", members = c("target_fast", "competitor"),
                attrs = list(component = "fast")) |>
      add_group("component:slow", members = c("target_slow", "competitor"),
                attrs = list(component = "slow")) |>
      set_metadata(mixture = list(
        mode = "sample",
        reference = "slow",
        components = list(
          component("fast", attrs = list(weight_param = "p_fast")),
          component("slow")
        )
      )) |>
      finalize_model()
  }
  params_example_7_mixture <- c(
    target_fast.meanlog = log(0.25),
    target_fast.sdlog = 0.15,
    target_slow.meanlog = log(0.45),
    target_slow.sdlog = 0.20,
    competitor.meanlog = log(0.35),
    competitor.sdlog = 0.18,
    p_fast = 0.2
  )

  example_10_exclusion <- function() {
    race_spec() |>
      add_accumulator("R1_acc", "lognormal") |>
      add_accumulator("R2_acc", "lognormal") |>
      add_accumulator("X_acc", "lognormal") |>
      add_pool("R1", "R1_acc") |>
      add_pool("R2", "R2_acc") |>
      add_pool("X", "X_acc") |>
      add_outcome("R1", inhibit("R1", by = "X")) |>
      add_outcome("R2", "R2") |>
      finalize_model()
  }
  params_example_10_exclusion <- c(
    R1_acc.meanlog = log(0.35),
    R1_acc.sdlog = 0.18,
    R2_acc.meanlog = log(0.45),
    R2_acc.sdlog = 0.18,
    X_acc.meanlog = log(0.35),
    X_acc.sdlog = 0.18
  )

  example_16_guard_tie_simple <- function() {
    race_spec() |>
      add_accumulator("go_fast", "lognormal") |>
      add_accumulator("go_slow", "lognormal") |>
      add_accumulator("gate_shared", "lognormal") |>
      add_accumulator("stop_control", "lognormal") |>
      add_pool("FAST", "go_fast") |>
      add_pool("SLOW", "go_slow") |>
      add_pool("GATE", "gate_shared") |>
      add_pool("STOP", "stop_control") |>
      add_outcome("Fast", inhibit(all_of("FAST", "GATE"), by = "STOP")) |>
      add_outcome("Slow", all_of("SLOW", "GATE")) |>
      finalize_model()
  }
  params_example_16_guard_tie_simple <- c(
    go_fast.meanlog = log(0.28),
    go_fast.sdlog = 0.18,
    go_slow.meanlog = log(0.34),
    go_slow.sdlog = 0.18,
    gate_shared.meanlog = log(0.30),
    gate_shared.sdlog = 0.16,
    stop_control.meanlog = log(0.27),
    stop_control.sdlog = 0.15
  )

  example_21_simple_q <- function() {
    race_spec() |>
      add_accumulator("go1", "lognormal") |>
      add_accumulator("go2", "lognormal") |>
      add_pool("R1", "go1") |>
      add_pool("R2", "go2") |>
      add_outcome("R1", "R1") |>
      add_outcome("R2", "R2") |>
    add_group(
      "shared_trigger",
      members = c("go1", "go2"),
      attrs = trigger(id = "go_trigger",
                      q = 0.10,
                      param = "go_trigger_q",
                      draw = "shared")
    ) |>
    finalize_model()
  }
  params_example_21_simple_q <- c(
    go1.meanlog = log(0.30),
    go1.sdlog   = 0.18,
    go2.meanlog = log(0.32),
    go2.sdlog   = 0.18,
    go_trigger_q = 0.10
  )

  example_22_shared_q <- function() {
    race_spec() |>
      add_accumulator("go_left", "lognormal") |>
      add_accumulator("go_right", "lognormal") |>
      add_pool("L", "go_left") |>
      add_pool("R", "go_right") |>
      add_outcome("Left", "L") |>
      add_outcome("Right", "R") |>
      add_group("par:q_shared", members = c("go_left", "go_right"),
                attrs = trigger(q = 0.10, param = "q_shared", draw = "independent")) |>
      finalize_model()
  }
  params_example_22_shared_q <- c(
    go_left.meanlog = log(0.30),
    go_left.sdlog = 0.18,
    go_right.meanlog = log(0.32),
    go_right.sdlog = 0.18,
    q_shared = 0.10
  )

  # Bundle target models -------------------------------------------------------
  models <- list(
    example_1_simple = list(spec = example_1_simple, params = params_example_1_simple),
    example_2_stop_mixture = list(spec = example_2_stop_mixture, params = params_example_2_stop_mixture),
    example_3_stop_na = list(spec = example_3_stop_na, params = params_example_3_stop_na),
    example_5_timeout_guess = list(spec = example_5_timeout_guess, params = params_example_5_timeout_guess),
    example_6_dual_path = list(spec = example_6_dual_path, params = params_example_6_dual_path),
    example_7_mixture = list(spec = example_7_mixture, params = params_example_7_mixture),
    example_10_exclusion = list(spec = example_10_exclusion, params = params_example_10_exclusion),
    example_16_guard_tie_simple = list(spec = example_16_guard_tie_simple, params = params_example_16_guard_tie_simple),
    example_21_simple_q = list(spec = example_21_simple_q, params = params_example_22_simple_q),
    example_22_shared_q = list(spec = example_22_shared_q, params = params_example_23_shared_q)
  )

  # Test each model -----------------------------------------------------------
  for (name in names(models)) {
    mod <- models[[name]]
    spec_obj <- mod$spec()
    structure <- finalize_model(spec_obj)
    params_vec <- mod$params

    params_df <- build_params_df(spec_obj, params_vec, n_trials = 500)
    data_df <- simulate(structure, params_df, seed = 123)

    probs <- response_probabilities(
      structure,
      build_params_df(spec_obj, params_vec, n_trials = 1L),
      include_na = TRUE
    )
    testthat::expect_equal(sum(probs), 1, tolerance = 5e-4, info = name)

    # Compare analytic probabilities to empirical outcome frequencies
    emp <- prop.table(table(data_df$outcome, useNA = "ifany"))
    emp_names <- names(emp)
    emp_names[is.na(emp_names)] <- "NA"
    emp <- as.numeric(emp)
    names(emp) <- emp_names
    common <- intersect(names(probs), names(emp))
    testthat::expect_true(length(common) > 0, info = paste(name, "no common outcomes"))
    testthat::expect_true(max(abs(probs[common] - emp[common])) < 0.35, info = name)

    # Likelihood equivalence across argument forms
    ctx <- build_likelihood_context(structure, params_df, data_df)
    ll_list <- as.numeric(log_likelihood(ctx, list(params_df)))
    ll_df <- as.numeric(log_likelihood(ctx, params_df))
    testthat::expect_equal(ll_list, ll_df, tolerance = 1e-8, info = name)
    testthat::expect_true(is.finite(ll_list), info = name)
  }
})
