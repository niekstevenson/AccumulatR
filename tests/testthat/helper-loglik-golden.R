`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

loglik_golden_case_registry <- function() {
  list(
    example_1_simple = list(
      build = function() {
        race_spec() |>
          add_accumulator("go1", "lognormal") |>
          add_accumulator("go2", "lognormal") |>
          add_outcome("R1", "go1") |>
          add_outcome("R2", "go2") |>
          finalize_model()
      },
      params = c(
        go1.m = log(0.30),
        go1.s = 0.18,
        go2.m = log(0.32),
        go2.s = 0.18
      ),
      n_trials = 12L,
      seed = 1001L,
      data_cols = c("trial", "R", "rt")
    ),
    example_2_stop_mixture = list(
      build = function() {
        race_spec() |>
          add_accumulator("go1", "lognormal") |>
          add_accumulator("stop", "exgauss", onset = 0.20) |>
          add_accumulator("go2", "lognormal", onset = 0.20) |>
          add_outcome("R1", inhibit("go1", by = "stop")) |>
          add_outcome("R2", all_of("go2", "stop")) |>
          add_component("go_only", members = c("go1"), weight = 0.5) |>
          add_component("go_stop", members = c("go1", "stop", "go2"), weight = 0.5) |>
          set_mixture_options(mode = "fixed") |>
          finalize_model()
      },
      params = c(
        go1.m = log(0.35),
        go1.s = 0.20,
        stop.mu = 0.10,
        stop.sigma = 0.04,
        stop.tau = 0.10,
        go2.m = log(0.60),
        go2.s = 0.18
      ),
      n_trials = 12L,
      seed = 1002L,
      simulate_args = list(keep_component = TRUE),
      data_cols = c("trial", "component", "R", "rt")
    ),
    example_3_stop_na = list(
      build = function() {
        race_spec() |>
          add_accumulator("go_left", "lognormal") |>
          add_accumulator("go_right", "lognormal") |>
          add_accumulator("stop", "lognormal", onset = 0.15) |>
          add_outcome("Left", "go_left") |>
          add_outcome("Right", "go_right") |>
          add_outcome("STOP", "stop", options = list(map_outcome_to = NA_character_)) |>
          finalize_model()
      },
      params = c(
        go_left.m = log(0.30),
        go_left.s = 0.20,
        go_right.m = log(0.32),
        go_right.s = 0.20,
        stop.m = log(0.15),
        stop.s = 0.18
      ),
      n_trials = 12L,
      seed = 1003L,
      data_cols = c("trial", "R", "rt")
    ),
    example_5_timeout_guess = list(
      build = function() {
        race_spec() |>
          add_accumulator("go_left", "lognormal") |>
          add_accumulator("go_right", "lognormal") |>
          add_accumulator("timeout", "lognormal", onset = 0.05) |>
          add_outcome("Left", "go_left") |>
          add_outcome("Right", "go_right") |>
          add_outcome("TIMEOUT", "timeout", options = list(
            guess = list(labels = c("Left", "Right"), weights = c(0.2, 0.8), rt_policy = "keep")
          )) |>
          finalize_model()
      },
      params = c(
        go_left.m = log(0.30),
        go_left.s = 0.18,
        go_right.m = log(0.325),
        go_right.s = 0.18,
        timeout.m = log(0.25),
        timeout.s = 0.10
      ),
      n_trials = 12L,
      seed = 1004L,
      data_cols = c("trial", "R", "rt")
    ),
    example_6_dual_path = list(
      build = function() {
        race_spec() |>
          add_accumulator("acc_taskA", "lognormal") |>
          add_accumulator("acc_taskB", "lognormal") |>
          add_accumulator("acc_gateC", "lognormal") |>
          add_outcome("Outcome_via_A", all_of("acc_taskA", "acc_gateC")) |>
          add_outcome("Outcome_via_B", all_of("acc_taskB", "acc_gateC")) |>
          finalize_model()
      },
      params = c(
        acc_taskA.m = log(0.28),
        acc_taskA.s = 0.18,
        acc_taskB.m = log(0.32),
        acc_taskB.s = 0.18,
        acc_gateC.m = log(0.30),
        acc_gateC.s = 0.18
      ),
      n_trials = 10L,
      seed = 1005L,
      data_cols = c("trial", "R", "rt")
    ),
    example_7_mixture = list(
      build = function() {
        race_spec() |>
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
      },
      params = c(
        target_fast.m = log(0.25),
        target_fast.s = 0.15,
        target_slow.m = log(0.45),
        target_slow.s = 0.20,
        competitor.m = log(0.35),
        competitor.s = 0.18,
        p_fast = 0.20
      ),
      n_trials = 12L,
      seed = 1006L,
      data_cols = c("trial", "R", "rt")
    ),
    example_10_exclusion = list(
      build = function() {
        race_spec() |>
          add_accumulator("R1_acc", "lognormal") |>
          add_accumulator("R2_acc", "lognormal") |>
          add_accumulator("X_acc", "lognormal") |>
          add_outcome("R1", inhibit("R1_acc", by = "X_acc")) |>
          add_outcome("R2", "R2_acc") |>
          finalize_model()
      },
      params = c(
        R1_acc.m = log(0.35),
        R1_acc.s = 0.18,
        R2_acc.m = log(0.45),
        R2_acc.s = 0.18,
        X_acc.m = log(0.35),
        X_acc.s = 0.18
      ),
      n_trials = 10L,
      seed = 1007L,
      data_cols = c("trial", "R", "rt")
    ),
    example_16_guard_tie_simple = list(
      build = function() {
        race_spec() |>
          add_accumulator("go_fast", "lognormal") |>
          add_accumulator("go_slow", "lognormal") |>
          add_accumulator("gate_shared", "lognormal") |>
          add_accumulator("stop_control", "lognormal") |>
          add_outcome("Fast", inhibit(all_of("go_fast", "gate_shared"), by = "stop_control")) |>
          add_outcome("Slow", all_of("go_slow", "gate_shared")) |>
          finalize_model()
      },
      params = c(
        go_fast.m = log(0.28),
        go_fast.s = 0.18,
        go_slow.m = log(0.34),
        go_slow.s = 0.18,
        gate_shared.m = log(0.30),
        gate_shared.s = 0.16,
        stop_control.m = log(0.27),
        stop_control.s = 0.15
      ),
      n_trials = 10L,
      seed = 1008L,
      data_cols = c("trial", "R", "rt")
    ),
    example_21_simple_q = list(
      build = function() {
        race_spec() |>
          add_accumulator("go1", "lognormal") |>
          add_accumulator("go2", "lognormal") |>
          add_outcome("R1", "go1") |>
          add_outcome("R2", "go2") |>
          add_trigger("shared_trigger", members = c("go1", "go2"), q = 0.10, draw = "shared") |>
          finalize_model()
      },
      params = c(
        go1.m = log(0.30),
        go1.s = 0.18,
        go2.m = log(0.32),
        go2.s = 0.18,
        shared_trigger = 0.10
      ),
      n_trials = 12L,
      seed = 1009L,
      data_cols = c("trial", "R", "rt")
    ),
    example_22_shared_q = list(
      build = function() {
        race_spec() |>
          add_accumulator("go_left", "lognormal") |>
          add_accumulator("go_right", "lognormal") |>
          add_outcome("Left", "go_left") |>
          add_outcome("Right", "go_right") |>
          add_trigger("q_shared", members = c("go_left", "go_right"), q = 0.10, draw = "independent") |>
          finalize_model()
      },
      params = c(
        go_left.m = log(0.30),
        go_left.s = 0.18,
        go_right.m = log(0.32),
        go_right.s = 0.18,
        q_shared = 0.10
      ),
      n_trials = 12L,
      seed = 1010L,
      data_cols = c("trial", "R", "rt")
    ),
    example_23_ranked_chain = list(
      build = function() {
        race_spec(n_outcomes = 2L) |>
          add_accumulator("a", "lognormal") |>
          add_accumulator("b", "lognormal", onset = after("a")) |>
          add_outcome("A", "a") |>
          add_outcome("B", "b") |>
          finalize_model()
      },
      params = c(
        a.m = log(0.30),
        a.s = 0.16,
        b.m = log(0.22),
        b.s = 0.16
      ),
      n_trials = 10L,
      seed = 1011L,
      data_cols = c("trial", "R", "rt", "R2", "rt2")
    )
  )
}

materialize_loglik_golden_case <- function(def) {
  structure <- def$build()
  params_sim <- build_param_matrix(structure, def$params, n_trials = def$n_trials)
  sim_args <- c(
    list(structure, params_sim, seed = def$seed),
    def$simulate_args %||% list()
  )
  sim <- do.call(simulate, sim_args)
  data_df <- as.data.frame(sim[, def$data_cols, drop = FALSE], stringsAsFactors = FALSE)
  list(
    structure = structure,
    data_df = data_df,
    params = def$params,
    seed = def$seed,
    tolerance = def$tolerance %||% 1e-4
  )
}

evaluate_loglik_golden_case <- function(structure, data_df, params) {
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  params_df <- build_param_matrix(structure$model_spec, params, trial_df = prepared)
  trial_ids <- unique(prepared$trial)
  list(
    total_loglik = as.numeric(log_likelihood(ctx, prepared, params_df)),
    trial_loglik = vapply(trial_ids, function(trial_id) {
      trial_df <- data_df[data_df$trial == trial_id, , drop = FALSE]
      trial_prepared <- prepare_data(structure, trial_df)
      trial_params_df <- build_param_matrix(
        structure$model_spec,
        params,
        trial_df = trial_prepared
      )
      as.numeric(log_likelihood(
        ctx,
        trial_prepared,
        trial_params_df
      ))
    }, numeric(1))
  )
}

build_loglik_golden_fixture <- function() {
  registry <- loglik_golden_case_registry()
  cases <- lapply(names(registry), function(label) {
    materialized <- materialize_loglik_golden_case(registry[[label]])
    values <- evaluate_loglik_golden_case(
      materialized$structure,
      materialized$data_df,
      materialized$params
    )
    list(
      label = label,
      seed = materialized$seed,
      tolerance = materialized$tolerance,
      data_df = materialized$data_df,
      params = materialized$params,
      total_loglik = values$total_loglik,
      trial_loglik = values$trial_loglik
    )
  })
  list(
    version = 1L,
    created_at = as.character(Sys.time()),
    cases = cases
  )
}
