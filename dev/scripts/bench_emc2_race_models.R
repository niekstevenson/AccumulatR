Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")

library(AccumulatR)
library(EMC2)

sizes <- c(100L, 1000L)
samples <- 7L
min_ll <- log(1e-10)
particle_offsets <- seq(-0.08, 0.08, length.out = 16L)
particle_counts <- c(`100` = 20000L, `1000` = 2000L)

models <- list(
  LBA = list(
    emc_model = EMC2::LBA,
    acc_dist = "lba",
    distinct = "v",
    source = c(t0 = "t0", p1 = "v", p2 = "B", p3 = "A", p4 = "sv"),
    acc_pars = c(A.v = 2.2, B.v = 2.6, B = 1.1, A = 0.3,
                 sv = 1, t0 = 0.15),
    emc_theta = c(v_lRA = 2.2, v_lRB = 2.6, sv = log(1),
                  B = log(0.8), A = log(0.3), t0 = log(0.15)),
    acc_theta = c(v_lRA = 2.2, v_lRB = 2.6, sv = log(1),
                  B = log(1.1), A = log(0.3), t0 = log(0.15))
  ),
  RDM = list(
    emc_model = EMC2::RDM,
    acc_dist = "rdm",
    distinct = "v",
    source = c(t0 = "t0", p1 = "v", p2 = "B", p3 = "A", p4 = "s"),
    acc_pars = c(A.v = 1.8, B.v = 2.2, B = 0.7, A = 0.2,
                 s = 1, t0 = 0.15),
    emc_theta = c(v_lRA = log(1.8), v_lRB = log(2.2), B = log(0.7),
                  A = log(0.2), t0 = log(0.15), s = log(1)),
    acc_theta = c(v_lRA = log(1.8), v_lRB = log(2.2), B = log(0.7),
                  A = log(0.2), t0 = log(0.15), s = log(1))
  ),
  LNR = list(
    emc_model = EMC2::LNR,
    acc_dist = "lognormal",
    distinct = "m",
    source = c(t0 = "t0", p1 = "m", p2 = "s"),
    acc_pars = c(A.m = -0.55, B.m = -0.35, s = 0.35, t0 = 0.15),
    emc_theta = c(m_lRA = -0.55, m_lRB = -0.35,
                  s = log(0.35), t0 = log(0.15)),
    acc_theta = c(m_lRA = -0.55, m_lRB = -0.35,
                  s = log(0.35), t0 = log(0.15))
  )
)

cpp_likelihood <- function(data, model, accumulatr_context = NULL) {
  designs <- EMC2:::get_designs_expanded(data, model)
  constants <- attr(data, "constants")
  if (is.null(constants)) constants <- NA_real_
  function(particles) {
    EMC2:::calc_ll(
      particles, data, constants, designs, model$c_name,
      model$bound, model$transform, model$pre_transform,
      names(model$p_types), min_ll, model$trend, accumulatr_context
    )
  }
}

trial_data <- function(n_trials) {
  data.frame(
    subjects = factor(rep(1L, n_trials)),
    trials = seq_len(n_trials),
    R = factor(rep(c("A", "B"), length.out = n_trials),
               levels = c("A", "B")),
    rt = seq(0.35, 1.35, length.out = n_trials)
  )
}

make_case <- function(definition, n_trials) {
  data <- trial_data(n_trials)
  emc_model <- definition$emc_model()
  formulas <- setNames(
    lapply(names(emc_model$p_types), function(parameter) {
      as.formula(paste(parameter, "~ 1"))
    }),
    names(emc_model$p_types)
  )
  formulas[[definition$distinct]] <- reformulate(
    "lR", response = definition$distinct, intercept = FALSE
  )
  design <- EMC2::design(
    data = data,
    model = definition$emc_model,
    formula = formulas,
    report_p_vector = FALSE
  )
  emc <- EMC2::make_emc(
    data, design, type = "single", compress = FALSE, rt_resolution = 1e-12
  )
  emc_data <- emc[[1]]$data[[1]]
  emc_model <- emc[[1]]$model()
  theta_names <- attr(emc_data, "p_names")
  theta <- list(
    EMC2 = matrix(definition$emc_theta[theta_names], 1L,
                  dimnames = list(NULL, theta_names)),
    AccumulatR = matrix(definition$acc_theta[theta_names], 1L,
                        dimnames = list(NULL, theta_names))
  )

  acc_model_spec <- AccumulatR::race_spec() |>
    AccumulatR::add_accumulator("A", definition$acc_dist) |>
    AccumulatR::add_accumulator("B", definition$acc_dist) |>
    AccumulatR::add_outcome("A", "A") |>
    AccumulatR::add_outcome("B", "B") |>
    AccumulatR::set_parameters(
      separate = setNames(list(TRUE), definition$distinct)
    ) |>
    AccumulatR::finalize_model()
  acc_data <- AccumulatR::prepare_data(
    acc_model_spec,
    data[setdiff(names(data), "subjects")],
    compress = FALSE
  )
  template <- AccumulatR::build_param_matrix(
    acc_model_spec, definition$acc_pars, n_trials = n_trials
  )
  source_names <- matrix(
    NA_character_, nrow(template), ncol(template), dimnames = dimnames(template)
  )
  for (target in names(definition$source)) {
    source_names[, target] <- definition$source[[target]]
  }
  attr(acc_data, "designs") <- attr(emc_data, "designs")
  attr(acc_data, "constants") <- attr(emc_data, "constants")
  attr(acc_data, "expand") <- seq_len(n_trials)
  acc_context <- list(
    native = AccumulatR::make_context(acc_model_spec)$cpp,
    bridge = list(defaults = template, source_names = source_names),
    trial_counts = rep.int(1L, n_trials)
  )
  acc_model <- emc_model
  acc_model$c_name <- "AccumulatR"

  list(
    EMC2 = cpp_likelihood(emc_data, emc_model),
    AccumulatR = cpp_likelihood(acc_data, acc_model, acc_context),
    theta = theta
  )
}

make_particles <- function(theta, n, parameter) {
  particles <- theta[rep.int(1L, n), , drop = FALSE]
  offset <- rep(particle_offsets, length.out = n)
  particles[, paste0(parameter, "_lRA")] <-
    particles[, paste0(parameter, "_lRA")] + offset
  particles[, paste0(parameter, "_lRB")] <-
    particles[, paste0(parameter, "_lRB")] - offset
  particles
}

timings <- list()
max_difference <- 0
# EMC2's AccumulatR branch has no native cens/trunc race likelihood.
for (n_trials in sizes) {
  for (model_name in names(models)) {
    definition <- models[[model_name]]
    benchmark <- make_case(definition, n_trials)
    n <- particle_counts[[as.character(n_trials)]]
    particles <- lapply(
      benchmark$theta,
      make_particles,
      n = n,
      parameter = definition$distinct
    )

    check_rows <- seq_len(min(16L, n))
    ll_emc <- drop(benchmark$EMC2(
      particles$EMC2[check_rows, , drop = FALSE]
    ))
    ll_acc <- drop(benchmark$AccumulatR(
      particles$AccumulatR[check_rows, , drop = FALSE]
    ))
    difference <- max(abs(ll_emc - ll_acc)) / n_trials
    max_difference <- max(max_difference, difference)
    stopifnot(difference < 1e-6)

    for (package in c("EMC2", "AccumulatR")) {
      for (warmup in 1:2) benchmark[[package]](particles[[package]])
    }
    for (sample in seq_len(samples)) {
      order <- if (sample %% 2L) {
        c("EMC2", "AccumulatR")
      } else {
        c("AccumulatR", "EMC2")
      }
      for (package in order) {
        elapsed <- system.time(
          benchmark[[package]](particles[[package]])
        )[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          scenario = "ordinary",
          model = model_name,
          n_trials = n_trials,
          package = package,
          sample = sample,
          n_particles = n,
          us_per_particle = elapsed * 1e6 / n
        )
      }
    }
  }
}

timings <- do.call(rbind, timings)
medians <- aggregate(
  us_per_particle ~ scenario + model + n_trials + package,
  timings,
  median
)
keys <- c("scenario", "model", "n_trials")
emc <- medians[medians$package == "EMC2", c(keys, "us_per_particle")]
acc <- medians[medians$package == "AccumulatR", c(keys, "us_per_particle")]
summary <- merge(
  emc, acc,
  by = keys,
  suffixes = c("_emc2", "_accumulatr")
)
summary$accumulatr_over_emc2 <-
  summary$us_per_particle_accumulatr / summary$us_per_particle_emc2

dir.create("dev/scripts/scratch_outputs", showWarnings = FALSE)
write.csv(
  timings,
  "dev/scripts/scratch_outputs/benchmark_emc2_race_models.csv",
  row.names = FALSE
)
print(summary, row.names = FALSE, digits = 4)
cat("maximum absolute log-likelihood difference per trial:",
    format(max_difference, digits = 8), "\n")
