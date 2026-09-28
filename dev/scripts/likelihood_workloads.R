# Source from the repository root, then pass likelihood_cases to benchmark_cumulative().
library(AccumulatR)
source('dev/examples/benchmark_models.R')

likelihood_case <- function(model, data = NULL, trials = 100L, particles = 20L) {
  parameters <- lapply(seq(.98, 1.02, length.out = particles), function(scale)
    build_param_matrix(model$structure, model$pars * scale, trials))
  if (is.null(data)) data <- simulate(model$structure, parameters[[1L]],
                                     seed = 831L, keep_component = TRUE)
  list(spec = model$structure, data = prepare_data(model$structure, data),
       parameters = parameters, expand = seq_len(trials))
}

likelihood_cases <- lapply(benchmark_models, likelihood_case)

# Equal, mixed and unique parameter rows exercise sharing and its overhead.
vary_rows <- function(case, mixed = FALSE) {
  trials <- length(case$expand)
  offsets <- seq(-.1, .1, length.out = trials)
  if (mixed) offsets[seq_len(trials %/% 2)] <- 0
  rows <- rep(offsets, each = length(case$spec$prep$accumulators))
  case$parameters <- lapply(case$parameters, function(p) {
    p[, 'p1'] <- p[, 'p1'] + rows
    p
  })
  case
}
for (name in c('example_6_dual_path', 'example_16_guard_tie_simple', 'stim_selective_stop')) {
  likelihood_cases[[paste0(name, '_mixed')]] <- vary_rows(likelihood_cases[[name]], TRUE)
  likelihood_cases[[paste0(name, '_unique')]] <- vary_rows(likelihood_cases[[name]])
}

for (name in c('example_1_simple', 'example_6_dual_path')) {
  model <- benchmark_models[[name]]
  labels <- names(model$structure$prep$outcomes)
  data <- data.frame(trials = 1:100, R = rep(labels, length.out = 100),
                     rt = seq(.22, .70, length.out = 100), LT = .15, UT = .90)
  likelihood_cases[[paste0(name, '_truncated')]] <- likelihood_case(model, data)
  data$rt <- NA_real_
  data$LC <- .28
  data$UC <- .48
  data$missingness <- rep(1:3, length.out = 100)
  likelihood_cases[[paste0(name, '_known_censored')]] <- likelihood_case(model, data)
  data$R <- NA_character_
  likelihood_cases[[paste0(name, '_unknown_censored')]] <- likelihood_case(model, data)
  for (type in c('truncated', 'known_censored', 'unknown_censored')) {
    label <- paste0(name, '_', type)
    likelihood_cases[[paste0(label, '_unique')]] <- vary_rows(likelihood_cases[[label]])
  }

  # Missing RT is an explicit observation mapping, not an invalid ordinary trial.
  spec <- model$structure$model_spec
  for (i in seq_along(spec$outcomes)) spec$outcomes[[i]]$options <- list(
    guess = list(labels = labels[i], weights = 1, rt_policy = 'na'))
  missing <- list(structure = finalize_model(spec), pars = model$pars)
  data <- data.frame(trials = 1:100, R = rep(labels, length.out = 100), rt = NA_real_)
  label <- paste0(name, '_missing_rt')
  likelihood_cases[[label]] <- likelihood_case(missing, data)
  likelihood_cases[[paste0(label, '_unique')]] <- vary_rows(likelihood_cases[[label]])
}

# Small/large k, including both boundaries; the evaluator sees ordinary pool models.
for (n in c(8L, 24L)) for (k in c(1L, 2L, n - 1L, n)) {
  members <- paste0('a', seq_len(n))
  spec <- race_spec()
  for (id in c(members, 'b')) spec <- add_accumulator(spec, id, 'lognormal')
  spec <- spec |> add_pool('pool', members, k = k) |>
    add_outcome('A', 'pool') |> add_outcome('B', 'b') |>
    set_parameters(separate = list(m = TRUE)) |> finalize_model()
  model <- list(structure = spec, pars = c(
    setNames(log(seq(.3, .5, length.out = n)), paste0(members, '.m')),
    b.m = log(.4), s = .2))
  likelihood_cases[[paste0('pool_', k, '_of_', n)]] <- likelihood_case(model)
}

# Same amount of numerical work, but interleaved contexts expose setup eviction.
likelihood_cases$alternating_contexts <- list(interleave = likelihood_cases[c(
  'example_1_simple', 'example_16_guard_tie_simple', 'example_23_ranked_chain')])
