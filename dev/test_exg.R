rm(list = ls())

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(AccumulatR)
  }
  library(EMC2)
})

source("tests/testthat/helper-beest.R")

load("~/Documents/2025/EMC2/23B_Spair_exg.RData")

emc_to_beest_par <- function(emc_fit) {
  post <- credint(emc_fit, probs = 0.5)[[1]]

  c(
    S.mu = unname(post[["p1_accumulatorS"]]),
    S.sigma = exp(unname(post[["p2_accumulatorS"]])),
    S.tau = exp(unname(post[["p3_accumulatorS"]])),
    stop.mu = unname(post[["p1_accumulatorstop"]]),
    stop.sigma = exp(unname(post[["p2_accumulatorstop"]])),
    stop.tau = exp(unname(post[["p3_accumulatorstop"]])),
    change.mu = unname(post[["p1_accumulatorchange"]]),
    change.sigma = exp(unname(post[["p2_accumulatorchange"]])),
    change.tau = exp(unname(post[["p3_accumulatorchange"]])),
    stop_trigger = pnorm(unname(post[["q_typeSTOP"]]))
  )
}

prepare_native_data <- function(data) {
  out <- data[, c("trial", "R", "rt", "component", "accumulator", "onset")]
  out$trial <- as.integer(out$trial)
  out$R <- factor(as.character(out$R), levels = c("S", "X"))
  out$rt <- as.numeric(out$rt)
  out$component <- factor(as.character(out$component), levels = c("go_only", "go_stop"))
  out$accumulator <- as.character(out$accumulator)
  out$onset <- as.numeric(out$onset)
  out
}

model <- race_spec() |>
  add_accumulator("S", "exgauss") |>
  add_accumulator("stop", "exgauss") |>
  add_accumulator("change", "exgauss") |>
  add_outcome("S", inhibit("S", by = "stop")) |>
  add_outcome("X", all_of("change", "stop")) |>
  add_component("go_only", members = "S", weight = 0.75) |>
  add_component("go_stop", members = c("S", "stop", "change"), weight = 0.25) |>
  add_trigger(
    "stop_trigger",
    members = c("stop", "change"),
    q = 0.05,
    param = "stop_trigger"
  ) |>
  set_mixture_options(mode = "fixed") |>
  finalize_model()

beest_par <- emc_to_beest_par(emc)
data_long <- prepare_native_data(emc[[1]]$data[[1]])
n_trials <- length(unique(data_long$trial))

ll_custom <- beest_loglik(beest_par, data_long)
ll_custom_conditional <- beest_loglik(
  beest_par,
  data_long,
  include_component_weight = FALSE
)

params_native <- build_param_matrix(model, beest_par, n_trials = n_trials)
ctx_native <- build_likelihood_context(model, data_long)
data_native <- prepare_likelihood_data(ctx_native, data_long)
ll_native <- as.numeric(log_likelihood(ctx_native, data_native, params_native))

print(beest_par)
cat("Standalone log-likelihood (with mixture weights): ", ll_custom, "\n")
cat("Standalone log-likelihood (conditional on component): ", ll_custom_conditional, "\n")
cat("Native log-likelihood (conditional on component):     ", ll_native, "\n")
cat("Difference (conditional):                              ", ll_custom_conditional - ll_native, "\n")

stopifnot(isTRUE(all.equal(ll_custom_conditional, ll_native, tolerance = 1e-3)))
