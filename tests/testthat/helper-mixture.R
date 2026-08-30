latent_sampled_mixture_spec <- function() {
  race_spec() |>
    add_accumulator("target_fast", "lognormal") |>
    add_accumulator("target_slow", "lognormal") |>
    add_accumulator("competitor", "lognormal") |>
    add_pool("TARGET", c("target_fast", "target_slow")) |>
    add_outcome("Target", "TARGET") |>
    add_outcome("Competitor", "competitor") |>
    add_component("fast", members = c("target_fast", "competitor")) |>
    add_component("slow", members = c("target_slow", "competitor")) |>
    set_mixture(mode = "sample", reference = "slow") |>
    test_separate_all_parameters() |>
    finalize_model()
}

latent_sampled_mixture_params <- function(p.fast) {
  c(
    target_fast.m = log(0.25), target_fast.s = 0.15,
    target_slow.m = log(0.45), target_slow.s = 0.20,
    competitor.m = log(0.35), competitor.s = 0.18,
    p.fast = p.fast
  )
}
