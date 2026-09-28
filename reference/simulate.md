# Simulate behavioral data from a model

Generate one trial for each accumulator block in `params_df`, using the
model's response rules, timing dependencies, triggers, and mixture
weights.

## Usage

``` r
simulate(
  structure,
  params_df,
  trial_df = NULL,
  seed = NULL,
  keep_detail = FALSE,
  keep_component = NULL
)
```

## Arguments

- structure:

  Finalized model structure.

- params_df:

  Parameter matrix from
  [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md),
  with one row per accumulator per trial, in the order accumulators were
  added to the model.

- trial_df:

  Optional trial-by-accumulator conditioning table containing complete
  `trials`/`racer` blocks in model order. It may supply per-trial
  `component` and per-accumulator `onset` values. Repeat a component
  label across all rows of its trial; `NA` draws a component from the
  mixture. An onset value replaces a fixed onset or adds a delay to a
  chained onset; `NA` uses the onset specified in the model.

- seed:

  Optional seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html). If `NULL`, use the
  current state of R's random-number generator.

- keep_detail:

  If `TRUE`, attach a `details` list containing latent source times and
  outcome candidates for each trial. Inactive or unreachable sources
  have infinite completion times.

- keep_component:

  Whether to keep the chosen mixture component in the output when the
  model has multiple components. If `NULL`, fixed mixtures keep the
  component label and sampled mixtures drop it.

## Value

A data frame with one row per trial and columns `trials`, `R`, and `rt`.
A trial with no observed response has `R = NA` and `rt = NA`. If
`n_outcomes > 1`, additional ordered response pairs such as `R2`/`rt2`
are included; unobserved later ranks are `NA`.

## Details

Use the same model to build `params_df` and simulate the data. Supply
valid parameter values and keep the matrix and conditioning table in
matching row order; simulation does not check their layout or domains.

## See also

[`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md),
[`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md),
[`set_mixture()`](https://niekstevenson.github.io/AccumulatR/reference/set_mixture.md)

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params <- c(m = 0, s = 0.1)
df <- build_param_matrix(structure, params, n_trials = 3)
simulate(structure, df, seed = 123)
#>   trials     R       rt
#> 1      1 A_win 1.083347
#> 2      2 A_win 1.168675
#> 3      3 A_win 1.131959
```
