# Evaluate log-likelihoods of behavioral data

Compute the summed log-likelihood by default, or trial-wise
log-likelihoods when `sum = FALSE`. Build the context and prepare the
observations once, then reuse them when evaluating candidate parameter
matrices for the same model.

## Usage

``` r
log_likelihood(
  context,
  data,
  parameters,
  ok = NULL,
  sum = TRUE,
  min_ll = log(1e-10)
)
```

## Arguments

- context:

  Context created with
  [`make_context()`](https://niekstevenson.github.io/AccumulatR/reference/make_context.md).

- data:

  Prepared data created with
  [`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md).

- parameters:

  Numeric matrix from
  [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md),
  with one accumulator block per prepared trial in matching order.

- ok:

  Optional logical vector with one value per prepared trial. `TRUE`
  evaluates that trial; `FALSE` assigns `min_ll`. These assigned values
  are included in the sum. For compressed data, use the retained trial
  order.

- sum:

  If `TRUE`, return the summed log-likelihood. If `FALSE`, return
  trial-wise log-likelihood values.

- min_ll:

  Minimum log-likelihood value used for excluded or impossible trials.

## Value

A summed log-likelihood by default, or a numeric vector of trial-wise
log-likelihood values when `sum = FALSE`.

## Details

Response-time observations contribute densities, so a log-likelihood can
be positive. Missing responses and censoring contribute probability
masses according to the observation rules.

Use
[`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md)
and
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md)
for the same model. Evaluation assumes matching layouts and valid
parameters and does not repeat preparation checks.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  structure,
  c(m = 0, s = 0.1),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
prepared <- prepare_data(structure, data_df)
ctx <- make_context(structure)
log_likelihood(ctx, prepared, params_df)
#> [1] 2.481128
```
