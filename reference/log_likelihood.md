# Evaluate log-likelihoods of behavioral data

Compute the summed log-likelihood by default, or trial-wise
log-likelihoods when \`sum = FALSE\`.

## Usage

``` r
log_likelihood(
  context,
  data,
  parameters,
  ok = NULL,
  sum = TRUE,
  min_ll = log(1e-10),
  ...
)

# S3 method for class 'accumulatr_context'
log_likelihood(
  context,
  data,
  parameters,
  ok = NULL,
  sum = TRUE,
  min_ll = log(1e-10),
  ...
)

# Default S3 method
log_likelihood(context, ...)
```

## Arguments

- context:

  Context created with \`make_context()\`.

- data:

  Prepared data created with \`prepare_data()\`.

- parameters:

  A parameter data frame, or a list of parameter data frames.

- ok:

  Logical vector marking which trials should contribute to the
  likelihood. Trials marked \`FALSE\` are assigned \`min_ll\`.

- sum:

  If \`TRUE\`, return the summed log-likelihood. If \`FALSE\`, return
  trial-wise log-likelihood values.

- min_ll:

  Minimum log-likelihood value used for excluded or impossible trials.

- ...:

  Unused; for S3 compatibility.

## Value

A summed log-likelihood by default, or a numeric vector of trial-wise
log-likelihood values when \`sum = FALSE\`.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  spec,
  c(m = 0, s = 0.1),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
prepared <- prepare_data(structure, data_df)
ctx <- make_context(structure)
log_likelihood(ctx, prepared, params_df)
#> [1] 2.481128
```
