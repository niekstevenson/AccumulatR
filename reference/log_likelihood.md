# Evaluate the log-likelihood of behavioral data

Compute the log-likelihood of prepared behavioral data under one or more
candidate parameter sets.

## Usage

``` r
log_likelihood(
  context,
  data,
  parameters,
  ok = NULL,
  expand = NULL,
  min_ll = log(1e-10),
  ...
)

# S3 method for class 'accumulatr_context'
log_likelihood(
  context,
  data,
  parameters,
  ok = NULL,
  expand = NULL,
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

- expand:

  Optional 1-based index vector mapping original trials to evaluated
  trial log-likelihood positions. If \`NULL\`, the \`expand\` attribute
  from \`data\` is used when present; otherwise no expansion is applied.

- min_ll:

  Minimum log-likelihood value used for excluded or impossible trials.

- ...:

  Unused; for S3 compatibility.

## Value

A numeric vector of log-likelihood values.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  spec,
  c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
prepared <- prepare_data(structure, data_df)
ctx <- make_context(structure)
log_likelihood(ctx, prepared, params_df)
#> [1] 2.481128
```
