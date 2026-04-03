# Evaluate the log-likelihood of behavioral data

Compute the log-likelihood of the behavioral data stored in a
\`likelihood_context\` under one or more candidate parameter sets.

## Usage

``` r
log_likelihood(
  likelihood_context,
  parameters,
  ok = NULL,
  expand = NULL,
  min_ll = log(1e-10),
  ...
)

# S3 method for class 'likelihood_context'
log_likelihood(
  likelihood_context,
  parameters,
  ok = NULL,
  expand = NULL,
  min_ll = log(1e-10),
  ...
)

# Default S3 method
log_likelihood(likelihood_context, ...)
```

## Arguments

- likelihood_context:

  Context created with \`build_likelihood_context()\`.

- parameters:

  A parameter data frame, or a list of parameter data frames.

- ok:

  Logical vector marking which trials should contribute to the
  likelihood. Trials marked \`FALSE\` are assigned \`min_ll\`.

- expand:

  Optional index vector used to expand compressed trial-level results
  back to the original trial count.

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
  c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
ctx <- build_likelihood_context(structure, data_df)
log_likelihood(ctx, list(params_df))
#> [1] 2.481128
```
