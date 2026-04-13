# Prepare behavioral data for likelihood evaluation

\`prepare_data()\` expands trial-level observations to the accumulator
layout expected by the compiled likelihood code and tags the result as
trusted likelihood input.

## Usage

``` r
prepare_data(structure, data_df, prep = NULL, native_bundle = NULL)
```

## Arguments

- structure:

  Finalized model structure.

- data_df:

  Behavioral data. In the simplest case this contains \`trial\`, \`R\`,
  and \`rt\`; for multi-outcome models it can also contain \`R2\`,
  \`rt2\`, and so on.

- prep:

  Optional preprocessed model bundle.

- native_bundle:

  Optional serialized native bundle.

## Value

An \`accumulatr_data\` object.

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
prepare_data(structure, data_df)
#>   trial     R        rt accumulator onset
#> 1     1 A_win 0.9679031           A     0
#> 2     2 A_win 0.9198333           A     0
```
