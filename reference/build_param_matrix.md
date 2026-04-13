# Create trial-level parameter values

This expands a named parameter vector into the trial-by-trial format
expected by \`simulate()\` and \`log_likelihood()\`.

## Usage

``` r
build_param_matrix(
  model,
  param_values,
  n_trials = 1L,
  component = NULL,
  trial_df = NULL,
  layout = NULL
)
```

## Arguments

- model:

  Model definition.

- param_values:

  Named numeric vector of parameter values.

- n_trials:

  Number of trials to generate.

- component:

  Optional component label or labels.

- trial_df:

  Optional trial/prepared data object. If it includes an \`accumulator\`
  column, parameter rows are built in that exact row order.

- layout:

  Optional storage layout.

## Value

A data frame or matrix of parameter values by trial.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
vals <- c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0)
build_param_matrix(spec, vals, n_trials = 2)
#>      q w t0 p1  p2 p3
#> [1,] 0 1  0  0 0.1  0
#> [2,] 0 1  0  0 0.1  0
```
