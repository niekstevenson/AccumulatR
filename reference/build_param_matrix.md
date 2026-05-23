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

  Optional trials/prepared data object. If it includes a \`racer\`
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
vals <- c(m = 0, s = 0.1)
build_param_matrix(spec, vals, n_trials = 2)
#>      q t0 p1  p2 p3 w
#> [1,] 0  0  0 0.1  0 1
#> [2,] 0  0  0 0.1  0 1
```
