# Simulate behavioral data from a model

Simulate behavioral data from a model

## Usage

``` r
simulate(
  structure,
  params_df,
  trial_df = NULL,
  seed = NULL,
  keep_detail = FALSE,
  keep_component = NULL,
  layout = c("auto", "rectangular", "long"),
  ...
)

# S3 method for class 'model_structure'
simulate(
  structure,
  params_df,
  trial_df = NULL,
  seed = NULL,
  keep_detail = FALSE,
  keep_component = NULL,
  layout = c("auto", "rectangular", "long"),
  ...
)

# Default S3 method
simulate(structure, ...)
```

## Arguments

- structure:

  Finalized model structure.

- params_df:

  Trial-level parameter values.

- trial_df:

  Optional data frame used to condition the simulation. When it includes
  an \`accumulator\` column, \`onset\` and \`component\` values are
  matched by trial and accumulator. Otherwise, values apply at the trial
  level.

- seed:

  Optional random-number seed.

- keep_detail:

  If \`TRUE\`, keep additional simulation detail.

- keep_component:

  Whether to keep the chosen mixture component in the output when the
  model has multiple components. If \`NULL\`, fixed mixtures keep the
  component label and sampled mixtures drop it.

- layout:

  Parameter layout. \`"rectangular"\` expects rows ordered by trial and
  accumulator. \`"long"\` expects rows aligned to \`trial_df\` or
  inferable from \`component\`. \`"auto"\` chooses automatically.

- ...:

  Unused; for S3 compatibility.

## Value

A data frame of simulated behavioral data. For standard models this
includes \`trial\`, \`R\`, and \`rt\`. If \`n_outcomes \> 1\`,
additional ordered response columns such as \`R2\`/\`rt2\` are included.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params <- c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0)
df <- build_param_matrix(spec, params, n_trials = 3)
simulate(structure, df, seed = 123)
#>   trial     R       rt
#> 1     1 A_win 1.083347
#> 2     2 A_win 1.168675
#> 3     3 A_win 1.131959
```
