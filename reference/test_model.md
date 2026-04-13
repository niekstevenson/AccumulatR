# Test Model Implementation

\`test_model()\` is a quick diagnostic for checking whether a model and
a parameter vector behave sensibly. It does three things:

## Usage

``` r
test_model(
  model,
  param_values,
  n_trials,
  seed = 123,
  include_na = TRUE,
  profile_points = 21L,
  profile_span = 0.5,
  plot = TRUE
)
```

## Arguments

- model:

  Model specification or finalized model.

- param_values:

  Named numeric vector of parameter values. Parameters not supplied are
  assumed to be \`0\`. Names should follow \`sampled_pars(model)\`, so
  custom names from \`set_parameters()\` are supported.

- n_trials:

  Number of trials to simulate.

- seed:

  Random seed used for simulation.

- include_na:

  If \`TRUE\`, include missing-response probability mass in the
  analytical and simulated comparison.

- profile_points:

  Number of grid points per parameter profile.

- profile_span:

  Relative profiling span. Positive and non-negative parameters are
  profiled on a multiplicative scale around their supplied value;
  unrestricted parameters are profiled on a symmetric additive scale.

- plot:

  If \`TRUE\`, draw one profile plot per supplied parameter.

## Value

Invisibly returns a list with the completed parameter vector, the
analytical and simulated probability comparison, the simulated data, the
prepared data, the model context, and the profile tables.

## Details

It compares analytical response probabilities with response proportions
from simulated behavioral data and profiles the log-likelihood for each
parameter supplied, holding the others fixed

This is especially useful when you are building a new model, checking a
custom \`set_parameters()\` mapping, or verifying that a parameter
vector gives sensible behavior before fitting real behavioral data.

## Examples

``` r
spec <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("go", "go") |>
  add_outcome("stop", "stop") |>
  set_parameters(list(
    drift = c("go.m", "stop.m")
  ))

res <- test_model(
  spec,
  c(drift = log(0.30), go.s = 0.16, stop.s = 0.18),
  n_trials = 50,
  plot = FALSE
)
#> Assuming 0 for missing parameter(s): go.q, go.t0, stop.q, stop.t0
#>  response analytical simulated abs_diff
#>        go        0.5       0.4      0.1
#>      stop        0.5       0.6      0.1

res$comparison
#>   response analytical simulated abs_diff
#> 1       go        0.5       0.4      0.1
#> 2     stop        0.5       0.6      0.1
```
