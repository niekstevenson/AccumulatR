# Define the external parameter names for a model

\`set_parameters()\` lets you rename parameters and share them across
accumulators using one simple mapping. Each list name is the external
parameter name users will supply; each value is one or more internal
model parameter names such as \`go.m\` or \`stop.t0\`.

## Usage

``` r
set_parameters(spec, parameters)
```

## Arguments

- spec:

  A \`race_spec\` object.

- parameters:

  A named list mapping external names to one or more model parameter
  names.

## Value

The updated \`race_spec\`.

## Details

Any parameters not mentioned in \`parameters\` keep their default names.

## Examples

``` r
spec <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("go", "go") |>
  add_outcome("stop", "stop") |>
  set_parameters(list(
    drift = c("go.m", "stop.m"),
    spread = c("go.s", "stop.s"),
    onset = "go.t0"
  ))

sampled_pars(spec)
#> [1] "drift"   "spread"  "go.q"    "onset"   "stop.q"  "stop.t0"
```
