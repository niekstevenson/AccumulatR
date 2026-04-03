# List the free parameters implied by a model

List the free parameters implied by a model

## Usage

``` r
sampled_pars(model)
```

## Arguments

- model:

  A \`race_spec\` or related model object.

## Value

A character vector of parameter names.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
sampled_pars(spec)
#> [1] "A.meanlog" "A.sdlog"   "A.q"       "A.t0"     
```
