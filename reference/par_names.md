# List the free parameters implied by a model

List the free parameters implied by a model

## Usage

``` r
par_names(model)
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
par_names(spec)
#> [1] "m"  "s"  "t0"
```
