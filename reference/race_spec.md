# Start a race-model specification

Create an empty specification. Add accumulators and response rules, then
call
[`finalize_model()`](https://niekstevenson.github.io/AccumulatR/reference/finalize_model.md)
before simulation or likelihood evaluation.

## Usage

``` r
race_spec(n_outcomes = 1L)
```

## Arguments

- n_outcomes:

  Number of ordered observed responses to retain per trial. Use `1` for
  standard choice/RT data, `2` when you also observe the second
  finishing response, and so on.

## Value

A `race_spec` object.

## Examples

``` r
race_spec()
#> $accumulators
#> list()
#> 
#> $pools
#> list()
#> 
#> $outcomes
#> list()
#> 
#> $triggers
#> list()
#> 
#> $parameters
#> $parameters$separate
#> list()
#> 
#> $parameters$share
#> list()
#> 
#> $parameters$rename
#> character(0)
#> 
#> 
#> $components
#> list()
#> 
#> $mixture_options
#> list()
#> 
#> $observation
#> $observation$n_outcomes
#> [1] 1
#> 
#> 
#> attr(,"class")
#> [1] "race_spec"
```
