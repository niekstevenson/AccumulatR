# Start a race-model specification

Start a race-model specification

## Usage

``` r
race_spec(n_outcomes = 1L)
```

## Arguments

- n_outcomes:

  Number of ordered observed responses to retain per trial. Use \`1\`
  for standard choice/RT data, \`2\` when you also observe the second
  finishing response, and so on.

## Value

A \`race_spec\` object.

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
#> list()
#> 
#> $components
#> list()
#> 
#> $mixture_options
#> list()
#> 
#> $metadata
#> $metadata$observation
#> $metadata$observation$mode
#> [1] "top_k"
#> 
#> $metadata$observation$n_outcomes
#> [1] 1
#> 
#> 
#> 
#> attr(,"class")
#> [1] "race_spec"
```
