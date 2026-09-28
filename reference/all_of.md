# Define a response that requires several processes to finish

Return a rule that finishes when all required events have completed. Any
[`none_of()`](https://niekstevenson.github.io/AccumulatR/reference/none_of.md)
conditions are checked at that finishing time.

## Usage

``` r
all_of(...)
```

## Arguments

- ...:

  Accumulator labels, pool labels, or response expressions.

## Value

An expression object.

## Examples

``` r
all_of("A", "B")
#> $kind
#> [1] "and"
#> 
#> $args
#> $args[[1]]
#> $args[[1]]$kind
#> [1] "event"
#> 
#> $args[[1]]$source
#> [1] "A"
#> 
#> 
#> $args[[2]]
#> $args[[2]]$kind
#> [1] "event"
#> 
#> $args[[2]]$source
#> [1] "B"
#> 
#> 
#> 
```
