# Define a response that occurs when the first listed process finishes

Return a rule whose finishing time is the earliest completion among its
arguments. Each argument must contain an event that can generate a
response.

## Usage

``` r
first_of(...)
```

## Arguments

- ...:

  Accumulator labels, pool labels, or response expressions.

## Value

An expression object.

## Examples

``` r
first_of("A", "B")
#> $kind
#> [1] "or"
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
