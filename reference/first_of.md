# Define a response that occurs when the first listed process finishes

Define a response that occurs when the first listed process finishes

## Usage

``` r
first_of(...)
```

## Arguments

- ...:

  Accumulator labels or expression objects to combine with OR.

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
#> $args[[1]]$k
#> NULL
#> 
#> 
#> $args[[2]]
#> $args[[2]]$kind
#> [1] "event"
#> 
#> $args[[2]]$source
#> [1] "B"
#> 
#> $args[[2]]$k
#> NULL
#> 
#> 
#> 
```
