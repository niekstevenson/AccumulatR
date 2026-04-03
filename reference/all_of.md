# Define a response that requires several processes to finish

Define a response that requires several processes to finish

## Usage

``` r
all_of(...)
```

## Arguments

- ...:

  Accumulator labels or expression objects to combine with AND.

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
