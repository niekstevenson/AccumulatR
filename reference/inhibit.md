# Define a response that is blocked by another process

Define a response that is blocked by another process

## Usage

``` r
inhibit(reference, by)
```

## Arguments

- reference:

  Response rule or accumulator label to be blocked.

- by:

  Blocking process or expression.

## Value

A guarded expression object.

## Examples

``` r
inhibit("A", "B")
#> $kind
#> [1] "guard"
#> 
#> $blocker
#> $blocker$kind
#> [1] "event"
#> 
#> $blocker$source
#> [1] "B"
#> 
#> $blocker$k
#> NULL
#> 
#> 
#> $reference
#> $reference$kind
#> [1] "event"
#> 
#> $reference$source
#> [1] "A"
#> 
#> $reference$k
#> NULL
#> 
#> 
```
