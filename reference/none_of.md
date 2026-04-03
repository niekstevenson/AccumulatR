# Define the absence of an event

Define the absence of an event

## Usage

``` r
none_of(expr)

exclude(expr)
```

## Arguments

- expr:

  Accumulator label or expression to negate.

## Value

An expression object.

## Examples

``` r
none_of("A")
#> $kind
#> [1] "not"
#> 
#> $arg
#> $arg$kind
#> [1] "event"
#> 
#> $arg$source
#> [1] "A"
#> 
#> $arg$k
#> NULL
#> 
#> 
```
