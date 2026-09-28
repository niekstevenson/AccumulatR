# Define the absence of an event

Use inside
[`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md)
to require that `expr` has not finished before the conjunction
completes. Absence supplies a condition, not a response time, so it
cannot be an outcome on its own or a standalone
[`first_of()`](https://niekstevenson.github.io/AccumulatR/reference/first_of.md)
branch.

## Usage

``` r
none_of(expr)
```

## Arguments

- expr:

  Accumulator label, pool label, or expression to negate.

## Value

An expression object.

## Examples

``` r
all_of("go", none_of("stop"))
#> $kind
#> [1] "and"
#> 
#> $args
#> $args[[1]]
#> $args[[1]]$kind
#> [1] "event"
#> 
#> $args[[1]]$source
#> [1] "go"
#> 
#> 
#> $args[[2]]
#> $args[[2]]$kind
#> [1] "not"
#> 
#> $args[[2]]$arg
#> $args[[2]]$arg$kind
#> [1] "event"
#> 
#> $args[[2]]$arg$source
#> [1] "stop"
#> 
#> 
#> 
#> 
```
