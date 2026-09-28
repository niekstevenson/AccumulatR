# Build a response rule from a quoted expression

Convert a source label or quoted R expression into a response rule for
[`add_outcome()`](https://niekstevenson.github.io/AccumulatR/reference/add_outcome.md).
Within a quoted expression, `&` combines requirements, `|` allows
alternative routes, and `!` specifies an absence condition. These
correspond to
[`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md),
[`first_of()`](https://niekstevenson.github.io/AccumulatR/reference/first_of.md),
and
[`none_of()`](https://niekstevenson.github.io/AccumulatR/reference/none_of.md).

## Usage

``` r
build_outcome_expr(expr)
```

## Arguments

- expr:

  Accumulator or pool label, symbol, quoted logical expression, or an
  expression object returned by an outcome helper.

## Value

An expression object used inside model specifications.

## Examples

``` r
build_outcome_expr(quote(A & !B))
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
#> [1] "not"
#> 
#> $args[[2]]$arg
#> $args[[2]]$arg$kind
#> [1] "event"
#> 
#> $args[[2]]$arg$source
#> [1] "B"
#> 
#> 
#> 
#> 
```
