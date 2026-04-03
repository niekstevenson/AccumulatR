# Turn a response rule into an internal expression

Use this when you want to write a response rule programmatically rather
than through the helper functions such as \`all_of()\` or \`inhibit()\`.

## Usage

``` r
build_outcome_expr(expr)
```

## Arguments

- expr:

  Expression or symbol describing an event or blocking rule.

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
#> $args[[1]]$k
#> NULL
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
#> $args[[2]]$arg$k
#> NULL
#> 
#> 
#> 
#> 
```
