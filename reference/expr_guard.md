# Build a blocking rule explicitly

Most users will prefer \`inhibit()\`, but \`expr_guard()\` is available
when you want to construct the guarded expression directly.

## Usage

``` r
expr_guard(blocker, reference)
```

## Arguments

- blocker:

  Blocking expression or accumulator label.

- reference:

  Target expression or accumulator label.

## Value

A guarded expression object.

## Examples

``` r
expr_guard("B", "A")
#> $kind
#> [1] "guard"
#> 
#> $blocker
#> [1] "B"
#> 
#> $reference
#> [1] "A"
#> 
```
