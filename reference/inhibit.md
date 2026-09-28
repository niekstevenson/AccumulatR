# Define a response that is blocked by another process

The response occurs when `reference` finishes, provided `by` has not
finished strictly earlier. A blocker that finishes later does not cancel
a response that has already occurred.

## Usage

``` r
inhibit(reference, by)
```

## Arguments

- reference:

  Response rule, accumulator label, or pool label to be blocked.

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
#> 
#> $reference
#> $reference$kind
#> [1] "event"
#> 
#> $reference$source
#> [1] "A"
#> 
#> 
```
