# Start one accumulator after another process finishes

This is useful for staged or contingent architectures, where one process
can only begin after an earlier accumulator or pool has finished.

## Usage

``` r
after(source, lag = 0)
```

## Arguments

- source:

  Accumulator or pool label that must finish first.

- lag:

  Optional non-negative delay added after \`source\` finishes.

## Value

A chained-onset specification for \`add_accumulator(onset = ...)\`.

## Examples

``` r
after("A")
#> $kind
#> [1] "after"
#> 
#> $source
#> [1] "A"
#> 
#> $lag
#> [1] 0
#> 
#> attr(,"class")
#> [1] "race_onset_after" "list"            
after("pool1", lag = 0.05)
#> $kind
#> [1] "after"
#> 
#> $source
#> [1] "pool1"
#> 
#> $lag
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "race_onset_after" "list"            
```
