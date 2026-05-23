# Build a compiled likelihood context from a model

A context stores compiled model/runtime state only. Behavioral data are
prepared separately with \`prepare_data()\` and supplied to
\`log_likelihood()\`.

## Usage

``` r
make_context(structure, prep = NULL, diagnostics = FALSE)
```

## Arguments

- structure:

  Finalized model structure.

- prep:

  Optional preprocessed model bundle.

- diagnostics:

  If \`TRUE\`, collect symbolic/compiled complexity metrics.

## Value

An \`accumulatr_context\` object.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
make_context(structure)
#> $cpp
#> $cpp$native
#> <pointer: 0x555eb1cb4b90>
#> 
#> $cpp$has_complexity_metrics
#> [1] FALSE
#> 
#> 
#> $required_p_slots
#> [1] 2
#> 
#> attr(,"class")
#> [1] "accumulatr_context"
```
