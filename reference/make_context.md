# Build a compiled likelihood context from a model

A context stores compiled model/runtime state only. Behavioral data are
prepared separately with \`prepare_data()\` and supplied to
\`log_likelihood()\`.

## Usage

``` r
make_context(structure, prep = NULL)
```

## Arguments

- structure:

  Finalized model structure.

- prep:

  Optional preprocessed model bundle.

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
#> <pointer: 0x5653d41e3b90>
#> 
#> $cpp$observed_identity
#> [1] TRUE
#> 
#> $cpp$identity_backend
#> [1] "direct"
#> 
#> $cpp$ranked_supported
#> [1] TRUE
#> 
#> 
#> $required_p_slots
#> [1] 2
#> 
#> attr(,"class")
#> [1] "accumulatr_context"
```
