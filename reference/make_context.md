# Build a compiled likelihood context from a model

A context stores compiled model/runtime state only. Behavioral data are
prepared separately with \`prepare_data()\` and supplied to
\`log_likelihood()\`.

## Usage

``` r
make_context(structure, prep = NULL, native_bundle = NULL)
```

## Arguments

- structure:

  Finalized model structure.

- prep:

  Optional preprocessed model bundle.

- native_bundle:

  Optional serialized native bundle.

## Value

An \`accumulatr_context\` object.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
make_context(structure)
#> $native_ctx
#> <pointer: 0x556f96ff8df0>
#> 
#> attr(,"class")
#> [1] "accumulatr_context"
```
