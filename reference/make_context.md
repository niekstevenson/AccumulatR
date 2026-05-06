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
#> <pointer: 0x55f381addf00>
#> 
#> $cpp$observed_identity
#> [1] TRUE
#> 
#> $cpp$identity_backend
#> [1] "exact"
#> 
#> $cpp$ranked_supported
#> [1] TRUE
#> 
#> $cpp$batch_coverage
#> $cpp$batch_coverage$summary
#>                         category count
#> 1                  BatchComplete     3
#> 2  BatchGroupedButScalarLeafMath     0
#> 3    BatchGroupedButScalarBounds     0
#> 4 BatchGroupedButScalarExprUpper     0
#> 5         ScalarIdentityShortcut     0
#> 6                    Unsupported     0
#> 
#> $cpp$batch_coverage$programs
#>     scope variant_index program_index root_id      category
#> 1 program             0             0      NA BatchComplete
#> 2 program             0             1      NA BatchComplete
#> 3 program             0             2      NA BatchComplete
#>                                            reason
#> 1 no unsupported or scalar leaf-math marker found
#> 2 no unsupported or scalar leaf-math marker found
#> 3 no unsupported or scalar leaf-math marker found
#> 
#> 
#> 
#> $required_p_slots
#> [1] 2
#> 
#> attr(,"class")
#> [1] "accumulatr_context"
```
