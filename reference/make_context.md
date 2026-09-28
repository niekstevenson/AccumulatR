# Build a compiled likelihood context from a model

Compile the model's response rules and dependencies for likelihood
evaluation. Reuse the context across candidate parameter values and
datasets for the same model. Prepare each dataset with
[`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md).

## Usage

``` r
make_context(structure, diagnostics = FALSE)
```

## Arguments

- structure:

  Finalized model structure.

- diagnostics:

  If `TRUE`, collect model compilation statistics for
  [`complexity_metrics()`](https://niekstevenson.github.io/AccumulatR/reference/complexity_metrics.md).

## Value

An `accumulatr_context` object.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
make_context(structure)
#> $cpp
#> <pointer: 0x5630d0b06d60>
#> 
#> $outcome_labels
#> [1] "A_win"
#> 
#> $observed_outcome_labels
#> [1] "A_win"
#> 
#> attr(,"class")
#> [1] "accumulatr_context"
```
