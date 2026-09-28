# List the parameter names used by a model

Return the names accepted by
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md),
after applying any grouping or renaming in
[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md).
The list includes nondecision times, trigger probabilities, and sampled
mixture weights where applicable. These names do not determine which
parameters an optimizer must estimate; parameters can be held fixed by
supplying constant values.

## Usage

``` r
par_names(model)
```

## Arguments

- model:

  A `race_spec` or finalized `model_structure` object.

## Value

A character vector of parameter names.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
par_names(spec)
#> [1] "m"  "s"  "t0"
```
