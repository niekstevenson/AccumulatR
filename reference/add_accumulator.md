# Add an accumulator to a model

An accumulator has a random finishing time drawn from `dist`, shifted by
its onset and nondecision time `t0`. Use
[`add_outcome()`](https://niekstevenson.github.io/AccumulatR/reference/add_outcome.md)
to connect its completion to an observed response.

## Usage

``` r
add_accumulator(spec, id, dist, onset = 0)
```

## Arguments

- spec:

  A `race_spec` object.

- id:

  Label for the accumulator.

- dist:

  Distribution family: `"lognormal"`, `"gamma"`, `"exgauss"`, `"LBA"`,
  or `"RDM"`. Names are case-insensitive.

- onset:

  Start time for the accumulator. This can be a fixed numeric onset or a
  chained onset created with
  [`after()`](https://niekstevenson.github.io/AccumulatR/reference/after.md).

## Value

The updated `race_spec`.

## See also

[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md),
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md),
[`after()`](https://niekstevenson.github.io/AccumulatR/reference/after.md)

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
```
