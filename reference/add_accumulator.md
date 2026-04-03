# Add an accumulator to a model

Add an accumulator to a model

## Usage

``` r
add_accumulator(spec, id, dist, onset = 0)
```

## Arguments

- spec:

  A \`race_spec\` object.

- id:

  Label for the accumulator.

- dist:

  Distribution family used for that accumulator.

- onset:

  Start time for the accumulator. This can be a fixed numeric onset or a
  chained onset created with \`after()\`.

## Value

The updated \`race_spec\`.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
```
