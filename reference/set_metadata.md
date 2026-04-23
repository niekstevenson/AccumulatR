# Store model-level metadata

\`set_metadata()\` stores additional metadata on a model specification.
At present, these values are carried through model finalization but do
not directly change simulation or likelihood evaluation.

## Usage

``` r
set_metadata(spec, ...)
```

## Arguments

- spec:

  A \`race_spec\` object.

- ...:

  Named metadata entries.

## Value

The updated \`race_spec\`.

## Examples

``` r
spec <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_outcome("go", "go")

spec <- set_metadata(spec, note = "example")
finalize_model(spec)$model_spec$metadata$note
#> [1] "example"
```
