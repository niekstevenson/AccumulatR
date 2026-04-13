# Store model-level metadata

\`set_metadata()\` stores additional metadata on a model specification.
At present, these values are carried through model finalization but do
not directly change simulation or likelihood evaluation. This is mainly
useful for annotating a model or passing labels such as
\`special_outcomes\` to downstream tools.

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
  add_accumulator("watch", "lognormal") |>
  add_outcome("go", "go") |>
  add_outcome("NR_CENSOR", "watch", options = list(class = "censor"))

spec <- set_metadata(spec, special_outcomes = list(censor = "NR_CENSOR"))
finalize_model(spec)$special_outcomes$censor
#> [1] "NR_CENSOR"
```
