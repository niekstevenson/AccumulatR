# Store model-level settings

Store model-level settings

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
spec <- race_spec()
spec <- set_metadata(spec, rel_tol = 1e-4)
```
