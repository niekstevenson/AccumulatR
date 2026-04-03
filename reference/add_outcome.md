# Define an observed response

Define an observed response

## Usage

``` r
add_outcome(spec, label, expr, options = list())
```

## Arguments

- spec:

  A \`race_spec\` object.

- label:

  Response label that should appear in the behavioral data.

- expr:

  Rule describing when that response is observed.

- options:

  Optional response settings.

## Value

The updated \`race_spec\`.

## Examples

``` r
spec <- race_spec()
spec <- add_outcome(spec, "A_win", "A")
```
