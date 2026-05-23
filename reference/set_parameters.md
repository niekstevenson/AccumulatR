# Control parameter grouping and names

Parameters are grouped by compatible type by default. For example, two
lognormal accumulators expose one \`m\`, one \`s\`, and one \`t0\`
parameter unless you ask for specific parameters to be separate.
Triggers expose their trigger name directly, and sampled mixtures expose
automatic \`p.\<component\>\` parameters.

## Usage

``` r
set_parameters(spec, separate = NULL, share = NULL, rename = NULL)
```

## Arguments

- spec:

  A \`race_spec\` object.

- separate:

  Named list. Each name is a grouped public parameter, and each value is
  one or more accumulator ids to split from that group. Use \`TRUE\` to
  split every member of a group.

- share:

  Named list mapping a new public name to default public names or
  internal parameter names that should share one value.

- rename:

  Named character vector mapping current public names to new public
  names.

## Value

The updated \`race_spec\`.

## Examples

``` r
spec <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("go", "go") |>
  add_outcome("stop", "stop") |>
  set_parameters(
    separate = list(m = c("go", "stop")),
    rename = c(s = "spread", t0 = "onset")
  )

par_names(spec)
#> [1] "go.m"   "spread" "onset"  "stop.m"
```
