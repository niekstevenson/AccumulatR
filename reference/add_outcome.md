# Define an observed response

Attach a response label to a finishing rule. In a single-response model,
the first available outcome determines the observed response and time.

## Usage

``` r
add_outcome(spec, label, expr, options = list())
```

## Arguments

- spec:

  A `race_spec` object.

- label:

  Response label that should appear in the behavioral data.

- expr:

  Accumulator or pool label, or a rule built with
  [`first_of()`](https://niekstevenson.github.io/AccumulatR/reference/first_of.md),
  [`all_of()`](https://niekstevenson.github.io/AccumulatR/reference/all_of.md),
  [`none_of()`](https://niekstevenson.github.io/AccumulatR/reference/none_of.md),
  [`inhibit()`](https://niekstevenson.github.io/AccumulatR/reference/inhibit.md),
  or
  [`build_outcome_expr()`](https://niekstevenson.github.io/AccumulatR/reference/build_outcome_expr.md).

- options:

  Named list of response settings:

  - `component`: component labels in which this outcome is active. Omit
    to use it in every component.

  - `map_outcome_to`: another declared outcome label, or `NA_character_`
    to record no response when this outcome wins.

  - `guess`: a list with declared outcome `labels`, their probability
    `weights` (summing to one), and `rt_policy = "keep"` or `"na"`. When
    this outcome wins, draw its recorded label using these weights;
    `"na"` discards the response time. The default policy is `"keep"`.

  Guessing and remapping require `n_outcomes = 1`.

## Value

The updated `race_spec`.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
```
