# Define a mixture component

Components are useful when trials can come from qualitatively different
processing modes, such as fast versus slow processing.

## Usage

``` r
add_component(
  spec,
  id,
  members,
  weight = NULL,
  weight_param = NULL,
  n_outcomes = NULL,
  attrs = list()
)
```

## Arguments

- spec:

  A \`race_spec\` object.

- id:

  Component label.

- members:

  Accumulator labels that belong to this component.

- weight:

  Optional fixed mixture weight.

- weight_param:

  Optional parameter name for a fitted mixture weight.

- n_outcomes:

  Optional component-specific override for the number of observed
  ordered responses.

- attrs:

  Additional component attributes.

## Value

The updated \`race_spec\`.
