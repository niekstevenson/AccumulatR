# Define a mixture component

Components are useful when trials can come from qualitatively different
processing modes, such as fast versus slow processing. Component
declarations define membership only; component probabilities are
configured with \`set_mixture()\`.

## Usage

``` r
add_component(spec, id, members, n_outcomes = NULL)
```

## Arguments

- spec:

  A \`race_spec\` object.

- id:

  Component label.

- members:

  Accumulator labels that belong to this component.

- n_outcomes:

  Optional component-specific override for the number of observed
  ordered responses.

## Value

The updated \`race_spec\`.
