# Add a shared trigger or gate

Triggers let several accumulators share the same stochastic gating
event, which can be useful in stop-signal or contingent-processing
models.

## Usage

``` r
add_trigger(spec, id, members, q = NULL, draw = c("shared", "independent"))
```

## Arguments

- spec:

  A \`race_spec\` object.

- id:

  Trigger label. If \`NULL\`, the first member label is used.

- members:

  Accumulator labels controlled by the trigger.

- q:

  Failure probability for the trigger. The trigger \`id\` is the
  parameter name for this quantity.

- draw:

  Whether trigger failures are shared across members or drawn
  independently.

## Value

The updated \`race_spec\`.
