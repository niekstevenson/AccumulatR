# Add a shared absence trigger

A trigger is a named absence-probability parameter. All members in one
trigger call share the same absence draw. Use separate trigger calls for
independent absence draws.

## Usage

``` r
add_trigger(spec, name, members)
```

## Arguments

- spec:

  A \`race_spec\` object.

- name:

  Trigger parameter name.

- members:

  Accumulator labels controlled by the shared absence draw.

## Value

The updated \`race_spec\`.
