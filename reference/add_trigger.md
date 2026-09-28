# Add a shared absence trigger

With probability given by `name`, all member accumulators are absent on
a trial. Otherwise they follow their specified onsets and finishing-time
distributions. Use separate triggers for independent absence events;
[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md)
can give those events a common probability.

## Usage

``` r
add_trigger(spec, name, members)
```

## Arguments

- spec:

  A `race_spec` object.

- name:

  Trigger parameter name. Supply its probability in `[0, 1]` to
  [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md).

- members:

  Accumulator labels controlled by the shared absence draw.

## Value

The updated `race_spec`.
