# Control how mixture components are combined

Fixed mixtures use known component probabilities. Sampled mixtures
expose automatic `p.<component>` parameters for every non-reference
component; the reference component receives the residual probability.
Both modes draw a component during simulation. In likelihood evaluation,
a supplied component label conditions on that component; an absent or
`NA` label averages over components using their probabilities.

## Usage

``` r
set_mixture(
  spec,
  mode = c("fixed", "sample"),
  weights = NULL,
  reference = NULL
)
```

## Arguments

- spec:

  A `race_spec` object.

- mode:

  `"fixed"` stores probabilities in the model. `"sample"` makes
  probabilities parameters supplied to
  [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md).

- weights:

  Named numeric component probabilities for fixed mixtures. If `NULL`,
  fixed mixtures use uniform component probabilities.

- reference:

  Reference component for sampled mixtures. Its probability is one minus
  the sum of non-reference probabilities. Defaults to the last component
  added to the model.

## Value

The updated `race_spec`.
