# Control how mixture components are combined

Fixed mixtures use known component probabilities. Sampled mixtures
expose automatic \`p.\<component\>\` parameters for every non-reference
component; the reference component receives the residual probability.

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

  A \`race_spec\` object.

- mode:

  Mixture mode. Fixed mixtures use known component probabilities;
  sampled mixtures estimate probabilities for all non-reference
  components.

- weights:

  Named numeric component probabilities for fixed mixtures. If \`NULL\`,
  fixed mixtures use uniform component probabilities.

- reference:

  Reference component for sampled mixtures. Its probability is the
  residual probability after non-reference component probabilities.

## Value

The updated \`race_spec\`.
