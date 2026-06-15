# Evaluate marginal response probabilities

\`response_probabilities()\` evaluates the model-implied marginal
probability of each observed response label for a finalized model and
parameter set. Component labels in row-form parameters condition the
calculation on those observed components; latent or \`NA\` components
are marginalized according to the model's mixture specification.

## Usage

``` r
response_probabilities(structure, params_df, include_na = TRUE, ...)

# S3 method for class 'model_structure'
response_probabilities(structure, params_df, include_na = TRUE, ...)

# Default S3 method
response_probabilities(structure, params_df, include_na = TRUE, ...)
```

## Arguments

- structure:

  Finalized model structure.

- params_df:

  A parameter data frame or rectangular parameter matrix.

- include_na:

  If \`TRUE\`, include residual mass as \`"NA"\`.

- ...:

  Unused; for S3 compatibility.

## Value

A named numeric vector of marginal response probabilities. Names are
observed outcome labels. When \`include_na = TRUE\`, a residual \`"NA"\`
entry is included if the model assigns probability mass to unobserved or
\`NA\`-mapped outcomes.

## Examples

``` r
spec <- race_spec() |>
  add_accumulator("left", "lognormal") |>
  add_accumulator("right", "lognormal") |>
  add_outcome("left", "left") |>
  add_outcome("right", "right") |>
  set_parameters(separate = list(m = TRUE))

model <- finalize_model(spec)
params <- build_param_matrix(
  spec,
  c(left.m = log(0.25), right.m = log(0.40), s = 0.20),
  n_trials = 1
)

response_probabilities(model, params)
#>       left      right 
#> 0.95171491 0.04828509 
```
