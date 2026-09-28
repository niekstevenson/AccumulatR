# Evaluate marginal response probabilities

Calculate the probability of each first response, integrated over
response time and averaged over mixture components. With multiple trial
blocks in `parameters`, return the mean probabilities across those
blocks.

## Usage

``` r
response_probabilities(context, parameters, include_na = TRUE)
```

## Arguments

- context:

  Context created with
  [`make_context()`](https://niekstevenson.github.io/AccumulatR/reference/make_context.md).

- parameters:

  Parameter matrix from
  [`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md).
  Use `n_trials = 1` for one set of response probabilities.

- include_na:

  If `TRUE`, include a `"NA"` entry when there is residual probability
  of no observed response.

## Value

A named numeric vector of marginal response probabilities. Names are
observed outcome labels. When `include_na = TRUE`, a residual `"NA"`
entry is included if the model assigns probability mass to unobserved or
`NA`-mapped outcomes.

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
  model,
  c(left.m = log(0.25), right.m = log(0.40), s = 0.20),
  n_trials = 1
)

response_probabilities(make_context(model), params)
#>       left      right 
#> 0.95171491 0.04828509 
```
