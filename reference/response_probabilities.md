# Compute predicted response probabilities

This function returns the model-implied response probabilities for a
given parameter set. It is useful when you want the predicted response
distribution without evaluating a full trial-by-trial likelihood.

## Usage

``` r
# Default S3 method
response_probabilities(structure, ...)

response_probabilities(structure, params_df, include_na = TRUE, ...)

# S3 method for class 'model_structure'
response_probabilities(structure, params_df, include_na = TRUE, ...)
```

## Arguments

- structure:

  Finalized model structure.

- ...:

  Unused; for S3 compatibility.

- params_df:

  Parameter data frame for one or more trials.

- include_na:

  If \`TRUE\`, include any leftover probability mass assigned to missing
  responses.

## Value

A named numeric vector of response probabilities.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  spec,
  c(A.m = 0, A.s = 0.1, A.q = 0, A.t0 = 0),
  n_trials = 1
)
response_probabilities(structure, params_df)
#> A_win 
#>     1 
```
