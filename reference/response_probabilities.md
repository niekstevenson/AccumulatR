# Evaluate marginal response probabilities

Evaluate marginal response probabilities

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
