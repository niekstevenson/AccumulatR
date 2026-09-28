# Create trial-level parameter values

Expand a named parameter vector into the matrix used by
[`simulate()`](https://niekstevenson.github.io/AccumulatR/reference/simulate.md)
and
[`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md).
Each trial receives the same parameter values; rows are grouped by
trial, with accumulators in their model declaration order.

## Usage

``` r
build_param_matrix(model, param_values, n_trials = 1L)
```

## Arguments

- model:

  Finalized model structure.

- param_values:

  Named numeric vector using the names from
  [`par_names()`](https://niekstevenson.github.io/AccumulatR/reference/par_names.md).
  Omitted nondecision-time parameters (`t0`) default to zero. Supply all
  other parameters on their natural scale, after any fitting
  transformations.

- n_trials:

  Positive integer number of trial blocks to create.

## Value

A numeric matrix with one row per trial/accumulator pair. Columns
contain trigger absence probability `q`, nondecision time `t0`,
distribution parameters `p1`, `p2`, and so on, followed by sampled
mixture weights. Distribution slots follow the order listed in the
Supported Distributions vignette. Accumulators without a trigger have
`q = 0`.

## Details

Parameter domains are checked before expansion. All values must be
finite, `t0` must be nonnegative, and distribution scales must be
positive. Trigger probabilities must lie in `[0, 1]`; sampled mixture
weights must be nonnegative and sum to at most one. See
[`vignette("distributions", package = "AccumulatR")`](https://niekstevenson.github.io/AccumulatR/articles/distributions.md)
for distribution-specific constraints.

For trial-specific parameters, modify the relevant matrix rows while
preserving their order and column layout. Likelihood and simulation
calls use those values directly, so modified values must satisfy the
same domains.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
vals <- c(m = 0, s = 0.1)
build_param_matrix(finalize_model(spec), vals, n_trials = 2)
#>      q t0 p1  p2
#> [1,] 0  0  0 0.1
#> [2,] 0  0  0 0.1
#> attr(,"class")
#> [1] "accumulatr_parameters" "matrix"                "array"                
```
