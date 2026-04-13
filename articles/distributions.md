# Supported Distributions

This vignette lists the accumulator distributions currently available in
`AccumulatR` and the parameter names they use.

``` r
library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

## Naming convention

Parameter names are built as `accumulator.parameter`. For example,
`go.m`, `stop.shape`, and `choice.v` refer to parameters for the
accumulators `go`, `stop`, and `choice`.

## Available distributions

``` r
data.frame(
  distribution = c("lognormal", "gamma", "exgauss", "LBA", "RDM"),
  parameters = c(
    "m, s",
    "shape, rate",
    "mu, sigma, tau",
    "v, B, A, sv",
    "v, B, A, s"
  ),
  stringsAsFactors = FALSE
)
```

    ##   distribution     parameters
    ## 1    lognormal           m, s
    ## 2        gamma    shape, rate
    ## 3      exgauss mu, sigma, tau
    ## 4          LBA    v, B, A, sv
    ## 5          RDM     v, B, A, s

- `lognormal`: `m` and `s`, the usual log-scale location and spread
  parameters.
- `gamma`: `shape` and `rate`.
- `exgauss`: `mu`, `sigma`, and `tau`.
- `LBA`: `v`, `B`, `A`, and `sv`.
- `RDM`: `v`, `B`, `A`, and `s`.

Outside of these all accumulators can have a non-decision time `t0` and
a trigger probability `q`. For the `exgauss` `t0` and `mu` are redundant
(they both constitute a shift in rt) and one must be fixed.

## Example

The same model can combine different accumulator distributions. The only
thing that changes is the parameter suffix used in the parameter vector.

``` r
model <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "exgauss") |>
  add_outcome("go", "go") |>
  add_outcome("stop", "stop") |>
  finalize_model()

params <- c(
  go.m = log(0.30),
  go.s = 0.18,
  stop.mu = 0.10,
  stop.sigma = 0.04,
  stop.tau = 0.08
)

build_param_matrix(model, params, n_trials = 2)
```

    ##      q w t0        p1   p2   p3
    ## [1,] 0 1  0 -1.203973 0.18 0.00
    ## [2,] 0 1  0  0.100000 0.04 0.08
    ## [3,] 0 1  0 -1.203973 0.18 0.00
    ## [4,] 0 1  0  0.100000 0.04 0.08

Choose the distribution that matches the accumulator you want to
specify, then use its parameter names consistently in
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md),
[`simulate()`](https://niekstevenson.github.io/AccumulatR/reference/simulate.md),
and likelihood evaluation.
