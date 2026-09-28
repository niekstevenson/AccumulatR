# Supported Distributions

Each accumulator defines a distribution of finishing times after its
onset. This guide describes the available families, their parameter
names, and the values accepted by
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md).

``` r

library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

## Naming convention

By default, parameters are grouped by parameter name across compatible
accumulators. Use
[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md)
when an accumulator should have its own parameter, such as `go.m`,
`stop.shape`, or `choice.v`.

## Available distributions

| Family | Parameters, in matrix-slot order | Meaning and constraints |
|----|----|----|
| `lognormal` | `m`, `s` | Mean and standard deviation of log finishing time; `s > 0`. |
| `gamma` | `shape`, `rate` | Gamma shape and rate; both positive. The mean is `shape / rate`. |
| `exgauss` | `mu`, `sigma`, `tau` | Mean and standard deviation of a Gaussian component, and mean of an exponential component; `sigma > 0`, `tau > 0`. Their sum is conditioned to be positive. |
| `LBA` | `v`, `B`, `A`, `sv` | Mean drift, threshold gap, start range, and drift standard deviation; `sv > 0`. Drift is drawn from a normal distribution conditioned to be positive. |
| `RDM` | `v`, `B`, `A`, `s` | Drift, threshold gap, start range, and diffusion noise scale; `v >= 0`, `s > 0`. |

For both `LBA` and `RDM`, `A` is the full width of the uniform
starting-point distribution on `[0, A]`, and `B` is the gap from its
upper end to the threshold. The absolute threshold is `B + A`, giving a
uniform initial distance to threshold on `[B, B + A]`. Both `A` and `B`
must be non-negative, with `B + A > 0`; `A = 0` gives a fixed starting
point.

All distribution parameters must be finite. The builder also requires
finite ratios `mu / sigma` and `sigma / tau` for exgauss, `v / sv` for
LBA, and `v / s` and `(B + A) / s` for RDM.

Every accumulator has a nonnegative nondecision time `t0`, added to its
finishing time after onset. It defaults to zero when omitted. For
exgauss, `mu` belongs to the Gaussian component before conditioning on a
positive sum; `t0` shifts the resulting positive finishing time.

Use
[`add_trigger()`](https://niekstevenson.github.io/AccumulatR/reference/add_trigger.md)
to assign probability to an accumulator being absent on a trial. A
trigger probability lies in `[0, 1]`.

## Example

A model can combine different distributions. Use the parameter names
appropriate to each family;
[`par_names()`](https://niekstevenson.github.io/AccumulatR/reference/par_names.md)
shows the names required by the model’s parameter grouping.

``` r

model <- race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "exgauss") |>
  add_outcome("go", "go") |>
  add_outcome("stop", "stop") |>
  set_parameters(separate = list(m = TRUE, s = TRUE, mu = TRUE, sigma = TRUE, tau = TRUE)) |>
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

    ##      q t0        p1   p2   p3
    ## [1,] 0  0 -1.203973 0.18 0.00
    ## [2,] 0  0  0.100000 0.04 0.08
    ## [3,] 0  0 -1.203973 0.18 0.00
    ## [4,] 0  0  0.100000 0.04 0.08
    ## attr(,"class")
    ## [1] "accumulatr_parameters" "matrix"                "array"
