# AccumulatR

<img src="man/figures/logo.png" align="right" alt="AccumulatR logo" width="180" />

`AccumulatR` is an R/C++ toolkit for race-model simulation and likelihood evaluation.
Define a model in R, then run simulation and native likelihood kernels with minimal overhead.

## Installation

```r
# from a local checkout
install.packages(".", repos = NULL, type = "source")
```

## Tiny Example

```r
library(AccumulatR)

spec <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_outcome("A", "A") |>
  add_outcome("B", "B")

model <- finalize_model(spec)

pars <- c(
  A.meanlog = log(0.28), A.sdlog = 0.16, A.q = 0, A.t0 = 0,
  B.meanlog = log(0.35), B.sdlog = 0.18, B.q = 0, B.t0 = 0
)

param_df <- build_param_matrix(spec, pars, n_trials = 8)
sim <- simulate(model, param_df, seed = 123)
head(sim[c("trial", "R", "rt")])

ctx <- build_likelihood_context(model, sim[c("trial", "R", "rt")])
log_likelihood(ctx, param_df)
```

`simulate()` returns trial-level choice/RT data; `log_likelihood()` scores those data under a parameter set.

## More Examples

For larger simulation + fitting walkthroughs, see `dev/scripts/vignettes/`.
