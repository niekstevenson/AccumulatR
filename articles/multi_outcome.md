# Multiple Outcomes

This vignette shows how to work with data in which more than one ordered
response is observed on each trial. Here each trial records the first
response (`R`, `rt`) and the second response (`R2`, `rt2`), and so on
for any number of responses (`R3`, `rt3` etc.).

``` r
library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

**Define the model** We use a simple two-accumulator race with direct
responses `A` and `B`. Setting `n_outcomes = 2` tells the model to
retain the first and second finishing responses.

``` r
model <- race_spec(n_outcomes = 2L) |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_outcome("A", "A") |>
  add_outcome("B", "B") |>
  finalize_model()

true_params <- c(
  A.m = log(0.30),
  A.s = 0.18,
  B.m = log(0.38),
  B.s = 0.22
)
```

**Simulate data** We generate data and keep both ordered responses for
each trial.

``` r
set.seed(123456)

n_trials <- 500
params_df <- build_param_matrix(model, true_params, n_trials = n_trials)

sim <- simulate(model, params_df)

data_df <- data.frame(
  trial = sim$trial,
  R = factor(sim$R),
  rt = sim$rt,
  R2 = factor(sim$R2),
  rt2 = sim$rt2,
  stringsAsFactors = FALSE
)

head(data_df)
```

    ##   trial R        rt R2       rt2
    ## 1     1 A 0.3394130  B 0.3514512
    ## 2     2 A 0.2373401  B 0.4565748
    ## 3     3 A 0.3709420  B 0.4913625
    ## 4     4 B 0.2974072  A 0.3441526
    ## 5     5 A 0.3297834  B 0.4790856
    ## 6     6 A 0.3671582  B 0.4663272

**Estimate parameters with
[`optim()`](https://rdrr.io/r/stats/optim.html)** We estimate `A.m`,
`A.s`, `B.m`, and `B.s`. The variance parameters are optimized on the
log scale and transformed back inside the function.

``` r
prepared <- prepare_data(model, data_df)
ctx <- make_context(model)
neg_loglik <- function(theta) {
  est <- true_params
  est[c("A.m", "A.s", "B.m", "B.s")] <- theta[c("A.m", "A.s", "B.m", "B.s")]
  est[c("A.s", "B.s")] <- exp(est[c("A.s", "B.s")])
  params_df <- build_param_matrix(
    model,
    est,
    trial_df = prepared
  )
  ll <- log_likelihood(ctx, prepared, params_df)
  -as.numeric(ll)
}

start <- c(
  A.m = log(0.24),
  A.s = log(0.12),
  B.m = log(0.48),
  B.s = log(0.12)
)

set.seed(123456)
fit <- optim(start, neg_loglik, method = "Nelder-Mead")
```

``` r
fit_params <- fit$par
fit_params[c("A.s", "B.s")] <- exp(fit_params[c("A.s", "B.s")])
target <- true_params[c("A.m", "A.s", "B.m", "B.s")]

data.frame(
  true = target,
  recovered = fit_params,
  miss = abs(target - fit_params)
)
```

    ##          true  recovered        miss
    ## A.m -1.203973 -1.2097763 0.005803478
    ## A.s  0.180000  0.1770936 0.002906442
    ## B.m -0.967584 -0.9585490 0.009035000
    ## B.s  0.220000  0.2148618 0.005138242

The workflow is the same as in the single-response case, but the
likelihood now uses the second ranked response as additional
information.
