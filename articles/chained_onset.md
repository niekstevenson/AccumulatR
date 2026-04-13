# Chained Onsets

This vignette shows how to define staged processing with
[`after()`](https://niekstevenson.github.io/AccumulatR/reference/after.md).
Here accumulator `C` starts only after accumulator `B` has finished.

``` r
library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

**Define the model** `A` is directly observed. `B` is latent. `C` is
observed and starts only after `B` has finished.

``` r
model <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_accumulator("C", "lognormal", onset = after("B")) |>
  add_outcome("A", "A") |>
  add_outcome("C", "C") |>
  finalize_model()

true_params <- c(
  A.m = log(0.28),
  A.s = 0.14,
  B.m = log(0.1),
  B.s = 0.1,
  C.m = log(0.15),
  C.s = 0.1
)
```

**Simulate data** Each trial contributes an observed response label and
response time.

``` r
set.seed(123456)

n_trials <- 2000
params_df <- build_param_matrix(model, true_params, n_trials = n_trials)

sim <- simulate(model, params_df)

data_df <- data.frame(
  trial = sim$trial,
  R = factor(sim$R),
  rt = sim$rt,
  stringsAsFactors = FALSE
)

table(data_df$R)
```

    ## 
    ##    A    C 
    ##  446 1554

**Estimate parameters with
[`optim()`](https://rdrr.io/r/stats/optim.html)** We estimate `A.m`,
`A.s`, `B.m`, `B.s`, `C.m`, and `C.s`. The spread parameters are
optimized on the log scale.

``` r
prepared <- prepare_data(model, data_df)
ctx <- make_context(model)

neg_loglik <- function(theta) {
  est <- true_params
  est[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")] <- theta[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")]
  est[c("A.s", "B.s", "C.s")] <- exp(est[c("A.s", "B.s", "C.s")])
  params_df <- build_param_matrix(
    model,
    est,
    trial_df = prepared
  )
  ll <- log_likelihood(ctx, prepared, params_df)
  -as.numeric(ll)
}

start <- c(
  A.m = log(0.22),
  A.s = log(0.10),
  B.m = log(0.28),
  B.s = log(0.10),
  C.m = log(0.28),
  C.s = log(0.10)
)

set.seed(123456)
fit <- optim(start, neg_loglik, method = "Nelder-Mead")
```

``` r
fit_params <- fit$par
fit_params[c("A.s", "B.s", "C.s")] <- exp(fit_params[c("A.s", "B.s", "C.s")])
target <- true_params[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")]

data.frame(
  true = target,
  recovered = fit_params,
  miss = abs(target - fit_params)
)
```

    ##          true   recovered        miss
    ## A.m -1.272966 -1.25003011 0.022935562
    ## A.s  0.140000  0.09451568 0.045484321
    ## B.m -2.302585 -2.84037595 0.537790857
    ## B.s  0.100000  0.02842790 0.071572096
    ## C.m -1.897120 -1.64098410 0.256135882
    ## C.s  0.100000  0.10959761 0.009597612

Use chained onsets when the model requires a staged dependency between
processes. We do note that these models can suffer from weak
identifiability. So use with care!
