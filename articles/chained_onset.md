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

## Define the model

`A` is directly observed. `B` is latent. `C` is observed and starts only
after `B` has finished.

``` r

model <- race_spec() |>
  add_accumulator("A", "lognormal") |>
  add_accumulator("B", "lognormal") |>
  add_accumulator("C", "lognormal", onset = after("B")) |>
  add_outcome("A", "A") |>
  add_outcome("C", "C") |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
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

## Simulate data

Each trial contributes an observed response label and response time.

``` r

set.seed(123456)

n_trials <- 2000
params_df <- build_param_matrix(model, true_params, n_trials = n_trials)

sim <- simulate(model, params_df)

data_df <- data.frame(
  trials = sim$trials,
  R = factor(sim$R),
  rt = sim$rt,
  stringsAsFactors = FALSE
)

table(data_df$R)
```

    ## 
    ##    A    C 
    ##  446 1554

## Estimate parameters with `optim()`

We estimate `A.m`, `A.s`, `B.m`, `B.s`, `C.m`, and `C.s`. The spread
parameters are optimized on the log scale.

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
    n_trials = n_trials
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

fit <- optim(start, neg_loglik, method = "Nelder-Mead")

fit_params <- fit$par
fit_params[c("A.s", "B.s", "C.s")] <- exp(fit_params[c("A.s", "B.s", "C.s")])
target <- true_params[c("A.m", "A.s", "B.m", "B.s", "C.m", "C.s")]

data.frame(
  true = target,
  recovered = fit_params,
  miss = abs(target - fit_params)
)
```

Only the combined finishing time of `B` and `C` is observed through
response `C`. Separating their timing parameters can therefore be
difficult. Additional observations of the stages or constraints on their
parameters can help identify the individual contributions.
