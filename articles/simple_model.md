# A Simple Race Model

This vignette shows the basic simulate-and-fit workflow for a race
model. The example has two accumulators feeding response `R1` and one
accumulator feeding response `R2`.

``` r
library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

**Define the model** We use three accumulators. `R1_A` and `R1_B` both
feed response `R1`, while `R2` feeds response `R2`. The
[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md)
call shares `A`, `B`, and a variability parameter across the LBA and RDM
accumulators, and shares `t0` across all three accumulators.

``` r
model <- race_spec() |>
  add_accumulator("R1_A", "LBA") |>
  add_accumulator("R1_B", "RDM") |>
  add_accumulator("R2", "lognormal") |>
  add_pool("R1", c("R1_A", "R1_B")) |>
  add_outcome("R1", "R1") |>
  add_outcome("R2", "R2") |>
  set_parameters(list(
    B_shared = c("R1_A.B", "R1_B.B"),
    A_shared = c("R1_A.A", "R1_B.A"),
    noise_shared = c("R1_A.sv", "R1_B.s"),
    t0_shared = c("R1_A.t0", "R1_B.t0", "R2.t0")
  )) |>
  finalize_model()

true_params <- c(
  R1_A.v = 2,
  R1_B.v = 3,
  B_shared = 1,
  A_shared = 0.3,
  noise_shared = 1,
  R2.m = log(0.4),
  R2.s = 0.18,
  t0_shared = 0.05
)
```

**Simulate data** We generate response outcomes and response times for
each trial.

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
    ##   R1   R2 
    ## 1572  428

**Evaluate the likelihood** We prepare the data, build a model context,
and evaluate the log-likelihood at the true parameter values.

``` r
prepared <- prepare_data(model, data_df)
ctx <- make_context(model)

params_df_true <- build_param_matrix(
  model,
  true_params,
  trial_df = prepared
)

ll_true <- as.numeric(log_likelihood(ctx, prepared, params_df_true))
ll_true
```

    ## [1] 1281.428

For comparison, we can evaluate a clearly misspecified parameter set.

``` r
wrong_params <- true_params
wrong_params["t0_shared"] <- 0.08

params_df_wrong <- build_param_matrix(
  model,
  wrong_params,
  trial_df = prepared
)

ll_wrong <- as.numeric(log_likelihood(ctx, prepared, params_df_wrong))
ll_wrong
```

    ## [1] 1045.375

**Estimate parameters with
[`optim()`](https://rdrr.io/r/stats/optim.html)** We estimate six
parameters: `R1_A.v`, `R1_B.v`, shared `B`, `R2.m`, `R2.s`, and shared
`t0`. `B`, `R2.s`, and `t0` are estimated on the log scale. The shared
`A` and variability parameters are held fixed for scaling constraints.

``` r
neg_loglik <- function(theta) {
  est <- true_params
  est["R1_A.v"] <- theta[["R1_A.v"]]
  est["R1_B.v"] <- theta[["R1_B.v"]]
  est["B_shared"] <- exp(theta[["log_B_shared"]])
  est["R2.m"] <- theta[["R2.m"]]
  est["R2.s"] <- exp(theta[["log_R2.s"]])
  est["t0_shared"] <- exp(theta[["log_t0_shared"]])
  params_df <- build_param_matrix(
    model,
    est,
    trial_df = prepared
  )
  ll <- log_likelihood(ctx, prepared, params_df)
  -as.numeric(ll)
}

start <- c(
  R1_A.v = 1.5,
  R1_B.v = 1.5,
  log_B_shared = log(1.0),
  R2.m = log(0.32),
  log_R2.s = log(0.15),
  log_t0_shared = log(0.03)
)

set.seed(123456)
fit <- optim(start, neg_loglik, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-9))
```

``` r
fit_params <- c(
  R1_A.v = fit$par[["R1_A.v"]],
  R1_B.v = fit$par[["R1_B.v"]],
  B_shared = exp(fit$par[["log_B_shared"]]),
  R2.m = fit$par[["R2.m"]],
  R2.s = exp(fit$par[["log_R2.s"]]),
  t0_shared = exp(fit$par[["log_t0_shared"]])
)
target <- c(
  R1_A.v = true_params[["R1_A.v"]],
  R1_B.v = true_params[["R1_B.v"]],
  B_shared = true_params[["B_shared"]],
  R2.m = true_params[["R2.m"]],
  R2.s = true_params[["R2.s"]],
  t0_shared = true_params[["t0_shared"]]
)
data.frame(true = target, recovered = fit_params, miss = abs(target - fit_params))
```

    ##                 true   recovered       miss
    ## R1_A.v     2.0000000  1.77020920 0.22979080
    ## R1_B.v     3.0000000  2.80726564 0.19273436
    ## B_shared   1.0000000  0.88548265 0.11451735
    ## R2.m      -0.9162907 -0.94747448 0.03118375
    ## R2.s       0.1800000  0.19522887 0.01522887
    ## t0_shared  0.0500000  0.06410211 0.01410211

This is the core usage pattern of the package: define a model, simulate
or prepare data, evaluate the likelihood, and fit the parameters of
interest.
