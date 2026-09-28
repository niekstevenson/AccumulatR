# Getting Started with AccumulatR

`AccumulatR` describes responses as races between latent processes. This
guide introduces the model objects, data columns, and steps needed to
simulate observations and evaluate their likelihood.

``` r

library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

## The basic workflow

Most workflows in `AccumulatR` follow the same four steps:

1.  Define a model with accumulators, pools, and observed responses.
2.  Finalize that model with
    [`finalize_model()`](https://niekstevenson.github.io/AccumulatR/reference/finalize_model.md).
3.  Work with behavioral data, either simulated with
    [`simulate()`](https://niekstevenson.github.io/AccumulatR/reference/simulate.md)
    or collected in an experiment.
4.  Prepare the data, build a model context, and evaluate candidate
    parameter values.

## Core ideas

### `trials`

`trials` identifies the behavioral observation. In the simplest case,
each row of your data corresponds to one trial.

### `R`

`R` is the observed response label on a trial. For a two-choice task,
typical values might be `"left"` and `"right"`, or `"red"` and
`"green"`, etc.

### `rt`

`rt` is the observed response time for `R`. The examples use seconds;
onsets and time parameters use the same units. In a single-response
model, `R = NA` and `rt = NA` record a trial with no observed response.

### Component

A component specifies which accumulators are active on a trial. It can
represent an observed trial type or an unobserved processing mode. A
`component` column in the data conditions on the named component;
omitting it or using `NA` averages over components. Each trial belongs
to one component. See [Working with
Mixtures](https://niekstevenson.github.io/AccumulatR/articles/mixtures.md)
for examples.

### Onset

An `onset` is the time at which an accumulator becomes active. Most
accumulators start at time `0`, while in some cases, extra information
appears during the trial. In `AccumulatR`, an onset can be:

- a fixed delay, such as `onset = 0.15`
- a chained onset, such as `onset = after("A")`, meaning the process
  starts only after another accumulator (in this case `A`) or pool has
  finished

Both onset and non-decision time (`t0`) shift response times. The
difference is where they are supplied: `t0` belongs to the parameter
vector, while onset is declared in the model or supplied with trial
conditions. Omitted `t0` values default to zero.

## A minimal model

The example below defines a basic two-choice race model with two
lognormal accumulators.

``` r

model <- race_spec() |>
  add_accumulator("left", "lognormal") |>
  add_accumulator("right", "lognormal") |>
  add_outcome("left", "left") |>
  add_outcome("right", "right") |>
  set_parameters(separate = list(m = TRUE, s = TRUE)) |>
  finalize_model()
```

## Simulated behavioral data

[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md)
repeats a named parameter vector across trials. It creates one row per
accumulator per trial, with each trial’s accumulators in model order.
Here six trials and two accumulators give twelve parameter rows.
`par_names(model)` lists the accepted parameter names;
[`set_parameters()`](https://niekstevenson.github.io/AccumulatR/reference/set_parameters.md)
determines which values are shared or separate.

``` r

pars <- c(
  left.m = log(0.28), left.s = 0.16,
  right.m = log(0.35), right.s = 0.18
)

param_df <- build_param_matrix(model, pars, n_trials = 6)
sim <- simulate(model, param_df, seed = 123)

head(sim)
```

    ##   trials     R        rt
    ## 1      1  left 0.3182631
    ## 2      2  left 0.3414185
    ## 3      3  left 0.2883234
    ## 4      4  left 0.3669476
    ## 5      5  left 0.3057125
    ## 6      6 right 0.2456584

Each row is one trial, `R` is the observed response, and `rt` is the
response time.

## From behavioral data to likelihood

To fit a model, the usual pattern is:

``` r

prepared <- prepare_data(model, sim[c("trials", "R", "rt")])
ctx <- make_context(model)
log_likelihood(ctx, prepared, param_df)
```

    ## [1] 5.291506

[`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md)
checks and arranges the observations, and
[`make_context()`](https://niekstevenson.github.io/AccumulatR/reference/make_context.md)
compiles the model’s response rules. Reuse these objects while varying
parameters for the same model. Supply one parameter block per prepared
trial, in matching order.

[`log_likelihood()`](https://niekstevenson.github.io/AccumulatR/reference/log_likelihood.md)
returns the summed log-likelihood; use `sum = FALSE` for one value per
trial. Response-time observations contribute densities, so the result
can be positive. See [A Simple Race
Model](https://niekstevenson.github.io/AccumulatR/articles/simple_model.md)
for an example using [`optim()`](https://rdrr.io/r/stats/optim.html) to
estimate parameters.
