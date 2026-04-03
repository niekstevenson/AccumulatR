# Getting Started with AccumulatR

This article is a short orientation to the main objects and column names
used throughout `AccumulatR`. It is meant to answer the first questions
most users have before they start fitting specific models.

``` r
library(AccumulatR)
```

    ## 
    ## Attaching package: 'AccumulatR'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

## The basic workflow

Most analyses in `AccumulatR` follow the same four steps:

1.  Define a model with accumulators, pools, and observed responses.
2.  Finalize that model with
    [`finalize_model()`](../reference/finalize_model.md).
3.  Work with behavioral data, either simulated with
    [`simulate()`](../reference/simulate.md) or collected in an
    experiment.
4.  Build a likelihood context and evaluate candidate parameter values.

## Core ideas

### Trial

A `trial` is a single behavioral observation. In the simplest case, each
row of your data corresponds to one trial.

### `R`

`R` is the observed response label on a trial. For a two-choice task,
typical values might be `"left"` and `"right"`, or `"go"` and `"stop"`.

### `rt`

`rt` is the observed response time for `R`.

### `R2` and `rt2`

If you use `race_spec(n_outcomes = 2)`, the model keeps the first two
ordered responses on each trial. In that case:

- `R` and `rt` are the first observed response and its response time.
- `R2` and `rt2` are the second observed response and its response time.

The same logic extends to `R3`, `rt3`, and beyond when more ordered
outcomes are retained.

### Component

A `component` is a latent submodel inside a mixture model. This is
useful when you think behavior comes from qualitatively different trial
types, such as:

- a fast-guess process on some trials
- a more controlled decision process on others

If your model has only one component, you usually do not need to think
about this column at all.

### Onset

An `onset` is the time at which an accumulator becomes active. Some
accumulators start at time 0, while others begin later. In `AccumulatR`,
an onset can be:

- a fixed delay, such as `onset = 0.15`
- a chained onset, such as `onset = after("cue")`, meaning the process
  starts only after another accumulator or pool has finished

## A minimal model

The example below defines a basic two-choice race model.

``` r
spec <- race_spec() |>
  add_accumulator("left", "lognormal") |>
  add_accumulator("right", "lognormal") |>
  add_outcome("left", "left") |>
  add_outcome("right", "right")

model <- finalize_model(spec)
```

## Simulated behavioral data

We can generate a small behavioral dataset from this model.

``` r
pars <- c(
  left.meanlog = log(0.28), left.sdlog = 0.16, left.q = 0, left.t0 = 0,
  right.meanlog = log(0.35), right.sdlog = 0.18, right.q = 0, right.t0 = 0
)

param_df <- build_param_matrix(spec, pars, n_trials = 6)
sim <- simulate(model, param_df, seed = 123)

sim[c("trial", "R", "rt")]
```

    ##   trial     R        rt
    ## 1     1  left 0.3182631
    ## 2     2  left 0.3414185
    ## 3     3  left 0.2883234
    ## 4     4  left 0.3669476
    ## 5     5  left 0.3057125
    ## 6     6 right 0.2456584

Here each row is one trial, `R` is the observed response, and `rt` is
the response time.

## From behavioral data to likelihood

To fit a model, the usual pattern is:

``` r
ctx <- build_likelihood_context(model, sim[c("trial", "R", "rt")])
log_likelihood(ctx, param_df)
```

    ## [1] 5.291506

[`build_likelihood_context()`](../reference/build_likelihood_context.md)
prepares the behavioral data for repeated likelihood evaluation.
[`log_likelihood()`](../reference/log_likelihood.md) then tells you how
well a parameter set explains those observed data.

## Where to go next

- If you want the basic fitting workflow, read the simple race-model
  vignette.
- If your experiment records ordered responses, read the multi-outcome
  vignette.
- If your model has staged processing, read the chained-onset vignette.
