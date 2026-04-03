# Expand behavioral data to one row per accumulator

This is mainly a helper for advanced workflows where you want
trial-level behavioral data repeated across the accumulators that are
active on each trial.

## Usage

``` r
nest_accumulators(model_spec, data)
```

## Arguments

- model_spec:

  Finalized model returned by \`finalize_model()\`.

- data:

  Behavioral data with one row per trial.

## Value

A data frame with rows repeated per accumulator and an \`accumulator\`
column.
