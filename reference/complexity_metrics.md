# Inspect the size of a compiled likelihood plan

Report how many symbolic regions, numerical operations, and integration
kernels the model requires. Use these counts to investigate models whose
context construction or likelihood evaluation is expensive.

## Usage

``` r
complexity_metrics(context)
```

## Arguments

- context:

  Context created with `make_context(diagnostics = TRUE)`.

## Value

A list containing a `variants` data frame and a `total` list. Each
variant is a compiled component plan. Columns count symbolic regions and
cells, compiled roots and nodes, and integral kernels. Fields beginning
with `max_` report maxima; other total fields sum across variants.
