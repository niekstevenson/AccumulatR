# Finalize a model for simulation and fitting

Check source references, timing dependencies, response rules, and
parameter grouping, and create the model object used by
[`simulate()`](https://niekstevenson.github.io/AccumulatR/reference/simulate.md),
[`prepare_data()`](https://niekstevenson.github.io/AccumulatR/reference/prepare_data.md),
and
[`build_param_matrix()`](https://niekstevenson.github.io/AccumulatR/reference/build_param_matrix.md).
Build its likelihood context with
[`make_context()`](https://niekstevenson.github.io/AccumulatR/reference/make_context.md).

## Usage

``` r
finalize_model(model)
```

## Arguments

- model:

  Model specification.

## Value

A `model_structure` object.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
finalize_model(spec)
#> $model_spec
#> $accumulators
#> $accumulators[[1]]
#> $accumulators[[1]]$id
#> [1] "A"
#> 
#> $accumulators[[1]]$dist
#> [1] "lognormal"
#> 
#> $accumulators[[1]]$onset
#> [1] 0
#> 
#> 
#> 
#> $pools
#> list()
#> 
#> $outcomes
#> $outcomes[[1]]
#> $outcomes[[1]]$label
#> [1] "A_win"
#> 
#> $outcomes[[1]]$expr
#> $outcomes[[1]]$expr$kind
#> [1] "event"
#> 
#> $outcomes[[1]]$expr$source
#> [1] "A"
#> 
#> 
#> $outcomes[[1]]$options
#> list()
#> 
#> 
#> 
#> $triggers
#> list()
#> 
#> $parameters
#> $parameters$separate
#> list()
#> 
#> $parameters$share
#> list()
#> 
#> $parameters$rename
#> character(0)
#> 
#> 
#> $components
#> list()
#> 
#> $mixture_options
#> list()
#> 
#> $observation
#> $observation$n_outcomes
#> [1] 1
#> 
#> 
#> attr(,"class")
#> [1] "race_spec"
#> 
#> $prep
#> $prep$accumulators
#> $prep$accumulators$A
#> $prep$accumulators$A$id
#> [1] "A"
#> 
#> $prep$accumulators$A$dist
#> [1] "lognormal"
#> 
#> $prep$accumulators$A$onset
#> [1] 0
#> 
#> $prep$accumulators$A$onset_spec
#> $prep$accumulators$A$onset_spec$kind
#> [1] "absolute"
#> 
#> $prep$accumulators$A$onset_spec$value
#> [1] 0
#> 
#> 
#> $prep$accumulators$A$components
#> character(0)
#> 
#> $prep$accumulators$A$shared_trigger_id
#> NULL
#> 
#> 
#> 
#> $prep$pools
#> named list()
#> 
#> $prep$outcomes
#> $prep$outcomes$A_win
#> $prep$outcomes$A_win$label
#> [1] "A_win"
#> 
#> $prep$outcomes$A_win$expr
#> $prep$outcomes$A_win$expr$kind
#> [1] "event"
#> 
#> $prep$outcomes$A_win$expr$source
#> [1] "A"
#> 
#> 
#> $prep$outcomes$A_win$options
#> list()
#> 
#> 
#> 
#> $prep$components
#> $prep$components$ids
#> [1] "__default__"
#> 
#> $prep$components$weights
#> [1] 1
#> 
#> $prep$components$attrs
#> $prep$components$attrs$`__default__`
#> list()
#> 
#> 
#> $prep$components$mode
#> [1] "fixed"
#> 
#> $prep$components$reference
#> [1] "__default__"
#> 
#> 
#> $prep$observation
#> $prep$observation$n_outcomes
#> [1] 1
#> 
#> $prep$observation$global_n_outcomes
#> [1] 1
#> 
#> $prep$observation$component_n_outcomes
#> named list()
#> 
#> 
#> $prep$shared_triggers
#> named list()
#> 
#> $prep$parameter_lookup
#>  A.m  A.s A.t0 
#>  "m"  "s" "t0" 
#> 
#> $prep$outcomes_by_component
#> $prep$outcomes_by_component$`__default__`
#> [1] "A_win"
#> 
#> 
#> $prep$observed_outcomes_by_component
#> $prep$observed_outcomes_by_component$`__default__`
#> [1] "A_win"
#> 
#> 
#> 
#> $simulation
#> <pointer: (nil)>
#> 
#> attr(,"class")
#> [1] "model_structure" "list"           
```
