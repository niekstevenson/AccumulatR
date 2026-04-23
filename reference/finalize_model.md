# Compile a model for simulation and fitting

This converts a human-readable model specification into the finalized
object used by \`simulate()\`, \`prepare_data()\`, \`make_context()\`,
and related functions.

## Usage

``` r
finalize_model(model)
```

## Arguments

- model:

  Model specification.

## Value

A \`model_structure\` object.

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
#> $outcomes[[1]]$expr$k
#> NULL
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
#> list()
#> 
#> $components
#> list()
#> 
#> $mixture_options
#> list()
#> 
#> $metadata
#> $metadata$observation
#> $metadata$observation$mode
#> [1] "top_k"
#> 
#> $metadata$observation$n_outcomes
#> [1] 1
#> 
#> 
#> 
#> attr(,"class")
#> [1] "race_model_spec"
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
#> $prep$accumulators$A$q
#> [1] 0
#> 
#> $prep$accumulators$A$params
#> $prep$accumulators$A$params$t0
#> [1] 0
#> 
#> 
#> $prep$accumulators$A$components
#> character(0)
#> 
#> $prep$accumulators$A$shared_trigger_id
#> NULL
#> 
#> $prep$accumulators$A$shared_trigger_q
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
#> $prep$outcomes$A_win$expr$k
#> NULL
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
#> $prep$components$attrs$`__default__`$weight_param
#> NULL
#> 
#> 
#> 
#> $prep$components$has_weight_param
#> [1] FALSE
#> 
#> $prep$components$mode
#> [1] "fixed"
#> 
#> $prep$components$reference
#> [1] "__default__"
#> 
#> 
#> $prep$observation
#> $prep$observation$mode
#> [1] "top_k"
#> 
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
#> list()
#> 
#> $prep$onset_specs
#> $prep$onset_specs$A
#> $prep$onset_specs$A$kind
#> [1] "absolute"
#> 
#> $prep$onset_specs$A$value
#> [1] 0
#> 
#> 
#> 
#> $prep$onset_dependencies
#> $prep$onset_dependencies$A
#> character(0)
#> 
#> 
#> $prep$onset_sources
#> named list()
#> 
#> $prep$onset_topology
#> [1] "A"
#> 
#> $prep$onset_has_dependencies
#> [1] FALSE
#> 
#> 
#> $accumulators
#>        dist onset q shared_trigger_id shared_trigger_q params components
#> A lognormal     0 0              <NA>               NA      0           
#> 
#> $components
#>   component_id weight has_weight_param attrs  mode   reference
#> 1  __default__      1            FALSE  NULL fixed __default__
#> 
#> $shared_triggers
#> list()
#> 
#> attr(,"class")
#> [1] "model_structure"     "generator_structure" "list"               
```
