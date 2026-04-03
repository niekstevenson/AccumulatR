# Prepare behavioral data for repeated likelihood evaluation

A likelihood context stores the processed behavioral data, model
structure, and compiled backend needed for fast repeated calls to
\`log_likelihood()\`. Build it once, then reuse it while searching over
parameter values.

## Usage

``` r
build_likelihood_context(structure, data_df, prep = NULL, native_bundle = NULL)
```

## Arguments

- structure:

  Finalized model structure.

- data_df:

  Behavioral data. In the simplest case this contains \`trial\`, \`R\`,
  and \`rt\`; for multi-outcome models it can also contain \`R2\`,
  \`rt2\`, and so on.

- prep:

  Optional preprocessed model bundle.

- native_bundle:

  Optional serialized native bundle.

## Value

A \`likelihood_context\` object.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_outcome(spec, "A_win", "A")
structure <- finalize_model(spec)
params_df <- build_param_matrix(
  spec,
  c(A.meanlog = 0, A.sdlog = 0.1, A.q = 0, A.t0 = 0),
  n_trials = 2
)
data_df <- simulate(structure, params_df, seed = 1)
build_likelihood_context(structure, data_df)
#> $structure
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
#> $shared_params
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
#> $prep$components$attrs$`__default__`$guess
#> NULL
#> 
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
#> $prep$special_outcomes
#> list()
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
#>   component_id weight has_weight_param      attrs  mode   reference
#> 1  __default__      1            FALSE NULL, NULL fixed __default__
#> 
#> $special_outcomes
#> list()
#> 
#> $shared_triggers
#> list()
#> 
#> attr(,"class")
#> [1] "model_structure"     "generator_structure" "list"               
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
#> attr(,".lik_id")
#> [1] 1
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
#> $prep$components$attrs$`__default__`$guess
#> NULL
#> 
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
#> $prep$special_outcomes
#> list()
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
#> $prep$.id_index
#>     A A_win 
#>     1     2 
#> 
#> $prep$.label_cache
#> <environment: 0x557bc3cdf6d0>
#> 
#> $prep$.expr_compiled
#> $prep$.expr_compiled$nodes
#> $prep$.expr_compiled$nodes[[1]]
#> $prep$.expr_compiled$nodes[[1]]$id
#> [1] 1
#> 
#> $prep$.expr_compiled$nodes[[1]]$kind
#> [1] "event"
#> 
#> $prep$.expr_compiled$nodes[[1]]$expr
#> $prep$.expr_compiled$nodes[[1]]$expr$kind
#> [1] "event"
#> 
#> $prep$.expr_compiled$nodes[[1]]$expr$source
#> [1] "A"
#> 
#> $prep$.expr_compiled$nodes[[1]]$expr$k
#> NULL
#> 
#> 
#> $prep$.expr_compiled$nodes[[1]]$sources
#> [1] 1
#> 
#> $prep$.expr_compiled$nodes[[1]]$args
#> NULL
#> 
#> $prep$.expr_compiled$nodes[[1]]$reference_id
#> [1] NA
#> 
#> $prep$.expr_compiled$nodes[[1]]$blocker_id
#> [1] NA
#> 
#> $prep$.expr_compiled$nodes[[1]]$unless_ids
#> NULL
#> 
#> $prep$.expr_compiled$nodes[[1]]$arg
#> [1] NA
#> 
#> $prep$.expr_compiled$nodes[[1]]$source
#> [1] "A"
#> 
#> $prep$.expr_compiled$nodes[[1]]$needs_forced
#> [1] TRUE
#> 
#> $prep$.expr_compiled$nodes[[1]]$scenario_sensitive
#> [1] FALSE
#> 
#> 
#> 
#> $prep$.expr_compiled$signatures
#> <environment: 0x557bc2eea060>
#> 
#> 
#> $prep$.competitors
#> $prep$.competitors$A_win
#> list()
#> 
#> attr(,"by_index")
#> [1] TRUE
#> 
#> $prep$.runtime
#> $prep$.runtime$expr_compiled
#> $prep$.runtime$expr_compiled$nodes
#> $prep$.runtime$expr_compiled$nodes[[1]]
#> $prep$.runtime$expr_compiled$nodes[[1]]$id
#> [1] 1
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$kind
#> [1] "event"
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$expr
#> $prep$.runtime$expr_compiled$nodes[[1]]$expr$kind
#> [1] "event"
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$expr$source
#> [1] "A"
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$expr$k
#> NULL
#> 
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$sources
#> [1] 1
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$args
#> NULL
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$reference_id
#> [1] NA
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$blocker_id
#> [1] NA
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$unless_ids
#> NULL
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$arg
#> [1] NA
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$source
#> [1] "A"
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$needs_forced
#> [1] TRUE
#> 
#> $prep$.runtime$expr_compiled$nodes[[1]]$scenario_sensitive
#> [1] FALSE
#> 
#> 
#> 
#> $prep$.runtime$expr_compiled$signatures
#> <environment: 0x557bc2eea060>
#> 
#> 
#> $prep$.runtime$label_cache
#> <environment: 0x557bc3cdf6d0>
#> 
#> $prep$.runtime$competitor_map
#> $prep$.runtime$competitor_map$A_win
#> list()
#> 
#> attr(,"by_index")
#> [1] TRUE
#> 
#> $prep$.runtime$id_index
#>     A A_win 
#>     1     2 
#> 
#> $prep$.runtime$pool_members_cache
#> <environment: 0x557bc3abf928>
#> 
#> $prep$.runtime$cache_bundle
#> $node_plan
#> $node_plan$nodes
#> $node_plan$nodes[[1]]
#> $node_plan$nodes[[1]]$id
#> [1] 1
#> 
#> $node_plan$nodes[[1]]$kind
#> [1] "event"
#> 
#> $node_plan$nodes[[1]]$expr
#> $node_plan$nodes[[1]]$expr$kind
#> [1] "event"
#> 
#> $node_plan$nodes[[1]]$expr$source
#> [1] "A"
#> 
#> $node_plan$nodes[[1]]$expr$k
#> NULL
#> 
#> 
#> $node_plan$nodes[[1]]$sources
#> [1] 1
#> 
#> $node_plan$nodes[[1]]$args
#> NULL
#> 
#> $node_plan$nodes[[1]]$reference_id
#> [1] NA
#> 
#> $node_plan$nodes[[1]]$blocker_id
#> [1] NA
#> 
#> $node_plan$nodes[[1]]$unless_ids
#> NULL
#> 
#> $node_plan$nodes[[1]]$arg
#> [1] NA
#> 
#> $node_plan$nodes[[1]]$source
#> [1] "A"
#> 
#> $node_plan$nodes[[1]]$needs_forced
#> [1] TRUE
#> 
#> $node_plan$nodes[[1]]$scenario_sensitive
#> [1] FALSE
#> 
#> 
#> 
#> $node_plan$signatures
#> <environment: 0x557bc2eea060>
#> 
#> 
#> $precomputed_values
#> <environment: 0x557bc3ac5b30>
#> 
#> $pool_templates
#> <environment: 0x557bc3ac5858>
#> 
#> $guard_quadrature
#> <environment: 0x557bc3ac93b0>
#> 
#> $guard_quadrature_meta
#> <environment: 0x557bc3ac90d8>
#> 
#> $guard_quadrature_limit
#> [1] 128
#> 
#> $native_ctx
#> <environment: 0x557bc3ac0818>
#> 
#> $version
#> [1] "2026-04-03 17:04:59 UTC"
#> 
#> attr(,"class")
#> [1] "likelihood_cache_bundle"
#> 
#> 
#> 
#> $native_ctx
#> <pointer: 0x557bc1ad0100>
#> 
#> $data_df
#>   trial     R        rt accumulator onset
#> 1     1 A_win 0.9679031           A     0
#> 2     2 A_win 0.9198333           A     0
#> 
#> $data_df_cpp
#>   trial R_id        rt onset
#> 1     1    2 0.9679031     0
#> 2     2    2 0.9198333     0
#> 
#> $param_layout
#> $param_layout$row_trial
#> [1] 1 2
#> 
#> $param_layout$row_acc
#> [1] 1 1
#> 
#> $param_layout$row_component
#> [1] NA NA
#> 
#> $param_layout$n_trials
#> [1] 2
#> 
#> $param_layout$trial_ids
#> [1] 1 2
#> 
#> $param_layout$rectangular
#> [1] TRUE
#> 
#> 
#> $n_trials
#> [1] 2
#> 
#> $trial_ids
#> [1] 1 2
#> 
#> $rank_width
#> [1] 1
#> 
#> $outcome_label_id_map
#> A_win 
#>     2 
#> 
#> $component_label_idx_map
#> __default__ 
#>           0 
#> 
#> $rel_tol
#> [1] 1e-05
#> 
#> $abs_tol
#> [1] 1e-06
#> 
#> $max_depth
#> [1] 12
#> 
#> attr(,"class")
#> [1] "likelihood_context"
```
