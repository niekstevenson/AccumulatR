# Draw a processing tree for the responses in a model

Show the observed responses, the accumulators or pools that feed them,
and any blocking relationships as a directed graph.

## Usage

``` r
processing_tree(model, outcome_label = NULL, return_dot = FALSE)
```

## Arguments

- model:

  A finalized model structure.

- outcome_label:

  Optional response label. If supplied, only that response is shown.

- return_dot:

  If `TRUE`, return the graph description as a list instead of rendering
  it with DiagrammeR.

## Value

If `DiagrammeR` is available and `return_dot = FALSE`, a `DiagrammeR`
graph. Otherwise, a list with `dot`, `nodes`, and `edges`.
