# Draw a processing tree for the responses in a model

This is a compact visual summary of the response rules in a model. It
avoids equations and instead shows the observed responses, the
accumulators or pools that feed them, and any blocking relationships.

## Usage

``` r
processing_tree(model, outcome_label = NULL, return_dot = FALSE)
```

## Arguments

- model:

  A race model or finalized model structure.

- outcome_label:

  Optional response label. If supplied, only that response is shown.

- return_dot:

  If \`TRUE\`, return the Graphviz DOT string instead of a plot.

## Value

If \`DiagrammeR\` is available and \`return_dot = FALSE\`, a
\`DiagrammeR\` graph. Otherwise, a list with \`dot\`, \`nodes\`, and
\`edges\`.
