# Refresh the compiled backend inside a context

If the stored native pointer is no longer valid, this function rebuilds
it from the model specification. Most users will not need to call this
directly, but it can be useful when caching contexts across sessions or
long workflows.

## Usage

``` r
ensure_native_ctx(context, model_spec = NULL)
```

## Arguments

- context:

  A context-like object containing \`native_ctx\` plus either \`prep\`
  or \`structure\$prep\`.

- model_spec:

  Model specification from \`finalize_model()\`. Usually inferred from
  \`context\` when available.

## Value

The same object, with \`native_ctx\` refreshed if needed.
