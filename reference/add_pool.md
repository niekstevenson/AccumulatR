# Pool several accumulators under a shared label

Pools let you talk about several accumulators as one source when
defining observed responses.

## Usage

``` r
add_pool(spec, id, members, k = 1L)
```

## Arguments

- spec:

  A \`race_spec\` object.

- id:

  Label for the pool.

- members:

  Accumulator labels included in the pool.

- k:

  Threshold for a \`k\`-of-\`n\` pool rule.

## Value

The updated \`race_spec\`.

## Examples

``` r
spec <- race_spec()
spec <- add_pool(spec, "P1", members = c("A", "B"), k = 1L)
```
