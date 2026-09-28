# Pool several accumulators under a shared label

A pool finishes when its `k`th member finishes. Use the pool label in
[`add_outcome()`](https://niekstevenson.github.io/AccumulatR/reference/add_outcome.md),
another pool, or
[`after()`](https://niekstevenson.github.io/AccumulatR/reference/after.md).

## Usage

``` r
add_pool(spec, id, members, k = 1L)
```

## Arguments

- spec:

  A `race_spec` object.

- id:

  Label for the pool.

- members:

  Accumulator or pool labels included in the pool.

- k:

  Number of members that must finish, from `1` to `length(members)`. The
  default `1` gives the first member's finishing time.

## Value

The updated `race_spec`.

## Examples

``` r
spec <- race_spec()
spec <- add_accumulator(spec, "A", "lognormal")
spec <- add_accumulator(spec, "B", "lognormal")
spec <- add_pool(spec, "P1", members = c("A", "B"), k = 1L)
```
