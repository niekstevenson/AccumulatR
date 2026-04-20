# Not Supported

This file lists model and data shapes that one could try to construct, but that
are not supported.

## 1. `none_of(...)` as a branch inside `first_of(...)` / `or`

Example model:

```r
race_spec() |>
  add_accumulator("go", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("X", first_of("go", none_of("stop")))
```

Status:

- not supported

Reason:

- `none_of(stop)` is a passive truth condition, not an event-time arrival.
- `first_of(...)` / `or` requires a branch to become true first at a specific
  event time.
- Mixing those semantics would require treating absence conditions as
  event-producing branches, which the current exact transition planner does not
  do cleanly.

Boundary:

- [src/eval/exact_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_kernel.hpp:687)

## 2. Ranked (`n_outcomes > 1`) models that are not direct event outcomes

Example models:

```r
race_spec(n_outcomes = 2L) |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_outcome("A", "a") |>
  add_outcome("B_guard", inhibit("b", by = "stop"))
```

```r
race_spec(n_outcomes = 2L) |>
  add_accumulator("a1", "lognormal") |>
  add_accumulator("a2", "lognormal") |>
  add_pool("A", c("a1", "a2")) |>
  add_accumulator("b", "lognormal") |>
  add_outcome("R1", "A") |>
  add_outcome("R2", "b")
```

Status:

- not supported

Reason:

- Ranked likelihood is currently restricted to direct event outcomes so that
  each rank corresponds to one realized event source.
- Guarded or onset-dependent ranked outcomes would require carrying latent
  prerequisite-time distributions across later ranks.
- Pooled, overlapping, or wrapper-based ranked outcomes add extra latent
  structure or ambiguous observation semantics.
- Keeping ranked outcomes direct-only avoids reintroducing ad hoc ranked
  evaluator branches.

Boundary:

- [R/model_definition.R](/Users/nstevenson/Documents/2025/AccumulatR/R/model_definition.R:347)
- [src/eval/exact_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_kernel.hpp:2409)
- [src/eval/exact_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_kernel.hpp:2482)

## 3. Exact conditioning on overlapping forced pool states

Example models:

```r
race_spec() |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal") |>
  add_pool("P", c("a", "b")) |>
  add_outcome("X", all_of("a", "P"))
```

```r
race_spec() |>
  add_accumulator("a", "lognormal") |>
  add_accumulator("b", "lognormal") |>
  add_accumulator("c", "lognormal") |>
  add_pool("P1", c("a", "b")) |>
  add_pool("P2", c("b", "c")) |>
  add_outcome("X", all_of("P1", "P2"))
```

Status:

- not supported

Reason:

- In these shapes, the exact kernel would have to condition simultaneously on a
  pool state and an overlapping substate built from the same members.
- The current forced-state representation is set-based and does not encode a
  coherent joint state for overlapping pools and members.
- Allowing this would mix incompatible constraints instead of computing a
  clean exact conditioning law.

Boundary:

- [src/eval/exact_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_kernel.hpp:1832)
- [src/eval/exact_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/exact_kernel.hpp:1838)

## 4. Finite RT with missing response label

Example:

```r
R = NA
rt = 0.43
```

Status:

- not supported

Boundary:

- [src/eval/observed_kernel.hpp](/Users/nstevenson/Documents/2025/AccumulatR/src/eval/observed_kernel.hpp:456)
