# Leaf Reduction Plan

## Goal

Replace the current "full triple everywhere" kernel evaluation with one demand-aware architecture that computes only the channels a query actually needs.

This is aimed at the current bottleneck on the compressed realistic workload for the stop/ignore model:

- total compressed `log_likelihood()`: about `16.2 ms`
- compressed `go` observed: about `1.9 ms`
- compressed `stop` observed: about `12.1 ms`
- compressed `stop NA`: about `2.1 ms`

The target is to make the compressed observed-stop path materially faster than the current code for this model.

## Non-goals

- No model-specific fast paths.
- No second execution architecture kept alive.
- No compatibility branch that preserves the old global-triple kernel path.
- No partial refactor that only changes one call site while the rest still assumes full triples.

## Current Problem

The current kernel machinery is structurally overcomputing.

- `eval_node_recursive_dense()` forces `density`, `survival`, and `cdf` all to `true` in `src/distributions.cpp`.
- `kernel_event_eval_callback()` ignores demand and calls `eval_event_ref_idx(..., kEvalAll, ...)`.
- `KernelRuntimeState::slots` stores a plain `KernelNodeValues` triple with no validity mask.
- Fused observed-density evaluation in `node_density_with_competitors_from_state()` still requests full triples for target and competitors, then only consumes target density and competitor survival.

So the architecture cannot represent "this slot has only survival" or "this leaf only needs cdf". The code compensates by computing everything.

## Single Architecture

The clean design is:

- one structural `KernelProgram`
- one demand-aware `KernelQueryPlan`
- one runtime cache with per-slot validity masks
- one executor that uses the query plan to decide which channels each slot must produce

Nothing else stays alive.

## Concrete Data Structure Changes

### 1. Replace `KernelEvalNeed`

Use `EvalNeed` everywhere. Remove the duplicate kernel-specific demand type.

Current split:

- `EvalNeed` in `src/distributions.cpp`
- `KernelEvalNeed` in `src/kernel_executor.h`

That split is pointless now. One demand type is enough.

### 2. Add per-slot validity

Replace the runtime slot payload with:

```cpp
struct KernelSlotState {
  KernelNodeValues values;
  std::uint8_t valid_mask{0};
};
```

Then change:

- `KernelRuntimeState::slots` in `src/kernel_executor.h`

from:

```cpp
std::vector<KernelNodeValues> slots;
```

to:

```cpp
std::vector<KernelSlotState> slots;
```

This is the minimum change that makes partial evaluation sound.

### 3. Add compiled query plans

Add:

```cpp
struct KernelQueryPlan {
  std::vector<std::uint8_t> slot_need_mask;
  std::vector<int> target_slots;
  int max_slot{-1};
  bool valid{false};
};
```

to `src/context.h` near the existing kernel structs.

This plan is compiled from:

- the structural `KernelProgram`
- the query roots
- the root channel demand

## Backward Demand Rules

These rules define the new architecture.

- `Event(mask)`: terminal
- `Guard(D) -> ref(D), blocker(S)`
- `Guard(S|C) -> guard survival/CDF routine`
- `And(D) -> all children get D|C`
- `And(S|C) -> all children get C`
- `Or(D) -> all children get D|S`
- `Or(S|C) -> all children get S`
- `Not(C|S) -> child gets C`
- `Not(D) -> invalid request`

These rules match the current composition math in `src/kernel_executor.cpp`.

## Exact Motivating Query

For observed stop-component `A` in the motivating model:

- target root `A_out` needs `D`
- competitor root `B_out` needs `S`

This implies mixed, but not full, demand through the tree.

That is the key point of the refactor:

- mixed demand is real
- full triple everywhere is not

## File-Level Cut Sequence

### Step 1. Move demand type to one place

Files:

- `src/distributions.cpp`
- `src/kernel_executor.h`

Work:

- move `EvalNeed` out of `src/distributions.cpp` into a shared header, likely `src/kernel_executor.h` or a new tiny header
- delete `KernelEvalNeed`
- add helper conversions:
  - `need_mask(EvalNeed)`
  - `needs_density(EvalNeed)`
  - `needs_survival(EvalNeed)`
  - `needs_cdf(EvalNeed)`

Done when:

- all demand-carrying APIs use `EvalNeed`
- there is no second demand type left

### Step 2. Make event evaluation demand-aware

Files:

- `src/kernel_executor.h`
- `src/distributions.cpp`

Work:

- change `KernelEventEvalFnPtr` to accept `EvalNeed`
- change `KernelEventEvaluator::operator()` to pass demand
- change `kernel_event_eval_callback()` to call:

```cpp
eval_event_ref_idx(..., need, ...)
```

instead of `kEvalAll`

This is mandatory. Without it, leaves still overcompute and the refactor is mostly fake.

### Step 3. Add slot validity masks

Files:

- `src/kernel_executor.h`
- `src/kernel_executor.cpp`

Work:

- replace `KernelNodeValues` slots with `KernelSlotState`
- teach `reset_kernel_runtime()` to initialize `valid_mask = 0`
- teach invalidation to clear validity from the invalidated slot onward
- teach cache hits to verify that the requested channels are already valid before reusing a slot

Done when:

- executor can safely reuse a slot computed for `S` and later extend it to `D|S` without lying about `C`

### Step 4. Add the query-plan compiler

Files:

- `src/kernel_program.h`
- `src/kernel_program.cpp`

Work:

- keep `compile_kernel_program()` as the structural topo compiler
- add `compile_kernel_query_plan(...)`
- query-plan inputs should be:
  - structural `KernelProgram`
  - root node indices or root slots
  - root demand masks
- propagate backward over topo slots using the demand rules above
- store one aggregated demand mask per slot
- compute `max_slot`

This remains one architecture because the executor always runs a structural program plus a compiled demand plan.

### Step 5. Replace executor API with plan-based execution

Files:

- `src/kernel_executor.h`
- `src/kernel_executor.cpp`

Work:

- delete the API shape that accepts one global demand triple for the whole execution
- replace it with:

```cpp
bool eval_kernel_plan_incremental(
    const KernelProgram& program,
    const KernelQueryPlan& plan,
    KernelRuntimeState& runtime,
    const KernelEventEvaluator& event_eval,
    const KernelGuardEvaluator& guard_eval,
    std::vector<KernelNodeValues>& out_values);
```

- inside the executor, each slot is evaluated against `plan.slot_need_mask[slot]`
- child composition uses the per-slot mask, not a single global mask

Done when:

- there is only one kernel execution API left
- it is plan-based

### Step 6. Update guard evaluation to use real demand

Files:

- `src/distributions.cpp`

Work:

- `kernel_guard_eval_callback()` already branches on demand, but now it must trust the real mask
- keep:
  - `D` computes `ref_density * blocker_survival`
  - `S|C` computes `guard_cdf_internal(...)`
- do not compute `guard_cdf_internal(...)` if only `D` is requested

This step should save real work on observed stop, because many guard nodes in the target branch only need `D`, and many competitor branches only need `S`.

### Step 7. Compile one observed-density plan

Files:

- `src/distributions.cpp`

Work:

- in `node_density_with_competitors_from_state()`, stop requesting full triples for target and competitors
- build one query plan:
  - target root gets `D`
  - each competitor root gets `S`
- execute that plan
- consume:
  - target density
  - competitor survivals

This directly attacks the current observed-stop bottleneck.

### Step 8. Use the same plan architecture for competitor clusters

Files:

- `src/distributions.cpp`

Work:

- replace the ad hoc batch kernel call in `compute_competitor_cluster_value()`
- compile a plan where every cluster member root gets `S`
- execute via the same plan-based executor

That keeps one architecture. No special case remains.

### Step 9. Use the same plan architecture for generic single-node evaluation

Files:

- `src/distributions.cpp`

Work:

- replace `eval_node_recursive_dense()` and any direct kernel call sites with tiny wrappers that compile or fetch the right query plan
- remove the old hard-coded:

```cpp
kernel_need.density = true;
kernel_need.survival = true;
kernel_need.cdf = true;
```

This is the actual deletion of the old architecture.

## Caching Strategy

Use one runtime cache model only:

- values per slot
- valid-mask per slot

Use one query-plan cache only:

- keyed by:
  - context
  - query kind
  - root slot set
  - root demand masks

No legacy slot cache and no legacy full-triple plan cache.

## What Must Be Deleted

These are explicit cleanup requirements:

- `KernelEvalNeed`
- any kernel entry point that accepts one global full-triple request
- any call site that forces `density = survival = cdf = true`
- any event callback path that calls `eval_event_ref_idx(..., kEvalAll, ...)`

If any of those remain, the refactor is incomplete.

## Expected Performance Outcome

The win should come from the compressed realistic workload, not the fake duplicated-row benchmark.

Current reference numbers on the motivating model:

- compressed total: `~16.2 ms`
- compressed stop observed: `~12.1 ms`
- compressed stop `NA`: `~2.1 ms`

Success means:

- compressed total beats `~16.2 ms`
- compressed stop observed beats `~12.1 ms`
- profile share for:
  - `kernel_guard_eval_callback`
  - `eval_pdf_single`
  - `eval_cdf_single`
  drops on the observed-stop slice

Honest target:

- `15%` to `30%` end-to-end improvement on the compressed realistic benchmark is plausible
- less than `10%` means the churn likely was not worth it

## Verification

Minimum verification set:

1. Correctness:
- full `devtools::test()`

2. Structural:
- no remaining references to `KernelEvalNeed`
- no remaining `kEvalAll` call from the event callback
- no remaining kernel path that globally requests full triples

3. Performance:
- rerun the simulated `2000`-trial compressed benchmark for the motivating model
- rerun stop-observed profile

4. Regression check:
- benchmark a simple `go` model and ensure this refactor does not undo the earlier component-projection win

## Order Constraint

This should be implemented in this order:

1. unify demand type
2. demand-aware event callback
3. slot validity masks
4. query-plan compiler
5. plan-based executor
6. observed-density path on plans
7. competitor-cluster path on plans
8. delete old global-triple execution path

Do not keep both systems in parallel for long. The only acceptable temporary overlap is during the edit sequence before the old entry points are deleted in the same refactor branch.
