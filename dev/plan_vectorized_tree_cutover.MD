# Vectorized Tree Cutover Contract (No Spin)

Last updated: 2026-03-23

## Brutal Status Right Now

1.  The codebase still does not have one vectorized runtime end-to-end. It has one shared lane container and one shared coupling evaluator, but ranked and exact still retain their own scheduler shells.
2.  The old scalar fallback machinery is mostly gone from the `cpp_loglik()` path. The remaining cost is not obvious scalar fallback branches; it is scheduler/orchestration work that still re-groups, re-materializes state, and re-walks logic.
3.  Coupling no longer dispatches by the old payload families, but it still carries transitional metadata in [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h) and a scalar adaptive segment scheduler in [src/coupling_eval.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_eval.cpp).
4.  Ranked no longer uses `RankedTransitionCompiler`, `RankedTransitionTemplate`, or `RankedTransitionStep`, but it still executes through `RankedBatchPlan`, `RankedProgramSliceRef`, and `sequence_prefix_density_batch_resolved(...)` in [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h) and [src/ranked_transitions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.cpp), not the canonical `VectorProgram`.
5.  Exact now uses the shared lane container, but guard evaluation and scenario orchestration are still structurally exact-specific around [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp), especially `exact_fill_guard_prepared_batch(...)` and `exact_eval_guard_batch_common(...)`.
6.  Native mixture and parameter-set evaluation still retain scalar outer orchestration in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp).
7.  Honest summary: current code is a successful scalar fallback deletion plus partial runtime unification. The remaining work is scheduler/runtime consolidation, not deleting more obvious scalar helper names.

## Primary Objective

1.  Build one canonical lowered execution program that both coupling and ranked execution consume.
2.  Make vectorization a property of the runtime itself, not a helper layer around scalar semantics.
3.  Delete old runtime splits instead of wrapping them.
4.  Treat temporary regressions as acceptable only when they buy real structural deletion in the same batch.
5.  Recover speed only after the runtime is genuinely consolidated.

## Non-Negotiables

1.  No adapter shells that preserve the old evaluator bodies under a cleaner API.
2.  No "unified" entrypoint that still dispatches to old `Pair`, `NWay`, and `GenericNodeIntegral` execution bodies.
3.  No ranked "batching" that still loops state-by-state through scalar transition logic.
4.  No tests or benchmark scripts will be allowed to block deletion of dead code purely because they assert stale counters or stale symbol names.
5.  When an old runtime family is replaced, its stats, tests, and dev scripts must be rewritten or deleted in the same phase.

## What Counts As Success

1.  The hot path walks one lowered program over many active lanes.
2.  Scheduler code only groups lanes, applies state deltas, and reduces results.
3.  Runtime op evaluation does not care whether the work came from coupling or ranked execution.
4.  SIMD becomes possible because the leaves operate over contiguous lane arrays with low divergence.

## What Does Not Count

1.  Keeping `OutcomeCouplingProgram` as a permanent runtime alias and calling the result "one tree."
2.  Keeping `RankedBatchPlan` and `RankedProgramSliceRef` as a permanent alternate runtime instead of lowering ranked work into the canonical `VectorProgram`.
3.  Adding a `Vectorized*` wrapper that still re-enters old scalar control flow per state, per segment, or per transition.
4.  Preserving old runtime stats and tests by keeping dead code reachable.
5.  Declaring victory because correctness is intact while the runtime still has multiple evaluator families.

## Current Split Inventory

### Coupling Split

1.  Runtime core that is already shared:
    1.  `uuber::VectorProgram`
    2.  `uuber::VectorLaneBatch`
    3.  `CouplingProgramInterpreter`
2.  Remaining coupling-specific residue:
    1.  `OutcomeCouplingProgram` alias in [src/coupling_program.h](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.h)
    2.  `VectorProgramPattern` transitional metadata in [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h)
    3.  adaptive segment scheduling in [src/coupling_eval.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_eval.cpp)
3.  Honest read:
    1.  evaluator-family split is gone
    2.  scheduler split and transitional naming are not

### Ranked Split

1.  Ranked scheduler/runtime family that still remains:
    1.  `RankedBatchPlan`
    2.  `RankedProgramSliceRef`
    3.  `RankedBatchPlanner`
2.  Runtime hot path:
    1.  `sequence_prefix_density_batch_resolved(...)`
3.  Honest read:
    1.  the old ranked template/compiler runtime is gone
    2.  lower-level exact pieces are batched
    3.  ranked still does not execute the canonical `VectorProgram` directly
    4.  frontier grouping, scenario collection, and collapse are still ranked-specific scheduler work

### Test And Tooling Split

These files currently assert or expose the old split vocabulary and must be rewritten or deleted when the architecture moves:

1.  [tests/testthat/test_exact_scalar_density.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_exact_scalar_density.R)
2.  [tests/testthat/test_shared_trigger_param_matrix.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_shared_trigger_param_matrix.R)
3.  [tests/testthat/test_step5_vectorized_coupling.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_step5_vectorized_coupling.R)
4.  [tests/testthat/test_ranked_cutover_unified.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_ranked_cutover_unified.R)
5.  [src/runtime_stats.h](/Users/nstevenson/Documents/2025/AccumulatR/src/runtime_stats.h)
6.  [src/runtime_stats.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/runtime_stats.cpp)

These are not blockers to deletion. They are cleanup work that follows deletion.

## Target Runtime Architecture

### Core Principle

1.  Lower once.
2.  Walk once per compatible lane group.
3.  Only leaf numeric kernels do real math.
4.  Ranked and coupling differ only in scheduler inputs and reducer logic.

### Proposed Runtime Types

These names are concrete suggestions, not sacred names. The structure is what matters.

#### 1. Canonical Lowered Program

File target: [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)

```cpp
enum class VectorOpCode : std::uint8_t {
  EventDensity = 0,
  EventCDF = 1,
  EventSurvival = 2,
  And = 3,
  Or = 4,
  Not = 5,
  Guard = 6,
  Multiply = 7,
  Complement = 8,
  Accumulate = 9
};

struct VectorOp {
  VectorOpCode code;
  int dst_slot;
  int child_begin;
  int child_count;
  int ref_slot;
  int blocker_slot;
  int payload_idx;
  std::uint32_t flags;
};

struct VectorProgramOutputs {
  std::vector<int> node_idx_to_slot;
  std::vector<int> outcome_idx_to_slot;
  std::vector<int> coupling_op_idx_to_slot;
};

struct VectorProgram {
  std::vector<VectorOp> ops;
  std::vector<int> children;
  VectorProgramOutputs outputs;
  bool valid{false};
};
```

Rules:

1.  This is the only runtime program representation.
2.  `Pair`, `NWay`, and `GenericNodeIntegral` may exist only as builder-time lowering cases.
3.  Ranked transition plans may reference program slices or output slots, but not alternate runtime program families.

#### 2. Lane Batch State

File target: [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h) or a new [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h)

```cpp
struct VectorLaneBatch {
  std::vector<double> t;
  std::vector<double> upper;
  std::vector<double> weight;
  std::vector<const uuber::TrialParamsSoA *> trial_params_soa;
  std::vector<uuber::BitsetState> forced_complete_bits;
  std::vector<uuber::BitsetState> forced_survive_bits;
  std::vector<std::uint8_t> forced_complete_valid;
  std::vector<std::uint8_t> forced_survive_valid;
  std::vector<TimeConstraintMap> time_constraints;
  std::vector<int> group_id;
  std::vector<std::uint8_t> active_mask;
};
```

Rules:

1.  This is the thing that walks the program.
2.  `ExactScenarioBatch` in [src/exact_outcome_density.h](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.h) should be absorbed into, or made layout-compatible with, this lane batch.
3.  Do not keep a second near-identical batch container for ranked execution.

#### 3. Runtime Scratch

File target: new [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h)

```cpp
struct VectorEvalScratch {
  std::vector<double> slot_values;
  std::vector<double> tmp_values;
  std::vector<int> partition_offsets;
  std::vector<int> partition_indices;
  std::vector<std::uint8_t> partition_mask;
};
```

Rules:

1.  Scratch belongs to the executor, not to coupling or ranked specifically.
2.  Reuse scratch aggressively. Avoid per-op allocation.

#### 4. Ranked Scheduler Outputs

File target: [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h)

```cpp
struct RankedStateDelta {
  int trigger_bit{-1};
  int invalidate_slot{0};
  int source_mask_begin{-1};
  int source_mask_count{0};
};

struct RankedProgramSliceRef {
  int density_slot{-1};
  int cdf_slot{-1};
  int survival_slot{-1};
};

struct RankedBatchPlan {
  std::vector<RankedStateDelta> deltas;
  RankedProgramSliceRef slice;
  bool deterministic{false};
  bool valid{false};
};
```

Rules:

1.  Ranked compilation emits state deltas plus references into `VectorProgram`.
2.  Ranked compilation must stop emitting execution templates that imply their own runtime evaluator.

## Existing Types: Reuse, Demote, Delete

### Reuse

1.  `uuber::IrNode`, `uuber::IrEvent`, `uuber::IrOutcome`, and `uuber::IrContext` in [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)
2.  `KernelStateGraph` and competitor-cache machinery in [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)
3.  `TrialParamsSoA` in [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)
4.  `ExactScenarioBatch` layout ideas in [src/exact_outcome_density.h](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.h)

### Demote To Builder-Only

1.  `IrOutcomeCouplingKind`
2.  any "shape" enums used only while lowering source semantics to `VectorProgram`
3.  `for_each_sequence_node_transition(...)` if it remains as a compiler helper only

### Delete From Runtime

1.  `OutcomeCouplingProgram`
2.  `OutcomeCouplingOpKind`
3.  `OutcomeCouplingPairPayload`
4.  `OutcomeCouplingNWayPayload`
5.  `OutcomeCouplingGenericNodeIntegralPayload`
6.  `RankedTransitionStep`
7.  `RankedTransitionTemplate`
8.  `RankedNodeTransitionPlan`
9.  `RankedTransitionCompiler`

The rule is simple: if a type exists mainly to preserve a separate coupling or ranked runtime, it must die.

## File-By-File Destination

### [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)

1.  Add the canonical `VectorProgram` definitions.
2.  Keep `KernelProgram` temporarily only if needed for migration.
3.  Delete or demote `IrOutcomeCouplingKind` once lowering is complete.

### [src/coupling_program.h](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.h)

1.  Stop defining runtime payload families.
2.  Become a builder header that exposes lowering into `VectorProgram`.
3.  End state: no runtime `OutcomeCouplingProgram` alias either.

### [src/coupling_program.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.cpp)

1.  Lower `Pair`, `NWay`, and `GenericNodeIntegral` into `VectorProgram` ops.
2.  Builder-only code may branch by source shape.
3.  Runtime-visible payload families must disappear.

### [src/coupling_eval.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_eval.cpp)

1.  Replace shape dispatch with one evaluator over `VectorProgram`.
2.  Keep adaptive scheduling if necessary during migration.
3.  Adaptive scheduling may choose lane groups, but it may not choose different evaluator families.

### [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h)

1.  The old ranked templates are already gone. The next step is to replace `RankedProgramSliceRef` with canonical `VectorProgram` slot refs or equivalent direct program references.
2.  Keep ranked-specific metadata only for scheduler decisions and state deltas.
3.  If a compiler remains, it must compile scheduler metadata only, not an alternate execution contract.

### [src/ranked_transitions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.cpp)

1.  Rewrite `sequence_prefix_density_batch_resolved(...)` so the hot loop becomes:
    1.  group compatible frontier states
    2.  build lane batch
    3.  evaluate the canonical `VectorProgram` over those lanes
    4.  reduce back into ranked state collapse
2.  Ranked may keep frontier/scenario logic, but it may not keep a separate runtime evaluation contract.
3.  Accept that ranked state expansion remains ranked-specific. That is scheduler logic, not evaluator logic.

### [src/exact_outcome_density.h](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.h)

1.  Make `ExactScenarioBatch` layout-compatible with the new lane batch, or merge them.
2.  Do not preserve a second scenario container with the same fields under a different name.
3.  Exact helpers should fill lane batches and slot arrays, not invent alternate execution contracts.

### [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp)

1.  Keep exact scenario construction only as data preparation.
2.  Feed prepared scenarios into the shared vector executor.
3.  Reduce exact-private orchestration around guards and shared-state partitioning. In particular, `exact_fill_guard_prepared_batch(...)` and `exact_eval_guard_batch_common(...)` must stop acting like a private execution shell.
4.  Do not let exact retain a side runtime that coupling and ranked do not share.

### [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp)

1.  `cpp_loglik()` and native helpers should become scheduler/dispatcher code only.
2.  No more private runtime families hidden behind "entry" helpers.
3.  Entry points should route to the shared vector executor plus thin reducers.
4.  `native_trial_mixture_internal_idx(...)` and `cpp_loglik_multiple(...)` are still valid cleanup targets after the ranked/exact scheduler work is done.

### [src/runtime_stats.h](/Users/nstevenson/Documents/2025/AccumulatR/src/runtime_stats.h) and [src/runtime_stats.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/runtime_stats.cpp)

Replace current split-specific counters with generic ones:

1.  `vector_program_eval_calls`
2.  `vector_program_lane_total`
3.  `vector_program_partition_calls`
4.  `ranked_frontier_group_calls`
5.  `ranked_state_collapse_calls`
6.  `coupling_segment_group_calls`

Do not keep counters whose main purpose is proving old split paths still exist.

## Hard Milestones

### VT0: Freeze The Baseline

Deliverables:

1.  Run [dev/scripts/check_loglik_golden.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_loglik_golden.R)
2.  Run [dev/scripts/bench_centralized.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/bench_centralized.R)
3.  Record current split inventory in the plan and in the batch report

Counts:

1.  the benchmark and correctness baselines are frozen
2.  the split inventory is explicit

Does not count:

1.  saying "we know roughly where the split is"

### VT1: Define The Canonical Runtime

Main code changes:

1.  Add `VectorProgram`, `VectorLaneBatch`, and runtime scratch.
2.  Lower coupling shapes into the canonical program.
3.  Define ranked program-slice refs and state-delta outputs.
4.  Make `ExactScenarioBatch` layout-compatible with the lane batch.
5.  It is acceptable for transitional metadata to remain here, but only if it is explicitly scheduled for deletion later.

Files:

1.  [src/context.h](/Users/nstevenson/Documents/2025/AccumulatR/src/context.h)
2.  [src/coupling_program.h](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.h)
3.  [src/coupling_program.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.cpp)
4.  [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h)
5.  [src/exact_outcome_density.h](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.h)

Deletion requirement:

1.  no runtime code may be added that depends on the old payload families

Counts:

1.  one canonical runtime contract exists
2.  ranked and coupling can target it
3.  this phase defines the contract; it does not prove that all hot paths truly execute it yet

Does not count:

1.  adding the new contract while writing all new runtime code against the old one

### VT2: Cut Coupling Over To One Evaluator

Main code changes:

1.  Rewrite `evaluate_outcome_coupling_unified(...)` and `evaluate_outcome_coupling_unified_batch(...)` to consume only `VectorProgram`.
2.  Delete runtime branches over `OutcomeCouplingOpKind`.
3.  Delete the old payload-specific evaluator bodies.
4.  Keep adaptive segment grouping only as scheduler logic.
5.  Transitional naming or metadata may remain temporarily, but only if the evaluator body split is actually gone.

Files:

1.  [src/coupling_eval.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_eval.cpp)
2.  [src/coupling_program.h](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.h)
3.  [src/coupling_program.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.cpp)

Deletion requirement:

1.  `OutcomeCouplingProgram`
2.  `OutcomeCouplingOpKind`
3.  `OutcomeCouplingPairPayload`
4.  `OutcomeCouplingNWayPayload`
5.  `OutcomeCouplingGenericNodeIntegralPayload`

Honest note:

1.  Deleting the payload families and `OutcomeCouplingOpKind` is VT2-complete.
2.  Deleting the `OutcomeCouplingProgram` alias and `VectorProgramPattern` metadata is still VT5 work if they survive as transitional names only.

Counts:

1.  one coupling evaluator remains
2.  runtime shape dispatch is gone

Does not count:

1.  a renamed dispatcher
2.  a `switch` hidden inside a helper

### VT3: Lower Ranked Execution To Program Slices

Main code changes:

1.  Replace ranked transition templates with `RankedBatchPlan` plus `RankedProgramSliceRef`.
2.  Remove scalar node-evaluation semantics from the ranked compiler.
3.  Make the ranked compiler emit only:
    1.  admissible state metadata
    2.  state deltas
    3.  slice refs into `VectorProgram`
4.  This phase is about deleting the old ranked runtime contract, not proving the ranked hot loop is fully unified yet.

Files:

1.  [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h)
2.  [src/ranked_transitions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.cpp)
3.  [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp)
4.  [src/exact_outcome_density.h](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.h)

Deletion requirement:

1.  `RankedTransitionStep`
2.  `RankedTransitionTemplate`
3.  `RankedNodeTransitionPlan`
4.  `RankedTransitionCompiler`

Counts:

1.  ranked plans no longer define their own runtime evaluator
2.  ranked now targets the shared program

Does not count:

1.  compiling a new object that still drives the same scalar hot loop

### VT4: Replace Ranked Hot Execution

Main code changes:

1.  Rewrite `sequence_prefix_density_batch_resolved(...)`.
2.  The new hot path must be:
    1.  build frontier
    2.  group compatible states
    3.  materialize a lane batch
    4.  evaluate the canonical `VectorProgram` directly once per group
    5.  reduce results back into next-state collapse
3.  Scenario collection may remain ranked-specific, but evaluation cannot.
4.  `RankedProgramSliceRef` is allowed only as temporary scheduler metadata. It does not count as the final runtime contract.

Files:

1.  [src/ranked_transitions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.cpp)
2.  [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp)
3.  [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp)
4.  [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h)

Deletion requirement:

1.  no hot ranked path may execute per-transition scalar evaluation semantics
2.  no hot ranked path may retain a separate execution vocabulary that the coupling evaluator does not understand

Counts:

1.  ranked hot execution shares the runtime with coupling at the evaluator level, not just at the container level
2.  grouped state execution is visible and real
3.  ranked scheduler metadata no longer acts like a private execution shell

Does not count:

1.  looping over states and calling the same old evaluator behind a batch API
2.  lane-batching around `RankedProgramSliceRef` plus exact-private orchestration while still avoiding direct `VectorProgram` execution

### VT5: Delete The Old Architecture

Main code changes:

1.  Remove all transitional runtime aliases and metadata that only preserve the old split vocabulary.
2.  Collapse exact-private execution shells that still sit on top of the shared lane runtime.
3.  Remove all dead split-specific runtime stats.
4.  Remove all tests whose only purpose is to assert the old split vocabulary.
5.  Remove stale dev-script references to old counters or old runtime families.

Primary remaining targets:

1.  `OutcomeCouplingProgram` alias in [src/coupling_program.h](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_program.h)
2.  `VectorProgramPattern` in [src/vector_runtime.h](/Users/nstevenson/Documents/2025/AccumulatR/src/vector_runtime.h)
3.  `RankedProgramSliceRef` as a runtime-visible execution concept in [src/ranked_transitions.h](/Users/nstevenson/Documents/2025/AccumulatR/src/ranked_transitions.h)
4.  `exact_fill_guard_prepared_batch(...)` in [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp)
5.  `exact_eval_guard_batch_common(...)` in [src/exact_outcome_density.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/exact_outcome_density.cpp)
6.  `guard_cdf_internal_prepared(...)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp)

Required proof:

1.  grep inventory of deleted symbols
2.  hot-path call inventory
3.  proof that coupling and ranked both target `VectorProgram`

Counts:

1.  the codebase shrinks or stays flat through deletion

Does not count:

1.  marking old code deprecated but still reachable

### VT6: Speed Recovery On The Real Runtime

Allowed work:

1.  improve ranked scheduler grouping after VT4 is structurally true
2.  improve exact guard batching after VT5 removes the private execution shell
3.  improve adaptive coupling segment scheduling
4.  improve lane partitioning
5.  improve scratch reuse
6.  improve leaf kernels
7.  add SIMD-friendly loops

Forbidden work:

1.  resurrecting shape-specific evaluators
2.  reintroducing scalar ranked engines for speed
3.  preserving transitional architecture because it benchmarks well on a microcase

### VT7: Closure

Required:

1.  correctness stable
2.  centralized benchmark stable
3.  structural deletion proof still true
4.  remaining scheduler-only code is honestly scheduler-only, not a hidden evaluator family

## Test And Script Policy

### Tests To Keep

These stay because they guard behavior:

1.  [dev/scripts/check_loglik_golden.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_loglik_golden.R)
2.  [dev/scripts/bench_centralized.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/bench_centralized.R)
3.  behavior-level ranked, coupling, and multi-outcome tests that do not require the old runtime split vocabulary

### Tests To Rewrite Or Delete During Cutover

These are allowed to change or disappear if they block deletion of stale architecture:

1.  [tests/testthat/test_exact_scalar_density.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_exact_scalar_density.R)
    1.  stop asserting `exact_batch_density_calls`
    2.  rewrite around correctness plus new generic vector stats if still useful
2.  [tests/testthat/test_shared_trigger_param_matrix.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_shared_trigger_param_matrix.R)
    1.  stop asserting `shared_trigger_mask_batch_calls` and direct-group path names if those names become stale
3.  [tests/testthat/test_step5_vectorized_coupling.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_step5_vectorized_coupling.R)
    1.  stop asserting old `cpp_loglik_*_group_batch_*` counters if those counters preserve dead runtime categories
4.  [tests/testthat/test_ranked_cutover_unified.R](/Users/nstevenson/Documents/2025/AccumulatR/tests/testthat/test_ranked_cutover_unified.R)
    1.  stop asserting the current ranked batch vocabulary if ranked no longer has a separate runtime family

Rule:

1.  No test gets veto power over deletion of dead architecture because it asserted the dead architecture.

### Runtime Stats Policy

1.  During VT1-VT4, add only generic vector-runtime stats if genuinely needed.
2.  Delete split-specific stats as soon as the split they describe is gone.
3.  Do not keep counters as historical monuments.

## Verification Rules

### Correctness Gate

1.  Run [dev/scripts/check_loglik_golden.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/check_loglik_golden.R)
2.  Required: `max_abs_diff <= 1e-4`
3.  Any failure blocks the phase.

### Benchmark Gate

1.  Run [dev/scripts/bench_centralized.R](/Users/nstevenson/Documents/2025/AccumulatR/dev/scripts/bench_centralized.R)
2.  Use fresh installs for both baseline and candidate.
3.  Temporary regressions are allowed only through VT1-VT4 and only when they accompany real deletion.
4.  VT6 and VT7 require non-regression relative to the cutover lineage baseline.

### Structural Gate

At closure, all must be true:

1.  coupling runtime no longer dispatches by shape family
2.  ranked runtime no longer owns a separate transition-evaluation engine
3.  one canonical lowered program is shared
4.  remaining ranked-specific logic is scheduler and collapse logic only
5.  no tests or scripts require dead split code to remain

## Remaining Real Work

1.  Ranked still needs to execute the canonical `VectorProgram` directly from `sequence_prefix_density_batch_resolved(...)`, not through `RankedProgramSliceRef` plus exact-node orchestration.
2.  Exact guard evaluation still behaves like a private execution shell around `exact_fill_guard_prepared_batch(...)`, `exact_eval_guard_batch_common(...)`, and `guard_cdf_internal_prepared(...)`.
3.  Coupling still retains scalar adaptive segment scheduling in [src/coupling_eval.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/coupling_eval.cpp).
4.  `native_trial_mixture_internal_idx(...)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp) still has serial outer orchestration.
5.  `cpp_loglik_multiple(...)` in [src/distributions.cpp](/Users/nstevenson/Documents/2025/AccumulatR/src/distributions.cpp) is still scalar over parameter sets.

## Reporting Format For Every Meaningful Batch

1.  commit SHA
2.  UTC timestamp
3.  milestone label `VT0` through `VT7`
4.  correctness result
5.  benchmark summary
6.  deleted symbols
7.  what is still not unified
8.  whether any tests/scripts were rewritten or removed to permit cleanup

## Current Honest Status Board

1.  `VT0`: done enough for lineage, rerun only when a new benchmark baseline is needed
2.  `VT1`: mostly complete; canonical runtime types exist, but transitional metadata still remains
3.  `VT2`: evaluator-body cutover complete for coupling; scheduler-level cleanup still remains
4.  `VT3`: complete for deleting the old ranked template/compiler runtime
5.  `VT4`: not complete; ranked still runs through `RankedBatchPlan` / `RankedProgramSliceRef` and exact-node orchestration rather than direct canonical `VectorProgram` execution
6.  `VT5`: partially complete; stale helper/stat cleanup is largely done, but exact-private execution shells and transitional metadata still remain
7.  `VT6`: not formally started; benchmark-proven speed work is allowed only when it does not preserve dead architecture
8.  `VT7`: not started

## Immediate Next Step

1.  Replace `RankedProgramSliceRef`-driven hot execution inside `sequence_prefix_density_batch_resolved(...)` with direct canonical `VectorProgram` evaluation.
2.  After ranked is truly on the canonical evaluator, collapse the exact guard shell around `exact_fill_guard_prepared_batch(...)`, `exact_eval_guard_batch_common(...)`, and `guard_cdf_internal_prepared(...)`.
3.  Then flatten adaptive coupling scheduling and only after that spend time on `native_trial_mixture_internal_idx(...)` and `cpp_loglik_multiple(...)`.
