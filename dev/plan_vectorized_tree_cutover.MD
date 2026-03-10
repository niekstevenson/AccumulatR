# Vectorized Tree Cutover Contract (No Spin)

Last updated: 2026-03-10

## Brutal Status Right Now

1.  The codebase has one consolidated engine boundary, not one lowered vectorized tree.
2.  `evaluate_outcome_coupling_unified(...)` is unified only at the entry surface. It still branches by `OutcomeCouplingOpKind` and still carries separate semantic payloads for `Pair`, `NWay`, and `GenericNodeIntegral`.
3.  Ranked evaluation still owns a separate sequence-state scheduler built around `SequenceState`, `SequenceStateKey`, `RankedTransitionCompiler`, `for_each_sequence_node_transition(...)`, and `sequence_prefix_density_resolved(...)`.
4.  The kernel executor is vectorized over time batches, but the whole system is not yet lowered into one vectorized IR that both coupling and ranked execution consume.
5.  Honest summary: this is a successful lower-layer consolidation, not a final vectorized-tree architecture.

## Primary Objective

1.  Build one real lowered execution representation for outcome coupling and ranked execution.
2.  Make vectorization a property of the execution model itself, not a helper wrapped around scalar semantics.
3.  Keep correctness mandatory at every phase.
4.  Accept temporary performance regressions during structural cutover if, and only if, the cutover deletes real architectural splits.
5.  Recover speed only after the lowered tree is real.

## What “Single Vectorized Tree” Actually Means

1.  `Pair`, `NWay`, and `GenericNodeIntegral` no longer exist as separate runtime programs inside `evaluate_outcome_coupling_unified(...)`.
2.  Outcome coupling executes one lowered node program shape, with operation differences expressed as node/op data, not separate evaluator bodies.
3.  Ranked execution no longer evaluates transition semantics by repeated scalar `eval_node_recursive_dense(...)` calls inside `for_each_sequence_node_transition(...)`.
4.  Ranked transitions are compiled into the same lowered execution representation used elsewhere, then evaluated in batches over times and state groups.
5.  The remaining semantic split is only at the scheduler boundary:
    1.  non-ranked work schedules one execution slice
    2.  ranked work schedules many admissible state slices
    3.  both consume the same lowered node program and batch executor

## What Does Not Count

1.  Leaving `evaluate_outcome_coupling_unified(...)` as a dispatch shell over `Pair`, `NWay`, and `GenericNodeIntegral` and calling that “one tree”.
2.  Keeping `OutcomeCouplingPairPayload`, `OutcomeCouplingNWayPayload`, and `OutcomeCouplingGenericNodeIntegralPayload` as parallel runtime program shapes and only renaming the wrapper.
3.  Making ranked call a batch helper while `for_each_sequence_node_transition(...)` still performs scalar state-by-state `eval_node_recursive_dense(...)` work.
4.  Moving duplicate logic to another file without deleting the semantic split.
5.  Adding adapter layers or provider objects that preserve old control flow under a new interface.
6.  Claiming success because `kernel_executor.cpp` is vectorized over time while the main runtime still branches by coupling shape and ranked state path.
7.  Claiming success because benchmarks are green while the structural split remains.

## Automatic Fail Conditions

1.  `OutcomeCouplingOpKind` remains the active runtime split for coupling execution.
2.  `evaluate_outcome_coupling_unified(...)` still contains separate live evaluator bodies for `Pair`, `NWay`, and `GenericNodeIntegral`.
3.  `for_each_sequence_node_transition(...)` still uses scalar `eval_node_recursive_dense(...)` as the main transition-evaluation mechanism.
4.  Ranked batching is implemented by repeatedly calling the same scalar transition evaluator inside a loop over states or times.
5.  New wrappers are added while the old semantic split remains reachable.
6.  The hot file grows because of duplicated replacement bodies rather than shrinking or staying flat through real deletion.

## Required End State

### A. One Lowered Coupling IR

1.  Replace the current coupling-specific runtime payload family with one canonical lowered program representation.
2.  `build_outcome_coupling_program(...)` must produce one node/op graph shape, not one of several runtime payload families.
3.  `Pair`, `NWay`, and `GenericNodeIntegral` become construction-time lowering cases only.
4.  Runtime evaluation consumes only the lowered representation.

### B. One Batch-Native Coupling Executor

1.  The executor consumes batches of times against that lowered coupling IR directly.
2.  Direct CDF / term evaluation becomes one loop over lowered ops, not separate bodies for old coupling shapes.
3.  Non-coupling helper orchestration must stay thin and data-oriented.

### C. Ranked Lowering To The Same IR

1.  `RankedTransitionCompiler` must stop producing a control-flow shape that is then evaluated by scalar recursive node calls.
2.  It must instead emit state deltas and references into the same lowered node program used by non-ranked execution.
3.  `sequence_prefix_density_resolved(...)` may still own ranked scheduling, but it must schedule batches of lowered program evaluations, not scalar transition callbacks.

### D. Batched Ranked Execution

1.  Ranked execution must batch over grouped states and time vectors where semantics permit.
2.  State collapse may remain ranked-specific.
3.  Transition evaluation must stop being “one state, one scalar tree-walk, repeat”.
4.  The main hot path must become “group states, apply state deltas, evaluate lowered program batch, reduce”.

## Functions And Structures That Must Change Meaningfully

### Coupling Surface

1.  `evaluate_outcome_coupling_unified(...)`
   1.  Must stop branching by separate coupling runtime program kinds.
   2.  Must become a thin entry into one lowered-program evaluator.
2.  `build_outcome_coupling_program(...)`
   1.  Must stop returning a variant-like runtime program family.
   2.  Must lower all coupling shapes into one canonical representation.
3.  `OutcomeCouplingOpKind`
   1.  Must disappear as a runtime evaluator split.
   2.  Construction-time op kinds inside a canonical IR are allowed; top-level program-kind dispatch is not.
4.  `OutcomeCouplingPairPayload`, `OutcomeCouplingNWayPayload`, `OutcomeCouplingGenericNodeIntegralPayload`
   1.  Must disappear as active runtime payload types.
   2.  Temporary builder-local lowering structs are allowed only if they do not survive into runtime execution.

### Ranked Surface

1.  `SequenceState`, `SequenceStateKey`, and ranked state collapse remain valid concepts.
2.  `RankedTransitionCompiler`
   1.  Must stop compiling work that can only be evaluated by scalar recursive node calls.
   2.  Must compile state deltas plus lowered-program references.
3.  `for_each_sequence_node_transition(...)`
   1.  Must stop being the main scalar transition executor.
   2.  It may remain only as a builder/scheduler helper if it no longer performs the hot scalar node evaluation.
4.  `sequence_prefix_density_resolved(...)`
   1.  Must stop driving transition evaluation through scalar `eval_node_recursive_dense(...)` calls.
   2.  It may remain as the ranked scheduler if it dispatches grouped lowered-program batches instead.
5.  `eval_node_recursive_dense(...)`
   1.  May remain for legacy-scope or non-ranked scalar edge cases during migration.
   2.  It must no longer be the main ranked transition-evaluation engine at closure.

## Hard Milestones

### VT0: Baseline Freeze

1.  Freeze current behavior with:
    1.  `Rscript dev/scripts/check_loglik_golden.R`
    2.  `Rscript dev/scripts/bench_centralized.R`
2.  Record current structural split inventory:
    1.  coupling program-kind split
    2.  ranked scalar transition evaluation split
3.  This step is complete only when the split inventory is explicit.

### VT1: Canonical Coupling IR Definition

1.  Introduce one canonical lowered coupling program representation.
2.  Lower `Pair`, `NWay`, and `GenericNodeIntegral` into that representation.
3.  Do not yet claim runtime unification until the evaluator stops branching by old program kind.

Counts:
1.  old runtime payload family no longer defines the active evaluator contract
2.  one program representation is passed into the runtime

Does not count:
1.  adding a canonical IR alongside the old runtime program family while both remain live

### VT2: Coupling Runtime Cutover

1.  Make `evaluate_outcome_coupling_unified(...)` evaluate only the canonical lowered program.
2.  Delete old runtime program-kind branches.
3.  Delete old runtime payload-family evaluator bodies.

Counts:
1.  runtime branch by `OutcomeCouplingOpKind` is gone
2.  old payload-specific evaluator code is deleted or unreachable by proof

Does not count:
1.  one wrapper calling three old evaluator bodies
2.  one `switch` over old payload kinds inside a renamed helper

### VT3: Ranked Program Lowering

1.  Redefine ranked transition compilation so ranked transitions reference the canonical lowered program.
2.  Separate ranked-only responsibilities clearly:
    1.  admissible-state expansion
    2.  state delta application
    3.  state collapse
3.  Remove scalar tree evaluation from ranked transition compilation.

Counts:
1.  ranked transition objects carry lowered-program references and state deltas
2.  ranked hot work is no longer defined in terms of scalar recursive node calls

Does not count:
1.  building ranked transition lists that still call scalar node evaluation later

### VT4: Batched Ranked Execution Cutover

1.  Replace scalar transition evaluation in ranked execution with grouped lowered-program batch execution.
2.  Group compatible ranked states before evaluation.
3.  Evaluate the lowered program over time vectors per state group.
4.  Reduce results back into ranked state collapse.

Counts:
1.  the main ranked hot path no longer calls `eval_node_recursive_dense(...)` per transition as its core mechanism
2.  grouped state execution is visible in the runtime structure

Does not count:
1.  looping over states and calling the same scalar routine through a batch wrapper
2.  batching times but not states while the state path stays scalar-per-transition

### VT5: Structural Deletion Proof

1.  Prove the old coupling program-kind split is deleted.
2.  Prove ranked scalar transition evaluation is no longer the hot path.
3.  Prove no wrapper-only duplicate engine remains.

Required proof:
1.  grep/call-graph inventory
2.  hot-path function inventory
3.  deleted symbol inventory

### VT6: Speed Recovery On The Real Tree

1.  Only after VT5 is achieved.
2.  Optimize inside the lowered-tree executor and ranked grouped execution only.
3.  No reintroduction of shape-specific runtime engines.

Counts:
1.  speed wins inside the same lowered execution path

Does not count:
1.  restoring `Pair` / `NWay` / scalar-ranked special engines for speed

### VT7: Final Closure

1.  Correctness passes on the cutover commit lineage.
2.  Canonical benchmark is stable on the cutover commit lineage.
3.  Structural proofs still hold after speed recovery.
4.  Only then can the codebase honestly claim one lowered vectorized tree.

## Verification Rules

### Correctness Gate

1.  `Rscript dev/scripts/check_loglik_golden.R`
2.  Required: `max_abs_diff <= 1e-4`
3.  Any correctness failure blocks the phase.

### Benchmark Gate

1.  `Rscript dev/scripts/bench_centralized.R`
2.  The benchmark file of record is `dev/scripts/scratch_outputs/benchmark_centralized.csv`.
3.  Temporary regressions are allowed during VT1 through VT4 only if structural deletion occurs in the same batch.
4.  VT6 and VT7 require stable non-regression relative to the frozen cutover baseline for the same lineage.

### Structural Gate

At closure, all must be true:

1.  `evaluate_outcome_coupling_unified(...)` is not a shape-dispatch shell.
2.  Active runtime execution no longer depends on separate `Pair`, `NWay`, and `GenericNodeIntegral` program bodies.
3.  Ranked hot execution no longer depends on scalar `eval_node_recursive_dense(...)` transition evaluation.
4.  The canonical lowered program is shared between coupling and ranked execution.
5.  Remaining ranked-specific logic is only scheduler/state-collapse logic, not a separate node-evaluation engine.

## Reporting Format (Required Every Meaningful Batch)

1.  Commit SHA
2.  UTC timestamp
3.  Milestone: `VT0` through `VT7`
4.  Correctness result (`max_abs_diff`)
5.  Benchmark summary from the canonical benchmark
6.  Explicit deleted symbols / deleted split points
7.  Explicit statement of what still does not count as unified

## Current Honest Status Board

1.  `VT0`: not frozen yet under this contract
2.  `VT1`: not started
3.  `VT2`: not started
4.  `VT3`: not started
5.  `VT4`: not started
6.  `VT5`: not started
7.  `VT6`: not started
8.  `VT7`: not started

## Immediate Next Step

1.  Execute `VT0`.
2.  Freeze the current split inventory around:
    1.  `evaluate_outcome_coupling_unified(...)`
    2.  `build_outcome_coupling_program(...)`
    3.  `OutcomeCouplingOpKind`
    4.  `OutcomeCouplingPairPayload`
    5.  `OutcomeCouplingNWayPayload`
    6.  `OutcomeCouplingGenericNodeIntegralPayload`
    7.  `RankedTransitionCompiler`
    8.  `for_each_sequence_node_transition(...)`
    9.  `sequence_prefix_density_resolved(...)`
    10. `eval_node_recursive_dense(...)`
3.  Do not claim any vectorized-tree progress before that freeze exists.
