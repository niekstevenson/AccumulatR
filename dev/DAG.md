# DAG Evaluator Contract

This file is the contract for the likelihood evaluator. It is intentionally strict.
The target is not "a faster generic evaluator". The target is a coherent compiled
DAG/IR evaluator where complexity is opt-in, structure is discovered once during
preparation/compilation, and `log_likelihood()` spends its time on math.

## Non-Negotiable Rules

### 1. The Evaluation Boundary Is Trusted

The evaluation boundary starts at `log_likelihood()`.

After that point, inputs are trusted:

- no data-frame reconstruction
- no parameter reconstruction
- no column-name lookup
- no class/type validation
- no shape validation
- no semantic validation
- no "helpful" recovery
- no compatibility fallback

Invalid inputs at this boundary are a caller bug. The evaluator does not diagnose
them. It may return nonsense or crash. Correctness checks belong in model
construction, `prepare_data()`, `build_param_matrix()`, and `make_context()`.

Allowed at evaluation:

- reading already compiled native pointers
- reading already compiled column indices
- reading already canonical parameter slots
- executing compiled math/observation/exact plans
- applying mathematically required domain behavior that is part of the model
  semantics, for example probability zero when an observed response has no mass

Not allowed at evaluation:

- `as.matrix()`, `as.data.frame()`, `Rcpp::as<...>()` for inputs
- `containsElementNamed()`, `colnames()`, `dimnames()`, `names()`, string lookups
- `TYPEOF`, `Rf_isNull`, class checks, length checks, nrow checks as validation
- `std::unordered_map` construction from names or labels
- `.at()` or bounds checks used to protect trusted compiled indices
- fallback to a generic old evaluator when a compiled plan is incomplete

### 2. Evaluation Must Not Discover Structure

The evaluator must not infer what kind of model, observation, outcome mapping,
parameter layout, component layout, trigger layout, or data layout it is handling.
It must read that from precompiled context.

Structure discovery is allowed only in these phases:

- model/spec construction
- preparation of data
- preparation of parameter matrices
- native context compilation
- native trial-layout compilation
- debug-only validation before the evaluation boundary

The correct runtime shape is:

```text
trusted log_likelihood input
  -> native context pointer
  -> prepared layout pointer
  -> canonical parameter matrix pointer
  -> precompiled observation/exact/math programs
  -> numeric result
```

Anything that asks "what columns exist?", "which component does this name mean?",
"does this plan support this observation?", or "which backend should handle this?"
during evaluation is a violation.

### 3. Complexity Is Opt-In

A simple model must not pay for generality it did not request.

This does not mean adding a separate direct backend. It means the compiled plan
for a simple model contains a small canonical DAG, and the compiled plan for a
complex model contains the extra nodes needed for that complexity.

Allowed:

- a compiled opcode such as `LogDensity`, `NoResponseProbability`,
  `WeightedSum`, or `IntegralZeroToCurrent`
- loops over precompiled child spans, source spans, trigger states, or op spans
- dispatch on a precompiled opcode when every model lowers into the same plan
  representation
- a canonical algebraic lowering, if it is represented in the compiled plan and
  not re-discovered during evaluation

Not allowed:

- "try the simple case, else fall back" logic that inspects semantic shape at
  evaluation
- a second semantic engine for easy models
- a generic evaluator that all observations enter by default
- runtime allocation of branch records to ask the evaluator what to do next
- keeping an old path as backup after a compiled path exists

### 4. No Runtime Checking In The Hot Path

The hot path should not defend itself against broken compiled state. Compile-time
validation must prove the state is valid before evaluation.

Forbidden in hot evaluation code:

- null-pointer checks for required workspaces/plans
- invalid-index checks for compiled IDs
- bounds checks for compiled spans
- type checks for R objects
- support checks such as "does this exact plan support the observation?"
- finite/probability checks whose only purpose is validation
- fallback on failed checks

Allowed only if semantically necessary:

- probability zero for impossible events
- lower/upper integration bound checks that define an empty integral
- math-domain transforms that are part of distribution semantics
- debug assertions compiled out of release/hot benchmark builds

Current code still violates this rule in several places. Those violations are
listed below and should be removed, not normalized.

### 5. No Backups

Once a compiled path replaces an old path, the old path must be deleted. Keeping
parallel record-based or generic paths creates two semantic engines and makes
speed regressions likely.

If a model cannot lower into the compiled framework, the fix is a missing lowering
rule. It is not a reason to preserve a runtime fallback.

### 6. Workspaces Are Part Of The Plan Contract

Evaluation must not allocate structure repeatedly.

Allowed:

- workspace reuse sized from compiled plans
- per-thread or per-call numeric buffers when required for thread safety
- fixed-size local arrays for distribution parameters

Not allowed:

- per-trial maps
- per-trial branch records
- per-trial vector construction for structural layout
- per-evaluation name maps
- per-evaluation data/parameter canonicalization
- growing workspace vectors because the plan did not own its layout

### 7. Profiling Guides Every Step

Benchmarks answer whether speed changed. Profiles answer why.

Every nontrivial evaluator change must be guided by profiles for the affected
workload class. This is mandatory for hard cases, because their remaining
slowness is not automatically "math". The profile must show which cost is real
distribution/quadrature math and which cost is evaluator overhead.

Required for performance work:

- a before profile when changing an existing hot path
- an after profile from the same workload and build mode
- the exact profile command
- the workload/model label
- top frames or symbols before and after
- a classification of time spent in math versus evaluator overhead
- the next overhead target identified by the profile

Acceptable profiling targets include:

- source-channel dispatch and cache checks
- compiled math node interpreter dispatch
- quadrature callbacks into exact density roots
- trigger/state enumeration
- ranked frontier/state comparison
- observation request/record execution
- bridge/control construction
- workspace growth/checking

Performance work without a profile is speculation. A change is not complete if
it only improves the simple benchmark while leaving hard-case overhead
unmeasured.

## What Counts As A Real DAG Fix

A change counts as a real DAG fix only if all of these are true:

1. The decision is made before evaluation.
2. The evaluator reads a compiled opcode, span, index, or kernel slot.
3. The simple case has fewer compiled ops, not an evaluator-side escape hatch.
4. The complex case adds explicit compiled ops for the complexity it needs.
5. There is one semantic framework.
6. The old branch/record/fallback code becomes removable.
7. The speed benchmark improves or stays flat.
8. A profile explains the speed result and identifies remaining overhead.
9. Validation and golden likelihood tests pass.

Examples that count:

- lowering no-response into `NoResponseProbability`
- lowering missing RT into `Log(Complement(NoResponseProbability))` when that is
  the canonical observation semantics
- storing component weight parameter indices in the native context
- storing prepared data column indices in `PreparedTrialLayout`
- storing variant leaf row offsets in `ExactVariantPlan`

Examples that do not count:

- wrapping the old branch machine in a new function
- adding an `if (simple)` branch that reimplements semantics beside the plan
- preserving a generic record evaluator as the default path
- scanning names faster
- caching a result of runtime discovery instead of compiling the layout
- checking more carefully before falling back

## Current State

This is the honest map of the current code after the latest speed work.

### Compliant Or Mostly Compliant

`log_likelihood()` no longer reconstructs inputs.

- `R/likelihood_param_interface.R` no longer validates context/data, wraps
  parameters into a list, coerces row-form parameters, canonicalizes columns, or
  checks `nrow()` before calling native likelihood.
- `R/semantic_bridge.R` no longer calls `as.matrix(params)` or
  `as.data.frame(data)`.
- `R/model_definition.R` now emits canonical evaluator column order:
  `q`, `t0`, `p*`, `w`, then weight parameters.

Prepared data layout is now native metadata.

- `PreparedTrialLayout` owns trial/component/label/time column indices.
- `PreparedTrialLayout` owns compressed-trial expansion weights compiled from
  the prepared data `expand` attribute
- C++ evaluation reads columns by compiled index rather than repeated name scan.
- `ok` remains dynamic, but it is now read as a trusted raw logical pointer and
  checked before trial math; no selected mask is rebuilt per call

Parameter component layout is now native metadata.

- sampled-mixture weight parameters are assigned indices during context
  compilation
- evaluation reads them by ordinal slot rather than a name map

Latent component row mapping is now native metadata.

- exact plans own leaf row offsets for their projected component variant
- latent component evaluation uses these offsets instead of a per-trial
  accumulator-name map

Observation probability has a real compiled plan.

- finite density, finite probability, no-response probability, complement,
  weighted sum, and log are explicit `ObservationProbabilityPlan` ops
- `log_likelihood()` evaluates these ops directly in
  `evaluate_observation_plan_direct()`
- finite response, missing-RT response, and missing-all response dispatch now
  lowers to compiled observation state codes that index precompiled log/probability
  plans
- `response_probabilities()` uses the same direct compiled observation-plan
  executor; the old request/record executor has been deleted
- planning-only observation branch vectors are released after context compilation

Exact finite density and finite probability now execute through
`ExactProbabilityProgram`.

- simple terminal top-1 races lower to
  `WeightedTriggerSum(Top1LeafRaceDensity)`
- generic finite density lowers to
  `WeightedTriggerSum(GenericTransitionDensity)`
- finite outcome probability lowers into exact probability programs as
  `Integral(WeightedTriggerSum(...Density))`; the observation executor no longer
  owns the quadrature callback for missing-RT or finite-complement observation
  plans
- terminal no-response lowers to
  `WeightedTriggerSum(TerminalNoResponseProbability)`
- single-child weighted trigger programs own a precompiled `root_child` slot
- programs that do not need shared-trigger variation use a workspace-owned
  default trigger state and do not enumerate trigger-state vectors
- the old `exact_simple_race_target_density()` and
  `exact_no_response_probability()` backup functions have been deleted
- `exact_unranked_target_density()` no longer contains a simple-race
  try-else-generic branch; it executes the compiled exact probability program

### Still Violating The Contract

The C++ `.Call` bridge still has some boundary overhead.

- external pointer access still uses checked pointer wrappers.
- `min_ll` is still passed through the `.Call` boundary as an R object, though
  native code now reads it directly as a trusted scalar

Some of this is small, but it is still not the final contract.

The observation layer still has some transitional shape, but the old native
request/record executor is gone.

- `ComponentObservationPlan` still contains planning fields used while compiling
  observation plans; those fields are now released before runtime evaluation.
- observed-state selection still computes a state code from prepared label/RT
  columns at runtime, then indexes compiled plan tables. That is acceptable
  dispatch, but the remaining defensive checks around invalid components/plans
  should still be removed or made debug-only.

The log-likelihood observation loop still contains defensive checks.

- component code range checks
- plan presence checks
- variant invalid checks
- finite-value checks before accumulation

With trusted prepared data and compiled plans, these should either disappear or
be debug-only assertions.

The exact evaluator still has runtime structure checking.

- many invalid-index and bounds checks remain in `exact_sequence.hpp`
- ranked evaluation still builds dynamic frontiers and compares sequence states
  at runtime, but frontier entries now index workspace-owned sequence-state
  slots instead of owning/copying vector-backed states
- shared-trigger state enumeration is now compiled into
  `ExactCompiledTriggerStateTable`; evaluation reads pre-expanded started masks
  and multiplies q terms. This is no longer a runtime structural expansion.
- the old `exact_plan_supports_observation()` runtime support check has been
  deleted

The compiled math evaluator is still a node-schedule interpreter.

- evaluation dispatches on `CompiledMathNodeKind`
- source-product integral factors are now flattened into compiled factor spans,
  so quadrature samples no longer walk source-value node IDs or call the generic
  source-node interpreter for every factor
- source-product channels now own their source id, effective static source-view
  relation for every valid source, source kernel slot, direct leaf distribution
  kind, required component mask, and bounds; `source_view_id == 0` means an
  explicitly compiled `Unknown` relation, not inherited runtime structure
- source-product integrals now own `CompiledMathSourceProductProgram` objects
  for conditioned sources, exact gates, leaf-absolute kernels, onset
  convolution, and pool k-of-n composition. The residual
  `evaluate_source_product_channel_fill()` / `source_product_base_fill()` loader
  API has been deleted.
- source-product pool programs use compiled scratch offsets instead of allocating
  member/prefix/suffix vectors during evaluation
- cache layouts and condition ids are compiled, but cache-valid/evaluator/time
  checks still occur during evaluation
- some workspace growth/checking remains through `ensure_size()` and time-slot
  guards
- quadrature still calls back into compiled roots instead of always using
  canonical closed forms when the plan can provide them

The old source-channel load path has been deleted.

- ordinary `SourcePdf`, `SourceCdf`, and `SourceSurvival` roots now carry a
  compiled source arithmetic program id
- those ordinary roots execute the same typed source arithmetic program machinery
  used by source-product integrals
- `SourceChannelLoad`, `CompiledSourceChannelPlan`,
  `compiled_math_store_source_channels()`, `compiled_math_source_channels()`,
  source-channel workspace arrays, and
  `CompiledSourceView::source_channels_for_slot()` have been deleted
- the remaining `ExactSourceChannels` object is now a runtime holder for loaded
  leaf inputs, compiled source-view identity, direct-leaf availability, and
  bound-time reads; it no longer owns the old source-channel evaluator/cache path
- source-channel condition bounds are now compiled into
  `CompiledSourceBoundPlan` plus flat `CompiledSourceBoundTerm` spans; evaluation
  no longer walks condition fact arrays to discover exact/lower/upper bounds
- source-view relations are now dense compiled `(source_view_id, source_id)`
  tables; source-channel evaluation no longer walks overlay relation views
- source-channel plans are compiled, but non-source-product compiled math loads
  still execute through source-channel objects, cache checks, and distribution
  calls

The response-probability API now uses the same compiled observation plans, but
it is not a hot-path API.

- it still performs more R-side normalization and parameter-row alignment
- this is acceptable for now because `log_likelihood()` speed is the only hot
  performance target, but response probabilities must not reintroduce a second
  semantic evaluator

Harder cases are not yet proven math-dominated.

- simple cases now avoid most known checking and reconstruction overhead
- harder cases still contain evaluator overhead mixed with real numerical work
- future decisions about these cases must come from profiles, not intuition
- any claim that a hard case is slow because of math must be backed by a profile
  showing that evaluator overhead is no longer the dominant cost

Latest accepted profiles:

- ranked-only profile
  `dev/scripts/scratch_outputs/profile_ranked_state_pool.sample.txt`:
  `example_23_ranked_chain` rose from about 150k loops in the previous ranked
  profile to about 246k loops in the same 6s sample after ranked frontier entries
  stopped owning vector-backed states. Allocator frames dropped from a dominant
  ranked cost to residual setup/reuse noise.
- hard mixed profile
  `dev/scripts/scratch_outputs/profile_source_view_materialized_pruned_hard.sample.txt`:
  the old source-product loader symbols are absent
  (`evaluate_source_product_channel_fill`,
  `evaluate_source_product_conditioned_fill`, `source_product_base_fill`,
  `load_source_product_base_fill`, `load_source_product_conditioned_fill`).
  Source-product source-view inheritance symbols are also absent
  (`source_product_relation`, `parent_source_view_id`).
  The remaining source-product profile is the fused source-product integral
  loop, direct leaf distribution math (`Rf_pnorm_both`, `Rf_plnorm`,
  `Rf_dlnorm`, `log`, `exp`), and small compiled program frames such as
  `compiled_math_source_product_program_conditioned_fill`.

Latest accepted installed benchmark:

- `dev/scripts/scratch_outputs/benchmark_speed_source_arithmetic_cleanup_repeat.csv`
- versus `benchmark_speed_source_product_arithmetic_program_rerun.csv`, the
  clean rerun is effectively flat to slightly slower at benchmark granularity:
  `example_6_dual_path` 1.310 ms -> 1.330 ms,
  `example_16_guard_tie_simple` 4.120 ms -> 4.240 ms,
  `stop_change_shared_trigger` 0.155 ms -> 0.155 ms,
  `stim_selective_stop` 338 ms -> 345 ms.
- simple rows are unchanged at timer granularity. This cleanup is accepted as an
  architectural pruning step, not as a hard-case speed win. The small hard-row
  slowdown is not explained by old source-channel dispatch, which is absent from
  the profile; it is consistent with the remaining source-product/distribution
  wall and benchmark noise at this granularity.

Latest hard mixed profile:

- `dev/scripts/scratch_outputs/profile_source_arithmetic_cleanup_hard.sample.txt`
- workload:
  `example_2_stop_mixture,example_6_dual_path,example_10_exclusion,example_16_guard_tie_simple,example_23_ranked_chain,stop_change_shared_trigger,stim_selective_stop`
- old ordinary source-channel symbols are absent:
  `SourceChannelLoad`, `compiled_math_store_source_channels`,
  `compiled_math_source_channels`, `source_channels_for_slot`,
  `evaluate_source_channel_plan_at`
- top collapsed frames are now:
  `compiled_math_source_product_direct_leaf_scalar` 47 samples,
  `compiled_math_source_product_value_for_ops` 42,
  `evaluate_compiled_source_product_integral_kernel` 42,
  `Rf_pnorm_both` 40, `exp` 36, `log` 33, `Rf_dlnorm` 15,
  `Rf_pnorm5` 15, and
  `CompiledSourceChannels::compiled_source_bounds` 9.
- conclusion: the requested old source-channel cleanup is real; the remaining
  hard-case wall is source-product integral arithmetic, distribution kernels,
  and small compiled bound/source-program frames.

## Required Pruning

These are not optional cleanups. They are required to keep one coherent framework.

1. Delete row-form/log-likelihood compatibility paths permanently.
2. Delete runtime support checks such as exact observation support once
   preparation/compilation guarantees valid observations.
3. Delete fallback logic that returns `min_ll` for malformed compiled state.
   Malformed compiled state must fail before evaluation.
4. Delete source/name maps from any path reachable after `log_likelihood()`.
5. Delete per-trial structural vectors and replace them with workspace-owned
   buffers or compiled spans.
6. Delete unused enums/structs as soon as their old paths are removed.
7. Delete "old direct backend" concepts if they reappear. Canonical algebra must
   lower into the one exact framework.

## Required Development Work

### A. Keep Observation Dispatch Fully Compiled

Current status: completed for the first-rank observation probability layer.
Finite, missing-RT, and missing-all observations now map to compiled state codes
that index precompiled observation probability plans.

Target:

```text
PreparedObservationState
  -> component code
  -> observed label code / missing state code
  -> precompiled ObservationProbabilityPlan pointer/index
```

No plan-null checks. No branch discovery. No request records. No fallback.

Remaining work under this heading is cleanup, not a new semantic path:

- remove or debug-assert the remaining trusted-state defensive checks
- keep planning branch vectors out of the runtime context
- extend the same state-code discipline if ranked/partial observation dispatch
  grows beyond first-rank observations

### B. Keep Response Probabilities On Compiled Plans

Current status: completed for native evaluation. `response_probabilities()` now
reuses the same `ObservationProbabilityPlan` op executor as `log_likelihood()`.

Maintained target:

- build query rows in prepared/canonical form before evaluation
- execute direct plan ops
- do not reintroduce `ObservationPlanRequest`
- do not reintroduce `ObservedRecord`
- do not reintroduce record-based exact query batches

### C. Fold Exact Probability Into A Unified Exact Probability Program

Current status: active for finite density, finite probability mass, and
canonical terminal no-response. `ExactVariantPlan` owns an
`ExactProbabilityProgramSet`, and `log_likelihood()` finite-density and
finite-probability evaluation enter that program rather than branching between
simple and generic exact implementations.

Target exact probability op examples:

- `Top1LeafRaceDensity`
- `TerminalNoResponseProbability`
- `GenericTransitionDensity`
- `GenericTransitionProbability`
- `RankedTransitionSequence`
- `Integral`
- `WeightedTriggerSum`

Every exact observation lowers into a sequence of these ops. The evaluator executes
the program. It does not try one shape and fall back to another.

Implemented now:

- `Top1LeafRaceDensity`
- `TerminalNoResponseProbability`
- `GenericTransitionDensity`
- `GenericTransitionProbability` executor support
- `Integral`
- `WeightedTriggerSum`

Lowered exact probability program shapes now include:

- finite density: `WeightedTriggerSum(Top1LeafRaceDensity)` or
  `WeightedTriggerSum(GenericTransitionDensity)`
- finite probability: `Integral(WeightedTriggerSum(Top1LeafRaceDensity))` or
  `Integral(WeightedTriggerSum(GenericTransitionDensity))`
- terminal no-response:
  `WeightedTriggerSum(TerminalNoResponseProbability)`

The finite probability shape intentionally keeps `Integral` outside
`WeightedTriggerSum` for shared-trigger cases. Hoisting trigger enumeration
outside the integral changed shared-trigger finite-complement numerics; the
current shape preserves the existing mixture-density quadrature semantics while
still making the exact program own the integration.

Remaining work:

- `RankedTransitionSequence`
- non-terminal no-response/missing-response semantics beyond the terminal
  survival-product canonical case
- canonical finite-probability forms that can safely avoid quadrature; these
  must be represented in the exact probability program and proven by validation,
  not rediscovered during evaluation

### D. Compile Trigger And State Enumeration

Current status: completed for shared-trigger table expansion; improved but not
fully eliminated for ranked frontier state handling.

Shared-trigger variants now compile `ExactCompiledTriggerStateTable` into the
exact variant plan. Each state owns:

- a fixed probability multiplier
- flat q-weight terms for dynamic trigger probabilities
- a pre-expanded `shared_started` mask

Evaluation reads these states directly. The old runtime
`enumerate_trigger_states_into()` expansion has been deleted.

Target:

- trigger states are compiled into arrays of q slot indices, fixed probabilities,
  and started masks
- shared-trigger cases are pre-expanded into a plan-specific trigger-state table
- evaluation reads the table and multiplies terms
- no runtime allocation or reconstruction of trigger state buffers

Remaining work under this heading:

- ranked sequence state/frontier handling is still dynamic, but it now reuses
  workspace-owned state pools and moves frontier entries by state index. The
  remaining dynamic work is frontier merging/state comparison, not per-entry
  vector allocation.
- trigger-state probability is still computed per evaluation from q values, which
  is expected; the structure of the state table is no longer discovered there

### D2. Compile Source-Channel Condition And Source-View Structure

Current status: completed for source bounds, static source-view relation lookup,
condition-cache layouts, guard/expr upper-bound term lowering, source-product
channel component requirements, and source-product arithmetic lowering.

Implemented:

- `CompiledSourceBoundPlan` now points at flat `CompiledSourceBoundTerm` spans
  containing time-slot ids
- source exact/lower/upper bound structure is built during context compilation
- source-channel evaluation no longer calls the old `resolve_source_bounds()`
  or `condition_has_fact()` machinery
- source-view relations are compiled into dense relation tables, and the old
  runtime `RelationView` overlay chain has been deleted
- condition source relations used by source ordering are compiled into dense
  condition/source relation tables
- source-channel condition-cache ids are compiled as static ids or dynamic
  time-dependency spans; evaluation no longer asks the workspace to discover
  condition time dependencies
- compiled math guard/expr upper-bound nodes now read flat timed upper-bound
  term spans; evaluation no longer calls `compiled_math_condition_by_id()` or
  walks guard/expr condition fact arrays
- source-product integral kernels now compile unique source-channel slots and
  the required component mask (`pdf`, `cdf`, `survival`) for each slot
- direct leaf-absolute source-channel plans can evaluate only the components
  requested by a source-product channel; survival-only competitors no longer
  compute unused densities, and pdf-only factors no longer compute unused
  CDF/survival terms
- source-product factors now reference a compiled
  `CompiledMathSourceProductProgram`; conditioned lower/between/exact handling,
  base exact gates, onset convolution, pool k-of-n composition, and direct leaf
  kernels are represented as typed source-product program kinds
- pool source-product programs own member program spans and compiled scratch
  offsets; evaluation reuses workspace scratch instead of allocating DP vectors
- the old source-product scalar loader/cache API has been deleted:
  `evaluate_source_product_channel_fill()`,
  `evaluate_source_product_conditioned_fill()`,
  `source_product_base_fill()`, `load_source_product_base_fill()`, and
  `load_source_product_conditioned_fill()` no longer exist
- source-product relation-sensitive programs now carry a compiled final
  source-view relation. `source_product_relation()` and the evaluator's
  `parent_source_view_id` inheritance path have been deleted, and compile-time
  validation rejects relation-sensitive source-product programs without a
  materialized relation.
- the unused source-product channel cache arrays and store helper were pruned;
  source-product sample caching now uses only compiled program ids.

Remaining work:

- source-view evaluator objects still exist as cache/evaluator identities, though
  they now wrap a compiled source-view id rather than reconstructing overlays
- ordinary source roots still dispatch through the compiled math node interpreter
  before entering the source arithmetic program
- the source arithmetic program still has small program/cache/bounds frames; hard
  profiles show these are residual compared with source-product integral work and
  distribution kernels, but they are not zero

### E. Flatten Compiled Math Execution

The current compiled math schedule is still an interpreter.

Target:

- root schedules own compact op arrays
- op inputs are offsets/spans, not semantic nodes
- source-channel loads are fused where possible
- cache use is opt-in per compiled op, not checked globally
- time slots are known and pre-sized
- workspace is fully sized by context/plan

### F. Make Checks Debug-Only Or Move Them Before Evaluation

Runtime checks that validate compiled state should become:

- compile-time validation
- preparation-time validation
- debug assertions compiled out of benchmark/release paths

No hot-path code should branch to recover from impossible compiled state.

### G. Keep Evaluation Control Trusted

Current status: completed for `ok` and `expand`.

- compressed-trial expansion weights are compiled into `PreparedTrialLayout`
- `ok` is dynamic and is read as a trusted raw logical pointer
- false `ok` trials keep their `min_ll` contribution but are skipped before
  exact/observation math

Maintained target:

- do not rebuild selected masks per call
- do not rebuild expand weights per call
- do not reinterpret expansion structure after layout preparation

## Forbidden Patterns For Future Changes

Do not add code matching these patterns after the evaluation boundary:

```cpp
Rcpp::DataFrame(...)
Rcpp::as<...>(...)
containsElementNamed(...)
Rf_isNull(...)
TYPEOF(...)
std::unordered_map<std::string, ...>
colnames / dimnames / names
if (plan is simple) do direct semantic calculation else generic semantic calculation
records.push_back(...)
requests.push_back(...)
try compiled; if not, use old path
```

Exceptions must be explicitly justified in this file before they are added. The
default answer is no.

## Success Gates

Correctness gates:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE, helpers=FALSE); testthat::test_local(reporter = "summary")'
Rscript dev/validation/run_validation.R
```

Speed gate:

```sh
R CMD INSTALL --preclean .
ACCUMULATR_BENCH_INSTALLED=true \
ACCUMULATR_BENCH_TRIALS=50 \
ACCUMULATR_BENCH_N_REP=5 \
ACCUMULATR_BENCH_TARGET_SEC=0.10 \
ACCUMULATR_BENCH_MAX_INNER_REPS=2000 \
Rscript dev/scripts/benchmark_speed.R
```

The default benchmark script still supports `pkgload::load_all()` for development
smoke runs. Use `ACCUMULATR_BENCH_INSTALLED=true` for performance claims; debug
`load_all()` timings are not comparable to optimized installed-package baselines.

Focused regression gate for the cases that exposed the failure:

```sh
Rscript - <<'RS'
suppressPackageStartupMessages(library(pkgload))
load_all('.', quiet = TRUE, helpers = FALSE)
source(file.path('dev', 'examples', 'new_API.R'))
labels <- c('example_1_simple', 'example_5_timeout_guess',
            'example_21_simple_q', 'example_22_shared_q')
make_case <- function(label) {
  structure <- new_api_examples[[label]]
  pars <- new_api_example_params[[label]]
  params_df <- build_param_matrix(structure, pars, n_trials = 50)
  data_df <- simulate(structure, params_df, seed = 123, keep_component = TRUE)
  prepared <- prepare_data(structure, data_df)
  ctx <- make_context(structure)
  params_slim <- build_param_matrix(structure, pars, trial_df = prepared)
  invisible(log_likelihood(ctx, prepared, params_slim))
  list(ctx = ctx, prepared = prepared, params = params_slim)
}
bench <- function(label, inner = 5000L, reps = 5L) {
  case <- make_case(label)
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc(FALSE)
    times[[i]] <- system.time({
      for (j in seq_len(inner)) {
        invisible(log_likelihood(case$ctx, case$prepared, case$params))
      }
    })[['elapsed']] * 1000 / inner
  }
  data.frame(label = label,
             inner_reps = inner,
             median_per_eval_ms = median(times),
             min_per_eval_ms = min(times),
             max_per_eval_ms = max(times))
}
print(do.call(rbind, lapply(labels, bench)), row.names = FALSE)
RS
```

Current focused benchmark reference from the latest run:

```text
example_1_simple          0.0094 ms/eval
example_5_timeout_guess   0.0210 ms/eval
example_21_simple_q       0.0187 ms/eval
example_22_shared_q       0.0125 ms/eval
```

Any change that worsens these without moving cost into real mathematical work is
a regression.

Profiling gate:

```sh
ACCUMULATR_PROFILE_R_SCRIPT=dev/scripts/profile_workload_nested.R \
ACCUMULATR_PROFILE_FILE=dev/accumulatr_profile_cpp.txt \
ACCUMULATR_PROFILE_DURATION=20 \
ACCUMULATR_PROFILE_INTERVAL=1 \
bash dev/scripts/profile_cpp_simple.sh
```

Mixed profiling gate for exact evaluator work:

```sh
ACCUMULATR_PROFILE_TARGET_SIGNALS=1 \
ACCUMULATR_PROFILE_R_SCRIPT=dev/scripts/profile_workload_mixed.R \
ACCUMULATR_PROFILE_FILE=dev/scripts/scratch_outputs/profile_mixed.sample.txt \
ACCUMULATR_PROFILE_DURATION=20 \
ACCUMULATR_PROFILE_WORKLOAD_SECONDS=20 \
ACCUMULATR_PROFILE_INTERVAL=1 \
bash dev/scripts/profile_cpp_simple.sh
```

The mixed workload defaults to a balanced set of simple race, q/shared-q,
observation-wrapper, mixture, exclusion, ranked, and shared-trigger cases. Use
`ACCUMULATR_PROFILE_CASES` for a comma-separated override, especially when a
single pathological case such as `stim_selective_stop` should be isolated rather
than allowed to dominate a mixed profile.

Latest exact probability program profile reference:

```text
focused simple/q cases:
  before exact probability program:       104,962 loops/case
  recursive first implementation:          96,738 loops/case
  final flattened/trusted program:        119,116 loops/case

mixed exact workload:
  before exact probability program:        10,036 loops/case
  recursive first implementation:           9,550 loops/case
  final flattened/trusted program:         11,122 loops/case
```

Latest hard-case C profile reference:

```text
stim_selective_stop:
  before C finite-probability lowering: 12 loops / 10s profile
  after C finite-probability lowering:  10 loops / 8s profile
  dominant stack after: ObservationPlan -> ExactProbabilityProgram ->
    GenericTransitionDensity -> compiled source-product integral
  top self frames after: Rf_pnorm_both, standard_leaf_channels,
    CompiledSourceChannels::evaluate_direct_leaf_absolute_at,
    compiled_math_store_source_channels,
    evaluate_compiled_source_product_integral_sample

example_6_dual_path:
  before single-rank dispatch cleanup: 3303 loops / 10s profile
  after single-rank dispatch cleanup:  1846 loops / 6s profile
  exact_ranked_loglik_for_trial no longer appears in the single-rank profile;
  remaining cost is compiled source-product integral/source-channel math

example_16_guard_tie_simple:
  after C profile: 638 loops / 6s profile
  top self frames: Rf_pnorm_both,
    evaluate_compiled_source_product_integral_sample,
    compiled_math_store_source_channels,
    CompiledSourceChannels::evaluate_direct_leaf_absolute_at

stop_change_shared_trigger:
  after C profile: 13711 loops / 6s profile
  top self frames: Rf_pnorm_both,
    CompiledSourceChannels::evaluate_direct_leaf_absolute_at,
    compiled_math_store_source_channels,
    standard_leaf_channels
```

Conclusion from these profiles: for these hard cases, section C moved ownership
of finite probability semantics into the exact probability program, but it did
not make the hard cases materially faster. The remaining dominant cost is not
observation branching or ranked wrapper dispatch. It is section E/source-channel
execution: compiled source-product integrals, source-channel object/cache
dispatch, and distribution kernel calls.

Latest section-E source-channel cleanup profile reference:

```text
Change:
  direct leaf-absolute source-channel candidacy is compiled into
  CompiledSourceChannelPlan; direct evaluation no longer rediscoveries kernel
  shape or revalidates the trusted kernel slot during evaluation; compiled math
  root execution also skips redundant workspace size checks.

Benchmarked against the post-C installed build:
  example_6_dual_path:          2.9254 ms -> 2.8806 ms  (-1.5%)
  example_16_guard_tie_simple:  8.0435 ms -> 7.8462 ms  (-2.5%)
  stop_change_shared_trigger:   0.3850 ms -> 0.3775 ms  (-1.9%)
  stim_selective_stop:        826.0000 ms -> 811.0000 ms (-1.8%)

8s isolated sample loop counts:
  example_6_dual_path:          2600 -> 2643 loops
  example_16_guard_tie_simple:   941 ->  968 loops
  stop_change_shared_trigger:  19035 -> 19891 loops
  stim_selective_stop:             9 ->    10 loops

Conclusion:
  this is a real section-E cleanup, not a semantic fast path, but the payoff is
  small. The dominant hard-case stack remains source-product integral execution,
  source-channel cache/object dispatch, and distribution kernels. The next
  material target is to flatten/fuse the source-product integral program itself.
```

Latest source-product integral flattening profile reference:

```text
Change:
  Source-product integral kernels now own flat CompiledMathSourceValueFactor
  spans. Quadrature samples evaluate those factors directly with a fused
  source-channel cache and trusted time-slot access. The old source-product
  integral loop no longer walks source-value node IDs or calls
  compiled_math_source_value_for_node()/compiled_math_store_source_channels()
  for every factor. Those old source-product helper calls were deleted.

Benchmarked against the post-C installed build:
  example_6_dual_path:          2.9254 ms -> 2.7761 ms  (-5.1%)
  example_16_guard_tie_simple:  8.0435 ms -> 7.2414 ms (-10.0%)
  stop_change_shared_trigger:   0.3850 ms -> 0.3375 ms (-12.3%)
  stim_selective_stop:        826.0000 ms -> 774.0000 ms (-6.3%)

8s isolated sample loop counts:
  example_6_dual_path:          2600 ->  2707 loops
  example_16_guard_tie_simple:   941 ->  1013 loops
  stop_change_shared_trigger:  19035 -> 22119 loops
  stim_selective_stop:             9 ->    10 loops

Conclusion:
  this is a real section-E improvement. It removes an interpreter layer from
  the hot quadrature loop and produces measurable hard-case gains. It is still
  not the final payoff: profiles remain dominated by distribution kernels
  (`Rf_pnorm_both`, `Rf_dlnorm`, `exp`, `log`) plus fused source-factor/channel
  cache execution. The next material target is not more observation wrapping; it
  is further specializing source-channel factor programs and identifying closed
  forms or lower-order quadrature where the compiled plan proves they are valid.
```

Latest source-product component-pruning profile reference:

```text
Change:
  Source-product integral kernels now compile unique source-product channels.
  Each channel owns a required component mask, so direct leaf-absolute channel
  evaluation computes only pdf, cdf, and/or survival components that downstream
  source factors actually read. Full-channel callers keep the straight original
  leaf arithmetic path.

Benchmark against benchmark_speed_final_dag_consolidation.csv:
  example_1_simple:             0.0090 ms ->   0.0090 ms (flat)
  example_2_stop_mixture:       0.1150 ms ->   0.1000 ms (-13.0%)
  example_3_stop_na:            0.1150 ms ->   0.1125 ms ( -2.2%)
  example_5_timeout_guess:      0.0210 ms ->   0.0210 ms (flat)
  example_6_dual_path:          2.6700 ms ->   2.1300 ms (-20.2%)
  example_7_mixture:            0.0130 ms ->   0.0120 ms ( -7.7%)
  example_10_exclusion:         0.0710 ms ->   0.0575 ms (-19.0%)
  example_16_guard_tie_simple:  6.6897 ms ->   6.4483 ms ( -3.6%)
  example_21_simple_q:          0.0160 ms ->   0.0160 ms (flat)
  example_22_shared_q:          0.0140 ms ->   0.0130 ms ( -7.1%)
  example_23_ranked_chain:      0.0220 ms ->   0.0199 ms ( -9.5%)
  stop_change_shared_trigger:   0.3184 ms ->   0.2584 ms (-18.8%)
  stim_selective_stop:        748.0000 ms -> 558.0000 ms (-25.4%)
  geomean: -10.2%

Hard mixed profile:
  profile_hard_current_state_pool.sample.txt ->
  profile_hard_source_product_channel_mask_full.sample.txt

  The workload completed more hard-case loops in the same 8s sample. Despite
  that, Rf_pnorm_both top samples fell from 932 to 643, Rf_dlnorm from 187 to
  133, exp from 624 to 505, and log from 636 to 556. The old
  compiled_math_store_source_factor_channels frame is gone. The new top
  evaluator frame is standard_leaf_channels_mask, which is expected: it is the
  component-pruned leaf math requested by the compiled source-product channel.

Conclusion:
  this is the first source-channel cleanup in this series that materially
  reduces hard-case math work rather than only moving cache symbols. Remaining
  work is still section E: the source-product integral loop itself is still a
  compact interpreter, source-view evaluator objects still exist, and generic
  conditioned/onset/pool source channels need to become compiled component
  programs rather than full EventChannels machinery.
```

Latest source-channel component-program reference:

```text
Change:
  Required component masks now live on ordinary CompiledSourceChannelPlan slots,
  not only on source-product integral channels. SourceChannelLoad evaluation
  requests the slot's compiled mask. Timed source caches track valid components
  per source/time/context, so a later request only computes missing lanes.

  Leaf-absolute, conditioned, onset-convolution, and pool k-of-n source kernels
  now propagate component masks internally:
    - leaf-absolute calls standard_leaf_channels_mask()
    - lower-bound conditioning requests survival or cdf only as required
    - finite upper conditioning requests cdf for bound mass and only requested
      output lanes for the current time
    - onset convolution integrates source pdf and only requested shifted lanes
    - pool k-of-n skips member pdfs and output pdf accumulation unless pdf is
      requested

  Dead source_channel_use_counts bookkeeping was deleted. It was no longer read
  by the compiled evaluator and did not belong in the trusted evaluation path.

Benchmark against benchmark_speed_source_product_channel_mask_full.csv:
  Two same-config runs were taken after the full-mask cache fix.
  Geomean vs previous: -0.4% on run 1, +3.8% on run 2, +1.6% using the
  geometric mean of the two current runs. This is not a reliable speed win over
  the previous source-product component-pruning state.

  Averaged hard rows vs previous:
    example_6_dual_path:          2.1300 ms ->   2.1866 ms (+2.7%)
    example_16_guard_tie_simple:  6.4483 ms ->   6.5688 ms (+1.9%)
    stop_change_shared_trigger:   0.2584 ms ->   0.2606 ms (+0.8%)
    stim_selective_stop:        558.0000 ms -> 566.9400 ms (+1.6%)

  Averaged table vs benchmark_speed_final_dag_consolidation.csv remains faster:
    geomean -8.7%, stim_selective_stop -24.2%, stop_change_shared_trigger
    -18.1%, example_6_dual_path -18.1%.

Hard mixed profile:
  profile_hard_source_product_channel_mask_full.sample.txt ->
  profile_hard_source_channel_component_programs_cachefull.sample.txt

  The profile moved in the intended structural direction but did not expose a
  new hard-case wall-clock payoff. Source-plan samples fell, standard leaf mask
  samples fell slightly, source-product sample cost stayed essentially flat,
  and distribution math remained dominant. This confirms that the change removes
  remaining full-channel work, but the current hard rows are already dominated by
  quadrature plus pnorm/dlnorm/exp/log work after source-product pruning.

Rejected attempted shortcut:
  Precompiling source-view relation directly onto CompiledSourceChannelPlan was
  invalid. A slot with source_view_id = 0 can be evaluated inside a nonzero
  parent source-view context, so the effective relation is not determined by the
  slot alone. Validation caught this in the overlapping-guard test. Do not add a
  fallback around that; the correct future fix is to make effective source-view
  context part of the compiled program shape where it is truly static.

Conclusion:
  This is real DAG progress, not a second semantic engine: source-channel
  component demand is compiled and executed by one path. It is not the next speed
  breakthrough. The remaining section-E overhead is source-view context handling
  and the source-product quadrature program itself; most visible hard-case cost
  is now mathematical distribution evaluation under quadrature.
```

Latest source-product typed arithmetic program reference:

```text
Change:
  Source-product integral channels now carry executable typed fields copied from
  the compiled source-channel plan:
    - source id and source kernel slot
    - direct leaf index, distribution kind, parameter count, and absolute onset
    - source-condition bound plan and condition cache ids
    - required pdf/cdf/survival component mask
    - static source-view relation when the source-product channel has an
      explicit source_view_id

  The hot quadrature loop calls a source-product channel evaluator that reads
  these fields directly. Direct leaf-absolute source-product channels call
  distribution-specific scalar/fill functions instead of
  CompiledSourceChannels::evaluate_source_channel_plan_at().

  Deleted from the source-product hot path:
    - compiled_math_source_product_channel_channels()
    - evaluate_compiled_source_product_integral_sample()
    - source-value node walking inside source-product factors

  A scalar direct-factor attempt was rejected and removed. It looked more
  specialized, but benchmarked worse for example_6, example_16, and
  stop_change_shared_trigger because it duplicated component work and added
  another evaluator branch. It was not the clean path.

Benchmark against benchmark_speed_dag_audit_current.csv:
  example_6_dual_path:          2.1333 ms ->   1.2000 ms (-43.8%)
  example_16_guard_tie_simple:  6.8000 ms ->   4.5000 ms (-33.8%)
  stop_change_shared_trigger:   0.2667 ms ->   0.1800 ms (-32.5%)
  stim_selective_stop:        569.0000 ms -> 398.0000 ms (-30.1%)

Final profile artifacts:
  profile_hard_example_6_dual_path_source_product_static_view_final.sample.txt
  profile_hard_example_16_guard_tie_simple_source_product_static_view_final.sample.txt
  profile_hard_stop_change_shared_trigger_source_product_static_view_final.sample.txt
  profile_hard_stim_selective_stop_source_product_static_view_final.sample.txt

Top-frame counts, before -> final:
  compiled_math_source_product_channel_channels:
    example_6 740 -> 0; example_16 665 -> 0; stop_change 829 -> 0;
    stim 546 -> 0
  evaluate_compiled_source_product_integral_sample:
    example_6 1176 -> 0; example_16 598 -> 0; stop_change 1121 -> 0;
    stim 387 -> 0
  evaluate_source_channel_plan_at:
    example_6 206 -> 0; example_16 319 -> 0; stop_change 243 -> 15;
    stim 320 -> 0
  standard_leaf_channels_mask:
    example_6 451 -> 5; example_16 570 -> 0; stop_change 492 -> 29;
    stim 927 -> 68

Honest conclusion:
  This achieves the requested source-product flattening target for the
  source-product integral path. The source-product integral now owns
  source-channel loads, final source-view relation data, component requirements,
  source bounds, typed base source programs, and leaf distribution kernel choice
  up front. The old generic source-product lookup/sample frames disappear from
  hard profiles.

  This is still not pure math. The remaining top evaluator frame is
  evaluate_compiled_source_product_integral_kernel(), because the fused
  quadrature/product arithmetic now lives there. Generic source-product
  conditioned/base/onset/pool cases lower into
  CompiledMathSourceProductProgram objects rather than EventChannels or residual
  source-product loader functions. The remaining non-math source-product frames
  are small compiled-program interpreter frames and cache fills. The
  mathematical wall is now clear: Rf_pnorm_both, Rf_dlnorm, Rf_pnorm5, exp, and
  log are top frames in the hard profiles.

Latest source-product arithmetic-program cleanup:
  Source-product generic channels now lower into typed
  CompiledMathSourceProductProgram objects. The old scalar residual loader API
  in src/eval/exact_source_channels.hpp has been deleted; conditioned bounds,
  exact gates, onset convolution, pool k-of-n composition, and direct leaf
  kernels are evaluated from compiled program ids and workspace scratch.

  Benchmark against benchmark_speed_source_product_scalar_generic_final.csv:
    example_1_simple:             0.0100 ms ->   0.0100 ms
    example_2_stop_mixture:       0.0650 ms ->   0.0700 ms
    example_3_stop_na:            0.1100 ms ->   0.1150 ms
    example_5_timeout_guess:      0.0200 ms ->   0.0200 ms
    example_6_dual_path:          1.2700 ms ->   1.3300 ms
    example_16_guard_tie_simple:  4.4600 ms ->   4.1400 ms
    stop_change_shared_trigger:   0.1500 ms ->   0.1550 ms
    stim_selective_stop:        347.0000 ms -> 343.0000 ms

  Source-view materialization rerun against
  benchmark_speed_source_product_arithmetic_program_rerun.csv:
    example_1_simple:             0.0100 ms ->   0.0100 ms
    example_2_stop_mixture:       0.0700 ms ->   0.0650 ms
    example_3_stop_na:            0.1150 ms ->   0.1150 ms
    example_5_timeout_guess:      0.0200 ms ->   0.0200 ms
    example_6_dual_path:          1.3300 ms ->   1.3100 ms
    example_16_guard_tie_simple:  4.1400 ms ->   4.1200 ms
    stop_change_shared_trigger:   0.1550 ms ->   0.1550 ms
    stim_selective_stop:        343.0000 ms -> 338.0000 ms

  Profile artifact:
    profile_source_view_materialized_pruned_hard.sample.txt

  Profile check:
    evaluate_source_product_channel_fill,
    evaluate_source_product_conditioned_fill, source_product_base_fill,
    load_source_product_base_fill, and load_source_product_conditioned_fill are
    absent because those functions no longer exist. source_product_relation and
    parent_source_view_id are absent too. Remaining source-product time is the
    fused integral kernel, direct leaf distribution math, and small
    compiled-program frames such as
    compiled_math_source_product_program_conditioned_fill.
```

The profile workload must match the change. Use simple workloads for simple
regressions and hard workloads for hard-case work. The profile report must name
the top evaluator-overhead frames and separate them from genuine mathematical
cost such as distribution kernels and quadrature.

## Review Checklist

Before merging evaluator work, answer these questions in the PR or commit note:

1. What structure is newly compiled?
2. What runtime discovery was removed?
3. What checks moved out of evaluation?
4. What old path was deleted?
5. Which benchmark rows changed?
6. Which profile command and artifact were used?
7. What were the top overhead frames before and after?
8. What cost is now real math, and what evaluator overhead remains?
9. Which validation/golden tests prove semantic equivalence?
10. Is there now one semantic path, or did this introduce another engine?

If the answer to question 4 is "none", the change is suspect. If the answer to
question 10 is "another engine", reject it.
