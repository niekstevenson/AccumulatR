# C++ Port Plan

## Overview
- Goal: migrate the hottest likelihood paths backed by `super_large_likelihood.R` into C++17 with Rcpp interfaces and RcppParallel workers.
- Strategy: peel the stack in layers—distributions → accumulators → pool fast paths → guard/competitor integrations—while keeping the R boundary thin and maintaining current caching semantics.
- Guardrails: never touch R APIs from worker threads, mirror NA/edge handling from the R originals, and regression-test every slice against the existing code before moving up a layer.

## Immediate Hotspots
- `R/dist.R` (`exgauss_pdf`, `exgauss_cdf`, `dist_registry`) underpin every accumulator evaluation; moving them to C++ unlocks native math for downstream calls.
- Accumulator helpers in `R/likelihood_primitives.R:1` (`.acc_density`, `.acc_survival`, `.acc_cdf_success`, etc.) are thin wrappers over the distributions and a natural second target once the kernels live in C++.
- Pool fast paths (`pool_coeffs` in `R/pool_math.R`, `.pool_density_fast_value` / `.pool_survival_fast_value` in `R/likelihood_primitives.R:115`) are tight loops that benefit immediately from native code and threaded evaluation over members/time points.
- Forced-scenario pool logic (`.build_pool_templates`, `.pool_density`, `.pool_survival` in `R/likelihood_primitives.R:202+`) still allocates per-call structures; a C++ port can precompute templates once and reuse them safely across trials.
- Guard density/integration (`.make_guard_density_fn`, `.make_guard_cdf_fn` at `R/likelihood_kernels.R:831`) and `.scenario_density_with_competitors` (`R/likelihood_kernels.R:2663`) still rely on `stats::integrate`; replacing them with deterministic C++ quadrature removes the dominant latency spike.
- Competitor survival aggregation (`.fast_survival_product`, `.compute_survival_product` in `R/likelihood_kernels.R:142` and `:2733`) runs on every guarded outcome and should leverage the same native accumulators/pool results.

## Phased Migration
1. **Distributions Kernel**
   - Implement lognormal, gamma, and ex-Gaussian PDF/CDF/rng functions in C++17 (using `<random>` and `<cmath>`), replicating all NA checks.
   - Expose them via Rcpp modules and update `dist_registry()` to return lightweight R wrappers over the native calls.
   - Add unit tests comparing R and C++ results for random parameter grids.
2. **Accumulator Layer**
   - Port `.acc_density`, `.acc_density_success`, `.acc_survival`, `.acc_cdf_success` so they call the C++ distributions directly, applying onset and `q` logic in native code.
   - Return plain doubles; keep scenario attribute management on the R side for now.
   - Bench against the R implementations inside `speed_profile.R`.
3. **Pool Fast Paths**
   - Reimplement `pool_coeffs`, `.pool_density_fast_value`, `.pool_survival_fast_value` in C++, threading across pools or time points with RcppParallel.
   - Ensure shared caches (pool member maps) stay on the R side; pass precomputed member indices/parameters into the C++ worker structs.
4. **Forced Templates & Scenarios**
   - ✅ `.build_pool_templates`, `.pool_density`, `.pool_survival` now live in C++, templates cached in the bundle.
   - ✅ Scenario weights/forced IDs returned as structs (converted to R lists on the main thread).
   - **Next extension**: port guarded-expression evaluation (density/survival/scenario) to C++ so R only dispatches calls.

5. **Native Guard Evaluator Roadmap**
   1. Extend the native context to include compiled-node metadata (node kind, child IDs, fast ops) so C++ code can resolve references/blockers/unless expressions without touching R closures.
      - ✅ `native_context_build()` now emits accumulator, pool, and compiled-node tables (IDs, kinds, child links) so native code can look up references/blockers/unless nodes directly.
   2. Build reusable native event/pool evaluators (density/survival/cdf with forced IDs) on top of the accumulator/pool kernels to support guards and other expression types.
      - ✅ `native_node_eval_cpp` now walks AND/OR/NOT nodes recursively using the accumulator/pool kernels, memoizes per-node evaluations, and mirrors forced-ID/component handling (validated in `tests/test_native_node_eval.R`).
   3. Port `guard_density`, `guard_effective_survival`, and guard scenario assembly to C++ (using the Boost Gauss–Kronrod integrator). Wire R guard helpers to call into the native module.
      - ✅ Guard density/effective survival live in C++ (`native_guard_eval_cpp` / `native_guard_effective_survival_cpp`), and scenario lists now come from `native_guard_scenarios_cpp` so R guard helpers only marshal arguments.
   4. After guards, port competitor clustering/survival products to C++ so `.scenario_density_with_competitors()` is fully native.
      - ✅ `native_competitor_survival_cpp` clusters expressions by shared sources and mirrors the R guard-aware chaining logic; `.compute_survival_product()` now delegates to the native path (validated in `tests/test_competitor_cpp.R`).
      - ✅ `native_density_with_competitors_cpp` evaluates compiled outcomes (including guards) and multiplies by the native competitor survival, so `.scenario_density_with_competitors()` now reduces to a thin wrapper when compiled nodes are available.
   5. Move outcome integrals to C++ so R only hands off arguments.
      - ✅ `native_outcome_probability_cpp` integrates the native density (via Gauss–Kronrod) and `.integrate_outcome_probability()` now just tries the native path before falling back to the legacy R integrator (see `tests/test_outcome_probability_cpp.R`).
5. **Guard Integrals**
   - Introduce a C++ guard evaluator that consumes the accumulator/pool interfaces, computes densities, and integrates with a Gauss–Kronrod or adaptive Simpson routine honoring `.integrate_rel_tol()` / `.integrate_abs_tol()`.
   - Replace `.make_guard_*` closures with thin R wrappers around the native functions; thread the setup via compiled node metadata.
6. **Competitor Aggregation**
   - Port `.fast_survival_product` and `.compute_survival_product` so competitor clusters reuse the C++ survival evaluators, enabling parallel evaluation and eliminating redundant R recursion.
   - Keep cluster discovery in R initially; pass clustered node IDs to C++ for execution.
7. **Scenario Aggregation**
   - Mirror `.node_scenarios_at` / `.expr_scenarios_at` in C++ so compiled nodes (event/pool/AND/OR/guard) emit scenario records directly from the native context. The R helpers then become thin wrappers, removing the last heavy R-only recursion.

6. **Prep Serialization & Context Bootstrap**
   - ✅ Defined a `NativePrepProto` struct-of-arrays plus a compact binary serializer so compiled preps can be snapshotted once and reloaded without re-running the R prep (see `native_prep_serialize_cpp` / `native_context_from_proto_cpp`).
   - ✅ `build_native_context()` now routes through the proto builder, so both live and serialized preps produce identical native contexts and share the same cache wiring.
   - ✅ Cache bundles now stash the proto blob alongside the native pointer, letting `.prep_native_context()` rebuild contexts lazily after serialization/resume without re-running `.precompile_likelihood_expressions()`; ad-hoc callers still get the blob via the pointer attribute when no bundle exists.
   - ✅ Proto parameters carry scalar/vector numeric and logical values, and native accumulators now store the resolved distribution coefficients directly, so recursive evaluators no longer rebuild `Rcpp::List` arguments on every density/survival call.
   - ✅ Trial parameter rows are now treated as canonical per-accumulator specs: `.build_trial_overrides()` simply hydrates the supplied rows into accumulator definitions (optional fallback to the legacy defaults behind `options(uuber.require.param.table = FALSE)`), which becomes the default estimation flow both in R and when we move the orchestrator to C++.
   - **Next extension**: allow the serialized blob to be generated outside of R (e.g., from a generator spec dumped to disk) so pure C++ clients can bootstrap contexts without linking against the R prep helpers.
7. **Native Trial Driver**
   - ✅ `native_trial_mixture_cpp` evaluates compiled outcome nodes for every component in a trial using the shared native context, competitor IDs, and the per-trial `TrialOverrides*`, returning the weighted mixture density entirely in C++.
   - ✅ `.likelihood_mixture_likelihood()` now calls the native driver whenever trial rows were materialized via `.likelihood_build_trial_plan()`, so `log_likelihood_from_params()` issues a single `.Call` per trial and falls back to the R loop only when RTs are missing or nodes lack compiled fast paths.
   - **Next extension**: teach the driver about NA/alias outcomes, donor/guess redistribution, and guard-integrated probabilities so R no longer needs to special-case those paths before accumulating log-likelihoods.

## Operational Considerations
- Cache compatibility: continue to populate and clone `likelihood_cache_bundle` (`R/likelihood_common.R:133`) so overrides still get copy-on-write bundles.
- Thread safety: restrict R interaction to the main thread; share data across workers via immutable structures or thread-safe containers.
- Testing: extend `tests/` with regression checks for each migrated layer and compare full likelihood runs to `likelihood_old.R` where feasible.
- Profiling: repeat `speed_profile.R` after each phase to confirm gains and catch regressions early.

## Next Actions
1. Capture current hotspot timings (profvis or `speed_profile.R`) to validate priorities.
2. Scaffold an Rcpp package skeleton (e.g., `src/likelihood_core.cpp`) and wire in the distribution kernel.
3. Write regression tests for distribution/accumulator parity before touching pools.

## Scenario Aggregation Contract

### Shared entry semantics
- `.node_scenarios_at()` always canonicalizes `forced_complete` / `forced_survive` via `.coerce_forced_ids()` before doing anything else, so all downstream helpers receive sorted, de-duplicated integer vectors tagged as `forced_ids` (R/likelihood_kernels.R:1381-1392, R/likelihood_prep.R:78-115).  
- If a compiled node exposes a prebuilt `scenario` function (`node$scenario_fn` or `.node_ops_get(node, "scenario", …)`), that delegate is invoked immediately and the recursive logic is skipped; current examples are guard helper closures that call `native_guard_scenarios_cpp` (R/likelihood_kernels.R:1382-1387, 1848-1879).  
- Scenario results are cached inside the evaluation state entry keyed by `(node_id, component, t, forced_complete, forced_survive)`; misses create an entry via `.node_eval_entry()` and `store()` writes the list back with `.state_entry_set(entry, "scenarios", …)` (R/likelihood_kernels.R:1393-1408).  
- Every record is created through `.make_scenario_record()` which rejects non-positive weights and re-coerces forced sets before returning `list(weight, forced_complete, forced_survive)` (R/likelihood_common.R:103-123). Native code must match both validation and ordering semantics.

### Compiled node behaviors
- **Event nodes** (`expr$kind == "event"`):  
  - Evaluate `.event_density_at()` to get the finisher density; bail out if the returned scalar is non-finite or ≤ 0 (R/likelihood_kernels.R:1409-1419).  
  - If the density carries `attr(weight, "scenarios")` (populated by pool templates), treat each nested record as a child scenario: multiply the event’s forced-complete set by the source ID plus each pool finisher’s forced IDs, union forced-survive sets, and discard any invalid results (R/likelihood_kernels.R:1419-1428).  
  - Without pool scenarios, emit exactly one record with `forced_complete = forced_complete ∪ {source_id}` and `forced_survive` unchanged (R/likelihood_kernels.R:1429-1434).
- **AND nodes**:  
  - Recursively fetch each child’s scenarios, skip empty children, and for every scenario multiply its weight by the conditional CDF of every sibling evaluated at the same `(t, component)` but under the scenario’s forced sets (R/likelihood_kernels.R:1436-1467).  
  - Abort a scenario as soon as any sibling’s CDF is non-finite or ≤ 0, or if the evolving weight becomes non-positive.  
  - Every time we incorporate a sibling we expand `forced_complete` with that sibling’s `sources` so later branches see the combined completion set (R/likelihood_kernels.R:1462-1467).  
  - Final records keep the scenario’s original `forced_survive`.
- **OR nodes**:  
  - Recursively collect child scenarios, then validate each one against the remaining siblings’ `sources` (R/likelihood_kernels.R:1469-1518).  
  - For each sibling, if all required source IDs are already in the scenario’s `forced_complete`, the scenario is invalid because another branch would have already completed; otherwise pick the first “witness” ID that is not yet forced complete and add it to `forced_survive`.  
  - Weight stays untouched (just the child weight), but scenarios fail fast if they cannot produce a witness for some sibling (R/likelihood_kernels.R:1502-1518).
- **Guard nodes**:  
  - Resolve compiled IDs for the reference, blocker, and unless members up front, then pull reference scenarios (or an empty list if the reference is missing) (R/likelihood_kernels.R:1519-1546).  
  - For each reference scenario, call `.guard_effective_survival()` using the scenario’s forced sets; discard any result with `S_eff ≤ 0` (R/likelihood_kernels.R:1547-1566).  
  - Multiply the reference scenario weight by `S_eff` and union blocker source IDs into `forced_survive` (R/likelihood_kernels.R:1560-1571).  
  - `forced_complete` is left exactly as produced by the reference scenario. Guard scope filtering happens inside `.guard_effective_survival()` before integration, so native code has to replicate the same scoping rules.
- **Other compiled kinds**: after exhausting explicit branches, the current code falls back to an empty list. Non-standard nodes are expected to provide their own `scenario` delegate if they can emit scenarios (R/likelihood_kernels.R:1572-1575).

### Guard helper behavior
- `.guard_effective_survival()` first filters forced IDs down to the guard’s scope, caches results by `(guard signature, component, t, forced sets)`, and prefers native guard evaluation when the guard is compiled; otherwise it integrates blocker density × unless-survival on the R side (R/likelihood_kernels.R:1577-1708).  
- `.guard_scenarios_slow()` mirrors the compiled guard path but works over raw expressions, again applying blocker-source unions into `forced_survive` (R/likelihood_kernels.R:1709-1778). Any mismatch here would surface when `.expr_scenarios_at()` handles uncompiled guards.

### Uncompiled expression fallback
- `.expr_scenarios_at()` coerces forced sets, looks for a compiled node, and caches results in `.eval_state_entry()`; compiled hits delegate to `.node_scenarios_at()`, misses drop into `.expr_scenarios_at_slow()` (R/likelihood_kernels.R:2330-2357).  
- `.expr_scenarios_at_slow()` reimplements the same per-kind logic as the compiled path using expression helpers:  
  - **Event**: identical to the compiled version, consuming `.event_density_at()` and pool scenario attributes (R/likelihood_kernels.R:2369-2422).  
  - **AND**: multiplies each child scenario by the unconditional CDF of the remaining expressions (because state caching already captures forced sets) and unions source IDs (R/likelihood_kernels.R:2423-2463).  
  - **OR**: applies the same witness logic as nodes, using `.expr_sources()` to discover required IDs (R/likelihood_kernels.R:2464-2506).  
  - **Guard**: delegates to `.guard_scenarios_slow()` (R/likelihood_kernels.R:2507-2515).  
  - **Fallback**: evaluate the expression density directly and force all sources complete in one record (R/likelihood_kernels.R:2516-2533).  
- Because both compiled and uncompiled paths populate the same cache slot, native code must conserve this contract so mixed graphs (compiled + interpreted nodes) keep sharing memoized scenario lists.

## Full Native Likelihood Roadmap

1. **Native Orchestration Layer**
   - Implement a top-level `native_likelihood_cpp` that mirrors `.likelihood_at()` / `.super_large_likelihood`: iterate trials, components, and time points directly in C++, invoking the native node evaluators (density/survival/scenarios) without returning to R between nodes.
2. **State & Cache Mirroring**
   - Port `.eval_state_entry` / `likelihood_cache_bundle` semantics into C++ (per-node/component/time/forced-set caches, pool template cache, guard quadrature cache) so repeated native calls reuse memoized results exactly as the R path does today.
3. **Prep Serialization / Builder**
   - Either serialize the compiled prep from R into a native-friendly structure or add a native builder that reproduces `.precompile_likelihood_expressions()` so C++ can bootstrap contexts without relying on R helpers.
4. **Trial & Override Plumbing**
   - Recreate R’s per-trial adjustments (shared triggers, forced IDs, competitor overrides) inside the native layer so each run can apply overrides before evaluation. This includes shared-trigger state, guard blockers, and competitor clustering utilities.
5. **Unified API Surface**
   - Expose the native orchestration through a clean C++ function (and corresponding Rcpp wrapper) that accepts the prep context, data, and runtime options. This becomes the single entry point both R and any future direct C++ consumers call, paving the way for a fully native likelihood loop.
