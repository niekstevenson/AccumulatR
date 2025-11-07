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
5. **Guard Integrals**
   - Introduce a C++ guard evaluator that consumes the accumulator/pool interfaces, computes densities, and integrates with a Gauss–Kronrod or adaptive Simpson routine honoring `.integrate_rel_tol()` / `.integrate_abs_tol()`.
   - Replace `.make_guard_*` closures with thin R wrappers around the native functions; thread the setup via compiled node metadata.
6. **Competitor Aggregation**
   - Port `.fast_survival_product` and `.compute_survival_product` so competitor clusters reuse the C++ survival evaluators, enabling parallel evaluation and eliminating redundant R recursion.
   - Keep cluster discovery in R initially; pass clustered node IDs to C++ for execution.

## Operational Considerations
- Cache compatibility: continue to populate and clone `likelihood_cache_bundle` (`R/likelihood_common.R:133`) so overrides still get copy-on-write bundles.
- Thread safety: restrict R interaction to the main thread; share data across workers via immutable structures or thread-safe containers.
- Testing: extend `tests/` with regression checks for each migrated layer and compare full likelihood runs to `likelihood_old.R` where feasible.
- Profiling: repeat `speed_profile.R` after each phase to confirm gains and catch regressions early.

## Next Actions
1. Capture current hotspot timings (profvis or `speed_profile.R`) to validate priorities.
2. Scaffold an Rcpp package skeleton (e.g., `src/likelihood_core.cpp`) and wire in the distribution kernel.
3. Write regression tests for distribution/accumulator parity before touching pools.
