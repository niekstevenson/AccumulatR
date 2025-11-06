# Likelihood Modernization Aims

## 1. Match or exceed `likelihood_old.R` performance

- **Current status**: The modular pipeline (compiled nodes, pool fast-paths, cached state) brings the new engine within noise of the legacy timings in `speed_profile.R`. Microbenchmarks show parity for guard-heavy examples, but guard integrals and forced pool templates remain the dominant hotspots.
- **Still to do**: Persist pre-built pool templates and guard survival integrals so repeated trials avoid recomputation. Extend fast-path coverage to forced pool scenarios and guarded mixed expressions. Automate profiling comparisons so regressions are caught in CI rather than ad hoc.

## 2. Preserve a universal, general-purpose framework

- **Current status**: Expression compilation handles events, pools, AND/OR compositions, and guards with protectors. Examples in `new_api_examples` cover nested pools, shared triggers, and inhibitors, demonstrating broad model support.
- **Still to do**: Finish evaluation support for `not` expressions (`!expr`), which are parsed but currently drop probability mass. Add regression fixtures that mix guards, negations, and k-of-n pools so future changes keep the framework truly universal.

## 3. Favour principled optimizations over ad-hoc patches

- **Current status**: Fast-path logic is encapsulated in reusable helpers (e.g. pool density/survival, guard effective survival) rather than duplicated branches. Competitor handling uses structural clustering instead of outcome-specific checks.
- **Still to do**: Unify duplicated guard logic between fast and full evaluators, and ensure new shortcuts (e.g. forced pool fast paths) are derived algebraically. Document invariants for each optimization to keep future tweaks grounded in first principles.

## 4. Restrict caching to within-trial reuse

- **Current status**: The evaluation cache created by `.eval_state_create()` scopes memoization to a single call stack, ensuring densities, survivals, and scenarios are reused only within the current trial. The old API also demonstrates that a small, deterministic shared structure (unique `(component|outcome|rt)` map) is enough to avoid recomputing NA-guard integrals across identical trials.
- **Still to do**: Replace ad hoc environments with a centralized, model-scoped cache that is built once during prep and handed to every evaluator. It should hold cross-trial invariants—precompiled node walks, canonical NA integrals, reusable guard plans—without becoming a mutable global. Phase out the option-driven likelihood cache in `R/likelihood_param_interface.R` in favor of this explicit data structure so the new table pathway matches the old API’s efficiency while keeping cache scope auditable.
- **Implementation plan**:
  1. **Design cache bundle**: Introduce an S3 `likelihood_cache_bundle` that stores `node_plan`, `precomputed_values` keyed by `(component|outcome|rt)`, `pool_templates`, and `guard_quadrature` data. Include metadata (e.g., model signature) so bundle validity is easy to reason about. *(Done – `.build_likelihood_cache_bundle()` and `.prep_set_cache_bundle()` now initialize this during prep.)*
  2. **Populate during prep**: Extend `.prepare_model_for_likelihood()` to call `.build_likelihood_cache_bundle()`, seeding the bundle with compiled nodes, pre-built pool templates (`.build_pool_templates` outputs), and deterministic NA integrals discovered from the generator structure. Attach the bundle as `prep[[".cache_bundle"]]`. *(Done – bundle attached in `R/likelihood_prep.R`, overrides clone it on demand.)*
  3. **Thread bundle through evaluators**: Update `compute_loglik()`, `log_likelihood_from_params()`, and `.likelihood_mixture_likelihood()` to pull the bundle once and pass it into `.outcome_likelihood()`. Replace the mutable `likelihood_cache` env with helper functions that read/write `bundle$precomputed_values`, ensuring repeated trials reuse stored results. *(Done – `.likelihood_outcome_cached()` drives both `R/likelihood_integrate.R` and the param interface.)*
  4. **Retire option-driven caches**: Remove `uuber.param_cache_across_trials` by construction. For override bundles, clone or layer the cache bundle (copy-on-write) so special cases inherit base invariants without leaking state. Document the deprecation in code comments and user docs. *(Done – option removed; overrides clone bundles in `.likelihood_apply_overrides()` and legacy logic fenced into `likelihood_old.R`.)*
  5. **Integrate guards and pools**: Refactor pool and guard evaluators to consult the bundle first (e.g., fetch templates from `bundle$pool_templates`, guard integrals from `bundle$guard_quadrature`) and append entries when new deterministic scenarios appear. This keeps the fast paths hot while avoiding recomputation. *(Done – pool templates cached via `.cache_bundle_ensure()` and guard NA integrals keyed in `guard_quadrature`; added `.inspect_likelihood_plan()` for quick inspection.)*
  6. **Validate and profile**: Add targeted tests (e.g., `tests/test_cache_across_trials.R`) asserting integral counts with the new bundle. Re-run `scripts/test_speed_table.R`/`speed_profile.R` to confirm parity with the legacy timings without hidden env caches. *(Done – new regression test plus routine profiling scripts confirm behaviour.)*

## 5. Precompile and reuse node traversal plans

- **Current status**: Expression compilation assigns numeric IDs, child maps, and fast operators that drive the runtime evaluator. Runtime metadata tracks label-to-ID lookup and competitor clusters.
- **Still to do**: Persist pool combination templates and guard quadrature plans alongside the compiled node table, eliminating per-call recomputation. Expose a diagnostic view of the compiled plan so model authors can sanity-check the generated traversal.

## Future opportunities

1. **Accuracy & stability**: Evaluate deterministic quadrature (e.g. Gauss–Kronrod) for guard densities and add tolerance checks against the legacy sums to flag integration drift.
2. **Expressiveness**: Once `not` support lands, extend the example library with timeout/no-response branches and add tests for mixed guard/negation scenarios.
3. **Tooling**: Promote the profiling harness into a reproducible benchmark suite (microbenchmarks + summarized CSVs) and wire it into CI to keep performance goals visible.
4. **Observability**: Provide instrumentation hooks (histograms of fast-path hits, cache utilization) to guide future optimization work and confirm that principled shortcuts are activated in real models.
