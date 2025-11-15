## Likelihood Cache Strategy

### 1. Structural Cache (per trial-type)
- Build once per unique combination of model topology and component assignment.
- Contents: compiled outcome expressions, pool templates, mixture metadata, and any other purely structural artefact.
- No parameter overrides stored here; the cache is keyed only by structure (e.g., hash of component/pool layout).
- Reused across all trials of the same trial-type so prep stays small and deterministic.

### 2. Parameter-Dependent Integral Cache (shared across trials)
- Stores only integrals that **do not** depend on the observed time `t` (e.g., NA/timeout mass, guard normalizers).
- Keyed by `(trial-type, params-hash)` where the hash is a deterministic summary of the `params_df` rows used for that trial.
- Allows multiple trials sharing identical parameter rows to reuse expensive calculations without growing unbounded.
- These caches live alongside the structural cache but never mix different parameter slices.

### 3. Within-Trial Scratch Cache
- Lightweight map that exists only for the duration of a single likelihood evaluation.
- Holds expensive intermediates that *do* depend on `t` (e.g., shared-gate pieces when evaluating multiple outcomes of the same trial).
- Reset on every trial to avoid memory growth and to keep the logic easy to reason about.

### Reconciliation With Current Implementation
- Today’s `.runtime$cache_bundle` mixes structural data with serialized override bundles; this plan splits the layers so parameter data never pollutes the structural cache.
- `.likelihood_outcome_cached()` keys must incorporate the `(trial-type, params-hash)` so entries are reused only when both structure and parameter slice match.
- Native kernels should mirror the same structure: compiled DAG per trial-type, hash map for parameter-only integrals, and per-call scratch buffers.
- No reliance on global environments—each cache lives inside the prep/nativ context for a specific structure, keeping memory predictable and making eviction trivial (drop the context).
