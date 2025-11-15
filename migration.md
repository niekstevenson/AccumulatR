# Parameter-Table Migration Plan

## Objectives

- Move all simulators and likelihood functions to accept only explicit parameter tables (`params_df`) together with the structural model definition.
- Remove hidden parameter overrides inside model specs; the race specification should describe topology only.
- Retain `sampled_pars()` (and `sampled_parameters`) as the authoritative description of the core parameter slots so downstream tools know which base parameters exist before design-matrix expansion.
- Lay groundwork for distinguishing *core* parameters from trial-by-trial variations driven by design matrices or condition-level modifiers.

### Canonical Parameter Tables Are Mandatory

- Every downstream consumer (simulation, response-probability utilities, likelihood evaluation, tests) must receive a fully populated `params_df`. There are no implicit defaults or auto-populated overrides.
- Specs produced via `race_spec()` are purely structural; the DSL no longer carries numeric defaults.
- The canonical table layout includes: `trial`, `accumulator_id` (and/or index), `component`, `role`, per-accumulator distribution parameters (as reported by `sampled_pars()`), and timing slots (`onset`, `q`, `t0`).
- Helper `build_params_df(model, core_params, n_trials, component)` converts core vectors into this table; design-matrix logic will later expand core parameters across trials/conditions.

## Implementation Steps

1. **Define Canonical Parameter Tables**
   - Specify required columns (trial id, accumulator id/index, component, onset, q, all distribution params, etc.).
   - Update DSL builders so they never store numeric defaults; only structural metadata remains in `race_spec`.
   - Provide helpers to convert a `sampled_pars()` vector into a single-row parameter table (for baseline fits) and to expand via design matrices for experiments.

2. **Refactor Caching Layers**
   - *Topology Cache*: keep existing `prepare_model()` output (pools, guards, components); no change.
   - *Parameter Block Cache*: replace override bundles with hashes/keys derived from rows of `params_df`, enabling memoisation of density kernels per (accumulator, parameter set).
   - *Trial Plan Cache*: simplify to reference parameter-table row indices; remove logic that merges per-trial overrides into accumulator definitions.

3. **Simplify Likelihood Evaluation**
   - Remove code that mutates `prep$accumulators` with overrides; each evaluator reads parameters directly from the indexed row.
   - Delete helper paths that synthesise default parameter vectors (since parameter tables supply every value).
   - Keep `sampled_pars()` so users can inspect the “core” parameter structure before expansion; document how design matrices map these cores to per-trial values.
   - Revise response-probability helpers to work off parameter tables (single-row tables for analytic probabilities; expanded tables for trial-wise likelihoods).

4. **Design-Matrix Integration (Future)**
   - Introduce a layer that maps core parameters (as reported by `sampled_pars()`) through design matrices to produce per-trial tables.
   - Allow multiple experimental conditions/components to share the same core parameter names but with different design-matrix rows, providing clean compartmentalisation between model structure and experimental manipulations.

## Expected Simplifications

- No need for `populate_model_defaults()` or other override utilities.
- Shared-parameter groups become design-matrix constraints instead of runtime overrides.
- Likelihood caching logic becomes purely structural + table-index-based, reducing complexity in both R and C++ paths.
