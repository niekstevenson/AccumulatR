## Parameter-Table Likelihood Paths

We now have two native entrypoints for evaluating the table-driven likelihood:

1. `native_loglik_from_params_cpp` – the original implementation that expects every trial entry to carry an `override_ptr` built from the parameter rows. Plans constructed this way must be rebuilt whenever the parameter table changes because the override pointer embeds the concrete parameter values.

2. `native_loglik_from_buffer_cpp` – a new entrypoint intended for samplers that already produce a dense parameter table every iteration. Instead of requiring per-trial `override_ptr`s, the plan only needs to record which row indices belong to each trial. At evaluation time the C++ code rebuilds the `TrialOverrides` directly from the parameter data frame by slicing those indices, so the same plan can be reused across thousands of parameter draws.

### R Surface

`log_likelihood_from_params_buffer()` mirrors the existing `log_likelihood_from_params()` signature but:

- requires a buffer-style trial plan (`.likelihood_build_trial_plan(..., build_override_ptr = FALSE, keep_trial_rows = FALSE, keep_component_rows = FALSE)`),
- forwards the parameter data frame, row indices, and trial metadata straight to `native_loglik_from_buffer_cpp`, and
- never falls back to the R path—if the native call cannot be evaluated the function errors so the caller can diagnose the plan/data mismatch.

Both helpers share the same fallback path: if the native batch call fails for any reason we drop back to the legacy R implementation so behaviour stays identical.

### When To Use Which Path

| Scenario | Recommended Path |
| --- | --- |
| Existing `model_tables` callers or code that mutates prep per trial | `log_likelihood_from_params()` / `native_loglik_from_params_cpp` |
| MCMC / optimisation loops that regenerate full parameter frames every iteration | `log_likelihood_from_params_buffer()` / `native_loglik_from_buffer_cpp` |

The buffer path is designed for samplers that already perform steps 1–2 below:

1. Draw a parameter vector.
2. Use design matrices to expand it into a per-accumulator parameter table.
3. Call the likelihood with `structure`, the freshly built parameter table, and the observed data.

Only step 3 needed a slimmer native path; now the sampler can reuse a single trial plan and simply ship the matrix-like parameter frame each evaluation.

### Row Indices And Overrides

`trial_plan$trials[[tid]]$row_index` stores zero-based indices into the parameter table. The new C++ helper `build_trial_overrides_from_indices()` reads the parameter columns in-place (no per-trial data frame copies) and constructs the `TrialOverrides` bundle just-in-time. For backwards compatibility we still honour `override_ptr` if it is present; otherwise the native batch function rebuilds it from `row_index`.

### Caching Expectations

- Row indices are per-trial and deterministic, so a single plan can be reused across an entire MCMC run even though the numeric parameter values change every iteration.
- Cross-trial caches (compiled nodes, guard quadrature, shared-gate templates) still live in the native prep context, so only the per-row overrides are rebuilt on each call.
- NA-integral caching is still per evaluation; future work can hash `(task_id, component, param_signature)` if we want cross-trial reuse within a single likelihood call.

### Next Steps

This is the first step in the “native buffer” roadmap:

1. Introduce the buffer driver (done here).
2. Provide a slim trial plan (`build_override_ptr = FALSE`, `keep_trial_rows = keep_component_rows = FALSE`) so buffer callers do not ship per-trial data frames.
3. Allow samplers to pass a contiguous parameter matrix (instead of an R data frame) and read it directly in C++.
4. Layer additional caching for repeated NA/guard integrals keyed by parameter signatures.

For now the API is deliberately conservative so existing callers can flip between the two helpers. As we shave more R-side scaffolding, the buffer path will become the default for high-volume likelihood calls (e.g., samplers, optimisers, large grids). 
