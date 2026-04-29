# Exact likelihood precompilation

The first useful step is not a special-case kernel and not an
observed-data cache. It is to store a flatter exact runtime plan in
`likelihood_context()`.

The current exact plans already contain the semantic result: outcome
scenarios, competitor subsets, source ids, expression ids, relation
templates, and support information. The slow part is that likelihood
evaluation still interprets that structure repeatedly through
`readiness_density()`, `after_survival()`, `scenario_truth_cdf()`, and
`competitor_subset_win_mass()`.

This precompile layer should be model-specific and data-agnostic:

- store source/expression value requests, not numeric values;
- store readiness and survival formulas as flat products or sums of
  products;
- store competitor subsets with same-active scenario lists already
  filtered per target scenario;
- keep mutable numeric storage in per-call/per-step scratch, not in the
  context.

The initial implementation should be brave but small:

1.  Add flat runtime formulas to each `ExactVariantPlan`.
2.  Compile each exact transition scenario into:
    - active source id;
    - relation template;
    - `readiness_cdf`;
    - `readiness_density`;
    - `after_survival`.
3.  Compile competitor subsets once per target outcome, and compile
    per-target scenario views that only contain the same-active scenario
    indices.
4.  At evaluation time, precompute the observed-time competitor mass
    once per target scenario, then only evaluate the readiness-varying
    same-active term inside the quadrature loop.

This avoids observed-data cache growth while cutting repeated
model-structure discovery. It intentionally leaves expression evaluation
itself recursive for now; if the next profile still shows expression
tree walking as material, that can be flattened in a second step.
