# Fast Path Reference

This note summarizes the algebraic shortcuts used by the modern likelihood in
the places where we bypass the fully general evaluation.  The goal is to keep
the optimisations transparent so we can reason about them (and extend them) in
a principled way.

## Guard density and survival

For a guard expression `guard(reference, blocker)`, with no protectors, the
likelihood of completion at time `t` factors into the reference density and the
blocker survival:

```
f_guard(t) = f_ref(t) * S_blocker(t)
```

where

* `f_ref(t)` is the density of the reference event at `t`, and
* `S_blocker(t)` is the probability that the blocker has not yet completed by
  `t` (conditional on any currently forced completions/survivals).

The fast path implemented in `R/likelihood_kernels.R` uses exactly this
expression.  The helper `.guard_density_fast_builder()` pulls the reference
density and blocker survival from the already compiled node table.  The
associated CDF/survival are then obtained by integrating that density (see
`.guard_fast_survival_from_density()`), matching the general path that integrates
`f_guard`.

When protectors are present we fall back to the general integral because
`S_blocker(t)` depends on the protector ordering.

## Pool templates (k-of-n races)

For a pool with members `m_1 … m_n` and completion threshold `k`, the density of
the `k`-th finisher can be written as the sum over all templates:

```
f_pool(t) = Σ_{i} f_{m_i}(t) * Π_{j in complete_i} F_{m_j}(t) * Π_{l in survive_i} S_{m_l}(t)
```

where each template encodes which members must have completed (`complete_i`) and
which must still be surviving (`survive_i`) when `m_i` finishes at time `t`.

`.build_pool_templates()` enumerates these combinations once and the resulting
templates are cached per `(pool, component, k)` in the bundle.  During
evaluation `.pool_density()` and `.pool_survival()` simply walk those cached
templates, multiplying the member densities/CDFs/survivals that are already
accessible via the node table.  The fast-path helpers
`.pool_density_fast_value()` and `.pool_survival_fast_value()` use the same
algebra, but avoid building templates when no forced sets are present by relying
on the binomial coefficient identities (`pool_coeffs`).

## Competitor clusters

Competitor expressions are clustered by shared sources.  When a cluster contains
no guards we can reuse the product of the component survivals:

```
S_cluster(t) = Π_i S_{expr_i}(t)
```

The helper `.compute_survival_product_cluster()` memoises this product inside
the evaluation state so repeated visits with the same `(component, t)` reuse the
cached value.  If guards are present the cluster falls back to the
scenario-based evaluator (`.joint_survival_for_cluster`) that honours the guard
ordering.

## Cache bundle invariants

The per-prep cache bundle is the only place where we persist fast-path artefacts
across trials.  `bundle$precomputed_values` stores `(component|outcome|rt)`
results (e.g. NA guard integrals), `bundle$pool_templates` stores the
`k`-of-`n` templates mentioned above, and `bundle$guard_quadrature` stores guard
integrals keyed by guard signature, component, upper limit, and competitor set.

Everything else that lives in `likelihood_cache.R` is confined to a single
evaluation and is recreated for each call to `.outcome_likelihood()`.
