General racing architecture for RT and choice

Overview
- Likelihood engine now uses a lean, numeric-ID node evaluator with per-trial caching.
- Each expression is compiled into node functions (scenario density, conditional CDF, conditional survival).
- A per-trial evaluation context (ctx) memoizes intermediate results within a trial; no cross-trial caching.

Key points
- Precomputed registry assigns integer IDs to all nodes and stores metadata (kind, child_ids, sources, guard wiring).
- Event and pool evaluators support scenario attribution; pools use prebuilt templates and are cached via ctx.
- Guards integrate on demand using compiled node functions; only sub-calls are memoized within-trial.

Public API
- compute_loglik(model, data): unchanged; one ctx per trial is created internally.
- response_probability(model, outcome_label, component=NULL): unchanged; uses fresh ctx per evaluation.
- observed_response_probabilities(model, component=NULL, include_na=TRUE): unchanged semantics.

Performance/caching
- ctx keys: prefix|node_or_label|component|t_key|forced_complete|forced_survive.
- Integral results are not cached across trials; only sub-calls memoize within a trial.

Migration notes
- Legacy recursive walkers have been removed; all evaluation routes through ID-compiled functions.
- Outcome expression IDs are precomputed in prep under .outcome_expr_ids for O(1) lookup.

Files of interest
- R/likelihood_new.R: compiled node evaluator, ctx, and public APIs.
