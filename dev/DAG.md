# DAG Evaluator Contract

The likelihood evaluator executes a compiled model. Model structure is resolved
once during construction and context compilation; evaluation performs numeric
work over that compiled representation.

## Trusted evaluation boundary

The boundary starts at `log_likelihood()`. Before it, model construction,
`prepare_data()`, `build_param_matrix()`, and `make_context()` validate and
canonicalize their inputs. After it, evaluation trusts:

- native context and layout pointers;
- compiled integer ids, spans, offsets, and opcodes;
- canonical data and parameter columns;
- workspace sizes derived from the compiled plans.

The hot path must not reconstruct data frames or matrices, inspect names or
classes, resolve string ids, infer model shape, or repair malformed compiled
state. Invalid compiled state is a construction bug, not an evaluation case.

Runtime branches define mathematical support, empty integration domains,
impossible events, censoring, truncation, and numerical convergence.
Parameter-domain validation belongs to parameter preparation.

## Compilation owns structure

Compilation must decide:

- component and observation variants;
- source, pool, trigger, and expression dependencies;
- condition and time relations;
- transition and ranked-state schedules;
- distribution operations and their parameter offsets;
- quadrature kernels, bounds, and scratch layout;
- homogeneous trial groups and lane programs.

Evaluation reads this information. It must not rediscover it through semantic
tree walks, maps, dynamic condition records, or try-one-path-then-fallback
control flow.

## Generality is encoded in plans

There is one semantic framework. A simple model is fast because its compiled
plan contains fewer operations. A complex model adds the operations required by
its declared structure. Generality must not be implemented as a second generic
runtime engine or model-specific likelihood formulas.

New functionality belongs in semantic lowering or the numerical kernel layer.

## Workspace and lane execution

Plans own immutable structure and required capacities. Reusable workspaces own
mutable numeric buffers. Evaluation may resize a workspace when lane capacity
grows, but it must not allocate per-trial structural vectors, maps, or branch
records.

Homogeneous operations execute over contiguous lanes. Distribution kernels may
classify genuinely different numerical regimes, but model identity must not
select special kernels. Vector math may change the final bit relative to scalar
libm; it must preserve the equations and the documented numerical accuracy.

## Performance standard

Evaluation cost should be explainable by the compiled schedule and dominated by
its mathematical work:

- distribution PDF/CDF/survival operations;
- vector transcendental operations;
- algebra nodes;
- quadrature evaluations;
- transition, inclusion, and ranked-state terms.

Performance work requires a matching benchmark and profile. A wall-clock change
without a profile is not evidence about its cause. The retained tools are:

```sh
Rscript dev/scripts/benchmark_speed.R

bash dev/scripts/profile_cpp_simple.sh
```

Benchmark and profile output belongs in the ignored
`dev/scripts/scratch_outputs/` directory, not in this contract.

## Verification

Validate evaluator changes against independent analytic and adversarial
references. Benchmark simple and complex compiled plans, including small races
where call overhead matters. A performance report should state the workloads,
timing boundaries, numerical differences, and profile evidence for the measured
costs. See [validation](validation/README.md) for the available checks and
[integral evaluation](scripts/cumulative_integrals.md) for benchmark usage.
