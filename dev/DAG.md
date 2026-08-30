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

Checks that define model semantics remain numeric runtime work. These include
empty integration domains, impossible events, distribution domains, finite
observations, censoring, truncation, and numerical convergence.

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

New functionality belongs in semantic lowering or the generic numerical kernel
layer. Once a compiled path replaces an old path, the old path is deleted.

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

ACCUMULATR_PROFILE_R_SCRIPT=$PWD/dev/scripts/profile_workload_mixed.R \
  bash dev/scripts/profile_cpp_simple.sh
```

Benchmark and profile output belongs in the ignored
`dev/scripts/scratch_outputs/` directory, not in this contract.

## Verification

Evaluator changes must preserve the mathematical model, not merely reproduce an
older engine. Relevant analytic and adversarial validation cases must pass.
Benchmarks must cover both simple and complex compiled plans so complexity does
not leak into ordinary models.

Before accepting evaluator work, answer:

1. What structure moved into compilation?
2. What runtime discovery, check, allocation, or old path was deleted?
3. Does the change apply through a general lowering or numerical kernel?
4. Which benchmark rows changed?
5. What did the matching profile identify before and after?
6. Which validation cases establish semantic equivalence?
7. Is there still exactly one semantic execution framework?
