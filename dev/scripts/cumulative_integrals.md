# Compiled integrals and likelihood benchmarks

This guide describes integral evaluation and the tools for measuring likelihood
cost. Public model and data conventions are documented in the package reference
and vignettes.

## Integral evaluation

The compiler moves factors that do not depend on the integration variable
outside an integral and derives lower limits from source support where possible.
The resulting kernels execute in batches using adaptive 15-point
Gauss–Kronrod quadrature. Composite integrals use absolute tolerance `1e-12`
and relative tolerance `1e-3`, with a limit of 6,000 evaluations per lane.
Adaptive response-interval calculations use `1e-8` and `1e-6`.

These tolerances control estimated quadrature error. They are not bounds on
error in a complete log-likelihood, particularly when integrals are nested or
densities are close to zero. Accuracy checks should compare independent formulas
or tighter, independently subdivided numerical references.

## Cumulative reuse

For eligible integrals, requests with exactly matching source parameters,
trigger probabilities, and onsets are grouped and sorted by their upper limits.
Groups with at least four requests integrate consecutive intervals and recover
cumulative values by summing those intervals. Each group's tolerance budget is
divided among its nonempty intervals. Accumulated error is checked at every
requested limit, including when signed contributions cancel.

Eligibility is determined by the compiled dependencies. Integrals that depend
on external times or ranked-response history use direct evaluation. Different
parameter groups do not share numerical results. Cumulative reuse can also
apply to inner integrals within a nested calculation.

Integral estimates are cached within an evaluation. Numerical results are
invalidated between parameter evaluations. Contexts own reusable workspaces
and retain the latest prepared-data schedule, so repeated calls can reuse
allocated storage and structural information.

The relevant implementation is in:

- [`compiled_integral_planning.hpp`](../../src/eval/compiled_integral_planning.hpp)
  for integral construction and algebraic simplification;
- [`compiled_math_kernel_planning.hpp`](../../src/eval/compiled_math_kernel_planning.hpp)
  for kernel dependency planning;
- [`exact_compiled_lane_eval.hpp`](../../src/eval/exact_compiled_lane_eval.hpp)
  for execution and support limits;
- [`cumulative_lane.hpp`](../../src/eval/cumulative_lane.hpp) for grouping,
  prefix sums, and accumulated error;
- [`exact_adaptive.hpp`](../../src/eval/exact_adaptive.hpp) for quadrature.

## Comparing source versions

Run the benchmark from the repository root. The harness uses Apple's
Accelerate framework and requires macOS, Rcpp, and a C++17 compiler. Load the
checkout whose R interface will prepare the cases:

```r
pkgload::load_all(".", quiet = TRUE, helpers = FALSE, debug = FALSE)
source("dev/scripts/likelihood_workloads.R")
source("dev/scripts/bench_cumulative.R")

results <- benchmark_cumulative(
  likelihood_cases,
  baseline = "path/to/baseline-checkout",
  repetitions = 7
)
```

`baseline` accepts a Git revision or a directory containing a `src/` snapshot.
`candidate` defaults to the working checkout. Both versions must accept the
same native context interface and prepared data and parameter conventions.
The harness compiles each evaluator in a temporary directory and constructs
its own contexts. It does not install either package. Compilation, context
construction, and parameter preparation are excluded from the timers.

Each named case contains a finalized model in `spec`, prepared observations
in `data`, a list of parameter matrices in `parameters`, and an `expand` index
mapping original trials to prepared trials. An `interleave` list alternates
cases to measure context switching. Results contain median `seconds`,
`total_ll_difference` with original-trial multiplicities, and raw `results`.

Use enough repeated calls for stable timings. Cover small races as well as
composite models, shared and distinct trial parameters, censoring, missing
responses, and context switching. Report numerical differences separately from
runtime differences. Comparing two implementations establishes agreement;
independent mathematical references establish correctness.

For the public R-call benchmark and a native profile, use:

```sh
Rscript dev/scripts/benchmark_speed.R
bash dev/scripts/profile_cpp_simple.sh
```

Keep generated results in `dev/scripts/scratch_outputs/`. See the
[validation guide](../validation/README.md) for analytic and compiler checks.
