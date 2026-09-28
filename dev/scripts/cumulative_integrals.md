# Compiled integral reuse — 13 September 2026

Comparison against `cb62e1f79e350208e318ba4ce340fd8051fcdee7`, on macOS arm64,
R 4.5, Apple Clang 17, C++17 `-O2`, Accelerate, one numerical-library thread.
Hart and the distribution kernels are unchanged. These are native likelihood
timings with prepared parameters, not end-to-end EMC2 sampler timings.

## Implementation

The compiler removes integration-domain gates and moves factors independent
of the integration variable outside the integral, respecting nested bindings
and source-view context. It also identifies mandatory source PDF/CDF factors
whose support supplies an exact lower integration limit.

For a compiled univariate cumulative integral, requests with identical relevant
leaf parameters, trigger state and onsets are sorted by upper limit. Consecutive
intervals are integrated once and their prefix sums answer the requests. This
also applies to eligible inner integrals at different outer quadrature points.
There is no interpolation, distribution switch or model-name recognition.

Groups need at least four requests. Integrands depending on external times,
sequence history, and nonmatching parameter sets retain direct integration.
Leaf-only races never enter this code. Scratch allocation is reused, but numeric
results are invalidated between particles, including independent likelihood calls.

The existing adaptive GK15 remains: `abs=1e-12`, `rel=1e-3`. Shared intervals
divide the error budget by the largest active group size; accumulated error estimates are checked
at each requested prefix, including signed cancellation. Absolute error budgets
account for external multipliers, and cached estimates include their errors.
This is numerical error estimation, not a rigorous uniform accuracy guarantee.
Outer censoring/truncation settings and the public C++ API are unchanged.

Duplicate integral-construction logic and the old late source-factor evaluation
were replaced. There is one integrator, not competing fixed/adaptive/cache modes.

## Posterior comparison

Six Weber subjects: 105, 108, 113, 115, 116, 118. Each has 563–600 observations;
100 evenly spaced saved particles per subject, five timed repetitions. Timings
below are medians; the mean across subjects is **0.459 → 0.041 seconds per 100
particles**, about **11.2× faster**. A subsequent cache-layout cleanup retimed
at 0.466 → 0.041 seconds and produced identical likelihoods.

| Subject | Before (s / 100 particles) | After | Speedup |
|---|---:|---:|---:|
| 105 | 0.3054 | 0.0400 | 7.6× |
| 108 | 0.5859 | 0.0417 | 14.0× |
| 113 | 0.6190 | 0.0394 | 15.7× |
| 115 | 0.2505 | 0.0423 | 5.9× |
| 116 | 0.4979 | 0.0402 | 12.4× |
| 118 | 0.4929 | 0.0410 | 12.0× |

Against the independently tightened, subdivided reference described in
`integrator_comparison.md`, across all 600 particles:

| Absolute total log-likelihood error | Before | After |
|---|---:|---:|
| 99th percentile | 0.03090 | 0.0001434 |
| Maximum | 0.07085 | 0.0001434 |

Errors include original-trial multiplicities. The largest after-change absolute
trial-density error was 0.0001762. Mathematical equations are unchanged; floating
point summation and quadrature subdivisions are not bit-identical to baseline.

An ablation with cumulative sharing disabled cost 0.0436 seconds per 100
particles versus 0.0391 with sharing in that run. Thus **most of the Weber gain
is algebra and support-aware integration**, with roughly another 10% from
sharing. In particular, starting at a known support boundary avoids wasting
panels before onset and missing the onset discontinuity. No ex-Gaussian-specific
code was added to achieve this.

## Other models and non-reuse controls

Twenty perturbed parameter sets and 32 broad RT/outcome combinations per model.
These are stress grids, not additional posterior samples.

| Model | Before (ms / 20 particles) | After | Max absolute density error before → after |
|---|---:|---:|---:|
| Shared gate | 0.801 | 0.198 | 4.96e-5 → 8.62e-10 |
| Inhibition | 0.314 | 0.086 | 1.54e-6 → 1.04e-13 |
| Guarded shared gate | 1.893 | 0.522 | 1.45e-5 → 1.00e-10 |
| Three-deep inhibition | 354.687 | 23.423 | 9.00e-5 → 2.45e-13 |
| k=2 pool/shared gate | 1.073 | 0.578 | 6.48e-7 → 1.00e-10 |
| Gamma inhibition | 12.950 | 2.350 | 1.08e-3 → 5.60e-5 |

Without sharing, the deep model cost 244.7 ms: reuse is a major gain here,
not merely a small addition to compiler simplification. Tight independent
reference rules and the common reporting floor are documented in the earlier
integrator comparison. Sub-floor tails can have large raw log differences
despite negligible absolute density differences; the table uses density error.

Giving every trial distinct parameters removes cross-trial reuse. Current/before
time ratios were 0.61, 0.83, 0.63, 0.19, 0.85 and 1.00 for the six models above.
The deep model can still reuse inner integrals within a trial.

Longer interleaved controls used 1,000 independent particle calls per repetition,
30 randomized rounds. Ordinary 600-trial LBA/RDM/LNR races were within about 1%
of baseline, with identical likelihoods. Their shorter initial measurements had
suggested a 2–5% regression; this did not persist in the interleaved measurement.

This is **not universally regression-free**: a single very early (`rt=0.04`)
deep-model observation cost 4.22 rather than 3.76 microseconds per call, about
12% slower. At `rt=0.6`, it cost 14.6 rather than 28.2 microseconds. Grouping can
lose when the underlying nested calculations are already extremely cheap.
The other one-trial controls were unchanged or faster within timing variation.

## Reproduction and verification

Source `dev/scripts/bench_cumulative.R` from the repository root, then call
`benchmark_cumulative(cases, baseline = "path/to/baseline")`. Both snapshots
must use the current native context interface; historical measurements below
used the harness from their respective revisions. Each case supplies a finalized `spec`, prepared
`data`, a list of particle `parameters`, and original-trial `expand` indices.
Use enough particles per repetition for millisecond-scale timings. The script
compiles baseline and current evaluators separately; each constructs its own
native context with matching headers. Compilation, parameter mapping and context
construction are outside the timer. No package is installed by the benchmark.

Prepared input snapshots and strong references are temporarily retained in
`/tmp/accumulatr-rules.NSXOuF/`; experiments, ablations and raw timings are in
`/tmp/accumulatr-cumulative.tsmdVf/`. Important files are `reproduced.rds`,
`full.rds`, `no_sharing.rds`, `final_controls.rds`, and `layout_controls.rds`.
The earlier integrator report records the reference construction and cross-checks.

The existing validation suite is retained: 73/73 checks across 31 cases passed,
without a case exceeding the 110-second limit. Existing testthat tests and all
30 compiler-architecture acceptance checks passed. No tests were added or removed.
The installed-package likelihood and EMC2's normal `calc_ll_manager` path matched
the native benchmark exactly on all six subjects and 600 particles.

## Cleanup pass

Integral and top-level evaluation now share execution dispatch and signed-result
cleanup. Error convergence has one definition. Removed the two redundant integral
constructor wrappers, cumulative-bound backup/restoration, and setup for inactive
sharing requests. Completed quadrature panels are summed linearly instead of
repeatedly removing them from a heap. The obsolete `bench_integrators` code was
deleted; its historical report remains. Two unused test calculations were removed,
but no test cases were deleted or added.

Against the pre-cleanup implementation, Weber and ordinary-race timings remained
within 1%; the additional composite controls ranged from 5.4% faster to 2.4%
slower in the final run. Maximum total log-likelihood difference was `9.8e-12`.
In 18 shared/grouped/unique-parameter checks, reordering trials changed densities
by at most `3.8e-13`, and returning to an earlier parameter matrix reproduced its
values exactly. Batch versus individual calls agreed within numerical integration
tolerance, not bitwise equality. Raw cleanup results are in
`/tmp/accumulatr-cleanup.9eoSJd/`.

## Likelihood-wide audit — 17 September 2026

No likelihood implementation or installed package was changed in this audit.
The existing `bench_cumulative.R/.cpp`, example models, EMC2 comparison and
mixed likelihood profiler were reused. Temporary inputs/results are in
`/tmp/accumulatr-audit.aQaXkC/`. The validation suite is unchanged; it was not
rerun because this pass changed no production code.

### Coverage beyond the original workloads

Thirty synthetic cases, 100 trials and 20 perturbed particles each, versus the
same `cb62e1f` baseline. Initial timings used three repetitions. Sixteen cheap
or apparently unchanged cases were then repeated with 4,000 independent
particle calls and seven repetitions to resolve timing noise. Context creation
and parameter preparation are outside these C++ timers.

| Likelihood workload | Measured result |
|---|---:|
| Shared-gate race | 2.46× faster |
| Inhibition | 2.11× faster |
| Guarded shared gate | 3.30× faster |
| Observed/latent stop mixtures | 2.92–3.94× faster |
| Selective-stop composite | 15.91× faster |
| Shared-gate missing RT | 7.37× faster |
| Shared-gate truncation, with/without censoring | 7.93–9.00× faster |
| Shared-gate/guarded models with unique trial parameters | 1.47–1.59× faster |
| Ranked chain, random onset, 2-of-12 pool, simple bounded races | Approximately unchanged |

The longer controls ranged from 3.5% faster to 5.0% slower. The cheap ordinary
two-runner lognormal case was 5.0% slower: 1.76 versus 1.67 microseconds per
100-trial particle. This does not establish a universally regression-free change.
It is also not the initial noisy 52% slowdown from a 33-microsecond batch.

The general improvement applies to eligible compiled integrals, not every
likelihood operation. Ranked history disables cumulative sharing/support
clipping; random-onset convolution and pool arithmetic are separate operations.
Outer censor/trunc and missing-RT calculations benefit from faster inner
integrands, but do not deduplicate complete requests across equal parameter rows.

No tested likelihood hit the reporting floor. Maximum baseline/current total
log-likelihood difference was 0.0007086 per 100 trials; the largest individual
log-probability difference was 0.0006991. The largest change in a likelihood
ratio between tested particles was approximately 0.0704%. These are differences
from the previous implementation, **not errors against an independent oracle**.
There were no timeouts. Raw files: `broad_inputs.rds`, `broad_results.rds`,
`broad_summary.rds`, `controls_results.rds`, `controls_summary.rds`.

### Existing EMC2 comparison

Ran `Rscript dev/scripts/bench_emc2_race_models.R`. The only API repair passes
the context as the final `calc_ll` argument instead of an obsolete data attribute.
The run selects ordinary observations; bounded-case definitions remain available.

EMC2: installed `AccumulatR` branch, source commit `47bb4f41`; AccumulatR:
installed current `optimize_lane` implementation. Both use Hart, Accelerate,
single-thread settings and no fast-math. Seven repetitions; 100/1,000 trials;
20,000/2,000 particles per C++ call, respectively. The repeated 16-offset grid is
not posterior sampling. Design mapping, transformations and bounds are inside
both timed C++ calls; R preparation and compilation are outside. LBA absolute
thresholds are matched despite the packages' different B parameterizations.

| Model | Trials | EMC2 µs/particle | AccumulatR µs/particle | AccumulatR / EMC2 |
|---|---:|---:|---:|---:|
| LBA | 100 | 6.65 | 5.30 | 0.80 |
| LBA | 1,000 | 67.00 | 52.00 | 0.78 |
| RDM | 100 | 8.15 | 7.50 | 0.92 |
| RDM | 1,000 | 81.50 | 74.00 | 0.91 |
| LNR | 100 | 3.45 | 3.05 | 0.88 |
| LNR | 1,000 | 33.50 | 28.00 | 0.84 |

Maximum absolute log-likelihood difference per trial: 3.53e-12. This is a
large-particle-batch result, not a claim that independent one-particle calls
are faster. Ordinary races contain no newly optimized integral, so this lead
must not be attributed to the recent cumulative-integral change.

Branch inspection found existing EMC2 optimizations in `origin/cens_trunc2-SS`
(`e9ac0b84`): stop-success memoization and support-window clipping. The descendant
`origin/cens_trunc2-SS-dEXG3mu` retains them. These should not be reimplemented.
No inspected branch combines these changes with the current AccumulatR bridge.
The installed bridge branch lacks the native race censor/trunc likelihood, so
this audit does not report a bounded-data EMC2 comparison. No branch was checked
out, modified or installed; the inspected remote-tracking refs were not fetched.

### Remaining general likelihood opportunities

1. **Keep error budgets local to each reuse group.**
   `exact_compiled_lane_eval.hpp:504` divides both tolerances for every lane by
   the largest matching group. An unrelated singleton can therefore receive
   unnecessarily strict tolerances. Use per-group/per-lane budgets while retaining
   cumulative-prefix error checks. This is a concrete structural inefficiency;
   its cost in mixed shared/unique batches has not been isolated here.

2. **Reuse complete probability/interval requests.**
   Equivalent outer censor/trunc normalizers and missing-RT probabilities still
   repeat work across parameter-identical rows (`exact_interval.hpp:451`,
   `exact_sequence.hpp:182`). Evaluate representatives and scatter/reduce using
   exact equality of all dependencies: measure/outcome, bounds, source parameters,
   onsets, triggers and context. Numeric values remain particle-local. Existing
   inner reuse already helps: a guarded-model marginal-probability check took
   0.4 ms for 25 identical parameter rows and 1.3 ms for 100; it nevertheless
   repeats the same complete requests. These measurements motivate an extension,
   not a claimed speedup from an unimplemented cache.

3. **Persist preparation across independent calls and context switches.**
   The native workspace/schedule already persists for the same context/data;
   do not add another cache for that case. It retains only the last context and
   data object. The existing nine-model, 50-trial mixed profiler consequently
   shows substantial allocation/freeing and schedule rebuilding when contexts
   alternate (`profile_cpp.txt`). That workload is not a same-subject MCMC loop.
   EMC2 separately reconstructs string-to-cell bridge bindings once per
   `calc_ll` call, already amortized across particles. A worker-owned prepared
   evaluator could retain immutable binding/schedule metadata plus its scratch,
   without retaining parameter-dependent numerical results across particles.

4. **Reduce pool and ranked-source repeated work.**
   Pool prefix/suffix tables currently compute all counts although only counts
   below k are consumed (`exact_source_lane_eval.hpp:606`). Truncating their
   coefficient range gives O(nk), instead of O(n²), work/storage without changing
   the equations. Ranked history-bound normalizers are also recomputed inside
   quadrature calls (`:287`); retain them per source/state/particle. Simply enabling
   current cumulative keys for ranked history would be incorrect.

Smaller candidates are constant mixture log-weights, vector exp/log in
observation mixtures, native weighted reduction for standalone summed likelihood,
and support-bound propagation into random-onset convolution. Convolution depends
on the outer time throughout its integrand, so cumulative prefix addition is not
applicable there. These were opportunities identified by the audit; the measured
implementation follows below.

## Further likelihood optimizations — 17 September 2026

This comparison starts from the implementation audited above, **including** its
previous cumulative-integral optimization. The before-source snapshot is
`/tmp/accumulatr-likelihood.VgDqxT/before/`; it is not the older `cb62e1f` baseline.

### Changes kept

- Cumulative error budgets are local to each parameter group and count actual
  nonempty intervals, not duplicate requests. Unrelated lanes keep their normal
  tolerance. Signed-prefix checks and external-factor error scaling remain.
- Numerical probability operations group requests with exactly equal parameters,
  triggers, onsets and bounds, evaluate representatives, then scatter results.
  Outcome/measure/variant are homogeneous within each call. Grouping is rebuilt
  each call; no numerical result is reused across particles. Existing cheap
  closed-form survival endpoints bypass grouping.
- Evaluation scratch belongs to its context, replacing the single thread-local
  last-context slot. Independent borrowers receive separate scratch; locks cover
  ownership transfer only. Dataset references are traced by R's garbage collector,
  including cycles. Each workspace still retains only its latest dataset schedule.
- EMC2's existing `AccumulatR` branch caches immutable bridge bindings and bound
  ranges on the recipe, with call-local numeric buffers. Changed metadata or a
  cleared serialized pointer rebuilds the plan. No R recipe schema/API changed.
- Pool evaluation streams positive probability/density recurrences: O(nk) work
  and O(k × lanes) scratch. Removed stored member buffers and quadratic
  prefix/suffix tables: 26 lines added, 96 removed in the pool evaluator.

The ranked-normalizer cache was implemented and compared with a cache-disabled
version. It produced no material gain on the supported ranked workloads, so it
was removed entirely. No simulator changes or new permanent tests were added.

### Measurements against the pre-change implementation

The existing C++ timer now also supports alternating contexts and source-directory
snapshots. `likelihood_workloads.R` supplies 44 cases: existing examples, mixed and
unique parameter rows, missing RT, censor/truncation, pool boundaries and context
switches. Data/model/parameter preparation stays outside timing; every particle
is an independent native call. Cheap workloads were repeated to reach at least
approximately 15 ms per sample; seven repetitions, one numerical-library thread.

| Workload | Before → after, ms per 20 particles | Improvement |
|---|---:|---:|
| Mixed-parameter shared gate | 1.041 → 0.899 | 14% less time |
| Mixed-parameter guarded shared gate | 5.025 → 3.947 | 21% less time |
| Mixed-parameter selective stop | 12.166 → 9.143 | 25% less time |
| Composite truncation, repeated parameters | 42.656 → 2.049 | 20.8× |
| Composite known censoring/truncation | 54.077 → 2.523 | 21.4× |
| Composite unknown censoring/truncation | 79.940 → 2.581 | 31.0× |
| Simple missing RT, repeated parameters | 1.351 → 0.064 | 21.0× |
| Composite missing RT, repeated parameters | 7.824 → 0.574 | 13.6× |
| 2-of-24 pool | 1.402 → 0.478 | 2.9× |
| Alternating three contexts (20 calls each) | 5.724 → 3.720 | 35% less time |

An isolated pool comparison used 22 workloads and nine randomized rounds,
including defective/shared-trigger, unique-parameter and ranked/onset cases.
Speedups were 1.25–1.56× for eight members and 1.43–2.80× for 32 members,
depending on k. Differences were at floating-point rounding scale.

There is **not** a universal zero-overhead guarantee. Fully unique simple
selected-censor requests were 4.1% slower and unique simple missing-RT requests
2.0% slower. Composite unique-request controls ranged from 2.3% faster to 1.3% slower.
Ordinary controls were within about 3%. The initially larger regressions on
closed-form endpoints were removed by keeping their existing direct operation.

The isolated bridge benchmark reused `bench_emc2_race_models.R` constructors and
timed C++ calls, with the same installed AccumulatR evaluator on both sides.
Warm one-particle calls improved 3.3–3.4× at 100 trials and 5.6–6.7× at 1,000;
100-particle batches took approximately 5–11% less time. Cold calls still pay
plan construction (up to 6% overhead on single-particle calls). All compared
likelihood vectors were bit-identical.

Default equations, Hart kernels and base quadrature settings are unchanged.
Largest before/after total-log-likelihood difference in the 44-case sweep was
4.99e-5 per 100 trials, from removing excess refinement in mixed groups.
Tightened adaptive tolerances (`abs=1e-13`, `rel=1e-9`) gave a largest discrepancy
of 0.001145 per 100 trials on a unique-parameter truncated workload; this existing
quadrature error was essentially unchanged by the optimization (before/after
difference below 3e-8 there). Tightened rules are a numerical cross-check, not a
rigorous error bound. Existing validation retains its independent formula checks.
The reference tightens inner and complete-response interval integration; selected
response intervals retain `abs=1e-8`, `rel=1e-6`, and missing-RT outer integration
retains its existing fixed rule.

The six current Weber fits (saved September 14) provide a separate posterior
check: 100 particles each, 557–597 trials, five timed repetitions.

| Subject | Before → after, ms per 100 particles | Time saved |
|---|---:|---:|
| 105 | 42.020 → 39.696 | 5.5% |
| 108 | 41.878 → 39.119 | 6.6% |
| 113 | 39.659 → 38.864 | 2.0% |
| 115 | 43.828 → 42.094 | 4.0% |
| 116 | 40.748 → 40.032 | 1.8% |
| 118 | 50.097 → 48.146 | 3.9% |

For each subject, ten particles were checked against tighter inner quadrature:
the five largest before/after discrepancies plus five evenly spaced others.
The largest total-log-likelihood error was **0.0007063**, versus **0.000000721**
before on that particle (subject 118); largest single-trial log-probability error
was 0.0001178. This is about 0.0707% in the full likelihood. Removing accidental
extra refinement does reduce numerical accuracy here; unchanged equations and
configured tolerances do **not** imply unchanged numerical results. This modest
posterior speedup must not be confused with the much larger repeated-normalizer
speedups above. These current fits differ from the older posterior snapshots
used in the earlier sections.

Reproduction from the repository root:

```r
source("dev/scripts/likelihood_workloads.R")
source("dev/scripts/bench_cumulative.R")
results <- benchmark_cumulative(likelihood_cases,
  baseline = "/tmp/accumulatr-likelihood.VgDqxT/before")
```

The baseline argument also accepts a git revision. Raw final measurements and
accuracy comparisons are in that temporary directory: `final.rds`, `summary.rds`,
`accuracy.rds`, `pool_timing.rds`, `emc_timing.rds`, `posterior.rds` and
`posterior_accuracy.rds`. Context-owned scratch
retains memory until the context is released; it is not an unbounded global cache.

Final verification: all 73 checks across the existing 31 validation cases passed,
all existing testthat tests passed without failures/errors/warnings/skips, and all
30 compiler-architecture acceptance checks passed. Temporary ownership probes
also checked garbage-collection cycles, concurrent leases (not foreign-thread R
calls), changed bridge metadata and serialized-pointer reconstruction.

Both packages were reinstalled: AccumulatR `optimize_lane` and EMC2 `AccumulatR`
(base commit `47bb4f41`, only the existing bridge implementation changed).
The normal `calc_ll_manager` path and serialized reconstruction in forked workers
matched direct installed AccumulatR likelihoods exactly for all 600 posterior
particles. Restart existing R sessions to load the rebuilt libraries.

The repository's existing ordinary EMC2 race benchmark was also rerun, with
Hart on both sides and the design/transform/bounds pipeline included for both:

| Model | Trials | EMC2, µs/particle | AccumulatR bridge, µs/particle |
|---|---:|---:|---:|
| LBA | 100 | 6.60 | 5.20 |
| LBA | 1,000 | 68.50 | 52.50 |
| RDM | 100 | 8.00 | 7.50 |
| RDM | 1,000 | 81.00 | 72.50 |
| LNR | 100 | 3.40 | 3.05 |
| LNR | 1,000 | 32.50 | 27.00 |

These are warm batched timings, not single independent calls. Maximum absolute
log-likelihood difference per trial was 3.53e-12. No EMC2 censor/truncation
comparison was added or inferred from these ordinary-race measurements.
