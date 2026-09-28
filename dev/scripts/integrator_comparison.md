# Integration-rule comparison — 12 September 2026

Historical results at `cb62e1f79e350208e318ba4ce340fd8051fcdee7`.
See `cumulative_integrals.md` for subsequent executor changes.

## Conclusion

No tested rule was a broadly faster replacement at comparable accuracy.
The strongest new candidate is a mixed rule: try GK15 once, accept converged
integrals, and send the remaining integrals to globally adaptive GK7 starting
on four equal subintervals. It makes no model-specific decisions.

On the posterior workload it saves 25% of the current GK15 runtime with similar
99th-percentile total log-likelihood error. However, it is slower on five of six
additional model shapes, including 42% slower on the deep-guard example.
I would not make it the general default on this evidence.

## Posterior confirmation

100 evenly spaced saved particles per subject, subjects 105, 108, 113, 115,
116 and 118; 563–600 trials per subject. Three timing repetitions, median per
subject, then mean across subjects. Subjects 108, 115 and 118 were not used in
the initial rule screen. These are the updated September 11 fits, not the older
posterior workload used to choose the existing default.

All timings execute the actual likelihood in C++, including nested integrals.
Parameters are prepared outside the timer. Hart, data, parameters and workspace
reuse are shared. This is not an end-to-end sampler timing.

| Rule | Seconds / 100 particles | Time versus current GK15 | p99 absolute total LL error | Maximum absolute total LL error |
|---|---:|---:|---:|---:|
| Fixed 4 × GL7 | 0.251 | 0.53 | 0.6506 | 0.6755 |
| Fixed GL31 | 0.270 | 0.57 | 2.1401 | 2.3879 |
| GK7, four initial panels, rel=0.003 | 0.351 | 0.74 | 0.0562 | 0.1074 |
| GK15 → GK7, four initial panels, rel=0.001 | 0.353 | 0.75 | 0.0336 | 0.0478 |
| GK7, four initial panels, rel=0.001 | 0.388 | 0.82 | 0.0217 | 0.0422 |
| Current GK15, rel=0.001 | 0.471 | 1.00 | 0.0309 | 0.0709 |

Absolute tolerance is `1e-12` for adaptive candidates. Errors aggregate the
actual trial multiplicities, not unweighted compressed rows. None of these
final candidate runs reached the evaluation limit.

Simply loosening GK15 to `rel=0.01` was also tested: approximately 0.322 seconds,
p99 error 0.2850 and maximum error 0.5585. It was inferior to the small-rule
compromises on accuracy at similar speed. A faster rule's small total LL error
does not guarantee small errors on every trial.

## Generality check

Six additional model shapes, 20 perturbed parameter sets each, 32 deliberately
broad RT/outcome combinations per set. These are stress grids, not simulated
or posterior-predictive data. They cover shared gates, inhibition, guarded
shared gates, deep nested inhibition, a k=2 pool with a shared gate, and gamma
accumulators with an endpoint singularity.

| Model | Current GK15 (ms / 20 particles) | Mixed GK15 → GK7 (ms / 20 particles) | Mixed/current | Current maximum absolute density error | Mixed maximum absolute density error |
|---|---:|---:|---:|---:|---:|
| Shared gate | 0.783 | 0.879 | 1.12 | 4.96e-5 | 4.96e-5 |
| Inhibition | 0.303 | 0.409 | 1.35 | 1.54e-6 | 2.99e-5 |
| Guarded shared gate | 1.758 | 2.094 | 1.19 | 1.45e-5 | 1.55e-5 |
| Deep guard | 334.6 | 473.5 | 1.42 | 9.00e-5 | 4.48e-4 |
| Pool/shared gate | 0.963 | 1.257 | 1.30 | 6.48e-7 | 1.43e-5 |
| Gamma guard | 12.20 | 7.927 | 0.65 | 1.08e-3 | 1.16e-3 |

Fixed rules were particularly poor for the gamma example: GL31's maximum
absolute density error was 0.0639, versus 0.00108 for the current GK15.
The fixed-rule error did not improve monotonically with order or panel count.

For synthetic tails, the additional log-error summaries apply a common
reporting floor of `1e-10`. The production simple-observation path substitutes
`min_ll` for impossible/nonfinite values but can return positive densities below
that floor. Comparing an approximately zero density with a sub-floor positive
density otherwise creates large log differences from negligible absolute
differences. The density errors above are unfloored. No production floor
behavior was changed.

## Scope of the screen

92 distinct configurations were evaluated, including:

- Fixed GL23/31/47/63 and combinations of GL7/11/15 on 2–8 equal panels.
- Globally adaptive GK7/11/15/21/31/41; relative tolerances 0.001–0.03,
  selected 0.0001 runs, and 1/2/4 initial panels.
- Nested Gauss–Patterson, Fejér-II and Clenshaw–Curtis rules. Previously
  evaluated nodes are reused when increasing order; unresolved integrals
  switch to the existing adaptive integrator.
- Fixed nested rules, and polynomial changes of variable (`u²`, `u³`,
  `3u²−2u³`) with the appropriate Jacobian.
- Boost GK31 and tanh-sinh; Cubature's vectorized h- and p-adaptive APIs,
  with one or 16 trial lanes sharing a subdivision scheme.
- Four mixed GK15-screen/GK7-refinement configurations.

Screening used 20 particles on subjects 113, 105 and 116. Methods exceeding
five times the fixed baseline were dropped after the first subject. Individual
screening calls had a ten-second limit, checked inside integrand callbacks.
Confirmation/reference calls used larger limits, never more than 60 seconds.

Boost GK31 and tanh-sinh cost approximately 20× and 43× the fixed baseline
on the first screen. Cubature h with 16 shared lanes cost approximately 19×.
Cubature p cost approximately 8×/59× for one/16 lanes and also exhausted
evaluation budgets. These were discarded. The lane-native nested rules
generally needed more evaluations than GK15 and did not improve the trade-off.

## Reference reliability

An unsplit GK15 run at `rel=1e-10` was not a sufficiently reliable reference:
independent tight GK21 evaluations differed by up to 0.0134 log-density units
on the posterior screen. Tight tolerance alone does not detect missed features.

The final posterior reference instead uses GK31 with 16 mandatory initial
panels, `abs=1e-12, rel=1e-10`. On 20 particles from each of difficult subjects
108 and 113, doubling to 32 panels and tightening to `rel=1e-11` changed total
LLs by at most `5.25e-10` and `1.03e-9`, respectively. No final reference run
reached its evaluation limit. This is an empirical cross-check, not a rigorous
uniform error bound over the parameter space.

The synthetic references use tight GK15 and were cross-checked against tight
GK21. Maximum absolute density disagreement was approximately `1e-10`;
large raw log disagreements were confined to sub-floor tails.

## Sources and experiment files

- [Boost Gauss–Kronrod](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/gauss_kronrod.html): embedded rules, available nodes/weights, and conservative error estimation.
- [Burkardt Gauss–Patterson tables](https://people.sc.fsu.edu/~jburkardt/cpp_src/patterson_rule/patterson_rule.html): nested orders and the tables used here.
- [Trefethen, Gauss versus Clenshaw–Curtis](https://ora.ox.ac.uk/objects/uuid%3A4d851558-1a72-4b74-b7b7-8a02cb22b542): why nominal polynomial exactness does not determine practical performance.
- [Trefethen, Exactness of quadrature formulas, 2022](https://ora.ox.ac.uk/objects/uuid%3A09759af5-daaa-4fce-a3ad-257fe70d90e2): further limitations of choosing a rule by polynomial degree alone.
- [Cubature](https://github.com/stevengj/cubature): vectorized h/p interfaces, smoothness trade-offs, and shared-component error control.
- [GSL integration documentation](https://www.gnu.org/software/gsl/doc/html/integration.html): QNG and CQUAD were researched; GSL/CQUAD was not installed or timed.
- [Endpoint-transform research, 2024](https://www.sciencedirect.com/science/article/pii/S001046552400047X): transformations can help singular integrands, but its tailored-exponent rule was not implemented here; the tested polynomial maps were generic.

Experimental code, method definitions, prepared input snapshots and raw results
are retained temporarily in `/tmp/accumulatr-rules.NSXOuF/`. Key files are
`methods.rds`, `retimed.rds`, `finalists_strong.rds`, `models_final.rds`,
`mixed.rds`, `posterior_inputs.rds` and `model_inputs.rds`. Experimental code is
deliberately outside the package.
