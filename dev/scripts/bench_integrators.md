# Historical composite quadrature comparison

The old benchmark was removed because its executor replacement no longer matches
the compiled-integral interface. The measurements below predate cumulative reuse.
See `cumulative_integrals.md` and `bench_cumulative.R` for the current comparison.

## Selected defaults and historical results

Composite integrals: GK15, `abs=1e-12, rel=1e-3`.
Outer censoring/truncation integrals retain `abs=1e-8, rel=1e-6`.

These measurements used 750 saved draws from each of 24 Weber participants
(18,000 draws), with one timed repetition per method/participant. Subject 113
was excluded because its original fit contained the separate ex-Gaussian bug.
The reference was GK15 with `abs=1e-12, rel=1e-10`, not an analytic error bound.

| Method | Seconds / 100 particles | 99th percentile absolute total log-likelihood error | Largest estimated parameter-mean shift (SD) |
|---|---:|---:|---:|
| Fixed GL31 | 0.258 | 0.9314 | 0.6075 |
| GK15, abs=1e-12, rel=1e-3 | 0.342 | 0.0832 | 0.0090 |
| GK15, abs=1e-12, rel=1e-4 | 0.436 | 0.0209 | 0.0079 |
| Previous GK15, abs=1e-8, rel=1e-6 | 0.618 | 0.0156 | 0.0072 |

Parameter-mean shifts came from importance reweighting saved draws, not refits.
The selected setting's worst total log-likelihood error was 1.16, and worst
individual-trial relative error was 68.4%. These results are not a uniform
accuracy guarantee.

## Discarded experiments

Higher fixed orders, independently adaptive half-intervals and a batched GK31
prototype did not improve the trade-off. Scalar
[Boost GK15/GK31](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/gauss_kronrod.html)
cost approximately 12/20 times GL31 in this executor;
[Boost tanh-sinh](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/double_exponential/de_tanh_sinh.html)
exceeded two minutes on one case. Their code and the one-off analysis scripts
were removed after selecting the default. This does not establish that Boost
is intrinsically slow: scalar callbacks forfeit the executor's lane batching.

[GSL CQUAD](https://www.gnu.org/software/gsl/doc/html/integration.html) and
[recent SIMD-oriented tanh-sinh work](https://joss.theoj.org/papers/84c42780ad01b47c0c3b5f1ab5c26260)
were researched but not benchmarked here.
