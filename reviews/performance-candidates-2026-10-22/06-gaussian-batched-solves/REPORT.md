# 06: combine design and response covariance solves

**Decision: conditional; consider a dense-path-only implementation, and do not
adopt the unconditional replacement as tested.** Complete dense profile fits
showed a modest local benefit, but structured paths showed no clear benefit,
diagonal temporary allocation increased, and the low-rank optimizer reached a
different successful stopping criterion. No production code was edited.

## Basis and prototype

`regress._profile_covariance_fit` already reuses each covariance factor. This
proposal does not eliminate a factorization. In each likelihood evaluation it
currently solves separately for the predictor design, response, and fitted
residual. The prototype stacks `[design, response]`, solves those columns in
one call, and keeps the residual solve separate. This reduces solve calls from
three to two per evaluation, at the cost of stacking and retaining a combined
array. It deliberately avoids computing the inverse residual by subtraction,
which could amplify cancellation for nearly exact fits.

Official [SciPy 1.17 `cho_solve`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.cho_solve.html)
supports reuse of the supplied Cholesky factor. The underlying
[LAPACK `dpotrs`](https://netlib.org/lapack/explore-html/d3/dc8/group__potrs_ga70a04d13ff2123745a26b1e236212cf7.html)
accepts multiple right-hand-side columns, and its implementation calls
[BLAS `dtrsm`](https://github.com/Reference-LAPACK/lapack/blob/master/BLAS/SRC/dtrsm.f).
These interfaces motivate reducing call/scanning overhead; they do not imply
that every shape or structured solver will become faster. [Ho & Ané
(2014)](https://academic.oup.com/sysbio/article/63/3/397/1649891) provide the
phylogenetic Gaussian likelihood context. The experiment retains the existing
model, ML/REML normalization, variance starts/bounds, optimizers, finite checks,
separate residual quadratic and output schema.

The scoped prototype compiles the existing function with two local changes in
each of its dense and structured branches. It restores the original function
and any already imported NWKIT aliases on exit. `candidate06_plugin.py` applies
that exact prototype before existing tests are collected.

## Measurements

All timings cover complete `_profile_covariance_fit` calls, excluding imports,
fixture creation, output comparison and prototype installation. There is one
warmup per implementation, seven alternating paired measurements, and batches
targeting 40 ms (capped at 100 calls). All runtime thread limits are one and
the shared measurement lock serializes this study. Other user-owned jobs were
active outside that lock, so the ratios are exploratory local observations.

Fixtures use seed 606, an intercept plus random predictors, and random
Gaussian residuals. Dense covariance is `0.7**abs(i-j) + 0.2*I`; the fixed-noise
case adds known diagonal variance 0.2. The diagonal component has heterogeneous
variances. The low-rank case has rank 16 and heterogeneous known diagonal noise.
Only the fixed-noise and low-rank cases require iterative variance optimization;
the other cases use NWKIT's existing closed-form single-scale fit.

| Profile workload | N / predictors | Baseline median | Candidate median | Ratio | Solve calls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense, ML | 128 / 3 | 0.492 ms | 0.378 ms | 1.30x | 6 → 4 |
| Dense, REML | 512 / 8 | 8.149 ms | 7.294 ms | 1.12x | 6 → 4 |
| Dense with fixed noise, REML | 128 / 3 | 4.792 ms | 3.880 ms | 1.24x | 63 → 42 |
| Diagonal, REML | 2,000 / 5 | 0.347 ms | 0.363 ms | 0.96x | 6 → 4 |
| Low rank, REML | 512 / 4 | 8.431 ms | 8.179 ms | 1.03x | 123 → 82 |

The 128-observation dense case has nonoverlapping raw ranges: baseline
0.451–0.546 ms and candidate 0.357–0.436 ms. All seven fixed-noise pairs favored
the candidate, but their overall ranges overlap. For dense 512, one pair was
slower (paired ratios 0.97–1.24x); its modest median gain needs a quiet-host
confirmation. Diagonal and low-rank paired ratios span both directions and
provide no useful speedup evidence. These are profile-call results, not CLI or
total-analysis speedups. These fixtures are not measurements over real trees.

Peak Python allocations were, respectively, 932,128 → 934,470 bytes;
14,743,913 → 14,775,609; 696,646 → 699,125; 326,867 → 439,139;
and 474,204 → 490,922. Thus the diagonal case uses about 34% more temporary
allocation despite no measured speed benefit. Tracemalloc excludes some native
allocations and is not total RSS. Raw timings, batch sizes, solve counts/shapes,
allocations and numerical deltas are in [results.json](results.json).

## Equivalence and checks

Dense workloads differed by at most 5.68e-14 absolute across the numeric output
fields; the diagonal workload matched exactly. Convergence, REML and boundary
flags matched in all workloads. The low-rank iterative fit is numerically close
but is **not an exact metadata match**:

- Objective difference: 1.12e-10, about 7.71e-14 scaled.
- Coefficient maximum absolute difference: 6.74e-9.
- Coefficient covariance difference: 5.67e-12.
- Fitted covariance difference: 4.24e-7; inverse residual difference: 1.61e-6.
- Quadratic difference: 4.64e-5, about 1.69e-8 scaled.
- Component variance difference: 1.74e-7; log-variance difference: 3.21e-6.
- Baseline optimizer message: `CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH`.
  Candidate: `CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL`.

The initial harness stopped at this message mismatch. The checker was then
changed to retain both messages explicitly while still requiring the same
success/boundary flags and comparing all numeric fields (`rtol=1e-4`,
`atol=1e-6`). Production and existing test tolerances were never changed. This
failure history is retained in `exploration-initial.log`. It demonstrates
optimizer sensitivity to arithmetic ordering and is a reason to avoid a
blanket structured-path change.

Additional checks covered ML/REML, response scaling by 1e-6 and 1e6 with
sampling variances scaled consistently, nearly exact fits with 1e-8 residual
noise, and covariance condition number 3.2e9. The ill-conditioned and near-exact
dense comparisons matched numerically exactly in these fixtures. Singular
covariance, rank-deficient design, and negative known variance produced the
same exception type and message. No random stream is used inside this helper.

**36 existing tests passed, no skips**, with the prototype installed through
`tools/check.py test`: the complete numerical-invariance file plus six selected
regression nodes (one parameterized), including independent dense/diagonal GLS,
low-rank and nested factors, size limits, and a reconciled REML consumer.
[tests.log](tests.log) contains the exact command. An earlier narrower selection
also passed 11 tests; it is retained separately and is not added to the final
count. Ruff lint and formatting checks passed on the experimental scripts.

## Reproduction and adoption boundary

Baseline commit `d7116f80b9979739c69a5431a5e34b822795408b`, NWKIT 0.43.54;
macOS 27.0.1 ARM64; CPython 3.12.14; NumPy 2.5.4; SciPy 1.17.1. The parent's
temporary environment passed the documented preflight after a same-version
ETE 4.4.0 source rebuild. The original `.venv` was preserved.

From the repository root, with a compatible Python environment:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-candidates-2026-10-22/06-gaussian-batched-solves/benchmark.py
$PY reviews/performance-candidates-2026-10-22/06-gaussian-batched-solves/benchmark.py --tests
```

The script sets thread limits before numerical imports and acquires the shared
lock internally. `--resume` retains completed timing cases and reruns edge checks.

The most defensible follow-up is a small dense-branch-only patch with a
quiet-host repetition on larger and real-tree workloads. Such a narrowed
production implementation was not integrated or timed here. The unconditional
candidate is not recommended: stacking is counterproductive on diagonal data
and offers little evidence of value for low-rank data. No full/release suite,
population calibration, other-platform runs or end-to-end CLI timing was run
for this isolated, unshipped experiment.
