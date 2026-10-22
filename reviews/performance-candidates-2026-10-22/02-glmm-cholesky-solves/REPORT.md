# Candidate 02: use the Cholesky factor in dense GLMM solves

**Judgment: adoptable within the tested finite, dense-helper scope; whole-fit
performance remains conditional.** This is a small, plausible production
optimization, but the experiment does not constitute an integrated release.
No production files were edited.

## Rationale and sources

The dense GLMM Newton and marginal-covariance helpers first compute a lower
Cholesky factor `L`, then call `np.linalg.solve(L.T,
np.linalg.solve(L, rhs))`. NumPy documents that `solve` uses LAPACK `gesv`;
that interface treats the factors as general matrices. The prototype replaces
only these nested calls with `scipy.linalg.cho_solve((L, True), rhs)`.
LAPACK `potrs` uses the existing factor with two triangular solves. This
avoids unnecessary general factorizations; the Newton convergence criteria,
backtracking, optimization bounds, fallback optimizer, precision matrix,
likelihood, and existing bounded Cholesky jitter policy are untouched.

Primary/official references consulted:

- [NumPy `linalg.solve`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html):
  documents `gesv` and vector/matrix right-hand sides.
- [SciPy 1.17 `cho_solve`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.cho_solve.html):
  specifies solving from a supplied Cholesky factor. The prototype keeps
  `check_finite=True`, the documented default.
- [LAPACK `dpotrs` implementation](https://netlib.org/lapack/explore-html/d3/dc8/group__potrs_ga70a04d13ff2123745a26b1e236212cf7.html):
  the lower-factor implementation calls `dtrsm` twice rather than refactoring.
- [Ives and Helmus (2011), PGLMMs](https://ives.labs.wisc.edu/pdf/Ives%26Helmus_2011.pdf),
  DOI [10.1890/10-1264.1](https://doi.org/10.1890/10-1264.1):
  model context for preserving phylogenetic covariance, distinct from the
  numerical implementation proposal.
- Repository guide: `docs/guides/PHYLOGENETIC_REGRESSION.md`, especially the
  Laplace likelihood, finite converged mode, and dense/sparse contracts.

The existing sparse path already factors and solves structured precision
systems. This proposal targets the remaining dense path, including dense
custom covariance inputs. It offers no demonstrated sparse-path benefit.

## Prototype and measured workloads

The scoped AST transformation changes the nested Cholesky solve in exactly
five existing helpers: `_random_mode`, `_scalar_random_mode`,
`_polish_scalar_random_mode`, `_marginal_coefficient_covariance`, and
`_marginal_scalar_coefficient_covariance`. The functions retain the original
module globals, so existing regression-test monkeypatches continue to work.
The context manager restores every original function.

Baseline: `d7116f80b9979739c69a5431a5e34b822795408b`, NWKIT 0.43.54.
Host: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4,
SciPy 1.17.1, pandas 3.0.6, ETE 4.4.0 (same-version source build).
BLAS/OpenMP thread limits were all one. The temporary environment passed
runtime imports; the original, incompatible `.venv` was preserved.

Fixtures use `C[i,j] = 0.8**abs(i-j) + 0.2*I[i,j]` and deterministic
seeds. This is a dense correlated-covariance workload, not a claim about
all phylogenetic tree shapes. Poisson modes have one observation per tip;
multinomial modes have three categories and two random coordinates per tip.
Complete fits use the public `fit_phylogenetic_glmm` function, two fixed
coefficients, the same dense covariance factory, and the default Student-t
penalty. Fixture construction and prototype installation are excluded from
timing; the entire helper or fit call is included.

Each workload has one warmup per implementation followed by seven alternating
paired trials, alternating which implementation runs first. Short calls are
batched to target 50 ms per trial; JSON retains batch sizes, all per-call
times, medians, ranges, and equivalence differences. Memory is separately
measured peak Python allocation, not process RSS or a complete accounting
of native allocations.

| Operation | Tips | Baseline median (ms) | Candidate median (ms) | Median ratio |
| --- | ---: | ---: | ---: | ---: |
| Poisson mode | 64 | 0.842 | 0.638 | 1.32x |
| Poisson mode | 256 | 7.574 | 4.623 | 1.64x |
| Poisson mode | 512 | 31.182 | 21.995 | 1.42x |
| Scalar marginal covariance | 64 | 0.196 | 0.148 | 1.32x |
| Scalar marginal covariance | 256 | 1.393 | 1.007 | 1.38x |
| Scalar marginal covariance | 512 | 6.015 | 3.766 | 1.60x |
| Multinomial mode | 64 | 4.480 | 4.100 | 1.09x |
| Multinomial mode | 256 | 27.056 | 22.109 | 1.22x |
| Multinomial mode | 512 | 160.352 | 116.050 | 1.38x |
| Complete penalized Poisson fit | 24 | 51.963 | 48.472 | 1.07x |
| Complete penalized Poisson fit | 64 | 88.601 | 68.191 | 1.30x |

At 512 tips the baseline/candidate ranges did not overlap for the three
helper workloads. Their paired ratios were respectively 1.20–1.70x,
1.44–1.74x and 1.11–1.50x. Smaller gains often overlap noise.
The 24-tip complete-fit effect is inconclusive. Although all seven
64-tip complete-fit pairs favored the candidate, its absolute ranges
overlap, so 1.30x is exploratory rather than a deployable performance promise.

Our measurement lock serialized this NWKIT study. Other user-owned benchmark
jobs were observed using CPU on the same host; they were preserved. This
contention limits all timing conclusions, especially modest gains. Repeat
the same harness on an otherwise idle host before a release performance claim.

No material allocation saving was demonstrated: the scalar 512-tip helper's
peak was identical at 10,512,204 bytes, and the 512-tip multinomial mode
was approximately 58.80 MB for both versions. Cholesky/matrix storage
dominates these peaks. The proposal is about redundant arithmetic.

## Correctness and failure checks

- Helper output differences were at most `1.14e-13` across tested modes,
  objective values and weights; scalar covariance differences were at most
  `8.68e-19`.
- Complete fits retained convergence and boundary-warning status. The
  largest coefficient difference was `1.50e-6`, covariance difference
  `2.19e-8`, random-mode difference `4.37e-11`, and log-likelihood difference
  `1.01e-7`. The outer finite-difference optimizer can follow a slightly
  different floating-point path even when individual solves agree closely.
- The harness comparison uses `rtol=2e-4`, `atol=2e-6`; it changes no
  production or existing regression-test tolerances. Exact differences
  are recorded rather than hidden behind a boolean.
- Constructed SPD systems with condition numbers from 1 through about
  `1e12` had relative solution differences at most `4.21e-15`. At condition
  `1e12` the candidate relative residual is `5.48e-6`; the problem's
  conditioning remains a limitation, not an accuracy improvement.
- Two materially indefinite precision cases raised `LinAlgError` under
  both versions. NaN/Inf precision probes also raised `LinAlgError` under
  both versions. The existing jitter policy can accept nearly singular
  systems; this experiment neither removes nor broadens that policy.
- **26 existing affected regressions passed in 8.03 seconds with the exact
  prototype installed** through a pytest plugin and `tools/check.py test`.
  They cover supported scalar/categorical families, dense/sparse reference
  fits, group/predictor effects, failed-mode score validation, small-penalty
  comparisons, penalized inference and a real null-bootstrap refit.
  There were no skips in this selection. `validation.txt` records the full
  command and test results.
- Local Ruff lint/format checks passed on the experiment scripts. No full
  release suite, full study calibration, cross-platform timing or real-data
  CLI throughput benchmark was run, because production code was not changed.

## Reproduce

From the repository root with a compatible environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python \
  reviews/performance-candidates-2026-10-22/02-glmm-cholesky-solves/benchmark.py
/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python \
  reviews/performance-candidates-2026-10-22/02-glmm-cholesky-solves/validate.py
```

Both entrypoints take the study measurement lock. The temporary interpreter
path and lock path are machine-local and can be replaced together on a new
host. `results.json` is the machine-readable evidence; `validation.txt`
contains prototype regression results.

For production adoption, replace the five explicit expressions directly,
retain factor validation and solver settings, then run affected numerical,
CLI-consumer and release checks on the integrated change. There is no need
to introduce an AST transformation or monkeypatch into the library.
