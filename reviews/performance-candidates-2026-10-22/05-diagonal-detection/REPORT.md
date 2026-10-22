# 05: Exact bulk diagonal detection with an early-exit guard

**Decision: adoptable within the tested numeric-array scope, using the guarded
prototype. Reject an unconditional full-array `count_nonzero` replacement.**
Production files are unchanged. The strong effects concern diagonal detection
and `factor_covariance`; they are not measurements of whole CLI analyses.
External CPU workloads were active, so the exact ratios need a quiet-host
confirmation before a release performance claim.

## Sources and hypothesis

The baseline `gaussian.is_diagonal` performs two `np.any` reductions for each
row until it finds an off-diagonal nonzero. For diagonal matrices this entails
`2*n` Python-to-NumPy calls. The exact identity

```text
off_diagonal_nonzero_count = count_nonzero(A) - count_nonzero(diag(A))
```

allows bulk reductions without an `n*n` mask, an approximate tolerance, or an
allocated dense off-diagonal copy. This inference follows from disjoint
diagonal/off-diagonal entries and the existing exact-zero contract.

- [NumPy `count_nonzero`](https://numpy.org/doc/stable/reference/generated/numpy.count_nonzero.html)
  documents truth-based counting and the scalar whole-array result.
- [NumPy `any`](https://numpy.org/doc/stable/reference/generated/numpy.any.html)
  documents the logical reduction, including truth of NaN and infinities.
- [NumPy SIMD implementation overview](https://numpy.org/doc/stable/reference/simd/index.html)
  explains compiled kernels and runtime CPU dispatch. It motivates measuring
  bulk operations, but does not promise a particular speedup for these kernels.

The naive full count has two drawbacks. It always reads the entire dense
matrix, losing baseline's first-row exit, and it evaluates diagonal object
truth values that the baseline never examines. The guarded prototype first
retains a first-row nonzero exit for numeric arrays (`dtype.kind` in `biufc`),
then counts the remaining rows and remaining diagonal. Other dtypes delegate
to the baseline. No cache or changed matrix representation is introduced.

## Comparable measurements

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
NWKIT 0.43.54 at `d7116f80b9979739c69a5431a5e34b822795408b`.
BLAS/OpenMP thread controls are all one. An exclusive advisory lock coordinates
the ten candidate agents; unrelated user-owned host workloads remain outside
that lock and were not modified.

Every case uses identical matrices for baseline, naive bulk, and guarded
versions. Each version warms up and calibrates a loop count to at least 30 ms;
there are seven recorded measurements per version, reversing version order
on alternating repetitions. Values below are per-call medians. Fixture creation,
imports and calibration are excluded. `factor_covariance` rows include its
complete factorization, diagonal positivity validation and dispatch; the small
scoped monkeypatch assignment overhead is shared across versions. Raw samples,
strides, loop counts and allocations are in [results.json](results.json).

| Operation / matrix | Baseline | Guarded | Baseline / guarded |
| --- | ---: | ---: | ---: |
| `factor_covariance`, diagonal 64 | 155.44 µs | 6.48 µs | 24.00x |
| `factor_covariance`, diagonal 512 | 1327.39 µs | 103.76 µs | 12.79x |
| `factor_covariance`, dense SPD 64 | 15.15 µs | 13.90 µs | 1.09x, inconclusive |
| `factor_covariance`, dense SPD 512 | 1709.81 µs | 1747.15 µs | 0.98x, inconclusive |
| `factor_covariance`, late off-diagonal SPD 64 | 206.10 µs | 21.74 µs | 9.48x |
| `factor_covariance`, late off-diagonal SPD 512 | 3109.98 µs | 1799.28 µs | 1.73x |
| Detector only, diagonal 2048 | 6582.77 µs | 1598.68 µs | 4.12x |

At 64 diagonal entries, full-factor relative timing ranges were 2.4% baseline
and 3.6% guarded; at 512 they were 30.2% and 5.1%. Their baseline and candidate
ranges do not overlap. Late off-diagonal factors had wider spreads (up to 118%
baseline and 53% guarded), but their ranges also do not overlap. Dense factor
differences are modest relative to timing spread and do not establish a gain
or a regression. No high-level fit or end-to-end latency claim is made.

The naive bulk detector took **1604.47 µs** for an early-nonzero dense 2048
matrix, compared with **3.30 µs** baseline, a roughly **486x slowdown**. This
is a decisive reason not to replace the current helper unconditionally.
The guarded detector took 2.12 µs on that case; its improvement is small in
absolute terms and not used as the main adoption evidence.

Guarded diagonal detector gains also appeared for a 512 strided view (3.12x)
and a Fortran-contiguous 512 matrix (3.53x). The Fortran case is slower than
naive bulk counting because excluding the first row creates a less favorable
layout, but it still beats the original row loop. Additional layout-specific
dispatch was not implemented or needed to establish the proposed scope.

Peak traced Python allocation stayed negligible for detection (at most about
3.5 KiB). Whole dense factor peaks remain approximately 4.00 MiB at 512 in all
versions. These readings exclude some native allocations and input storage;
the optimization does not reduce the quadratic size of a dense covariance.

## Exactness and regression checks

- All factors in the recorded workloads are **exactly array-equal** to baseline.
  Detector booleans also match exactly.
- Explicit cases cover empty/scalar/vector/non-square arrays, all-zero and
  signed-zero matrices, diagonal and off-diagonal NaN/Inf, the smallest
  representable positive float off diagonal, complex/int/bool/string/object
  dtypes, C/F layouts, strided slices and negative strides.
- A custom object whose `__bool__` raises was placed on the diagonal.
  Baseline and guarded versions return `True` without inspecting it; naive
  bulk counting raises. This illustrates the compatibility cost of applying
  the numeric count identity to arbitrary truth-bearing objects.
- The guarded prototype was applied to `gaussian`, `regress`,
  `measurement_error`, and `ordinary_regression` aliases in the test process.
  `tools/check.py test` ran Gaussian factor tests, all numerical invariance
  tests, diagonal-profile independent-reference cases and two ordinary
  regression consumers: **50 passed in 0.41 s, no skips**. See [tests.log](tests.log).
- Ruff lint and format checks passed for both experimental Python files.

## Reproduction and remaining boundaries

From the repository root, after recreating a compatible environment and the
protocol's lock directory if necessary:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python \
  reviews/performance-candidates-2026-10-22/05-diagonal-detection/benchmark.py
```

The benchmark holds `/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock` internally.
For the regression command, hold the same `fcntl.flock(..., LOCK_EX)` in a
wrapper while invoking:

```sh
PYTHONPATH="$PWD/reviews/performance-candidates-2026-10-22/05-diagonal-detection" \
PYTEST_PLUGINS=candidate05_plugin NWKIT_CI_DIAGNOSTICS=0 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python tools/check.py test -- \
  tests/test_gaussian.py tests/test_numerical_invariance.py \
  tests/test_regress.py::test_diagonal_profile_fit_matches_dense_reference \
  tests/test_ordinary_regression.py::test_independent_correlation_matches_ordinary_least_squares \
  tests/test_ordinary_regression.py::test_custom_covariance_and_model_comparison_are_recorded_in_audit -rs
```

No production source, scientific model, numerical tolerance, random stream or
saved-data schema was changed. Full/release, distribution builds, other
platforms and high-level fit timings were not run for this isolated prototype.
Integration should preserve the first-row exit and nonnumeric fallback, add
direct detector contract cases, and use normal affected source verification.
There is no evidence here for replacing dense arrays with a new sparse format
or raising existing dense-work limits.
