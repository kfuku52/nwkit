# 04: Batch grouped predictor covariance projections

**Decision: adoptable within the tested dense grouped-uncertainty scope.**
The isolated prototype preserves the calculation and passed numerical and
existing regression checks. It removes a large Python-loop cost from the
projection helper, and improved one complete small Poisson GLMM fit. Its
temporary-memory increase and quiet-host confirmation need attention during
integration. No production code was changed.

## Evidence and implementation

[Ives, Midford & Garland (2007)](https://doi.org/10.1080/10635150701313830)
derive phylogenetic analyses that retain within-species measurement error;
dropping covariance terms can bias estimates. This supports preserving the
full projected covariance, rather than approximating it by its diagonal.
It does not establish a speedup for NWKIT.

At baseline `d7116f80b9979739c69a5431a5e34b822795408b`,
`phylogenetic_glmm._project_predictor_uncertainty` computes each grouped block
as `S.T @ C[:, :, i, j] @ S` in a nested Python loop. For `n` observations,
`k` selected predictor columns and `d` linear predictors per response,
this requires `n²` Python iterations, small matrix products and slice updates.
Independent within-observation uncertainty uses a separate `n`-iteration loop.

NumPy's official [`matmul` documentation](https://numpy.org/doc/stable/reference/generated/numpy.matmul.html)
defines broadcasting of matrix stacks. The prototype transposes the covariance
view into observation-pair/matrix axes, then applies the same two products to
the stack. It preserves the multiplication order of each block and writes
back through an observation-major four-dimensional view. The joint case runs
in chunks of at most 64 first-observation rows, limiting intermediate arrays;
the local case batches the diagonal blocks. This changes Python overhead,
not the mathematical complexity or model.

The [official `einsum` API](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html)
and [`einsum_path`](https://numpy.org/doc/stable/reference/generated/numpy.einsum_path.html)
were also reviewed: contraction optimization can trade memory for speed and
change contraction order. Batched `matmul` was selected here to retain the
baseline block multiplication order. No speed claim is made for an untested
`einsum` variant.

Continuous single-column, GMRF, grouped low-rank and joint low-rank dispatch
remain unchanged. All covariance terms, observation/trait order, float64
operations, final symmetrization, shape checks, optimizer settings and
scientific limits are retained. The prototype uses source-guarded replacement
of the two loop bodies inside a scoped monkeypatch; it fails explicitly if
those baseline bodies change.

## Measurements

All values are complete projection-helper wall times, including output and
temporary allocation, final symmetrization and per-call chunk setup. Shared
input construction and imports are excluded. There is one warmup and seven
alternating measurements; short calls are batched (6–1,000 calls per sample).
The selected variant is the bounded, chunked version.

| Uncertainty | n | k | d | Baseline median | Chunked median | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cross-observation | 16 | 2 | 1 | 0.619 ms | 0.0089 ms | 69.40x |
| Cross-observation | 64 | 2 | 2 | 10.84 ms | 0.427 ms | 25.39x |
| Cross-observation | 128 | 4 | 3 | 48.48 ms | 3.522 ms | 13.76x |
| Cross-observation | 256 | 3 | 2 | 197.80 ms | 8.239 ms | 24.01x |
| Within-observation | 64 | 2 | 2 | 0.166 ms | 0.0343 ms | 4.85x |
| Within-observation | 256 | 4 | 3 | 1.370 ms | 0.789 ms | 1.74x |

For every helper case, the slowest chunked measurement was faster than the
fastest baseline measurement. External user-owned workloads were active
outside our shared advisory lock. Relative timing ranges were substantial
(baseline 16–61%; chunked 30–155%), so the direction of these effects is better
established than their exact magnitude.

A complete **12-observation Poisson GLMM fit**, with two correlated predictor
uncertainty columns, measured **126.79 ms → 55.92 ms (2.27x)**. It includes
model setup, likelihood optimization, mode finding and coefficient inference;
imports and shared fixture construction are excluded. Every paired candidate
fit was faster, but the ranges overlap (baseline 110–221 ms, candidate 53–140
ms) under external contention. This is one small dense fit, not an overall
CLI speedup or a result for sparse/large multinomial fits.

The unbounded batched alternative was measured as a memory tradeoff. At
`n=256,k=3,d=2`, peak traced Python allocation was baseline **4.06 MiB**,
unbounded **7.00 MiB**, chunked **4.56 MiB**. Across larger joint cases,
chunking still increased the peak by about **12–33%**, but avoided the larger
unbounded temporary. The output remains dense and quadratic in `n*d`.
`tracemalloc` is not total process RSS and may exclude native allocations.

Raw timings, per-case memory and all numerical differences are saved in
[results.json](results.json). No random draws occur in the projection helper;
the fit fixtures use the same seeded data in both variants.

## Equivalence and validation

- All six timed helpers returned **bitwise-identical** full matrices, including
  every cross-observation block and final symmetry.
- Direct edge checks cover single observations, three response dimensions,
  a 65-observation chunk boundary, zero and `1e-200`/`1e200` covariance scales,
  non-contiguous reversed views, duplicate selected columns, multiple terms,
  continuous and structured uncertainties, empty terms and NaN propagation.
  All baseline comparisons were bitwise identical.
- An independent dense linear operator, `I_n ⊗ S.T`, verified both local and
  joint observation-major covariance projection. Its results were also
  bitwise identical in the checked fixtures. See
  [additional-checks.json](additional-checks.json).
- Malformed term counts and covariance shapes raised the same exception type
  and message. No PSD projection, tolerance relaxation or extra fallback was
  introduced.
- In the complete Poisson fit, coefficients, coefficient covariance, modes,
  likelihood and other compared numerical fields had maximum difference **0**;
  convergence, boundary flags, messages and remaining metadata matched.
- Existing tests ran with the prototype installed: **39 passed, 177
  deselected, no skips**, in 4.85 s. These included the complete numerical
  invariance file, dense/sparse grouped-uncertainty comparisons, scalar and
  categorical GLMM consumers, latent categorical predictors and cross-column
  uncertainty. Full output: [prototype-tests.log](prototype-tests.log).
- Ruff lint and formatting passed for the two new Python files.

Cross-platform bitwise identity, dense fits close to singular boundaries,
all optimizer families, large multinomial fits and full release gates remain
unverified. Adoption should retain the chunked memory bound, run ordinary
integration/static checks, and confirm representative fits on a quiet host.
The broad kernel gain should not be substituted for that integration evidence.

## Reproduction

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
NWKIT 0.43.54. The temporary environment passed the documented preflight after
ETE 4.4.0 was rebuilt from source at the same version. BLAS/OpenMP thread counts
are one. The original `.venv` was preserved.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python \
  reviews/performance-candidates-2026-10-22/04-uncertainty-contraction/benchmark.py
```

The script acquires the shared `fcntl.flock` lock internally. `--checks-only`
reproduces the additional direct checks without changing the saved timings.
For existing tests, hold the protocol's shared lock while running:

```sh
PYTHONPATH="$PWD/reviews/performance-candidates-2026-10-22/04-uncertainty-contraction" \
PYTEST_PLUGINS=candidate04_plugin NWKIT_CI_DIAGNOSTICS=0 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python tools/check.py test -- \
  tests/test_ordinary_regression.py tests/test_regress.py tests/test_numerical_invariance.py \
  -k 'grouped_predictor_uncertainty or predictor_posterior or categorical_response_and_predictor or phylogenetic_glmm_response_families or numerical_invariance or latent_categorical_predictor or cross_column_uncertainty' -rs
```

No full/release, R backend, population-calibration study or external service
test was run: this experiment changes only an isolated linear-algebra helper
and preserves the statistical procedure, while production remains unchanged.
