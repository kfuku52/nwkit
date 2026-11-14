# 01: Reuse fixed covariance factors for zero-error K permutations

**Decision: adoptable within the tested zero-sampling-error permutation scope.**
This is an isolated prototype; the shipped library was not edited. The measured
effect is large and the numerical/RNG checks passed. Unrelated CPU workloads on
this host could not be stopped, so the timings are exploratory, not release-wide
performance claims; confirm the magnitude on an otherwise idle machine before
publishing a speed claim.

## Basis and current bottleneck

Blomberg, Garland & Ives (2003) formulate K through phylogenetic GLS and test
label exchangeability with randomization. The tree covariance does not change
when the tip values are permuted. This makes its Cholesky factor, the GLS weight
vector, weight sum, and K normalization invariant during a zero-error test.
These are mathematical consequences of the existing formula, not an altered
statistical test. See the [primary paper](https://doi.org/10.1111/j.0014-3820.2003.tb00285.x)
and its [full text](https://pages.cs.wisc.edu/~larget/botany940/Blomberg_etal2003.pdf).
SciPy explicitly provides [cho_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)
to solve repeatedly using a given Cholesky factor.

At baseline commit `d7116f80b9979739c69a5431a5e34b822795408b`,
`signal_stats.permutation_test` calls `k_fit` and `k_statistic` for every draw;
each call reaches `gls` and repeats `np.linalg.cholesky(covariance)` plus the
constant-vector solve. The prototype prepares those quantities once inside
each call and preserves the baseline residual solve, summation order, tolerance,
permutation draws, plus-one p-value correction, and per-draw checks. Its work
changes from approximately `O(B*n^3)` to `O(n^3 + B*n^2)`; dense storage remains
quadratic.

Positive or mixed sampling errors take the original path, because their
rate-adjusted covariance changes with each permutation and must be refitted.
Nonpositive draw counts delegate to the original helper to preserve its
existing return/error behavior. The prepared factor is local to one invocation;
there is no global cache keyed by mutable covariance arrays.

## Measurements

Every timing covers the complete `permutation_test`, including preparation of
the candidate factor and RNG initialization. Both versions exclude shared input
construction and the observed K computation, along with CLI imports and tree I/O.
There is one warmup per version and seven paired measurements, reversing order
on alternating repetitions. All raw values, relative ranges, memory readings,
and numerical checks are in [results.json](results.json).

| Covariance | Tips | Draws | Baseline median | Candidate median | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Balanced BM | 64 | 999 | 105.24 ms | 26.38 ms | 3.99x |
| Balanced BM | 256 | 999 | 528.21 ms | 115.13 ms | 4.59x |
| Balanced BM | 512 | 199 | 556.89 ms | 105.19 ms | 5.29x |
| Ladder BM | 128 | 199 | 22.46 ms | 6.86 ms | 3.28x |
| Star | 128 | 199 | 23.00 ms | 6.94 ms | 3.31x |

The first two cases show substantial timing spread: full baseline relative
ranges are 137% and 83%, respectively; candidate ranges are 189% and 48%.
Nevertheless, the slowest candidate measurement is faster than the fastest
baseline measurement in each of the five cases. The other cases have baseline
relative ranges of 8–22% and candidate ranges of 8–16%. This supports the
direction of the effect while limiting precision of its magnitude.

Peak traced Python allocations are essentially unchanged: about 0.07, 1.01,
4.03, 0.26, and 0.26 MiB for the five cases, respectively. This is not total RSS
and excludes untraced native allocations; no memory reduction is claimed.

## Equivalence and checks

- Every timed call returned the exact same p-value as baseline; RNG states after
  the calls also matched exactly, including the star-tree tie case.
- For 100 additional permutations per covariance case, the maximum absolute K
  difference was **0.0**. Another 100 small-case comparisons exercised ties.
- All-positive and mixed sampling errors preserved exact p-values and RNG state
  through the original path. Singular covariance and constant traits raised the
  same exception type and message. Zero draws returned 1 even for an invalid
  covariance, matching baseline.
- Existing signal formula, independent-reference, error-aware, CLI reproducibility,
  input-order, and calibration-audit tests ran with the prototype installed by
  `prototype_plugin.py`: **43 passed in 46.25 s**, no skips. Full output is in
  [prototype-tests.log](prototype-tests.log).

Near-threshold Monte Carlo comparisons are discontinuous: a changed statistic
can change the integer exceedance count. The prototype deliberately keeps the
same residual solve and operation order; bitwise K identity was observed here.
It should not be combined with a new inverse, batched reduction, or relaxed
floating-point tolerance without repeating these checks. Cross-platform bitwise
identity is unverified.

## Reproduction

Environment: macOS 27.0.1 ARM64; CPython 3.12.14; NumPy 2.5.4; SciPy 1.17.1;
NWKIT 0.43.54. Runtime threads are fixed at one. The original `.venv` was
preserved; a fresh temporary environment passed the parent preflight after ETE
4.4.0 was rebuilt from source at the same version.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python \
  reviews/performance-candidates-2026-10-22/01-signal-factor-reuse/benchmark.py
```

The benchmark acquires the protocol's `fcntl.flock` measurement lock internally.
For tests, hold that same lock in a wrapper while invoking:

```sh
PYTHONPATH="$PWD/reviews/performance-candidates-2026-10-22/01-signal-factor-reuse" \
PYTEST_PLUGINS=prototype_plugin NWKIT_CI_DIAGNOSTICS=0 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python tools/check.py test -- \
  tests/test_signal.py tests/test_signal_calibration.py -rs
```

No full/release gate, live R comparison, or population calibration experiment
was run: production and model behavior were not changed, existing independent
reference tests passed, and the prototype leaves the statistical procedure
unchanged. Integration would still need normal source checks and a focused
regression for factorization count and lifetime. The documented dense size/work
limits are retained; this does not establish linear-time tree signal inference
or a speedup for lambda bootstrap or error-aware K.
