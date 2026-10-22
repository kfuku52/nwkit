# 10: Construct only observed covariance during partial-response PGLS optimization

**Decision: conditional; consider a highly missing dense-path implementation,
not an unconditional replacement for every partial-response fit.** The proposed
arithmetic is compatible with the tested outputs, and a highly missing complete
fit improved locally. Moderate missingness can make the construction slower,
however, and complete-fit memory remains constrained by the full covariance
returned to callers. Production code was not changed.

## Basis and implementation

[Clavel, Escarguel & Merceron (2015), mvMORPH](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12420)
describe multivariate phylogenetic likelihoods, evolutionary covariance matrices,
missing values and measurement errors. Their paper motivates avoiding unnecessary
covariance construction; it does not establish a NWKIT performance result. The
repository's `docs/guides/PHYLOGENETIC_REGRESSION.md` defines the actual contract:
the likelihood uses the observed joint-response submatrix, retains all
cross-trait covariance terms, and does not impute missing values or add
undocumented diagonal noise.

The [official NumPy `kron` indexing identity](https://numpy.org/doc/stable/reference/generated/numpy.kron.html)
gives the direct construction. With Fortran-stacked responses, observed index
`i` corresponds to trait `i // n_tips` and tip `i % n_tips`. Thus the selected
entry of `Sigma ⊗ C` is `Sigma[trait_i, trait_j] * C[tip_i, tip_j]`.

At baseline `d7116f80b9979739c69a5431a5e34b822795408b`,
`multivariate_pgls.fit_multivariate_pgls.state` allocates and fills every
`(n_tips*n_traits)²` covariance element for each likelihood evaluation, and then
selects the observed entries. The scoped prototype constructs the observed
covariance directly for intermediate partial-response evaluations. Complete-data
and sparse likelihood branches retain their existing calculation. Component
order, float64 products, fixed covariance addition, Cholesky validation, ML/REML
normalization, parameter bounds and optimizer settings are retained.

The returned `MultivariatePglsFit.fitted_covariance` is a public output containing
all trait-tip combinations. The prototype therefore invokes the original full
construction once at the final parameter vector. It does **not** substitute an
observed-only array in that field. This requirement substantially limits the
possible whole-fit peak-memory reduction. A source-guarded scoped replacement
fails explicitly if the baseline code changes, and the pytest plugin restores
the original functions after the test process.

## Complete-fit measurements

The same seeded fixture, settings and optimizer are used in each pair. All
times cover the complete `fit_multivariate_pgls`, including validation, index
preparation, every likelihood evaluation and construction of the final full
covariance. Imports, fixture creation, prototype compilation and output comparison
are excluded. One warmup per version precedes seven alternating paired trials;
short fits are batched to target at least 25 ms. The seed is `1010+n_tips+n_traits`.
The covariance is `0.65**abs(i-j)`; these are dense correlated synthetic inputs,
not a survey of real phylogenetic tree shapes.

Every retained tip is forced to retain at least one response, so the actual
missing fraction differs from the requested random mask probability. The table
reports actual missingness, not that requested probability.

| Tips | Traits | Missing | Model | Baseline median | Candidate median | Ratio |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 32 | 2 | 0% | REML | 5.92 ms | 6.35 ms | 0.93x, inconclusive |
| 64 | 2 | 26.6% | REML | 9.88 ms | 10.08 ms | 0.98x, inconclusive |
| 64 | 3 | 53.6% | REML | 23.49 ms | 21.74 ms | 1.08x, inconclusive |
| 192 | 2 | 37.2% | REML | 49.39 ms | 43.58 ms | 1.13x, exploratory |
| 24 | 2 | 29.2% | ML, two components plus fixed noise | 33.08 ms | 36.90 ms | 0.90x, inconclusive |
| 128 | 4 | 71.1% | REML | 110.72 ms | 62.55 ms | 1.77x, exploratory |

The highly missing 128-tip case is the promising whole-fit result. Six of seven
paired trials favored the candidate, with paired ratios 0.88–2.17x. Its ranges
overlap: baseline 89.08–141.10 ms, candidate 60.06–106.84 ms. Other user-owned
CPU workloads were active outside this study's advisory lock. The exact ratios
need an idle-host confirmation; modest effects cannot establish a reliable
speedup. No CLI-wide or ordinary complete-data speed claim is made.

Whole-fit traced Python allocation is essentially unchanged. For the high-missing
case it was 4,537,802 → 4,530,457 bytes, because the final full 512×512 covariance
is still built. Other cases had peaks from roughly 0.13 to 3.39 MB, with no
material reduction. `tracemalloc` excludes some native allocations and is not
process RSS.

## Covariance-construction measurements

These separate kernel timings include output allocation and selection, but
exclude factorization, coefficient estimation, optimization and final covariance
publication. Index-pair setup is excluded here because it is reused across
iterations; the complete-fit timing above includes that setup.

| Tips / traits | Missing | Observed elements | Baseline | Candidate | Ratio | Peak Python allocation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 256 / 4 | 35.6% | 659 | 2.20 ms | 2.87 ms | 0.77x | 16.91 → 10.55 MB |
| 256 / 4 | 71.7% | 290 | 1.43 ms | 0.636 ms | 2.25x | 16.91 → 2.15 MB |
| 512 / 4 | 71.3% | 588 | 4.82 ms | 2.30 ms | 2.10x | 67.24 → 8.42 MB |

At 512 tips the baseline and candidate ranges do not overlap: 4.74–5.00 ms
versus 2.24–2.40 ms. The approximately 87% temporary-allocation reduction in
the highly missing kernels is real within this traced scope, but it is **not**
an 87% reduction in complete-fit memory. At moderate missingness, advanced
indexing and multiplication cost more than full Kronecker construction followed
by slicing; the 256-tip moderate case has nonoverlapping timing ranges in the
unfavorable direction. These controls rule out blanket adoption.

## Equivalence and validation

- Every field of all six complete fits matched **exactly**, including
  coefficients, coefficient covariance, all component trait covariances,
  full fitted covariance, likelihood, optimizer message/convergence, boundary
  warning, response count and ML/REML metadata. The largest absolute numeric
  difference was zero.
- Direct covariance-construction comparisons were bitwise identical. The
  response index maps preserve Fortran stacking and all cross-trait blocks.
- An additional callable covariance with an estimated shape parameter and
  fixed sampling noise matched exactly, including the estimated shape and
  stopping message. No component factory cache behavior was changed.
- An independent observed-data ML calculation using the returned full covariance,
  selected response/design, `slogdet` and `solve` matched the reported likelihood
  with absolute difference zero in the checked missing-response fixture.
- All-missing tips, insufficient observations for a response, and singular
  covariance raised the same exception type and message as baseline. No
  tolerance, search bound or failure check was weakened.
- The exact prototype was installed through `candidate10_plugin.py` before
  test collection. `tools/check.py test` ran the complete numerical-invariance
  and multivariate-numerics files, plus five explicit dense/sparse/missing-response
  and CLI consumer nodes: **58 passed in 3.12 seconds, no skips**. The selection
  includes dense/sparse reference agreement, fixed low-rank covariance,
  rank-deficient partial-response rejection, returned-covariance likelihood,
  retained partial tips and biological-replicate CLI handling. The exact command
  and output are in [tests.log](tests.log).
- Ruff lint and formatting checks passed for both experiment scripts.

## Reproduction and adoption limits

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
NWKIT 0.43.54. The fresh shared temporary environment passed preflight after
a same-version ETE 4.4.0 source rebuild; the incompatible original `.venv` was
preserved. All BLAS/OpenMP thread limits are one. Each command acquires the
common exclusive `fcntl.flock` lock internally.

From the repository root:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-candidates-2026-10-22/10-observed-covariance/benchmark.py
$PY reviews/performance-candidates-2026-10-22/10-observed-covariance/benchmark.py --additional
$PY reviews/performance-candidates-2026-10-22/10-observed-covariance/benchmark.py --tests
```

[results.json](results.json) retains all seven timing samples, batch sizes,
memory readings, exact field comparisons and exception probes. The additional
mode appends the high-missing whole-fit case and independent likelihood check
without rerunning successful earlier timings.

A production implementation should choose a validated high-missing dispatch
criterion and preserve full covariance publication; a universal threshold has
not been established by these few fixtures. An alternative blockwise indexing
implementation may improve moderate cases but was not implemented or claimed.
Cross-platform checks, extreme dynamic ranges, real-data workload prevalence,
full release gates and calibration studies were not run for this isolated
numerically equivalent prototype. A quiet-host whole-fit confirmation and
normal source integration checks remain before adoption.
