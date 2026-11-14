# Candidate 01: spectral zero-error lambda profiles

**Verdict: adoptable within the tested zero-error, well-conditioned scope.**
The complete lambda fit is substantially faster, including spectral setup.
Positive or mixed sampling-error fits retain their original implementation.
This is an experimental implementation; production code is unchanged.
Integration should retain the demonstrated guards and add focused interior-fit
and fallback regressions. It must not promise bitwise-identical optimized lambda
coordinates or validity for every ill-conditioned covariance.

## Hypothesis and implementation

For `D = diag(diag(C))`, write the existing lambda covariance as
`D^(1/2) [(1-lambda) I + lambda R] D^(1/2)`, where
`R = D^(-1/2) C D^(-1/2)`. One symmetric eigendecomposition of `R` gives a
diagonal profile for every later lambda. Projecting the response and intercept
once makes each profiled mean, diffusion rate and determinant calculation O(n),
after O(n^3) setup, instead of a Cholesky factorization at every search point.
This algebra is an inference from the documented covariance formula; it does
not introduce a different model or likelihood.

The prototype keeps the exact existing grid, optimizer, `xatol`, null test,
unidentifiability rule, profile interval search and tail tolerance. Both lambda
endpoints are evaluated by the original `rate_fit`, retaining their numerical
values and the positive-definiteness failure of a singular supplied `C`.
The interior calculation subtracts the first response before projection to
avoid a large-offset cancellation. All factors/projections are local to one fit.
No cache survives input mutation or another fit.

Positive/mixed SEs and constant traits call the original function directly.
For a positive-definite matrix whose normalized eigenvalues cannot be resolved
above `64*eps*n*max(abs(eigenvalues))`, the prototype also retains the original
Cholesky profile. This criterion selects the computational method; it neither
reclassifies the matrix rank nor relaxes the original positive-definiteness or
scientific tolerances. A singular original matrix fails before that dispatch.
Finite checking remains enabled for `scipy.linalg.eigh`.

## Complete-operation measurements

All times include constructing the fit-local spectral basis, the complete
existing search and, except the explicitly labelled 512-tip case, profile
intervals. Inputs and imports are outside the timed call. Each measurement
uses warmups, seven alternating paired repetitions and the shared exclusive
measurement lock, with one numerical thread. Every candidate timing is below
every baseline timing in its corresponding fit workload. These are exploratory
host-local measurements with unrelated CPU jobs active, not a general CLI
speedup or an uncontended performance guarantee.

| Workload | Baseline median ms | Prototype median ms | Ratio |
| --- | ---: | ---: | ---: |
| Balanced, 8 tips | 20.294 | 6.593 | 3.08x |
| Balanced, 64 tips | 30.949 | 4.705 | 6.58x |
| Balanced, 256 tips | 285.630 | 8.766 | 32.59x |
| Balanced, 512 tips, no interval | 1485.251 | 31.506 | 47.14x |
| Pectinate, 128 tips, non-ultrametric | 52.137 | 5.431 | 9.60x |

The positive-LR 64-tip, 19-replicate **complete bootstrap** improved **5.81x**
in the final additional run; every fit constructs its own spectral basis.
Its observed LR is 6.172442679748826, and both methods give Monte Carlo P=0.05
with exactly equal final RNG states. The initial 32-tip fixture had LR=0 and
therefore exercised the existing no-draw shortcut. It remains explicitly
labelled as an early-exit control in results.json; its near-unity timing ratio
is not bootstrap acceleration evidence.

Peak traced allocations, measured separately from timing, were
2,734,382→3,172,594 bytes at 256 tips and
10,580,872→12,631,424 bytes at 512 tips: approximately 16% and 19% increases.
The extra eigenvectors/workspace are O(n²). Smaller measured cases used less
or slightly more traced allocation. These figures exclude retained inputs,
some native allocations and total RSS; this proposal is a CPU optimization
with a bounded per-fit memory tradeoff, not a general memory reduction.

## Numerical evidence and existing acceptance criteria

- Every returned fit field was compared, including status, lambda, ML rate,
  mean, likelihoods, LR, P and interval bounds. All measured statuses agreed.
  The ordinary-workload maximum scaled difference was `4.94e-8`, entirely
  dominated by the optimized lambda coordinate. The 8-tip boundary fit matched
  every output exactly; boundary=1 and its profile limits also passed the
  existing recorded-reference regression.
- The non-ultrametric 64-tip audit evaluates eleven identical fixed lambda
  points. An independent general solve obtains the GLS mean and ML rate, then
  `scipy.stats.multivariate_normal.logpdf` evaluates their Gaussian density.
  Maximum candidate differences from that direct MLE were likelihood
  `2.84e-14`, rate `2.22e-16`, and mean `6.67e-16`.
  The relevant unchanged existing direct-optimization test in
  `tests/test_signal.py` uses likelihood absolute tolerance `1e-7`, rate
  relative tolerance `1e-5` and mean absolute tolerance `1e-5`.
  These fixed-profile differences are substantially below those existing
  reference criteria; no looser replacement criterion was introduced.
- Optimized lambda itself can move along the likelihood's numerical plateau:
  the largest observed absolute coordinate difference was `7.13e-8` in the
  extreme-unit interior case. That exceeds the optimizer's `xatol=1e-9`, which
  is its stopping setting rather than a guarantee of agreement across different
  floating-point objective evaluations. Neither that setting nor a scientific
  tolerance was changed. The extreme-unit case's LR difference was `1.82e-12`
  and CI lower-bound difference `1.31e-13`; ordinary-workload LR differences
  were at most `2.28e-13`. These observations do not guarantee identical
  selection on every flat/multimodal profile or every platform.
- A second set of nineteen seeded null refits had LR differences at most
  `2.85e-14`. Exceedance counts were identical at all 24 audited cutoffs,
  including the baseline replicate LRs themselves, using the unchanged
  `1e-12*max(1,abs(observed))` tie rule. This is a numerical tail/RNG check,
  not an empirical population calibration study.
- Positive/mixed SE fits, star trees, constant traits and invalid covariance
  failures retained their original outcomes. A `1e-14` near-singular diagonal
  perturbation triggered the guarded original profile and matched exactly.
  A separately retained `1e-12` perturbation was spectrally resolved, with
  maximum scaled difference `2.67e-15`. Interior fits with a `1e9` response
  offset and power-of-two extreme time/trait units retained status and finite
  results; original-unit absolute parameter errors must be interpreted with
  their recorded scale, not as unit-free errors.

The prototype pytest plugin patched `nwkit.signal_stats.lambda_fit` before test
collection and the real CLI consumer alias. The immutable baseline checkout's
existing `tests/test_signal.py` and `tests/test_signal_calibration.py` then
passed **49 tests, zero failures/skips**, in 50.97 seconds through
`tools/check.py test`. These include independent references, lambda boundary
and profile checks, known/mixed SEs, CLI output/seed/order contracts,
bootstrap reproducibility and small study replay/invariance checks. They are
prototype-exercising tests, although the unaffected K/error-aware cases naturally
retain their baseline algorithms.

## Sources and limits

- [Pagel (1999), historical patterns of biological evolution](https://www.nature.com/articles/44766)
  provides the phylogenetic transformation/ML context. NWKIT's SIGNAL guide
  defines the particular `[0,1]`, ML, sampling-error and interval contracts.
- [Kang et al. (2008), efficient control of population structure](https://pubmed.ncbi.nlm.nih.gov/18385116/)
  is primary background on efficiently fitting covariance-based mixed models;
  this implementation's lambda identity is derived above, rather than claimed
  as an NWKIT result established by that paper. Full PMC text was unavailable
  to this session's browser because of its verification gate.
- [SciPy `eigh`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html)
  documents the symmetric eigensystem, finite checks, driver and workspace
  tradeoffs. The experiment uses the installed SciPy 1.17.1 `evd` driver;
  the current online manual may describe a newer release.
- [SciPy multivariate normal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html)
  supplies the independent density evaluator in the additional audit.

No production integration, complete release gate, new R/phytools execution,
cross-platform LAPACK check, exhaustive near-singular/flat-profile audit or
large population calibration study was performed. The existing scientific
model, search bounds, seeds and acceptance criteria remain unchanged. Larger
than 512-tip or high-memory behavior is unmeasured. Error-aware fits have no
demonstrated speedup and deliberately retain their prior algorithm.

## Reproduction

Baseline is immutable commit `86c8f2c537769d852c17e2c789340f6c9de4ef54`
(NWKIT 0.43.55), archived under the path below. The JSON records the imported
`nwkit/signal_stats.py` path from that snapshot. The environment was CPython
3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6 on macOS ARM64; the shared
protocol records the validated environment and external-load limitation.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export NWKIT_CI_DIAGNOSTICS=0
export PYTHONPATH=/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
EXPERIMENT=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/01-lambda-spectrum/benchmark.py
$PY "$EXPERIMENT"
$PY "$EXPERIMENT" --additional
$PY "$EXPERIMENT" --tests
```

Each command takes the common measurement lock internally. The main
[results](results.json), [additional audit](additional-results.json),
[prototype](benchmark.py), [pytest plugin](candidate01_plugin.py) and
[test log](prototype-tests.log) retain the raw evidence.
