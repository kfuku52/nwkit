# 08: Reuse integer log-gamma terms in count transitions

**Decision: adoptable within the tested log-domain count-kernel scope.**
The improvement applies to rare-family likelihoods that enter the existing
log-domain fallback. Ordinary probability-domain likelihoods do not call this
kernel and receive no demonstrated benefit. Production code is unchanged.
Timing ratios are exploratory local observations because unrelated host jobs
were active; quiet-host confirmation remains necessary for release claims.

## Evidence and proposed change

At baseline `d7116f80b9979739c69a5431a5e34b822795408b`,
`wgd_count_model._birth_death_log_transition` evaluates the exact branching
probability by summing over surviving ancestral lineages. Its nested array
formula calls `gammaln` six times in each of M ancestral-count rows, repeatedly
evaluating integer arguments independent of duplication/loss rates or time.

The prototype evaluates `gammaln(arange(M+2))` once per call. Simple slices
replace five repeated terms, and one local M-by-M lookup replaces the repeated
descendant-minus-survivor term. This changes neither the formula's addition
order nor its `logsumexp` reduction. Nonpositive integer gamma arguments only
occur in impossible entries that the original kernel subsequently masks; they
are represented by the table's `gammaln(0) = inf` before the same mask.
There is no cross-call cache, mutable cached model, or cubic cache. The extra
retained table uses O(M²) storage and is released at return; the existing
summation work remains O(M³), although gamma evaluation falls to O(M).

Primary sources and official APIs consulted:

- [Rabier, Ta & Ané (2014), detecting and locating WGDs](https://doi.org/10.1093/molbev/mst263)
  provides the birth/death plus retained-copy pulse model context. The paper's
  branching and pulse model motivates preserving exact transitions and original
  copies; it does not establish this prototype's performance.
- [NIST DLMF 5.4.1](https://dlmf.nist.gov/5.4.E1) gives the integer gamma/factorial
  relationship behind the log-combinatorial terms.
- [SciPy `gammaln`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaln.html)
  specifies log-absolute-gamma evaluation and its nonpositive-integer poles.
- [SciPy `logsumexp`](https://cache-docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html)
  documents stable log-domain summation, retained unchanged here.
- Repository [WGD count guide](../../../docs/guides/WGD_COUNT.md) and
  [scientific audit](../../../docs/validation/WGD_SCIENTIFIC_VALIDATION.md)
  define the authoritative omitted-tail, ascertainment and state-convergence
  contracts. Existing study data and protocols were not regenerated.

The current web manuals describe SciPy 1.18; this experiment runs SciPy 1.17.1
and verifies its actual outputs. No dependency upgrade was part of the proposal.

## Comparable measurements

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
NWKIT 0.43.54. All BLAS/OpenMP controls were one. The fresh temporary
environment passed the documented preflight after the same-version ETE 4.4.0
source rebuild; the original `.venv` was preserved.

Each case includes the complete kernel or complete fixed-parameter per-family
likelihood operation, including local table construction. Inputs, rates, state
bounds and ascertainment are identical. Imports, fixture construction and
prototype compilation are outside timing. Both tree versions pay the same
scoped function assignment/restoration overhead. There is a warmup per version
and seven alternating paired trials, with short calls batched to at least
30 ms. [results.json](results.json) retains raw samples, batches, spread,
equivalence checks and separately measured peak Python allocations.

| Complete operation | Baseline median | Candidate median | Ratio |
| --- | ---: | ---: | ---: |
| Log kernel, M=20, birth .2 / loss .3 / time .9 | 1.276 ms | 1.111 ms | 1.15x, modest/noisy |
| Log kernel, M=60, same rates/time | 6.364 ms | 4.598 ms | 1.38x |
| Log kernel, M=120, same rates/time | 29.392 ms | 17.784 ms | 1.65x |
| Ordinary 3-tip likelihood, 3 families, M=40 | 0.645 ms | 0.655 ms | 0.98x, no demonstrated benefit |
| Rare 3-tip likelihood with missingness, 3 families, M=86 | 55.059 ms | 36.314 ms | 1.52x, broad outliers |
| Rare 3-tip likelihood with pulse/detection/missingness, 2 families, M=86 | 103.716 ms | 74.820 ms | 1.39x |

The ordinary control invoked the log kernel zero times under both versions;
its small timing difference is noise. The two rare cases invoked it four and
six times, respectively. They use the same birth rate `1e-8`, loss zero and
root mean one, with observed high count 84. The pulse case uses retention .3,
fraction .5, multiplicity two and tip detection (.8,.9,1). These are synthetic
stress cases for the real fallback, not an estimate of its frequency in real
datasets.

M=60 kernel ranges were 6.334–6.655 ms baseline and 4.520–4.749 ms candidate;
M=120 ranges were 28.138–30.676 and 17.017–20.037 ms. Neither pair overlaps.
The complete pulse likelihood's ranges were 101.694–104.855 and
74.190–76.874 ms, also disjoint. The missing-only case has severe host outliers
(relative ranges 91% and 267%), so its exact ratio is less reliable. No claim
is made for total optimization/bootstrap/CLI runtime or all rate regimes.

Peak tracemalloc allocation increases from 948,106 to 1,065,310 bytes for the
M=120 kernel (about 12%). Complete rare likelihood peaks increase from
926,097 to 986,986 bytes and 1,167,286 to 1,227,843 bytes (about 7% and 5%).
These are traced Python/NumPy allocations, not process RSS or complete native
memory. The proposal exchanges a bounded extra table for less repeated work.

## Equivalence, failures and verification

- Every recorded kernel and per-family log likelihood is exactly array-equal
  to baseline; maximum finite absolute difference is **0**.
- Thirty additional kernel combinations (bounds 1,20,86) cover zero rates/time,
  pure birth, pure death, equal and nearly equal rates, underflowing high-copy
  probabilities, long branches and reciprocal `1e±200` rate/time rescaling.
  Each is exactly equal to baseline and also agrees with the independent
  probability-convolution implementation within existing numerical precision.
- Negative duplication, NaN loss and infinite time preserve exception type and
  message. No probability tail is renormalized; zero-event and missing-subtree
  behavior remains covered by existing regressions.
- With the exact prototype installed by `candidate08_plugin.py`,
  `tools/check.py test -- -p candidate08_plugin tests/test_wgd_count_model.py -rs`
  passed **77 tests in 1.25 s**, no skips. This includes independent CTMC,
  convolution/log equivalence, brute-force pruning, Yule-law, pulse, tail,
  missingness, ascertainment, units, limiting and deterministic simulation checks.
- The real command consumer
  `tests/test_wgd_count.py::test_real_cli_and_metadata_no_uncalibrated_probability`
  passed **1 test in 1.64 s**, no skips, through the same prototype and checker.
  It exercises actual fitting, staged TSV/JSON outputs and state diagnostics.
- Ruff lint and format checks passed for both experiment scripts.

No full/release suite, large replicated scientific study, real-data frequency
study, cross-platform benchmark or complete bootstrap timing was run. This is
an isolated implementation evaluation with no inferential model change. A
production integration still needs ordinary source checks and full applicable
delivery validation. Large bounds can make the extra M² table material; no
claim is made beyond the measured bounds or existing supported work limits.

## Reproduce

From the repository root, with the environment described above:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-candidates-2026-10-22/08-wgd-combinatorial-reuse/benchmark.py
$PY reviews/performance-candidates-2026-10-22/08-wgd-combinatorial-reuse/benchmark.py --tests
$PY reviews/performance-candidates-2026-10-22/08-wgd-combinatorial-reuse/benchmark.py --consumer-tests
```

The script sets one runtime thread and acquires the shared exclusive advisory
lock before numerical imports, fixtures, measurement and tests. Existing users'
unrelated CPU jobs remain outside that lock and were not touched.
