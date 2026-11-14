# Candidate 09: ER diagonal-plus-constant log pruning

**Verdict: adoptable within the tested fixed-rate, 32–64-state ER scope;
conditional for free-rate fitting.** Complete MK inference, including the
outside pass and all posterior construction, improved 1.17–1.39x at fixed
rates. Seven of seven paired measurements improved and their raw ranges
were separated in those cases. The small-state and freely fitted cases do
not justify an unrestricted replacement. Production code was not changed.

## Mechanism and applicability

For an actual transition matrix with constant diagonal `a`, constant
off-diagonal `b`, and `a >= b >= 0`, `P = (a-b) I + b 11'`. Therefore its
log-message is

```text
logaddexp(log(a-b) + child_log_inside,
          log(b) + logsumexp(child_log_inside))
```

This changes the transition/message contraction from O(k²) to O(k), after
checking the matrix structure. No probabilities are exponentiated back out of
the log domain and no signed sum is used. Zero diagonal difference or zero
off-diagonal is represented by negative infinity. Both sides remain negative
infinity for impossible child likelihoods.

The descriptor is cached by retained matrix identity only within a single
inside pass. The structure test uses **exact equality**, not the approximate
generator classifier: even a one-ULP off-diagonal change is rejected. A
supplied fixed P is checked directly, regardless of the generator's label.
NonER, unequal diagonal, asymmetric, ordered and mixed-regime matrices retain
the original logsumexp path. The generic fixed-P path still computes its logs
in the original way; this is not the prior round's rejected generic log-P
cache proposal. The original outside algorithm is unchanged.

The measurement prototype uses a 16-state crossover guard to examine that
boundary. Results support a more conservative **32-state minimum for adoption**;
the 16-state result is exploratory evidence, not a recommendation. The
32/64-state timed cases take exactly the same path under either threshold.
Integration should explicitly limit the fast method to fixed-rate inference
until free-rate objective/selection sensitivity has been addressed. The
experimental generic inside patch intentionally exposes that sensitivity.

## Complete consuming operation

Every timed call runs `compute_mk_marginals`, including rate-matrix construction,
the likelihood fit, a fresh inside pass, outside likelihoods, posterior
normalization and all fit fields. Fixtures, imports and CLI parsing are
excluded; algorithm preparation and descriptor classification are included.
Inputs contain missing and ambiguous tips. Fixed rate is `0.15/k`.

| Tips/states | Rate/lengths | Baseline ms | Prototype ms | Ratio |
| --- | --- | ---: | ---: | ---: |
| 32/2 | Fixed, repeated | 12.274 | 12.519 | 0.98x, unchanged-path control |
| 64/8 | Fixed, repeated | 23.611 | 23.717 | 1.00x, unchanged-path control |
| 128/16 | Fixed, repeated | 49.322 | 47.078 | 1.05x, overlapping ranges |
| 128/32 | Fixed, repeated | 58.420 | 49.981 | **1.17x** |
| 128/64 | Fixed, repeated | 82.026 | 59.085 | **1.39x** |
| 128/64 | Fixed, all unique | 101.741 | 80.008 | **1.27x** |
| 32/16 | Free, repeated | 126.168 | 113.671 | 1.11x, numerical selection caveat |

For 32 states, baseline/candidate raw ranges were 55.58–71.74 / 48.15–54.88 ms.
For 64 repeated states they were 81.55–84.78 / 58.21–61.29 ms. All-unique
64-state branches were 101.26–104.22 / 79.37–81.22 ms. The 16-state fixed
case's ranges overlap, so its modest change is not treated as a demonstrated
benefit. The cached 64-state single-message kernel improved 2.43x; this is
separate from the full-inference effects above.

Peak traced allocations were 688,476→645,883 bytes at 32 states and
1,063,759→1,020,874 bytes at 64 repeated states. With all branch lengths unique,
unchanged stored transition/output arrays dominate: 9,384,645→9,380,424 bytes.
These are small reductions in tracemalloc peak, not total RSS or a general
memory optimization. Descriptor records are O(unique matrices); the temporary
classification mask is O(k²), and no dense copy is retained in the cache.

## Numerical, reference and behavior evidence

All fields were compared, including likelihood/rate metadata, transition
matrices, log-inside and inside arrays, normalized posteriors, keys and
nonfinite masks. Fixed-rate maximum scaled difference was `2.66e-15`, where
scaling is `abs(delta)/max(1,abs(baseline))`; posterior differences are thus
ordinary absolute differences. The low-state controls matched exactly.
This algebra reassociates the contraction, so bitwise equality is not claimed.

The unrestricted **free-rate** 16-state example moved the fitted rate by
`1.42e-6` (approximately the same relative scale), the largest posterior by
`2.11e-8`, and log-likelihood by `1.51e-10`. Status and metadata keys matched.
The optimizer, starts, tolerances and bounds were unchanged. Passing existing
tests does not establish identical optimum selection for all flat profiles;
free-rate fitting is consequently excluded from the adoption recommendation.
Free 32–64-state fitting and cross-platform optimizer behavior were not measured.

- Eighty direct log-message cases cover 2/8/16/64 states; zero, `1e-300`,
  ordinary and saturated rates; one/ambiguous/all-missing/impossible and
  dynamic-range child likelihoods. They match the original at `1e-13`
  relative/absolute checks; maximum scaled difference is `2.22e-16`.
- Five exact-structure rejection checks cover one-ULP changes, unequal
  diagonal, negative entries, diagonal below off-diagonal, and NaN.
- An independent four-state, four-tip full enumeration sums **16,384** joint
  assignments without pruning. Its maximum posterior difference is
  `6.66e-16`, log-likelihood difference `8.88e-16`.
- Five additional complete inside-pass controls cover original fixed P,
  one-ULP nonER fixed P, ordered asymmetric P, mixed ER/nonER branch regimes,
  and mutation followed by a new call. NonER fallback and mutation results
  match exactly. Four invalid/zero-rate failure comparisons retain identical
  exception type and message.
- A 32-state, 20-draw stochastic-map comparison at seed 971 retains exactly
  equal output counts. This is a selected seed/serial check, not a guarantee
  of equal draws for every reassociated posterior or all thread settings.
- Prototype-injected existing ASR/discrete-model tests passed **79 tests**,
  with 14 explicitly deselected and zero skips, in 10.69 seconds. The plugin
  forces the same ER method even at two states to exercise the existing tiny
  independent reference/CLI/tiny-rate cases; generic model paths retain their
  original algorithm. No test tolerance was changed.

The test selection excludes broad stochastic-map, RSC and origin diagnostics;
the additional seeded map check above does not replace their full suites.
Full release, optional external backends, empirical datasets, population
calibration, other platforms, and combined-candidate integration were not run.

## Reproduction and primary sources

See [results.json](results.json), [additional-results.json](additional-results.json),
[prototype](benchmark.py), [additional checks](additional_checks.py),
[pytest plugin](candidate09_plugin.py) and [test log](tests.log).

The source is the immutable 0.43.55 archive
`86c8f2c537769d852c17e2c789340f6c9de4ef54`; JSON records the actual imported
snapshot path and baseline function SHA-256. All substantive runs held the
shared exclusive lock, used one BLAS/OpenMP thread, warmups and seven
alternating pairs. Unrelated user CPU jobs were active and preserved, so
measurements remain exploratory host-local evidence. CPython 3.12.14, NumPy
2.5.4, SciPy 1.17.1, pandas 3.0.6, macOS ARM64; complete metadata is in JSON.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/09-er-pruning
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests
$PY "$DIR/additional_checks.py"
```

Scripts set numerical thread variables before imports and take the shared
lock internally. Existing tests use the snapshot's `tools/check.py test`.
The shared [protocol](../PROTOCOL.md) explains source/environment recreation.

- [Felsenstein (1981)](https://link.springer.com/article/10.1007/BF01734359)
  supplies primary likelihood-pruning context; only the publisher's summary
  was accessible, and no claim is made that this proposal implements a new
  algorithm from that paper. The ER identity above is an algebraic inference
  from NWKIT's existing transition formula.
- [NumPy logaddexp](https://numpy.org/doc/stable/reference/generated/numpy.logaddexp.html)
  documents stable two-term logarithmic addition.
- [SciPy logsumexp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html)
  documents stable reduction over the child-state axis. Current web manuals
  can describe newer versions than the measured SciPy 1.17.1 implementation.
- [ASR model guide](../../../docs/guides/ASR.md) remains authoritative for
  complete versus ordered graphs, root priors and scientific contracts.
