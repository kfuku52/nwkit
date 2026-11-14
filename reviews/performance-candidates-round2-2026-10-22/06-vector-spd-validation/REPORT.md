# Candidate 06: remove discarded vector-covariance inverse work

**Verdict: conditional; insufficient evidence to prioritize adoption as a
performance change.** The minimal prototype preserves all compared numerical
and failure behavior, but the complete small MV-BM fit did not demonstrate a
speedup. Known-parameter conditioning at 128/256 traits improved by 6–7% in
all seven pairs, with overlapping raw timing ranges. The measured overall
memory differences are only a few kilobytes. A quiet-host confirmation on an
actual consuming workload is needed before making a useful-performance claim.

Production source, existing tests, dependencies, and Git history were not
changed. All numerical runs used immutable baseline
`86c8f2c537769d852c17e2c789340f6c9de4ef54` (NWKIT 0.43.55), isolating this
study from unrelated edits in the live checkout.

## Mechanism and retained checks

`vector_gaussian._propagate` calls `_inverse(covariance, "Branch covariance")`
only for validation and discards both returned values. `_inverse` checks
symmetry, computes a Cholesky factor, solves against an identity matrix to
construct an inverse, and computes its log determinant. The latter two
calculations do not contribute to any posterior or likelihood calculation.
The prototype replaces that one call with the same symmetry check and
Cholesky factorization. It retains the factor-finiteness check formerly made
inside `cho_solve` through `numpy.asarray_chkfinite`. The preceding `_array`
still checks the complete input's shape and finiteness. No finite check,
symmetry tolerance, SPD check, or exception conversion was weakened.

The necessary inverses for observation covariance, exact-coordinate
innovations and root posterior precision remain unchanged. The precision-form
update, its general solve, log determinant, symmetrization and all downstream
arithmetic are also unchanged. There is no new cache or retained storage.
The hypothesis is reduced discarded work, not a new Gaussian pruning method.

[SciPy `cho_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html)
documents the SPD factorization and warns that only the chosen triangle is
used; this is why the separate full symmetry check remains.
[SciPy `cho_solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)
documents the additional system solve and its default finite check.
[NumPy `asarray_chkfinite`](https://numpy.org/doc/stable/reference/generated/numpy.asarray_chkfinite.html)
supports retaining that check without solving against an identity matrix.
[Bastide et al. (2020)](https://arxiv.org/abs/2003.10336) provides primary
context for general Gaussian tree models with missing data and errors. These
sources motivate the proposal; they do not establish its NWKIT performance.
Current web documentation describes newer versions than the installed SciPy
1.17.1; the experiment ran the recorded installed implementation.

## Complete-operation measurements

Each timing includes the complete `condition_vector_tree` call: compiled-tree
construction, observation preparation, pruning and all-node smoothing. Input
fixture construction and imports are excluded. Results include every returned
field and seeded joint posterior draws. Each experiment warms up both paths,
uses seven alternating pairs and one numeric thread, and holds the study's
exclusive measurement lock. Unrelated user CPU jobs remained active.

| Workload | Baseline median (ms) | Prototype median (ms) | Observed ratio |
| --- | ---: | ---: | ---: |
| 8 tips / 2 traits, exact, shared transitions | 2.428 | 1.917 | 1.267x |
| 64 tips / 4 traits, partial/mixed errors, unique transitions | 19.900 | 18.533 | 1.074x |
| 64 tips / 16 traits, noisy, shared transitions | 19.408 | 19.114 | 1.015x |
| 128 tips / 32 traits, partial/mixed errors | 57.481 | 55.136 | 1.043x |
| 32 tips / 64 traits, noisy, flat root | 23.455 | 21.717 | 1.080x |
| 64 tips / 4 traits, zero-branch control | 4.967 | 5.162 | 0.962x |

All raw ranges overlap. The small 8-tip case had particularly large variation
(baseline 1.51–4.44 ms, prototype 1.45–5.28 ms), so its 1.267x ratio is not a
stable speedup claim. The zero-branch control largely bypasses propagation;
its variation is another reason to avoid interpreting modest ratios alone.
See [results.json](results.json) for all raw samples and exact field comparisons.

Because discarded inverse construction is cubic in trait dimension, a bounded
follow-up tested two larger dimensions and an actual complete MV-BM fit.
The follow-up used longer target batches (150 ms where a call was shorter).

| Follow-up operation | Baseline median (ms) | Prototype median (ms) | Ratio |
| --- | ---: | ---: | ---: |
| 16 tips / 128 traits, noisy known-parameter conditioner | 35.071 | 33.047 | 1.061x |
| 16 tips / 256 traits, noisy known-parameter conditioner | 142.593 | 133.608 | 1.067x |
| Complete 8-tip / 2-trait partial/noisy MV-BM fit | 1,823.326 | 1,879.173 | 0.970x |

The high-dimensional conditioner improved in seven of seven pairs, but ranges
overlap (128 traits: 31.91–58.02 vs 30.88–43.81 ms; 256 traits:
137.76–166.63 vs 128.97–152.47 ms). This is suggestive evidence for a small
known-parameter conditioning benefit, not an established whole-fit gain.
The complete fit improved in only four of seven pairs with substantial
overlap. Its result does not establish a regression either.

`condition_vector_tree` accepts these known-parameter dimensions. However,
16 tips do **not** provide enough observed coordinates to estimate a full
128/256-dimensional covariance under the existing MV-BM fit guard. These rows
must not be advertised as CLI parameter-estimation speedups. The actual
8-tip fit explicitly calls `fit_pruning_mvbm`; an untimed audit counted
**12,893 `_propagate` calls**, confirming that it exercises the changed code
rather than a separate complete-case scalar/whitened path. The 64-tip
conditioner audit counted 127 calls. See [call-audit.json](call-audit.json).

## Memory, equivalence and checks

Tracemalloc peaks were essentially unchanged. At 128/256 traits they were
17,173,320→17,171,460 and 68,390,744→68,387,892 bytes. The complete fit was
452,524→449,337 bytes. These small differences do not establish a useful
overall memory reduction. Tracemalloc excludes some native allocations and
is not total RSS.

All nine timed operations returned exactly equal fields. The complete fit
comparison includes every optimizer/status field, fitted covariance and every
node's posterior mean/covariance. The first six conditioner cases also retain
exactly equal 19-draw joint samples at the same seed.

Twenty-nine additional checks cover 1/2/8 traits, exact/noisy/partial
observations, flat/proper roots, wrong covariance shape, asymmetry, indefinite and
singular covariance, NaN/Inf, tiny/huge finite covariance, nearly singular SPD
and accepted symmetry roundoff. Exception types/messages and every valid
result match exactly. A separate fault injection makes the factor contain
NaN and verifies that both paths still raise the same finite-array exception.

Prototype-injected existing checks ran through the repository entrypoint:

```text
python tools/check.py test -- -p candidate06_plugin \
  tests/test_vector_gaussian.py tests/test_vector_fit.py \
  tests/test_asr_multivariate_diagnostics.py -rs
20 passed in 10.13s
```

These include independent dense-joint likelihood and posterior comparisons,
noisy and exact coordinates, flat-root scalar references, exact zero-edge
contraction, general OU stationarity, nearly singular roots, BM/OU fit
comparison against dense implementations, disparate units, pruning routing,
and reproducible multivariate sample CLI behavior. No skips were reported.
Ruff lint and formatting passed for the experimental Python files.

The full release suite, other platforms, replicated scientific studies and
large identifiable high-dimensional parameter-estimation fits were not run.
No model, root semantics, uncertainty interpretation, tolerance, optimizer
bound or random schedule was changed. Normal integration and delivery gates
would still be required if this small simplification is later selected.

## Reproduction

The scripts set numeric-thread variables before imports, enforce the shared
measurement lock and import the pinned source snapshot. JSON records the
loaded source path and original `_propagate` SHA-256. Python 3.12.14,
NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6, macOS 27.0.1 ARM64 were used.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/06-vector-spd-validation
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests
$PY "$DIR/additional_review.py"
$PY "$DIR/audit_calls.py"
```

See [additional-results.json](additional-results.json), [tests.log](tests.log)
and the corresponding scripts/logs. The original environment and all unrelated
checkout changes remain untouched.
