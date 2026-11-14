# Candidate 03: aggregate cluster scores once

**Judgment: adoptable within the tested numerical scope; integration and an idle-host performance confirmation remain.** This is a promising change to the explicit `cluster-hc1` regression path. It does not accelerate the default REML model, ordinary PGLS, or the whole CLI by the factors below.

## Evidence and proposed change

At baseline `d7116f80b9979739c69a5431a5e34b822795408b`, `nwkit/regress.py::_fit_model` constructs one `inverse == event_index` mask per species event, then sums each cluster's score outer product. With N observations, G events and P predictors, the repeated masks cost O(NG), in addition to the score calculation.

The prototype computes `scores[g,j] = sum_i(design[i,j] * weights[i] * residual[i])` for the rows in event g using one weighted `np.bincount` per predictor. `scores.T @ scores` is the same sandwich middle matrix. Its aggregation work is O(NP + GP²), with O(GP) additional score storage. The fitted coefficients, event weights, variance normalization, existing G/(G-P) correction, t inference and output schema remain as implemented by NWKIT.

[Cameron and Miller (2015), equations 10–11 and the score-matrix identity on page 10](https://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf) justify assembling the sandwich matrix from cluster scores. Their discussion also distinguishes finite-cluster corrections; this experiment preserves NWKIT's existing correction rather than selecting a different one. [NumPy's weighted bincount documentation](https://numpy.org/doc/stable/reference/generated/numpy.bincount.html) specifies the required repeated-index weighted summation. [NumPy's add.at documentation](https://numpy.org/doc/stable/reference/generated/numpy.ufunc.at.html) supports a second equivalent implementation, also benchmarked; bincount was faster here.

## Measured results

Whole `_fit_model` time includes pandas extraction, event identification, fitting and all coefficient rows. Fixture creation and imports are excluded. Synthetic seeded fixtures have shuffled repeated event IDs, heterogeneous contrast variances, correlated within-event residuals and 3–8 predictors. All use seed 442. Seven alternating measurements follow warmup; the 200-row cases were remeasured with batches targeting 40 ms because their initial single-call timings were noisy. `results-initial.json` retains that initial evidence; `results.json` replaces only those small cases.

| N / G / P | Weighting | Baseline median ms | Candidate median ms | Ratio |
|---|---|---:|---:|---:|
| 200 / 50 / 3 | event | 0.633 | 0.375 | 1.69× |
| 200 / 50 / 3 | contrast | 0.643 | 0.369 | 1.74× |
| 3,000 / 600 / 8 | event | 9.917 | 2.161 | 4.59× |
| 3,000 / 600 / 8 | contrast | 10.345 | 2.131 | 4.85× |
| 12,000 / 3,000 / 4 | event | 122.589 | 6.999 | 17.52× |
| 12,000 / 3,000 / 4 | contrast | 138.637 | 9.526 | 14.55× |

The large event-weighted case ranges were 116.091–127.982 ms before and 6.655–8.560 ms after. The contrast-weighted candidate had wider spread, 6.903–18.558 ms, versus baseline 129.069–166.383 ms. The parent's process inspection found unrelated CPU-heavy jobs outside our advisory lock. Consequently these ratios are **exploratory host-local observations**, not uncontended release claims. The consistent order-of-magnitude separation on large cases and removal of O(NG) scanning support further adoption work; exact ratios need idle-host confirmation.

Isolated score-kernel ratios are much larger (approximately 38–750×) and must not be quoted as whole-fit gains. Peak tracemalloc allocation for the whole helper was approximately 2.87 MB for both implementations at N=12,000, dominated by existing input/design work. This is Python allocation, not RSS or a complete account of native allocations. The score matrix still adds O(GP) temporary storage.

## Correctness and checks

Every returned coefficient-row field was compared, including standard errors, test statistics, p-values, confidence intervals and metadata. The six workloads' maximum scaled difference was 1.18e-15 (largest unscaled difference 7.11e-14). Categorical omnibus rows also matched. Independent response rescalings by 1e-12 and 1e12 matched to floating precision; the latter had 1.22e-4 maximum absolute difference at its very large units, 3.97e-16 scaled. The same exception types/messages were retained for insufficient events, rank deficiency and zero variance. The existing reference-formula regression also passed with the prototype installed.

The scoped pytest plugin replaces only `_fit_model` in the pytest process, with no shipped source edits:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-candidates-2026-10-22/03-cluster-score-aggregation/benchmark.py
$PY reviews/performance-candidates-2026-10-22/03-cluster-score-aggregation/benchmark.py --short-only
$PY reviews/performance-candidates-2026-10-22/03-cluster-score-aggregation/benchmark.py --tests
```

The script sets one BLAS/OpenMP thread and holds the common exclusive measurement lock. The last command invokes `tools/check.py test -- -p candidate03_plugin tests/test_regress.py tests/test_numerical_invariance.py tests/test_regression_inference_contract.py -k 'cluster or event_weighting or gaussian_fit_is_invariant' -rs`: **8 passed, 141 deselected, no skips**. Ruff lint and format checks passed for the two experimental Python files.

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NWKIT 0.43.54, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6. The parent's fresh environment passed import preflight following a same-version ETE 4.4.0 source rebuild; the original `.venv` was preserved.

## Limits before integration

No production implementation was changed. A production patch should keep source-level regression coverage for unbalanced clusters, near-cancelling scores and predictor-rich designs, retain the existing inferential correction and failure behavior, and confirm the performance result without external CPU contention. These synthetic timing fixtures use balanced cluster counts; they do not establish field-data frequency or total workflow speedup. The change reorders floating-point sums, so bitwise output identity is not promised. The full regression/release suite and cross-platform checks were not run for this isolated experiment.
