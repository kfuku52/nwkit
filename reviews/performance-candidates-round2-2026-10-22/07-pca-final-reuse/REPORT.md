# Candidate 07: reuse the final PCA factor and whitened observations

**Verdict: adoptable within the tested BM/fixed-lambda scope.** Complete
`fit_pca` calls improved about 1.3–1.7x for the principal BM/fixed-lambda
workloads, retaining every returned field exactly. The measured large cases
also reduced peak traced allocations. Estimated-lambda fits showed no
demonstrated speed benefit because their unchanged search dominates. These
are exploratory host-local results, not a general CLI speedup.

Production sources, existing tests and Git history were not changed.

## Mechanism and scientific boundary

The baseline calculates a Cholesky factor, stable GLS centering and whitened
traits for the PCA decomposition, then repeats all three in the final
`_profile_likelihood` call. The prototype retains the unscaled whitened array
before `white / scale` creates the PCA-specific array, and passes it with the
existing factor into an extracted likelihood tail.

That tail is copied directly from the baseline: column-unit normalization,
`svd(..., compute_uv=False)`, full-dimensional rank criterion, ML divisor `n`,
tree log determinant and likelihood expression remain the same. The separate
PCA decomposition, its `n-1` scaling, reproducible component signs and repeated
eigenvalue rule are unchanged. Rank-deficient fixed models still report
likelihood `None`; free lambda still rejects a deficient likelihood. No
whitened values are recovered by reversing the correlation scaling, and no
inverse, changed SVD/rank approximation or new numerical tolerance is used.

The lambda optimizer and every search-point likelihood retain their original
implementation. Only the final repeated calculation is removed. Factors and
whitened arrays stay local to one fit, without a persistent tree/input cache.
See the [PCA guide](../../../docs/guides/PCA.md) for the model and output contract.

## Complete-operation measurements

Fixtures use a balanced hierarchical Brownian covariance, and deterministic
continuous traits generated from it. The non-ultrametric fixture adds unequal
tip variances. Every timed call includes validation, model selection/search,
final factorization, PCA decomposition, scores, loadings, status and reported
likelihood. Fixture construction and imports are excluded. These are complete
Python fit operations, not end-to-end CLI/tree-reading measurements.

| Whole `fit_pca` operation | Baseline median ms | Candidate median ms | Ratio |
| --- | ---: | ---: | ---: |
| BM covariance, 8 tips / 2 traits | 0.1546 | 0.1154 | 1.34x |
| BM covariance, 128 / 4 | 0.3568 | 0.2648 | 1.35x |
| BM covariance, 512 / 4 | 7.8755 | 4.5771 | 1.72x |
| BM correlation, 512 / 16 | 6.4870 | 4.4327 | 1.46x |
| BM correlation, 128 / 32 | 0.6535 | 0.5643 | 1.16x |
| Wide/rank-deficient BM correlation, 32 / 64 | 0.4692 | 0.4360 | 1.08x |
| Fixed lambda=0.4 correlation, 128 / 4 | 0.3550 | 0.2678 | 1.33x |
| Fixed lambda=0 covariance, 512 / 4 | 5.3075 | 3.5094 | 1.51x |
| Estimated lambda correlation, 64 / 4 | 8.8825 | 8.8507 | 1.004x |
| Non-ultrametric estimated lambda covariance, 128 / 4 | 16.7504 | 16.7329 | 1.001x |

Every BM/fixed-lambda pair improved in seven alternating measurements. Most
raw ranges were separated: at 512/16 correlation, baseline 5.679–7.213 ms
versus candidate 3.621–4.764 ms. At 512/4 covariance, however, the baseline's
6.303–9.090 ms and candidate's 3.864–6.388 ms overlap slightly. Host contention
therefore limits precision of that 1.72x estimate. The small wide-case gain is
not grounds for claiming a broad high-dimensional speedup. Both estimated
lambda cases have overlapping ranges and only three of seven faster pairs;
their ratios around one are treated as noise.

Peak traced allocations at 512/4 BM fell **8,467,824→6,302,374 bytes**, and
at 512/16 correlation **8,668,752→6,302,374 bytes**. The retained unscaled
white array costs O(n*p), while avoiding a second factor and centered/whitened
working arrays removes more storage in these cases. At 128/4 BM the peak fell
548,092→397,894 bytes; estimated lambda's search peak was effectively unchanged.
These are separate `tracemalloc` measurements, not complete native memory/RSS.
The dense covariance output remains O(n²); no general linear-memory claim is made.

## Equivalence and existing tests

All fields of all ten timed fits matched exactly: centers, scales, rotations,
loadings, eigenvalues, scores, fitted covariance, lambda values/flags,
likelihood, status and repeated-eigenvalue flags. Thirty-one additional cases
matched the same values or exception class/message. They cover covariance and
correlation modes; fixed lambda 0/0.4/1; non-ultrametric free lambda;
rank deficiency, nearly dependent columns and more traits than tips; repeated
eigenvalues; power-of-two extreme tree/trait units; large trait offsets;
constant/nonfinite input; asymmetric/singular covariance; invalid mode, model
and fixed-lambda combinations.

An independent explicit inverse/GLS plus matrix-normal determinant calculation
differed in log likelihood by **2.84e-14**, below the unchanged existing
independent matrix-normal test's `1e-7` absolute tolerance. This is numerical
reference evidence, not a statistical population calibration experiment.

The scoped plugin applies the prototype before the tests and real CLI module
import their `fit_pca` aliases. The immutable baseline's check entrypoint ran:

```text
python tools/check.py test -- -p candidate07_plugin \
  tests/test_pca.py tests/test_numerical_invariance.py -rs
57 passed in 2.31s
```

There were no skipped cases. PCA checks include recorded phytools 2.3.0
references, independent matrix-normal optimization, GLS orthogonality,
transformation reconstruction, units/rank contracts, CLI artifacts and
failure rollback. The additional shared-invariance cases naturally use their
unaffected baseline algorithms. Ruff lint and formatting pass for both
experimental Python files; results.json uses strict standard JSON.

No full release gate, new live R/phytools comparison, other-platform numerical
backend, empirical-data benchmark or combined-candidate measurement was run.
Integration still needs the normal affected checks and delivery gates. Existing
scientific bounds, seeds, tolerances and search algorithms remain unchanged.

## Sources and reproduction

- [Revell (2009), Size-correction and principal components for interspecific
  comparative studies](https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2009.00804.x)
  establishes the phylogenetic centering/PCA context. The reuse proposal is an
  inference from NWKIT's duplicate calculations, not a result asserted by that
  paper.
- [SciPy `solve_triangular`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_triangular.html)
  documents solving with a supplied triangular factor, matrix right-hand sides
  and finite checking. Checks are retained; this proposal skips a duplicate
  evaluation rather than weakening them.
- [NumPy `linalg.svd`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html)
  documents `full_matrices` and `compute_uv`; both original PCA and
  singular-values-only likelihood procedures remain intact.

Baseline is `86c8f2c537769d852c17e2c789340f6c9de4ef54`, NWKIT 0.43.55.
All imports/tests used its fixed archive, isolated from unrelated live-worktree
edits; results.json records `nwkit.__file__` and target source hashes. Runtime:
CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6, macOS ARM64. Online
manuals may describe newer releases. Other user-owned CPU jobs were active and
were preserved. The shared study lock serialized numerical work, with warmups,
seven alternating pairs, calibrated short-call batches and one numerical thread.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/07-pca-final-reuse
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests
```

The script sets numerical thread variables before imports, forces the snapshot
path, and obtains the shared lock internally. The tests also use the existing
ARM64 Cairo library for actual PNG/PDF/SVG output. [results.json](results.json),
[benchmark.py](benchmark.py), [candidate07_plugin.py](candidate07_plugin.py) and
[tests.log](tests.log) retain the raw data and exact procedure.
