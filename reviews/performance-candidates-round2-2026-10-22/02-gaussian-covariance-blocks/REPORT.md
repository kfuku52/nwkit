# Candidate 02: descendant-block scalar Gaussian covariance

**Verdict: adoptable within the tested scope.** Building dense covariance by
LCA descendant blocks gave a 3.40x median improvement in a complete 256-tip
Brownian PGLS fit, with every compared fit field exactly equal. The seven raw
baseline/candidate timing ranges do not overlap for that workload. Large
covariance calls improved 5.00–9.72x across balanced, pectinate and star trees.
These are exploratory local observations under unrelated host CPU load, not a
whole-application or uncontended-machine performance guarantee.

Production code and tests were not edited. The baseline is immutable archive
`86c8f2c537769d852c17e2c789340f6c9de4ef54` (0.43.55) at
`/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline`. `results.json` records both
the loaded `nwkit.__file__` and baseline method SHA-256. This avoids an
unrelated concurrent task's edits in the shared checkout.

## Mechanism and boundary

The existing method already uses binary-lifting LCA and carefully handles
zero, signed, underflowed and extreme affine slopes. This proposal changes
the remaining nested Python pair loop, rather than suggesting a new LCA
implementation. Each requested output position belongs to a node. In reverse
preorder, its node's own positions and child-subtree positions are merged.
Pairs in different groups have that node as their LCA, so an entire cross
block can be filled at once. Incremental group merging avoids quadratic
numbers of tiny NumPy calls at a high-degree polytomy. Child membership arrays
are released when merged; no dense ancestor-index matrix is retained.

The covariance arithmetic is retained: unit slopes use the ancestor's
variance; direct finite slopes multiply the two cumulative loadings and the
same ancestor factor in the original order. Off-diagonal exceptional products
use the original stable logarithmic calculation and output-position order.
If cumulative loadings/factors cannot safely use that path, the original
implementation remains the fallback. A flat root, an outside-tree request,
and overflowed moments retain their original failure messages.

The final proposal dispatches to blocks for **at least 64 requested nodes**.
A 32-request experiment did not establish a benefit (0.849x, overlapping raw
ranges), so it is retained only as a dispatch control. The 16-request workload
uses the original implementation. All principal timing cases use exactly the
same block/original path before and after choosing this final threshold.
The fallback with zero slopes was 0.918x with overlapping ranges: no fallback
speedup is claimed, and it does incur an additional linear setup pass.

Two product arrays plus masks are tiled against a 4 MiB workspace estimate;
the previous block is explicitly released before the next one. Required dense
output storage is still quadratic. Existing sparse tree precision remains
the preferred representation when consumers can use it; this proposal does
not replace it or claim linear complexity for a dense output matrix.

## Comparable measurements

Every row includes complete method/fit setup, seven alternating pairs,
warmups, one numerical thread and an exclusive study lock. Short calls are
calibrated into batches. Input fixture construction and imports are excluded;
the method's own topology/index construction is included. Peak allocations
are measured separately from timing, using tracemalloc, not total RSS.

| Operation | Baseline median (ms) | Candidate median (ms) | Ratio |
| --- | ---: | ---: | ---: |
| 16-tip balanced BM covariance, original-path control | 0.120 | 0.120 | 0.994x |
| 64-tip balanced BM covariance | 1.238 | 0.840 | 1.47x |
| 128-tip balanced BM covariance | 4.545 | 1.720 | 2.64x |
| 512-tip balanced BM covariance | 78.844 | 8.108 | 9.72x |
| 512-tip balanced OU covariance, alpha=0.4 | 104.193 | 13.452 | 7.75x |
| 256-tip pectinate BM covariance | 25.873 | 3.810 | 6.79x |
| 256-tip pectinate OU covariance, alpha=0.01 | 30.644 | 6.123 | 5.00x |
| 512-tip star BM covariance | 52.819 | 6.644 | 7.95x |
| Complete 64-tip Brownian REML PGLS fit | 1.728 | 1.391 | 1.24x |
| Complete 256-tip Brownian REML PGLS fit | 22.531 | 6.628 | 3.40x |

For the complete 256-tip fit, raw baseline times are 22.03–22.77 ms and
candidate times 6.37–6.84 ms. Larger effects are clear on this host; smaller
cases still warrant a quiet-host confirmation. Raw timings, loop calibration,
all output comparisons and allocation peaks are retained in `results.json`.

At 512 tips, traced peaks are 2,511,224→2,515,668 bytes for BM and
2,519,648→3,085,576 bytes for OU. The complete 256-tip fit is
4,223,901→4,228,222 bytes. The 64/128-tip covariance cases add about 29/78 KB,
respectively, so no general memory reduction is claimed.

## Equivalence and checks

All timed matrices and every returned key/field of both complete fits are
exactly equal. Twenty-one direct checks cover arbitrary internal nodes,
shuffled order, repeats, repeated root, empty requests, signed affine slopes,
zero slopes, extreme variances (1e±280), cumulative-loading underflow, flat
roots, outside-tree nodes and moment overflow. Eight additional random affine
processes contain signed slopes and shuffled duplicate requests; an OU process
also uses all nodes twice in shuffled order.

A review identified a real trap: filling a repeated node's entire own block
with its marginal variance skips the original off-diagonal loading product.
`additional-results.json` reproduces 16 differing entries, maximum 8.88e-16,
for that initial shortcut. The final prototype keeps the loading-product
arithmetic for off-diagonal repeats and assigns marginal variance only to the
true output diagonal. It matches exactly. Tiny 24-byte and 192-byte workspace
controls and the default 4 MiB control also match exactly.

The scoped pytest plugin forces the block algorithm even for tiny existing
fixtures, so the focused tests do not merely exercise the size fallback:

```text
python tools/check.py test -- tests/test_gaussian_tree.py tests/test_evolution.py \
  tests/test_numerical_invariance.py -p candidate02_plugin -rs
90 passed in 1.01s
```

This was run after the duplicate-node correction, from the immutable archive,
under the common exclusive lock, with the study directory on `PYTHONPATH`.
Ruff checks and formatting also pass for this candidate's Python files.
No full release gate, optional external backend, or replicated scientific
study was run for this isolated prototype. The scientific model, optimizer,
root treatment, likelihood, tolerances, bounds and seeds are unchanged.

## Reproduction and primary references

The script takes the same one-CPU measurement lock described in the shared
protocol. From the immutable archive's root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 NWKIT_CI_DIAGNOSTICS=0
export PYTHONPATH=/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/02-gaussian-covariance-blocks
$PY "$DIR/benchmark.py"
$PY "$DIR/additional_checks.py"
```

The actual runtime is CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
pandas 3.0.6 on macOS 27.0.1 ARM64. The original `.venv` is untouched.

- [Ho and Ané (2014), A Linear-Time Algorithm for Gaussian and Non-Gaussian
  Trait Evolution Models](https://doi.org/10.1093/sysbio/syu005), supports
  exploiting tree structure while preserving the model. The particular block
  implementation here is an inference from NWKIT's covariance identity, not
  a claim to implement that paper's linear-time likelihood algorithm.
- [Bastide et al. (2020), Efficient Bayesian Inference of General Gaussian
  Models on Large Phylogenetic Trees](https://arxiv.org/abs/2003.10336), provides
  primary context for general Gaussian tree structure and retained root/model
  semantics. It does not establish the measured speedup here.
- [NumPy `ix_`](https://numpy.org/doc/stable/reference/generated/numpy.ix_.html)
  documents open-mesh indexing for the Cartesian product of position groups.
- [NumPy indexing](https://numpy.org/doc/stable/user/basics.indexing.html)
  documents advanced-index assignments, copies and broadcasting; these inform
  the extra-memory accounting and tiled products.
