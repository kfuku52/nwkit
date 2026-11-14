# Candidate 07: cache logarithms of shared ASR transition arrays

**Decision: conditional / defer adoption.** The prototype removes redundant
matrix logarithms and preserves tested numerical behavior exactly, but the
complete fixed-generator marginal calculation has no demonstrated speedup.
Do not adopt this optimization on a performance claim from this experiment.
An unconditional cache also increases memory on unique-length trees. No
production files were edited.

## Algorithm and literature

The baseline already computes transition matrices once per distinct branch
length for a common generator. Its symmetric-generator eigendecomposition
and ER closed-form paths are existing optimizations, not this proposal.
Inside `_compute_log_inside_likelihoods`, however, the fixed-transition branch
calls `_log_probabilities(matrix)` again for each edge even when several edges
reference the same NumPy array. The regime-specific path already caches both
the matrix and its logarithm by generator identity and branch length.

The proposed change memoizes these logarithms by array identity within one
likelihood call. It retains the original `logsumexp`, scaling, traversal,
probability validation, child terms and returned transition arrays. The
shared-only variant first counts array identities and retains logs only for
arrays referenced more than once. Its counter pass costs O(edges), but avoids
retaining O(edges × states²) logarithms when all matrices are distinct.

Relevant primary and official sources:

- [Felsenstein (1981), likelihood pruning](https://doi.org/10.1007/BF01734359),
  [primary paper PDF](https://gen-files.princeton.edu/sysbio/Felsenstein81.pdf):
  the child-to-parent likelihood recursion. The recurrence and stochastic
  process are unchanged here; caching is an inference from NWKIT's source.
- [SciPy 1.17 `logsumexp`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.special.logsumexp.html):
  stable log-domain summation is preserved rather than replaced by an
  unscaled probability product or `logaddexp.reduce`.
- [NumPy `log`](https://numpy.org/doc/stable/reference/generated/numpy.log.html):
  an elementwise operation returning a freshly allocated array by default;
  zero probabilities yield negative infinity, which must retain its meaning.
- [Python 3.12 `id`](https://docs.python.org/3.12/library/functions.html#id):
  identities are unique and stable during an object's lifetime. The fixed
  transition dictionary retains its arrays for the duration of the call.
- Repository [ASR guide](../../../docs/guides/ASR.md): root/polytomy,
  nonnegative-length, model and fixed-generator contracts.

No cache of numerical values survives a call. Mutating transition arrays
between successive calls therefore recomputes logs. Concurrent mutation
during one likelihood calculation was not tested. Identity caching does not
coalesce different array objects that happen to contain equal values.

## Measurements

Baseline commit: `d7116f80b9979739c69a5431a5e34b822795408b`, NWKIT 0.43.54.
Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
ETE 4.4.0 rebuilt at the same version. The parent's documented environment
preflight passed. One BLAS/OpenMP thread is used. Fixtures are seeded balanced
or star trees with symmetric, asymmetric or disconnected generators. The
repeated-length fixture uses length 0.2 on every edge; unique lengths are
0.02 + edge_index × 0.001. The full results include every raw timing.

Both implementations have a warmup followed by seven paired measurements,
with order reversed on alternating repetitions. Short helpers are batched
to target at least 40 ms for the fastest warmup estimate; larger calls use
one invocation per sample. Prototype compilation is outside measured calls.
The shared advisory lock serializes this study, but unrelated user-owned jobs
were active on the same host. All ratios below are exploratory. Their ranges
overlap; small effects cannot be separated reliably from host contention.

Inside-helper times include the complete pruning call with already computed
transition matrices, excluding tree/Q/transition construction. Selected
shared-only measurements are:

| Tips / states | Lengths / structure | Baseline median ms | Candidate median ms | Ratio |
| --- | --- | ---: | ---: | ---: |
| 128 / 4 | repeated, balanced | 10.098 | 10.252 | 0.985× |
| 512 / 4 | repeated, balanced | 43.190 | 40.778 | 1.059× |
| 256 / 16 | repeated, balanced | 26.950 | 21.937 | 1.229× |
| 128 / 32 | repeated, balanced | 13.562 | 12.803 | 1.059× |
| 256 / 8 | unique, balanced | 20.844 | 20.897 | 0.997× |
| 128 / 4 | repeated, star | 5.859 | 5.049 | 1.160× |
| 128 / 4 | regime path, already cached | 10.331 | 13.524 | 0.764× |

The 256/16 helper result has baseline 21.550–49.916 ms and candidate
20.986–48.166 ms ranges. Its 1.229× median is not a reliable independent
performance claim. The regime path removes no existing computation; its
apparent slowdown also has overlapping ranges and substantial host noise.

The higher-level measurement calls the complete `compute_mk_marginals` with
CUSTOM fixed Q, including generator validation, likelihood, inside/outside
passes, posterior normalization and fit metadata. Imports and fixture creation
are excluded. These are full API calls, not full CLI runs:

| Tips / states | Lengths | Baseline median ms | Candidate median ms | Ratio |
| --- | --- | ---: | ---: | ---: |
| 128 / 4 | repeated | 48.529 | 50.636 | 0.958× |
| 128 / 16 | repeated | 54.830 | 55.439 | 0.989× |
| 128 / 8 | unique | 60.350 | 60.605 | 0.996× |

All three complete-call timing ranges overlap. Reducing the count of logs
does not establish a whole-operation benefit: the stable reductions, outside
pass, posterior work and validation remain. An optimizing Mk fit was not
timed, so no fit-search or CLI throughput claim follows.

The mechanism itself is confirmed: on the 128-tip repeated-length tree,
matrix-log calls fall from **254 to 1**; leaf-vector logs remain 128. At
512 tips they fall from **1022 to 1**. Regime and unique-length cases retain
their original counts (2 and 510 matrix logs respectively).

Peak traced Python allocations for the 256-tip/8-state unique-length inside
call are baseline **397,589 bytes**, unconditional cache **761,582 bytes**,
shared-only **401,539 bytes**. The guard avoids this retention cost. Other
shared-only peaks and complete-call peaks are essentially unchanged. These
are Python allocation peaks, not RSS or complete native-memory accounting.

## Equivalence and verification

- Every nested inside array, scale, child-term array and transition output
  matched exactly in all nine main cases. Symmetric/asymmetric, disconnected,
  star, regime and unique-length structures are included.
- Every posterior and every returned fixed-Q fit field matched exactly in all
  three complete-call cases; no tolerances were relaxed.
- Additional zero-length/identity transitions, zero-probability impossible
  observations and negative-infinity propagation matched exactly for three
  generator structures. Arrays mutated between calls also matched baseline,
  demonstrating that the prototype does not reuse stale logs across calls.
- With the shared-only prototype installed through `candidate07_plugin.py`,
  existing ASR, polytomy, regime, and discrete-model tests ran via
  `tools/check.py test`: **112 passed in 17.37 seconds**, no skips.
  This includes independent transition/reference comparisons, impossible
  state handling, stochastic maps/seeds, fixed models, free-rate models and
  CLI consumers. These free-rate tests establish regression behavior in
  their fixtures; they are not timing evidence for optimizing fits.
- Ruff lint and format checks passed on the two experiment scripts.

## Reproduce and adoption boundary

From the repository root:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-candidates-2026-10-22/07-asr-log-transition-cache/benchmark.py
$PY reviews/performance-candidates-2026-10-22/07-asr-log-transition-cache/benchmark.py --tests
```

The script sets the four thread limits to one and holds the common exclusive
`fcntl.flock` for both benchmarking and tests. `--tests` invokes:

```sh
python tools/check.py test -- tests/test_asr.py tests/test_asr_polytomy.py \
  tests/test_asr_regimes.py tests/test_discrete_asr_models.py -rs
```

The subprocess receives `PYTEST_PLUGINS=candidate07_plugin` and this directory
on `PYTHONPATH`; [prototype-tests.log](prototype-tests.log) records the command
and results. [results.json](results.json) contains all raw samples, memory
readings, logarithm call counts and edge-check records.

Current recommendation is to defer this change. Consider it again only if
an otherwise idle-host representative workload shows a consistent complete-fit
benefit large enough to justify the additional cache/guard complexity.
Unconditional retention should remain rejected on unique-length workloads.
Cross-platform measurements, full/release checks and large scientific studies
were not run because this is an isolated experiment with unchanged models and
production sources. They would be required as appropriate during integration.
