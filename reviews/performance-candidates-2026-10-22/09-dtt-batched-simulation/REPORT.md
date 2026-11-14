# 09: Bounded replica batches for Brownian DTT simulations

**Decision: adoptable within the tested low-dimensional simulation scope, with
an explicit memory tradeoff.** The compatible prototype improves the complete
serial simulation operation and passes the existing exact seed/thread contracts.
Integration should retain a small batch bound and confirm representative larger
trait dimensions and platforms. No production sources or tests were edited.

## Rationale and sources

[Harmon et al. (2003), Tempo and mode of evolutionary radiation in iguanian
lizards](https://doi.org/10.1126/science.1084786) describe disparity-through-time
comparisons against Brownian simulations. That motivates repeated simulation
work, not a particular optimization or performance ratio. NWKIT's authoritative
[DTT guide](../../../docs/guides/DTT.md) specifies its squared-distance disparity,
stable clade moments, conditional simulation envelope, indexed random streams,
and work limits. All are retained here; no simpler disparity or independent-trait
approximation replaces the existing joint Brownian model.

[NumPy broadcasting](https://numpy.org/doc/stable/user/basics.broadcasting.html)
supports performing the same deterministic elementwise operations across replica
axes. [NumPy's SeedSequence documentation](https://numpy.org/doc/stable/reference/random/bit_generators/generated/numpy.random.SeedSequence.html)
and [parallel RNG documentation](https://numpy.org/doc/stable/reference/random/parallel.html)
describe spawned independent streams. Each replica therefore still creates the
same `default_rng(seed)` and draws the same `(nodes-1, traits)` normal matrix;
its multiplication by the same rate factor is also kept as an individual call.
The [matmul documentation](https://numpy.org/doc/stable/reference/generated/numpy.matmul.html)
explains batched products, but bitwise-equivalence measurements determined which
operations were safe to batch here.

At baseline `d7116f80b9979739c69a5431a5e34b822795408b`, `_simulate_chunk`
calls `_simulate_one` for every replica, repeating node traversal, square roots,
count bookkeeping and stable subtree-moment setup. The prototype replaces only
`_simulate_chunk`. It processes at most 32 replicas together, computes fixed
counts/square roots once per chunk, and propagates values and merges stable
means/sums across a replica axis. It preserves child order. Singleton/empty
chunks use the original per-replica helper. The existing `simulate_curves`
dispatch, spawned process pool, indexed chunks, output sorting, model fitting,
work guards and public output schema remain unchanged.

## Exactness required a restricted batching strategy

The first fully batched prototype changed 13 of 165 values in a 5-trait
polytomy check, by at most **2.22e-16**. Weakening an `array_equal` check would
have broken the existing thread/prefix contract. The failed checks are retained
in `checks-initial.log`, `checks-sparse-matvec.log`, `checks-scalar-dot.log` and
`checks-array-division.log`.

Diagnostics found branch values and subtree means identical while some subtree
sums differed (`rounding-diagnostic.json`). A fresh contiguous copy of each
replica's small delta vector before the original `np.dot(row, row)` removed the
observed discrepancy. Offset rows of a batch buffer and fresh vectors can follow
different low-level dot reductions; the exact backend cause was not established.
The final prototype conservatively retains these per-replica dot calls, the
original sparse matrix-vector curve projection, and explicit array denominators.
It does not claim that arbitrary batched reductions are bitwise interchangeable.
The compatibility copies have a cost, included in every timing below.

## Comparable whole-operation measurements

Environment: macOS 27.0.1 ARM64, CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
NWKIT 0.43.54. ETE 4.4.0 uses the same-version source build after the original
wheel import failed; the original `.venv` was preserved. BLAS/OpenMP limits are
one. The common exclusive lock serialized our candidate jobs. Other user-owned
CPU workloads remained active outside this lock, so all timing magnitudes are
exploratory and need quiet-host confirmation before release claims.

Timings include the complete serial `simulate_curves` call: seed spawning,
generators, normal draws, rate-factor products, local batch setup, propagation,
disparity computation, indexed collection, sorting and final array construction.
Imports, fixture construction and prototype installation are excluded. Fixtures
are rooted ultrametric balanced, ladder, star and three-way polytomy trees,
with the same correlated, seeded rate factor and seed 42 for both versions.
There is a warmup per implementation followed by seven paired trials, reversing
order on alternate trials. Short operations are calibrated into batches of at
least approximately 25 ms. Raw timings and full output comparisons are in
[results.json](results.json).

| Shape | Tips | Traits | Draws | Baseline median | Candidate median | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Balanced, small control | 8 | 2 | 4 | 0.320 ms | 0.208 ms | 1.53x |
| Balanced | 64 | 2 | 128 | 66.57 ms | 16.04 ms | 4.15x |
| Balanced | 256 | 3 | 128 | 257.85 ms | 59.11 ms | 4.36x |
| Ladder | 128 | 1 | 96 | 119.54 ms | 25.40 ms | 4.71x |
| Star | 128 | 2 | 96 | 64.26 ms | 12.93 ms | 4.97x |
| Polytomy | 81 | 4 | 96 | 48.24 ms | 11.44 ms | 4.22x |
| Balanced, single-draw control | 64 | 2 | 1 | 0.543 ms | 0.522 ms | 1.04x, inconclusive |

For each multi-draw case the slowest candidate trial was faster than the fastest
baseline trial. Minimum paired ratios for the five larger cases ranged from
2.93x to 4.20x. Their baseline relative timing ranges were 20–41% and candidate
ranges 16–51%; the separation supports the direction, not an exact stable ratio.
Single-draw behavior delegates to the original kernel and shows no established
benefit. These are full simulation measurements, not full CLI/model-fit timings.
The star case is a throughput/control workload; its timing is biologically
uninformative according to the existing DTT guide.

## Memory tradeoff

Temporary arrays scale as `O(batch * nodes * traits)` instead of a single
replica. With 32-replica batches, measured peak traced Python allocation rose
from **134,003 bytes to 1,977,778 bytes** in the 256-tip, 3-trait, 128-draw case
(about 14.8x, or 1.84 MiB extra). Other larger cases increased approximately
3.1–9.0x. Input arrays are excluded; `tracemalloc` is not total process RSS
and may miss native allocations. The 32 bound limits this extra work space,
but no claim is made about every supported 64-trait input.

A separately paired seven-trial **8-replica** control of that same 256-tip case
retained exact curves and measured **2.53x**, with peak allocation
**134,196 → 570,044 bytes**. See [results-batch8.json](results-batch8.json).
This demonstrates a useful speed/memory choice, not that 32 is optimal.
An integrated implementation should choose/document a conservative memory
budget rather than allocating all requested simulations together.

## Numerical and contract checks

- Every timed curve is **array-equal** to baseline. Means, medians, 2.5/97.5%
  envelopes, complete areas and interpolated partial areas also match exactly.
- Six direct cases cover balanced/ladder/star/polytomy shapes, 1–5 traits,
  4–65 replicas and batch boundaries. The prototype's serial and two-worker
  outputs match exactly; one- and four-draw prefixes match longer runs, and
  selected curves match the untouched `_simulate_one` oracle exactly.
- Additional checks cover a valid tree with zero edges, a unary node and
  polytomies; reversed tip indexing; non-contiguous factors; and factors scaled
  by `1e-140` and `1e140`. All curves and downstream summaries remain identical.
- Generator probes for 35 replicas confirmed the same normal-draw shape and
  exact post-draw bit-generator state for every replica. This checks the actual
  draw sequence, not merely statistical similarity.
- Zero/single draws match baseline. Constant/nonfinite factors and a negative
  seed preserve exception type and message. Public validation and work limits
  are not replaced or relaxed. Additional evidence is in
  [direct-checks.json](direct-checks.json).
- With the exact prototype installed through `candidate09_plugin.py`,
  `tools/check.py test -- tests/test_dtt.py tests/test_dtt_audit.py -rs` completed
  **51 passed in 6.91 s, no skips**. This includes the independent geiger and
  pairwise/GLS references, dense BM distribution comparison, exact two-worker
  seed behavior, CLI tables/figures, missingness, units, unusual topologies and
  transactional failures. Full command/output: [tests.log](tests.log).
- Ruff lint and formatting checks passed for the experimental Python files.

## Reproduction and remaining scope

From the repository root with the temporary environment:

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
  $PY reviews/performance-candidates-2026-10-22/09-dtt-batched-simulation/benchmark.py
$PY reviews/performance-candidates-2026-10-22/09-dtt-batched-simulation/benchmark.py --checks-only
$PY reviews/performance-candidates-2026-10-22/09-dtt-batched-simulation/benchmark.py --tests
```

The script sets the same thread limits and acquires the shared lock internally.
For the batch-8 control, import the benchmark module in a lock-holding wrapper,
set `BATCH_SIZE = 8`, and invoke `benchmark("balanced", 256, 3, 128)`; the output
record contains all parameters and seven raw baseline/candidate measurements.

Cross-platform bitwise identity, high-dimensional traits, throughput with more
than two processes, peak RSS, and full CLI wall time remain unmeasured. No live R
run, replicated calibration study, full/release gate or distribution build was
run because this is an isolated implementation experiment, without a shipped
model/interface change. Production integration still requires its normal static
checks and focused regression coverage, plus quiet-host performance confirmation.
