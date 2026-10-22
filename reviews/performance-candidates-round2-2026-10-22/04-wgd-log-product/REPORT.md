# Candidate 04: bounded row batches for WGD log matrix products

**Verdict: conditional; not a priority for production adoption.** The matrix
kernel improves, but the complete consuming likelihood improves only modestly
in the stronger repeat, while requiring about one additional MiB of traced
allocation. No production code, scientific parameter, test tolerance, optimizer
bound or random seed was changed. Candidate 03's family batching is not active.

## Proposal and model boundary

`_log_matrix_product(left, right)` currently calls `logsumexp` once per output
row. The prototype broadcasts a bounded number of independent rows and reduces
the same latent-state axis in the same order. It retains the original function
when the estimated row batch is one or the input has at most one row.
The selector uses a 1 MiB working-array estimate and at most 32 rows. This is
an estimate, **not a strict peak-memory or RSS limit**: SciPy's stabilization
and array bookkeeping create additional bounded temporaries.

The only patched function is the log-semiring product used to compose the
before-event birth/death transition, retained-copy pulse, and after-event
transition. Unbounded branching probabilities, omitted state/root tails,
ascertainment, missing-data integration, detection, rate-category mixing,
adaptive state-convergence checks and zero-retention split-edge checks retain
their baseline implementations. The ordinary probability path is unchanged.

Primary sources and official documentation:

- [SciPy `logsumexp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html)
  documents multi-axis array reductions and stable log summation. Batching
  independent rows is an inference from this API and the source, not a documented
  NWKIT speed guarantee.
- [Blanchard, Higham & Higham, log-sum-exp accuracy](https://arxiv.org/abs/1909.03469)
  explains why algebraic rearrangements require floating-point validation.
  The prototype keeps SciPy's stabilized operation instead of substituting a
  probability-domain product that could underflow.
- [Rabier, Ta & Ané (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4038794/)
  provides the birth/death and retained-WGD model context. The current
  [NWKIT guide](../../../docs/guides/WGD_COUNT.md) remains authoritative for the
  more limited native command and interpretation.

The current SciPy web manual describes 1.18.0; actual measurements used 1.17.1.

## Evidence

Both implementations were warmed up, measured in alternating order seven
times, and compared on identical inputs. Short kernels were calibrated into
approximately 30 ms batches. The full-operation initial run experienced large
host-load changes, so one limited follow-up used approximately 0.3 s samples.
All raw times are preserved in [results.json](results.json) and
[followup-results.json](followup-results.json). Other agents' numerical work was
serialized through the common lock; unrelated user-owned CPU jobs remained
active.

| Operation | Baseline median | Prototype median | Ratio | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Kernel, count bound 20 | 1.059 ms | 0.141 ms | 7.51x | Prepared input kernel only |
| Kernel, count bound 40 | 2.860 ms | 0.844 ms | 3.39x | Prepared input kernel only |
| Kernel, count bound 86 | 9.729 ms | 5.787 ms | 1.68x | Prepared input kernel only |
| Complete rare-family likelihood, 1 family/bound 86, longer samples | 81.962 ms | 69.627 ms | 1.18x | Ranges overlap; one of seven pairs worsened |
| Complete rare-family likelihood, 16 families/bound 86, longer samples | 80.163 ms | 74.855 ms | 1.07x | Seven of seven pairs improve, ranges narrowly separate |
| Complete log-transition preparation, bound 40, longer samples | 37.166 ms | 35.654 ms | 1.04x | Large range overlap; inconclusive |

For the stronger 16-family comparison the baseline range was 78.356–86.642 ms
and the prototype range 72.611–77.619 ms. This establishes a small local effect
for that input, not a stable release-level gain. Ordinary-family and zero-pulse
controls make zero calls to the patched product and have no claimed speedup.
At bound 128 the workspace guard selects the original row loop; no gain is
claimed there. Initial large apparent full-operation gains or regressions are
not adoption evidence because their raw values overlap widely.

All 14 initial operation comparisons, four follow-up operation comparisons,
and 26 additional helper checks matched **exactly**, including dtype and
nonfinite masks. Helper checks include float32/float64, single/empty input rows,
rectangular products, noncontiguous arrays, impossible rows/columns, extremely
negative logs, retention 0/0.4/1, multiplicity 2/3, and tiny/default workspace
selectors. Finite numerical differences were zero. These are tested cases,
not a proof for every array or floating-point implementation.

The measured peak traced allocation for the complete bound-86 rare likelihood
was 1.200→2.201 MB for one family and 1.234→2.231 MB for 16 families (decimal MB).
The bound-40 kernel was 0.121→1.608 MB. These are separate `tracemalloc` runs,
not total RSS; inputs are prepared before tracing. This extra allocation is
material relative to a roughly 7% complete-operation improvement.

The scoped prototype was installed in the immutable baseline for existing
`tools/check.py test` regressions: **47 passed, 35 deselected**, followed by the
specific zero-retention split-edge convergence test: **1 passed**. These cover
the independent high-copy Yule pulse sum, log/probability pruning agreement,
missingness, pulse limits, unnormalized tails, rate mixture/selection, invalid
event nodes and convergence. See [tests.log](tests.log) and
[split-edge-test.log](split-edge-test.log). Ruff lint and formatting passed for
the experimental Python files.

## Reproduction and limits

Source baseline is `86c8f2c537769d852c17e2c789340f6c9de4ef54`, NWKIT 0.43.55.
Imports were asserted against the immutable snapshot at
`/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline`; the live workspace had unrelated
changes from another task. Environment: CPython 3.12.14, NumPy 2.5.4,
SciPy 1.17.1, pandas 3.0.6, macOS 27.0.1 ARM64. The metadata and baseline
function hash are retained in the JSON. The original `.venv` was preserved.

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
STUDY=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/04-wgd-log-product
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export NWKIT_CI_DIAGNOSTICS=0 PYTHONPATH="$PWD"
$PY "$STUDY/benchmark.py"
$PY "$STUDY/followup.py"
$PY "$STUDY/run_checks.py"
$PY "$STUDY/run_checks.py" \
  tests/test_wgd_count_model.py::test_zero_retention_state_check_cannot_hide_split_edge_truncation -rs
```

Each script acquires the shared advisory measurement lock. No full optimization
or scan/bootstrap timing was run: those repeat many ordinary-probability paths
and cannot inherit a kernel speedup by assumption. No scientific studies or
optional external backends were run because this prototype changes only exact
array batching. Integration would require confirming worthwhile whole-search
benefit on a quiet host and a workspace policy that accounts for SciPy's actual
temporaries. The present evidence supports retaining the prototype for review,
not prioritizing another production path.
