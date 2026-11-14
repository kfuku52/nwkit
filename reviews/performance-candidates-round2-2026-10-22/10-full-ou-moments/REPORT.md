# Candidate 10: consecutive full-OU moment-product reuse

**Verdict: adoptable within the tested native sorted-design diagnostic scope.**
The complete local-identifiability diagnostic improved **4.27x for four traits**
and **7.61x for eight traits**, including every finite-difference evaluation,
matrix exponential and final SVD. All moment entries, unscaled Jacobian entries
and returned diagnostic fields were exactly equal. Traced peak allocations
were essentially unchanged. This is not a measured speedup of the full OU
optimizer, a complete fit, a CLI invocation or population calibration.

Production code, existing tests, scientific settings and Git history were not
changed. The prototype and evidence are isolated in this directory.

## Mechanism and numerical boundary

The original `_moment_vector` computes
`(slopes[left] @ covariance @ slopes[right].T)[first, second]` independently
for every selected observed covariance coordinate. `observation_moment_design`
returns a sorted tuple, so different trait coordinates with the same two
distances occur consecutively. The prototype retains only the immediately
preceding distance pair and its complete matrix product, then extracts each
requested cell in the original design order.

Both matrix multiplications remain in their original order and shape. The
prototype does not substitute row dot products, reorder design entries,
diagonalize the attraction matrix, or introduce approximate transitions.
The parameterization, `decode_full_ou`, `scipy.linalg.expm`, distance set,
finite-difference step, Jacobian column scaling, rank tolerance, SVD and
complete/subset interpretation remain unchanged. Only one extra `d x d`
result is retained locally; it is replaced at the next pair and released
when the function returns. Computing the replacement also needs the normal
matrix-product temporaries, so this is an O(d²) cache bound rather than a
total-memory or RSS cap. No cache survives another parameter evaluation.

An externally supplied unsorted design is still evaluated in its own order.
It benefits only when an identical pair happens to be consecutive. Unique
pairs receive no saved multiplication, and their tuple comparison adds some
overhead. The native sorted-design scope is therefore important to the
recommendation.

## Comparable complete-diagnostic measurements

All times are median seconds per complete `full_ou_identifiability` call.
Input construction and `observation_moment_design` are excluded as fixture
setup; every diagnostic call includes its own parameter decoding, slopes,
product reuse, finite differences and SVD. This is the complete consuming
diagnostic API, not merely one `_moment_vector` kernel. Each case uses warmups,
seven alternating baseline/candidate samples, calibrated short calls, one
numeric thread and the study's exclusive measurement lock.

| Workload | Entries / consecutive pairs | Baseline | Prototype | Ratio |
| --- | ---: | ---: | ---: | ---: |
| Native non-ultrametric star, 2 traits / 12 tips | 219 / 55 | 0.007034 | 0.004198 | 1.68x |
| Native non-ultrametric star, 4 traits / 16 tips | 1,034 / 65 | 0.086061 | 0.020162 | **4.27x** |
| Native non-ultrametric star, 8 traits / 32 tips | 4,096 / 65 | 1.799186 | 0.236524 | **7.61x** |
| Native ultrametric deficient design, 2 traits | 7 / 2 | 0.000515 | 0.000424 | 1.21x |
| Externally shuffled design, 4 traits | 1,034 / 1,017 | 0.090354 | 0.090696 | 1.00x |
| Unique distance-pair control, 4 traits | 64 / 64 | 0.027466 | 0.028008 | 0.98x |
| Capped incomplete design, 2 traits | 2 / 1 | 0.000318 | 0.000301 | 1.06x |

The four-trait raw ranges were 0.08545–0.08822 versus 0.01945–0.02041 s;
the eight-trait ranges were 1.77475–1.84723 versus 0.23187–0.27685 s.
Every one of seven pairs improved in both principal cases. The two-trait
native case had substantial external-load noise and overlapping ranges, so
its larger median ratio is not a firm performance claim. The shuffled and
unique-pair controls also overlap; no improvement is claimed for either.
The unique-pair case's 0.98x ratio is retained rather than hidden.

The eight-trait design is capped/incomplete, but its full column rank is
sufficient for the unchanged `local_full_rank` diagnostic. This is not a
claim of global identifiability or precise estimation. Both four/eight-trait
fixtures have enough observed coordinates for the model's parameter-count
guard, but the outer fitted optimization was not benchmarked.

Traced whole-diagnostic peaks were 672,492→672,446 bytes for four traits and
9,913,662→9,913,294 bytes for eight traits. These negligible differences do
not establish a memory reduction; the retained Jacobian/SVD storage dominates.
Tracemalloc excludes some native allocations and total RSS. Raw timings,
loop counts, inputs, source hash and exact comparisons are in
[results.json](results.json).

## Correctness and verification

Every benchmark compares all moment values and all entries of the central-
difference Jacobian before column scaling, followed by exact comparison of
the status, rank and smallest/largest singular-value ratio. All comparisons
were exact. Seventeen additional controls cover ultrametric rotation,
near-ultrametric trees, complete and capped designs, shuffled/duplicate entries,
2/4/8 traits, time scales 1e-5 and 1e5, a nonnormal attraction with complex
eigenvalues and a negative diagonal, and an empty moment design.

The [additional audit](additional-results.json) checks the existing SVD rank
threshold more directly. With branch lengths near 1e-6, small differences
between tips give singular ratios on both sides of `1e-6`; baseline and
prototype retain exactly the same moment/Jacobian arrays and rank decisions.
Three failure controls also retain the exact exception type/message for an
out-of-range trait coordinate, too-short parameter vector and empty
identifiability design. No tolerance or failure check was relaxed.

The scoped pytest plugin replaces `_moment_vector` before collection while
retaining the original complete diagnostic and fit implementation:

```text
python tools/check.py test -- -p candidate10_plugin tests/test_full_ou.py -rs
9 passed in 6.63s
```

There were no skips. The selection includes the real free full-OU optimizer,
ultrametric/local-rank/subset classification, model-ranking exclusion,
independent dense GLS and original-unit references, CLI comparison and
posterior/bootstrap/predictive output checks. Fixed-matrix tests bypass
identifiability and are explicitly consumer controls rather than direct
exercises of the changed calculation. See [tests.log](tests.log).
Ruff lint/format and strict JSON parsing pass for the experimental files.

No full release gate, cross-platform BLAS comparison, complete-fit speed test,
large empirical dataset or replicated scientific study was run for this
isolated prototype. The unchanged normal affected checks and release gates
would still be required during production integration.

## Primary sources and reproduction

- [Bartoszek et al. (2012), multivariate adaptation](https://pubmed.ncbi.nlm.nih.gov/22940235/)
  is the primary multivariate-OU model context linked by NWKIT's guide. The
  [repository full-OU contract](../../../docs/guides/ASR.md#full-attraction-multivariate-ou)
  defines the particular parameterization and local-rank interpretation here.
- [NumPy `matmul`](https://numpy.org/doc/stable/reference/generated/numpy.matmul.html)
  documents the two-dimensional matrix multiplication and BLAS-backed `@`
  operation. Reusing an identical full product rather than selecting a
  differently associated row expression is the NWKIT-specific inference.
- [SciPy `expm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html)
  describes the matrix-exponential implementation and references Al-Mohy &
  Higham's scaling-and-squaring method. This prototype leaves that operation
  unchanged, including for nonnormal attraction matrices. The current web
  manual describes SciPy 1.18; the installed version measured here is 1.17.1.

The immutable baseline is `86c8f2c537769d852c17e2c789340f6c9de4ef54`
(NWKIT 0.43.55). Both result files record the imported snapshot path.
The validated environment is CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
pandas 3.0.6, macOS ARM64. Other user-owned CPU jobs remained active and were
not modified; these are exploratory host-local observations.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/10-full-ou-moments
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests
$PY "$DIR/additional_checks.py"
```

The scripts set thread limits before numerical imports, force snapshot imports,
and acquire the shared exclusive measurement lock internally. The shared
[protocol](../PROTOCOL.md) explains how to recreate the exact source archive.
