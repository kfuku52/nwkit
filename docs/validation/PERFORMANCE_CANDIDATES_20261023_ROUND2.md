# Further performance candidates, 2026-10-23

Six further candidates were investigated against
`efee7cb2098fe87de5c7ea3270ac6a6e4e023316`, after the
[first review](PERFORMANCE_CANDIDATES_20261023.md). Five were retained; the
clade-ID encoding cache was rejected after measuring its actual consumers.

## Protocol and decisions

The development preflight passed with conda-forge Python 3.14.7, NumPy 2.5.3
and SciPy 1.17.1 on macOS 27.0.1 ARM64. Benchmarks used one BLAS/OpenMP thread,
prebuilt identical fixtures, a warmup, alternating before/after order, and
five or seven measured repetitions. The final cladefreq microbenchmark used
nine repetitions. Process CPU and wall time were measured separately from
tracemalloc peak allocation; peaks exclude prebuilt fixtures and are not RSS.
Reported speedups apply to the named operations and inputs.

| # | Candidate | Decision and measured result |
|---|---|---|
| 1 | Read descendant leaf-interval bounds from the first/last child | Retained. At 1,024 balanced tips, helper CPU decreased 1.831 to 1.375 ms; the real edge-vector consumer improved 1.07–1.12× in the tested trees. |
| 2 | Enumerate sparse cladefreq masks without scanning every unset bit | Retained after revising the first prototype. Final handler CPU improved 1.20–2.48× across balanced/ladder trees with 256/1,024 tips. |
| 3 | Cache UTF-8 length-prefixed clade-ID fields | Rejected. Hashing all IDs alone improved 1.56–2.45× on some shapes, but reconciliation/contrast consumers were essentially unchanged or slightly slower, with higher memory use. |
| 4 | Reuse descendant-label strings across annotation aggregations | Retained. The 256-tip ladder/12-aggregation helper improved 9.78× CPU and peak allocation decreased 3,984,789 to 1,336,631 bytes. Real annotate handlers improved 1.31–2.91×. |
| 5 | Derive union cardinality from set sizes and intersection size | Retained. With precomputed species sets, a 512-tip distinct-species ladder score pass improved 4.09× CPU; peak allocation decreased 34,256 to 1,816 bytes. |
| 6 | Classify diagonal sampling covariance without a dense comparison matrix | Retained for allocation reduction. At 2,048 dimensions, the predicate peak decreased 37,749,697 to 3,929 bytes; diagonal-case CPU decreased 2.649 to 1.739 ms. |

The leaf-interval optimization relies on DFS child order: the first child's
start and last child's end bound the parent's contiguous leaf interval.
Duplicate-label inputs retain the original explicit-set implementation.
There were 2,189 node comparisons against direct descendant walks, including
singletons, unary nodes, polytomies, reordered children and duplicate names.
Membership and iteration order also matched where interval views apply.

The first cladefreq prototype visited only set bits but regressed on dense
large masks. The retained implementation uses a direct scan at densities of
at least one quarter, and otherwise shifts past consumed positions. This
changes iteration only; displayed clades, weights and ordering are unchanged.
The 24-case density sweep covered 128/512/2,048 bits and 1–100% density,
including 24/25/26% around the switch. CPU improved 1.14–26.0× in that sweep.
Intermediate sparse cases used up to about 22 KB more traced allocation;
the four real-handler peak measurements were essentially unchanged. Exhaustive
small-mask and invalid-mask comparisons retained results and exceptions.

The rejected encoding cache preserved every tested SHA256 identifier, but
added approximately 9–37 KB in actual consumers. Their seven-round CPU
comparisons ranged from approximately 2% faster to 6% slower. That did not
justify a persistent cache. The original stable-ID encoding remains intact.

Annotation computes each internal node's label string lazily for its first
report row, then reuses it for subsequent aggregations on that node. No report
or no aggregation means no label formatting. All eight helper cases and the
two real handlers matched exactly. Another 54 command cases covered singleton,
unary and polytomous trees, missing/nonnumeric values, absent reports, empty
aggregations, exceptions and previous output preservation. The one-aggregation
helper control had noisy timing differences and no allocation improvement;
the optimization targets repeated report metadata.

Duplication confidence still uses the identical integer numerator/denominator:
`intersection / (len(first) + len(second) - intersection)`. Empty unions and
nonbinary nodes retain zero scores. All 72 additional score comparisons were
exact, using both calculated and supplied species sets. The large speedup is
for the score pass with precomputed sets, as supplied by `subtree`; it does not
include species-set construction or promise the same whole-command speedup.

The sampling predicate uses the existing exact off-diagonal check and retains
`array_equal`'s rejection of NaNs on the diagonal. There were 720 reference
comparisons, including zero/rectangular sizes, tiny nonzero off-diagonals,
signed zero, infinity, NaNs, numeric dtypes and array layouts. At 1,024
dimensions, C/F/transposed/strided cases used about 2.9 KB instead of
9.4–9.5 MB. The 32-dimensional timing case increased roughly 39 to 54 microseconds;
the main justification is removing quadratic temporary storage, while keeping
the same numerical branches and failure checks.

## Integrated verification

The following check passed Ruff lint/format, mypy over 260 source files, and
654 tests in 71.20 seconds, with four slow regression cases deselected:

```sh
python tools/check.py quick -- \
  tests/test_util_tree_helpers.py tests/test_cladefreq.py tests/test_annotate.py \
  tests/test_subtree.py tests/test_dist.py tests/test_diff.py tests/test_skim.py \
  tests/test_transfer.py tests/test_mapping_tables.py tests/test_clade_mapping.py \
  tests/test_rooting_state.py tests/test_util_tree_io.py tests/test_tree_outputs.py \
  tests/test_regress.py tests/test_gaussian.py tests/test_numerical_invariance.py \
  tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_output_transaction.py -rs
```

The four deselected tests were then selected explicitly:

```sh
python tools/check.py test -- tests/test_regress.py -m slow -rs
```

All four passed in 64.16 seconds: large factor-loading predictor uncertainty,
both bootstrap variants refitting the automatic gene parameter, and the
multiple-response/predictor known-SE workflow. The 106 non-slow cases in that
file were deselected in this second run. Together the two runs passed 658
tests without skips or failures.

`python tools/check_maintainability.py` also passed all hard limits using the
existing isolated Radon 6.0.1 tooling environment. Against the actual baseline,
leaf-interval complexity decreased 13 to 10, annotation increased 11 to 12,
cladefreq increased 3 to 6, and the new sampling predicate scored 2. The shared
comparison baseline and exception limits were preserved.

The existing `tools/benchmark_core.py` was also run against both source trees
with one BLAS/OpenMP thread. Tree invariants and Gaussian regression objective,
coefficients and coefficient covariance matched exactly. Its unaffected
workload timings are recorded without claiming additional speedups.

Full/release checks, replicated studies, other Python versions and other
operating systems are outside this local validation run.

The local `performance-candidates-round2-20261023.tar.gz` evidence archive
contains replayable benchmarks, source patches, JSON results and logs. Its
README records how to reconstruct the archived baseline and individual
candidate modules. Exploratory outputs stay outside the repository.
