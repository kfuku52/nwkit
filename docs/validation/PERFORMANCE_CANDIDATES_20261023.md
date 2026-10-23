# Performance candidate review, 2026-10-23

Ten candidates were assigned to ten separate agents. Experiments started from
`9a229fc1c12d95828e9e905709fa1fc29ecb3435`; accepted changes are integrated into
the existing `master` checkout. Eight were accepted and two rejected. These are
local implementation optimizations, with the existing model settings, random
seeds and output contracts retained.

## Measurement protocol

Environment: macOS 27.0.1 ARM64, conda-forge Python 3.14.7, NumPy 2.5.3 and
SciPy 1.17.1. Interpreter startup, dependency metadata, compiled SciPy imports,
ETE tree parsing and the checkout version passed the development preflight.
The existing `.venv` could not start because of its CPU architecture and was
preserved.

This task serialized benchmarks and tests with an exclusive lock and set
BLAS/OpenMP thread counts to one. Unrelated jobs were also running on the host;
wall-time distributions were consequently noisy. CPU times, allocation counts
and alternating before/after runs support the decisions below. Ratios describe
the specified workloads, rather than general end-to-end command speedups.
Peak allocation means traced Python allocation during the operation, excluding
prebuilt fixtures; it is not total process RSS.

Each candidate used identical fixtures, a warmup and at least three measured
runs. Promising candidates received further equivalence comparisons, boundary
cases and focused regression checks. No scientific tolerances or search budgets
were relaxed.

## Decisions

| # | Candidate | Decision and evidence |
|---|---|---|
| 1 | Shared integer bitmasks for weighted tree distances | Rejected. Deep 1,200-tip trees used 97.6–97.9% less traced memory, but a 20,000-tip star increased from 20.34 to 62.07 MB and a balanced tree from 74.40 to 123.34 MB. |
| 2 | LCA depth alignment visits only set bits | Accepted. D+L cost CPU time decreased 8.9–11.9%; small Gaussian covariance consumers decreased 18–20%. Full reconciliation-table improvement was inconclusive. |
| 3 | Assemble consensus clades by current components | Accepted. At 1,024 tips, builder CPU speedups were 29.24× for nested clades, 49.98× for noncontiguous nested clades, 1.91× for balanced clades and 1.26× for a star. |
| 4 | Broadcast DTT denominators | Accepted for allocation reduction: 1.8–12.1% less peak allocation. CPU differences were −4.0% to +1.2%; no speedup claim. |
| 5 | Reduce Gaussian covariance temporary arrays | Accepted. Covariance construction CPU time decreased 22.0–45.3% across twelve cases; peak allocation decreased 9.5% in the balanced 512-tip signed-slope case. |
| 6 | Reuse discrete-ASR transition/log calculations | Rejected. Shared-matrix posterior CPU gains were only 0.7–1.9%, while a unique-matrix 32-state case regressed 8.5%. |
| 7 | Reuse monophyly group membership | Accepted. At 2,048 tips/two polyphyletic groups, handler CPU decreased 27.8% rooted and 7.3% unrooted. Small cases were inconclusive. |
| 8 | Remove redundant tree-to-table copies/sorting | Accepted. Median handler CPU reductions were 11.1%, 3.3% and 1.6% at 16, 256 and 4,096 tips; individual-condition effects remain noisy. |
| 9 | Reduce duplicate TSV parsing/copies | Accepted: remove only the redundant fresh-DataFrame copy. Peak allocation decreased 33.0–36.5% for large/wide inputs and 22–27% in the full tip-table consumer. |
| 10 | Reuse contrast metadata across traits | Accepted. At 512 tips/16 traits, CPU decreased 44.6% and peak allocation 17.3%; at 256 deep-tree tips/16 traits, CPU decreased 79.9% and peak allocation 58.4%. |

The rejected distance representation allocates an integer through its highest
set bit even for singleton clades. Large shallow trees therefore replace
linear singleton-set storage with quadratic bit storage. Its experiment passed
258 focused tests and preserved weighted RF exactly; that did not compensate
for the memory regression. No shape-dependent threshold was introduced.

LCA validation compared each implementation with a parent-walk oracle on
64,051 node pairs, including depth differences around powers of two. Query
indices, reconciliation CSV, D+L components and Gaussian covariance bytes
matched exactly. Alternating nine-run consumer comparisons supported the CPU
reductions; the complete reconciliation table remained dominated by other work.

Consensus construction skips whole previously assembled components instead of
rewriting an owner for every descendant tip at each ancestor. All 640 small
equivalence cases and eight benchmark snapshots preserved topology, child
order, branch lengths, support and properties exactly. Cases included partial
parents, noncontiguous masks, singleton/duplicate selections and missing
annotations. Balanced 1,024-tip peak allocation increased by about 18 KB (1.3%);
deep-tree allocation was essentially unchanged. A 4,096-tip baseline warmup
was stopped after exceeding the exploratory budget and supports no claim.

DTT replaces repeated denominator arrays with scalar/column broadcasting. It
preserves the dot reductions and row copies required by seeded replay, as well
as the original batching/workspace limits. All 87 additional comparisons
matched both the old implementation and the single-replica reference bit for
bit, including zero-length edges, unary/polytomous trees, 1–64 traits, extreme
scales and batch boundaries. DTT/audit tests passed 64 cases without skips.

Gaussian covariance uses equivalent broadcast integer indices instead of
repeated `np.ix_` calls and scales its loading product in place, retaining the
same operation order and tile dimensions. All twelve benchmark byte hashes
and 83 additional array comparisons matched exactly. These covered arbitrary
duplicate requests, small/default workspaces, root-prior variants and
overflow/underflow paths; both versions invoked the overflow fallback 490
times. Gaussian, numerical-invariance and branch-Gaussian CLI tests passed
102 cases without skips.

Discrete ASR already reuses transition matrices. A prototype additionally
cached logarithms only for matrix identities occurring more than once, to
avoid retaining dense log matrices for unique branches. Twenty-four timing
cases and 19 boundary comparisons retained exact results, but longer
five-by-five repeated measurements did not justify the added bookkeeping.
The rejected prototype received no additional pytest/static run.

Tree-to-table conversion builds rows in the insertion order of the branch-ID
dictionary, which already records level order, and avoids a second traversal
and redundant DataFrame sort/reset. Thirty-six timing conditions and 864
boundary comparisons preserved output bytes, errors and diagnostics exactly;
27 `table2nwk` roundtrips also matched. Focused table/rooting tests passed 101
cases without skips. Peak allocation decreased by only about 2.4–3.9 KB; this
is primarily a small-handler optimization, with no uniform whole-CLI speedup
claim.

Monophyly materializes the target-name set once per group instead of once per
foreign leaf. Fourteen benchmark conditions and 180 boundary comparisons
preserved TSV, messages and exceptions exactly; all nine existing command
tests passed. Eleven alternating rounds, each averaging three handler calls,
confirmed the large polyphyletic-group CPU reductions and found no material
regression for empty foreign-leaf sets. The large unrooted case used about
33 KB more peak allocation (0.83%), while rooted peak allocation was unchanged.

TSV loading retains raw-name reparsing, duplicate-header detection and blank-row
alignment validation. It simply assigns the preserved names to the already
new DataFrame. Thirty-six boundary inputs matched values, dtypes, serialized
bytes and exception types/messages exactly, including BOM, leading zeroes,
missing literals, multiline quoting and malformed input. Stdin and independent
subsequent reads also matched. The interface/ASR-review tests passed 52 cases
without skips. Raw-table CPU medians improved modestly, but host variability
and the initially reversed large-consumer measurement make memory reduction
the primary justification.

Contrast tables reuse first-trait rows only when a second trait actually needs
the same metadata/CSV strings. Every trait replaces its own label and all four
numerical contrast fields; coefficient and sampling-covariance work remains
trait-specific. Nine benchmark workloads and 52 additional boundary comparisons
matched independent single-trait results exactly, including model transforms,
orientation, filtering, quoting and dense/structured/sparse sampling. The agent's
contrast/replicate/evolution tests passed 98 cases without skips. Single-trait,
empty-trait and fully excluded cases do not allocate the cache. A fully
excluded control showed +0.581 ms CPU (+6.9%), with overlapping measurement
ranges; no speedup is claimed for that case.

## Verification and reproduction

The seven initially accepted changes passed this integrated check:

```sh
python tools/check.py quick -- \
  tests/test_clade_index.py tests/test_consensus.py \
  tests/test_gaussian_tree.py tests/test_numerical_invariance.py \
  tests/test_branch_gaussian_cli.py tests/test_dtt.py tests/test_dtt_audit.py \
  tests/test_monophyly.py tests/test_nwk2table.py tests/test_table2nwk.py \
  tests/test_rooting_state.py tests/test_tree_outputs.py \
  tests/test_util_tree_io.py tests/test_util_tree_helpers.py \
  tests/test_asr_review_regressions.py tests/test_cli.py \
  tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_reconcile.py -rs
```

Result: Ruff lint/format and mypy (260 source files) passed; **582 tests passed**
in 25.75 seconds, with no skips or failures. After integrating the last
candidate, the following checks also passed:

```sh
python tools/check.py quick -- \
  tests/test_contrast.py tests/test_replicates.py tests/test_evolution.py \
  tests/test_regress.py::test_pgls_raw_mode_supports_mixed_response_replication_depth -rs
python tools/check.py test -- \
  tests/test_regress.py::test_pgls_raw_known_se_supports_multiple_responses_and_predictors -rs
python tools/check_maintainability.py
```

The final quick check passed Ruff lint/format, mypy and **99 tests** in 10.75
seconds. Its initial mypy run requested an explicit row-list type, which was
added before the successful rerun. The explicitly selected slow/integration
regression passed in 29.34 seconds. Both test runs had no skips or failures.

Maintainability passed all hard limits. Radon 6.0.1 was installed in a fresh
temporary tooling environment using `constraints-dev.txt`, because it was
absent from the numerical environment. Against the actual source baseline,
LCA complexity decreased 7 to 6 and contrast-table complexity increased 30 to
34, below the limit of 40. Other reported growth warnings predate this change;
the shared baseline and exception limits were preserved.

Full/release checks, replicated scientific studies, other Python versions and
other operating systems were not run for this local commit. The earlier LCA
quick run deselected four slow GeneRax integrations and one free full-OU
optimizer case. Focused successes do not imply those backends or the complete
delivery gate were validated.

Individual exploration scripts, patches, raw JSON and reports are supplied
together in the local `performance-candidates-20261023.tar.gz` evidence archive.
Its README explains reconstruction of baseline and candidate source trees;
the exploratory outputs are not bundled into the installed package.

The existing shared core harness also exercises consensus construction and
unaffected tree/Gaussian operations on either source checkout:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python tools/benchmark_core.py --checkout /path/to/checkout --output /tmp/core.json
```

The baseline/integrated core runs matched RF and consensus invariants, sister
IDs, interval cardinalities, and Gaussian regression objective, coefficients
and coefficient covariance exactly. Their noisy wall-time differences are not
used as performance claims.
