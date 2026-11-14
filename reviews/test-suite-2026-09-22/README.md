# Test suite review, 2026-09-22

Baseline: `8dddd3e`, 191 test modules and 4,140 collected cases. Review used the
suite inventory, assertion/call inspection, duplicate-body checks, and comparison
with the exercised product paths. Test count and coverage were not retention goals.

## Removed or consolidated

| Suites (`tests/test_*.py`) | Decision and remaining failure detection |
| --- | --- |
| `label`, `rescale`, `printlabel` | Remove weaker repetitions of exact-output tests. Keep target selection, forced replacement of existing names, uniqueness, lengths, and invalid inputs. |
| `prune`, `subtree` | Combine retained-tip assertions with branch-length checks. Keep root promotion, partial selection, whitespace parsing, metadata and failure cases. |
| `sanitize`, `util_tree_helpers` | Remove private quote-state tests and repeated singleton examples. Keep serialized quote round trips, summed paths, root stems, zero-length resolution edges, and metadata preservation. |
| `mark` | Remove intermediate annotation flags and repeated suffix/selection checks. Assert the resulting named clades, leaf names and lengths together. Keep prefix/separator, no-match and single-target cases. |
| `drop`, `transfer`, `nhx2nwk` | Remove duplicate smoke tests. Existing scenarios now check actual support removal, unmatched-name fill and preserved tip names instead of file existence alone. |
| `shuffle` | Remove helper-list checks and repeated tip-count smoke tests. Keep output reproducibility, branch-length multiset, changed topology/RF, and unnamed/deep/unrooted inputs. |
| `dist`, `cladefreq` | Keep exact RF rows and combine frequency values with reference-clade checks on the same input. Weighted and parallel execution contracts remain. |
| `intersection` | Remove equality/startswith/endswith wrapper tests. Test both partial-match directions through actual removal selection instead. |
| `util_tree_io` | Remove plain file/explicit-format reads subsumed by annotated, named and format-specific input tests. |
| `asr_comparison` | Check numeric alignment in the existing ranked/unranked overview, avoiding another independent figure fixture. |
| `gaussian`, `gaussian_tree` | Remove checks of SuperLU's internal factor representation and comparison of two views of the same covariance factory. Independent solves, likelihoods, indefinite-matrix rejection and consumer agreement remain. |
| `image_helpers` | Remove exact worker counts, fetch-buffer constants and trivial scoring-helper checks. Keep candidate ordering, fallback, cache reuse/invalidation and thread-bound resource cleanup. Merge GET-only retry policy into the real session test. |
| `shift_native_limits` | Remove a snapshot of default budget constants. Keep caps under explicit budgets, traversal limits, dispatch agreement and large-tree numerical checks. |
| `transformed_continuous_asr` | Remove finite-fit smoke cases and a minimum optimizer-grid count. Keep Brownian reductions, model domains, estimated parameter bounds and profile intervals. |

## Reduced repeated cases

- Reserved NHX properties are rejected before target-format and keep/drop handling.
  Three cases exercise every reserved key and output format instead of all 18
  combinations; output-preservation tests remain separate.
- Random whitening uses four fixed seeds covering both root-variance modes and
  zero/nonzero slopes instead of 20 seeds. Explicit observation-order, singleton,
  independent/Brownian and structure-reuse cases remain.
- Wide-tree GLS retains partially observed and fully observed 512-tip trees,
  using fixed and random roots respectively. Small dense-oracle tests already
  cross root modes with observation patterns; intermediate size 448 adds no mode.
- The noise-shift oracle uses 15 rather than 240 random examples: all three
  topologies with exact-anchor and all-noisy observations. Dedicated extreme-scale,
  multimodal-likelihood and zero-edge regressions remain unchanged.

This removes 75 test functions and 114 collected cases from the baseline suite,
plus 225 repetitions inside the noise-shift test. No runtime speedup is claimed.

## Retention decisions

Independent dense/quadrature/closed-form numerical oracles, reference fixtures,
identifiability and unit-invariance checks still expose plausible wrong scientific
results. Mocks at network, process and output-failure boundaries test real parsing,
fallback or rollback behavior and are not merely mock self-checks. File integrity,
atomic publication, unsafe media rejection and worker cleanup retain distinct
failure paths. Those tests, external-tool integrations and the delivery checks were
kept; assertion tolerances and repository coverage gates were not relaxed.

## Verification

- `tools/check.py full`: 3,976 passed, 83 skipped; lint, formatting, types,
  dependency consistency, security audit and maintainability checks passed.
  Reported coverage was 85%; existing thresholds were unchanged.
- `tools/check.py quick -- tests/test_label.py tests/test_mark.py
  tests/test_discrete_asr_models.py tests/test_optimization.py`: 74 passed,
  covering final assertion edits and the concurrently completed audit fixes.
- `tools/check.py dist`: version 0.43.26 wheel/sdist contents and reproducibility
  passed. These full and distribution checks together cover the release entrypoint.

Validation used an isolated Python 3.14 environment, with ETE4 4.4.0 rebuilt from
source after its downloaded wheel failed to import. The existing local virtual
environment was not modified. Concurrent audit fixes were committed separately as
`e5c46da`; their 33 additional cases explain the difference from the initial
4,140-case inventory: 4,140 + 33 - 114 = 4,059 final cases, including skips.
