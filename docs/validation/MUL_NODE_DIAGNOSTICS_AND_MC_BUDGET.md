# Node Diagnostics And Fixed-Data MC Budget Follow-Up

Executed 2026-10-03 on local NWKIT 0.43.37. These are optional diagnostics,
not a change to the default D+L model, a production WGD test, or biological
node-origin posterior probabilities. Previous studies and their failures
remain unchanged. The budget experiment follows the
[predeclared protocol](../../examples/mul-locus/CONVERGENCE.md).

## Auditable Node Output And GeneGalleon Join

D+L `--node-out` retains every optimal assignment for every globally tied
best candidate, including every tip/internal node and the costs that reconstruct
each per-gene score. Rooted topology and descendant-tip clade IDs ignore child
order, branch lengths, support and internal labels. Legacy `--report` still
exports only the first best candidate. Caps abort the whole staged bundle;
assignments are never truncated or treated as probabilities. Duplicate literal
MUL occurrences cannot supply distinct clade identities and are rejected only
for this optional output. Other score models reject the option explicitly.

GeneGalleon optionally runs a separate per-family search on its original full
species tree and unmodified gene tree. Its join verifies tree/clade identity,
species mappings, candidate metadata, complete mappings, uniqueness and score
totals before adding any node attributes. Original WGD/SSD/unresolved evidence
and TSV columns remain unchanged, including when native MUL duplication labels
disagree. Consistency means agreement of mapped occurrence clades/species/local
costs across the complete co-optimal set, not biological certainty.

The output includes `node_diagnostics.tsv`, raw `mul/nodes.tsv`, model/scores,
NHX properties and `duplication_origins_mul.pdf`. Colors encode original origin
evidence; circles/triangles encode consistent/ambiguous MUL states. Node labels
join to the table. A reproduced deep-tree label collision was fixed by measuring
actual glyph widths and allocating space per depth level. Rendered examples
were visually inspected; 80-tip comb trees and 200-character wide labels pass
artist-bound checks for overlap and clipping. Disabled reruns clear stage-owned
attributes and invalidate enabled-mode cached bundles. Failure preserves the
previous family ZIP and provenance, including after output-store archiving.

See the [NWKIT command guide](../guides/MUL_RECONCILE.md) and
[GeneGalleon settings/output guide](https://github.com/kfuku52/genegalleon/blob/main/docs/wgd-ssd.md#optional-mul-node-diagnostics).

## Fixed-Data Budget Experiment

| Draws Per Bank | Histogram Completed / 21 | Histogram Support Failures | Detection / Hybrid Completed |
| ---: | ---: | ---: | ---: |
| 2,000 | 1 | 20 | 21 / 21 each |
| 4,000 | 5 | 16 | 21 / 21 each |
| 8,000 | 10 | 11 | 21 / 21 each |

All 21 observed dataset files are byte-identical across budgets. Models,
selection, seeds, 10 families/dataset, 19 null replicates per each of four
generating grid points, exact hidden-tip limit four and work cap 100,000
are unchanged. There are no reference-generation failures. Both follow-up
runners exit 1 because histogram support failures remain, retaining all
planned analyses rather than dropping failed cases.

| Method / Measure | 2,000 | 4,000 | 8,000 |
| --- | ---: | ---: | ---: |
| Detection: correct point parent among six alternatives | 5 / 6 | 6 / 6 | 6 / 6 |
| Hybrid: correct point parent among six alternatives | 5 / 6 | 6 / 6 | 6 / 6 |
| Detection: median contrast-MC interval width, 21 datasets | 139.256 | 96.013 | 64.360 |
| Hybrid: median contrast-MC interval width, 21 datasets | 140.228 | 94.728 | 64.662 |
| Detection: alternatives passing upper-MC P <= 0.05 | 0 / 6 | 0 / 6 | 0 / 6 |
| Hybrid: alternatives passing upper-MC P <= 0.05 | 0 / 6 | 0 / 6 | 0 / 6 |

No completed observed contrast interval has infinite bounds; failed analyses
have no valid interval and are not assigned zero width. Histogram interval
medians cannot be compared as a fixed panel because completion changes.
Both weighted methods have point P = 0.05 for all six alternatives at every
budget. At 8,000, their two baseline alternatives have upper-MC P = 0.25;
turnover and missing-ILS alternatives still have upper-MC P = 1. Neither
weighted method reports an event for the 12 on-grid null datasets at any budget.

Median probability-interval widths use the same dependent pattern/bank cells
at all budgets (820 baseline, 980 turnover, 960 missing-ILS per method):

| Scenario / Method | 2,000 | 4,000 | 8,000 |
| --- | ---: | ---: | ---: |
| Baseline / Detection | 0.081268 | 0.057657 | 0.039174 |
| Baseline / Hybrid | 0.081221 | 0.056881 | 0.039509 |
| Turnover / Detection | 0.083388 | 0.055821 | 0.038164 |
| Turnover / Hybrid | 0.081660 | 0.053779 | 0.037755 |
| Missing-ILS / Detection | 0.107970 | 0.074491 | 0.051836 |
| Missing-ILS / Hybrid | 0.107229 | 0.074028 | 0.051556 |

Increasing the budget repairs the point-parent error in this particular
missing-ILS/allop-B dataset, but this is one dataset per truth case and only
two budget doublings. It proves neither convergence, parent identifiability,
empirical power, a 5% error guarantee, nor robustness to sequence/tree inference
error. Histogram support remains unreliable, and conservative MC uncertainty
still blocks event decisions. D+L remains the workflow default; DL+ILS stays
opt-in research functionality.

## Reproducibility And Audit Repair

Directories are `/tmp/nwkit-locus-conditional-probe-v2-20261003`,
`/tmp/nwkit-locus-convergence-4000-20261003` and
`/tmp/nwkit-locus-convergence-8000-20261003`. The same pre-run source archive
is copied into both follow-ups. Its SHA-256 is
`b50053b9159913bdb698d31ef89bb3330db827b8330b022777ed80122f21a3e2`.
All scientific source entries match their stored protocol hashes, and the
follow-up source snapshots are identical. Relative to the 2,000 control, only
the previously documented integer-count guard in `mul_locus_integral.py`
differs; every stored observed fit replays exactly.

The independent arithmetic audit recomputes point P, both MC-overlap endpoints
and the four-grid supremum for **10,792 complete null searches** across the
three budgets (3,268 / 3,572 / 3,952). All **142 completed observed fits**
replay exactly from saved banks (43 / 47 / 52), including histograms.

This broader replay exposed a research helper defect: `integration.read_banks`
restored pooled counts but omitted saved histogram strata. The loader now
restores stratum counts, weights, selected counts and budgets. Flat and
stratified round trips both preserve exact probabilities. This loader repair
occurred **after** sampling and is recorded as a post-freeze change; it does
not alter the stored studies or their source archives. The probe's sampling
path with no `--source` never calls this loader. The audit includes the repair
when replaying saved fits, not when retroactively generating new results.

Follow-up commands differ only in the budget and fresh output path:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/integration.py \
  --output /tmp/nwkit-locus-convergence-4000-20261003 \
  --samples 4000 --families 10 --bootstrap 19 --replicates 1 \
  --seed 20261111 --bank-seed 20261112 --calibration-seed 20261113 \
  --reference conditional --interval chernoff-kl
```

The second command substitutes `8000` in both budget and directory name.
Source archives, protocols, bank moments, all failures and audits are retained.

## Executed Checks

NWKIT: isolated CPython 3.12.14, macOS 27 ARM64. Import/pip preflight passes.

```sh
PY=/tmp/nwkit-wgd-dev-20261001/bin/python
"$PY" tools/check.py test -- tests/test_mul_reconcile_nodes.py \
  tests/test_mul_reconcile.py tests/test_cli.py tests/test_cli_contracts.py \
  tests/test_interface_conventions.py tests/test_provenance.py \
  tests/test_output_transaction.py tests/test_rooting_state.py \
  tests/test_util_tree_io.py tests/test_tree_outputs.py
"$PY" tools/check.py quick -- tests/test_mul_reconcile_nodes.py \
  tests/test_mul_locus_integration_study.py -x -rs
"$PY" -m ruff check examples/mul-locus/integration.py
"$PY" -m ruff format --check examples/mul-locus/integration.py
"$PY" tools/check_maintainability.py
```

Broad tests: 350 passed before the additional multi-gene case. Final focused
quick: 20 passed, no skips; it includes that case and both loader round trips.
Ruff/format and incremental mypy (251 sources) pass. Maintainability respects
all hard limits; existing comparison warnings are not suppressed.

GeneGalleon: Docker Linux ARM64/CPython 3.12.14,
`local/genegalleon:mul-node-diagnostics-20261003`, immutable image
`sha256:304bd0724a2d9b1a1d063da53b878488b3e9d0be71ec4f26c260e3dd58c0a6f5`.
This is a local NWKIT package-source overlay on the previously verified locus
integration runtime, not a pinned upstream default. Daily freshness passes.
For each following command, set `GG_TEST_RUNTIME=docker` and
`GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:mul-node-diagnostics-20261003`:

```sh
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_wgd_ssd_contracts.py \
  workflow/tests/test_wgd_ssd_gene_stage.py -x -rs
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_wgd_ssd_runtime.py \
  workflow/tests/test_native_mul_reconcile_runtime.py \
  workflow/tests/test_wgd_evidence.py -x -rs
bash ./dev check static
bash workflow/tests/run_in_runtime.sh env GIT_CONFIG_COUNT=1 \
  GIT_CONFIG_KEY_0=safe.directory \
  GIT_CONFIG_VALUE_0=/Users/kf/repos/genegalleon bash ./dev lint
```

Contracts/stage: 159 passed. Broader real native count/synteny/dS, origin and
MUL workflow checks: 45 passed. Static suite: 240 passed. No skips in these
selections. Lint/config check passes (8 entrypoints, 23 common parameters).
Final diff/whitespace and feature-list ordering are checked separately.

SIF, minimum Python 3.10, complete release/coverage/security/distribution lanes
and empirical real-data inference-error calibration were not run. Host evidence
does not establish SIF compatibility. No commit, push, Wiki deployment,
threshold retuning or curated-input changes are part of this local follow-up.
