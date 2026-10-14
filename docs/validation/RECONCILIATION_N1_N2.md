# Common contract, N1 import and N2 static drawing — 2026-10-14

This increment is isolated from the clean primary nwkit checkout at
`d1e7f97a84fc069dd336b6e2661ecdb5cc9b80e9`. The implementation checkout is
`/Users/kf/.codex/worktrees/reconciliation-viewer-pnxjotzn/nwkit`; the aioi peer
is `/Users/kf/.codex/worktrees/reconciliation-contract/aioi`. Existing primary
work, historical runs and the original aioi specification bundle are preserved.
Nothing was pushed or published.

## Delivered stages

| Stage | Implementation | Evidence and limits |
|---|---|---|
| Contract | Versioned document/schema, independent namespaces, explicit histories/events/selection/support/units, X1 boundary records, reserved view/comparison files | Byte-identical aioi schemas and SHA manifest; syntax and semantic checks; comparison/view execution is deferred |
| N1 | Canonical JSON, single Newick/NHX, phyloXML and declared recPhyloXML profile; saved TSV adapter; JSON/XML export | Round trips, reordered children/internal naming, metadata preservation, explicit ambiguous/unsupported-profile errors; no generic GeneRax-run directory adapter |
| N2 | Neutral scene, standalone/one upper level with multiple lower histories, explicit transfer arrows and losses, display-only clade selection, SVG/PDF/PNG | Real vector files, selectable TrueType PDF text, tiled A4/A3 pagination, continuation labels, transactional related outputs |
| A0 peer | Representative derivation definition and exact development oracle | 60 conditions / 1,500 topology conditions; native reconstruction/output is not implemented |

The document schema is packaged under `nwkit/data_schema`; the authoritative
hashes and peer bases are in [the common manifest](../guides/reconciliation_contract_manifest.json).
Known representative selection rules, support conditioning and display-only
selection metadata are checked structurally as well as semantically. Cross-tree
references, connectedness, event placement and X1 boundary identity require the
semantic validator; passing JSON Schema alone does not establish a drawable
or scientifically complete history.

## Actual validation

ARM64 macOS, Python **3.14.7**, Matplotlib **3.10.9**, jsonschema **4.26.0** and
pypdf **6.20.1**. A dedicated venv lives alongside this worktree. The primary
checkout's existing Intel Python venv failed with `Bad CPU type in executable`;
it was preserved. Native bundled Poppler was used for PDF QA because the
existing `/usr/local/bin/pdftoppm` also targets Intel. No primary environment
or installation was modified.

The affected regression selection below passes **586 tests** (including the
45 new contract/scene test cases). Ruff check and formatting pass, and mypy
reports no issues in **267 source files**. This is an affected-suite regression,
not a complete nwkit release or multi-Python/OS test matrix. Maintainability
checks meet the existing hard limits; complexity warnings remain in the log.
[The quick-check log](reconciliation-n1-n2-quick.log) records the actual commands.

```sh
python tools/check.py quick -- \
  tests/test_reconciliation_document.py tests/test_reconciliation_scene.py \
  tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_draw.py tests/test_draw_robustness.py tests/test_draw_branch_ids.py \
  tests/test_result_plot.py tests/test_output_transaction.py tests/test_reconcile.py \
  tests/test_reconciliation_exports.py tests/test_util_tree_io.py \
  tests/test_tree_outputs.py tests/test_rooting_state.py -rs
python tools/check_maintainability.py
python tools/check.py dist
```

Distribution validation checks direct/sdist-built wheels and independently
normalized sdists, their byte equality, and packaged schemas/examples. Its log
is stored **outside** `docs`, at
`/Users/kf/.codex/worktrees/reconciliation-viewer-pnxjotzn/validation/distribution.log`.
A first attempt logged inside the `graft docs` input and failed reproducibility
solely because the included log grew between builds. That failure is retained
outside the package inputs as `distribution-initial-self-log-failure.log`.
The external distribution log and wheel-smoke JSON are the actual final evidence,
not a claim that the entire release gate ran.

PDF/SVG tests check all recorded event IDs, distinct lower-tree lanes,
child-order invariance, directed donor/recipient endpoints, XML source confidence
without a fabricated posterior, long labels, repeated page context and
cross-page continuation text. Continuation positions are independently checked
against sampled Matplotlib-rendered curves (including a curve crossing a page
twice) and bent branches. Visual QA found and corrected the earlier endpoint-based
continuation label positions. Input/output/audit alias tests and injected
second-output failures verify preservation of preexisting files. CLI tests also
prove document drawing does not read implicit stdin or call reconciliation
inference. Selected incoming transfers retain external donor context; normalized
selection documents retain all original trees/events/summaries.

A synthetic balanced **1,000 species / 1,000 genes** geometry check verifies
finite placement of all **1,999 recorded events**. It is not a benchmark of
readable dense output, peak memory, HTML latency or native inference.

The toy PDF is synthetic, with duplication, speciation, explicit loss and C→A
transfer. It was rendered with Poppler and visually inspected for readable
legend, labels, separated lanes, transfer direction and clipping. Source fixture:
[recPhyloXML](../../examples/reconciliation/toy.recphyloxml).
A three-page synthetic transfer test was also rendered and inspected after
correcting continuation labels. Review output: `output/pdf/reconciliation-toy.pdf` (ignored local artifact).

## Remaining gates and order

1. **A1** native DL traceback with independent event-weight re-scoring, fixed-state
   fidelity, numerical certificates and deterministic tie rules.
2. **A2** native DTL recipient/donor paths; **A3** native DTLX latent H/X1 mapping;
   **A4** saved-run restoration/common-document producer and rollback/resource tests.
3. **N3** arbitrary-panel comparison and model-aware differences; five methods
   become an aioi wrapper with explicit missing/failed/partial statuses and
   separately measured inference/traceback costs.
4. **N4** offline interactive HTML, browser PDF parity and reproducible state;
   three-layer drawing follows the agreed next-stage boundary.

Expectation maps, topology-transform rendering, physical-time layout, unsupported
recPhyloXML ghost events and incomplete transfer placements fail explicitly.
No placement, loss, support, ILS label or history is inferred from missing data.
Global crossing minimization and large dense readability are not certified.
See [the user guide](../guides/RECONCILIATION_DOCUMENT.md) for CLI and profile limits.
