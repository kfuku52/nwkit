# Reconciliation drawing audit, 2026-10-15

Preserved the existing isolated N1–N4 checkout and all prior build/distribution
artifacts. No pushes, publication or version changes. Native inference remains in
aioi; this viewer only consumes saved files.

## Confirmed fixes

- Panel/tree identity scopes HTML details, linked selection and neutral event
  coordinates; repeated reconciliation/event IDs cannot overwrite other panels.
- Initial HTML collapse keeps expandable geometry. Expand/reset work in-browser;
  downloadable PDF preserves the initial view. Embedded upper context collapse is
  disabled. Missing-history tree panels honor display state; ambiguous clade IDs
  across source documents are rejected.
- Comparison canvas retains the width of nested source scenes, preventing clipping.
- Hash verification, parsed document and ZIP use one captured source-byte snapshot,
  even if an external process edits the file while drawing.
- Identity/X1 transformations reconstruct the whole target topology. Structural
  scalar container corruption raises ValueError, not unhandled TypeError. Recorded
  event comparison includes observed arity and marginalized side when available.

## Executed evidence

Interpreter: dedicated Python 3.14.7 ARM64; optional shipped-JavaScript tests use
`NWKIT_TEST_NODE` with the bundled ARM64 Node.js. The default Intel Node executable
on this host cannot execute (Bad CPU type); this failure was preserved and corrected
by selecting the already installed runtime, not installing another dependency.

`tools/check.py quick -- tests/test_reconciliation_document.py
 tests/test_reconciliation_view.py tests/test_reconciliation_scene.py
 tests/test_reconciliation_exports.py tests/test_reconciliation_html.py
 tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py
 tests/test_provenance.py tests/test_output_transaction.py tests/test_draw*.py
 tests/test_draw.py -rs`: **487 passed**, no skips. Ruff lint/format and mypy
**272 source files** passed. Maintainability hard limits passed (3877 functions).
This is the affected offline selection, not the entire repository suite.

Systematic corruption covers 254 required scalar/reference mutations across the
synthetic saved document. Whole-target negative fixtures preserve a locally valid
X1 block but change a distant topology. All **336** previously exported native
history documents passed the strengthened whole-target verifier.

`tools/check.py dist`: direct/sdist wheel inventories and byte reproducibility,
normalized source distributions and metadata validation passed. Earlier build,
`dist` and `direct-dist` were moved to an audit backup before this command.

Real Codex in-app browser, local loopback server only: an initial paired lower
root collapse hid 20 primitives; selecting its root and toggling restored all 20,
including gA1/gA2/gA3. The downloaded comparison ZIP retained all five exact source
byte streams and hashes and loaded through read_comparison after extraction.
Standalone SVG keyboard automation could not focus a zero-height root group in
this browser backend; coordinate input succeeded. Shipped JavaScript keyboard,
collapse/reset and panel-detail handlers passed the separate Node DOM tests.
PDF output retained selectable vector content; all nine pages were rendered and
visually checked. Temporary server and tab were closed.

## Limits

Native inference/posterior inference, dated event placement and inference-input
rehashing are outside the viewer. A shared input hash is a producer declaration;
actual document bytes are verified, and fixed-gene mode additionally checks actual
serialized G topology/lengths/units. Browser Print depends on the host dialog;
printing completion was not automated. Optional Node tests explicitly skip when
Node is absent/unusable unless an explicit NWKIT_TEST_NODE is supplied. No live
external services, release/publication gate or full repository test claim.

Raw audit logs and temporary rendered previews are under
`/Users/kf/.codex/worktrees/reconciliation-viewer-pnxjotzn/validation/audit-20261015/`.
