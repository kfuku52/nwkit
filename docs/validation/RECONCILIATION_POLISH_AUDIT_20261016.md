# Reconciliation drawing polish audit, 2026-10-16

Local follow-up to the drawing commits through `50a217f5`, covering the native
aioi examples and the reusable nwkit renderer. No inference/model code, source
event schema, package version, or previously delivered slide PDF is changed.

## Confirmed corrections

- Long wrapped labels could extend past the export canvas; the reproducer had
  six clipped labels. Shared export bounds now include text, thick corridor
  strokes and transfer-curve control hulls. A translated/expanded display copy
  keeps event lookup coordinates aligned without mutating the supplied scene.
- The nested offset algorithm averaged corner directions, creating diagonal
  segments from rectangular hosts. Segment-normal intersections now retain
  right angles. Nested histories reuse frontier lanes so serial duplication
  ancestors do not consume extra concurrent-copy lanes. Concurrent speciation
  forks are staggered within the junction, with both daughters sharing their
  parent's fork x coordinate. Recorded curved middle hosts remain curved.
- Marginalized terminal sides were omitted from the bottom nested layer.
  Nested and two-layer views now share the same classification and annotation:
  a red loss cross only for the documented pure-loss/full-detection boundary;
  otherwise a gray E diamond. Stub endpoints are display midpoints, with no
  inferred time, posterior support, or additional recorded event.
- Pure-loss interpretation also checks actual recorded leaf placement, rather
  than trusting mapping metadata that contradicts an observed terminal copy.
- Long wide labels in three-layer views could overlap neighboring layer
  columns: four overlapping pairs in the reproducer. Conservative column widths
  remove those overlaps; multiline lane clearance is retained.
- Helvetica lacks the arrow glyph used in generated level/transformation
  headings on this host. Generated headings now use ASCII words. Source labels
  are retained verbatim; arbitrary missing source glyphs are not substituted.

## Executed verification

Runtime: macOS ARM64, Python 3.14.7, nwkit 0.43.50. Compiled ETE/NumPy/SciPy,
pandas/Matplotlib/Pillow import preflight and `pip check` passed. HTML JavaScript
tests used the bundled ARM64 Node explicitly through `NWKIT_TEST_NODE`.

From the worktree root, with that environment's Python and the checkout on
`PYTHONPATH`:

```sh
python tools/check.py quick -- \
  tests/test_reconciliation_bounds.py tests/test_reconciliation_corridors.py \
  tests/test_reconciliation_document.py tests/test_reconciliation_scene.py \
  tests/test_reconciliation_view.py tests/test_reconciliation_html.py \
  tests/test_reconciliation_exports.py tests/test_draw.py \
  tests/test_draw_robustness.py tests/test_draw_branch_ids.py \
  tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_output_transaction.py tests/test_provenance.py \
  tests/test_util_tree_io.py tests/test_tree_outputs.py tests/test_rooting_state.py -rs
python tools/check_maintainability.py
```

Final selection: **659 passed, no skips**, in 22.78 s. Ruff lint and formatting
passed (609 files); mypy passed (275 source files). Maintainability hard limits
passed: 3891 functions, average 6.96 (informational), maximum 50 (the existing
draw exception). Existing baseline-increase warnings remain informational; no
complexity limit or tolerance was loosened.

Actual CLI exports: **15 fixtures x PDF/SVG/PNG/HTML = 60 exports**, plus a
Helvetica PDF. The fixtures include all seven native aioi records used for the
four-scenario and DL/DTL/DTLX slide figures, recPhyloXML, nested 2/3/8-copy trees,
nested loss and incomplete detection, wide labels, and a recorded middle
transfer curve. Native source hashes were verified before use. Every input byte
snapshot remained unchanged; each normalized SVG companion retained the full
validated source record. PDFs parse, SVGs are valid XML, and HTML retains its
offline embedded PDF without external scripts/fetches. The Helvetica export
embeds `/EDQXKT+Helvetica` and emits no missing-glyph warning. System-font
subsetting reports nonfatal font-table/timestamp metadata warnings.

Regression checks inspect actual Matplotlib text boxes for clipping and nested
column overlap, rectangular copy paths, distinct fork positions, nested
loss/E annotations and exact source-record preservation. Representative
DL, nested-copy, loss, transfer and wide-label raster outputs were visually
inspected, including the final versions.

`python tools/check.py dist` passed on a separate source snapshot, with
`SOURCE_DATE_EPOCH` taken from the starting commit. Direct and sdist-built wheels
are byte-identical, independently normalized source archives match, and
contents/metadata checks pass, including the new `reconciliation_bounds.py`.
All reconciliation modules in the wheel match the working source bytes; an
extracted-wheel CLI smoke check rendered the nested-copy SVG successfully.
Existing worktree build artifacts were not altered. This packaging check does
not imply that the full regression/release gate was executed.

Raw logs, export inventory/hashes, local audit driver, sources and previews:
`/Users/kf/.codex/visualizations/2026/10/16/01a1431c-36a1-78e2-865a-68864ecca85f/nwkit-audit/`.

## Scope and limits

This is affected verification, not the complete repository regression/release
gate. Unrelated numerical/study suites, optional external inference backends,
live browser interaction and cross-platform font rendering were not run.
Dense arbitrary crossing minimization remains uncertified. Export text bounds
are conservative character estimates, not a proof for every font/script.
No cause-specific ILS interpretation is introduced for X1, and no integrated E
history is expanded into a reconstructed loss count/time.
