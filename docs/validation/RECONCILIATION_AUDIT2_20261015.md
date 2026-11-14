# Reconciliation audit II, 2026-10-15

This extends the previous audit with numerical boundary cases, deterministic
filesystem faults, concurrent writers, rendering semantics and deep XML. Previous
work, reports, wheels and native results are preserved; no push/publication.

## Confirmed corrections

- Near-critical finite DL/DTL passed summed inference but failed independent
  replay's bounded Neumann witness. Independent DL replay now solves triangular
  mean/extinction equations directly. DTL can build and eliminate its own mean
  system and checks the resulting positive witness. It never reads a backend's
  factorization or stored scores. Weighted extinction residuals are converted
  back to an absolute bound. Original likelihood/search equations are retained.
- A pure-observation single lineage with detection1e-14 previously returned log
  reconciliation0.000799597430194865, although the conditioned probability is
  exactly1. Stable Bernoulli observation unions correct this to log0. Tiny
  nontrivial normalizers (Z<1e-10) whose estimated error exceeds1e-8 relative
  precision explicitly fail with numerical exit4. This guard is not a proof of
  accuracy for every parameter range. It does not change observation data/model.
- Fixed-state serialization validates tree structure before indexed traversal.
- X1 highlighting now marks branches only; biological D/S/T symbols and legends
  retain their meanings. Selected lower clades keep upper-context event positions;
  collapsed descendants no longer leave stale displayed-event coordinates.
- A valid recPhyloXML tree deeper than Python's recursion limit imported but could
  not export. An iterative XML serialization path handles deep clades/opaque
  metadata without changing the global recursion limit. Ordinary ElementTree
  formatting remains identical. Expanded namespaces, attributes and text/tails
  survive the deep path.

## Tests and evidence

Native Release **30/30**, affected ASan/UBSan **5/5** including the Darwin fault
interposer. LeakSanitizer remains unsupported here and is not claimed. Exact
Fraction forward history/channel oracle **1500/1500**. Existing16176 positive
history conditions remain covered; **300** heterogeneous detection/zero-detection,
uneven species trees, D/T/L rates and X-channel endpoint conditions were added
(seed2026101502; three species shapes). Maximum added-case replay discrepancy is
9.23706e-14. Near-critical DL and DTL with respective propensities0.99999 now
trace successfully; the old binary exits4. Pure-observation detection1e-6/1e-14/
1e-20 and undefined zero-detection conditioning are explicit tests.

Real filesystem fault injection, test-only Darwin dylib not linked into the
product: install rename failure restores all five original hashes; rollback
rename failure preserves the exact original events backup and reports its path;
cleanup failure after deleting one old backup preserves the fully committed
bundle. All release and sanitizer fault cases passed. An initial MODULE build
could not be loaded by dyld; SHARED fixed the harness. An accidentally overlapping
fixture execution was discarded and rerun serially. Normal fixtures use separate
build directories. Twenty two-writer trials produced20 success/20 resource
refusals, each with a complete verified inventory and no remaining owned lock.
This is not exhaustive power-failure/OS fault coverage.

The affected nwkit `tools/check.py quick` selection (document/view/scene/XML/export/
HTML/CLI/conventions/provenance/output transactions/legacy draw) passed **492**,
no skips. Ruff lint/format, mypy **273 source files**, and maintainability hard
limits passed. Deep XML exercises1050 levels; full trees/reconciliation records
round-trip exactly, and nested opaque namespaces/text/tails retain their expanded
meaning. A1100-copy comb Newick import also completed with2199 nodes and absent
lengths remaining null. X1 source colors and the legend were visually checked in
PNG/SVG; no new physical X1 event or time is introduced.

Distribution inventories, byte-identical direct/sdist wheels, normalized source
archives and metadata validation passed, after backing up all previous build,
dist and direct-dist directories. Optional JS tests explicitly use the bundled
ARM64 Node. These are affected tests, not all optional/live backends or the entire
repository suite. Native runs remain independent of Python/Rust/nwkit.

The final five-method run/equality audit is recorded by aioi, with fresh timings
kept separate from previous audits. Existing host jobs are preserved; there is no
idle-host or speed superiority claim.

Raw logs/previews: task-local `validation/audit2-20261015/`. This is a follow-up to
[the first audit](RECONCILIATION_AUDIT_20261015.md), not a replacement for its evidence.
