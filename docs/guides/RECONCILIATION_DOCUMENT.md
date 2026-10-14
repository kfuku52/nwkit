# Generic reconciliation documents and embedded drawing

The file contract is `reconciliation_document` version **0.1.0**, independent
of any inference program's run-schema version. The authoritative schema is
`nwkit/data_schema/reconciliation_document.schema.json`; aioi keeps an identical
copy for its future C++ producer. Schema bytes and their SHA-256 are frozen in
the stage readout. Changes require a coordinated producer/consumer update and
a version change when meaning or compatibility changes.

## Contract

Documents contain flat `trees`, `reconciliations`, optional topology
`transformations`, and `provenance`. Every tree has its own node namespace,
explicit root, labels, child references and nullable branch lengths, with a
declared `branch_length_unit`. Parent references, roots, connectivity, cycles,
duplicate IDs and observed leaf identities are validated. Node IDs are opaque;
they are never equated across trees by array index or name alone.

Each reconciliation names an upper and lower tree and records:

- `history_kind`: representative, sample, expectation, or unspecified.
- `selection`: representatives must declare `method`. This is a choice rule,
  not a posterior probability or marginal-tree MAP claim.
- `model`: source model identity and conditions. Unknown conditions stay unknown.
- An ordered `events` array, allowing multiple events on the same lower node's
  incoming branch. Every event names its lower node and nullable upper placement.
- Explicit donor, recipient and recipient-entry lower node for a full transfer.
- Optional time intervals with units and support with its kind and conditioning.
- `summary` for aggregate quantities. A loss count is not a list of placed losses.

Event types include speciation, duplication, loss, leaf, origin, unresolved,
empty observation, full transfer, and source-format transfer out/in events.
Absent support stays absent/null. A representative's presence is never encoded
as support 1. JSON NaN/Inf and duplicate keys are rejected, including metadata.
Unknown core fields are rejected; extension data belongs in `annotations`.

X1 uses a separate transformation connecting source H and target G. Each block
records adjacent consumed source/target node pairs and three ordered boundary
clade mappings A,B,C for `((A,B),C)`. Resolution 1 is `((A,C),B)` and resolution 2
is `((B,C),A)`. Latent consumed blocks must not overlap; boundary leaf identities
and the target local split are validated. X1 has no fabricated species-branch
placement. N2 does **not** render transformation documents yet and refuses to
omit them. Arbitrary other transformation models are preserved as extension
blocks but do not acquire X1 semantics.

View and comparison schemas are reserved contracts for N3/N4. Their files are
not yet accepted as executable draw configurations. Comparison documents have
an arbitrary number of panels, separate inference/history statuses, source
hashes and distinct inference/traceback timings; they do not fix five methods.

## N1 adapters

`--document` reads canonical JSON, rooted phyloXML, the supported recPhyloXML
profile, or a single Newick/NHX tree. Newick/NHX has no inferred reconciliation.
Single-tree branch lengths and annotations are retained; the new scene is
currently topological. Existing ordinary `nwkit draw` is unchanged.

XML uses defusedxml and rejects DTDs/entities/external references. The supported
recPhyloXML profile requires one `spTree` and named, unique upper references.
It supports leaf, speciation, duplication, loss, branchingOut and transferBack.
An outgoing transfer must have exactly one immediate recipient-entry child.
Ambiguous/unpaired transfers, bifurcationOut/ghost events, and unsupported
top-level elements are rejected. timeSlice and source confidence remain source
attributes: time slices do not become physical times, and confidence does not
become a calibrated posterior. Non-event XML node annotations are kept as
`opaque_xml`, rather than silently discarded; they are not all rendered.

Import IDs derive from canonical subtree signatures and collision checks, not
traversal positions. Identical unobserved subtrees have occurrence suffixes in
canonical order; they have no invented distinct biological identities.
Internal renaming and child order do not change structural IDs.

Canonical JSON round trips every contract field. XML round trips the supported
imported profile. XML export refuses model conditions, selection rules,
summaries, X1 or extensions it cannot represent; use canonical JSON for these.

Legacy saved TSV can use `--view embedded`. Existing validation checks its
gene/parent clade IDs and species topology. `implied_losses` remains a summary.
A legacy mapped transfer without a donor/recipient pair cannot be embedded;
the existing paired report remains available without this new option.

## N2 examples

```sh
nwkit draw --document examples/reconciliation/toy.recphyloxml \
  --outfile toy.pdf --document-out normalized.json
nwkit draw --document normalized.json --outfile toy.svg
nwkit draw --document examples/reconciliation/toy.recphyloxml \
  --outfile toy.svg --document-out roundtrip.recphyloxml
nwkit draw --infile gene.nwk --species-tree species.nwk \
  --reconciliation events.tsv --view embedded --outfile embedded.pdf
```

The toy input is a synthetic display fixture with D, T, speciation and an
explicit loss, not an aioi or GeneRax inference result. No fitting, reconciliation
inference, new event classification or numerical-score calculation is run.

The neutral scene contains paths, directed transfers, event symbols and labels
with stable semantic IDs. All recorded events must be displayed. It currently
supports standalone trees or one upper tree with one reconciliation level and
multiple lower histories. Unmapped/incomplete placements, multi-level embedding
and expectation maps fail explicitly instead of becoming fabricated histories.
Upper and lower branch lengths are retained but not mixed as a time scale.

Display-only subtree selection uses `--document-clade TREE_ID/NODE_ID`; IDs are
available in normalized JSON. It retains the complete upper context, records
selection in provenance, retains all original scientific trees/events and
full-source summaries, and does not
recompute inference or invent a new origin. Upper-tree pruning is rejected.
Incoming transfer arrows retain their recorded donor context outside the selected
lower clade. The selected document remains a complete source document with a
display selection, rather than a new inferred family.

SVG and PDF use the same scene. PDF embeds TrueType fonts and selectable text.
`--pdf-pages auto|single|tiled` and `--pdf-page-size A4|A3` control pagination;
paper sizes are landscape. Tall PDFs paginate by default, repeat context, and
mark crossing branches/transfers with continuation page references. Single-page
outputs over 200 inches are rejected. PDF options with SVG/PNG are rejected.
SVG uses text and semantic artist IDs. `--font-family`, `--font-size`,
`--figure-width`, `--legend`, and `--transparent` are supported. Other nondefault
ordinary-tree/report options are rejected rather than ignored.
The new document mode defaults to DejaVu Sans; an explicit font-family option is
honored. Ordinary draw defaults are unchanged.

Figure and normalized-document writes are staged as a recoverable bundle.
Validation/write failures preserve previous files. Input/output aliases,
including audit and normalized-output paths, are rejected. Document drawing
does not consume the default implicit stdin; actual inputs appear in `--audit`.

## Remaining stages

N3 arbitrary comparison layouts and model-aware diffs, N4 interactive/offline
HTML and browser PDF export, and three-level drawing are not implemented.
Global crossing minimization and dense-event layout quality are not certified.
1000-tip geometry is a separate check from readable 1000-tip vector output or
interactive performance. See [the stage readout](../validation/RECONCILIATION_N1_N2.md) for actual checks and limits.
