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
placement. N3 renders source/target trees and highlights the recorded consumed clades. Arbitrary other transformation models are preserved as extension
blocks but do not acquire X1 semantics.

View and comparison schemas are executable N3/N4 contracts. Use `--view-config`,
`--view-out` and `--comparison` with `nwkit draw`. Comparison documents have
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
supports standalone trees, a shared upper tree with multiple lower histories,
and the three-layer chains described below. Broader hierarchies use adjacent-level
panels. Unmapped/incomplete placements and expectation maps fail explicitly.
Upper and lower branch lengths are retained but not mixed as a time scale.

Species corridors have visible outlines, aligned tips and one uniform width
throughout each upper tree, sufficient for its widest occupied lane envelope.
Both daughters turn at the same parental x coordinate, including offset gene
lanes and duplications within a species branch. Concurrent copy lineages use
slightly staggered x coordinates at a shared species junction, keeping their
vertical fork segments separate inside its corridor. These offsets are display
spacing, not differences in speciation time. An ancestor and its descendant
copies share lanes;
serializing a duplication node does not count it as a third simultaneous copy.
This remains a topology layout, including when species lengths are available.

An aioi representative's unary speciation can record a marginalized terminal
side without a placed loss event. When its recorded undated model has D=T=0,
positive loss propensity, detection=1 for that terminal species, and no observed
copy in either mapping or recorded leaf placements there, draw shows a short
copy branch ending in a loss cross. This marks required terminal extinction
under those conditions. Its midpoint is a
display convention, with no inferred event time or posterior support. The source
event list and integrated E history are unchanged; the HTML primitive's details
identify this display interpretation. Unknown/incomplete conditions, incomplete
detection, internal empty sides and positive D/T do not acquire a loss cross.
Other recorded marginalized terminal sides use a separate empty-side stub ending
in a gray diamond (E), preserving the distinction from a placed loss. These sides
have independent lanes even when a transferred observed copy reaches the same
species. The source event list and support fields remain unchanged. Figure
titles, legends, footers and continuation labels honor the requested base font
size along with ordinary labels.

Duplication nodes are red squares, speciation nodes blue circles, and losses red
crosses. Transfers use an orange hollow circle at the donor and an orange
directed arc with an arrowhead at the recipient. Imported explicit transfer-entry
markers remain orange; X1 consumed branches stay purple, separate from transfer.

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
Export bounds include corridor strokes, transfer curves and wrapped labels;
the renderer expands and translates a display copy when necessary. The source
scene and its recorded events remain unchanged. Generated level/transformation
headings use ASCII text so fonts such as Helvetica do not need an arrow glyph.

Figure and normalized-document writes are staged as a recoverable bundle.
Validation/write failures preserve previous files. Input/output aliases,
including audit and normalized-output paths, are rejected. Document drawing
does not consume the default implicit stdin; actual inputs appear in `--audit`.

## Remaining stages

N3 arbitrary comparison manifests, compatible-model event differences, H/G/X1
views, N4 offline HTML and PDF downloads, and recorded three-layer nesting are
implemented. Larger branching hierarchies use adjacent-level overview panels.
Global crossing minimization and dense-event layout quality are not certified.
1000-tip geometry is a separate check from readable 1000-tip vector output or
interactive performance. See [the stage readout](../validation/RECONCILIATION_N1_N2.md) for actual checks and limits.


## Executable comparison and display state (N3)

```sh
nwkit draw --comparison comparison.json --outfile comparison.pdf
nwkit draw --document history.json --view paired --outfile paired.svg
nwkit draw --document history.json --view-config saved.view.json \
  --view-out effective.view.json --outfile revised.pdf
```

A comparison contains any number of panels. Each document reference resolves
relative to the manifest and must match its SHA-256. Inference and history
statuses/timing dates/conditions remain distinct; missing, failed, timeout and
unsupported panels stay visible. An inferred standalone tree can be present
while its history is missing. Event differences are reported only for matching
history kind, selection and model conditions, and include recipient placement,
recorded observed arity and the marginalized side when supplied.
A fixed-gene-tree comparison also verifies actual G topology, lengths and units,
in addition to the recorded shared context hashes.

Display state holds layout, explicit leaf-ID order, collapsed clades, event
colors, font and output profile. A comparison applies a stable clade selection
where that clade exists; absent panels retain their context. If an opaque tree/node
ID names different clades in different panels, a shared collapse is rejected. Explicit CLI font
and PDF page settings override saved values. The output path determines export
format. Scientific trees, events, parameters and original summaries remain
unchanged; `--document-out` retains the full source rather than a pruned history.
Zoom/pan are browser viewport controls; their position is not part of this view
schema. Collapsing an upper context in an embedded view is not supported.

H and G are separate recorded trees, with X1 consumed clades highlighted and
source/target correspondence retained. No gene substitution length becomes an
H duration. A three-layer chain nests lower histories in the recorded middle
incoming corridors and those corridors in the outer tree. Middle transfer
curves are retained; lower transfers use their explicit donor/recipient. Arc
length and lane positions are display coordinates, not inferred time. Other
hierarchies use linked adjacent-level panels with unchanged middle-tree IDs.
Nested histories share the two-layer frontier-lane calculation, including
multiline label clearance, staggered copy forks and marginalized terminal-side
stubs. Rectilinear host corridors preserve right angles through lane offsets;
recorded curved hosts retain their curves. Separate label columns reserve room
for wide wrapped text across all three layers.
Expectation maps require a separate quantitative representation: the renderer
does not fabricate events from expected counts.

## Self-contained interactive export (N4)

```sh
nwkit draw --document history.json --outfile history.html
nwkit draw --comparison comparison.json --outfile comparison.html
```

HTML contains SVG geometry, complete source records, controls and a selectable
vector PDF. It loads no external scripts, fonts, images or services. Search
labels/IDs/event types, select branches/events for recorded model details,
toggle clades, zoom, pan, fit the overview, and save display state. Downloading
the embedded PDF preserves the initial scientific display/pagination; printing
uses the current browser viewport and display selection. Saved view JSON can be
rendered by the CLI for vector output with the same font/collapse choices.
Single-source export downloads canonical JSON. Comparison source export
downloads a single ZIP containing its manifest and original document bytes, with
adjusted relative filenames and unchanged document hashes.

Unreported support remains unreported. Imported raw support values do not
become posterior probabilities. In native aioi representative documents,
nonobserved sides are integrated in E; they are not placed biological losses.
The maximum partially marginalized derivation is not full R MAP or H MAP.

Vector PDF export is native Matplotlib, with selectable text, A4/A3 pagination
and continuation markers. Browser printing is a browser-controlled secondary
export. Scientific fitting/restoration belongs to aioi or an external producer;
nwkit performs no likelihood, root or parameter optimization during drawing.


## Audit corrections (2026-10-15)

Comparison primitives and event coordinates retain their panel and tree identity,
so reused source IDs cannot mix recorded details or overwrite another panel's
positions. The canvas also retains the expanded width required by nested layers.
Every referenced file is read once: its hash, parsed science and downloadable ZIP
use the same byte snapshot. Shared-context hashes remain producer-recorded input
claims; the viewer verifies their equality and the document-byte hashes. It cannot
rehash absent original inference inputs. Fixed-gene mode additionally checks the
actual serialized observed trees.

Identity and X1 transformations must reconstruct the complete target topology,
including parts outside each local block. Malformed structural references and
excessive JSON decoder nesting produce input errors. HTML keeps all geometry when
loading a collapsed view, allowing expand/reset; its embedded PDF retains that
initial collapsed view. Upper context collapse remains disabled in embedded mode.
Tree-only panels with missing histories honor the same display state.

See [the audit evidence](../validation/RECONCILIATION_AUDIT_20261015.md) for executed
checks, scope and limitations. No inference is performed by these corrections.


The [second audit](../validation/RECONCILIATION_AUDIT2_20261015.md) preserves event
colors when highlighting X1 branches, keeps display event coordinates consistent
with clade selection/collapse, and supports deep recPhyloXML output without a
process-wide recursion setting. Imported opaque XML retains expanded namespace,
attribute, text and tail semantics on the deep path; normal formatting is unchanged.

The [drawing polish audit](../validation/RECONCILIATION_POLISH_AUDIT_20261016.md)
covers clipping, wide nested labels, orthogonal copy forks and nested terminal
side annotations, with affected regression and CLI export evidence.
