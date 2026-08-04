![](https://raw.githubusercontent.com/kfuku52/nwkit/master/logo/logo_nwkit_large.png)

[![Run Tests](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/tag/kfuku52/nwkit?label=release)](https://github.com/kfuku52/nwkit/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/kfuku52/nwkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview
**NWKIT** ([/njuːkit/](https://ipa-reader.com/?text=nju%CB%90kit&voice=Joanna)) is a toolkit for manipulating phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format).

## Installation

Packaged releases of NWKIT are available from [Bioconda](https://anaconda.org/bioconda/nwkit). The repository can be ahead of the current Bioconda package; compare the release and Bioconda badges above when selecting a version. For users requiring a `conda` installation, please refer to [Miniforge](https://github.com/conda-forge/miniforge) for a lightweight conda environment.

#### Install from Bioconda

```
conda install bioconda::nwkit
```

#### Verify the installation by displaying the available options

```
nwkit -h
```

#### (For advanced users) Install the development version from GitHub

```
pip install git+https://github.com/kfuku52/nwkit
```

NWKIT requires Python 3.10 or newer.

#### Optional SVG rasterization support

Raster-image validation and processing use Pillow, which is installed with
NWKIT. Install the optional CairoSVG integration to normalize, trim, or
resize SVG files with `nwkit image`, or to render SVG tip images with
`nwkit draw`:

```
pip install "nwkit[image]"
```

## Subcommands
See [Wiki](https://github.com/kfuku52/nwkit/wiki) for usage.

Shared option naming, standard-input rules, TSV schemas, missing-value policy,
and output-column vocabulary are defined in
[CLI and TSV conventions](https://github.com/kfuku52/nwkit/blob/master/CLI_TSV_CONVENTIONS.md).

- [`annotate`](https://github.com/kfuku52/nwkit/wiki/nwkit-annotate): Attaching tip-table values and aggregating them as Newick properties
- [`asr`](https://github.com/kfuku52/nwkit/wiki/nwkit-asr): Inferring categorical ancestral states and imputing missing tip states under an Mk model
- [`constrain`](https://github.com/kfuku52/nwkit/wiki/nwkit-constrain): Generating a species-tree-like Newick file for topological constraint
- [`collapse`](https://github.com/kfuku52/nwkit/wiki/nwkit-collapse): Collapsing internal branches by support and/or branch length
- [`compose`](https://github.com/kfuku52/nwkit/wiki/nwkit-compose): Assembling compatible roots, values, and annotations from multiple trees
- [`cladefreq`](https://github.com/kfuku52/nwkit/wiki/nwkit-cladefreq): Summarizing clade frequencies across a tree collection
- [`consensus`](https://github.com/kfuku52/nwkit/wiki/nwkit-consensus): Generating a consensus tree or transferring consensus support to a reference tree
- [`diff`](https://github.com/kfuku52/nwkit/wiki/nwkit-diff): Reporting interpretable clade, root, value, and annotation differences between trees
- [`dist`](https://github.com/kfuku52/nwkit/wiki/nwkit-dist): Comparing tree topology and branch lengths with multiple distance metrics
- [`draw`](https://github.com/kfuku52/nwkit/wiki/nwkit-draw): Drawing phylogenetic trees with Cartesian, polar, equal-angle unrooted, spiral, or recursive fractal geometry; layouts can allocate space from measured labels and annotations, and compatible rooted layouts can use tidy subtree packing
- [`drop`](https://github.com/kfuku52/nwkit/wiki/nwkit-drop): Removing node and branch information
- [`image`](https://github.com/kfuku52/nwkit/wiki/nwkit-image): Retrieving representative species images with license-aware filtering
- [`info`](https://github.com/kfuku52/nwkit/wiki/nwkit-info): Printing tree information
- [`intersection`](https://github.com/kfuku52/nwkit/wiki/nwkit-intersection): Dropping non-overlapping leaves/sequences between two trees or between a tree and an alignment
- [`label`](https://github.com/kfuku52/nwkit/wiki/nwkit-label): Adding unique node labels
- [`mark`](https://github.com/kfuku52/nwkit/wiki/nwkit-mark): Adding texts to node labels by identifying the targets with a leaf name regex
- [`mcmctree`](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree): Introducing divergence time constraints for PAML's mcmctree
- [`monophyly`](https://github.com/kfuku52/nwkit/wiki/nwkit-monophyly): Assessing whether species or trait-defined groups are monophyletic
- [`nhx2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-nhx2nwk): Generating Newick from NHX
- [`nwk2table`](https://github.com/kfuku52/nwkit/wiki/nwkit-nwk2table): Converting a Newick tree into a parent-child table
- [`printlabel`](https://github.com/kfuku52/nwkit/wiki/nwkit-printlabel): Searching and printing node labels
- [`prune`](https://github.com/kfuku52/nwkit/wiki/nwkit-prune): Pruning leaves
- [`rename`](https://github.com/kfuku52/nwkit/wiki/nwkit-rename): Renaming nodes using a TSV mapping or regular expression
- [`rescale`](https://github.com/kfuku52/nwkit/wiki/nwkit-rescale): Rescale branch length with a given factor
- [`root`](https://github.com/kfuku52/nwkit/wiki/nwkit-root): Placing or transferring the tree root
- [`sanitize`](https://github.com/kfuku52/nwkit/wiki/nwkit-sanitize): Eliminating non-standard Newick flavors
- [`sample`](https://github.com/kfuku52/nwkit/wiki/nwkit-sample): Selecting a representative leaf subset by filters, ranks, and sampling method
- [`shuffle`](https://github.com/kfuku52/nwkit/wiki/nwkit-shuffle): Shuffling branches and/or labels
- [`skim`](https://github.com/kfuku52/nwkit/wiki/nwkit-skim): Sampling leaves from clades with shared traits
- [`subtree`](https://github.com/kfuku52/nwkit/wiki/nwkit-subtree): Generating a subtree Newick file
- [`table2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-table2nwk): Converting a parent-child table into a Newick tree
- [`transfer`](https://github.com/kfuku52/nwkit/wiki/nwkit-transfer): Transferring information between trees
- [`validate`](https://github.com/kfuku52/nwkit/wiki/nwkit-validate): Validating one or more Newick trees and reporting structural issues

## Tree drawing layouts

`nwkit draw --layout` provides eight deterministic geometries from the same
Newick input. Subtree packing is selected independently:

```sh
nwkit draw -i tree.nwk --layout rectangular -o rectangular.svg
nwkit draw -i tree.nwk --layout slanted     -o slanted.svg
nwkit draw -i tree.nwk --layout cladogram   -o cladogram.svg
nwkit draw -i tree.nwk --layout circular   -o circular.svg
nwkit draw -i tree.nwk --layout radial     -o radial.svg
nwkit draw -i tree.nwk --layout unrooted   -o unrooted.svg
nwkit draw -i tree.nwk --layout spiral     -o spiral.svg
nwkit draw -i tree.nwk --layout fractal    -o fractal.svg

# Restrict circular geometry to an upper semicircle (a fan drawing).
nwkit draw -i tree.nwk --layout circular --angular-span 180 \
  -o circular-sector.svg

# Compact the same rectangular geometry without changing topology or depth.
nwkit draw -i tree.nwk --layout rectangular --subtree-packing tidy \
  -o rectangular-tidy.svg
```

The 10-versus-100-tip gallery used in the Wiki is reproducible with
`python scripts/make_draw_layout_gallery.py --output-prefix PATH`; it writes
both PNG and editable SVG output at 8 pt Helvetica.

- `rectangular` is the conventional orthogonal phylogram and remains the
  default.
- `slanted` retains root-to-node branch-length depth but connects each parent
  and child directly.
- `cladogram` ignores branch lengths, aligns all tips, and displays topology
  with slanted branches.
- `circular` maps branch-length depth to radius and uses circular connectors
  followed by radial branch segments. `--angular-span` restricts it to a
  sector; 180 degrees draws the conventional fan in the upper half-plane.
- `radial` maps the rooted phylogram to polar coordinates and connects nodes
  with straight branches. It accepts the same angular-sector controls.
- `unrooted` suppresses a degree-two input root and uses a central-node,
  branch-length-aware equal-angle arrangement. Add
  `--unrooted-method equal-daylight` to iteratively distribute open angular
  space around internal nodes while preserving every edge length. Each pass is
  checked for branch crossings; unsafe rotations are reduced or rejected. The
  pass limit is controlled by `--daylight-iterations`. Equal-daylight is
  intended for detailed views of up to 2,000 displayed nodes; use equal-angle
  or collapse the drawing above that limit.
- `spiral` wraps the orthogonal phylogram around an Archimedean spiral. Use
  `--spiral-turns FLOAT` to override the tip-count-dependent turn count. The
  warped geometry uses a non-overlapping radial band narrower than the pitch
  between turns, retaining topology and relative depth without folding at the
  inner turn. Euclidean distances in the finished image are not branch lengths.
- `fractal` recursively assigns descendant-rich clades more space and fits the
  resulting crossing-free radial-fractal drawing to the requested rectangle. It is a
  topology-and-balance view: branch lengths intentionally do not determine
  geometry.

Subtree packing is independent of geometry. The default
`--subtree-packing standard` uses conventional tip rows or angles.
`--subtree-packing tidy` uses a branch-length-aware implementation of the
non-layered tidy tree algorithm of
[van der Ploeg (2014)](https://doi.org/10.1002/spe.2213), as introduced for
phylogenetic trees by
[de Vienne (2022)](https://doi.org/10.1093/molbev/msac204). It preserves
topology and root-to-node depth while moving non-overlapping subtrees closer
together; physical tip coordinates may differ from traversal order because
branches and terminal labels at disjoint depths can pass one another. Tidy
packing composes with `rectangular`, `circular`, and `spiral`; other geometries
reject it because their straight or recursively placed connectors do not yet
have the same non-overlap guarantee. Complete circles place their seam outside
the extrema of the compacted coordinates, retaining a one-to-one, planar polar
transform.

Angular extent is independent of circular versus straight radial connectors.
`--angular-span` accepts values greater than zero through 360 degrees for
`circular` and `radial`. `--angular-center` rotates the sector using 0 degrees
for right, 90 for up, 180 for left, and 270 for down. The defaults are 360 and
90 degrees; consequently, a 180-degree sector occupies the upper half-plane
from right to left. Spiral, unrooted, fractal, and Cartesian layouts reject
non-default angular-sector settings rather than silently changing their
meaning. When `--figure-height` is omitted, circular and radial sectors derive
their height from the occupied angular bounds and measured label extents; an
upper semicircle therefore uses approximately half the height of a complete
circle instead of leaving an unused lower half-plane.

Tip spacing is independent of tree layout. The default
`--tip-spacing uniform` gives each tip the same row or angular allocation.
`--tip-spacing label-aware` instead uses the measured height of wrapped labels,
badges, tracks, and tip images. It varies rows in Cartesian layouts, boxes used
by tidy packing, and angular sectors in circular, radial, unrooted, spiral,
and fractal layouts. Thus a label-aware phylogram is requested directly, for
example, as `--layout circular --tip-spacing label-aware`; it is not a separate
layout type.

Tip labels can be wrapped for display without changing names in the Newick
tree:

```sh
nwkit draw -i tree.nwk --layout circular --tip-spacing label-aware \
  --tip-label-wrap auto -o circular.svg
nwkit draw -i tree.nwk --layout rectangular --tip-spacing label-aware \
  --tip-label-wrap taxonomy -o taxa.svg
nwkit draw -i tree.nwk --layout fractal --tip-spacing label-aware \
  --tip-label-wrap 12 -o wrapped.svg
```

`--tip-label-wrap auto` measures wrapping targets corresponding to one through
five lines with the selected Matplotlib font and wraps an individual label
only when the reduction in layout congestion is meaningful. A positive integer sets the
maximum characters per displayed line; `taxonomy` keeps a conservatively
recognized underscore-delimited `Genus_species` prefix together. Wrapping
prefers whitespace and
`_`, `-`, `/`, or `|` boundaries before making a hard break. The underlying
tip name is never modified. `--tip-label-font-style taxonomy` italicizes only
labels that consist exactly of such a binomial; `italic` applies to every tip.

`--tip-label-position auto` aligns labels in a standard rectangular layout and
puts them at branch ends in other geometries and in tidy-packed rectangular
drawings. Branch-end labels are included in tidy collision geometry. Use
`--tip-labels no` for dense overview
figures where individual names would be illegible. Tip-image columns currently
apply to `rectangular` (with either packing strategy), `slanted`, and
`cladogram`; node markers,
support, pies, property labels, and categorical tip styling follow nodes in
every layout.

Branch support and ordinary Newick node labels can be selected independently.
A parsed node label is available through the `name` property, so the same
property-label controls also handle named clades:

```sh
nwkit draw -i named-and-supported.nhx -o annotated.svg \
  --support-labels yes --support-min 0.9 \
  --node-label-property name --node-label-target root,intnode
```

`--node-label-filter`, `--node-label-decimals`, and `--node-label-prefix`
provide selective formatting for any numeric or textual Newick/NHX property.
When node pies and labels are combined, their rendered extents are included in
fitting and collision reporting.

### Annotation-aware and auditable drawing

Rendered label extents, badges, support labels, node-property labels, node pies,
legends, tip tracks, and branches participate in a deterministic post-render
collision check. The default `--collision-policy resolve` moves annotations
that can be moved safely and retains the best non-worsening placement;
`warn`, `error`, and `ignore` make the remaining-collision policy explicit. A
machine-readable report records the layout, figure dimensions,
font settings, branch-length encoding, equal-daylight convergence, wrapping,
automatic collapsing, figure-boundary overflow, annotation occupancy, and a
complete branch/annotation collision audit. It also counts proper branch-pair
crossings for drawings within the bounded crossing-audit size and reports
explicitly when a larger drawing was skipped:

```sh
nwkit draw -i tree.nhx -o annotated.svg \
  --layout circular --tip-spacing label-aware \
  --tip-track habitat --tip-track temperature \
  --branch-color-property regime \
  --branch-width-property posterior \
  --layout-report annotated.layout.json
```

Tip tracks read categorical or numeric Newick/NHX properties; `auto` infers
the shared type separately for each property. Continuous tracks use
`--tip-track-palette` (default `viridis`), while categorical tracks use
`--trait-palette` and `--property-color` overrides. Branch color and width can
likewise be mapped from properties, with `--branch-width-range MIN,MAX` setting
the rendered width interval. Legends automatically use multiple columns and
move to the right when they become too dense; the position and column count
can also be fixed explicitly.

An exact scale can be added to rectangular, circular, and unrooted
drawings:

```sh
nwkit draw -i tree.nwk -o tree.svg \
  --scale-bar auto --branch-length-unit substitutions/site
```

The scale label is placed above the bar. A dedicated strip below the tree
panel keeps both elements separate from branches, tip labels, and annotation
layers; fixed-height figures reserve this space from the requested height.

Topology-only or geometrically warped layouts reject a scale bar rather than
implying a false distance interpretation. Slanted and straight radial layouts
retain root-to-node depth only as a projection, so they reject a bar that could
be mistaken for the length of a visible branch. Use `--depth-guide` instead:

```sh
nwkit draw -i tree.nwk --layout slanted --depth-guide auto \
  --branch-length-unit substitutions/site -o slanted.svg
nwkit draw -i tree.nwk --layout radial --depth-guide 0.5 \
  --branch-length-unit substitutions/site -o radial.svg
nwkit draw -i tree.nwk --layout spiral --depth-guide auto \
  --branch-length-unit substitutions/site -o spiral.svg
```

`slanted` receives a horizontal depth axis with vertical grid lines; `radial`
receives concentric rings for a complete circle or matching concentric arcs for
an angular sector, labelled in the largest empty direction within the displayed
sector; and
`spiral` receives a linear root-to-node distance key for the depth encoded
across its spiral band. These guides describe cumulative root-to-node distance,
not the Euclidean length of a slanted, radial, or warped edge. `auto` chooses a
readable interval and a positive number sets it
explicitly. On a dense radial tree with no safe sector for ring numbers, the
numbers are omitted and the interval remains explicit in the lower guide
caption. Both scale bars and depth guides reserve a lower annotation strip.
Cladogram and fractal layouts encode topology rather than branch-length depth
and therefore support neither guide. A requested scale or guide interval cannot
exceed the displayed tree-depth span.

For very large overview trees,
`--max-visible-tips INT` collapses clades in a drawing-only copy and marks the
resulting pseudo-tips; `--collapse-label` customizes their labels. The original
tree file is never changed, and the report lists every collapsed clade.
Collapsed properties are retained only when present and identical in every
descendant tip; incomplete and heterogeneous properties become `partial` and
`mixed`, respectively. `--collapse-property-aggregation mean` explicitly
requests a mean for complete numeric properties. The report records the
status and value of each summarized property.

## Drawing trees with species images

`nwkit image` writes a tip-keyed `manifest.tsv` that `nwkit draw` can consume
without contacting external services again. When `--max-per-species` writes
multiple ranked candidates for a tip, `nwkit draw` uses the first manifest row
for that tip:

```sh
nwkit image -i tree.nwk \
  --out-dir species-images \
  --style silhouette \
  --max-per-species 1 \
  --output-format png \
  --max-edge 256 \
  --canvas square \
  --background transparent \
  --trim transparent

nwkit draw -i tree.nwk \
  --tip-image-manifest species-images/manifest.tsv \
  --tip-image-size 18 \
  --tip-image-gap 4 \
  --figure-width 5 \
  -o tree-with-images.svg
```

For PhyloPic, `nwkit image` prefers the taxon's curated `primaryImage`. If that
image is unavailable or excluded by the license policy, candidates at the same
taxonomic rank are ranked by license openness, vector availability, drawable
aspect ratio, and resolution before genus or family fallbacks are considered.
The manifest records `is_primary`, `is_vector`, and `selection_reason` for an
auditable choice. This is a curated-representative policy, not a download-count
ranking; the public PhyloPic API does not expose per-image download counts.

Relative `local_path` values are resolved from the manifest directory.
Use `--tip-image-root PATH` when a manifest has been moved independently from
its image directory. Missing tip rows follow `--unmatched warn|error|ignore`;
for duplicated tip rows, the first row is used deterministically. Broken local
paths are rejected. The generated figure embeds the selected images, while
license and creator information remains in the `ATTRIBUTION.md` written by
`nwkit image`; distribute that file with the figure when its licenses require
attribution.

## Tree comparison, composition, and provenance

`diff` reports how two trees differ before information is combined. `compose`
then assembles compatible components from named sources and records every match
or conflict in a TSV report:

```sh
nwkit diff -i topology.nwk -i2 bootstrap.nwk -o differences.tsv
nwkit compose -i topology.nwk \
  --root-source rooted.nwk \
  --support-source bootstrap.nwk \
  --length-source chronogram.nwk \
  --property-source habitat=state@annotated.nhx \
  --report composition.tsv \
  -o combined.nwk
```

`transfer`, `compose`, `diff`, and transferred roots accept
`--taxon-mode intersection`. Reports distinguish `exact_match` from
`projected_match`: the latter is unique only in the trees induced by their
shared tips and does not establish that the original branches are identical.
Ambiguous projections are reported and left unchanged. Projected node names
and NHX properties can be transferred under the default `compatible-only`
policy, whereas projected support values and branch lengths require an explicit
`--allow-projected-values yes`. Strict policy accepts exact matches only.

Matching defaults to rooted descendant clades. For edge-associated values in
trees with different rootings, `--match-basis split` instead uses canonical,
root-independent splits:

```sh
nwkit transfer -i target.nwk -i2 rerooted-source.nwk \
  --support yes --match-basis split --report transfer.tsv \
  -o supported.nwk
```

Root alignment preserves support values, internal names, and NHX properties on
their original unrooted branch splits. A bifurcating root represents one
unrooted edge as two child branches, so `transfer` and `compose` resolve that
representation explicitly. By default, one source annotation is copied to
both target halves, equal source annotations are treated as one value, and
conflicting annotations follow the source half with the same projected
descendant taxa. Branch lengths use the source edge total and preserve the
target root-position ratio; a zero-length target root edge is split equally.
No result depends on traversal order.

Use repeatable `--root-edge-policy TARGET_PROPERTY=POLICY` options to override
the default per target property. `skip` leaves the edge unchanged,
`equal-only` accepts only a unique or equal annotation, `matching-side` follows
projected descendant taxa, `mean`, `min`, and `max` reduce numeric annotations,
and `edge-total` applies the branch-length rule. `auto` is the annotation
default; `edge-total` is the length default. `*` can set a fallback policy.
For example:

```sh
nwkit compose -i topology.nwk \
  --support-source bootstrap.nwk \
  --length-source chronogram.nwk \
  --match-basis split \
  --root-edge-policy support=mean \
  --root-edge-policy length=edge-total \
  --report composition.tsv \
  -o combined.nwk
```

A composition manifest accepts a top-level `root_edge_policies` object, and a
custom `properties` entry can contain `root_edge_policy`. Precedence is
top-level manifest, property entry, then CLI. Reports record the policy,
resolution, candidate counts, and all candidate values. Projected support and
length values still require `--allow-projected-values yes`, and strict mode
still rejects every projected match.

`compose` changes the target rooting only when `--root-source` is supplied.
The `midpoint`, `outgroup`, `mad`, `mv`, `taxonomy`, and `transfer` rooting
methods apply the same split-based annotation preservation. For tree
collections, `validate --require-same-rooting yes` compares the root
bipartition when leaf sets are identical.

Arbitrary NHX properties can be transferred or renamed directly:

```sh
nwkit transfer -i target.nwk -i2 source.nhx \
  --property color --property-map posterior=source_posterior \
  --report transfer.tsv -o annotated.nhx
```

`annotate` joins a TSV keyed by `leaf_name` to tree tips and can aggregate
values onto ancestral clades:

```sh
nwkit annotate -i tree.nwk --table traits.tsv \
  --columns habitat \
  --aggregate habitat:unique:shared_habitat \
  --aggregate body_mass:mean:mean_body_mass \
  -o annotated.nhx
```

Every functional command accepts `--audit PATH`. Each invocation appends one
JSON Lines record containing NWKIT version and arguments, input and output
SHA-256 hashes, input interpretation, random seeds, external-data settings,
warnings, runtime status, and captured command messages:

```sh
nwkit sanitize -i raw.nwk --audit workflow.audit.jsonl |
  nwkit root --method midpoint --audit workflow.audit.jsonl > rooted.nwk
```

## Citation
There is no published paper on NWKIT itself, but we used and cited NWKIT in several papers including [Fukushima & Pollock (2023, Nat Ecol Evol 7: 155-170)](https://www.nature.com/articles/s41559-022-01932-7).

## Development

Install the development and optional image dependencies, then run the same checks used by CI:

```
pip install -e ".[dev,image]"
ruff check nwkit tests setup.py
pytest tests/ -q
python -m build
```

See [CHANGELOG.md](https://github.com/kfuku52/nwkit/blob/master/CHANGELOG.md) for changes and
[RELEASING.md](https://github.com/kfuku52/nwkit/blob/master/RELEASING.md) for the release checklist.

# Licensing
This program is MIT-licensed. See
[LICENSE](https://github.com/kfuku52/nwkit/blob/master/LICENSE) for details.
