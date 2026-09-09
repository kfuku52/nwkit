# Reconciliation and dating figures

`nwkit reconcile` and `nwkit radte` can publish a PDF, SVG, or PNG report with
`--figure-out`. `nwkit draw` can render the same reports from saved results
without rerunning reconciliation, optimization, or interval estimation. The
reports use the rectangular coordinate engine shared with ordinary tree drawing.

## Reconciliation

```sh
nwkit reconcile --infile examples/radte/gene.nwk \
  --species-tree examples/radte/species.nwk \
  --species-map-tsv examples/radte/species-map.tsv \
  --outfile output/events.tsv --figure-out output/reconciliation.pdf

# Reuse the table; child order and internal node names may differ.
nwkit draw --infile examples/radte/gene.nwk \
  --species-tree examples/radte/species.nwk \
  --reconciliation output/events.tsv --outfile output/reconciliation.svg
```

The gene tree shows the recorded event types using both colors and shapes.
Internal labels have the form `gene node / mapped species node`; the second
tree identifies the species nodes. Long or duplicate internal names receive
short display identifiers. Complete tip names are retained and wrapped.
Unresolved events and unmapped species placements are marked, not silently
reclassified. The view is topological: it does not compare gene substitution
lengths with species times, and it does not invent loss events absent from the
input table.

Saved tables must cover exactly one supplied rooted gene tree. Node matching
uses clade IDs and parent clades, not traversal indices or internal node labels.
Species-event IDs must belong to the supplied species topology. An incompatible
tree or table fails before replacing the figure.

## RADTE

```sh
nwkit radte --gene-tree examples/radte/gene.nwk \
  --species-tree examples/radte/species.nwk \
  --species-map-tsv examples/radte/species-map.tsv \
  --reconcile lca --max-age 30 --rate-sd 0.3 --uncertainty profile \
  --out-prefix output/family --figure-out output/dating.pdf

# No original gene tree or alignment is needed for redrawing.
nwkit draw --radte-prefix output/family \
  --species-tree examples/radte/species.nwk --outfile output/dating.svg
```

The dated gene tree and species tree use the same age axis, with older ages
on the left. Dashed guides identify shared speciation ages. Species-tree
positions use the saved estimated species ages, including species ages fitted
within supplied calibration bounds. Outlined diamonds distinguish fixed hard
calibrations from estimated ages. PAML soft priors are not labeled fixed.
When PAML leaves species nodes unestimated, gray species branches show the input
chronogram; estimated ages and intervals are overlaid on the matching node rows.
Unestimated internal nodes are explicitly labeled `input only`, and their
missing estimates remain missing in the numerical tables.

Age intervals come directly from `.nodes.tsv`. Their estimator and level are
reported beside the tree; native conditional intervals are not labeled MCMC
posterior intervals. When intervals were not requested or are unavailable, the
report shows point estimates without making up uncertainty bars. The lower
panel shows duplication ages, their original calibration ranges, and boundary
estimates, alongside the saved inference diagnostics. PAML ranges are explicitly
identified as original ranges, since not every range is a selected PAML prior.

Units default to the input time units. Use `--branch-length-unit Ma` only when
the input chronogram is in millions of years; this changes a display label,
not any ages or branch lengths.

Redrawing reads `.dated.nwk`, `.nodes.tsv`, `.species.tsv`, and `.manifest.json`
from the prefix. The three numerical inputs must match the complete run
manifest's hashes. Gene lengths, parent clades, shared ages, interval endpoints,
species topology, and original species ages are also validated. The supplied
species tree can be reordered, but cannot introduce a different calibration
input. All numerical members of the bundle are protected against output aliases.
`--audit` includes the bundle input files and never consumes implicit STDIN in
this mode.

## Layout, files, and failure behavior

- The suffix of `--figure-out` selects PDF, SVG, or PNG. Saved-result drawing
  also accepts the existing `--image-format` option.
- `--figure-width`, `--figure-height`, and `--font-size` control report sizing.
  Automatic height expands for wrapped tip names, duplication rows, and
  diagnostics. Reports require at least 8 inches of width and sufficient
  height for their contents.
- Result reports have a fixed panel arrangement. Single-tree layout, trait,
  image, and annotation options are not report options; incompatible nondefault
  options such as `--layout radial` are rejected.
- Inference commands stage a requested figure together with their numerical
  outputs. A drawing error or a handled publication failure preserves the
  previous complete outputs. RADTE records the requested figure's hash in its
  run manifest. A separately redrawn figure does not change the numerical bundle.
- Ordinary `nwkit draw` without a result-table option retains its existing
  behavior and default width.

See [the reproducible figure examples](examples/radte/visualization/README.md)
for the small shared-age illustration and the bundled GeneRax example.

For external species intervals, separate hard-range and fitted-interval layers,
and three-mode comparison reports, see
[Species-age uncertainty](RADTE_SPECIES_UNCERTAINTY.md).
