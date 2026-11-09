# VertLife reference trees

Select `--vertlife-dataset amphibians-2018-consensus` to use the Jetz and Pyron
reference, a time-scaled consensus of 10,000 sampled trees. It includes 7,238
amphibian species, including species placed using taxonomic information, and the
outgroup `Homo_sapiens`. See the [study](https://doi.org/10.1038/s41559-018-0515-5),
[data citation](https://doi.org/10.5061/dryad.cc3n6j5), and
[VertLife downloads](https://vertlife.org/files_20170703/).

NWKIT downloads the registered file on first use and verifies its SHA-256 on
every load. The default location is
`~/.cache/nwkit/vertlife/<dataset-id>/`. Use `--download-dir PATH` to choose
another cache, or `--vertlife-file PATH` to read a local copy with the same
checksum. Concurrent first uses share a single download. Subsequent runs can
use the cache offline. If verification fails, provide a verified local copy
or remove the corrupt cache file and retry.

The commands accept labels with spaces or underscores and UTF-8 species lists
with an optional BOM. Exact reference names take precedence over species
parsing. Use `--species-parser`,
`--species-regex`, or a `--species-map-tsv` with an explicit `taxonomy_query`
or `species_label` for sequence labels and alternate names. An explicit map
entry takes precedence over exact matching. The reference is restricted to
matched input species, so the outgroup appears only when requested.

## Topology and dated branch lengths with `constrain`

```sh
nwkit constrain --backbone vertlife --vertlife-dataset amphibians-2018-consensus \
  --species-list species.txt --branch-length original --report matches.tsv \
  -o species.dated.nwk
```

The default `--branch-length none` produces a topology constraint.
`--branch-length original` retains durations in millions of years, sums edges
that become a single branch, and removes the stem above the retained species'
most recent common ancestor. It preserves the reference's branch precision and
any original deviations from exact ultrametricity. An explicit `--outformat`
in dated mode must retain internal and terminal lengths.

Dated mode requires at least two species and one label per species. Use
`--collapse yes` to emit one tip per matched species with its reference name.
In topology mode, repeated labels for one species can form a polytomy at that
species' position. The optional report has one row per input with `leaf_name`,
`taxonomy_query`, `vertlife_leaf_name`, `status`, `dataset`, `source_doi`, and
`source_sha256`. The matching statuses are `exact`, `parsed`, `mapped`,
`unmatched`, and `duplicate_species`. Unmatched labels and unresolved
duplicates in dated mode stop tree generation and still write the report.

## Root placement with `root`

```sh
nwkit root -i gene.nwk --method taxonomy --taxonomy-source vertlife \
  --vertlife-dataset amphibians-2018-consensus -o rooted.nwk
```

The root of the pruned VertLife reference supplies the split to locate in the
input topology. At least two distinct species must match. Multiple tips for a
species must occupy one side of an edge in the unrooted input topology.
If some tips are unmatched, the resolved split must identify a unique input
edge. Branch lengths are written at full precision on the input tree's scale.
Specify VertLife explicitly in the source list to enable it. The default
source order remains `ncbi,opentree,timetree`.

## Divergence age reports with `mcmctree`

```sh
nwkit mcmctree -i species.nwk --vertlife-dataset amphibians-2018-consensus \
  --report node-ages.tsv -o species.copy.nwk
```

This mode requires a rooted input tree with at least two species and one tip
per matched species. Rooted polytomies are supported. It writes the input tree
and a report with one row per internal node. `branch_id` is the input's
zero-based level-order node index, as
in `nwk2table`. A match requires the same descendant species and the same
partition among the node's children after subsetting the reference.
`clade_absent` and `child_clades_differ` rows have empty age fields.

The report's `age_ma` is the mean of the ages implied by the matched reference
node's child branches, calculated from the extant tips upward. The shortest
and longest paths to retained tips appear as `min_tip_path_myr` and
`max_tip_path_myr`, along with `dataset`, `source_doi`, and `source_sha256`.
The path spread describes the source tree's deviation from equal
root-to-tip times. Quantifying uncertainty across trees requires the
underlying sample of dated trees. This mode leaves the topology, branch
lengths, and existing node annotations intact.
