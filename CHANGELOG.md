# Changelog

All notable changes made after the `v0.21.1` tagged release are tracked here.

## [Unreleased]

## [0.31.0] - 2026-07-21

### Added

- A documented shared contract for CLI spelling, standard input, tip-keyed TSV
  validation, missing values, unmatched rows, branch IDs, node classes, taxon
  sets, and empty-result headers.

### Changed

- Canonical multiword long options now use kebab-case throughout. Historical
  snake_case and other replaced spellings remain compatibility aliases and
  emit a warning naming the canonical replacement when used.
- Tip-keyed tables used by annotation, ASR, drawing, monophyly checks, and
  sampling now share `leaf_name`, duplicate, missing-value, and unmatched-row
  handling through `--missing-values` and `--unmatched`.
- Node reports now use 0-based level-order `branch_id` values, `node_class`, and
  the `*_taxa` / `num_*_taxa` vocabulary. ASR uses `leaf` rather than `tip`,
  including for canonical `--target` values.
- Image-specific mappings are named `--species-name-tsv`; image output paths
  are `--manifest-out` and `--attribution-out`; sample tables use `--report`;
  and skim group tables use `--group-table-prefix`.
- Auxiliary outputs require file paths, and at most one input may read from
  standard input in a command.

### Fixed

- `image` now honors `--species-parser` and `--species-map-tsv`.
- `info` and `printlabel` now honor `--outfile`; non-Newick commands no longer
  advertise an unused `--outformat`.
- `table2nwk` preserves node names such as `NA`, `NaN`, and `null`.
- Empty TSV results retain stable headers.

## [0.30.0] - 2026-07-17

### Added

- `transfer` and `compose` now accept repeatable
  `--root-edge-policy TARGET_PROPERTY=POLICY` options for deterministic
  bifurcating-root ambiguity handling. Available policies are `auto`, `skip`,
  `equal-only`, `matching-side`, `mean`, `min`, `max`, and `edge-total`.
  Composition manifests support both top-level policy mappings and per-property
  overrides.
- Transfer reports now include the selected root-edge policy, resolution,
  target/source candidate counts, and a deterministic JSON record of every
  source candidate value.

### Changed

- Conflicting support, names, and NHX values on the two halves of a source root
  edge now follow the half with matching projected descendant taxa by default.
  Unique and equal values retain their edge-wide behavior. Numeric reducers and
  conservative rejection remain explicitly selectable.
- Ambiguous root-edge branch lengths now use the source edge total while
  retaining the target root-position ratio. If the target edge total is zero,
  the length is divided equally. An explicit `matching-side` policy can instead
  retain source-side lengths.

### Fixed

- Root transfer and the root alignment used by `compose` and `transfer` now
  preserve internal-node support, names, and NHX properties on their canonical
  unrooted branch splits instead of leaving them attached to nodes whose
  descendant taxa change during rerooting.
- `compose` annotation and value sources no longer change the target rooting
  when `--root-source` is absent.
- Failed root transfers no longer leave the target tree partially rerooted.
- All `root` methods now preserve internal support, names, and NHX properties
  on canonical branch splits, preserve root and tip metadata, and leave their
  input tree unchanged. `root` output also retains custom NHX properties.
- Collapsing a singleton root no longer loses names or custom properties.
- `validate --require-same-rooting yes` now compares the actual root
  bipartition for trees with identical leaf sets, rather than only comparing
  rooted/unrooted status.
- When branch lengths are transferred from a differently rooted tree, the
  source now supplies the total length of the aligned root edge while the
  target/root source retains its root-position ratio. This avoids replacing an
  intentional root position with the rerooting library's temporary 1:1 split.
- Split-based length transfer now recognizes the leaf and internal children of
  a root placed on a pendant edge as two halves of one physical edge. This
  prevents half-edge lengths from being lost when target and source rootings
  differ.
- `diff --comparison unrooted` now counts each canonical root edge once and
  selects terminal/internal splits independently of the displayed root, so a
  pure reroot no longer appears as an unrooted topology difference.
- Split-based transfer now treats equal support, names, and NHX values on the
  two children of a bifurcating root as one root-edge value.

## [0.29.0] - 2026-07-17

### Added

- Rooted-clade or root-independent canonical-split matching for `transfer` and
  `compose`, selected with `--match-basis` and recorded in TSV reports.

### Changed

- Partial-taxon correspondences are now reported as projected rather than
  exact matches. Projected support and branch-length transfer is blocked by
  default and requires `--allow-projected-values yes`; strict policy rejects
  all projected matches.
- `intersection` now uses a built-in FASTA reader and writer that preserves
  retained records verbatim; non-FASTA sequence formats are no longer accepted,
  and Biopython is no longer a runtime dependency.

## [0.28.0] - 2026-07-17

### Added

- `compose` for assembling a topology, root, node names, support values, branch
  lengths, and arbitrary NHX properties from compatible source trees, with an
  optional per-clade provenance report and JSON manifest.
- `diff` for rooted-clade or unrooted-split comparisons that report leaf-set,
  root, topology, value, and arbitrary-property differences as TSV.
- `annotate` for joining tip metadata to trees and aggregating categorical or
  numeric values onto internal nodes.
- A shared `--audit` option that appends command arguments, input and output
  hashes, detected input semantics, random seeds, warnings, and runtime status
  as JSON Lines.
- Auditable clade mapping across partially overlapping tip sets through
  `--taxon-mode intersection` for `transfer`, `compose`, `diff`, and transferred
  roots.
- Arbitrary NHX-property transfer and renaming with `--property` and
  `--property-map`, plus strict or compatible-only transfer policies and TSV
  reports.

### Changed

- Root alignment used internally for annotation and branch-value mapping now
  preserves source branch lengths instead of redistributing an already matching
  root edge.
- Named-tree parsing preserves support values explicitly stored as NHX
  properties, allowing names and support to coexist in composed output.

## [0.27.0] - 2026-07-15

### Added

- A reproducible `--seed` option for the `shuffle` and `skim` commands.
- Wheel-content, installed-wheel smoke, lint, and current-Python checks in CI.
- A tag-driven GitHub Release workflow that validates the tag against the package version.
- Standard `nwkit` and `python -m nwkit` entry points with `--version` support.

### Changed

- The minimum supported Python version is now 3.10; CI covers Python 3.10 through 3.14.
- Packaging metadata and command entry points are declared in `pyproject.toml`.
- Runtime modules use explicit imports instead of wildcard imports from `nwkit.util`.
- TimeTree requests and ETE documentation links use HTTPS.

### Fixed

- Test and legacy script packages are no longer included in wheels.
- Tests no longer silently skip when untracked legacy fixture data is absent.
- `mark` output format selection no longer depends on process-wide `sys.argv` state.
- The image-provider User-Agent now follows the package version.

[Unreleased]: https://github.com/kfuku52/nwkit/compare/v0.30.0...HEAD
[0.30.0]: https://github.com/kfuku52/nwkit/compare/v0.29.0...v0.30.0
[0.29.0]: https://github.com/kfuku52/nwkit/compare/v0.28.0...v0.29.0
[0.28.0]: https://github.com/kfuku52/nwkit/compare/v0.27.0...v0.28.0
[0.27.0]: https://github.com/kfuku52/nwkit/compare/v0.21.1...v0.27.0
