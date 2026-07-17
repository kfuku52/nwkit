# Changelog

All notable changes made after the `v0.21.1` tagged release are tracked here.

## [Unreleased]

### Fixed

- Root transfer and the root alignment used by `compose` and `transfer` now
  preserve internal-node support, names, and NHX properties on their canonical
  unrooted branch splits instead of leaving them attached to nodes whose
  descendant taxa change during rerooting.
- `compose` annotation and value sources no longer change the target rooting
  when `--root-source` is absent.
- Failed root transfers no longer leave the target tree partially rerooted.
- Split-based transfer now treats equal support, names, and NHX values on the
  two children of a bifurcating root as one root-edge value. Conflicting values
  remain ambiguous and are not transferred.

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

[Unreleased]: https://github.com/kfuku52/nwkit/compare/v0.29.0...HEAD
[0.29.0]: https://github.com/kfuku52/nwkit/compare/v0.28.0...v0.29.0
[0.28.0]: https://github.com/kfuku52/nwkit/compare/v0.27.0...v0.28.0
[0.27.0]: https://github.com/kfuku52/nwkit/compare/v0.21.1...v0.27.0
