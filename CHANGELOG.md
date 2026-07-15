# Changelog

All notable changes made after the `v0.21.1` tagged release are tracked here.

## [Unreleased]

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

[Unreleased]: https://github.com/kfuku52/nwkit/compare/v0.27.0...HEAD
[0.27.0]: https://github.com/kfuku52/nwkit/compare/v0.21.1...v0.27.0
