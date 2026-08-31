# Developing NWKIT

Use an isolated environment and the same development-tool constraints as CI:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -c constraints-dev.txt -e '.[dev,image]'
```

The image extra needs a native Cairo installation; see the installation notes
in the [README](README.md). Runtime dependency ranges are intentionally separate
from `constraints-dev.txt`; do not add upper bounds without a demonstrated
incompatibility.

## Choose the smallest useful check

| Command | Checks |
| --- | --- |
| `python tools/check.py quick` | Ruff lint/format, incremental mypy, all tests except `slow` |
| `python tools/check.py quick -- tests/test_asr.py -k time_units` | The same static checks, with a focused pytest selection |
| `python tools/check.py test -- tests/test_numerical_invariance.py` | Only the requested tests, with no default marker exclusion |
| `python tools/check.py quick -- -m slow tests/test_regress.py` | Explicitly selected slow tests, plus static checks |
| `python tools/check.py full` | Uncached mypy, lint/format, dependency/security checks, **all** tests with branch coverage, complexity checks |
| `python tools/check.py dist` | Independent wheel/sdist builds, archive contents, metadata and byte reproducibility |
| `python tools/check.py release` | `full` followed by `dist` |

Pass pytest paths and options after `--` for `quick` and `test`. `full`, `dist`,
and `release` reject test-selection arguments to avoid accidentally reporting a
partial run as complete. `dist` clears the derived `build/`, `dist/`, and
`direct-dist/` directories; use a separate source copy if those contain
artifacts you need to keep.

The `slow` marker identifies expensive numerical/bootstrap or concurrency
checks, not unreliable tests. Every full source CI run still executes them.
Keep small invariance tests and `tests/test_cli_contracts.py` in the quick
suite. The latter invokes the real parser and handler for every subcommand;
only external service boundaries are replaced with offline fixtures.

## CI coverage

`tools/ci_matrix.py` classifies both sides of renamed paths. All source changes
run full quality/security/coverage checks on Linux with the newest supported
Python, plus the complete test suite on the minimum supported Python.

- Documentation-only changes build and inspect reproducible distributions
  without installing the numerical runtime dependencies. A plain patch-version
  change in `nwkit/__init__.py` does not turn this into a source run.
- Numerical changes retain the minimum/newest Python pair. Filesystem, drawing,
  CLI, and otherwise unclassified source changes also run macOS/Windows tests
  and a clean macOS image installation.
- Dependency/build/workflow changes, weekly runs, manual runs, and major/minor
  release version changes exercise all supported Python versions and platforms.
- Pushes are limited to the default integration branches, avoiding a duplicate
  branch-push run for ordinary pull requests. Superseded runs are cancelled.

Windows currently needs a source-path fix when building ETE4 4.4.0.
`tools/build_ete_windows.py` verifies the upstream source checksum and checks
that the exact patch site still exists. Its wheel is cached by OS, architecture,
Python and the build script hash (including the source checksum and patch).
Remove the workaround and cache when upstream fixes extension-module paths or
provides a compatible Windows wheel; revalidate Windows before removing it.
This CI workaround is not a runtime ETE upper bound.

## Keep complexity from growing

`tools/check_maintainability.py` measures individual functions, methods and nested
functions. Existing functions may not exceed their recorded complexity; new
functions have a ceiling of 40. Average complexity is informational, so deleting
small functions cannot make a cleanup fail.

After a verified cleanup, run:

```sh
python tools/check_maintainability.py --update-baseline
```

The updater validates first, then records reductions/new functions and removes
deleted entries. It cannot approve a ceiling increase. For an unchanged function
moved to another module, move its baseline key with it and review that move.

Drawing stages exchange the typed records in `draw_types.py`. Input readers and
primitives live in `draw_helpers.py`, validation/measurement/layout in
`draw_setup.py`, ordered rendering and quality evaluation in `draw_render.py`,
and image/report serialization in `draw_output.py`. `draw.py` keeps the command
entrypoint and orchestration. Preserve artist order, coordinates, and existing
option meanings when changing these boundaries.

## Measure actual behavior before and after

The benchmark runner imports the selected checkout in a fresh process for each
repetition and uses one BLAS/OpenMP thread. `--baseline` can point at an archived
source tree; it does not need an installed baseline package or a Git checkout.

```sh
python tools/benchmark.py --case version --baseline /path/to/old/source
python tools/benchmark.py --case help --baseline /path/to/old/source
python tools/benchmark.py --case regress-help --baseline /path/to/old/source
python tools/benchmark.py --case gaussian --tips 512 --repeat 3 --baseline /path/to/old/source
```

JSON output includes wall time, process peak RSS, and comparable outputs. CLI
time includes command-module imports and parsing; Gaussian time covers fitting
after imports and fixture construction. RSS includes the entire child process;
it is reported as `null` on Windows, where `resource` is unavailable. The Gaussian
case uses a seeded balanced tree, one response/predictor, REML, and automatic
lambda estimation. It is not a benchmark of every family or tree shape.

CLI output hashes ignore only the package version. Gaussian coefficients,
coefficient covariance, likelihood, evolutionary parameter, variance components
and convergence status are compared (`rtol=1e-4`, `atol=1e-6` for numbers); a
mismatch fails the command. Inspect actual numerical differences as well as
the pass/fail result. Do not run other numerical jobs during timing, and report
the environment, repetitions and output-equivalence check with any speed claim.

Scalar likelihood searches retain all objective values but at most two complete
fits, releasing cached arrays after selecting the winner. Ordinary automatic
lambda fits also reuse one validated Brownian covariance: tip variances stay
fixed while shared covariances scale with lambda. Preserve resource-limit
exceptions instead of turning them into invalid likelihood candidates. Regression
tests explicitly cover that contract and the lifetime of evicted arrays.

For drawing changes, compare deterministic SVGs/reports and inspect representative
rendered layouts, including time annotations and tip images, in addition to running
the drawing tests. For release checks, see [RELEASING.md](RELEASING.md).
