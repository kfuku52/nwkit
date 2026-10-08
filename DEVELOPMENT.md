# Developing NWKIT

Run commands below from the repository root. Use Python 3.10 or newer and an
isolated environment with the same development-tool constraints as CI:

```sh
python --version
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -c constraints-dev.txt -e '.[dev,image]'
```

The image extra needs a native Cairo installation; see the installation notes
in the [README](README.md). Runtime dependency ranges are intentionally separate
from `constraints-dev.txt`; do not add upper bounds without a demonstrated
incompatibility.

Before reusing an environment, run its Python and check imports, not just
`pip check` (which checks metadata):

```sh
python -c 'import sys; assert sys.version_info >= (3, 10); print(sys.version)'
python -m pip check
python -c 'from ete4 import Tree; import numpy, scipy.linalg, scipy.sparse.linalg, pandas, matplotlib, PIL; assert len(list(Tree("(A:1,B:1);", parser=1).leaves())) == 2'
python -m nwkit --version
```

Success means all commands exit zero and the last prints the checkout version.
The SciPy submodule imports exercise compiled libraries that a top-level
`import scipy` can leave unloaded. Passing preflight is still not a substitute
for the affected command tests below.
Activate the environment first; a system `python3` may be older than supported.
If `.venv` belongs to another host/architecture or cannot start, preserve it and
create a separate environment with a working supported interpreter. On POSIX:

```sh
NWKIT_DEV_DIR=$(mktemp -d)
python -m venv "$NWKIT_DEV_DIR/venv"
. "$NWKIT_DEV_DIR/venv/bin/activate"
python -m pip install -c constraints-dev.txt -e '.[dev,image]'
```

If imports fail despite consistent package metadata, retain the traceback and
identify the failing binary/package before repairing only the isolated environment.
For an ETE wheel import failure, a same-version source rebuild is an option when
a compiler is available; do not change dependency versions to hide the problem.

### SciPy 1.18.0 / 1.18.1 matrix-exponential nontermination

These two versions are excluded from runtime metadata because a finite matrix
encountered by the full OU optimizer causes `scipy.linalg.expm` to stop making
progress inside `matrix_exponential_d` / `scipy_cblas_idamax`. This reproduces
outside NWKIT on Linux ARM64, Python 3.12.14, NumPy 2.5.3:

```python
import numpy as np
from scipy.linalg import expm

matrix = np.array([
    [-6.08572431674702e16, -1.1080979213680645e17],
    [-1.1080979240863336e17, -2.0176415639168192e17],
])
expm(matrix)
```

Run this reproducer in a disposable process with an external timeout: both
versions exceeded an eight-second timeout for this 2-by-2 input, and the full
suite remained in the same native call for over 35 minutes. SciPy 1.17.1
returns in under one millisecond on the same input and platform (its NaN
result is rejected by the existing OU objective). The incompatibility is in
the dependency; NWKIT does not replace its exponential or narrow the scientific
search. Remove an exclusion when that version's corrected upstream distribution
completes this reproducer and passes the full OU optimizer and release checks.
Future versions remain eligible; this is not a blanket upper bound.

### SciPy 1.15.3 wheel on macOS 27 ARM64 / Python 3.10

On macOS 27.0 (26A428), the PyPI `macosx_14_0_arm64` wheel for SciPy 1.15.3
fails to load `_spropack` with `__thread_bss` / `offset field is not zero`.
This also reproduces with `import scipy.sparse.linalg` outside NWKIT. A fresh
download matches the published wheel, so clearing pip's cache does not fix it.

For this specific failure, use the same SciPy version's compatible
`macosx_12_0_arm64` wheel in an isolated Python 3.10 environment. The platform
tag selects a different binary build; it does not change your macOS version.
After installing NWKIT in that environment, run:

```sh
python -c 'import platform, sys, scipy; assert sys.version_info[:2] == (3, 10) and platform.system() == "Darwin" and platform.machine() == "arm64" and scipy.__version__ == "1.15.3"'
NWKIT_SCIPY_WHEELS=$(mktemp -d)
python -m pip download --only-binary=:all: --no-deps \
  --platform macosx_12_0_arm64 --python-version 3.10 --implementation cp --abi cp310 \
  'scipy==1.15.3' --dest "$NWKIT_SCIPY_WHEELS"
python -m pip install --no-deps --force-reinstall \
  "$NWKIT_SCIPY_WHEELS/scipy-1.15.3-cp310-cp310-macosx_12_0_arm64.whl"
python -m pip check
python -c 'import scipy.linalg, scipy.sparse.linalg'
```

Then run [the quick start](docs/guides/QUICK_START.md) and the affected tests.
This recovery was verified with Python 3.10.21, including a PROPACK SVD and
the tree/ASR example. An unconstrained reinstall may select the failing wheel
again. This is a binary-build workaround, not a global SciPy version pin or a
claim that all macOS/Python 3.10 builds are broken. Retire it when the normally
selected upstream wheel passes the compiled imports and command tests on the
affected platform; do not patch installed binaries or suppress import errors.

### CPython 3.12.14 timed traceback hang on macOS ARM64

With the bundled CPython 3.12.14 runtime on Darwin 27 ARM64, an opt-in timed
traceback stalled inside CPython's `dump_frame` during the long heterogeneous
WGD bootstrap test under coverage. All four scientific replicates completed
and pytest reported the test as passed, but teardown then waited indefinitely
in `faulthandler.cancel_dump_traceback_later()`. Native process samples showed
the diagnostic thread spinning while the main thread waited for it.

For this observed runtime failure, use the standard local diagnostics setting:

```sh
NWKIT_CI_DIAGNOSTICS=0 python tools/check.py release
```

The release command still runs all tests, including slow cases, branch coverage,
security checks and distribution checks with their original scientific settings.
Re-enable optional timed traceback dumps when the Python runtime completes the
same long coverage test with `NWKIT_CI_DIAGNOSTICS=1`. Preserve the failed log
and native samples when diagnosing this condition; it is separate from a failed
scientific assertion.

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

Coverage collects both line and branch opportunities; the configured 80% gate
applies to coverage.py's combined percentage, not to branch-only coverage.

Pass pytest paths and options after `--` for `quick` and `test`. `full`, `dist`,
and `release` reject test-selection arguments to avoid accidentally reporting a
partial run as complete. `dist` clears the derived `build/`, `dist/`, and
`direct-dist/` directories; use a separate source copy if those contain
artifacts you need to keep.

The `slow` marker identifies expensive numerical/bootstrap or concurrency
checks, not unreliable tests. Every full source CI run still executes them.
The complete test and quality jobs have a two-hour execution budget; the
scientific/bootstrap suite exceeds the former 30-minute limit. This budget
does not exclude slow tests or reduce their scientific replicate counts.
Keep small invariance tests and `tests/test_cli_contracts.py` in the quick
suite. The latter invokes the real parser and handler for every subcommand;
only external service boundaries are replaced with offline fixtures.

CI sets `NWKIT_CI_DIAGNOSTICS=1` for the existing test and quality jobs. The
same check entrypoint then prints bounded, escaped, line-terminated test-start
markers and individual test progress, and requests a
diagnostic stack after 600 seconds in a test; it does not terminate or retry
that test. Job deadlines, complete suite selection and scientific settings
remain unchanged. `tools/ci_diagnostics.py` prints an allowlist of runtime,
package and BLAS metadata without a hostname or environment dump. An archived
replay mismatch reports its exact cell, replicate, lane and bootstrap seed.

### Select checks by change

Use the rows as starting points, then follow callers of changed shared helpers.
Pass selected paths to `python tools/check.py quick -- ...`; use `test` instead
when static checks already passed on the same code. Add a regression case for a
new behavior or verified bug rather than duplicating existing assertions.

| Changed surface | Minimum focused coverage |
| --- | --- |
| CLI options, command examples, shared TSV policy | `tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py`, plus the affected command's tests |
| Tree reading/writing or rootedness | `tests/test_util_tree_io.py tests/test_tree_outputs.py tests/test_rooting_state.py`, plus affected commands; staged writes also need `tests/test_output_transaction.py` |
| Numerical/model code | The affected model tests and its CLI consumer; shared Gaussian/regression changes also need `tests/test_numerical_invariance.py`. Preserve independent reference comparisons and include relevant `slow` cases explicitly |
| Drawing/media | Affected `tests/test_draw*.py` or `tests/test_image*.py`; inspect representative rendered outputs as described below |
| Check/build tooling | `tests/test_check_tools.py tests/test_distribution_reproducibility.py`, then the changed runner mode; package changes also need `dist` |
| Prose/agent instructions only | Verify linked paths and try changed commands; `dist` matches documentation CI. A plain patch-version bump does not require numerical tests by itself |

A small offline starting check exercises real command handlers, output schemas,
and CLI compatibility using pytest temporary directories:

```sh
python tools/check.py quick -- tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py
```

Expect a nonempty pytest selection, no failures, and successful lint/format/type
checks. Exit status 5 (no tests collected) is not success. Add `-rs` after `--`
to explain skipped tests. `quick` without targets still runs thousands of tests;
it is not a tiny smoke test. Explicit `-m` replaces its `not slow` selection.

### Cost and environment boundaries

- The focused CLI check above uses small local fixtures; image services and the
  legacy SHIFT backend are replaced at external boundaries. It does not validate
  those live services or the R backend itself.
- `slow` is a cost marker, not an offline/network marker. Optional IQ-TREE,
  IQ-TREE library worker, PAML/MCMCTree and R backend tests may skip when their
  runtimes are absent. Inspect the selected tests and report those skips; do not
  install/download external runtimes merely to make the summary green.
- `full`/`release` include `pip-audit`, which needs package/advisory service
  access. `dist`/`release` use isolated builds that may fetch build dependencies.
  A failed audit or unavailable network is an incomplete gate, not a pass.
- Calibration and benchmark programs under `tools/` and study scripts under
  `examples/` are separate from routine checks. Read the relevant guide and
  runtime/data requirements before running them; use small inputs and a fresh
  temporary output directory for exploratory work.

Before a requested push/release, use the complete delivery checks in
[RELEASING.md](RELEASING.md). Focused successes do not replace that gate.

## CI coverage

For core tree/Gaussian performance regressions, run the identical harness on
each source checkout (warmup plus three timed repetitions):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python tools/benchmark_core.py --checkout /path/to/checkout --output /tmp/core.json
```

It covers deep consensus construction, RF, star-tree sister IDs, interval
iteration, and a Gaussian profile fit. Compare numerical fields as well as
timings. Memory is peak Python allocation, not total process RSS. Run without
competing CPU workloads when making release performance claims.

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
functions. Both existing and new functions have a hard ceiling of 40. Increases
relative to `tools/complexity_baseline.json` produce review warnings rather than
failing the check. Average complexity is informational, so deleting small
functions cannot make a cleanup fail. Complexity is a branching heuristic, not
a correctness or readability score; split functions only along useful responsibilities.

After reviewing changes and running the relevant tests, update the comparison
baseline when appropriate:

```sh
python tools/check_maintainability.py --update-baseline
```

The updater checks hard limits first, then records current values (including
reviewed increases) and removes deleted entries. It never changes a hard limit.
Exceptions belong in `tools/complexity_exceptions.json` and require an explicit
limit above 40, a rationale, and existing test-file paths documenting the
relevant coverage. Those tests must be run; listing them is not proof of success.
The legacy `draw_main` limit of 50 is preserved as a documented exception.
Review exception changes explicitly; never raise a limit merely to pass a check.
For an unchanged function moved to another module, move its baseline and any
exception keys with it and review that move. Remove exceptions for deleted functions.

Ruff formatting remains mandatory. Apply `python -m ruff format` to edited files
before committing; CI continues to use `ruff format --check`.

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
