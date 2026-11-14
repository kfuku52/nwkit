# Second-round production adoption protocol

Baseline: `8664baa7e88afae7e384c0e4d8d16437846493a2` (0.43.56).
Candidate research: [ten independent evaluations](../performance-candidates-round2-2026-10-22/REPORT.md).
This follow-up measures the production implementations together, using fresh
baseline/current interpreter processes without loading prototype patches.

Use the supported isolated environment and import preflight in
[DEVELOPMENT.md](../../DEVELOPMENT.md). The recorded environment is CPython
3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6, ETE 4.4.0, macOS ARM64.
An existing ARM64 Cairo library is selected through
`DYLD_FALLBACK_LIBRARY_PATH=/Users/kf/miniforge3-arm64/lib`. Numerical threads
are fixed to one by the parent runner before each fresh process starts.

```sh
BASELINE_DIR=$(mktemp -d)
git archive 8664baa7e88afae7e384c0e4d8d16437846493a2 | tar -x -C "$BASELINE_DIR"
python reviews/performance-adoption-round2-2026-10-22/benchmark.py \
  --case lambda --baseline "$BASELINE_DIR" \
  --output /tmp/nwkit-lambda-adoption.json
```

Repeat for the `--case` choices in [benchmark.py](benchmark.py). The fixtures
in [fixtures.py](fixtures.py) are copied from the reviewed research scripts;
imports bind to the selected checkout. Each row stores the actual imported
NWKIT path and version. Some measurements precede the metadata-only patch
version update; numerical implementations are unchanged by that update.

Each operation is warmed, calibrated to approximately 50 ms when short, and
measured in seven alternating baseline/current pairs. Process startup and
fixture creation are outside timings. A separate invocation under tracemalloc
records allocation peaks; those are not total process RSS. Every timed result
is compared field by field using saved NumPy signatures. Only lambda uses the
existing scientific allowances (estimate absolute 1e-5, other fields 1e-7);
ER posterior values use absolute 1e-13. Other outputs, including free ER rate
estimation, final PCA fields, diagnostic status/rank and stochastic event times,
must match exactly. Tolerances in existing tests are unchanged.

The operation-local measurement lock serializes this experiment. Other
user-owned CPU work remains untouched. Raw timing ranges and every comparison
are retained, including controls without an expected speed benefit. Do not
infer effects for other tree sizes, error models, populations or BLAS builds.

Focused tests use `python tools/check.py quick -- ... -rs`; Full-OU also
explicitly includes its slow optimizer test with `-m ''`. Delivery uses
`NWKIT_TEST_WORKERS=3 NWKIT_CI_DIAGNOSTICS=0 python tools/check.py release`,
with all numerical thread variables set to one and the same Cairo path.
The disabled optional traceback diagnostic follows the documented CPython
runtime workaround; all slow regression tests and the coverage floor remain.
Existing build artifacts are copied aside before the release runner clears
build directories. External-backend skips and opt-in scientific studies are
reported separately from regression passes.
