# Ten independent performance candidate evaluations

Baseline: `d7116f80b9979739c69a5431a5e34b822795408b` (NWKIT 0.43.54).
These are experimental implementations, not changes to the shipped library.

Each of ten subagents owns one numbered candidate directory and evaluates one
specific optimization. Read `DEVELOPMENT.md` and
`.agents/skills/benchmark-performance/SKILL.md` before measuring. Investigate
primary papers and official documentation, link the exact relevant sources,
and inspect existing optimizations before claiming novelty.

Use `/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python`. The original `.venv`
cannot execute on this host and is preserved. The fresh environment's initial
ETE wheel failed to import (`PyInit_tree`); a same-version source rebuild is
being performed. Do not change dependency versions to hide this failure.

All CPU-heavy benchmarks and tests must use an exclusive advisory lock on
`/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock`, held with `fcntl.flock`.
Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
MKL_NUM_THREADS=1` before starting Python. Import/fixture setup goes inside
the lock if substantial. Do not spawn multiple timed workers. Other agents
may research or write their experimental scripts while the lock is occupied.

Compare the same input, settings, seed, API behavior, and output fields.
Use warmups and at least seven alternating baseline/candidate measurements
where feasible. Save all timings, median, spread, and equivalence residuals
to JSON. Separate setup cost from amortized cost and state whether a benchmark
is a kernel or whole operation. Measure peak Python allocation separately
with tracemalloc when useful; it excludes some native allocations. Do not
describe microbenchmark gains as CLI or whole-program gains. Inspect numerical
edge cases, exception behavior, random streams, and memory tradeoffs.

Write `benchmark.py`, `results.json`, and `REPORT.md` in your own directory.
Keep changes out of `nwkit/`, `tests/`, existing studies, and other candidates.
Prototype with a local function or a scoped monkeypatch; preserve scientific
search bounds, models, tolerances, and seeds. Run appropriate existing tests
via `tools/check.py test -- ...` under the same lock, with scoped prototype
patches if practical, and explicitly label baseline-only test runs.

Report one of: adoptable within tested scope; conditional (needs integration
or broader evidence); reject; already implemented. Include exact commands,
environment, workload sizes, references, benefits, risks, skipped checks, and
the boundary between demonstrated effects and unverified extrapolation.

## Observed host contention

Other user-owned numerical/build jobs were active on this host during the
experiment. The advisory lock coordinates our ten agents only. Do not stop or
modify other jobs. These are exploratory local measurements, not uncontended
release performance guarantees. Retain raw timings and spreads, amortize short
operations into at least 20–50 ms batches where feasible, and treat modest or
inconsistent gains as inconclusive. Large effects still require numerical
equivalence and a stated tested scope; a quiet-host confirmation is needed
before publishing a stable performance claim.

The temporary environment passed the complete documented preflight after the
same-version ETE source rebuild: Python 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
ETE 4.4.0, pandas 3.0.6, NWKIT 0.43.54, macOS 27.0.1 ARM64. No changes were
made to the original environment.
