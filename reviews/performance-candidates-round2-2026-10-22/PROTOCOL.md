# Second round: ten independent performance candidates

Baseline: `86c8f2c537769d852c17e2c789340f6c9de4ef54`, NWKIT 0.43.55.
This round evaluates prototypes only. Do not edit production code, tests,
existing studies, versions, or Git history. Each of ten distinct subagents
owns one numbered directory below this directory. At most three run at once.

Another task began editing `nwkit/time_tree.py` after our clean start. Preserve
that unrelated change. An immutable `git archive` of the baseline is available
at `/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline`. Run benchmarks and tests
with that directory as cwd and PYTHONPATH (plus your own plugin/helper directory
when necessary), and record `nwkit.__file__` to confirm the imported snapshot.
Use an absolute path to your experimental script. Run the snapshot's
`tools/check.py` with its own tests, rather than mixing live workspace sources.
No Git metadata is present in the archive; avoid unrelated Git/build tests.
On macOS `/tmp` resolves to `/private/tmp`; compare resolved paths rather than
rejecting that canonical spelling. The snapshot import preflight passed too.

To recreate the snapshot if the temporary directory has been removed, make an
empty directory at that path and extract this exact revision, not current HEAD:

```sh
mkdir -p /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
git archive 86c8f2c537769d852c17e2c789340f6c9de4ef54 | \
  tar -x -C /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
```

Use DEVELOPMENT.md to recreate a supported development environment if the
temporary interpreter is no longer available; retain the recorded versions
when attempting a numerical/timing reproduction.

Read AGENTS.md, DEVELOPMENT.md, the relevant scientific guide, and
`.agents/skills/benchmark-performance/SKILL.md`. Follow the previous study's
evidence boundaries without repeating its ten hypotheses. Investigate primary
papers or official API documentation and cite exact URLs in your report.
References motivate a hypothesis; measurements establish its local benefit.

Use `/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python`. Its complete documented
preflight passed again: CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas
3.0.6, ETE 4.4.0, macOS ARM64. Preserve the original incompatible `.venv`.

Before importing numerical modules set OPENBLAS_NUM_THREADS=1,
OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1, MKL_NUM_THREADS=1 and
NWKIT_CI_DIAGNOSTICS=0. For image tests the existing ARM64 Cairo is available
through DYLD_FALLBACK_LIBRARY_PATH=/Users/kf/miniforge3-arm64/lib.
Every substantial fixture, test or measurement run must hold an exclusive
`fcntl.flock` on `/tmp/nwkit-perf-round2-20261022-rkzoGe/measurement.lock`.
Use the shared `support.py` helper if convenient. Do not launch parallel timed
workers. Other user-owned CPU jobs are active; do not stop or alter them.

Measure representative small/control and larger workloads, warm up both
implementations, retain at least seven alternating baseline/candidate timings,
and calibrate short kernels into approximately 20-50 ms samples where feasible.
Keep setup-inclusive and amortized measurements distinct. Include at least one
complete consuming operation, or explicitly limit the conclusion to a kernel.
Use identical inputs, seeds, optimizer bounds and tolerances. Compare every
relevant numerical field and status, including edge/failure behavior. For
stochastic code verify random stream and threaded/serial behavior as relevant.
Do not weaken scientific checks to obtain agreement or faster execution.

Measure peak traced allocations separately when relevant, explaining that
tracemalloc is not total RSS. Bound any additional working arrays/cache lifetime.
Reject claims when timings are too noisy or correctness is insufficient. Report
raw timings, medians, ranges, residuals, environment, exact commands, sources,
tradeoffs, and skipped checks. Use one of: adoptable within tested scope;
conditional; reject; already implemented. Distinguish integration readiness
from exploratory evidence.

In your numbered directory save `benchmark.py`, `results.json`, `REPORT.md`,
and relevant small logs. Prototype with local functions or scoped monkeypatches.
Run affected existing tests through `tools/check.py test -- ...` under the same
lock, with a prototype pytest plugin when practical. Explicitly label tests
that exercise only baseline code. Avoid duplicating broad expensive suites;
independent-reference and focused failure/invariance checks are preferred.
Use fresh temporary outputs under `/tmp/nwkit-perf-round2-20261022-rkzoGe/`.
Do not create additional agents. Return your evidence and limitations to root.
