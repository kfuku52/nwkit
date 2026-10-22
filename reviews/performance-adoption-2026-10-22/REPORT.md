NWKIT 0.43.55 performance adoption — 2026-10-22

Eight candidates were integrated in separate implementation commits after
additional numerical, memory and full-operation checks. Candidate 06 remains
deferred; candidate 07 was not adopted. The earlier
[ten-candidate study](../performance-candidates-2026-10-22/REPORT.md) records the
initial prototypes; this document describes the production implementation.

| Candidate | Implementation commit | Additional result | Important implementation decision |
| --- | --- | --- | --- |
| 01 K permutations | `033224c8`, follow-up `80d315fa` | 4.70x, complete 256-tip/199-draw test | Local factor reuse; preserve error-aware refits and RNG consumption on failure |
| 02 dense GLMM | `1232e1ef` | 1.32x, complete 64-tip Poisson fit | Replace all seven redundant general solves of existing Cholesky factors, including precision construction |
| 03 cluster scores | `87a36d99` | 3.74x, complete 6,000-row/1,500-event regression helper | Stable grouped indices retain the original dot/outer-product order rather than using bincount |
| 04 grouped uncertainty | `f86a5a82` | 17.14x, complete 128-observation projection | Tile both observation axes against a 4 MiB product-workspace budget |
| 05 diagonal predicate | `07093112` | 13.58x, complete 512-diagonal covariance factor | Retain dense first-row exit and object-array off-diagonal semantics |
| 08 count transitions | `85f84907` | 1.56x, complete rare-family fixed-parameter likelihood | Read-only shifted gamma view needs O(M) extra storage, replacing the prototype's O(M²) table |
| 09 DTT simulation | `e8002009` | 3.12x, complete 256-tip/3-trait/128-draw simulation | At most 32 replicas with a 4 MiB workspace estimate; use single-replica processing when needed |
| 10 missing-response PGLS | `d09896e` | 1.44x, complete 128-tip/4-trait high-missing fit | Direct observed construction only for dimension >=512 and at least two-thirds missing; publish full covariance once at the final fit |

The ratios are exploratory host-local observations, not an overall NWKIT speedup.
Other user-owned CPU jobs were active and were preserved. These are complete
operations at the stated API/helper level, excluding imports and shared fixture
construction. Each of seven alternating baseline/current measurements uses a
fresh process, one numerical thread, a warmup and calibrated short-call batches.
The baseline is a clean archive of `d7116f80b9979739c69a5431a5e34b822795408b`.
Remote master was refreshed to `3da1df5724b432e59b3b1f04d60c649852cb55b3`
(policy-only update) before implementation; that is the fixed pre-push base.

All eight comparisons retain raw timings, equivalence checks and peak traced
allocation in the matching JSON files. Seven operations matched every compared
output exactly. GLMM's maximum scaled numeric difference was `1.24e-6`, with
convergence/status retained; the existing independent-reference regression
tolerances were unchanged. The core RF, consensus and Gaussian control harness
also retained exactly equal numerical fields.

The additional cancellation experiment found an important limitation of the
original candidate 03: a shuffled one-column score with terms near `1e16`
gave dot=4 and bincount=0. Stable grouped row indices eliminate the O(N*G)
repeated scans while retaining the same per-event BLAS products and event
accumulation order. The final large-case traced peak is effectively unchanged.
Unbalanced-event and cancelling-score regressions cover this decision.

Candidate 04's traced peak at the measured size was 2.43→3.15 MB. Candidate 09
was 0.134→1.844 MB; the workspace budget limits batch growth and retains the
single-replica method for larger dimensions. Candidate 08's revised linear
gamma workspace reduced the prototype's memory penalty to about 0.6% in the
measured complete likelihood (1.183→1.189 MB). These are tracemalloc readings,
not complete RSS, and exclude retained inputs/output storage. No blanket memory
reduction is claimed.

Re-evaluating conditional candidate 06 after the first seven integrations gave
dense complete-profile ratios 1.11–1.17x and a diagonal ratio 0.89x. Its small
benefit, added arrays and optimizer-order sensitivity did not justify another
path. Candidate 10 was stronger when restricted: 128/256-tip high-missing fits
gave 1.35/1.45x, with every pair faster and nonoverlapping raw ranges, whereas
moderate cases were 0.94/0.98x. This established the conservative dispatch used
in production. The signatures, all output fields and full final covariance
were exactly equal. See [conditional-results.json](conditional-results.json).
That experiment ran against the seven-candidate target sources, reproducible
from `80d315fa`; its guarded prototype script intentionally requires those
pre-adoption functions.

The individual affected checks passed: signal 47; dense GLMM 58; cluster 35;
projection 59; diagonal 58; count-model/CLI/reference 110; DTT 63; multivariate
66. Later integration checks passed 235, including all CLI contracts, image
artifacts, K failure-state handling, DTT exact seed/thread behavior, 1–64 traits,
extreme scaling, workspace boundaries and the cancellation regression. Counts
overlap. Each `*-check.log` records the exact `tools/check.py quick -- ...`
command and selections. Successful quick checks include Ruff lint/format and
mypy, not only pytest. Existing scientific test tolerances, seeds, search
bounds, work guards and inference methods are unchanged.

Reproduction uses Python >=3.10 with the project's development constraints.
The actual environment is CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas
3.0.6, ETE 4.4.0 source build, and CairoSVG 2.9.1 on macOS 27.0.1 ARM64.
The original incompatible `.venv` remains untouched. Existing build/dist
artifacts were backed up to
`/tmp/nwkit-performance-adoption-20261022/artifact-backup` before delivery gates.
An x86 Cairo on `/usr/local/lib` could not load; using the existing ARM64
library through `DYLD_FALLBACK_LIBRARY_PATH=/Users/kf/miniforge3-arm64/lib`
passed native PNG rendering and compiled scientific imports. No optional
R/IQ-TREE/PAML/GeneRax executable was installed to turn skips into passes.

```sh
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
$PY reviews/performance-adoption-2026-10-22/benchmark.py \
  --baseline /tmp/nwkit-performance-adoption-20261022/baseline \
  --case signal --output /tmp/signal-comparison.json
# Other cases: cluster, glmm, projection, diagonal, wgd, dtt, multivariate.

DYLD_FALLBACK_LIBRARY_PATH=/Users/kf/miniforge3-arm64/lib \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
NWKIT_TEST_WORKERS=3 NWKIT_CI_DIAGNOSTICS=0 $PY tools/check.py release
```

The release gate covers uncached mypy, lint/format, dependencies, Bandit,
pip-audit, the complete regression suite including slow tests, branch coverage,
complexity ceilings, independent wheel/sdist builds and byte reproducibility.
Replicated scientific studies remain separate: the statistical models and
uncertainty methods were not changed. A subsequent validation record records
the completed release outcome and actual optional skips.
