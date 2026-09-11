# OU shift response validation — 2026-09-10

This records the initial validation phase. Its results and original source hashes
are historical. The subsequent [bounded improvement attempt and integration](SHIFT_CALIBRATION_DECISION.md)
closed without a replacement calibration and supersedes the original no-commit scope.

The pBIC backend acceptance defect is closed by a numerical capability check.
The native grid-envelope method has independent numerical and simulation audits,
but **the complete scientific acceptance criteria have not all been met**.
In particular, the paired experiment does not demonstrate that power loss is
at most five percentage points. Continuous-alpha calibration, heterogeneous
measurement-error calibration and later-stage convergence decisions remain
outside the validated error-control scope. No new inference engine or automatic
criterion substitution was introduced by this response.

## Implementation and ownership

The existing shift implementation was imported read-only from the separate
working checkout at `/Users/kf/repos/nwkit`; its hashes and import history are in
`reviews/shift-implementation-base.json`. The native inference and candidate
sources used here have SHA-256 `5a26bf6eea4bb1db490479e4e9fd52282575ef533e687e568b6fd35b10dbc3ca`
and `bdf10b084ea6d9b8cda0ffe0a7758b617a41e694ad481aed01ff7c65f4e510a0`.
They were not changed during this response. The source checkout and external
kfl1ou checkout were not edited. The new changes concern backend acceptance,
output interpretation, independent verification and integration in this worktree.

`nwkit.shift_backend_probe` is shared by production IC inference and the alpha
and joint research drivers. It checks 22 quantities through public kfl1ou fits
against an independent dense Gaussian/QR calculation: fixed and estimated alpha,
fixed and stationary roots, unconstrained and singleton coordinates, and a
shared-optimum constrained fit. The tolerance is 2e-6. The installed original
3.0.9 fails; the separately installed corrected 3.0.9 passes. This identifies
behavior without inventing a release version or pinning all future packages to
one hash. The probe preserves R random-number state, runs in the fitting process,
and records the resolved package library plus installed-file SHA-256 hashes.
It is a bounded capability test, not a proof of all optimizer or pBIC properties.
Correcting coordinates does not repair the nonregular small-alpha behavior of
the pBIC approximation; see the [mathematical boundary diagnosis](SHIFT_PBIC.md).

A real CLI test confirms rejection before reading user data and preservation of
pre-existing outputs. Another accepts BIC with the original backend, and a real
corrected-backend pBIC fit exports the attestation. See the
[IC example](examples/shift/response-ic/model.json) and the
[calibrated example](examples/shift/response-calibrated/model.json).

Schema-7 output now explicitly distinguishes finite OU, the scaled-effect drift
limit, and the independent-tip limit; it records likelihood and simulation scale
rules, finite-grid scope, interpretation of shared-family selection and failure
policy. Missing finite optima remain missing. A regime map passed to ASR causes
a new fit and does not transfer these boundary models or selection uncertainty.

## Independent numerical checks

`tools/shift_continuous_reference.py` constructs means and covariance by branch
recursion, projects using an independently constructed orthonormal contrast basis,
and profiles a dense Gaussian likelihood without production fitting helpers.
Tests compare all candidates on 4-, 8- and 16-tip trees at alpha-height 0, 0.4,
30 and infinity, including a pectinate tree. Known-error tests cover zero and
positive process variance. Separate existing tests cover units, permutations,
constant offsets, candidate enumeration and numerical boundaries.

The continuous reference scouts log alpha from 1e-8 to 1e8 and refines each
interior local maximum, also evaluating exact 0 and infinity. Known-error
variance has its own continuous profile. It is not a certificate of a global
continuous optimum. A fixed numerical fixture improves the grid likelihood by
about 3.5e-5; this demonstrates that the grid is an approximation, not that its
model choices or calibrated probabilities match a continuous-search procedure.
The proposed 1% bound on selection disagreement has **not** been established.
Consequently no continuous inference mode is offered or claimed to be validated.

## Null experiment

The frozen primary protocol uses seed 20260923 and 1,000 independent datasets
for each combination of 4/8/16 tips, balanced/pectinate shape and true
alpha-height 0.2/2.1, without measurement error. Both values are between fitted
grid points; they test two off-grid truths, not the whole continuous parameter
space. The generator uses independent
branch-recursion covariance with process tip variance one. Both convergence-off
and convergence-on searches use each dataset; these are paired settings, not
additional independent data. The main B=199 experiment has 12,000 datasets and
24,000 fits. It uses the actual CLI default search/grid/level, but varies the
Monte Carlo seed by dataset; it is not a study of a permanently fixed seed=1.

Before seeing validation outcomes, the acceptance gate was fixed as a one-sided
95% exact binomial upper bound, Bonferroni adjusted over 24 cell/lane comparisons,
no greater than 0.075. With N=1,000 this permits at most 51 false selections
plus failures in each of the 24 comparisons. Every failed fit remains in the denominator and is counted
as a false selection for this worst-case bound. This is an engineering tolerance
of 2.5 percentage points above nominal 0.05, not a proof of uniform 5% control.
A separate B=999 check uses the same 2,000 observations in the two eight-tip
balanced cells; it is a sensitivity comparison, not 2,000 additional independent
datasets. Its four comparisons have their own simultaneous bound (at most 56 false
selections plus failures per N=1,000). These two families of bounds are reported
separately, not as a joint 95% guarantee across all 28 comparisons.

The protocols, source snapshots, observations, generating covariances, seeds,
every attempted fit and cell summaries are retained in
[primary evidence](examples/shift/null-contract-validation/protocol.json) and
[B=999 evidence](examples/shift/null-contract-b999/protocol.json).
The read-only verifier regenerates every input, independently checks the selected
likelihood and predicted means, checks stage probability metadata and denominators,
reconstructs all cellwise bounds, and optionally replays the complete search.

All 24 primary comparisons passed the prespecified gate, with zero failed fits.
The largest simultaneous one-sided 95% upper bound was **6.23%**, below 7.5%.
The two convergence lanes have the same first-stage decisions on these data;
they are still retained separately and the prespecified 24-comparison adjustment
was not reduced after observing this agreement.

| Tips | Tree shape | True alpha × H | B=199 false selections / 1,000, off / on | Simultaneous upper, off / on |
|---:|---|---:|---:|---:|
| 4 | balanced | 0.2 | 23 / 23 | 4.01% / 4.01% |
| 4 | balanced | 2.1 | 31 / 31 | 5.01% / 5.01% |
| 4 | pectinate | 0.2 | 24 / 24 | 4.14% / 4.14% |
| 4 | pectinate | 2.1 | 35 / 35 | 5.50% / 5.50% |
| 8 | balanced | 0.2 | 39 / 39 | 5.99% / 5.99% |
| 8 | balanced | 2.1 | 32 / 32 | 5.14% / 5.14% |
| 8 | pectinate | 0.2 | 35 / 35 | 5.50% / 5.50% |
| 8 | pectinate | 2.1 | 27 / 27 | 4.51% / 4.51% |
| 16 | balanced | 0.2 | 41 / 41 | 6.23% / 6.23% |
| 16 | balanced | 2.1 | 31 / 31 | 5.01% / 5.01% |
| 16 | pectinate | 0.2 | 41 / 41 | 6.23% / 6.23% |
| 16 | pectinate | 2.1 | 36 / 36 | 5.62% / 5.62% |

Both B=999 cells and convergence lanes also passed, with zero failures:

| True alpha × H | B=999 false selections / 1,000, off / on | Simultaneous upper, off / on |
|---:|---:|---:|
| 0.2 | 42 / 42 | 5.85% / 5.85% |
| 2.1 | 40 / 40 | 5.62% / 5.62% |

The [B=999 audit](examples/shift/response-audits/b999.json) regenerated every
input, independently checked 4,000 winning fits and replayed four complete
searches. Maximum independent likelihood error was 1.20e-14.
The [primary audit](examples/shift/response-audits/primary.json) regenerated all
12,000 inputs, independently checked all 24,000 winning fits and mean predictions,
and replayed 24 complete searches. Maximum independent likelihood error was
7.38e-12. All saved counts, probability metadata and simultaneous bounds matched.

These cells do not validate arbitrary tree shapes, shallow/long-branch extremes,
continuous alpha, exact endpoint generating distributions or measurement errors.
A separate [12-dataset known-error timing pilot](examples/shift/known-error-timing-pilot/summary.json)
completed 24 fits, but is not a
calibration validation. Known-error and later-stage procedures remain plug-in
methods. Selecting a shared-effect family does not prove biological convergence.
The overall scientific acceptance gate remains open despite any primary-null
cellwise success.

## Paired detection and recovery

The [paired replay](examples/shift/calibration-paired/summary.json) uses all 600
no-error datasets in the existing 700-dataset grid-envelope bundle. It reruns the
archived plug-in engine (`e30e9aab8c4886811c315817e0b730e7d3966bfddb743899ec18c6cded6809cf`)
on the same observations with the same B=199 and Monte Carlo seeds. This isolates
the software change on existing data; it is not a fresh held-out power study.
Convergence is enabled. Known-error cases are excluded because their first-stage
method did not change. All paired fits completed.

| Primary alternative / generating root | Plug-in detection | Envelope detection | Paired losses / gains | One-sided 95% upper bound on power loss |
|---|---:|---:|---:|---:|
| Single / fixed | 36/50 | 33/50 | 3 / 0 | 14.8% |
| Single / stationary | 23/50 | 23/50 | 0 / 0 | 5.8% |
| Distinct / fixed | 25/50 | 21/50 | 4 / 0 | 17.4% |
| Distinct / stationary | 30/50 | 26/50 | 4 / 0 | 17.4% |
| Shared / fixed | 19/50 | 19/50 | 0 / 0 | 5.8% |
| Shared / stationary | 29/50 | 28/50 | 1 / 0 | 9.1% |

The bound is conservatively the exact upper bound on the paired-loss probability;
gains are nonnegative. It is per-cell, without multiplicity adjustment. Even with
zero observed losses, N=50 cannot demonstrate a five-point noninferiority margin
with this bound. The distinct-shift cells show an observed eight-point reduction.
The planned power criterion is therefore **not demonstrated**, and point estimates
also show a real caution against treating conservatism as cost-free. The evidence
retains exact tip-partition recovery and mean RMSE separately from detection.
It does not establish recovery of each true shift location or correctness of
convergence claims. The [independent audit](examples/shift/response-audits/paired.json) checked all
1,200 selected fits with maximum likelihood error 2.85e-14.
No tuning or further sampling was performed after these
results to manufacture a pass.

The inherited 700-dataset envelope bundle was separately audited read-only:
all generating inputs and source snapshots match, and saved null selections are
7/250. This pooled number is descriptive and cannot replace the larger cellwise
experiment. Historical IC and plug-in bundles retain their historical meaning.

## Reproduction and software checks

Use one BLAS/OpenMP thread per process; these commands require a fresh output
path and do not overwrite the saved evidence:

```sh
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python tools/validate_shift_null_contract.py --phase validation \
  --replicates 1000 --workers 4 --output /tmp/shift-null-new
python tools/verify_shift_null_contract.py /tmp/shift-null-new --replay-stride 1000
python tools/validate_shift_null_contract.py --phase validation \
  --replicates 1000 --bootstrap-replicates 999 --cell-ids 4,5 \
  --workers 2 --output /tmp/shift-null-b999-new
python tools/compare_shift_null_calibration.py --workers 2 \
  --output /tmp/shift-paired-new
```

`verify_shift_calibration.py` now defaults to read-only auditing. Refit requests
must name a separate output directory outside the evidence bundle; source hash
mismatches are rejected by default. An explicit `--allow-source-revision`
requires full bootstrap replay and records changed hashes. Archived null-contract
bundles use `verify_shift_null_contract.py --frozen-engine` to audit their
hash-verified historical engine, rather than assert current-CLI validation. The input bundle is not rewritten to fit a new verifier.

The initial response-worktree repository check passed locally on macOS/Python 3.10:
3,326 tests passed, 27 skipped, branch coverage 85% (required minimum 80%).
Ruff, nonincremental mypy (173 modules), dependency consistency, Bandit,
`pip-audit` and complexity checks passed. The original and corrected R libraries
were both enabled for the real backend tests. The subsequent static check also
passed for all 368 files after the final audit/distribution utility additions.
See [software results](examples/shift/response-audits/software-checks.json) and
[the complete check log](examples/shift/response-audits/full-check.log).
Other Python versions and CI operating systems were not executed in this task.

Distribution validation uses `python tools/check.py dist`: separate wheel/sdist
builds, required archive contents and metadata, and byte reproducibility.
The final distribution run log and outcome are retained as
`reviews/shift-distribution-check.log` and `reviews/shift-distribution-check.json`
restored from the archived response worktree; these historical operational
records are outside the sdist and do not validate subsequent code changes.
Current integration checks are recorded in the decision report.
That initial phase made no commit, GitHub push or release. The final integration
commit is covered by the linked decision report.
