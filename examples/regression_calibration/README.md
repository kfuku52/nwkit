# Regression calibration evidence — 2026-09-10

See [protocol and usage](../../REGRESSION_CALIBRATION.md) for the estimands,
generators, model assumptions, acceptance bands and limitations. These are
fixed-model experiments; they do not validate post-selection inference.

## Completed experiments

| Directory | Generated datasets | Purpose |
|---|---:|---|
| `smoke-20260910` | 102 | 51 cases × 2; diagnostic only, bootstrap B=4 |
| `wald-pilot-20260910` | 10,200 | 51 cases × 200, original optimizer, Wald/oracle |
| `bootstrap-pilot-20260910` | 600 | RSC 2/5/20 events × 200, B=999, paired methods |
| `tip-pilot-20260910` | 1,400 | 7 raw-tip/reconciliation/replicate cases × 200 |
| `glmm-fixed-pilot-20260910` | 4,800 | 24 GLMM cases × 200, unpenalized multistart fix |
| `repeated-bootstrap-pilot-20260910` | 40 | Two repeated-contrast cases × 20, B=199 |
| `rsc-confirm-20260910` | 10,000 | Raw 8-tip RSC and 7-event contrast control × 5,000 |

All seven experiments completed. `audit-20260910.json` verifies every saved
input hash and seed, unique/complete replicate sets, applicable methods,
source archive hashes and exact summary regeneration for 27,142 records.
This count includes paired inputs across experiments, not 27,142 independent
samples. The 102 smoke records are not calibration evidence.

Each directory contains `protocol.json`, `records.jsonl.gz`, `summary.json`
and `source.tar.gz`. The protocol records the exact case definitions, seeds,
methods, runtime versions, source hashes and completion state. Source snapshots
are evidence for these runs, not upstream version defaults. Experiments used
GeneGalleon Docker `local/genegalleon:standard-iqtree-dev`, image
`deca052d1d71`, with frozen read-only NWKIT source overlays and single-threaded
BLAS. No SIF/Apptainer or rebuilt distribution image was validated.

## Findings

* With independent Gaussian errors and one contrast per event, the existing
  coefficient bootstrap rejected 28%, 11.5% and 7% at 2, 5 and 20 events;
  its 95% percentile coverage was 72.5%, 88.5% and 93.5%. Wald rejection was
  3%, 5% and 6%; validation-only null bootstrap was 2.5%, 5% and 6%.
  These use 200 outer datasets and B=999. At 2 events, the bootstrap rejection
  Wilson interval is 22.2–34.6%. The exact infinite-B limit in this simple
  case is `2*t.sf(norm.ppf(.975), events-1)`, about 30.03% for 2 events.
  Refitting variance alone does not studentize a coefficient bootstrap.
* A 9% raw-tip pilot rejection rate did not replicate: an independent seed
  with 5,000 datasets gave 243/5,000 = 4.86% rejection (MC95 4.30–5.49%)
  and 95.14% coverage (94.51–95.70%). The 7-event contrast control gave
  4.66% and 95.34%. Both meet the prespecified rate/coverage bands in these
  known Brownian, single-copy conditions. This does not establish validity
  with duplication, measurement errors or uncertain trees.
* Successful unpenalized GLMM optimizer results could have worse likelihood
  than a nested intercept-only fit. The production fix evaluates the existing
  multiple starts even after success. Two exact reproductions are regression
  tests. This fixes a numerical defect without proving a global optimum.
* `paired-glmm-20260910.json` checks all 4,800 before/after input hashes.
  Rare binary n=8 retains only 62/200 available P values. NB2 n=30 changes
  from 199 to 197 available values and 19/199 to 16/197 rejections; post-fix
  coverage is 181/197. The fix is not a general small-sample calibration cure.
* `laplace-reference-20260910.json` and
  `laplace-reference-fixed-20260910.json` compare independent importance-QMC
  integration at fitted parameters. Some stable NB2 likelihood differences
  are about 0.1–0.2. The extreme binary example fails reference convergence;
  its numerical discrepancy is not a reliable reference estimate. These
  checks do not establish reference optima, standard errors or coverage.

Alternative-effect cases, including copy-count measurement-error cases with
beta=0.5, are not Type I error experiments. Report availability separately
from conditional rejection and coverage. The repeated-contrast bootstrap
pilot has only 20 outer datasets per case and cannot establish calibration.

## Historical evidence caveats

The diagnostic smoke source used NB alpha as size in the validation-only
null bootstrap; this was corrected before later reference runs. Do not use
its NB null-bootstrap output as scientific evidence. The Wald pilot does not
call that reference generator and is unaffected by that error.

The raw-tip pilot predates preflight classification for leaves with no
biological observations. Its seven partial-missing input errors are stored
as `fit_failed`; they are ineligible input datasets, not optimizer failures.
Its protocol scope text also predates the raw-tip extension; saved cases and
source define what actually ran. Raw CLI inner-bootstrap attempts are not
instrumented, as recorded; direct-engine runs provide that accounting.

## Reproduce or audit

Restore an experiment's archived source in a separate directory and use a
GeneGalleon container. Reuse the saved `cases`, master seed, replicate range,
methods, B and timeout from `protocol.json` with
`tools/validate_regression_calibration.py`; pass case definitions through
`--case-file` and write to a new output directory. Confirmation used master
seed 20260912, pilots 20260911 and smoke 20260910. Do not replace historical
records with reruns. Use a fresh seed for additional confirmation.

`tools/verify_regression_calibration.py` accepts experiment directories and
`--output` for a new audit file. `tools/compare_regression_calibration.py`
accepts `--before`, `--after`, `--output` and refuses mismatched inputs.
`tools/plot_regression_calibration.py --results <this-directory> --output <stem>`
creates PNG/PDF summaries from the recorded estimates and MC intervals.

Related runtime tests passed: 234 existing non-slow tests, 6 slow tests,
and 20 calibration tests. Maintainability hard limits passed. This is not
an assertion that the full NWKIT suite or a release build was run. An isolated
source distribution build and inclusion check passed for the new documentation,
validation scripts, test and audit; independent wheel reproducibility was not run.
