# Independent OU selection and alpha-bound sensitivity

This research experiment compares BIC and corrected pBIC, two-stage and joint
selection, and two lower bounds on alpha using the same independently generated
data. This is historical research evidence, not the settings of the current
[calibrated CLI default](SHIFT_CALIBRATION.md). It used corrected pBIC and the
explicit bounds recorded below. Current reruns additionally use the shared
[backend capability check](SHIFT_PBIC.md).

The [technical report](examples/shift/alpha-validation/report.html) contains the
results, uncertainty, diagnostics, and interpretation. The
[protocol](examples/shift/alpha-validation/protocol.json) was written before any
primary fits. This study uses new datasets, unlike the earlier correction replay.

## Frozen design

There are 520 independent datasets and eight selected models per dataset:

- Primary: eight tips, four truths (null, single shift, distinct optima,
  convergent optima), two root treatments, 50 replicates per cell: 400 datasets.
- Extensions: weak attraction, known observation error, and sixteen tips;
  null and convergent truths under both roots, 10 replicates per cell: 120 datasets.

Trees are balanced and bifurcating with unit branches. With height `H`, generating
parameters are `alpha * H = 2.1`, `sigma2 * H = 0.25`, and nonzero optima `+2` or
`-2`. Weak attraction uses `alpha * H = 0.2`; known observation error uses tip
standard error 0.2. Other cells have zero observation error. The generator draws
independent branch innovations and uses a fixed root at zero or a stationary
random root, matching the fitted root treatment. It does not call the R backend
to simulate data.

Each dataset uses `SeedSequence([20260911, cell_id, replicate])`. The complete
input set is saved before inference. The bounds are `alpha * H >= 1e-7` or `0.1`,
with a common upper bound of 10 and starting value 1. The former backend default
upper bound excluded the previous generating alpha, so this experiment explicitly
uses a wider upper bound containing both generating values. Comparisons with
older pilots therefore change more than replication count.

All fits allow at most two shifts. Two-stage selection uses exhaustive
unconstrained branch search followed by backward convergence. Joint selection
enumerates the [identifiable candidate space](SHIFT_JOINT.md): 195 models for
eight tips and 899 for sixteen. Models are fitted through public kfl1ou APIs.
One fit per candidate provides pBIC and the public `stats::BIC()` result, avoiding
criterion-dependent duplicate numerical fits. The raw R object for a joint BIC
selection therefore still carries the pBIC criterion; the recorded `score` is
the selected criterion and `reported_score` preserves the object's original
score. Native R objects are retained in the raw run, not shipped in the evidence
bundle. No external inference implementation is copied into NWKIT.

## Definitions and audit

An effective shift changes the optimum group relative to its nearest selected
ancestor. A retained branch merged into the ancestral group does not count as
a false positive. Both retained and effective counts are saved. This definition
differs from historical pilots that counted retained branches.

Under the null, false positive means at least one effective shift. Shared-optimum
recovery requires an exact, label-invariant match of the complete tip partition.
Exact effective branch recovery is a separate metric. Mean-prediction RMSE is
computed against the generating expected tip means, excluding observation noise.
Root integration is included in the random-root generating expectation.

Rates retain their numerators, completed-fit denominators, Wilson 95% intervals,
and bounds obtained by treating all failed fits as either successes or failures.
Paired comparisons retain both attempted and completed pairs and count improvements
and losses separately. Results are stratified by design cell; the eight fits on
one dataset are not eight independent observations.

Every selected model has its tip mean and ordinary Gaussian log likelihood
independently recomputed with NWKIT's tree-process evaluator. BIC is additionally
checked from this likelihood and the declared parameter count. Mean, likelihood,
and BIC absolute tolerances are `1e-6`; equality of optima within each declared
group uses absolute `1e-6` or relative `1e-10` for very large optima. Candidate
ledgers verify complete unique IDs and the selected joint score minimum.
The audit records boundary proximity within a relative `1e-4` of either bound.

The ledger also refits the same selected baseline model where present in the
joint candidate set. Its score difference is stored separately from the gap to
the joint minimum. Redundant baseline histories may be absent from the joint
candidate set. Candidate failures, numerical refit differences, or excluded
baseline histories preclude a certified complete search comparison for that case.
No claim of a global continuous optimum is made even when these checks pass.

A separate, explicitly post-hoc diagnostic checks the null model in all 100
primary null datasets under both bounds. It profiles the intercept and variance
analytically, evaluates 161 log-spaced alpha values, and refines every grid-local
maximum. It does not change any selected model or the frozen design. The script
is `tools/shift_alpha_null_diagnostic.py`; its result is saved as
`null-profile-diagnostic.json` in the evidence bundle.

## Evidence and reproduction

The backend is the locally corrected, unreleased kfl1ou 3.0.9 described in
[SHIFT_PBIC.md](SHIFT_PBIC.md). A behavioral probe rejects the original coordinate
mismatch. The installed package file hashes identify the exact library used;
the version number alone is insufficient. The run saves a source snapshot before
fitting, and the evidence retains compressed input/truth, selected-model metrics,
and every candidate attempt, including failures.

From NWKIT, after installing the corrected backend in an isolated R library:

```bash
R_LIBS=/tmp/kfl1ou-pbic-library PYTHONPATH=. python tools/validate_shift_alpha.py \
  --output /tmp/shift-alpha-run --rscript Rscript
PYTHONPATH=.:tools python tools/summarize_shift_alpha.py \
  --input /tmp/shift-alpha-run --output /tmp/shift-alpha-evidence
python tools/report_shift_alpha.py /tmp/shift-alpha-evidence
```

The first command refuses an existing output directory. Raw runs retain all
selected RDS models, tip/effect tables, backend logs, and full candidate ledgers.
The report script creates canonical `artifact.json` and `report-data.json`.
Rendering `report.html` additionally requires the Data Analytics plugin's
portable report delivery tool; the statistical evidence and reproduction do not.
That renderer requires SQL provenance, so the report script materializes its
reviewed aggregates in an in-memory SQLite table and executes the exact read
recorded in report metadata. It does not alter the statistical results.

## Limits

This is a finite simulation experiment, not a calibrated significance test.
Fifty primary replicates give limited rate precision; ten extension replicates
are exploratory. Wilson intervals are pointwise, not simultaneous across cells.
The study does not cover nonbalanced trees, multiple traits, more than two shifts,
bootstrap coverage, or model-average calibration. A lower-bound sensitivity result
is not authorization to choose a bound after seeing real-data results. Any tuning
based on these data requires another independent validation set.
