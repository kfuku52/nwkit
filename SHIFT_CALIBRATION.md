# Calibrated small-tree OU shift selection

`nwkit shift` now defaults to `--selection calibrated`. This addresses the severe
small-sample overselection found in the [IC/alpha audit](SHIFT_ALPHA.md) by
calibrating the **entire search**, including nuisance parameters and shared
optimum candidates. Raising an arbitrary lower bound on alpha is not the fix.
The historical procedure remains available with explicit `--selection ic`.

```sh
nwkit shift -i examples/shift/tree.nwk \
  --trait examples/shift/traits.tsv --state-column value \
  --convergence --max-shifts 2 --calibration-replicates 199 --seed 1 \
  -o regimes.tsv --model-out calibrated.json \
  --effects-out effects.tsv --regime-parameters-out optima.tsv \
  --tip-summary-out predictions.tsv
```

The native Python implementation requires 4–16 tips, a rooted binary ultrametric
tree with positive branches, and at most two shifts. It exhaustively enumerates
identifiable location configurations; `--convergence` also enumerates their
nonredundant equality partitions. All enabled candidates count against
`--exhaustive-max-configurations` (default 5000). Inputs outside this scope fail
explicitly. There is no automatic fallback to IC selection.

## Selection rule

Let Q contain orthonormal contrasts with Q1=0. Each candidate is fitted to Qy
using its Gaussian contrast density. This removes the unknown overall mean and
the common root component. On an ultrametric tree, fixed-root and stationary-root
OU covariance differ by a constant rank-one matrix, which Q removes. Thus the
contrast selection rule is the same for either `--root-model`. This does not
make the original root models identical or estimate their root uncertainty.

The nested families are: no shifts; at most one shift; at most one nonbaseline
effect (including two locations sharing an effect, when convergence is enabled);
and the full candidate family. Duplicate families are removed. At each stage:

1. Fit the best candidate in the current family and compute twice the log-density
   improvement of the full family over that family.
2. Generate B datasets from the fitted current-family model. With no measurement
   error in later stages, bootstrap process variance uses residual degrees of freedom d−q,
   where d=n−1 and q is the number of independent mean contrasts. Likelihood
   evaluation itself uses its ML variance, RSS/d.
3. Refit **every location, equality and nuisance-grid candidate** on each dataset.
4. Compute p=(1 + number of simulated improvements at least as large as observed)
   /(B+1). Stop at the first family with p above `--calibration-level`; otherwise
   continue, ending at the full family.

For the first no-error test, NWKIT now evaluates the null distribution at every
alpha in the fitted grid, using common Gaussian draws. Baseline and process
scale cancel from the improvement statistic, so alpha is the only nuisance
parameter at this stage. Rejection requires p <= level at **every** grid point.
If an evaluated p already exceeds the level, acceptance is certain: NWKIT stops
and reports `p_value=1`, `p_value_kind=conservative_upper_bound`, and the partial
maximum as `p_value_lower_bound`. This is a computational bound on the grid
maximum, not a confidence interval or a posterior probability. If all grid
points are evaluated, `p_value_kind=grid_supremum` and the maximum is reported.
The evaluated points and probabilities are saved in `null_alpha_evaluations`.

At a true alpha belonging to the grid, this no-error construction is
conservative: the reported probability is at least the Monte Carlo p-value
computed under that true null covariance. A finite grid **does not** prove
uniform control for arbitrary continuous alpha. Known-error null tests and all
later stages retain explicitly labelled plug-in calibration. In particular,
the guarantee does not extend to selecting the exact number of true shifts.
The general nuisance-supremum principle is discussed by
[Berger and Boos (1994)](https://doi.org/10.1080/01621459.1994.10476836);
this implementation uses a finite grid, not their confidence-set procedure.

Defaults are B=199 and a stagewise level of 0.05. `--seed` determines calibration
randomness. At least 19 replicates and sufficient Monte Carlo resolution for the
requested level are required. The combination of finite-grid and plug-in
calibration is **not a proof of uniform 5% finite-sample error control**. Later tests and selected
parameters are also subject to selection uncertainty. `--bootstrap` is the old
IC selection-support feature, not this calibration; mixing it or `--fit-out`
with calibrated selection is rejected. `--criterion` requires `--selection ic`.

## Alpha limits and parameter meaning

Profiles use dimensionless a=alpha×tree height H, on the recorded grid containing
exact 0, exact infinity, and 25 log-spaced finite points from 0.001 to 1000.
This is a finite-grid approximation, not continuous optimization. Ties within
1e-10 log-density units prefer the exact limits, then the smaller finite value.
All bootstrap searches use exactly the same grid.

The fixed-root process covariance is parameterized by its marginal tip variance
v. At a=0 it is v times the shared-ancestry matrix divided by H (Brownian limit);
at infinity it is vI. For finite alpha,
v=sigma2×H×(1−exp(−2a))/(2a). Mean shift weights are normalized by 1−exp(−a),
so their limits remain finite. At zero they become the fraction of tree height
since each shift's parent. These are scaled-effect/drift limits: ordinary finite
OU optima are not identified there. At infinity the limiting tip means and v
are estimable, but finite alpha and diffusion rate are not separately estimated.

Schema-7 JSON records `parameters.alpha_status` as `finite`, `brownian_limit`
or `independent_limit`. Infinity is encoded as JSON null, with its explicit status;
no nonstandard NaN/Infinity literals are written. At the independent limit,
`sigma2` is also null. `process_tip_variance` and mean predictions remain available.

The entire selected-candidate alpha profile is exported. A drop of at most
1.920729410347062 log-density units defines a diagnostic support set. This is
**not** a confidence interval after model selection. If either alpha limit belongs
to this set, finite optima and optimum effects are withheld (`null`/TSV `NA`,
`optimum_identifiable=false`), even when the best grid point is finite. This is
an operational weak-identification flag, not a mathematical proof of
identifiability when it is true. Shared regime labels at the Brownian limit
encode constraints on scaled effects, not evidence for equal finite OU optima.

## Known observation errors

`--standard-error-column` supplies fixed independent variances D=diag(SE²).
They are never absorbed into an estimated overall scale. For each alpha,
calibration fits Q(vK+D)Qᵀ, including admissible v=0. Positive v values use 49
log-spaced points from 10⁻⁶ to 10⁶ times the median positive observation variance.
The grid is recorded in JSON and is independent of the trait observations.
Reaching its upper boundary in any fitted candidate or bootstrap dataset raises
an error instead of returning a clipped estimate. Zero v is skipped only when
its contrast covariance is singular. An all-zero error column exactly follows
the analytic no-error path. The grid is approximate; close comparisons warrant
sensitivity analysis, especially with very unequal errors.

Predicted absolute means use the fixed-root representative covariance to recover
the intercept. Adding a stationary-root common component leaves this point
estimate unchanged; root uncertainty is not exported. Exact zero residual
variance makes the no-error Gaussian likelihood undefined and is rejected.

Some mixtures of zero and positive errors do not admit a finite regular
likelihood maximum. With m exact observations and at least m−1 allowed shifts,
a candidate can fit every exact observation while process variance tends to
zero, making the contrast density diverge. These configurations are rejected
before bootstrap calibration. Other mixed-error fits are rejected if any
candidate reaches the smallest positive variance grid point at a singular
zero-variance boundary. An arbitrary variance or SE floor is not substituted.
One exact observation with all other errors positive still has a nonsingular
contrast error covariance and supports an admissible zero process variance.

Mean coefficients use QR least squares instead of normal equations. Bootstrap
residual tensors contain at most 64 replicates at a time; the smaller d×B normal
draw array is retained to preserve the seeded draw sequence.

## Outputs and validation

The branch/regime TSV retains the ASR interface. Auxiliary TSV columns retain
the existing names, with unidentified optima explicitly missing. ASR still
refits its own model and does not inherit the boundary model or selection
uncertainty. Schema 7 uses `contrast_log_likelihood`, not a BIC/pBIC score.
Schema 5 belongs to the explicit legacy IC backend.

The current [700-dataset validation](examples/shift/calibration-envelope/README.md)
uses seed 20260918 and a frozen implementation snapshot. It covers primary
8-tip null, single, distinct and shared effects, plus weak pull, known errors
and 16 tips. All datasets are generated branch by branch independently of the
fitting kernel. Separate studies contain [120 random-tree null datasets](examples/shift/calibration-envelope-stress/summary.json)
(seed 20260919) and [1,000 weak-pull null datasets](examples/shift/calibration-weak-null/README.md)
(seed 20260920). These studies use the alpha-envelope correction.

The preceding `calibration-validation`, `calibration-review` and
`calibration-stress` directories are historical snapshots of the earlier
plug-in method. Their results are retained, not overwritten or presented as
validation of the current correction. In particular, the review study's weak
null cell selected shifts in 7/50 datasets; that small-cell result alone cannot
distinguish nuisance-estimation bias from sampling variation.

| Current validation | False selections | Denominator |
|---|---:|---:|
| Main 700-dataset study: all null cells | 7 | 250 |
| Main study: primary 8-tip null | 3 | 100 |
| Dedicated weak-pull null study | 45 | 1,000 |
| Random-tree/unequal-error null study | 6 | 120 |

In the main study, correct tip partitions were recovered in 54/100 single-shift,
44/100 distinct-shift and 44/100 convergent-shift primary datasets. Weak-pull
convergence recovery remained 0/50. Empirical rates and Wilson intervals in the
saved summaries apply only to the tested cells. Pooling them does not establish
uniform error control. For example, the random-tree known-error subset returned
4/60 false selections; all four occurred among the 20 known-error 16-tip cases.
This small subgroup needs targeted follow-up and is not declared calibrated.

The evidence verifier reconstructs all selection/recovery flags and summaries,
refits every observed candidate and alpha profile, checks probability-bound
metadata, and optionally replays the complete seeded bootstrap traces and fit
parameters. Frozen snapshots are always hash-checked. `--allow-source-revision`
requires `--replay-bootstrap`; it records every changed source hash and must
reproduce the saved results. This was used to verify that the subsequent
mixed-zero-error guard leaves all supported simulation cases unchanged.

Reproduce into **new** directories:

```sh
OPENBLAS_NUM_THREADS=1 python tools/validate_shift_calibration.py \
  --output /tmp/shift-envelope --seed 20260918 --workers 4
OPENBLAS_NUM_THREADS=1 python tools/verify_shift_calibration.py \
  /tmp/shift-envelope --replay-bootstrap --workers 4
OPENBLAS_NUM_THREADS=1 python tools/stress_shift_calibration.py \
  --output /tmp/shift-stress --seed 20260919 --workers 2
OPENBLAS_NUM_THREADS=1 python tools/validate_shift_calibration.py \
  --output /tmp/shift-weak-null --focus weak-null \
  --extension-replicates 500 --seed 20260920 --workers 4
```

A [worked CLI example](examples/shift/calibrated/model.json) exports the regime
map, mean predictions, alpha diagnostics and missing optimum estimates.
The [review and next work plan](reviews/shift-calibration-review.md) records the
reproduced bugs, fixes, remaining limits and priorities.

## Follow-up known-error diagnostics

The [paired known-error study](examples/shift/known-error-diagnostic/README.md)
compares fitted-nuisance bootstrap calibration with an oracle using the true
generating mean and covariance, separately for no-shift and one-shift families.
It also documents a research-only finite alpha/variance envelope that reports
unresolved computations explicitly. Neither experiment changes the CLI default
or establishes continuous-parameter or later-stage error control.
