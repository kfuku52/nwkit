# OU optimum-shift discovery

The default is now **calibrated selection** for 4–16 tips and at most two shifts.
It repeats the full candidate search inside a parametric bootstrap and includes
both exact α limits. The no-error null test covers the entire fitted α grid. See [calibration, boundary semantics and validation](SHIFT_CALIBRATION.md).
No R installation is needed for this mode. Larger searches and the historical
IC methods require an explicit `--selection ic`; they are not calibrated by this change.

## Legacy IC backend

`nwkit shift --selection ic` is an experimental single-trait interface to the separately
installed **kfl1ou >= 3.0.9** R package. kfl1ou performs inference; NWKIT validates
inputs and converts inferred clades into its branch-regime convention. This is
not a Python reimplementation of l1ou and does not vendor GPL inference code.

Install R and kfl1ou using the [kfl1ou installation guide](https://github.com/kfuku52/kfl1ou#install).
`Rscript` must be on PATH, or supplied with `--rscript /path/to/Rscript`.
NWKIT's Python installation does not install R. Runs use `Rscript --vanilla`,
so dependencies must be available without user R startup files.

## Example

From a source checkout:

```sh
nwkit shift --selection ic -i examples/shift/tree.nwk \
  --trait examples/shift/traits.tsv --state-column value \
  --max-shifts 1 --criterion BIC --search-strategy exhaustive \
  -o regimes.tsv --model-out shift.json --fit-out shift.rds \
  --effects-out effects.tsv --regime-parameters-out optima.tsv \
  --tip-summary-out predictions.tsv
```

The initial interface supports `auto`, `lasso`, `ensemble`, and capped
`exhaustive` candidate searches, and kfl1ou's `pBIC` (default), `pBICess`, `mBIC`,
`BIC`, and `AICc` criteria. The default maximum is two shifts; specify a smaller
value for very small trees. The exhaustive configuration limit defaults to 5000.
The ensemble seed defaults to 1. Computation uses one kfl1ou worker.
These are kfl1ou candidate-search rules, not guarantees of finding the global
optimum; inspect the saved diagnostics and candidate scores.

## Input contract

- One rooted, strictly bifurcating, ultrametric tree with unique nonempty tip
  labels, at least three tips, and finite positive non-root branch lengths.
  Root-to-tip depths must agree to relative tolerance 1e-8.
- A UTF-8 TSV with `leaf_name` and the column selected by `--state-column`.
  Every tip must occur exactly once; extra rows are rejected. Row order is free.
- One finite numeric observation per tip, with at least two distinct values.
  Missing data, replicate rows and multivariate traits are not supported.
  Independent known observation errors can be supplied as described below.
- No pruning, branch repair, rerooting or tree-height normalization is applied.
  Alpha and diffusion variance rate retain the original branch-time units.
- `--max-shifts` must be nonnegative and smaller than the tip count minus one;
  kfl1ou can further restrict admissible configurations.

The tree or trait table can use `-` for stdin, but not both. Tree format and
rootedness follow the shared NWKIT options. Model averaging is not yet supported.

## Known observation errors

Add `--standard-error-column se` to read the **standard error of each tip
observation** from the trait TSV. Supply SEs in the same trait units as the
observations, not variances or the standard deviation of raw replicates. Each
SE must be finite and nonnegative, with a representable square; missing entries
and squares that overflow or underflow to zero are rejected before running R.
Use a different column from `leaf_name` and the selected trait column.

```sh
nwkit shift --selection ic -i examples/shift/tree.nwk \
  --trait examples/shift/traits-with-se.tsv --state-column value \
  --standard-error-column se --max-shifts 1 --criterion BIC \
  -o regimes.tsv --model-out shift.json --tip-summary-out predictions.tsv
```

Rows are matched by tip label, then supplied to kfl1ou as named known variances
`input_error = SE^2`, with `measurement_error = FALSE`. Thus the observation
covariance is evolutionary covariance plus `diag(SE^2)`; no additional shared
error variance is estimated. Correlated errors and replicate aggregation are
not implemented. Predictions remain fitted expected observations/latent means;
residuals include observation noise and are not standardized.

Omitting the option means zero observation variance. An explicit all-zero SE
column uses the same backend path as omission, so it cannot change the fit merely
by selecting a different optimizer. The supplied SEs and variances remain in the
tip outputs and fingerprint. JSON records `standard_error_column`, the requested
`observation_error` model and `effective_observation_error` after zero-error
canonicalization. Time-unit changes do not rescale observation variances;
trait-unit changes rescale SEs along with observations.

## Outputs

`-o/--outfile` (stdout by default) is exactly the ASR regime-map schema:

```tsv
branch_id	regime
0	baseline
```

There is one row for every input-tree branch, including root 0. A selected
branch starts `shift_<branch_id>` and descendants inherit that regime until a
nested shift starts another. By default, baseline and different shifts are
distinct regimes. With `--convergence`, shifts sharing an optimum share a regime
label, including a return to the background optimum. Branch IDs match `nwkit nwk2table` on the
same input tree, format and rootedness. Regenerate the map after editing the tree.

`--model-out` is a required JSON file (schema version 5). It records:

- NWKIT/kfl1ou versions, root model, selection criterion, search settings;
- alpha, sigma2, kfl1ou's intercept, log likelihood and selected criterion score;
- `shift_effects`, `regime_parameters`, and `tip_predictions` with original labels;
- structured `search` coverage, evaluated counts, ensemble success/failure counts
  and alpha bounds/boundary flags for unconstrained discovery
  (`search_stage`; unavailable diagnostics are null);
- selected branch IDs and unconstrained candidate configurations/scores (nonfinite
  candidate scores are JSON null, not model weights or bootstrap support);
- topology/lengths, branch regimes, tip-token-to-original-label mapping;
- a SHA-256 of the serialized analysis tree and trait table in backend tip order;
- search diagnostics as an R expression and captured backend stdout/stderr.

The intercept is kfl1ou's fitted baseline expected tip value. Reported optima use
an explicit **root mean = baseline optimum** convention. For fixed-root OU,
ultrametric tip data do not separately identify the root value and baseline
optimum; these outputs do not establish their historical equality.

The following tables are always included in JSON and optionally written as TSV:

| Option | Columns |
| --- | --- |
| `--effects-out` | `branch_id, regime, parent_regime, mean_effect, optimum_effect, optimum_identifiable` |
| `--regime-parameters-out` | `regime, branch_id, optimum, optimum_identifiable` |
| `--tip-summary-out` | `leaf_name, branch_id, regime, standard_error, observation_variance, observed, predicted, residual, optimum, optimum_identifiable` |

`mean_effect` is the additive contribution to expected extant-tip values for
all descendants of a selected branch. It is not an instantaneous ancestral jump.
`optimum_effect` is the change from the inherited optimum; effects from nested
shifts accumulate. If T is the tip height and t is the depth of the selected
branch's parent, their relation for positive alpha is
`mean_effect = optimum_effect * (1 - exp(-alpha * (T - t)))`.
`residual = observed - predicted`. Predictions are fitted tip means, not ASR
posterior means, and the tables contain no confidence intervals.

At **alpha = 0**, finite mean effects can remain defined, but finite OU optima
are not identifiable. Their values are JSON null / TSV `NA`, including the
baseline optimum, and `optimum_identifiable` is false. Positive-alpha flags
refer only to the fitted configuration and the root=baseline convention; they
are not evidence of precise historical localization. Very small positive alpha
can produce very large optima and should be interpreted using the finite mean
effects and bound diagnostics. Stationary/random-root OU at exactly zero alpha
is undefined and is rejected rather than silently relabeled as fixed-root BM.

`search.globally_optimal` is the backend's claim within its admissible search
space. Missing values remain null. `alpha_at_lower_bound` and
`alpha_at_upper_bound` compare the fitted alpha with reported bounds using
relative tolerance 1e-6 and absolute tolerance 1e-8 / tree height; an unavailable
bound yields null. Bounds and raw R diagnostics remain available for inspection.

Backend tip identities, predicted means, residuals, and the relationship between
mean and optimum effects are checked before any output is installed. Optional
TSV outputs retain headers even when zero shifts are selected. The regime-optimum
audit table is not directly an ASR `--regime-parameters` file: it includes audit
columns, and root semantics still need matching.

`--fit-out` optionally saves the complete fit as RDS for use in R. Its tree and
trait names use safe temporary tokens (`t0`, ... and `trait`); use the JSON mapping
and trait field to recover original labels. This avoids R/Newick quoting or edge
number assumptions. The R adapter returns descendant-tip sets; NWKIT maps them
back to its input branch IDs and checks the selected configuration against the
candidate profile.

Output files must be distinct and cannot replace the input tree or trait file.
They are installed together only after successful inference and validation;
handled installation failures restore previous files. This is not crash-atomic.
Output parent directories must exist. The backend's optional dependency is checked
at execution time, so other NWKIT commands do not require R.

## Convergent optima

Add `--convergence` to search for optima shared by different shift branches.
The public kfl1ou `estimate_convergent_regimes` API uses its backward heuristic,
with alpha refitted for each proposed grouping. Root treatment, original branch
units and known observation variances are retained. Convergence is disabled by
default. AICc, BIC and pBIC are supported; mBIC and pBICess are rejected before
fitting. Convergent pBIC is a heuristic extension of the unconstrained criterion.
No claim of a globally optimal convergence grouping is made.

A [joint-search reference](SHIFT_JOINT.md) found inconsistent pBIC penalties in
released kfl1ou 3.0.9. The [backend correction](SHIFT_PBIC.md) fixes coefficient
coordinates, fitted-alpha evaluation and fixed-alpha parameter counting in a
local unreleased kfl1ou checkout. Recompute pBIC analyses with that correction;
installing unmodified 3.0.9 does not include it. Statistical calibration and
convergence-model weights remain unvalidated.
An [independent sensitivity study](SHIFT_ALPHA.md) compares BIC/pBIC, two-stage
and joint selection, and alpha lower bounds on 520 new simulated datasets.

This is a refit under equality constraints on OU optima, not a relabeling of
unconstrained estimates. It does not force observed tip values to be equal.
The background may share its optimum with later shifts. The JSON `convergence`
object lists each group's NWKIT branch IDs (0 is background), its stable regime
label and the number of merges. Each background/shift identity appears exactly
once. A group containing 0 uses `baseline`; other groups use `shift_<smallest ID>`.
The regime map passes these shared labels directly to ASR.

`parameters`, effects, predictions and the RDS describe the constrained refit.
The RDS remains an `l1ou` object; `attr(fit, "nwkit.unconstrained.fit")` retains
the complete original unconstrained fit for later comparison in R.
`unconstrained_parameters`, `search`, and `candidates` preserve the original
shift-discovery stage. `regime_parameters` contains one row per shared optimum;
its `branch_id` is the smallest member ID, and all member IDs are available in
`convergence.groups`. Effects still contain one row per shift location.
The unconstrained fit or a constrained refit with alpha=0 is rejected when
convergence is requested because finite OU optima are then unidentifiable.
No positive-alpha fallback is silently substituted.

## Bootstrap selection support

Add `--bootstrap 100 --bootstrap-seed 1` to simulate traits under the fitted
OU model and repeat the configured shift search for every replicate. Known
observation variances are retained in simulation and refitting. Runs are
sequential; `--seed` still controls the search and `--bootstrap-seed` controls
bootstrap simulation. Bootstrap is disabled by default (`bootstrap: null`).

The model JSON `bootstrap` object records attempted, successful and failed
counts, summarized failure messages, `failures_without_message` for errors whose
message the backend omitted, and every successful configuration in
NWKIT branch IDs. It also reports edge inclusion frequencies, exact configuration
frequencies and tip-partition frequencies. Partitions group original tip names
by their most recent ancestral shift, ignoring regime labels: different branch
configurations can produce the same partition. Equal partitions do not by
themselves establish likelihood equivalence under every model or constraint.

All frequencies use **successful refits** as their denominator. Failures can
bias support and must be inspected; partial failure also produces a stderr
warning. If all replicates fail, the command fails
without installing any outputs. These frequencies measure repeated selection
under the fitted model and configured search, not posterior probabilities,
confidence intervals for optima or ASR uncertainty. The RDS remains the original
selected fit. A small replicate count is useful for checking execution, but gives coarse
frequencies and is insufficient for stable inference.

With `--convergence`, every replicate repeats **both** shift discovery and
backward convergence search. Simulation uses the constrained fit, and refits
retain the original search settings and SEs. The JSON records
`selection: "shift_and_convergence"`, `successful_convergence_groups` for each
successful replicate, and `shared_optimum_partition_frequencies` for tip groups
sharing an optimum. The existing `tip_partition_frequencies` still describes
ancestry through shifts, irrespective of equality constraints. Failure in either
stage, including an unidentifiable alpha=0 result, counts as a failed replicate.

## ASR handoff

```sh
nwkit asr -i examples/shift/tree.nwk \
  --trait examples/shift/traits.tsv --state-column value \
  --trait-type continuous --model OUM --regime-map regimes.tsv -o ancestors.tsv
```

This command performs a **new stationary-root OUM fit**. Shift's default is
kfl1ou `OUfixedRoot`; `--root-model OUrandomRoot` requests its random-root model.
A map alone does not transfer the fitted parameters or establish numerical
identity between engines. Match root semantics, bounds and units before comparing
fits. In particular, kfl1ou can reach the exact Brownian boundary, whereas OUM
uses positive alpha. ASR intervals conditional on the chosen map exclude
uncertainty in shift selection. No ASR refit is launched automatically.
When the shift fit used known errors, pass the same trait table and
`--standard-error-column` to ASR; the regime map alone does not transfer SEs.

## Validation and next steps

Offline adapter tests run without R. To execute actual inference and the ASR
handoff on the bundled example:

```sh
NWKIT_TEST_RSCRIPT=/path/to/Rscript python -m pytest -q tests/test_shift.py tests/test_shift_reference.py tests/test_shift_audit.py tests/test_shift_bootstrap.py tests/test_shift_convergence.py
```

The reference tests use NWKIT's shared Gaussian tree transitions and covariance
engine. On a balanced eight-tip ultrametric tree, they compare zero shifts and
two nested shifts, fixed-root alpha in {0, 1e-7, 0.5}, and stationary/random-root
alpha in {1e-7, 0.5}, both with and without heterogeneous known observation
variances (including zero-error tips). kfl1ou fits at the specified alpha/configuration; its fitted
intercept, mean effects and sigma2 are then held fixed during the NWKIT
calculation. There is no Python parameter optimization in this comparison.

The tests compare tip means, covariance (against kfl1ou's public
`sqrt_OU_covariance` plus the known diagonal observation variance), and ordinary ML log likelihood. Mean/covariance tolerances
are rtol 1e-7 and atol 1e-8. Absolute likelihood tolerance is 1e-7 except for
small-alpha stationary roots (1e-5), whose covariance is ill-conditioned.
Rescaling every branch by 1000 and alpha/sigma2 by 1/1000 must preserve the
same observables. A hand-calculated BM covariance and explicit rejection of the
stationary alpha=0 boundary provide additional independent checks. The installed
kfl1ou 3.0.9 also fails when asked to fit random-root OU with alpha fixed at zero.

These tests validate the covered fixed-parameter models; they do not demonstrate
equivalence of default ASR fits, optimizers, search paths or arbitrary datasets.
The dense reference evaluator is for small validation examples and is not invoked
during production shift inference. Model averaging remains a next step; a native search engine still requires
representative end-to-end performance evaluation.

Additional audit tests reject contradictory diagnostic counts/coverage and
nonfinite or unrepresentable error inputs, verify tip-order independence,
protect existing files on installation failure, and check that tiny trait units
do not defeat result consistency checks. Remaining branch times are computed
from descendants, avoiding subtraction of nearly equal large root depths.

Bootstrap audit checks compare all 22 configurations with up to two shifts on a
four-tip tree against `kfl1ou::shift_tip_partition`, and compare branch support
with a direct call to `l1ou_bootstrap_support` from the saved RDS under both root
models with known SEs. They also cover no-shift replicates, partial failure,
missing error messages, malformed diagnostic quoting and output preservation
when R aborts after fitting.

Convergence integration tests cover fixed and random roots with and without
known SEs, independent same-parameter tip-mean and log-likelihood evaluation,
shared regimes on disconnected branches, nested returns to background,
bootstrap reproducibility and ASR consumption of the shared-regime map.

Result checks allow rounding at the scale of source effects and intermediate
values along the relevant ancestral path. This covers cancellation of large
opposing shifts without letting a large effect on an unrelated branch hide an
incorrect equality constraint. Search alpha bounds are exported from the
unconstrained discovery fit, matching the stage named by `search_stage`.
After validating the shift-effect reconstruction, directly fitted tip optima
are retained and used consistently within each regime. This avoids replacing
a small fitted optimum with cancellation error from large opposing deltas.
Regimes without extant tips retain their reconstructed optimum; background
always uses the fitted intercept under the stated root convention.

A reproducible independent simulation harness and a small paired BIC/pBIC
pilot are documented in [SHIFT_VALIDATION.md](SHIFT_VALIDATION.md). It measures
false shift selections, recovery of shifts and shared-optimum groups, and
bootstrap failures under both root treatments and observation SEs. Broader
calibration across tree sizes and signal strengths should precede model averaging. kfl1ou's current `model_average_l1ou`
rejects convergent fits because an unconstrained candidate profile does not
preserve the equality constraints; convergence-aware averaging needs its own
validated candidate set.
