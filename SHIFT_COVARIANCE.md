# Evolutionary covariance in native shift analysis

`nwkit shift --selection native --trait-covariance full` estimates evolutionary
covariance across traits while sharing shift locations/regime assignments. One
covariance matrix applies across regimes. This is an experimental native model;
numerical validation is not a general guarantee of shift-detection calibration.

```sh
nwkit shift --selection native --input-rooted yes \
  --infile dated.nwk --trait expression.tsv --state-column root,leaf \
  --trait-covariance full --alpha-model shared --max-shifts 2 \
  --model-out shifts.json --outfile regimes.tsv
```

Omitting `--criterion` invokes fitted-model parametric bootstrap calibration,
with 199 complete searches per test by default. `--criterion AIC`, `AICc` and
`BIC` are available research scores. Joint-covariance `pBIC` is rejected because
its appropriate extension has not been implemented. The defaults remain
diagonal trait covariance and trait-specific alpha.

## Model and estimation

- `--trait-covariance diagonal|full`: diagonal diffusion or a full positive
  semidefinite evolutionary diffusion matrix. Full covariance adds p(p−1)/2
  correlation parameters to the p marginal process variances.
- `--alpha-model shared|trait-specific`: one attraction strength or one per
  trait. The attraction matrix itself is diagonal; cross-trait attraction is
  not estimated. `--alpha` fixes original-time rates instead of estimating them.
- `--standard-error-column se_root,se_leaf` supplies known sampling errors.
  `--estimate-measurement-error` adds estimated diagonal observation variances.
  Missing trait coordinates are allowed. Sampling and evolutionary covariance
  are distinct; the shift fitting interface currently accepts diagonal sampling
  errors, not an unrestricted estimated observation covariance.
- `--process-tip-variance` can fix marginal process variances while estimating
  correlations. `--measurement-variance` fixes additional diagonal errors and
  cannot be combined with `--estimate-measurement-error`.
- `--covariance-engine auto|pruning`: `auto` uses an exact separable covariance
  profile with shared alpha, complete observations, no observation error and
  free process variances. Otherwise, problems with at most 384 observed
  coordinates and a conservative 64 MiB workspace estimate use observed-coordinate
  dense GLS with analytic covariance and mean-design gradients. Larger inputs
  use vector-tree square-root pruning, reusing the observed subtree/index plan.
  Both paths retain checked multistart optimization and an independent numerical
  score check after convergence. Dense optimization uses direct measurement
  variance coordinates so a zero starting error variance can leave the boundary;
  its numerical score check cancels finite-difference truncation bias without
  changing the convergence tolerance. `pruning` forces tree GLS and numerical gradients
  for validation. Fixed and stationary roots, missing coordinates, known errors,
  estimated diagonal measurement variance, and common zero/infinite alpha limits
  retain their existing definitions. Engine metadata identifies `dense_observed_gls`,
  `vector_tree_pruning`, or `separable_profile`. No covariance jitter is added.
  Analytic gradients can change local optimization paths and heuristic branch
  selection; mathematical likelihood equivalence does not imply identical search
  results on every dataset.

With A=diag(alpha_i), diffusion D and branch duration t, the transition is
F=diag(exp(-alpha_i*t)) and Q_ij=D_ij*(1-exp(-(alpha_i+alpha_j)*t))/(alpha_i+alpha_j).
Zero rates use the corresponding integral limit. A fixed root and a stationary
root imply different tip covariance. Diffusion is parameterized to remain valid
when attraction rates differ; an arbitrary stationary covariance must not be
mistaken for a valid diffusion parameterization.

The common Brownian limit (fixed root) and common independent-tip limit are
evaluated explicitly. Estimated trait-specific alpha does not enumerate every
mixed zero/finite boundary face; diagnostics report that incomplete coverage.
Mixed finite/infinite rates with nonzero cross-trait diffusion are unsupported.
At the independent limit, free process and diagonal measurement variances are
not separately identifiable; output records the canonical all-process allocation
and its ambiguity. Small marginal identifiability diagnostics are conservative
and do not certify full joint covariance identifiability.

When a complete error-free dataset has too few residual dimensions for a full
covariance estimate, the corresponding layout has no finite unconstrained ML
optimum and is explicitly excluded. No ridge or numerical jitter is silently
added. Shrinkage and low-rank covariance estimators are future, separate models;
this implementation reports unpenalized covariance fits.

## Search and computation

The existing exhaustive, beam (`auto`/`lasso`) and covariance-updated
`native-path` strategies support joint covariance. Screening whitens the response
and candidate effects with their joint covariance; it never starts by dropping
branches using only a diagonal-covariance filter. Retained candidates refit
their means and covariance without the screening penalty. Finite candidate and
refit budgets still make heuristic search approximate.

The square-root tree factor stores O(nodes*p²) state, rather than forming the
(tips*p)-square observed covariance. For fixed p and bounded tree degree its
factorization uses O(nodes*p³) small-matrix work; mean profiling and numerical
optimization add costs. The shared-alpha fast path profiles covariance exactly
and reuses a bounded cache of tree factors within each search. Observation errors
must not be discarded to obtain that fast path.

Bootstrap generates correlated evolutionary innovations, keeps the supplied
sampling errors and missing mask, and repeats candidate generation, covariance
fitting and selection. Estimated nuisance parameters make this a plug-in
procedure, not an exact composite-null test. Stability bootstrap also uses the
correlated generator. There is no oracle knowledge of generating covariance in
this operational procedure.

## Outputs

Native model JSON retains its branch/effect/regime/tip tables and additionally
exports `joint_covariance` with:

- `process_tip_covariance` in original trait units squared;
- `diffusion_covariance` in original trait units squared per original tree time,
  or null at the independent limit;
- estimated/fixed additional `measurement_covariance` (diagonal);
- full mean-coefficient covariance, in trait-major then regime-coefficient order;
- shared versus trait-specific alpha, parameter counts, actual engine,
  convergence diagnostics and identifiability limitations.

The top-level likelihood is joint across all traits. Per-trait `log_likelihood`
is null, since these are not independent additive contributions. Covariance
options participate in provenance/resume fingerprints. The matrix estimates do
not supply confidence intervals or propagate model-selection uncertainty.

For joint AIC/BIC/AICc, shared alpha counts once and covariance parameters count
once. The joint AICc/BIC sample-size convention uses observed tip vectors, not
independent scalar traits. The AICc correction is a heuristic, especially with
missing data and estimated phylogenetic dependence; it is not a multivariate
finite-sample calibration theorem. Penalized estimators cannot reuse this
parameter count unchanged.

## Simulation

`nwkit shift-simulate` generates unconditional datasets from explicit JSON
parameters or a completed native model on the same tree. It writes observed
traits, generating truth and optionally latent states at every node.

Example `parameters.json`:

```json
{
  "trait_names": ["root", "leaf"],
  "alpha": 1.0,
  "root_model": "OUfixedRoot",
  "diffusion_covariance": [[2.0, 1.6], [1.6, 2.0]],
  "shift_branch_ids": [],
  "regime_optima": [[0.0, 0.0]],
  "sampling_standard_errors": [0.1, 0.1]
}
```

```sh
nwkit shift-simulate --input-rooted yes --infile dated.nwk \
  --parameters parameters.json --replicates 1 --seed 17 \
  --outfile simulated.tsv --truth-out truth.json --latent-out latent.tsv

nwkit shift-simulate --input-rooted yes --infile dated.nwk \
  --model-in shifts.json --replicates 100 --seed 19 \
  --outfile simulations.tsv --truth-out fitted-truth.json
```

For multiple replicates, the TSV has repeated tip names distinguished by the
zero-based `replicate` column: select one replicate before passing it to `shift`.
For one replicate the TSV is directly usable, with trait columns and known
sampling-SE columns `se_<trait>`. Additional measurement errors are recorded in
the truth/model, not falsely labeled as known sampling errors in those columns.

Explicit parameter fields:

| Field | Meaning |
|---|---|
| `trait_names` | Ordered distinct trait names |
| `alpha` | Physical-time scalar or vector; string `"inf"` for the independent limit |
| `root_model` | `OUfixedRoot` (default) or `OUrandomRoot` |
| `diffusion_covariance` | Symmetric PSD matrix in physical trait/time units |
| `process_tip_covariance` | Alternative to diffusion; required at the independent limit. Unequal finite alpha must imply valid PSD diffusion |
| `shift_branch_ids` | Distinct non-root NWKIT branch IDs; default no shifts |
| `groups` | Optional partition of root 0 and shifts into shared regimes |
| `regime_optima` | One row per canonical group (sorted by minimum branch ID), one column per trait |
| `scaled_regime_coefficients` | Alternative means: baseline row then scaled offsets; supports the zero-alpha drift limit |
| `sampling_standard_errors` | Scalar, per-trait vector, or tip-by-trait array |
| `measurement_covariance` | Optional common PSD observation covariance, including cross-trait errors for simulation/stress tests |
| `missing` | Optional tip-by-trait boolean mask |

Tip order follows the input tree's leaf order. The truth JSON records this order,
branch topology/lengths, generating covariance, means, errors, missingness and seed.
Keep the input tree and parameter/model JSON to reproduce the same draw.
An unresolved variance decomposition in a fitted model requires explicit
generating parameters rather than inventing separate variances.

The simulation API is `ShiftSimulation` / `simulate_shift`; internal fitted-model
bootstrap uses `simulation_from_fit`. Generation shares model semantics with
inference, while tests independently reconstruct Gaussian means/covariances and
compare dense likelihoods to the pruning and separable implementations.

## Validation

The [implementation validation report](reviews/trait-covariance-2026-09-11/README.md)
records container checks, actual native timings and the limited fitted-null pilot.
The [shared versus trait-specific alpha comparison](reviews/shared-versus-specific-alpha-2026-09-11/REPORT.md)
estimates alpha in both models and separates error-free timing, observation-error
timing and a limited AIC/BIC discovery pilot.

See `tests/test_vector_whitening.py`, `tests/test_shift_joint_covariance.py` and
`tests/test_shift_joint_simulation.py`. They cover dense likelihood agreement,
diagonal nesting, parameter counting, original-unit transformations, missingness,
correlated simulation moments, actual CLI round trips, safe output publication,
joint screening and full-search bootstrap.

`tools/benchmark_shift_covariance.py` records actual native timing and a 100-tip
fitted-null bootstrap pilot. Its fixed-alpha calibration experiment estimates
covariance/means but does not validate estimated-alpha selection. Compare
recorded configurations, failures and numerical equivalence before interpreting
timings. The earlier GeneGalleon prototype used generating-null oracle calibration
and is a different experiment.
