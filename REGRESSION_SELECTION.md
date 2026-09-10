# Exploratory phylogenetic predictor selection

`nwkit regress-select` fits lasso or elastic-net models and evaluates their
predictions using nested, user-supplied group cross-validation. It supports
numeric Gaussian, binary (0/1), Poisson and negative-binomial responses.
Use `nwkit regress` for prespecified coefficient tests. Selection results do
**not** contain coefficient P-values, adjusted P-values, or confidence intervals.

```bash
nwkit regress-select \
  --tree examples/regression_selection/tree.nwk \
  --data examples/regression_selection/data.tsv \
  --response trait --predictor-file examples/regression_selection/predictors.txt \
  --family gaussian \
  --folds examples/regression_selection/folds.tsv \
  --strengths 1,0.1,0.01 --l1-ratios 1,0.5 \
  --out-prefix results/selection
```

Use the kebab-case option names shown above. The common CLI also accepts the
compatibility aliases `--predictor_file`, `--evolution_model`, `--l1_ratios`,
and `--out_prefix`, with a deprecation warning.

## Data and scope

- `--data` has one row per species, a unique `leaf_name` key, one selected
  response, and finite numeric predictor columns. Missing values are rejected.
- Name predictors with `--predictors OG1,OG2,...` or a one-name-per-line
  `--predictor-file`; the latter avoids command-line limits on large matrices.
- `--unpenalized` names a subset of predictors, such as prespecified confounders,
  to retain without a penalty. The intercept is always unpenalized.
- `--family` is `gaussian` (default), `binomial`, `poisson`, or
  `negative-binomial`. Binomial means numeric 0/1, with 1 the modeled event.
  Count responses must be non-negative integers. Family is never inferred
  from integer storage. Copy number used as a predictor does not determine
  the response family.
- Brownian selection requires a rooted tree. Rooting declarations are respected;
  use `--input-rooted yes` for an intentionally rooted polytomy, and `no` to
  declare an unrooted input (which Brownian selection rejects).
- `--evolution-model` is `brownian` (default) or `independent`. The tree supplies
  covariance shape; evolutionary shape parameters are not tuned in this command.
  The covariance is normalized to mean diagonal one before CV, using tree
  information only, so changing branch-length units does not change the grid.
- The current backend is dense and supports at most 500 tips. Predictor count
  may exceed tip count. Replicate uncertainty, categorical predictors, offsets,
  structural-zero models and reconciled contrasts are not yet supported here;
  several of these are available in the unpenalized `regress` command.

## Objective and optimization

Every training fit centers and scales predictors using its training rows only.
Scaling avoids underflow/overflow from squared predictor units, and prediction
uses those centered coordinates to avoid cancellation from large offsets.
Constant penalized columns are excluded for that fit and receive coefficient
zero. Unpenalized columns must remain full-rank with the intercept and fewer
in number than training tips; constant unpenalized columns are an error.

For standardized coefficients, the penalty is

```
strength * [l1_ratio * sum(abs(beta)) + (1-l1_ratio)/2 * sum(beta^2)]
```

Only penalized coefficients enter these sums. `l1_ratio=1` gives lasso;
`0 < l1_ratio < 1` gives elastic net. Strengths must be strictly positive.
Gaussian loss is `(y-Xb)' C^-1 (y-Xb) / (2*n)`, in the supplied response units.
The Gaussian residual scale is estimated from the final residuals after fitting;
it does not change this convex GLS elastic-net objective. It is a descriptive
ML residual-scale estimate, not an unbiased post-selection variance estimate.
Gaussian fitting uses quasi-Newton initialization followed by exact coordinate
minimization of the same convex GLS objective. Unpenalized terms are projected
out after whitening and recovered by QR. Final KKT checks are independent of
the initializer's relative-objective stopping flag, which can differ across
BLAS/runtime environments. Gaussian fits require a final maximum KKT violation
of at most `1e-6`, with a default budget of 20,000 iterations (2,000 for
non-Gaussian fits). Exhausting that budget does not waive the KKT check.
Both phases share the iteration budget; storage is
linear in predictor count, without a dense predictor-by-predictor Gram matrix.
For pure Gaussian lasso, a rank-deficient equicorrelation design is flagged as
`lasso_equicorrelation_design_rank_deficient`: the coefficients may not be
unique. Training loss can agree while held-out predictions and CV choices
differ across numerical environments. This is not repaired by clipping
predictions or silently adding ridge regularization.

For non-Gaussian families, loss is the Laplace-approximated negative marginal
log likelihood divided by training tip count. The model has random effect
`u ~ N(0, variance*C)` and logit (binary) or log (count) link. Random variance
and, for negative binomial, NB2 dispersion are fitted jointly with coefficients
by ML. Negative-binomial variance conditional on `u` is `mu + dispersion*mu^2`.
The existing NWKIT GLMM mode solvers are reused. Log variance and dispersion
are bounded to [-12,6]; a solution at a bound is explicitly flagged.

Nonnegative positive/negative coefficient parts represent the L1 penalty without
smoothing it. Fixed-linear Laplace gradients include the derivative of the
mode and determinant. Only the one/two variance parameters use numerical
finite differences. Output requires a finite, converged solution and a checked
projected gradient. These non-Gaussian fits are local optima, not a guarantee
of the global minimum. Outer-CV stability and boundary diagnostics matter.
If relative objective change stops L-BFGS before stationarity, up to two
curvature resets are attempted within the original iteration budget. The
objective, bounds and projected-gradient acceptance threshold are unchanged.

## Nested group validation

`--folds` is a TSV with exactly one row per tree tip:

```text
leaf_name	fold
species_A	clade_1
species_B	clade_1
species_C	clade_2
```

Provide at least three phylogenetic groups. In each outer iteration, one group
is held out. On the remaining groups, leave-one-group-out inner CV chooses a
strength and L1 ratio. That choice is fitted on the outer training rows and
used to predict the untouched outer group. All training splits need at least
four tips and a variable response; binary splits must contain both 0 and 1.
Unusable splits fail explicitly. Groups need not have equal size.

The grid is scored by tip-weighted binary log loss or squared prediction error
(Gaussian and count responses). Failed candidates are not compared using a
subset of their folds: they are ineligible. If every candidate fails in any
required tuning stage, no bundle is published. Exact score ties prefer the
stronger penalty, then the larger L1 ratio.

Each outer fold also evaluates a baseline with the same tree and prespecified
unpenalized covariates but no selectable predictors. Its loss helps assess
predictive benefit over phylogeny alone. Failed baseline fits are explicitly
marked and suppress the aggregate baseline score.

After outer evaluation, a separate group-CV pass over all rows tunes the final
full-data model. Outer prediction errors never choose its hyperparameters.
The strength and ratio grids are explicit; no response-based screening outside
training folds is performed.

`--prediction conditional` (default) adds the conditional random-effect mode
`C_test,train C_train,train^-1 u_hat_train` to the fixed linear predictor.
Gaussian `u_hat_train` is the training residual. `--prediction fixed` uses fixed
effects only. Both are **plug-in predictions**: the inverse link is applied
without integrating uncertainty in latent effects or coefficients. Count
predictions are therefore not claimed to be integrated marginal means.

Clade folds evaluate transport to withheld clades within the supplied tree.
They do not make training and test species statistically independent. The
conditional mode uses known tree relationships but never held-out trait values;
fixed-effect prediction omits that contribution. Choose the mode according to
the intended prediction task and record it when reporting performance.

## Output bundle

`--audit audit.json` records all three input files (and the predictor-name file
when used), the tree summary and all six output artifacts. Audit paths cannot
alias an input or a bundle member.

The command publishes these files together, preserving previous outputs if a
fit or handled write fails:

| Suffix | Contents |
|---|---|
| `.coefficients.tsv` | Final full-data coefficients in original units, standardized coefficients, training scales, selected/constant flags, penalty settings and optimizer diagnostics |
| `.path.tsv` | Full-data coefficients at every requested grid point; exploratory, not CV evaluation |
| `.cv.tsv` | Inner and final-tuning fold scores, failures and boundary warnings |
| `.predictions.tsv` | Outer held-out observations, predictions, baseline predictions/errors and selected hyperparameters |
| `.stability.tsv` | Fraction of outer fitted models selecting each penalized predictor |
| `.metadata.json` | Response family, tree model, loss, prediction interpretation, dimensions, grid and nested-CV mean loss |

`selected` means a penalized standardized coefficient has absolute value
above `1e-7`; intercept and unpenalized covariates are always retained but are
not labeled selected. Selection frequency is descriptive and is not an FDR
bound or formal stability-selection guarantee. Inspect correlated predictor
sets, not just a single selected representative.

Every coefficient row records
`inference_status=exploratory_no_post_selection_inference`. Refitting the
selected variables with ordinary `regress` and applying BH to those P-values
does not account for selection. Existing fixed-model bootstrap inference also
does not do so. A future inferential extension must incorporate selection and
tuning, define the hypotheses being tested, and validate error control under
phylogenetic simulations before it emits inferential claims.

## Validation

`tests/test_regression_selection.py` checks analytic Laplace gradients against
numerical derivatives, exact orthogonal lasso solutions, the small-penalty GLS
and existing-GLMM limits, correlated predictors with more columns than tips,
training-only preprocessing, held-out-response isolation, invalid splits,
input/output alias protection and the output bundle. Runtime integration with
GeneGalleon additionally exercises mixed response families and publication
rollback when a later trait fails.
