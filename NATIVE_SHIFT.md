# Native multivariate OU shift inference

`nwkit shift --selection native` implements shared shift locations and shared
regime identities with independent OU processes for each trait. It needs no R
backend. Fixed layouts support parameter estimation; discovery and calibration
are **research-only** pending the [adoption protocol](NATIVE_SHIFT_VALIDATION.md).
The existing default `--selection calibrated` is a separate, restricted method.

## Inputs and fixed layouts

Supply a rooted, binary, ultrametric tree with positive branch lengths and unique
named tips. No branch flooring, arbitrary polytomy resolution, tree repair or
missing-tip pruning is performed. Traits are TSV columns; `NA` coordinates are
allowed with at least three observed tips per trait. Known sampling SE columns
must be aligned one-for-one with trait columns. Each layout needs positive
residual degrees of freedom and full observed mean-design rank for every trait.

```bash
nwkit shift --selection native --infile dated.nwk \
  --trait traits.tsv --state-column root,leaf \
  --standard-error-column root_se,leaf_se \
  --regime-map assigned-regimes.tsv --estimate-measurement-error \
  --model-out fit.json --outfile regime-map.tsv \
  --effects-out effects.tsv --regime-parameters-out optima.tsv \
  --tip-summary-out predictions.tsv
```

The map contains `branch_id` and `regime` for **every** branch including root 0.
IDs use NWKIT level-order traversal of the supplied tree. Regime names may be
reused on disconnected or nested branches. The root mean equals the background
optimum. `--root-model OUfixedRoot` fixes the centered root; `OUrandomRoot` uses
a stationary Gaussian root and excludes alpha zero.

Each trait has its own attraction, process variance, optional additional
observation variance and regime coefficients. `--estimate-measurement-error`
estimates the extra variance in addition to supplied SE squared. Fixed
`--alpha`, `--process-tip-variance` and `--measurement-variance` accept either one
number or a comma-separated value per trait. Alpha uses original time units;
variances use original trait units squared. Full cross-trait covariance and
non-ultrametric trees are not implemented.

## Parameterization and numerical scope

Internal time is divided by tree height H; trait centering/scaling is reversed
on export, including the ordinary likelihood Jacobian. Let a = alpha * H and
v be the process variance at a tip. For a branch of normalized length t, the
transition slope is exp(-a*t). Fixed-root innovation variance is
`v * (1-exp(-2*a*t))/(1-exp(-2*a))`; stationary-root innovation variance is
`v * (1-exp(-2*a*t))` with root variance v. Observation covariance adds supplied
SE squared plus the extra estimated variance. All trait likelihoods are added.

Regime offsets are scaled by `1-exp(-a)`. Alpha zero is the explicit
fixed-root Brownian/drift limit; infinity is the independent-tip limit. These
limits are not replaced with small/large finite estimates. Ordinary Gaussian
ML profiles means using square-root tree whitening and QR, without constructing
the tip covariance matrix. Finite a is numerically searched over `[1e-6,1000]`;
variance searches use normalized `[1e-10,1e4]` and separate exact zero components.
The JSON records convergence and bound diagnostics. No continuous global
optimum is certified. Selection aborts on unresolved covariance modes or a
selected numerical variance bound.

At infinite alpha, unknown process and unknown observation variances cannot be
separated: both exports are null, with their estimable sum in
`unstructured_tip_variance`. Optima are NA at alpha limits, zero process variance,
numerical bounds, or when an evaluated alpha limit is within 1.920729 log
likelihood units. This is a conditional fixed-layout diagnostic, **not** a
post-selection confidence interval or a general proof of identifiability.
An additional local covariance-Jacobian rank diagnostic uses the distinct
observed MRCA depths, with normalized columns and rank tolerance `1e-8`. It can
detect loss of variance/alpha separation in small or nearly independent trees.
Unsupported variance components and sigma squared are withheld, with
`fitted_total_tip_variance` retaining the fitted sum. This is conservative:
additional information in the mean is not assumed, and a passed local check
does not prove global identifiability.
For an estimated alpha unsupported by this diagnostic (including zero process
variance), the primary alpha fields are null and the numerical maximizer remains
in explicitly named `alpha_candidate` fields. Supplied fixed parameters are
distinguished from estimates.
Scaled offsets and fitted tip means remain available. Model-averaged inference
and full covariance-component identifiability analysis remain outside this release.

## Discovery, calibration and stability

Omit `--regime-map` to discover shifts. `--max-shifts` is an explicit cap, not an
automatic small-tree/large-tree rule. `--convergence` enables shared regime groups
and nested returns. `--search-strategy exhaustive` enumerates locations and
regime partitions subject to `--exhaustive-max-configurations` (default 5000).
`auto` uses exhaustive search if that budget admits it, otherwise `lasso`.

The `lasso` strategy whitens under a fitted null covariance and uses a group-lasso
path to propose locations. Cached covariance profiles rank additions, then
unpenalized likelihoods refit retained models. Beam search includes shared and
distinct regimes and a drop/move/merge/split refinement pass. Defaults are
`--candidate-pool 24 --refit-budget 48 --screening-budget 2000 --beam-width 2
--lasso-iterations 150 --search-memory-mb 512`. These are explicit budgets;
coverage is incomplete and lasso coefficients never become final effect estimates.

Nested families have complexity `locations + free regime offsets`. Each family
is compared against the largest family; a plug-in parametric bootstrap repeats
the complete configured search, including data-dependent screening and covariance
fitting. Testing stops at the first non-rejection. `--calibration-replicates 199`,
`--calibration-level 0.05` and `--seed 1` are defaults. P-values include ties and
the plus-one correction. Known SE values and missingness masks remain fixed.
Failed draws abort rather than disappear from the denominator. Plug-in calibration
has no proven uniform composite-null guarantee.

`--bootstrap N` additionally simulates from the selected fit and repeats the
**entire calibrated selection**, including its inner bootstrap. This can be much
more expensive than fitting or calibration alone. Reported branch frequencies
are stability under a fitted model, not posterior probabilities or test p-values.
The JSON also records exact-layout frequency and regime-pair co-selection/shared
regime frequencies among the selected model's branches, including separate
unconditional and conditional-on-co-selection denominators.
`--bootstrap-seed` controls an independent stream.

## Artifacts and replay

Schema-8 JSON includes original-unit parameters, likelihoods, branch assignments,
tip predictions, candidate scores, calibration draws, support and research status.
Optional TSVs are long-form by trait. Missing or unsupported estimates use NA.
The branch-map TSV supplies assignments only: passing it to ASR refits its model
and does not transmit alpha limits or selection uncertainty.

`--resume-model previous.json` regenerates requested artifacts only when input,
complete configuration and implementation fingerprints match. The implementation
fingerprint includes Python sources and numerical dependency versions. This reuses
a **completed** analysis; interrupted bootstrap checkpoints are not supported.
Different input coordinates, errors, layouts, configuration or source code reject
reuse. Publication is transactional and protects input files and old outputs.

`tools/benchmark_native_shift.py` compares equivalent dense and tree fixed-layout
GLS for balanced or pectinate trees. Those measurements do not establish an
end-to-end search or kfl1ou speedup.

The 1,000-tip/100-shift development target and reproducible measurements are
described in [the scaling guide](NATIVE_SHIFT_SCALING.md).
