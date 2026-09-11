# Native multivariate OU shift inference

`nwkit shift --selection native` implements shared shift locations and shared
regime identities, with independent OU processes per trait by default and optional
[full evolutionary covariance and shared alpha](SHIFT_COVARIANCE.md). It needs no R
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
variances use original trait units squared. Full cross-trait covariance is selected
with `--trait-covariance full`; see its [model, output and simulation guide](SHIFT_COVARIANCE.md).
The scalar parameterization below describes the default diagonal model.
Non-ultrametric trees are not implemented.

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

Omit `--regime-map` to discover shifts. `--max-shifts` accepts a nonnegative
integer (default 2) or native-only `auto`. `--convergence` enables shared regime groups
and nested returns. `--search-strategy exhaustive` enumerates locations and
regime partitions subject to `--exhaustive-max-configurations` (default 5000).
`auto` uses exhaustive search if that budget admits it, otherwise the branch-pool
beam search (`lasso`). The covariance-updated AIC path below is opt-in: independent
100-tip validation improved weak-shift recovery, but exploratory tree/count
changes did not consistently improve mean prediction error.

### Automatic maximum shift count

`--max-shifts auto` starts with NWKIT's structural search limit, **N−2** for
N tips. If the complete space fits the exhaustive budget, that limit is used
without applying heuristic budgets. Otherwise, beam search resolves the cap to
`min(N−2, candidate_pool, refit_budget−1)`. Explicit `native-path` uses
`min(N−2, refit_budget−1)` because its proposals do not use the candidate pool.
These bounds preserve the existing search requirements; they do not estimate
the true number of shifts or guarantee that every shift count is explored.
Layout-specific observed-design rank and residual degrees of freedom are still
checked, including missing coordinates. AICc eligibility is evaluated per
candidate, allowing shared regimes to reduce the mean-parameter count.

For example, with 1,000 tips and the default beam budgets the cap is 128.
`--candidate-pool 24 --refit-budget 48 --screening-budget 2000` restores the
former smaller budgets. Increase budgets if the default coverage is insufficient. An integer cap is never
silently reduced: incompatible heuristic budgets remain an error.
`--search-strategy exhaustive` still errors if enumeration cannot fit its budget.
The enumeration preflight stops counting once the traversal bound is exceeded,
so a large automatic cap does not require computing enormous regime partitions.

JSON `configuration.max_shifts` preserves the requested `auto` or integer.
`search.shift_limit` records `requested`, `resolved`, the applicable `constraints`
and `budget_limited`. Replay retains the same request and complete budget
configuration; changing either invalidates reuse. Fixed-layout fitting does not
apply a discovery cap.

```bash
nwkit shift --selection native --infile dated.nwk \
  --trait traits.tsv --state-column root,leaf \
  --criterion AICc --max-shifts auto --convergence --search-strategy auto \
  --model-out fit.json --outfile regime-map.tsv
```

### Covariance-updated path search

`--search-strategy native-path` explicitly selects the native AIC/AICc path search
(without `--convergence`). It retains joint configurations from two group-lasso
paths over all branches: a fixed-root Brownian covariance seed, followed by the
best model's covariance under the requested criterion. Fixed alpha is respected, and random-root models
start from their fitted null covariance. Path columns represent OU optimum
increments and are not individually standardized. All final layouts are refitted
without shrinkage, using the existing alpha search and requested AIC or AICc formula.
This candidate restriction can improve recovery even when it does not attain
the lowest criterion score found by a broader search; it is not a global IC optimizer.

Each path has up to 80 log-spaced penalty points, with at most
`--lasso-iterations 150` iterations per point. `--screening-budget` caps the total
path points, and `--refit-budget` includes the covariance seed fit as well as
scored fits. The path stops after its active set exceeds twice `--max-shifts`;
at points exceeding the cap, the strongest effects propose a capped layout.
Finite iteration limits and cap truncation are recorded in JSON. Pool and beam
budgets apply only to `lasso`; the memory limit also guards the all-branch path.

The `lasso` strategy whitens under a fitted null covariance and uses a group-lasso
path to propose locations. Cached covariance profiles rank additions, then
unpenalized likelihoods refit retained models. Beam search includes shared and
distinct regimes and a drop/move/merge/split refinement pass. Defaults are
`--candidate-pool 128 --refit-budget 256 --screening-budget 100000 --beam-width 2
--lasso-iterations 150 --search-memory-mb 512`. These are explicit budgets;
coverage is incomplete and lasso coefficients never become final effect estimates.
These defaults are shared by the CLI and Python API. They are counts, not a
wall-clock deadline; runtime depends on tree size, traits, covariance fitting
and resampling. See [scaling measurements](NATIVE_SHIFT_SCALING.md) for scope.

By default, nested families have complexity `locations + free regime offsets`. Each family
is compared against the largest family; a plug-in parametric bootstrap repeats
the complete configured search, including data-dependent screening and covariance
fitting. Testing stops at the first non-rejection. `--calibration-replicates 199`,
`--calibration-level 0.05` and `--seed 1` are defaults. P-values include ties and
the plus-one correction. Known SE values and missingness masks remain fixed.
Failed draws abort rather than disappear from the denominator. Plug-in calibration
has no proven uniform composite-null guarantee.

`--bootstrap N` additionally simulates from the selected fit and repeats the
**entire configured selection**, including the inner bootstrap when using
bootstrap selection. This can be much
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

## Information-criterion selection

Add `--criterion AIC`, `--criterion AICc`, `--criterion BIC`, or `--criterion pBIC` to
`--selection native` to select the lowest score among all refitted candidates.
Omitting `--criterion` retains full-search bootstrap selection. Information
criteria run no calibration draws unless `--global-null-gate` is requested; `--bootstrap` still optionally estimates
selection support by replaying the selected criterion. A fixed regime map can
also be scored without searching. JSON records the criterion, score, penalty,
parameter count and each candidate score; resume provenance includes the criterion.

For K shift locations, G mean parameters per trait, q estimated covariance
parameters per trait, N tips and n observed values for a trait, AIC is
`-2 logLik + 2(K + sum(G+q))`; BIC is
`-2 logLik + K log(N) + sum((G+q) log(n))`.
Fixed parameters are excluded. Estimated parameters retain their declared count
at boundaries. Known observation variances are not estimated parameters.

AICc adds `2p(p+1)/(n_total-p-1)` to AIC, where
`p = K + sum(G+q)` and `n_total` is the sum of observed scalar coordinates across
traits. This matches kfl1ou's independent-trait convention, including counting
shared shift locations once and excluding explicitly fixed covariance parameters.
If `n_total <= p+1`, the candidate is ineligible; JSON records a null score and
`insufficient_aicc_sample_size`. AICc records its sample size and correction.
This is a conventional small-sample adjustment, not a calibrated branch-search
test or a proven finite-sample correction for arbitrary phylogenetic covariance.
The covariance-updated path uses the requested criterion to choose the covariance
for its second path, so AIC and AICc need not generate identical candidate sets.
The optional `--global-null-gate` continues to require AIC.

pBIC uses `-2 logLik + 2K log(E-1) + sum(q log(n) + logdet(I_theta) + G log(var(y)))`,
where E is the number of edges and I_theta is conditional mean information in
absolute OU-optimum coordinates. This includes the optimum-coordinate Jacobian
used by the corrected kfl1ou backend. A zero-alpha shifted fit has singular
optimum information and is ineligible. Finite and independent-tip limits are
scored at their likelihood fit, without optimizing alpha for the penalty.

These scores are research criteria, not calibrated significance tests.
Boundary fits, shared regimes, missing data and measurement error require
separate operating-characteristic validation. Both budgeted generators can miss
the global criterion minimum: beam search ranks quick likelihoods, while the
native AIC path retains joint penalized configurations before unpenalized
refitting. pBIC is evaluated for every refitted layout because its penalty varies
within a shift count. Small-alpha pBIC behavior is nonregular;
see [pBIC audit](SHIFT_PBIC.md).


## Global no-shift gate for AIC

For experimental AIC search, add `--global-null-gate`:

```bash
nwkit shift --selection native --criterion AIC --search-strategy native-path \
  --global-null-gate --calibration-replicates 199 --calibration-level 0.05 \
  --seed 1 --max-shifts 10 --infile tree.nwk --trait traits.tsv \
  --state-column x --model-out model.json --outfile regime-map.tsv
```

The statistic is `max(0, AIC(null) - minimum AIC among searched candidates)`.
Each draw simulates from the fitted no-shift model, preserves errors and missing
coordinates, then repeats the complete configured candidate search, covariance
estimation and AIC selection. Explicitly fixed parameters remain fixed. At the
requested level, rejection retains the original AIC winner; non-rejection returns
the fitted no-shift model and its predictions. Ordinary AIC scores remain unchanged.
The gate also supports native exhaustive and beam AIC search. Fixed regime maps
and other criteria/backends reject this option.

P-values are `(1 + number of bootstrap statistics >= observed)/(B + 1)`, including
numerical ties. All B draws run; a failed draw aborts the analysis. At level 0.05,
at least 19 draws are required, but that minimum gives only 0.05 resolution;
199 is the default. JSON records the statistic, every draw, p-value, decision and
ungated branch IDs. Resume fingerprints include the gate option. `--bootstrap`
selection support repeats the gate inside each outer draw, with independent seeds.

This is a plug-in test of the **global no-shift null**, with no proven uniform
error guarantee over nuisance parameters. It does not control erroneous branches
when some shifts really exist, and does not correct an overly large AIC shift
count after rejection. Treat it as experimental until independently validated
across the intended tree shapes, covariance parameters and measurement designs.
Its cost is B additional complete searches, multiplied again for optional support.
