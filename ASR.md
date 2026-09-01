# Ancestral trait reconstruction

`nwkit asr` reconstructs traits at every node and imputes missing tip values.
Discrete inference includes ER, SYM, ARD, F81, GTR, fixed-generator CUSTOM,
branch-regime Mk, and hidden-rates models. Continuous inference includes BM,
branch-rate BMS, early burst, directional drift, multivariate BM, stationary OU,
and branch-optimum OUM. The tree must be rooted under the
[shared rootedness contract](CLI_TSV_CONVENTIONS.md#rootedness-and-root-polytomies).
Root and internal polytomies, non-ultrametric trees, and finite non-negative
non-root branch lengths are supported. Root stem length is not part of any
model. No rerooting, arbitrary binary resolution, or automatic trait
transformation is performed.

## Input and automatic type selection

The TSV requires `leaf_name` and the column selected by `--state-column`.
`MV-BM` instead accepts two or more comma-separated columns, for example
`--state-column height,mass`.
Tip names are matched exactly, including numeric-looking names and literal
names such as `NA`. Duplicate names or column headers are errors. The shared
`--missing-values` and `--unmatched warn|error|ignore` policies apply.

`--trait-type auto` is the default:

- Only non-missing values at tips present in the input tree are considered.
- If all considered values parse as numbers, the trait is continuous. Non-finite
  numbers such as `inf` or overflowing exponents are errors, not categorical
  fallbacks. Recognized missing markers are removed before detection.
- Otherwise the trait is discrete. For example, `red`, `blue`, or the ambiguous
  value `0|1` select the discrete mode. A mixed numeric/text column is also
  discrete unless a token parses as a non-finite number, which is an error;
  select explicit discrete mode to retain such a literal category, or explicit
  continuous mode when unexpected text should fail.
- An all-missing column cannot be classified. Explicit discrete mode with
  `--states` can perform a prior-only analysis; continuous mode requires at
  least one observed value even when the rate is fixed.

Override with `--trait-type discrete` or `--trait-type continuous`. In particular,
numeric category codes such as `0,1,2` **require explicit discrete mode**.
Category spelling is preserved: `001`, `01`, and `1` remain different discrete
states. Supplying `--states`, `--model ER`, or another discrete-only option does
not override numeric detection; incompatible combinations fail with guidance.
The selected/requested types are reported on STDERR and in `--model-out`.
Automatic continuous selection also prints a reminder that numeric category
codes require explicit discrete mode.

| Option | Discrete | Continuous |
|---|---|---|
| `--model` | ER (default), SYM, ARD, F81, GTR, MK-REGIME, HRM, CUSTOM | BM (default), BMS, EB, BM-DRIFT, MV-BM, OU, OUM |
| `--root-prior` | model default; equal, empirical, or stationary | flat for BM-family models; stationary for OU/OUM |
| `--output` | probabilities (default), map | summary (default) |
| `--tree-annotation` | map (default), state, probability, all | summary (default), mean, all |
| Rate controls | `--rate`, `--rate-bounds`, `--rate-matrix` | `--sigma2`; OU/OUM also `--alpha`, `--alpha-bounds`, `--theta`; EB uses `--eb-rate`; BM-DRIFT uses `--drift` |
| Structure controls | `--states`, `--ambiguous-separator`, `--transition-graph`, `--regime-map`, `--regime-model`, `--hidden-categories` | `--regime-map`, `--regime-parameters`; comma-separated trait columns for MV-BM |
| Observation uncertainty | Ambiguous states | `--standard-error-column` except MV-BM |
| Interval coverage | Not applicable | `--ci-level` (default 0.95) |
| Stochastic mapping | `--stochastic-map-out`, `--n-sim`, `--threads`, `--seed` | Not applicable |

Model/output/prior defaults are selected after trait detection. The continuous
root prior describes uncertainty about a numeric root value; it is independent
of `--input-rooted`, which describes the interpretation of the tree.

## Discrete Mk models and transition structure

ER shares one rate across all allowed directed transitions, SYM fits one rate
per allowed unordered pair, and ARD fits every allowed directed rate separately.
The default `--transition-graph complete` reproduces the original complete Mk
models. `--transition-graph ordered` requires an explicit `--states` order and
allows only bidirectional transitions between adjacent states. A path instead
reads a directed TSV edge list:

```tsv
from_state	to_state
juvenile	adult
adult	senescent
```

This supports irreversible ARD models by listing only forward edges. SYM
requires a symmetric edge set. For ER, `--rate` fixes the shared rate; for SYM
and ARD it remains the common optimizer starting value. All fitted rates use
`--rate-bounds`. Fitted models use deterministic homogeneous, patterned, and
coordinate multistarts, retain the best converged likelihood, and report start
counts, failures, and lower/upper-bound rates in `--model-out`.

F81 and GTR require the complete transition graph. F81 uses target-specific
rates `q_ij = r_j`; equivalently, its equilibrium frequencies are
`pi_j = r_j / sum(r)` and its overall scale is `sum(r)`. GTR uses symmetric
pair exchangeabilities `s_ij` and fitted equilibrium-frequency ratios, with
`q_ij = s_ij * pi_j`. The first frequency weight is fixed to one to remove the
otherwise redundant frequency scale. Exchangeabilities use `--rate-bounds`;
frequency ratios use fixed numerical bounds reported in `--model-out`. F81 and
GTR default to a stationary root prior, while equal and empirical priors remain
selectable.

`--model CUSTOM --rate-matrix Q.tsv` uses a fixed labelled generator instead:

```tsv
state	x	y
x	0	0.2
y	0.4	0
```

The `state` rows must exactly match the state columns in order. State order is
inferred from this matrix unless `--states` is supplied, in which case it must
match exactly. Off-diagonal entries must be finite and non-negative. If every
diagonal is zero, nwkit derives each diagonal as the negative off-diagonal row
sum; otherwise the supplied diagonals must already satisfy that generator
constraint. CUSTOM does not fit or rescale Q.
Generator residual tolerances scale with each Q row, so a malformed very small
Q is not accepted merely because its absolute entries are small. Matrix
exponentiation repairs only roundoff-sized probability errors, with a tolerance
that scales conservatively with the exponent norm for long branches; material
negative entries or invalid row sums fail explicitly.

### Branch regimes and hidden rates

`--model MK-REGIME --regime-map regimes.tsv` jointly fits one Q matrix per
named branch regime. `--regime-model ER|SYM|ARD|F81|GTR` selects the structure
shared by those independently parameterized matrices (default ER). The map must
assign every `branch_id`, including root 0, exactly once:

```tsv
branch_id	regime
0	background
1	background
2	foreground
```

Branch IDs are the same deterministic IDs used in normal ASR output. Every
estimated regime must occur on a non-root branch. The root assignment selects
which regime Q defines a stationary root prior; it has no stem branch. Marginal
inference and stochastic mapping use the Q assigned to each incoming branch.

`--model HRM --hidden-categories H` expands each observed state into `H` latent
rate classes. Observed-state changes occur within a class, hidden-class changes
occur without changing the observed state, and every allowed expanded transition
gets its own ARD rate. Tip likelihoods sum over hidden classes; normal output and
stochastic maps are projected back to observed states, so hidden-only changes do
not appear as observed transitions. The expanded fit is subject to hidden-class
label switching and can be parameter-rich; nwkit rejects configurations requiring
more than 256 free rates.

`--root-prior stationary` derives root frequencies from the current Q, including
inside each fitted-rate likelihood evaluation. It requires a unique valid
stationary distribution; reducible generators with multiple stationary
distributions fail explicitly. Equal and empirical root priors remain available.
Marginal reconstruction, missing-tip imputation, zero-length branches,
polytomies, and stochastic transition-count mapping all use the selected graph
or fixed Q without binary tree resolution.

Uniformization caches small branch calculations, but does not retain every dense
matrix power or every state-pair event-count distribution for high-rate branches.
It instead constructs the required endpoint-specific backward vectors. A single
branch requiring more than 2,000,000 potential events is rejected before a large
allocation; reduce the rate/time scale or fitted rate bounds in that case.

## Continuous BM model

For each branch of length `t`, `X_child | X_parent` is normally distributed
with mean `X_parent` and variance `sigma2 * t`. The root has a flat prior with
respect to trait units, and its value remains uncertain rather than being
estimated and then treated as known. All-node inference uses an upward Gaussian
integration pass and a downward smoothing pass, taking O(nodes) time and memory
per rate evaluation. Estimates at an internal node condition on **all** observed
tips, not only its descendants.

If a known SE `s` is supplied, the observation has distribution
`Y_tip | X_tip ~ Normal(X_tip, s^2)`. SEs must be finite and non-negative, and
every observed value needs an SE. Missing markers are matched against the
original SE text before numeric conversion, including custom numeric markers
such as `--missing-values 999`. Missing tips need no SE. Without an SE column,
observed tips are exact; their reconstructed values equal the observations and
their conditional variances are zero. With positive SEs, the output reconstructs
the latent trait and can differ from the observed value. Its interval is not
an interval for a new noisy measurement.

`--sigma2 FLOAT` fixes a non-negative variance rate. Otherwise the rate is fitted
by REML, integrating the latent internal states and the unknown root value.
Exact-data fits use a closed form. Known-SE likelihoods can have multiple local
maxima; their rate search uses likelihood bounds over the full feasible range,
including zero, to check competing maxima to numerical tolerance. Rate-independent
within-position measurement residuals are excluded from optimization but retained
in the reported likelihood. Branch lengths are used in their original
units: the rate has units of trait-squared per branch-length unit. Multiplying
all branch lengths by `c` and dividing a fixed rate by `c` preserves inference.
Multiplying trait values and SEs by `a` multiplies the fitted rate and conditional
variances by `a^2`. Numerical centering/scaling is undone before output; it does
not normalize the evolutionary covariance to unit tip variance.

Intervals are Gaussian, equal-tail, and conditional on the fixed/fitted rate and
input tree, with coverage controlled by `0 < --ci-level < 1`. They include root
uncertainty but **exclude rate-estimation and tree uncertainty**. Positive,
bounded, or transformed traits are not detected/transformed automatically;
the analyst must choose an appropriate scale before running BM.

### Exact constraints and boundary fits

- A zero-length edge identifies its endpoints as the same latent state. The
  original nodes remain in the output. Contradictory exact observations joined
  by zero-length edges are errors, including when a positive rate is fixed.
- Identical exact observations at one zero-edge position count once for fitting
  and likelihood dimension. Noisy measurements at that position remain separate
  observations and are combined with their known SEs.
- Rate estimation requires observations at at least two distinct positions
  separated by positive-length edges. Otherwise supply a fixed `--sigma2`.
- Identical exact data can give `sigma2=0` and zero-width plug-in intervals.
  This is a boundary estimate, not evidence that the evolutionary rate is known
  to be zero. A warning and `fit_status` make this explicit. If the best positive
  rate is likelihood-indistinguishable from zero at floating-point precision,
  the zero boundary is preferred and reported.
- With measurement error, a zero evolutionary rate can retain positive
  uncertainty in the common latent value. Different exact values are impossible
  under fixed `--sigma2 0`.
- Singular feasible zero-rate limits have `fit_status=singular_zero_boundary`
  and an empty `restricted_log_likelihood`. Regular zero-rate limits have
  `fit_status=zero_boundary`; interior fits use `ok`. No epsilon rate/branch
  is substituted, and optimizer failure is an error. Unresolvable input dynamic
  ranges or an exhausted global-search budget also fail explicitly; a numerical
  lower limit is never reported as a successful rate estimate.

`fit_status` describes the numerical rate value and likelihood support, whether
the rate was fitted or fixed. Use `sigma2_estimated` or `estimation_method` to
distinguish those cases.

## Continuous OU model

`--model OU` implements a single-optimum process with a stationary root
distribution:

```text
X_root ~ Normal(theta, sigma2 / (2 * alpha))
X_child | X_parent ~ Normal(
    theta + exp(-alpha * t) * (X_parent - theta),
    sigma2 * (1 - exp(-2 * alpha * t)) / (2 * alpha)
)
```

Here `alpha` is the attraction strength per branch-length unit, `theta` is the
trait optimum, and `sigma2` is the diffusion variance per branch-length unit.
Each parameter is fixed when its corresponding `--alpha`, `--theta`, or
`--sigma2` option is supplied and fitted by ordinary ML when omitted. Stationary
OU requires strictly positive alpha and sigma2. The default positive alpha
bounds are `1e-6 / max_root_to_tip_depth` through
`50 / max_root_to_tip_depth`; override them with `--alpha-bounds MIN,MAX`.
Theta is profiled exactly from its quadratic tree likelihood in one linear pass.
The remaining zero-, one-, or two-dimensional covariance likelihood is searched
on a deterministic log grid, including exact boundaries, and polished from
competing local/boundary starts. Boundary-adjacent grid intervals are polished
before a boundary is accepted. In two dimensions, Powell fallback is reserved
for cases where all primary starts fail or a failed endpoint remains competitive
with the best converged likelihood, avoiding expensive work in inferior basins.
The free stationary variance is optimized over a reported data-scaled range
spanning 24 orders of magnitude. Alpha or stationary-variance boundary solutions
are reported in `fit_status` rather than silently treated as interior optima.
Optimizer convergence/failure and grid/start counts are retained; a coarse-grid
fallback is explicitly marked rather than reported as `ok`.

The number of distinct observed tree positions must be at least the number of
free OU parameters. Replicate noisy observations at one zero-length-contracted
position remain separate likelihood observations, but do not create additional
phylogenetic positions for this identifiability check.

OU uses the same O(nodes)-time/O(nodes)-memory upward and downward Gaussian
passes as BM for each parameter evaluation. It supports rooted polytomies,
non-ultrametric trees, exact observations, known measurement SEs, missing tips,
and exact zero-length contraction. Reported Gaussian intervals condition on the
fixed/fitted alpha, theta, sigma2, and tree; parameter and tree uncertainty are
not integrated. Unlike BM's flat-root REML quantity, OU reports an ordinary ML
log likelihood under its finite stationary root prior. The two likelihoods
therefore do not share an interchangeable AIC convention.
Fitting requires many such linear passes for grid evaluation and local polishing;
fixing alpha, theta, or sigma2 reduces that work, and fixing all three requires
only the final pruning/smoothing passes.

## Continuous model extensions

All scalar extensions retain the Gaussian all-node smoothing and conditional
interval contract described above. Known measurement SEs are supported by BMS,
EB, BM-DRIFT, and OUM as well as BM and OU.

### BMS: branch-regime Brownian rates

`--model BMS --regime-map regimes.tsv` uses variance
`sigma2_regime * t` on each incoming branch. The regime-map schema and root-row
requirement are the same as for MK-REGIME; the root regime is recorded but has
no stem variance. `--sigma2` fixes one shared rate. To fix different rates, use
a complete parameter table:

```tsv
regime	sigma2
background	0.2
foreground	1.5
```

Without either fixed-rate input, all regime rates are estimated jointly by
restricted likelihood. There must be at least one more distinct observed tree
position than fitted regimes, and the regime-specific covariance components
must be linearly distinguishable after removing the flat root. Deterministic
positive multistarts, data-scaled bounds, convergence counts, and boundary
status are reported. A fixed regime
rate may be exactly zero; estimated zero-boundary mixtures are not currently
profiled, so a lower-bound result is explicitly retained as a boundary fit.
Model-induced exact equalities from fixed zero rates retain an explicit singular
support status and do not report a density on a reduced observation space.

### OUM: branch-regime OU optima

OUM shares one positive `alpha` and one positive `sigma2`, but assigns an optimum
`theta_regime` to each branch. The transition on a branch uses its incoming-branch
regime; the stationary root distribution is centered on the root regime's
theta. `--theta` fixes one shared optimum, while a complete table fixes different
optima:

```tsv
regime	theta
cold	-1.0
warm	2.5
```

When omitted, all regime optima are estimated with alpha and sigma2 unless those
parameters are fixed separately. OUM uses ordinary stationary-root ML,
deterministic covariance multistarts, and explicit alpha/root-variance bounds.
The number of distinct observed positions must be at least the total number of
free parameters, and the regime-optimum mean design must have full column rank.
Per-regime alpha or sigma2 values are not part of this model.

### EB and BM-DRIFT

EB allows the Brownian diffusion rate to change exponentially with root depth
`d`. For a branch from depth `d` to `d+t`, its effective Brownian length is

```text
exp(eb_rate * d) * expm1(eb_rate * t) / eb_rate
```

with the continuous limit `t` at `eb_rate=0`. Negative rates describe a decline
in diffusion away from the root; positive rates describe acceleration.
`--eb-rate` fixes the exponent. Otherwise nwkit profiles it on
`--eb-rate-bounds` (default `-10/max_depth,10/max_depth`) while fitting or using
the supplied sigma2. Use the equals form for a negative lower CLI bound, for
example `--eb-rate-bounds=-1,1`. A flat profile, such as an equal-depth star
where rate and sigma2 are confounded, is rejected; bound optima are reported.

BM-DRIFT uses transition mean `X_parent + drift * t` and Brownian variance
`sigma2 * t`. `--drift` fixes the directional trend. A free trend is profiled
from observed tips at different root depths; it is confounded with the unknown
flat-prior root value on contemporaneous tips, so ultrametric observations
require a fixed drift. Both models reduce exactly to BM when their extension
parameter is zero. An exactly linear, error-free fitted trend is retained as an
explicit `sigma2=0`, `singular_zero_boundary` result rather than a spurious tiny
positive diffusion estimate.

A free drift is profiled *after* integrating the flat root; it is not integrated
as a second fixed effect. Consequently, its reported flat-root likelihood and
`residual_df` retain the `n_effective - 1` convention, rather than the
`n_effective - 2` convention of a two-fixed-effect REML analysis. Output records
this treatment explicitly, and intervals condition on the fitted drift.

### MV-BM: correlated continuous traits

`--model MV-BM --state-column trait1,trait2,...` models each branch increment as
`MultivariateNormal(0, Sigma * t)`. It estimates the evolutionary covariance-rate
matrix `Sigma` from generalized independent contrasts, then performs separable
all-node smoothing in linear tree time per trait. Primary output is long-form,
with one row per selected node and trait. `--covariance-out` writes the upper
triangle of each selected node's conditional covariance and correlation matrix;
`--model-out` records the fitted Sigma.

An observed tip must provide either every selected trait or none; partial vectors
are rejected, and arbitrary measurement-error covariance is not yet supported.
At least two distinct observed positions are needed. If there are too few
contrasts for a full-rank Sigma, marginal reconstruction remains available with
`fit_status=singular_covariance`, but the ordinary multivariate likelihood is
left empty. Intervals condition on the fitted covariance and exclude covariance
and tree uncertainty.

## Output schemas

Both modes retain shared `branch_id`, `parent` (`-1` for root), `node_class`, and
`name` columns. `--target all|intnode|leaf|missing-leaf` accepts comma-separated
classes; `intnode` includes the root. `is_imputed` identifies only missing tips,
not internal nodes or denoised observed tips. An empty selection still produces
a table header.

Continuous summary columns additionally contain:

| Columns | Meaning |
|---|---|
| `trait` | Selected input column; MV-BM emits one row for each selected trait |
| `observed_value`, `observed_se` | Original observation and its SE; empty for internal/missing nodes, SE zero for exact observations |
| `is_imputed` | Whether this is an unobserved tip |
| `mean`, `variance`, `sd` | Conditional latent-trait moments in original units |
| `ci_lower`, `ci_upper`, `ci_level` | Equal-tail conditional interval and its coverage |

For BM, `--model-out` records `trait_type`, `trait_type_requested`, the
trait/model, `root_prior`, `sigma2`, `sigma2_estimated`, `estimation_method`
(`REML` or `fixed`), `restricted_log_likelihood`, `num_observed`,
`num_effective_observations`, `residual_df`, and `fit_status`. OU instead records
alpha/theta/sigma2 values and estimated flags, stationary root variance, alpha
bounds (including the stationary-variance fitting bounds), ordinary
`log_likelihood`, separate effective-observation and distinct-position counts,
optimizer status/message/grid/start/converged/failed counts, and
`likelihood_kind=stationary_root_ml`. Both report SE-column selection and the
interval conditioning contract; parameter and tree uncertainty flags are false.

BMS and OUM additionally report regime order, root regime, source paths, each
regime parameter, optimizer counts, and data-scaled bounds where applicable. EB
and BM-DRIFT report the extension parameter, whether it was estimated, profile
search details, and the underlying sigma2 fit. MV-BM reports covariance rank and
every Sigma element under collision-safe hex-encoded trait identifiers; its
optional covariance sidecar contains readable trait names.

Discrete `--model-out` records the transition graph, `q_source` (`estimated`,
`fixed:--rate`, or `fixed:PATH`), fit/boundary status, optimizer start/convergence
counts, and every directed Q entry. CUSTOM has an empty fitted-rate-bounds field.
F81/GTR also report equilibrium frequencies; MK-REGIME reports every regime Q;
HRM reports every expanded-state Q entry.

For nondegenerate observations `y` with covariance `V`, the reported residual
log-likelihood uses the flat-root integral convention:

```text
mu_hat = (1' V^-1 y) / (1' V^-1 1)
logL = -0.5 * [(n-1) log(2*pi) + log|V| + log(1' V^-1 1)
              + (y-mu_hat)' V^-1 (y-mu_hat)]
```

An optional `log(n)/2` contrast-normalization constant is not included. Exact
zero-edge duplicates use the reduced observation space. This quantity is not
an ordinary ML likelihood or a directly comparable likelihood for different
trait scalings, observation spaces, root priors, or discrete models; do not use
it as an interchangeable AIC score. A fixed-rate fit still reports this same
residual likelihood convention.

`--tree-out` uses the shared Newick/NHX writer, preserving rooting metadata,
quoted numeric internal names, and missing-support conventions. Continuous
annotations include `asr_trait_type=continuous`; every non-BM model additionally
includes `asr_model`. BM retains its prior NHX property set. Annotation levels are:

- `mean`: `asr_mean` only.
- `summary` (default): also variance, SD, interval limits/level, and
  `asr_interval_kind` (`conditional_on_sigma2`, `conditional_on_parameters`, or
  `conditional_on_covariance`). MV-BM suffixes per-trait properties with the
  UTF-8 hexadecimal trait identifier and includes cross-trait covariances.
- `all`: also tip `asr_observed_value`, `asr_observed_se`, and `asr_is_imputed`.

Discrete primary schemas and stochastic mapping remain unchanged. Optional
discrete model metadata additionally records selected/requested trait types.
Auxiliary output paths must be distinct and cannot be STDOUT; only the primary
TSV may be written to `-`.

## Examples

```tsv
leaf_name	body_mass	body_mass_se
A	1.2	0.1
B	2.8	0.2
C	NA	NA
```

Automatic continuous detection with known measurement SEs:

```sh
nwkit asr -i '[&R](A:1,B:1,C:1);' \
  --trait traits.tsv --state-column body_mass \
  --standard-error-column body_mass_se \
  --model-out model.tsv --tree-out ancestral.nwk -o ancestral.tsv
```

Explicit continuous mode, fixed rate, exact observations, and only missing tips:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column body_mass \
  --trait-type continuous --sigma2 0.5 --target missing-leaf -o imputed.tsv
```

Stationary OU with alpha fixed and theta/sigma2 fitted by ML:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column body_mass \
  --model OU --alpha 0.8 --model-out ou-model.tsv -o ou-ancestral.tsv
```

Numeric categories must explicitly select discrete inference:

```sh
nwkit asr -i tree.nwk --trait categories.tsv --state-column state \
  --trait-type discrete --states 0,1,2 --model ER -o discrete.tsv
```

Ordered discrete states can restrict the complete Mk graph:

```sh
nwkit asr -i tree.nwk --trait stages.tsv --state-column stage \
  --trait-type discrete --states juvenile,adult,senescent \
  --model ARD --transition-graph ordered -o stages-asr.tsv
```

Regime-specific Brownian rates, estimated from the regime map:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column body_mass \
  --model BMS --regime-map regimes.tsv --model-out bms-model.tsv -o bms.tsv
```

Correlated multivariate Brownian reconstruction:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column height,mass \
  --model MV-BM --covariance-out node-covariance.tsv -o mvbm.tsv
```

Multivariate OU, threshold/liability and Lévy/jump processes,
parameter/tree-uncertainty integration, and continuous conditional trajectory
sampling are not implemented. Discrete transition-count stochastic maps are not
repurposed as continuous trajectories.

## Method references

- [ape ancestral character estimation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html): BM and REML conventions.
- [phytools fastAnc](https://search.r-project.org/CRAN/refmans/phytools/html/fastAnc.html): continuous ancestral estimates and uncertainty.
- [Hansen 1997](https://doi.org/10.2307/2411186): OU comparative models and adaptive optima.
- [Butler and King 2004](https://doi.org/10.1086/426002): multi-optimum OU models.
- [Beaulieu et al. 2013](https://doi.org/10.1093/sysbio/syt034): hidden-rate models for discrete traits.
- [Harmon et al. 2010](https://doi.org/10.1111/j.1558-5646.2010.01025.x): early-burst comparative models.

Tests compare the Gaussian passes against an independently assembled full-tree
precision matrix and tip-covariance residual likelihood, plus analytic stars,
zero-edge equivalences, and unit/offset invariance. External-package comparisons
must match rate estimation, root treatment, and interval conventions first.
A recorded `phytools 2.3.0 fastAnc` fixture checks exact-data means and variances
without adding R as a runtime or test dependency. OU tests independently assemble
the stationary patristic covariance and compare ordinary likelihood and all-node
conditional moments. Extension tests additionally compare equal-regime limits to
BM/OU, OUM smoothing to independent dense conditioning, and MV-BM covariance to
generalized independent contrasts and analytic stars.
