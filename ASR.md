# Ancestral trait reconstruction

`nwkit asr` reconstructs traits at every node and imputes missing tip values.
Discrete inference includes ER, SYM, ARD, F81, GTR, fitted rate-class designs,
Pagel independent/dependent binary-trait models, fixed-generator CUSTOM,
branch-regime Mk, hidden/covarion rates, across-character rate mixtures, and a
Bayesian threshold/liability model. Continuous inference includes BM,
Pagel lambda/kappa/delta transforms, EB/ACDC, directional and regime drift,
multivariate BM/shared-OU/diagonal-attraction OU, flexible-root OU, and
OUM/OUMA/OUMV/OUMVA regime models.
The tree must be rooted under the
[shared rootedness contract](CLI_TSV_CONVENTIONS.md#rootedness-and-root-polytomies).
Root and internal polytomies, non-ultrametric trees, and finite non-negative
non-root branch lengths are supported. Root stem length is not part of any
model. No rerooting, arbitrary binary resolution, or automatic trait
transformation is performed.

## Input and automatic type selection

The TSV requires `leaf_name` and the column selected by `--state-column`.
Multivariate Gaussian models instead accept two or more comma-separated columns, for example
`--state-column height,mass`.
Pagel discrete models require exactly two comma-separated binary-trait columns.
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
  `--states` can perform a prior-only analysis only with a fully fixed process:
  ER plus `--rate`, or CUSTOM plus `--rate-matrix`. Fitted ER, SYM, ARD, F81,
  GTR, MK-REGIME, and HRM parameters have no likelihood information and are
  rejected. COVARION, MK-MIXTURE, THRESHOLD, and fitted continuous models are
  also rejected without informative observations. Continuous mode requires at
  least one observed value even when the rate is fixed.

Override with `--trait-type discrete` or `--trait-type continuous`. In particular,
numeric category codes such as `0,1,2` **require explicit discrete mode**.
Category spelling is preserved: `001`, `01`, and `1` remain different discrete
states. Supplying `--states`, `--model ER`, or another discrete-only option does
not override numeric detection; incompatible combinations fail with guidance.
Discrete ASR requires a model state space containing at least two states. If an
invariant sample observes only one category, provide the biologically valid
larger state space explicitly with `--states`; a one-state CTMC has no estimable
transition process and is rejected consistently by every discrete model.
The selected/requested types are reported on STDERR and in `--model-out`.
Automatic continuous selection also prints a reminder that numeric category
codes require explicit discrete mode.

| Option | Discrete | Continuous |
|---|---|---|
| `--model` | ER (default), SYM, ARD, F81, GTR, MK-DESIGN, PAGEL-INDEPENDENT/DEPENDENT, MK-REGIME, HRM, COVARION, MK-MIXTURE, THRESHOLD, CUSTOM | BM (default), BMS, BMS-DRIFT, LAMBDA, KAPPA, DELTA, EB, ACDC, BM-DRIFT, MV-BM, MV-OU, MV-OU-DIAG, OU, OUM/OUMA/OUMV/OUMVA |
| `--root-prior` | equal/empirical/stationary for CTMCs; identified Gaussian for THRESHOLD | flat for BM-family models; stationary for OU-family models; OU also supports fixed or Gaussian |
| `--output` | probabilities (default), map | summary (default) |
| `--tree-annotation` | map (default), state, probability, all | summary (default), mean, all |
| Rate controls | `--rate`, `--rate-bounds`, `--rate-design`, `--rate-matrix`; mixture/covarion controls | `--sigma2`, transform parameters, `--alpha`, `--alpha-by-trait`, `--theta`, and `--drift` as applicable |
| Structure controls | states/graph/design/regime/hidden controls; multiple columns for MK-MIXTURE or Pagel | regime map/parameters; multiple trait and SE columns for MV-BM/MV-OU/MV-OU-DIAG |
| Observation uncertainty | Ambiguous states; THRESHOLD MCMC | Known per-trait measurement SEs for every continuous model |
| Interval/parameter uncertainty | THRESHOLD posterior probabilities and liability moments | Conditional Gaussian intervals; optional transform profile CI, joint posterior samples, PPC, and bootstrap |
| Simulation | CTMC stochastic maps except MK-MIXTURE/THRESHOLD | `--posterior-samples-out`, `--posterior-predictive-out`, `--bootstrap-out`, and `--seed` |

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
Structured models are rejected before optimization when they would require more
than 256 free transition parameters; the same total-parameter guard applies to
MK-REGIME and the base generator of MK-MIXTURE. ER transition probabilities use
a cancellation-safe formula, so positive rates remain positive even far below
the precision at which `exp(-x)` can be distinguished from one.

F81 and GTR require the complete transition graph. F81 uses target-specific
rates `q_ij = r_j`; equivalently, its equilibrium frequencies are
`pi_j = r_j / sum(r)` and its overall scale is `sum(r)`. GTR uses symmetric
pair exchangeabilities `s_ij` and fitted equilibrium-frequency ratios, with
`q_ij = s_ij * pi_j`. The first frequency weight is fixed to one to remove the
otherwise redundant frequency scale. Exchangeabilities use `--rate-bounds`;
frequency ratios use fixed numerical bounds reported in `--model-out`. F81 and
GTR default to a stationary root prior, while equal and empirical priors remain
selectable.

`--model MK-DESIGN --rate-design design.tsv` fits an arbitrary partition of
direct transition rates. The TSV must contain exactly these columns in order:

```tsv
from_state	to_state	rate_class
juvenile	adult	forward
adult	senescent	forward
senescent	adult	reverse
adult	juvenile	reverse
```

Only listed directed edges are allowed, and every edge sharing a `rate_class`
uses one fitted rate. State names must belong to the inferred or explicit model
state space; self edges, duplicate edges, empty classes, and an empty design are
rejected. The design therefore expresses ordered, irreversible, symmetric, or
other biologically motivated rate-sharing hypotheses without supplying fixed
rate values. `--rate` is a common optimizer starting value, while
`--rate-bounds` constrains every class. A separate `--transition-graph` is
invalid because the design already specifies both edges and sharing. Model
metadata reports the design path, class order, fitted class rates, and full Q.

`PAGEL-INDEPENDENT` and `PAGEL-DEPENDENT` model the correlated evolution of two
binary traits supplied as `--state-column trait1,trait2`. The four joint states
form one CTMC in which exactly one trait may change per instantaneous event;
simultaneous two-trait jumps are excluded. The independent model has four rates
(two directions for each trait, shared across the other trait's background).
The dependent model has eight rates because each direction may depend on the
other trait's state. Each column must resolve to exactly two model states. To
include an unobserved state or fix ordering, use
`--states 'A0,A1;B0,B1'`. If one tip trait is missing, its two compatible joint
states remain in that tip's likelihood, so information from the observed trait
is retained. Standard output, tree annotations, and stochastic maps use four
unambiguous JSON-encoded joint-state labels. The two Pagel models share
`likelihood_kind=pagel_joint_ml` and can be compared to each other; this
joint-tip likelihood is kept separate from the across-character MK-MIXTURE
likelihood.

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
exponentiation repairs only roundoff-sized probability errors. Long branches
are exponentiated at a moderate time scale and repeatedly squared with
stochastic-row validation, avoiding generic `expm` row-sum drift while
retaining slow modes. Material negative entries or invalid row sums fail
explicitly.

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
estimated regime must occur on a positive-length non-root branch whose
descendant subtree contains an informative observation. A regime confined to
zero-length branches or wholly missing subtrees has no likelihood information
and is rejected. The root assignment selects which regime Q defines a
stationary root prior; it has no stem branch. Marginal inference and stochastic
mapping use the Q assigned to each incoming branch.

`--model HRM --hidden-categories H` expands each observed state into `H` latent
rate classes. Observed-state changes occur within a class, hidden-class changes
occur without changing the observed state, and every allowed expanded transition
gets its own ARD rate. Tip likelihoods sum over hidden classes; normal output and
stochastic maps are projected back to observed states, so hidden-only changes do
not appear as observed transitions. The expanded fit is subject to hidden-class
label switching and can be parameter-rich; nwkit rejects configurations requiring
more than 256 free rates or 64 expanded states.

`--model COVARION` uses the same expanded state space with an identifiable,
parsimonious parameterization: ordered log-spaced hidden rate multipliers have
geometric mean one, observed-state changes share a base rate, and hidden classes
switch at one fitted rate. This removes arbitrary class-label and scale
confounding. Every effective hidden-class observed-transition rate is constrained
to `--rate-bounds`; a fixed `--rate` must lie strictly inside those bounds so the
spread remains identifiable. Model metadata reports both multipliers and effective
rates, and boundary-saturated fits are excluded from regular AIC/AICc/BIC ranking.
Standard output marginalizes hidden classes and stochastic maps project away
hidden-only changes.
To bound dense matrix-exponential cost, observed states times hidden classes may
not exceed 64 expanded states in either hidden model; the HRM 256-rate guard is
likewise evaluated before allocating its expanded transition graph.

`--model MK-MIXTURE --state-column c1,c2,...` jointly fits at least two
characters under one ER/SYM/ARD/F81/GTR base generator. `--rate-mixture gamma`
(default) estimates a mean-one discrete-gamma shape using equal-probability
bins represented by their conditional means; `free` estimates ordered mean-one
category rates and weights. `--rate-categories` accepts 2 through 8.
Each character gets its own posterior rate-category probabilities before node
states are averaged over categories. Primary output is stacked by `trait`.
Because the mixture is across characters rather than along branches, a single
branchwise stochastic map is undefined and rejected. Transition matrices are
computed once per category/unique branch length and reused across characters.

`--model THRESHOLD --states low,medium,high` treats the ordered categories as
intervals of a latent Brownian liability. Scale and location are not jointly
identifiable from categories, so the process is fixed to unit diffusion with
`X_root ~ Normal(0,1)`. Binary data use threshold zero. Ordinal data fix the
first threshold at zero and sample the remaining ordered thresholds; all states
must then be observed. `--thresholds` instead fixes every threshold. Data-
augmentation MCMC samples all node liabilities, with retained draws, burn-in,
thinning and chains controlled by `--liability-samples`,
`--liability-burnin`, `--liability-thin`, and `--liability-chains`.
Ambiguous ordered observations constrain liabilities to the union of their
allowed category intervals and remain valid as estimated thresholds move.
`--liability-out` reports posterior liability moments; `--model-out` reports
maximum R-hat, minimum lag-one ESS, and diagnostic status. Threshold inference
requires positive branch lengths and does not report a marginal likelihood or
support CTMC stochastic mapping.
R-hat requires at least two chains and two retained draws per chain; otherwise
it is reported as unavailable rather than as apparent convergence. Burn-in may
be zero, while retained draws, thinning, and chain counts must be positive.

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
allocation, and one endpoint history is capped at 256 MiB after accounting for
the expanded state count; reduce the state/rate/time scale or fitted rate bounds
in that case. A requested simulation set is also rejected when its conservative
preparation-plus-sampling bound exceeds 2,000,000 uniformization steps; reduce
`--n-sim`, rates, or the rate/time scale instead of starting an unbounded job.

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
- A fitted `sigma2_lower_boundary` or `root_variance_lower_boundary` denotes a
  nonregular zero-variance-component limit. Its numerical likelihood remains a
  diagnostic, but `asrcompare` and legacy IC comparison exclude it because the
  result depends on an artificial positive optimizer bound and regular-model
  AIC/AICc/BIC assumptions do not apply.

`fit_status` describes the numerical rate value and likelihood support, whether
the rate was fitted or fixed. Use `sigma2_estimated` or `estimation_method` to
distinguish those cases.

## Continuous OU model

`--model OU` implements a single-optimum process. Its default stationary root
distribution is:

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

`--root-prior fixed --root-mean M` instead fixes the root state. A proper
nonstationary prior is selected with
`--root-prior gaussian --root-mean M --root-variance V`. These modes fit alpha,
sigma2 and theta by ordinary proper-root ML through the shared affine-Gaussian
engine; the user-supplied root variance is not multiplied by sigma2.
Stationary-root OU rejects root mean/variance options. The biological root
assumption is therefore explicit rather than silently forcing equilibrium at
the beginning of a short or nonstationary tree.
Before optimization, free proper-root OU parameters are checked for local rank
through their induced observed-tip mean vector and covariance matrix. For
example, free alpha and sigma2 on an equal-depth fixed-root star are rejected
because only one variance combination is identified.

## Continuous model extensions

All scalar extensions retain the Gaussian all-node smoothing and conditional
interval contract described above. Known per-tip measurement SEs are supported
by every scalar model; all multivariate models accept one comma-separated SE
column per trait and allow trait-level missingness.

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

### BMS-DRIFT: regime rates and directional trends

`--model BMS-DRIFT` uses
`X_child = X_parent + drift_regime*t + error` with
`Var(error)=sigma2_regime*t`. A complete `--regime-parameters` table has
`regime`, `sigma2`, and `drift` columns. Without it, `--sigma2` or `--drift`
can fix one shared component while the other component is estimated per regime;
omitting both estimates both maps. Before optimization, nwkit verifies full rank
of the root-to-tip regime-time design for drift and the pairwise path-regime
design for diffusion. This catches ultrametric/global-drift and unrepresented-
regime confounding explicitly. The model uses the same flat-root integrated
likelihood convention as BM/BMS.

### OUM/OUMA/OUMV/OUMVA: branch-regime OU

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
deterministic covariance multistarts, and explicit alpha/diffusion bounds.
The number of observed tips must exceed the total number of free parameters;
non-identifiable regime designs should instead use fixed parameters.
The related models vary additional branch parameters while retaining a
stationary root under the root regime:

- OUM: theta varies; alpha and sigma2 are shared.
- OUMA: theta and alpha vary; sigma2 is shared.
- OUMV: theta and sigma2 vary; alpha is shared.
- OUMVA: theta, alpha, and sigma2 all vary.

The complete parameter table therefore uses `theta`; `theta,alpha`;
`theta,sigma2`; or `theta,alpha,sigma2`, respectively. A shared CLI parameter
may be fixed only when its regime-specific counterpart is not in that table.
Multi-regime fits use the same generic affine-Gaussian pruning/smoothing engine,
positive parameterization, deterministic multistarts, and boundary reporting.
With exactly one mapped regime, every OUM-family variant is the same statistical
model as ordinary stationary OU, so nwkit instead delegates to the canonical OU fitter.
This gives exactly identical parameters, likelihood, status, and marginals while
avoiding a second generic optimization.
Each free non-root regime parameter must affect a positive-length branch
ancestral to an observation. In OUMVA, a regime appearing only at the root
cannot have both alpha and sigma2 free because only their stationary-variance
ratio is identified.

### LAMBDA, KAPPA, DELTA, EB/ACDC, and BM-DRIFT

The transformed-time Brownian models use one shared implementation and profile
their shape parameter jointly with sigma2:

- LAMBDA applies Pagel's lambda in `[0,1]` to shared/internal covariance while
  preserving each tip variance.
- KAPPA replaces every positive branch length `t` by `t^kappa`
  (`kappa >= 0`); zero-length contractions remain zero, including at kappa 0.
- DELTA transforms ultrametric node depth to
  `height * (depth/height)^delta` (`delta > 0`) and rejects non-ultrametric trees.
- EB and ACDC let diffusion change exponentially with root depth. EB is
  constrained to non-positive early-burst change; ACDC permits either decline
  or acceleration.

For EB/ACDC, a branch from depth `d` to `d+t` has effective Brownian length

```text
exp(eb_rate * d) * expm1(eb_rate * t) / eb_rate
```

with the continuous limit `t` at zero. `--evolution-parameter` and
`--evolution-parameter-bounds` are the common controls; `--eb-rate` and
`--eb-rate-bounds` remain aliases for EB/ACDC. A flat profile, such as an
equal-depth star where shape and sigma2 are confounded, is rejected; bound
optima are reported against the active user-supplied bounds. `--profile-ci-level`
optionally inverts the one-parameter
likelihood-ratio profile and records bounds plus boundary-limited flags in
`--model-out`.

BM-DRIFT uses transition mean `X_parent + drift * t` and Brownian variance
`sigma2 * t`. `--drift` fixes the directional trend. A free trend is profiled
from observed tips at different root depths. The search expands geometrically
until the likelihood optimum is bracketed instead of imposing artificial drift
bounds, which is important with strongly unequal observation errors. Drift is
confounded with the unknown flat-prior root value on contemporaneous tips, so
ultrametric observations require a fixed drift. Both models reduce exactly to
BM when their extension parameter is zero. An exactly linear, error-free fitted
trend is retained as an explicit `sigma2=0`, `singular_zero_boundary` result
rather than a spurious tiny positive diffusion estimate.

A free drift is profiled *after* integrating the flat root; it is not integrated
as a second fixed effect. Consequently, its reported flat-root likelihood and
`residual_df` retain the `n_effective - 1` convention, rather than the
`n_effective - 2` convention of a two-fixed-effect REML analysis. Output records
this treatment explicitly, and intervals condition on the fitted drift.

### MV-BM, MV-OU, and MV-OU-DIAG: correlated continuous traits

`--model MV-BM --state-column trait1,trait2,...` models each branch increment as
`MultivariateNormal(0, Sigma * t)`. Complete exact vectors use the linear-time
contrast/smoothing path. Trait-level missingness or known errors automatically
select a dense observed-coordinate likelihood, with one comma-separated
`--standard-error-column` per trait. Every trait needs enough observations to
identify its mean/covariance; an implementation guard rejects more
than 1,000 dense observed coordinates. This measured cap bounds cubic covariance
factorization and repeated all-node solves; complete error-free MV-BM remains on
the linear-time path and is not subject to the dense-coordinate cap. Explicit
all-zero SE columns are normalized to the same exact-observation path. All three
multivariate models report sample size as the number of distinct observed phylogenetic
positions, independent of trait dimension or dense/fast implementation path.
Even an all-zero SE mapping is validated for complete tip coverage, vector
dimension, finite numeric values, and non-negativity before the fast path is
selected. Dense covariance rank is calculated in normalized trait coordinates,
so changing trait/SE units cannot turn a full-rank fit into a singular one.
Consistent exact observations at one zero-length-contracted position count once;
conflicting values at that position have zero likelihood and are rejected. A
covariance component is also rejected explicitly when the observed trait/branch
overlap leaves it absent from every independent flat-root contrast.

`--model MV-OU` fits a stationary separable process with one shared positive
alpha, a full stationary trait covariance `Sigma`, and one optimum per trait:
`Cov(X_i,X_j) = Sigma * exp(-alpha * patristic_distance(i,j))`. Its equivalent
diffusion covariance is `2*alpha*Sigma`. Alpha may be fixed; Sigma and optima are
estimated by proper stationary-root ML. MV-OU supports the same partial vectors
and diagonal known measurement errors as dense MV-BM. A free alpha is rejected
when every observed pair within each trait combination has only one constant
phylogenetic distance, because alpha can then be absorbed into `Sigma`.

`--model MV-OU-DIAG` replaces the shared attraction rate by a positive rate for
each trait. With `A = diag(alpha_1,...,alpha_d)` and positive-definite diffusion
covariance `D`, its stationary covariance is
`C_ij = D_ij / (alpha_i + alpha_j)`. For nodes `p,q` with LCA depth `s`,
`Cov(X_p,i, X_q,j)` is
`exp[-alpha_i(depth_p-s)-alpha_j(depth_q-s)] * C_ij`. This permits traits to
lose ancestral signal at different speeds while retaining correlated process
noise. Omit alpha options to estimate all `d` rates, use `--alpha FLOAT` to fix
one shared rate, or use `--alpha-by-trait a1,a2,...` in state-column order to
fix distinct rates. These fixed forms are mutually exclusive and cannot be
combined with `--alpha-bounds`. Every freely estimated trait alpha requires at
least two observations at distinct phylogenetic positions. The optimizer
parameterizes `D` by its Cholesky factor, so invalid diffusion covariances are
not searched. With fixed shared alpha, MV-OU-DIAG is exactly MV-OU and
`asrcompare` retains the former as an equivalent row rather than assigning
duplicate IC weight.

Primary output is long-form with one row per node/trait. `--covariance-out`
writes each selected node's conditional covariance and correlation upper
triangle; `--model-out` records stationary Sigma, rank, optimizer diagnostics,
and alpha/theta for OU models. MV-OU-DIAG additionally records every trait alpha
and diffusion-covariance element. Intervals condition on fitted parameters and exclude parameter
and tree uncertainty. Scalar posterior-sampling/PPC/bootstrap outputs are
currently rejected for multivariate fits rather than silently approximated.

### Model comparison and simulation diagnostics

`nwkit asrcompare` is the batch model-selection interface corresponding to
`rootcompare`. It reads the tree and trait table once, resolves the trait type,
fits every applicable requested candidate, and writes one row per candidate:

```sh
nwkit asrcompare -i tree.nwk --trait traits.tsv \
  --state-column body_mass --models all \
  --exclude-models BM-DRIFT -o comparison.tsv \
  --figure-out comparison.pdf --criterion aic
```

`--models all` is the default. It includes every registered model for the
resolved trait type. Models needing unavailable inputs, such as BMS without a
`--regime-map`, are retained with `status=not_applicable`; failures in one
automatically selected model are retained with `status=failed` and do not hide
successful fits. Automatic comparison records HRM as `status=not_fitted` rather
than spending most of the run on a hidden-class fit that cannot participate in
regular-model IC ranking; naming HRM explicitly requests that diagnostic fit.
A named `--models M1,M2,...` request is strict and stops on an inapplicable
candidate or fit exception; numerical diagnostic outcomes such as
non-convergence remain visible as result rows. `--exclude-models` accepts
families or an individual OU root variant.

The command calculates IC values only from finite fits marked as rankable. It
retains non-converged fits, singular fits, HRM label-switching fits, nonregular
COVARION boundaries, zero-variance-component boundaries, and models without a
marginal likelihood, but does not rank them. Other finite boundary fits remain
rankable and are labeled for inspection.
THRESHOLD is reported as `no_likelihood` without running its MCMC because its
posterior diagnostics do not define AIC/AICc/BIC. Structurally identical
parameterizations are retained as `status=equivalent`, but only one
representative contributes IC weight. This includes binary ER/SYM and binary
ARD/F81 under matching transition/root contracts, neutral fixed
lambda/kappa/delta or zero EB/drift reductions to BM, and one-regime reductions
to the corresponding ordinary Mk, BM, BM-DRIFT, or stationary OU model. Binary
GTR is not merged with ARD/F81: its bounded exchangeability/frequency-ratio
parameterization defines a different feasible rate space. EB and ACDC are
merged only when a shared fixed rate or explicit nonpositive bounds make their
complete parameter spaces equal; a coincident point estimate alone does not
merge them. Criterion ties use a small absolute numerical tolerance, dense
ranks, and mark every tied minimum as best; the rule is invariant to a common IC
offset and cannot chain through adjacent near-ties.

Model-specific options are routed only to candidates that consume them and are
rejected if no applicable selected model does. A shared transform value or bound
must be valid for every selected transformed model. Because fixed regime tables
have model-specific exact column schemas, one `--regime-parameters` file cannot
be shared across selected regime models requiring different columns.

Comparisons are partitioned by trait dimensionality, likelihood convention, and
root prior. Equal, empirical, stationary, and Gaussian discrete roots therefore
remain distinct, as do continuous roots. Flat-root integrated BM-family fits,
stationary-root OU fits, and proper fixed/Gaussian-root OU fits never receive a
shared delta or weight. Proper-root variants are also kept separate from one
another. Use `--models 'OU[stationary],OU[fixed],OU[gaussian]'` to request
variants explicitly; fixed and Gaussian candidates require `--root-mean`, and
Gaussian additionally requires positive `--root-variance`. Each criterion needs
at least two finite values in a group before its deltas or weights are emitted;
the selected criterion independently controls ranks and winners.
The optional PDF is a single-page table containing every evaluated candidate.
Shaded section headers visibly separate compatible comparison sets, and models
within each set are sorted from `#1` onward. Rows report fit status, sample and
parameter counts, log likelihood, selected criterion, delta, weight, and notes;
no-likelihood, not-fitted, failed, or input-inapplicable candidates appear in a
final unassigned section.
No bar plot or cross-set ranking is shown. `--criterion aic|aicc|bic` selects
the displayed ranking while all three criteria remain in the TSV. Unicode text
uses an installed font with verified glyph coverage, and rendered-width fitting
prevents wide titles or diagnostics from clipping; if no installed font covers
the requested text, predictable labels are rejected before model fitting and PDF
creation fails explicitly instead of emitting tofu. The error lists the glyphs
actually missing from installed font coverage.

The TSV includes candidate identity and root provenance; compatibility group;
fit, optimizer, and parameter-contract diagnostics; log likelihood and all IC
values/deltas/weights; selected-criterion rank and joint-best flag; and
`equivalent_to` for retained aliases. Count columns are written as integers when
present. `num_comparable_models` is the number with a finite value for the
selected criterion in that row's group. `shared_preparation_seconds` reports
one-time input preparation, including tree and trait parsing plus shared
auxiliary inputs, repeated on each row; `elapsed_seconds` contains only that
candidate's fit/evaluation time. Alias rows have zero fit time.

Multiple discrete columns select MK-MIXTURE; exactly two columns also make the
two Pagel models applicable when both traits are binary. Their different
likelihood kinds form separate comparison sets rather than a false cross-model
ranking. Multiple continuous columns select MV-BM/MV-OU/MV-OU-DIAG candidates.
Numeric categorical columns still need explicit `--trait-type discrete`.

`--compare-models M1,M2,... --model-comparison-out FILE` fits compatible models
alongside a normal `nwkit asr` reconstruction and reports log likelihood,
parameter count, AIC/AICc/BIC, deltas, and weights. It remains available as the
compact backward-compatible interface when every requested model already shares
one likelihood convention.
For binary data it retains ER/SYM and ARD/F81 aliases as explicit
`status=equivalent` rows but fits each exact transition/root contract once and
leaves alias weights empty, preventing duplicate parameterizations from
receiving multiple shares of the IC weight.
Only models sharing one likelihood convention may be compared: discrete ML
models are separate from proper-root OU ML and flat-root integrated BM-family
fits. Discrete BIC/AICc use the number of informative observed tip-character
entries (fully missing or all-state observations contribute no sample). Singular,
non-finite, zero-variance-component, and structurally non-regular COVARION
boundary fits are rejected instead of being ranked under regular-model
information criteria; other boundary statuses remain visible in the comparison
table for interpretation.

For any scalar Gaussian model, `--posterior-samples-out` writes exact joint
all-node conditional draws using forward-filter/backward sampling.
`--posterior-predictive-out` simulates replicated observed tips and reports
mean/variance/range discrepancy summaries with tail probabilities.
`--bootstrap-out` simulates the fitted process, preserves the original missing
pattern and known SEs, refits the same fixed/free parameter specification, and
writes successful and failed replicates. Associated count options default to
1,000 posterior draws, 1,000 predictive replicates, and 100 bootstrap refits;
`--seed` makes all paths reproducible. For a flat root, bootstrap simulation
fixes the fitted posterior root mean; proper-root models draw from their root
prior.

## Output schemas

Both modes retain shared `branch_id`, `parent` (`-1` for root), `node_class`, and
`name` columns. `--target all|intnode|leaf|missing-leaf` accepts comma-separated
classes; `intnode` includes the root. `is_imputed` identifies only missing tips,
not internal nodes or denoised observed tips. An empty selection still produces
a table header.

Continuous summary columns additionally contain:

| Columns | Meaning |
|---|---|
| `trait` | Selected input column; multivariate Gaussian models emit one row for each selected trait |
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

BMS, BMS-DRIFT, and OUM-family fits additionally report regime order, root
regime, source paths, each regime parameter, optimizer counts, and data-scaled
bounds where applicable. Transformed BM and BM-DRIFT report extension
parameters, estimated flags, search details, and the underlying sigma2 fit.
Multivariate models report covariance rank and every stationary/diffusion Sigma
element under collision-safe hex-encoded trait identifiers; the optional
covariance sidecar retains readable trait names.

Discrete `--model-out` records the transition graph, `q_source` (`estimated`,
`fixed:--rate`, or `fixed:PATH`), fit/boundary status, optimizer start/convergence
counts, and every directed Q entry. CUSTOM has an empty fitted-rate-bounds field.
F81/GTR also report equilibrium frequencies; MK-REGIME reports every regime Q;
MK-DESIGN and Pagel report fitted rate classes; Pagel also reports both trait
columns and binary state orders. HRM/COVARION report expanded-state details;
MK-MIXTURE reports category rates, weights and gamma shape; THRESHOLD reports its
identified process, thresholds, MCMC settings and diagnostics.

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
  `conditional_on_covariance`). Multivariate models suffix per-trait properties with the
  UTF-8 hexadecimal trait identifier and includes cross-trait covariances.
- `all`: also tip `asr_observed_value`, `asr_observed_se`, and `asr_is_imputed`.

MK-MIXTURE additionally stacks a `trait` column; THRESHOLD otherwise uses the
normal discrete probability/MAP schema. Optional discrete model metadata records
selected/requested trait types. Every auxiliary output path (model, comparison,
tree, map, covariance, liability, samples, PPC, or bootstrap) must be distinct
and cannot be STDOUT; only the primary TSV may be written to `-`.

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

Trait-specific multivariate OU attraction rates:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column height,mass \
  --model MV-OU-DIAG --alpha-by-trait 0.2,0.8 \
  --model-out mvou-diag-model.tsv -o mvou-diag.tsv
```

Pagel correlated evolution for two binary traits:

```sh
nwkit asrcompare -i tree.nwk --trait traits.tsv \
  --state-column habitat,behavior --trait-type discrete \
  --models PAGEL-INDEPENDENT,PAGEL-DEPENDENT -o pagel-comparison.tsv
```

Custom fitted rate sharing:

```sh
nwkit asr -i tree.nwk --trait stages.tsv --state-column stage \
  --trait-type discrete --model MK-DESIGN --rate-design design.tsv \
  --model-out design-model.tsv -o design-asr.tsv
```

Ordinal threshold reconstruction with explicit MCMC size:

```sh
nwkit asr -i tree.nwk --trait stages.tsv --state-column stage \
  --trait-type discrete --model THRESHOLD --states juvenile,adult,senescent \
  --liability-samples 2000 --liability-out liability.tsv -o threshold.tsv
```

Compare flat-root models and write continuous simulation diagnostics:

```sh
nwkit asr -i tree.nwk --trait traits.tsv --state-column body_mass \
  --model LAMBDA --compare-models BM,LAMBDA,KAPPA,DELTA \
  --model-comparison-out comparison.tsv --profile-ci-level 0.95 \
  --posterior-samples-out draws.tsv --posterior-predictive-out ppc.tsv \
  --seed 42 -o lambda.tsv
```

Compare ordinary discrete likelihood models and retain a diagnostic PDF:

```sh
nwkit asrcompare -i tree.nwk --trait states.tsv --state-column state \
  --trait-type discrete --models ER,SYM,ARD,F81,GTR,COVARION \
  -o model-comparison.tsv --figure-out model-comparison.pdf
```

Lévy/jump processes, tree-uncertainty integration, correlated measurement-error
matrices, and branchwise continuous stochastic trajectories remain outside the
current model set. Joint Gaussian node samples are provided instead of
repurposing discrete transition-count maps.

## Method references

- [ape ancestral character estimation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html): BM and REML conventions.
- [phytools fastAnc](https://search.r-project.org/CRAN/refmans/phytools/html/fastAnc.html): continuous ancestral estimates and uncertainty.
- [Hansen 1997](https://doi.org/10.2307/2411186): OU comparative models and adaptive optima.
- [Butler and King 2004](https://doi.org/10.1086/426002): multi-optimum OU models.
- [Beaulieu et al. 2013](https://doi.org/10.1093/sysbio/syt034): hidden-rate models for discrete traits.
- [Harmon et al. 2010](https://doi.org/10.1111/j.1558-5646.2010.01025.x): early-burst comparative models.
- [Felsenstein 2012](https://doi.org/10.1086/664553): threshold-model comparative inference.

Tests compare the Gaussian passes against an independently assembled full-tree
precision matrix and tip-covariance residual likelihood, plus analytic stars,
zero-edge equivalences, and unit/offset invariance. External-package comparisons
must match rate estimation, root treatment, and interval conventions first.
A recorded `phytools 2.3.0 fastAnc` fixture checks exact-data means and variances
without adding R as a runtime or test dependency. OU tests independently assemble
the stationary patristic covariance and compare ordinary likelihood and all-node
conditional moments. Extension tests additionally compare equal-regime limits to
BM/OU, generic Gaussian smoothing to independent dense conditioning,
multivariate covariance to dense matrix or contrast oracles, shared-alpha
MV-OU-DIAG reduction to MV-OU, Pagel nested rate partitions, threshold constraints
and seeded sampling, and cached versus uncached discrete-mixture likelihoods.
