# Branch-specific scalar Gaussian processes

The Python API in `nwkit.branch_gaussian` builds a single process with fixed
BM, OU and prescribed Gaussian jump parameters assigned to individual branches.
It returns the existing `GaussianTreeProcess`, usable directly for likelihood,
ancestral-state conditioning, covariance and simulation. Parameters and branch
assignments are supplied by the caller. The `nwkit asr --model BRANCH-GAUSSIAN` CLI reads
models from TSV files. Adding `--branch-fit` estimates selected diffusion
parameters at a fixed regime assignment; automatic model search is not performed.
The process-building API itself continues to accept fully specified parameters.

## Model and units

Each incoming branch obeys `child = a * parent + b + independent error`, with
error variance `q`. For branch length `t`:

| Diffusion | a | b | q |
| --- | --- | --- | --- |
| `BrownianBranch(variance_rate=sigma2)` | 1 | 0 | sigma2 × t |
| `OUBranch(alpha, variance_rate=sigma2, optimum=theta)` | exp(−alpha × t) | (1−a) × theta | sigma2 × (1−a²) / (2 × alpha) |
| No diffusion (`None`) | 1 | 0 | 0 |

OU with `alpha=0` is exactly BM; small positive alpha uses a stable finite-time
formula. Alpha has inverse tree-time units, sigma2 has trait²/tree-time units,
and theta has trait units. Finite, nonnegative lengths, alpha and variance
parameters are required. Zero variance is allowed, including deterministic
branches. Unrepresentable positive diffusion variances raise an error.

`BranchGaussianModel(diffusion, jump=GaussianJump(mean=m, variance=v))`
adds an independent event **after diffusion, at the branch end**: `b += m`
and `q += v`. Jump variance is per event and is not multiplied by branch length.
A jump on a zero-length branch still occurs. `diffusion=None` gives a pure jump
with identity propagation between parent and child. A model must specify at
least a diffusion or a jump. This differs from the existing CLI `JUMP-BM`, which
integrates latent Poisson jump counts and is not marginally Gaussian.

## Branch assignments and root

Use IDs from `nwkit.util.assign_branch_ids(tree)`: level-order traversal with
root ID 0, matching the existing branch-ID convention. Supply **every non-root
ID exactly once**. Root ID 0 is excluded because the root prior is a separate,
required `GaussianRootPrior`. Integer IDs are required; strings and booleans
are rejected. Pass a mapping, not a list of pairs (which can contain duplicate
IDs). The supplied tree must be a root node; detach an attached subtree first.
Recompute assignments after topology or child-order changes.
Repeated node names do not affect ID-based assignment.

The root can be fixed (zero variance), Gaussian (positive variance), or flat
(`variance=None`). The `stationary` label is also accepted with an explicitly
supplied mean and positive variance; heterogeneous branches have no automatically
chosen stationary root. Root variance and jump variances are absolute trait²
values, independent of each branch's diffusion rate. Unlike the global
`variance_scale` in `build_evolutionary_process`, this builder does not rescale
the root. All innovations and jumps are mutually independent and independent
of the proper root draw.

## Example

Run the complete [mixed-model example](examples/branch_gaussian/mixed_process.py):

```sh
python examples/branch_gaussian/mixed_process.py
```

A minimal assignment is:

```python
from nwkit.branch_gaussian import (
    BranchGaussianModel, BrownianBranch, GaussianJump, OUBranch,
    build_branch_gaussian_process,
)
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids

ids = assign_branch_ids(tree)
background = BranchGaussianModel(BrownianBranch(variance_rate=0.8))
models = {branch_id: background for node, branch_id in ids.items() if not node.is_root}
selected_id = next(branch_id for node, branch_id in ids.items() if node.name == "A")
models[selected_id] = BranchGaussianModel(
    OUBranch(alpha=0.6, variance_rate=1.3, optimum=0.9),
    GaussianJump(mean=0.2, variance=0.15),
)
process = build_branch_gaussian_process(
    tree, models, root=GaussianRootPrior("gaussian", mean=0.3, variance=0.9),
)
```

Shared regime parameters can be expressed as a dictionary of named
`BranchGaussianModel` instances, then expanded using a branch-ID-to-regime map;
the runnable example shows this pattern. Model objects are immutable.

Pass `process` to `gaussian_tree_likelihood`, `condition_gaussian_tree`,
`simulate_gaussian_process`, or `sample_gaussian_posterior` in
`nwkit.gaussian_inference`. Likelihood and conditioning accept leaf-name values
and `standard_errors` (standard deviations, not variances); missing values use
the existing inference API's `None` convention. Conditional intervals describe
state uncertainty given the supplied parameters, excluding parameter-selection
uncertainty. Prior simulation returns latent node values, without observation
error; add independent measurement noise separately if needed. Flat-root prior
simulation requires explicit `root_values`, and unconditional covariance requires
a proper root. Exact deterministic constraints follow the existing inference
engine's singular-Gaussian likelihood conventions.

Numerical checks also compare 153 extreme parameter combinations with a
750-digit Decimal oracle, including weak selection whose unscaled alpha × time
product underflows.

Tests compare mixed-process covariance, likelihood and all-node conditioning
with an independently assembled structural-equation oracle, simulation moments,
uniform BM/OU reductions, the alpha→0 limit, and zero-duration end jumps.


## CLI and TSV inputs

`nwkit asr --model BRANCH-GAUSSIAN` provides `--output summary` (default),
`likelihood`, and `prior-samples`. See the [runnable CLI example](examples/branch_gaussian/README.md).
All modes require a rooted tree and one of these assignment forms:

- `--branch-models FILE`: direct `branch_id`, `model` and parameter columns.
- `--branch-regimes FILE --regime-models FILE`: a `branch_id, regime` assignment
  TSV and a `regime, model` parameter TSV. Every used regime must have exactly
  one definition, and unused definitions are rejected.

The comma-separated column lists above describe tab-separated files. Required
parameters depend on `model`:

| model | Required numeric columns | Other allowed parameters |
| --- | --- | --- |
| `BM` | `sigma2` | Optional `jump_mean` and `jump_variance` together |
| `OU` | `sigma2`, `alpha`, `theta` | Optional `jump_mean` and `jump_variance` together |
| `JUMP` | `jump_mean`, `jump_variance` | No diffusion parameters |

`model` is case-sensitive. `sigma2` is the variance rate, not a standard
deviation. In a mixed table, unused cells must be blank; zero is a specified
value, not a missing marker. Numeric fields must be finite, and rate, alpha and
variance fields must be nonnegative. Extra/duplicate columns, duplicate IDs,
empty/ragged rows, missing parameters and nonblank unused parameters are errors.
A UTF-8 BOM is accepted. Every non-root ID must occur exactly once. Root 0 is
excluded from **both** assignment forms; this differs from the existing ASR
`--regime-map`, which also assigns the root.

### Root and observation options

`--root-prior` is mandatory. `fixed` requires `--root-mean` and has zero
variance. `gaussian` and `stationary` require both `--root-mean` and a positive
`--root-variance`. The stationary label never estimates or selects a root from
the heterogeneous branches. `flat` accepts neither root parameter, and prior
simulation additionally requires `--prior-root-value`. That starting value is only
valid for flat-root simulation.

ASR and likelihood require `--trait FILE --state-column NAME`. An optional
`--standard-error-column NAME` supplies nonnegative measurement standard
deviations. Observed tips require a nonmissing SE when that column is selected.
Missing trait values and tree tips absent from the trait table are unobserved;
at least one observation is required. The shared `--missing-values` and
`--unmatched` policies apply to trait tables, not model parameter tables.
`--ci-level` defaults to 0.95 for node summary and prior distribution plots.

Prior simulation accepts `--prior-samples` (default 1000, 1–10000) and an
optional nonnegative `--seed` for reproducibility. It does not condition on
observations or add measurement error, and rejects trait and SE inputs.
`--state-column` optionally labels the prior plot (default `trait`). Output is
limited to two million node values, including roots. For joint posterior draws
use `--posterior-samples-out FILE --posterior-samples N` with summary output.

### Outputs and reproducibility

`-o/--outfile` defaults to stdout. ASR uses the scalar ASR summary columns:
`branch_id`, `parent`, `node_class`, `name`, `trait`, `observed_value`,
`observed_se`, `is_imputed`, `mean`, `variance`, `sd`, `ci_lower`, `ci_upper`,
`ci_level`. `--target` selects nodes (default all); figures always show every
node. `--tree-out` exports the standard ASR NHX annotations. Likelihood emits
one fixed-assignment model summary row, including `model`, `trait`, `root_prior`,
`log_likelihood`, `likelihood_rank`, `num_observed`, and
`num_observed_positions`. Prior sampling emits `simulation`, `branch_id`, `parent`,
`node_class`, `name`, `value`; each one-based simulation includes every latent
node value. Parents use original branch IDs, with −1 for the root.

`--branch-models-out FILE` exports normalized direct assignments accepted by
`--branch-models`. The shared `--model-out FILE` remains a TSV of model and
likelihood metadata. `num_parameters_estimated` counts the free groups (zero
without `--branch-fit`); `parameter_estimation` contains the group estimates,
bounds and diagnostics as JSON in one TSV cell.
`--process-out FILE` writes JSON schema version 2 containing
the root, topology/lengths, normalized branch parameters, optional regime names,
affine transitions, used observations/SEs, operation, seed, draw count, interval
level and inference summary. The `estimation` object records estimated groups
and optimizer/identifiability diagnostics, or is null for fixed parameters.
JSON is a run record; it is not a separate model
input format. Auxiliary outputs require file paths, not stdout. Normalized model TSVs
use 17 significant digits for float round trips.

One input may read from stdin (`-`), including a model/assignment table. Output
paths must be distinct from all inputs and from one another, including aliases.
Related file outputs are staged together; handled write failures restore the
previous set. Stdout is written after file publication. `--audit FILE` records
all model/assignment inputs and exported outputs alongside the shared provenance
record.


### Figures and diagnostics

`--figure-out FILE.png|FILE.svg|FILE.pdf` uses the common ASR figure: phylogeny,
node means and intervals, branch model colors, and OU optima. A star marks each
prescribed end jump; an open diamond marks an imputed tip. Add
`--figure-tip-heatmap yes` and `--figure-trait-tip-labels yes` to locate observed
and missing tips. These are scalar models with fixed assignments. With `--branch-fit`, node
intervals, optima and histories use the fitted diffusion parameters.

`--figure-simulations N --figure-simulation-mode conditional` adds sampled
histories conditioned on all observed tips and their standard errors.
`unconditional` (the default) draws new histories; a proper root follows its
specified distribution, while a flat root starts at the inferred root mean.
`--figure-simulation-steps` controls the finite diffusion grid. Each original
branch's jump is applied exactly once after diffusion and is drawn at the same
time as its pre-jump endpoint, including on a zero-length branch. Connecting
ASR node means is not a reconstruction of a branch trajectory.

With `--output prior-samples`, the figure instead shows analytic prior node
means/intervals and optional unconditional histories. It has no observations or
imputation markers. A flat root starts at `--prior-root-value`. The node samples
in the TSV and finite-grid histories are separate seeded realizations of the
same process; adding a grid changes the random draws. They need not have the
same numerical endpoints.

Summary output also supports posterior predictive checks, leave-one-tip/clade
cross-validation, independent measurement errors and replicate observations.
All intervals and draws condition on the supplied/fitted parameters, branch assignments
and tree. Bootstrap uncertainty, model comparison, correlated
measurement-error input and tree ensembles are rejected for this fixed branch
assignment interface. `asrcompare --models all` records it as not applicable;
use the explicit ASR command to evaluate a supplied assignment.

See the [eight-tip plotting example](examples/branch_gaussian/plot/README.md)
for rendered posterior and prior figures and their exact inputs.


## Estimation at a fixed assignment

Add `--branch-fit fit.tsv` to `asr --model BRANCH-GAUSSIAN` with
`--branch-regimes` and `--regime-models`. Each input regime model still contains
numeric values: these become starting values for estimated parameters and stay
fixed for all other parameters. Branch IDs, BM/OU model types, root parameters,
and Gaussian jump means/variances remain fixed. This is an ASR fitting option;
it does not invoke `shift` or search for regime locations.

The fit TSV has exactly these five columns:

```tsv
regime	parameter	group	lower	upper
Background	sigma2	diffusion	0.01	3
Adapted	sigma2	diffusion	0.01	3
Adapted	alpha	pull	0	8
Adapted	theta	optimum	-5	5
```

This example estimates **three** parameters: the shared BM/OU diffusion rate,
OU alpha, and OU optimum. Distinct group names permit different values;
repeating a group ties the same parameter across those regimes. Every tied row
must have identical initial values and bounds. Parameters omitted from this TSV
stay fixed. The regime names must exactly match the assignment tables.
Direct `--branch-models` input cannot be combined with fitting; use named regimes
to make parameter sharing explicit.

Only `sigma2` (BM or OU), `alpha` (OU), and `theta` (OU) may be estimated. Each
regime/parameter appears at most once. All bounds must be finite and increasing;
starting values must lie within them. Estimated sigma2 needs a strictly positive
lower bound; alpha permits zero. Fixed parameters, including fixed sigma2=0,
retain the ordinary process rules. Bounds are part of the requested statistical
model, not confidence intervals. An estimate at a bound is reported as
`fit_status=boundary`, with the affected groups and boundary sides; it must not
be interpreted as an unconstrained optimum.

Proper roots (`fixed`, `gaussian`, or explicitly specified `stationary`) use
maximum likelihood with their supplied root parameters unchanged. A flat root
uses the existing root-integrated likelihood, labeled `flat_root_integrated`;
for all-BM shared-rate fitting without known errors this gives the usual
residual-contrast rate estimate. **Alpha cannot be estimated with a flat root**:
changing alpha changes the root loading and therefore the scale of its improper
integral. Fix alpha or supply a proper root. A `stationary` root is a supplied
reference distribution, not an equilibrium recalculated from fitted branch
parameters. Proper-root and flat-root likelihoods are not directly comparable.

The implementation uses bounded deterministic multistart optimization and
independently checks projected finite-difference gradients before accepting
convergence. Positive rate/alpha bounds use logarithmic coordinates; alpha
bounds starting at zero use log(1 + alpha × tree-time scale). This is a numerical
search, with no global-optimum guarantee. JSON records start counts, failures,
termination messages, checked gradients, initial likelihood and group estimates.
The fit centers the trait origin before constructing OU transitions and keeps
objective resolution independent of the initial likelihood. Finite input
precision still limits values whose differences are tiny relative to their
absolute magnitudes.

Estimation requires positive-definite observed (or flat-root contrast)
covariance, more residual observation dimensions than free groups, at most
20 free groups and 512 observed tips. Numerical local identifiability is checked
using whitened mean/covariance derivatives on the actually observed tips;
flat-root checks remove the root-mean direction. Rank-deficient specifications
or fitted optima are rejected. Full local rank is not a guarantee of precise
estimation or global identifiability. In particular, fitting a separate variance
to every branch is not generally justified by one trait per tip.

Standard ASR missing-value handling, known SEs, replicate aggregation, tree
annotations, conditional posterior samples and figures remain available.
Cross-validation refits re-estimate the same free groups using only training
observations, while retaining their bounds and fixed layout. Posterior-predictive
draws, node intervals and histories condition on the fitted parameters;
they exclude parameter-estimation and assignment-selection uncertainty.

`--branch-models-out` writes the fitted numeric assignment, so it can be reused
without `--branch-fit` to reproduce the ASR or generate prior samples. Fitting
requires observations and is rejected with `--output prior-samples`. All fit
inputs participate in stdin, input-alias protection and audit recording; all
outputs retain transactional publication.

See the [fitted eight-tip example](examples/branch_gaussian/fit/README.md) for
complete inputs, fitted values, plots and a reproduction script.
