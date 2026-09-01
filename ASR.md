# Ancestral trait reconstruction

`nwkit asr` reconstructs a single trait at every node and imputes missing tip
values. It supports discrete Mk models (ER, SYM, ARD) and continuous Brownian
motion (BM). The tree must be rooted under the
[shared rootedness contract](CLI_TSV_CONVENTIONS.md#rootedness-and-root-polytomies).
Root and internal polytomies, non-ultrametric trees, and finite non-negative
non-root branch lengths are supported. Root stem length is not part of either
model. No rerooting, arbitrary binary resolution, or automatic trait
transformation is performed.

## Input and automatic type selection

The TSV requires `leaf_name` and the column selected by `--state-column`.
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
| `--model` | ER (default), SYM, ARD | BM (default) |
| `--root-prior` | equal (default), empirical | flat (default) |
| `--output` | probabilities (default), map | summary (default) |
| `--tree-annotation` | map (default), state, probability, all | summary (default), mean, all |
| Rate controls | `--rate`, `--rate-bounds` | `--sigma2` |
| State-specific controls | `--states`, `--ambiguous-separator` | Not applicable |
| Observation uncertainty | Ambiguous states | `--standard-error-column` |
| Interval coverage | Not applicable | `--ci-level` (default 0.95) |
| Stochastic mapping | `--stochastic-map-out`, `--n-sim`, `--threads`, `--seed` | Not applicable |

Model/output/prior defaults are selected after trait detection. The continuous
root prior describes uncertainty about a numeric root value; it is independent
of `--input-rooted`, which describes the interpretation of the tree.

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

## Output schemas

Both modes retain shared `branch_id`, `parent` (`-1` for root), `node_class`, and
`name` columns. `--target all|intnode|leaf|missing-leaf` accepts comma-separated
classes; `intnode` includes the root. `is_imputed` identifies only missing tips,
not internal nodes or denoised observed tips. An empty selection still produces
a table header.

Continuous summary columns additionally contain:

| Columns | Meaning |
|---|---|
| `trait` | Selected input column |
| `observed_value`, `observed_se` | Original observation and its SE; empty for internal/missing nodes, SE zero for exact observations |
| `is_imputed` | Whether this is an unobserved tip |
| `mean`, `variance`, `sd` | Conditional latent-trait moments in original units |
| `ci_lower`, `ci_upper`, `ci_level` | Equal-tail conditional interval and its coverage |

`--model-out` records `trait_type`, `trait_type_requested`, the trait/model,
`root_prior`, `sigma2`, `sigma2_estimated`, `estimation_method` (`REML` or `fixed`),
`restricted_log_likelihood`, `num_observed`, `num_effective_observations`,
`residual_df`, `fit_status`, SE-column selection, and the interval's conditioning
contract. Rate and tree uncertainty flags are false.

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
annotations include `asr_trait_type=continuous` and:

- `mean`: `asr_mean` only.
- `summary` (default): also `asr_variance`, `asr_sd`, `asr_ci_lower`,
  `asr_ci_upper`, `asr_ci_level`, and `asr_interval_kind=conditional_on_sigma2`.
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

Numeric categories must explicitly select discrete inference:

```sh
nwkit asr -i tree.nwk --trait categories.tsv --state-column state \
  --trait-type discrete --states 0,1,2 --model ER -o discrete.tsv
```

BM is currently the only continuous model. OU, correlated multiple traits,
parameter/tree uncertainty integration, and continuous conditional trajectory
sampling are not implemented. Discrete transition-count stochastic maps are
not repurposed as continuous trajectories.

## Method references

- [ape ancestral character estimation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html): BM and REML conventions.
- [phytools fastAnc](https://search.r-project.org/CRAN/refmans/phytools/html/fastAnc.html): continuous ancestral estimates and uncertainty.

Tests compare the Gaussian passes against an independently assembled full-tree
precision matrix and tip-covariance residual likelihood, plus analytic stars,
zero-edge equivalences, and unit/offset invariance. External-package comparisons
must match rate estimation, root treatment, and interval conventions first.
A recorded `phytools 2.3.0 fastAnc` fixture checks exact-data means and variances
without adding R as a runtime or test dependency.
