# Phylogenetic principal components

`nwkit pca` summarizes continuous traits using a phylogenetic covariance model.
It returns species scores and can also write rotation coefficients, evolutionary
loadings, eigenvalues, the reusable transformation, ancestral PC estimates and
a phylomorphospace figure.

```sh
nwkit pca -i examples/pca/tree.nwk --trait examples/pca/traits.tsv \
  --columns size,shape,performance --mode corr \
  -o scores.tsv --loadings-out loadings.tsv --eigenvalues-out eigenvalues.tsv \
  --model-out pca.json --ancestral-out ancestors.tsv --figure-out morphospace.png
```

![Synthetic phylogenetic PCA example](examples/pca/morphospace.png)

The figure connects observed tips to conditional ancestral means in PC space.
Its shaded ellipses condition on fitted parameters and axes; lines are not
reconstructed or sampled paths along branches. The right panel shows
phylogenetic trait–PC correlations, not rotation coefficients.

## Input and model choices

- `--trait` follows the shared tip-keyed TSV contract (`leaf_name`, UTF-8,
  configurable `--missing-values`, and `--unmatched warn|error|ignore`).
- `--columns` selects at least two distinct continuous numeric traits. At least
  three complete tips are required. Constant columns are rejected explicitly.
- `--missing error` is the default. `--missing drop` removes incomplete tips
  jointly across all selected columns, reporting excluded tips in model JSON.
  No automatic imputation or pairwise covariance estimation is performed.
- `--model BM` (default) uses shared root-to-MRCA distances. `--model LAMBDA`
  estimates one common lambda by multivariate ML on `[0,1]`; `--lambda-value`
  fixes it instead. A fixed value requires `--model LAMBDA`.
- `--mode cov` (default) uses evolutionary covariance; trait units affect the
  result. `--mode corr` standardizes by the fitted evolutionary standard
  deviations. It does not standardize by ordinary sample standard deviations.

Rooted polytomies and non-ultrametric trees are supported under the shared
`--input-rooted` policy. All original non-root branches must have finite,
nonnegative lengths. The retained tree starts at the complete tips' MRCA,
excluding any common stem. Tip covariance must be positive definite under the
selected model; singular matrices are not repaired with jitter.

This command treats complete tip observations as exact. Sampling SEs, unknown
within-species covariance, multiple imputation and uncertainty in estimated
PC axes are not implemented here. Use an appropriate joint model when those
sources of uncertainty determine the scientific conclusion.

## Calculation and interpretation

For an `n × p` trait matrix `X` and tip covariance `C`, estimate the root mean
by generalized least squares:

```
center = (1' C^-1 X) / (1' C^-1 1)
R = X - center
S = R' C^-1 R / (n - 1)
```

For BM, `C` is the phylogenetic covariance. For lambda:

```
C(lambda) = diag(C_BM) + lambda * (C_BM - diag(C_BM))
```

For covariance PCA, set `scale=1` and diagonalize `S`. For correlation PCA,
set `scale=sqrt(diag(S))` and diagonalize the corresponding correlation matrix.
The implementation uses an SVD of whitened, centered traits rather than
forming an explicit inverse. It calculates

```
scores = ((X - center) / scale) @ rotation
```

Scores remain associated with species; they are not whitened contrast rows.
The `loading` field is the evolutionary correlation between a trait and a PC.
`rotation` gives the coefficients used by the score transformation. Explained
variance ratios describe **evolutionary** covariance/correlation, not the
ordinary variance of the observed tip scores.

Components are sorted by decreasing eigenvalue. To choose a reproducible sign,
the largest absolute rotation coefficient in each component is made positive
(coefficients tied within relative tolerance `1e-12` follow input trait order). This does not make signs biologically intrinsic.
Repeated or nearly repeated eigenvalues permit rotations within the corresponding
subspace; model JSON flags adjacent eigenvalues within `1e-8` of the largest
eigenvalue's scale.

Only numerically positive-rank components are output. Rank is determined from
singular values with tolerance `machine_epsilon * max(n,p) * largest_singular`.
BM or fixed-lambda PCA can handle dependent traits or more traits than residual
degrees of freedom, with `status=rank_deficient`. Free-lambda estimation requires
full trait rank and more tips than traits because its unrestricted multivariate
Gaussian likelihood otherwise has no finite covariance MLE. A star covariance
cannot identify lambda; use BM or an explicit fixed value instead.

For lambda estimation, the matrix-normal ML likelihood profiles both the GLS
mean and trait covariance, using the **ML divisor `n`**. The PCA decomposition
itself uses `n-1`. Likelihoods are calculated in original trait/tree units.
A rank-deficient fixed-model PCA has no full-dimensional log likelihood and
writes JSON `null` for that field. Optimization uses NWKIT's shared bounded
scalar search, including both endpoints.

The [phytools PCA documentation](https://search.r-project.org/CRAN/refmans/phytools/html/phyl.pca.html)
and Revell (2009), *Size-correction and principal components for interspecific
comparative studies*, describe the method. Comparisons with phytools 2.3.0
confirm the decomposition up to component signs. That implementation's
`likMlambda` plugs its `n-1` covariance into the likelihood; NWKIT profiles the
`n` covariance. For identical `C`, `n` and `p`, NWKIT's log likelihood is higher
by `np/2 * log(n/(n-1)) - p/2`, a constant in lambda. Phytools may also search
lambda above 1; comparisons must use the same bounds or a fixed lambda.

## Outputs

Primary `-o/--outfile` output is a wide TSV with `leaf_name, PC1, PC2, ...`,
containing only complete retained tips in sorted name order. It defaults to
stdout. Either the tree or the trait table may use stdin, but not both.

| Option | Output |
|---|---|
| `--loadings-out` | Long TSV: `trait, component, rotation, loading, center, scale` |
| `--eigenvalues-out` | `component, eigenvalue, explained_variance_ratio, cumulative_variance_ratio` |
| `--model-out` | JSON with fitted transformation, lambda, likelihood, rank/status, used/excluded tips and original root branch ID |
| `--ancestral-out` | Long TSV: original `branch_id, parent_branch_id, name, node_class, component, mean, variance, ci_level, ci_lower, ci_upper` |
| `--figure-out` | PDF, SVG or PNG with phylomorphospace and loadings panels |

All auxiliary outputs require file paths. Destinations must be distinct and
must not alias either input. Tables, JSON and figures are staged together:
handled fitting, rendering or installation failures preserve previous outputs.
`--audit` uses the shared provenance machinery and records the auxiliary files.

The model JSON's `trait_names` fixes column order; `center`, `scale` and
`rotation` are sufficient to project a new complete numeric matrix using the
formula above. New observations do not refit the model automatically.

## Ancestors and figures

Ancestral projection reuses NWKIT's continuous BM ASR with each PC's eigenvalue
as a fixed diffusion rate, and the fitted lambda branch transformation when
applicable. The root mean has a flat prior. Variances include uncertainty in
that mean, conditional on the tree, fitted lambda, axes and diffusion rates.
They omit uncertainty in PCA estimation, model choice, missing-data selection
and the input phylogeny.

Original input-tree level-order branch IDs are preserved before excluding tips.
The retained MRCA is labelled `node_class=root` and has parent `-1`, even if its
original branch ID is not zero. Nodes ancestral only to omitted tips are absent;
unary nodes on retained paths remain. Observed tips are included with exact
scores and zero variance. The table contains marginal normal intervals;
`--ci-level` defaults to 0.95.

`--figure-components 1,2` selects two distinct available components. A figure
requires rank at least two; a rank-one analysis can still write tables and JSON.
`--figure-tip-labels yes|no` controls tip names (default `yes`). The figure's
ellipses are joint two-PC Gaussian regions at `--ci-level`, using the chi-square
quantile with two degrees of freedom. Their definition differs from two
independent marginal intervals.

The current numerical implementation uses dense tip covariance and SVD. It is
intended for moderate-dimensional comparative analyses; it does not claim
linear memory use or support for arbitrary extreme floating-point ranges.
