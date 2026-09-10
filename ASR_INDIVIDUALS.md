# Joint evolutionary and within-species covariance

`nwkit asr --model MV-BM --within-species-covariance full` jointly estimates
Brownian evolutionary covariance **Sigma** and a common within-species
individual covariance **W**. It reconstructs ancestral states and latent species
means, and predicts missing traits of identified individuals. Unequal numbers
of individuals and partially observed trait vectors are supported.

```bash
nwkit asr \
  --infile examples/individual_asr/tree.nwk \
  --trait examples/individual_asr/individuals.tsv \
  --state-column x,y --model MV-BM \
  --within-species-covariance full --covariance-method REML \
  --target all --outfile species-and-ancestors.tsv \
  --individual-out individuals.tsv --model-out model.tsv \
  --covariance-out node-covariances.tsv
```

Use `--within-species-covariance diagonal` to constrain individual trait
correlations to zero. Sigma remains a full covariance matrix. The default
estimation method is REML; `--covariance-method ML` profiles the root mean.
The rooted tree and its branch-length units determine Sigma's rate units.
Zero-length edges are allowed when the observation design remains identifiable.

## Individual input

With `--within-species-covariance`, `--trait` uses this **long** schema:

```tsv
leaf_name	individual_id	trait	value
A	1	x	2.8
A	1	y	5.9
A	2	x	3.2
A	2	y	NA
B	1	x	5.1
B	1	y	7.4
```

`--state-column x,y` names the selected traits in their model/output order;
at least two are required. Each `(leaf_name, individual_id, trait)` key must
be unique, including rows with missing values. IDs are strings scoped within
species, so `A/1` and `B/1` are different individuals. Leading zeros and literal
IDs such as `NA` are preserved. Empty IDs are invalid.

Missing cells may be explicit missing-value tokens or omitted trait rows.
The shared `--missing-values` policy applies to **values**, not IDs or trait
names. All traits in the input must be listed in `--state-column`. Standard
`--unmatched warn|error|ignore` behavior applies; table-only species are excluded
under `warn`/`ignore`, and tree-only species can be reconstructed. A listed
individual with all traits missing contributes no likelihood information but
receives predictions. A wholly absent individual has no output row.

Every row represents a biological individual, with traits measured on that
same individual linked by its ID. Technical replicates must not be assigned
new individual IDs. Known instrument SEs and estimated W are different
observation models: `--standard-error-column`, `--measurement-covariance`, and
`--replicate-observations` cannot be combined with this mode. The existing
wide-table ASR and known-SE replicate modes retain their existing behavior.

## Model and likelihood

For species `s`, individual `i`, and trait vector `y`:

```text
root mean beta → latent species value z_s → individual value y_si
                 BM covariance Sigma       common covariance W
```

`y_si = z_s + e_si`, where `e_si ~ N(0, W)` independently between individuals
and independently of evolution. For observed scalar coordinates `a,b`,

```text
V[a,b] = C[species_a,species_b] * Sigma[trait_a,trait_b]
       + same_individual(a,b) * W[trait_a,trait_b]
```

`C` is the Brownian shared-path covariance measured from the input root.
The design matrix `X` has one intercept per trait. Only observed coordinates
enter V and X; this is the marginal observed-data likelihood, without filling
missing values before fitting or treating sample means as exact observations.
The missingness pattern is conditioned on; the model does not account for
informative selection of individuals or traits.

For `N` observed coordinates, `p` traits, GLS root estimate `beta_hat`, and
residual `r = y-X beta_hat`, the reported ML log likelihood is

```text
-1/2 * [N log(2 pi) + log|V| + r' V^-1 r].
```

REML adds `log|X' V^-1 X|` inside the brackets and uses `(N-p) log(2 pi)`.
This is the flat-root integrated convention; it does not add an `|X'X|`
normalization. All constants, including individual replicate contributions,
are retained. Traits are centered at their observed minimum and scaled by
their observed range, so normalization does not depend on row order. This
normalization is reversed in both the matrices and reported likelihood. REML likelihoods from different observed
coordinates or trait sets should not be ranked against each other.

## Outputs and uncertainty

- `--outfile` is one row per selected node/trait, retaining original tree
  branch IDs. `estimand` is `ancestral_mean` or `latent_species_mean`.
  `sample_mean`, `sample_mean_se_estimate`, and `num_individuals` describe input
  data at tips. The SE is `sqrt(W[t,t]/n_t)` under the fitted model. It is an
  estimate, not a known measurement SE. `mean`, `variance`, `sd`, and the
  confidence limits describe the **latent species/ancestral value**. Even an
  observed species generally has nonzero posterior variance. `is_imputed`
  indicates that the species has no observed individual for that trait.
- `--individual-out` is one row per listed individual and every model trait.
  `estimand=individual_value`; `is_imputed` identifies missing coordinates.
  Observed individual values are conditioned on exactly and have zero
  conditional variance. Missing coordinates use W and other measured traits
  of that same individual, as well as the phylogeny and other individuals.
- `--covariance-out` contains the joint conditional trait covariance at each
  selected node. This is distinct from the fitted rate matrix Sigma and W.
- `--model-out` records the likelihood method, counts, convergence and
  identifiability diagnostics, fitted root means, and all ordered entries of
  Sigma and W. Matrix fields are `sigma_<hex-trait>_to_<hex-trait>` and
  `within_<hex-trait>_to_<hex-trait>`, using the existing ASR UTF-8 hex encoding.
- `--tree-out` supports the existing multivariate NHX annotations. Tip
  observed-value/SE annotations under `--tree-annotation all` describe the
  sample means and their fitted SEs; reconstructed means describe latent
  species values. The TSV outputs explicitly name these distinctions.

Intervals include uncertainty in the unknown root mean, conditional on the
fitted Sigma, W, and input tree. They **exclude covariance-parameter and tree
uncertainty** and are labeled `conditional_on_fitted_Sigma_W` in the model
output. They are not bootstrap intervals. Bootstrap, held-out diagnostics,
posterior samples, figures, model comparison, and tree ensembles are not yet
implemented for this observation mode; incompatible options are rejected.

All requested files are staged together. Validation or a handled export error
preserves existing files; outputs cannot overwrite an input or alias each other.

## Identifiability and numerical limits

Each trait needs at least two within-species replicate degrees of freedom
(`sum_s max(n_st-1,0) >= 2`). Full W additionally requires at least two paired
within-species degrees of freedom for **every** trait pair. With less paired
information, choose diagonal W or collect paired data. These are conservative
minimum data requirements, not guarantees of precise estimates.

The covariance design is checked after projecting out trait root means.
Confounded Sigma/W components, unobserved traits, zero within-species variation,
and singular complete-data within-species contrasts are rejected. Unresolved
numerical scales and nonconverged fits produce errors, not ordinary estimates.
Matrix symmetry and positive semidefiniteness are checked in correlation units;
a very long unobserved branch cannot conceal an invalid observed covariance.
Unit restoration preserves representable subnormal variances and rejects
positive variances that would round to zero.
Missingness can remove information even when the original table has many rows.

Fitting uses three deterministic starts, positive-definite Cholesky factors,
analytic gradients, SLSQP optimization, and a projected-gradient convergence
check. Converged-start counts use that same check, not just the optimizer
success flag. Conditional predictions use Cholesky solves; observed individual
coordinates have exactly zero conditional covariance. On normalized
trait/time scales, log diagonal factors are bounded to `[-18,8]`, and
lower-triangular off-diagonals to `[-1000,1000]`. A factor at a bound or a scaled
covariance eigenvalue ratio below `1e-6` yields `fit_status=boundary_covariance`.
A smallest eigenvalue below `1e-8` times the largest eigenvalue of normalized
Sigma+W is also flagged, including components that approach zero as a whole.
Such matrices are numerical approximations near the positive-semidefinite
boundary, not evidence of a well-resolved full-rank covariance. Their
conditional intervals do not account for boundary/parameter uncertainty.
Three starts do not guarantee a global maximum of this nonconvex likelihood.

This first implementation uses dense observed-coordinate likelihoods. It
rejects fits whose estimated matrix workspace exceeds 512 MiB, and node
covariance reconstruction exceeding 256 MiB. Runtime is cubic in the number
of observed coordinates per likelihood evaluation. Species-specific W,
technical-replicate layers, shared-lambda and OU extensions are deferred.

## Verification

[The runnable example](examples/individual_asr/README.md) has eight species,
three or four individuals per species, two correlated traits, and one missing
coordinate. [reference.R](examples/individual_asr/reference.R) uses Rphylopars
0.3.10 with matched BM, W structure, ML/REML and individual grouping.
Tests also use independently assembled dense covariance/GLS and conditional
prediction calculations, finite-difference gradient checks, unit and row-order
changes, malformed/degenerate designs, three-trait partial data, file rollback,
and repeated simulation recovery of both Sigma and W.
