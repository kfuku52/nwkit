# Conventional PGLS and reconciled speciation contrasts

`nwkit pgls` supports conventional tip-level PGLS as well as reconciled
gene-tree contrasts. The input mode is selected explicitly by its input paths.

## Conventional tip-level PGLS

Supply one rooted species tree and one table containing the response and
predictor values:

```sh
nwkit pgls \
  --tree species_tree.dated.nwk \
  --data species_data.tsv \
  --responses expression \
  --predictors body_size,temperature \
  --evolution-model lambda \
  --outfile ordinary-pgls.tsv
```

The table has one `leaf_name` row per species in the simple case. The model
includes an intercept by default; use `--intercept no` for a regression through
the origin. The following evolutionary residual covariance models are
available:

| `--evolution-model` | Parameter and tip covariance |
|---|---|
| `brownian` | No shape parameter. Covariance is shared root-to-MRCA path length; the default. |
| `lambda` | `lambda` in `[0, 1]`. Multiplies Brownian off-diagonal covariance while preserving each tip variance. |
| `ou` | `alpha > 0`. Fixed-root Ornstein-Uhlenbeck covariance. |
| `kappa` | `kappa >= 0`. Raises every branch length to the power `kappa`. |
| `delta` | `delta > 0`. Raises relative node depths to `delta`; requires an ultrametric tree. |
| `eb` | `rate_change <= 0`. Integrates an exponentially decreasing rate through time; zero is Brownian motion. |
| `acdc` | Finite signed `rate_change`. Negative values decelerate and positive values accelerate the evolutionary rate through time. |
| `independent` | No shape parameter. Identity covariance, equivalent to ordinary least squares. |
| `custom` | No shape parameter. A user-supplied positive-definite named covariance matrix. |

Parameterized models estimate their shape parameter by default. Use the generic
`--evolution-parameter FLOAT` option to fix it. Identity points are
`lambda=1`, `kappa=1`, `delta=1`, and `rate_change=0`; OU approaches Brownian
motion as `alpha` approaches zero. Branch lengths must be positive and finite
unless `--branch-length unit` is selected. Both ML and REML are available, as
are Wald and variance/model-refitted parametric-bootstrap inference.
`independent` and `custom` do not consume tree branch lengths and report their
branch-length mode as `not-applicable`.

Automatic estimation uses finite, scale-aware search regions: `[0, 1]` for
lambda, `[0, 3]` for kappa, `[1e-4, 1e4]` for delta,
`[1e-6 / tree_height, 1e3 / tree_height]` for OU alpha,
`[-50 / tree_height, 0]` for EB rate change, and
`[-50 / tree_height, 50 / tree_height]` for ACDC rate change. Fits at or near a
search boundary are flagged. A scientifically justified fixed value may use
the broader domains in the model table, provided the transformed covariance
remains finite and positive definite.

Biological and technical response replicates use the same
`--biological-id`, `--technical-id`, `--technical-aggregation`, `--batch`, and
`--within-variance` rules as reconciled expression input. Known response means
and standard errors use `--within-variance known-se`. Estimated species-mean
sampling covariance is added to the evolutionary covariance rather than being
treated as extra species.

Predictors independently support the corresponding
`--predictor-biological-id`, `--predictor-technical-id`,
`--predictor-technical-aggregation`, `--predictor-batch`, and
`--predictor-within-variance` options. Known predictor SE and sample-size
columns use `--predictor-standard-error-columns` and
`--predictor-sample-size-columns`. Thus either role or both roles may be
replicated. When only one role is replicated, values for the other role must be
identical within each tip.

Predictor uncertainty is not added to response covariance as though it were
response noise. NWKIT fits a phylogenetic latent predictor and integrates its
posterior covariance into the response likelihood; this is an
errors-in-variables PGLS and avoids the usual attenuation from treating noisy
predictor means as exact. `--sampling-covariance-out` and
`--tip-summary-out` audit the response calculation, while
`--predictor-sampling-covariance-out` and `--predictor-tip-summary-out` audit
the predictor calculation.

For example, a table with known uncertainty for both roles can be fitted with:

```sh
nwkit pgls \
  --tree species_tree.dated.nwk \
  --data species_means.tsv \
  --responses expression \
  --predictors body_size \
  --within-variance known-se \
  --standard-error-columns expression_se \
  --predictor-within-variance known-se \
  --predictor-standard-error-columns body_size_se \
  --outfile ordinary-pgls.tsv
```

### Custom evolutionary covariance

Select `--evolution-model custom` and supply `--evolution-covariance` as a wide
TSV. Its first column and row names are `leaf_name`; its remaining column names
must be exactly the same tree-tip names:

```tsv
leaf_name	A	B	C
A	2.0	0.7	0.2
B	0.7	1.8	0.4
C	0.2	0.4	1.5
```

Rows and columns may be in any order, but their names must match exactly. The
matrix must be numeric, finite, symmetric, and positive definite. NWKIT
estimates one evolutionary rate multiplying this fixed covariance. This
interface also permits externally constructed multi-rate or other covariance
models without hard-coding their parameterization into NWKIT. Tree branch
lengths are not used by `custom`, although the rooted tree still defines and
validates the species set and regression row order.

### Evolutionary-model comparison

Conventional PGLS can compare models on the same response data and design:

```sh
nwkit pgls \
  --tree species_tree.dated.nwk \
  --data species_data.tsv \
  --responses expression \
  --predictors body_size,temperature \
  --evolution-model lambda \
  --compare-evolution-models brownian,lambda,ou,kappa,delta,eb,acdc,independent \
  --model-comparison-out evolutionary-models.tsv \
  --outfile ordinary-pgls.tsv
```

Comparison fits always use ML, even when the selected coefficient model uses
REML. Every parameterized comparison model estimates its own shape parameter;
`--evolution-parameter` fixes only the selected main model. The comparison TSV
reports log likelihood, AIC, AICc, BIC, parameter/convergence diagnostics,
delta AIC and AIC weights, and delta AICc and AICc weights. The likelihood
parameter count includes all regression coefficients, the evolutionary rate,
and a model shape parameter when present. An individual AICc is empty when
`n_species <= n_likelihood_parameters + 1`; delta AICc and AICc weights are
left empty for the whole response comparison if any candidate lacks AICc. AIC
weights remain available.
Include `custom` in the list only when `--evolution-covariance` is also
supplied.

## Reconciled speciation-contrast PGLS

NWKIT separates tree reconciliation, contrast calculation, and statistical
inference. This keeps every mapping decision inspectable and prevents repeated
observations of one species-tree split from being treated as independent.
The exact reconciliation rule, PIC recursion, covariance propagation, and a
hand-calculated two-species/two-paralog example are given in
[`RECONCILED_SPECIATION_CONTRAST_MATH.md`](RECONCILED_SPECIATION_CONTRAST_MATH.md).

## Required trees and ordering of operations

Use the reconciler-annotated tree to identify events, but use a time-calibrated
version of the same rooted topology to calculate model-aware contrasts.
GeneRax branch lengths are sequence-substitution lengths and should not be used
as divergence times. In a GeneGalleon workflow, run GeneRax, date the complete
tree with RADTE, and only then prune the annotated and dated trees to the same
set of tips. Dating before pruning avoids unconstrained old duplication ages.

When pruning an NHX tree with NWKIT, preserve its properties:

```sh
nwkit prune \
  --infile gene_tree.nhx \
  --pattern '^(gene1|gene2|gene3)$' \
  --invert-match yes \
  --preserve-properties yes \
  --outfile gene_tree.pruned.nhx
```

If pruning collapses a duplication or transfer node with only one retained
child, NWKIT records `NWKIT_COLLAPSED_EVENT_BOUNDARY=Y` on the surviving
lineage. `reconcile` exposes this as `collapsed_event_boundary` and starts a new
lineage cluster there, while correctly omitting a contrast for the removed
binary event.

## Workflow

The recommended interface runs reconciliation, both PIC calculations, and the
final model with one command. Expression is the response and the species trait
is the predictor:

```sh
nwkit pgls \
  --gene-tree gene_tree.pruned.dated.nwk \
  --reconciliation-tree gene_tree.pruned.nhx \
  --species-tree species_tree.dated.nwk \
  --expression gene_expression_replicates.tsv \
  --species-traits species_traits.tsv \
  --responses expression \
  --predictors body_size \
  --biological-id sample_id \
  --batch sequencing_batch \
  --predictor-biological-id organism_id \
  --predictor-batch phenotyping_batch \
  --event-source nhx \
  --species-map-tsv gene_to_species.tsv \
  --tree-id OG000001 \
  --gene-evolution-model kappa \
  --gene-evolution-parameter auto \
  --species-evolution-model delta \
  --species-evolution-parameter auto \
  --out-prefix OG000001
```

`--gene-tree` supplies the dated tree whose branch lengths define expression
PICs. `--reconciliation-tree` is optional and defaults to `--gene-tree`; use it
when event annotations are retained on a separate GeneRax/NHX tree. The two
gene trees must have identical rooted topology and tip names, although their
branch lengths and annotations may differ. Intermediate tables are constructed
in memory.

`--gene-evolution-model` controls expression contrasts on the gene tree, and
`--species-evolution-model` independently controls predictor contrasts on the
species tree. Both default to `brownian`. In raw end-to-end mode, omitting a
parameter or selecting `auto` estimates one gene-tree parameter per response
and one species-tree parameter per predictor. A numeric value fixes the same
parameter across the selected traits.

Gene parameters maximize the final reconciled ML/REML likelihood. Every trial
rebuilds ancestral estimates, raw contrasts, evolutionary variances, and the
replicate-derived sampling covariance while reusing the topology-dependent
reconciliation. Predictor parameters are instead estimated by marginal ML from
the species-tip trait and species tree. This prevents the response association
from selecting the predictor transformation. The final results report fixed or
estimated status, convergence, boundary, and predictor marginal-likelihood
diagnostics. The contrast tables retain the exact selected transforms.

When parametric-bootstrap inference is requested with an automatically
estimated gene parameter, every bootstrap replicate simulates in tip space,
re-estimates that parameter, rebuilds the response contrasts and sampling
covariance, and refits the variance components. Species-trait transforms remain
conditional on the observed predictors. This can be substantially slower than
a fixed-parameter bootstrap.

`--out-prefix OG000001` writes an inspectable bundle:

| Path | Contents |
|---|---|
| `OG000001.reconciliation.tsv` | Reconciliation decisions and event mappings |
| `OG000001.gene-contrasts.tsv` | Expression contrasts with reconciliation context |
| `OG000001.species-contrasts.tsv` | Species-trait predictor contrasts |
| `OG000001.response-sampling-covariance.tsv` | Full replicate-derived response covariance, when applicable |
| `OG000001.response-tip-summary.tsv` | Biological/technical replicate audit, when applicable |
| `OG000001.predictor-sampling-covariance.tsv` | Full predictor covariance after the species-tree PIC transform, when applicable |
| `OG000001.predictor-tip-summary.tsv` | Predictor biological/technical replicate audit, when applicable |
| `OG000001.random-effects.tsv` | Conditional event and lineage effects; header-only when none are fitted |
| `OG000001.pgls.tsv` | Final coefficient and variance-component results |

NWKIT stages every applicable table before committing the bundle. A per-prefix
lock serializes concurrent writers. If ordinary staging or commit fails, newly
installed files are removed and pre-existing bundle members are restored. An
audited run records the planned paths even when no new bundle is committed.

Without `--out-prefix`, only `--outfile` (or standard output) and an explicitly
requested `--random-effects-out` are written. This is useful for scripted runs
that do not need retained intermediates.

### Inspecting or reusing individual stages

`reconcile` and `contrast` remain available as low-level commands. Use them to
inspect event assignments, reuse one reconciliation across several expression
datasets, or start from an externally prepared reconciliation. Precomputed
gene and species contrasts can still be fitted by supplying `pgls --infile`
and `--predictor-contrasts`; raw and precomputed inputs cannot be mixed.

First map the pruned, rooted GeneRax tree onto the rooted species tree. Supply a
stable family identifier whenever multiple trees will later be combined:

```sh
nwkit reconcile \
  --infile gene_tree.pruned.nhx \
  --species-tree species_tree.nwk \
  --event-source nhx \
  --species-map-tsv gene_to_species.tsv \
  --tree-id OG000001 \
  --outfile reconciliation.tsv
```

`--event-source lca` applies LCA reconciliation. `nhx` reads GeneRax-style
`S`, `D`, and `H` properties, including transfer annotations of the form
`H=Y@source@destination`. A valid `S` value is the primary species-tree
placement, which is essential when transfer makes descendant-species LCA
incorrect. Missing `S` falls back to descendant LCA; an unknown, ambiguous, or
species-map-conflicting `S` is reported as unmapped instead of being silently
replaced. Internal `S` values therefore require the corresponding GeneRax
species-node labels to be retained in `--species-tree`. `species-overlap` is
available as an explicit heuristic.

Calculate organismal-trait contrasts on a time-calibrated version of the same
species-tree topology used above. The reconciliation step uses topology, but
the PIC calculation must not reuse arbitrary GeneRax or substitution-scale
branch lengths:

```sh
nwkit contrast \
  --infile species_tree.dated.nwk \
  --trait species_traits.tsv \
  --columns body_size \
  --evolution-model delta \
  --evolution-parameter 1.2 \
  --tree-id species_tree \
  --outfile species_contrasts.tsv
```

Then calculate expression contrasts on the matching pruned, time-calibrated
gene tree—not on the raw GeneRax substitution-length tree. The recommended
input has one row per biological observation and gene-tree tip:

```sh
nwkit contrast \
  --infile gene_tree.pruned.dated.nwk \
  --trait gene_expression_replicates.tsv \
  --columns expression \
  --evolution-model kappa \
  --evolution-parameter 0.8 \
  --biological-id sample_id \
  --batch sequencing_batch \
  --within-variance pooled \
  --reconciliation reconciliation.tsv \
  --event-type speciation \
  --tree-id OG000001 \
  --sampling-covariance-out expression_sampling_covariance.tsv \
  --tip-summary-out expression_tip_summary.tsv \
  --outfile gene_contrasts.tsv
```

The low-level `contrast` command supports `brownian`, `lambda`, `ou`, `kappa`,
`delta`, `eb`, `acdc`, and `independent`. A parameter is required for every
parameterized model. NWKIT calculates PICs on the exactly equivalent
transformed tree and records the model, parameter, and branch-length mode in
each output row. `delta` and contrast-based `ou` require an ultrametric tree.
An arbitrary `custom` tip covariance is not accepted here because it need not
factor into independent local tree contrasts; use it only in conventional
tip-level PGLS.

If technical replicates are present, identify them with `--technical-id` and
request `--technical-aggregation mean`. This averages already transformed
continuous expression values within a biological observation; technical rows
never increase `n_biological`. The default `error` policy prevents accidental
pseudoreplication. A categorical `--batch` is fitted as a fixed observation-
level effect. Leaf and batch effects must be separable; confounded designs are
rejected.

`--within-variance pooled` estimates one residual biological variance shared
across gene-tree tips. `leaf` estimates one variance per tip and requires at
least two biological observations for every tip; it cannot be combined with
batch adjustment. If means and standard errors have already been estimated,
use `--within-variance known-se` and columns named `<trait>_se`, or name them
with `--standard-error-columns`. Optional `--sample-size-columns` add sample
sizes to the audit without changing the supplied standard errors.

The contrast output contains each contrast's sampling variance, while the
sidecar contains the full covariance. NWKIT computes the PIC coefficient
matrix `L` and propagates the tip-mean covariance `D` as `M = L D L'`.
Off-diagonal elements are retained because tip sampling error and batch
adjustment can correlate contrasts. The tip summary records biological and
technical counts, fitted means, within-tip SDs, standard errors, variance
method, and batch status.

Fit precomputed contrasts through the low-level PGLS input mode:

```sh
nwkit pgls \
  --infile gene_contrasts.tsv \
  --predictor-contrasts species_contrasts.tsv \
  --response-sampling-covariance expression_sampling_covariance.tsv \
  --predictor-sampling-covariance body_size_sampling_covariance.tsv \
  --responses expression \
  --predictors body_size \
  --model hierarchical \
  --random-effects-out random_effects.tsv \
  --outfile pgls.tsv
```

`pgls` performs and validates the join
`gene_contrasts.species_event_id = species_contrasts.branch_clade_id`. It also
requires exact agreement between `species_event_taxa` and `descendant_taxa`,
and between the numerator and denominator event IDs. These SHA-256 clade
identifiers depend only on sorted descendant tip labels, so they remain stable
when Newick child order changes. Numeric `species_branch_id` is retained for
inspection within one exact serialization but is not a cross-file join key.
Both input tables must be direct `nwkit contrast` outputs: selected rows require
`tree_id`, `evolution_model`, `evolution_parameter_name`,
`evolution_parameter`, and `branch_length_mode`. NWKIT rejects mixed or
internally inconsistent transforms within a response or predictor contrast
group, and all selected predictor rows must come from one non-empty species-tree
ID. The coefficient table records that predictor-tree ID and the exact response
and predictor transforms on every row.
When `--predictor-sampling-covariance` is supplied, every selected predictor
contrast pair must be present and the species-contrast table must retain its
positive `contrast_variance` column. The legacy clustered estimator rejects
both response and predictor covariance sidecars; use `hierarchical` or
`replicate-reml`.

## PGLS model and pseudoreplication control

The model is a regression through the origin on signed reconciled contrasts:

```text
expression contrast = species-trait contrasts * fixed slopes + error
```

The default `--model hierarchical` fits the Gaussian covariance

```text
V = sigma^2 G + M + tau_event^2 Z_event Z_event'
    + tau_lineage^2 sum_j diag(x_j) Z_lineage Z_lineage' diag(x_j)
```

`G` is the diagonal transformed-tree evolutionary variance of each selected
gene-tree contrast and `M` is the fixed sampling covariance propagated from
biological replicates.
`Z_event` gives paralog contrasts at the same `species_event_id` a shared
response deviation. The final term is a partially pooled lineage-specific
random slope for each predictor; it uses `lineage_clade_id` and does not add an
uninterpretable lineage intercept. Evolutionary, event, and lineage variances
are estimated by REML by default. `--reml no` selects ML.

With predictor replicates, the observed species contrast is modeled as

```text
x_hat_j = x*_j + eta_j,       eta_j ~ N(0, Mx_j)
x*_j ~ N(0, rho_j Kx_j)
```

Here `Kx_j` contains the predictor contrast variances from the transformed
species tree and `rho_j` is estimated by marginal likelihood. Gaussian
conditioning gives posterior mean `m_j` and covariance `S_j` for each latent
predictor. If `R` maps one value per species event onto repeated gene-tree rows,
the response likelihood is

```text
y | x_hat ~ N(sum_j beta_j R m_j,
              V + sum_j beta_j^2 R S_j R')
```

The added covariance depends on the fitted slope and is optimized jointly with
the response variance components. Repeated paralog rows therefore reuse the
same latent species-event value and covariance; they do not become independent
predictor measurements. Output distinguishes raw predictor sampling variance
from remaining posterior latent-predictor variance and reports the fitted
predictor evolutionary rate and its diagnostics.
For conventional tip-level PGLS, the same calculation is performed in tip
space with `R = I` and an estimated predictor ancestral mean.

`--event-random-effect auto` and `--lineage-random-slope auto` include only
identifiable covariance components. `yes` requires the component and fails if
the data cannot identify it; `no` removes it. `--model replicate-reml` retains
`sigma^2 G + M` but omits both random effects. `--model legacy` keeps the older
event-cluster HC1 estimator only as a sensitivity analysis.

Rows sharing a `species_event_id` have equal total working weight by default,
so a species split represented by ten paralogs does not dominate a split
represented once. `--event-weighting observation` is an explicit sensitivity
analysis. In every model, reported residual degrees of freedom are
`n_species_events - num_parameters`, never the number of gene-tree rows.

Wald inference is the default. `--inference parametric-bootstrap` simulates
from the fitted covariance, refits all variance components, reports bootstrap
standard errors and percentile intervals, and computes a centered empirical
p-value. `--bootstrap-replicates` and `--seed` make it reproducible. Models
with fewer than 20 unique species events and fits at a variance boundary are
flagged. Conditional event deviations and lineage slopes can be written with
`--random-effects-out`.

The command fits one model per non-empty `tree_id` and response trait;
accidental pooling of unnamed gene families is rejected. Multiple predictor
traits are supported, and each can have its own replicate design or known SE.
The current likelihood assumes zero measurement-error covariance between
different predictor traits and between a predictor and response. If one assay
induces material cross-trait measurement covariance, fit or propagate a joint
measurement model externally rather than assuming those errors are independent.

## Coverage policy and gene loss

An eligible speciation event must have both immediate gene-tree children map
unambiguously to different daughter clades of the mapped species-tree node.
Eligibility permits gene loss, but `coverage_status` records whether both full
species daughter clades are represented. Numerator and denominator observed
taxa, full clade taxa, counts, and coverage fractions are all emitted.

For safety, `--event-type speciation` uses `--eligible-only auto` and
`--speciation-coverage complete`, retaining only completely sampled daughter
clades. A sensitivity analysis that admits explicitly reported partial coverage
can be requested with:

```sh
nwkit contrast \
  --infile gene_tree.pruned.dated.nwk \
  --trait gene_expression.tsv \
  --columns expression \
  --reconciliation reconciliation.tsv \
  --event-type speciation \
  --speciation-coverage any
```

Whether the complete-coverage analysis, the partial-coverage sensitivity
analysis, or a matched-sampling trait estimator is the primary estimand is a
scientific decision. NWKIT reports the information needed to make that choice
explicit rather than treating partial and complete events as equivalent.

## Repeated species events and downstream inference

An ancient duplication can produce several paralog lineages that cross the
same species split. NWKIT intentionally writes one row per gene-tree speciation
node. It does not average those expression contrasts and does not claim they
are independent replicates of the organismal-trait contrast.

Use these identifiers downstream:

- `tree_id`: caller-supplied gene-family identifier;
- `gene_clade_id`: stable gene-tree node identifier;
- `species_event_id`: stable species-tree event identifier and cluster;
- `lineage_clade_id`: stable lineage segment descended from the most recent
  observed or collapsed non-speciation event.

Inference must account for repeated `species_event_id` values. `nwkit pgls`
does so with equal total event weighting and, when identifiable, a shared
species-event effect. The effective information for an organismal predictor
therefore comes from unique species events, not from the number of paralogs
carrying a repeated predictor contrast. `lineage_clade_id` supplies the
partially pooled lineage-specific slopes.

`--eligible-only auto` means yes for `--event-type speciation` and no for
`all`, `duplication`, `transfer`, or `unresolved`. Thus `--event-type all`
actually reports all internal event categories, and selecting duplication or
transfer no longer silently produces an empty table.

## Validation and limitations

- Both trees must be rooted and strictly bifurcating. Polytomies are not
  silently resolved.
- Standardized contrasts require positive finite branch lengths and finite
  intermediate calculations. Use `--branch-length unit` only when unit lengths
  are scientifically intended.
- Missing observations within a replicate table are allowed, but every tree
  tip needs an estimable mean. A completely unobserved tip is rejected; prune
  it explicitly and regenerate the reconciliation for the resulting topology.
- Reconciliation rows are matched to the contrast tree by `gene_clade_id`, not
  traversal-order branch numbers. Enum values, coverage state, orientations,
  tree identity, and exact clade coverage are validated before calculation.
- LCA reconciliation maps duplication/speciation events under a
  duplication-loss model but does not enumerate loss events. It does not model
  incomplete lineage sorting or horizontal transfer; use reconciler NHX
  annotations when those processes matter.
- Conventional PGLS can estimate supported evolutionary shape parameters and
  compare models by ML. Reconciled raw PGLS estimates gene parameters within
  each response model and predictor parameters from species-trait marginal ML.
  Precomputed mode and the low-level `contrast` command require fixed
  parameters because they do not have both trees and original tip values.
- `delta` and reconciled/contrast `ou` require ultrametric trees. Conventional
  fixed-root OU can also use a non-ultrametric tree. Arbitrary `custom`
  covariance is available only for conventional tip-level PGLS because it need
  not admit independent local contrasts.
- Count likelihoods, cross-trait measurement-error covariance, and joint
  multivariate responses are not fitted.
- `pgls` deliberately fits separate models by `tree_id` rather than silently
  pooling gene families. Branch lengths must be comparable within each model.
- Small-sample inference remains uncertain. Results with fewer than 20 unique
  species events are flagged and should be interpreted cautiously.

The PIC calculation follows Felsenstein's independent-contrasts recursion.
Applying PICs to internal speciation and duplication annotations follows Dunn
et al. (2018). Begum and Robinson-Rechavi (2021) show why gene-tree sampling,
event ages, branch scaling, pruning order, and contrast dependence require
special care.

References:

- Felsenstein J. 1985. Phylogenies and the comparative method. *The American
  Naturalist* 125:1–15. https://doi.org/10.1086/284325
- Dunn CW et al. 2018. Pairwise comparisons across species are problematic
  when analyzing functional genomic data. *PNAS* 115:E409–E417.
  https://doi.org/10.1073/pnas.1707515115
- Begum T, Robinson-Rechavi M. 2021. Special care is needed in applying
  phylogenetic comparative methods to gene trees with speciation and
  duplication nodes. *Molecular Biology and Evolution* 38:1614–1626.
  https://doi.org/10.1093/molbev/msaa288
- Ives AR, Midford PE, Garland T Jr. 2007. Within-species variation and
  measurement error in phylogenetic comparative methods. *Systematic Biology*
  56:252–270. https://doi.org/10.1080/10635150701313830
