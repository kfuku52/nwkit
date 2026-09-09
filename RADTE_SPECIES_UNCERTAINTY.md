# Species-age evidence, sensitivity, and uncertainty in RADTE

These workflows analyze one gene family. They do not jointly estimate a species
chronogram across multiple families. Input species trees retain one fixed,
rooted, binary topology and contemporaneous tips.

## Three distinct uses of an age interval

| Input | Role | Effect on inference |
| --- | --- | --- |
| `--species-node-intervals-tsv` | External confidence/credible/HPD/percentile interval | Display and provenance only; no constraint or prior is created |
| `--species-node-bounds-tsv` | Hard allowed age range | Makes specified species ages variable within that range |
| `--species-tree-ensemble` | Joint species-chronogram samples | Conditions a separate gene-family fit on each entire sampled chronogram |

A 95% interval is not a hard bound: replacing one with the other changes its
meaning. A lower/upper pair also does not specify a probability distribution or
the dependence between species node ages. RADTE does not infer these missing
quantities or automatically turn external intervals into priors. For an input
posterior ensemble, account for overlapping sequence evidence before interpreting
a downstream analysis; using the same data again is not a new independent update.

### External interval TSV

```tsv
node	lower	upper	level	kind	source
AB	9.05	10.95	0.95	percentile	Synthetic chronogram sample quantiles
R	18.1	21.9	0.95	percentile	Synthetic chronogram sample quantiles
```

These numbers belong to the synthetic example below, not an empirical analysis.
Use actual tab separators. `node` must identify one unique species-tree node;
`species_event_id` can be supplied instead for stable clade-based matching. If
both are supplied, they must agree. The other five columns are required. `kind`
is `confidence`, `credible`, `hpd`, or `percentile`; `source` identifies the
external analysis. `level` is a probability between zero and one, not a percentage.
Endpoints must be finite, nonnegative and ordered. Duplicate/unknown nodes fail.
Tip intervals must be zero. Omitted nodes have no external interval. Internal
marginal intervals may overlap between parent and child; they are not sampled
independently or treated as simultaneous constraints.

The reference point age still comes from the species Newick branch lengths.
Embedded dating-tool annotations are not automatically interpreted as intervals;
export their endpoints and interval metadata into this TSV explicitly.

```sh
nwkit radte --gene-tree gene.nwk --species-tree species.nwk \
  --species-map-tsv species-map.tsv --reconcile lca --max-age 30 \
  --species-node-intervals-tsv species-intervals.tsv \
  --uncertainty profile --out-prefix fixed --figure-out fixed.pdf
```

The TSV is hashed and protected against output aliases. Its values are retained
as `input_interval_lower`, `input_interval_upper`, `input_interval_level`,
`input_interval_kind`, and `input_interval_source` in `.species.tsv`. The species
output also contains `interval_lower`, `interval_upper`, and `interval_status`
for the current fit, separately from `age_min`/`age_max` hard ranges.

The species panel shows external intervals in gray, hard ranges as dotted lines,
and species result intervals in blue. Gray points use original input ages; blue
nodes use result ages. The legend and diagnostics preserve the interval type and
source. Fixed calibration nodes remain outlined diamonds. Nodes represented only
by constraints are marked as such when variable; they are not claimed to be
independently estimated from a gene speciation event. In a species ensemble,
intervals for unrepresented species events describe the supplied chronogram
samples (`input-chronogram-samples`), not new gene-based estimates. PAML ranges are identified
as prior ranges. Saved `nwkit draw --radte-prefix` reads this evidence from the
verified result bundle without needing the original interval TSV.

## Range-constrained re-estimation

Use a separate `node, age_min, age_max` TSV with
`--species-node-bounds-tsv species-bounds.tsv`. Species nodes not listed remain
fixed. All gene speciation nodes mapped to the same species event share its
estimated age and interval. Species age movement and active calibration limits
should be interpreted together. Retain sufficient information to determine the
absolute time scale; making every positive age freely scalable cannot identify
absolute ages independently of rates.

## Separate input variation from conditional uncertainty

```sh
nwkit radte --gene-tree gene.nwk --species-tree species.nwk \
  --species-map-tsv species-map.tsv --reconcile lca --max-age 30 \
  --species-node-intervals-tsv species-intervals.tsv \
  --species-tree-ensemble species-samples.nwk --uncertainty input-ensemble \
  --ensemble-within-uncertainty profile --out-prefix ensemble
```

Supply retained samples after any external burn-in/thinning; samples receive
equal weight. Each input chronogram is preserved as a joint sample, including dependence
between its node ages. The reference point fit remains the main dated tree.
Ensemble percentiles in the node/species tables summarize point fits across
input samples; conditional intervals do not replace or widen these percentiles.

`--ensemble-within-uncertainty` defaults to `none`. `profile` computes conditional
profile intervals within each chronogram. `bootstrap` uses site resampling with
an alignment or rate simulation without one, with `--bootstrap-replicates` per
chronogram. These nested analyses can multiply the computation substantially.
The reference `--interval-level`, `--starts`, and `--maxiter` apply; sample seeds
are deterministic. A conditional failure does not discard a successful input
point fit. Its error and missing interval are retained separately. The manifest also retains
per-sample conditional diagnostics and profile optimizer attempts under
`input_ensemble.within_fit_diagnostics`.

Additional output tables are always published, empty when not requested:

* `.conditional-intervals.tsv`: original input sample number, shared age ID,
  point age, conditional interval endpoints/status, estimator method, error,
  and conditional bootstrap variance when available.
* `.uncertainty-components.tsv`: valid input/conditional counts, SD and sample
  variance (`ddof=1`) across input-conditioned point fits, mean conditional
  interval width, and mean conditional bootstrap variance when available.

An input summary needs at least 20 valid samples and 90% coverage, including
per-event presence checks. Conditional summaries additionally require at least
20 and 90% valid intervals for that event. Missing values are `NA`, not zero.
Mixed input estimators or conditional methods are flagged and withheld from
the corresponding component summaries. Fixed species ages within one chronogram may have zero conditional
width but nonzero variation between chronograms.

These statistics are reported separately. Profile widths are not converted into
variances; conditional bootstrap variance is not automatically added to the
variance of input-conditioned point estimates. This is not a combined Bayesian
posterior or a claim of total uncertainty coverage.

## Compare the three saved analyses

```sh
nwkit radte-compare --fixed-prefix fixed --bounded-prefix bounded \
  --ensemble-prefix ensemble --species-tree species.nwk --out-prefix comparison
```

This command runs no inference. All runs must use the same reference gene and
species inputs, reconciliation, external interval evidence, estimator/clock
settings, and interval level. The fixed and ensemble reference fits use fixed
species ages; the bounded fit must have at least one variable species age. The
ensemble must contain species chronograms and no gene-tree ensemble, so this
report isolates the species-age analysis. Native hard-bound runs are currently
supported. Mixed ensemble estimators are rejected.

The first PDF page shows three columns (fixed, bounded, ensemble), with the
gene tree above the species tree in each column. All six panels use the same
age scale, with older ages on the left. Intervals are drawn directly at their
nodes; each column reports its interval method and how many internal gene nodes
have intervals available. The second page places a reference tree beside the input-refit SD and
mean conditional interval width; each summary is aligned vertically with its
node. Each summary column shares its scale across gene and species rows, but
these two metrics have distinct meanings and must not be added.

Outputs are `.pdf` (tree-based age comparison and aligned uncertainty summaries),
`.comparison.tsv`, `.uncertainty-components.tsv`, and `.manifest.json` with
source/output hashes. Missing intervals remain visibly unavailable. Source
bundles are validated and protected against overwriting, including by audit
output. A failed render/publication preserves existing comparison outputs.

Run the complete [synthetic example](examples/radte/species-uncertainty/README.md)
to generate all three analyses and the two-page comparison report.
