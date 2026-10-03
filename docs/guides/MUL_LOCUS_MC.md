# Experimental Locus DL + ILS Comparison

`mul-reconcile --score-model locus-mc` compares no polyploidy plus ordinary
duplication/loss against one disomic direct-parent allopolyploid event.
It uses a generative locus-history model, not a product of D+L and MSC scores.
This is an explicit small-family research mode: finite-grid estimation and
Monte Carlo integration, not continuous MLE or a production WGD detector.
Default D+L and [conditional MSC](MUL_MSC.md) retain their contracts.

## Model and Observations

The three layers are population/subgenome history, locus birth/death history,
and gene genealogy. The locus process is continuous-time linear duplication
and loss, with rates per locus per generation on every unfolded population
branch, including an explicit ancestral DL stem. One already-existing locus
is present at the beginning of that stem. The induced species-root copy count
is random, including extinction; it is not forced to one at the species root.
The original ancestral coalescent population is infinite above the stem.

SSD marks a mother locus and a new daughter. All extant sampled lineages of
that daughter must coalesce by its birth time. Count-distribution dynamic
programming normalizes this bound before backward genealogy sampling;
there is no rejection of rare daughter histories or an ILS penalty added
to parsimony. This follows the bounded-coalescent construction in
[DLCoal's three-tree model](https://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html).
This implementation's direct-parent MUL extension is experimental, not a
claim to reproduce all DLCoal/AlloppNET features or their reconstruction tools.

Each present copy represents a **distinct locus**, with one sampled genome
per locus, not an allelic replicate. Different copies can be ordinary
paralogs or homoeologs; their hidden histories are generated and integrated.
The allopolyploid model unfolds H1's first-parent stem and a supplied H2
attachment, with fixed species ages and shared diploid/subgenome Ne.
WGD attachment itself is not an SSD daughter-founder constraint.
No subgenome exchange, copy-number hemiplasy, HGT, population-rate
heterogeneity, autopolyploidy or multiple polyploid events are modeled.

The full gene genealogy is sampled before independent per-species detection.
Undetected extant loci remain in daughter conditioning; they are not pruned
from the coalescent beforehand. Loss and non-detection are distinct processes.
The observation is the species-colored rooted binary topology and its copy
counts, not arbitrary gene IDs or a best homoeolog assignment. Probabilities
sum over all compatible hidden histories; they do not condition separately
on each family's exact observed copy counts.

Selection is explicitly `2 <= observed tips <= max_observed_tips`. Both
hypotheses use this same event and the same complete input family set.
Rejecting simulated families outside that event is its normalization, not
history truncation. An input outside it fails; no input family is silently
filtered. Other curation rules (for example requiring particular species,
single-copy outgroups or support thresholds) are **not** included in this
selection model and can invalidate the comparison.

## Configuration and CLI

Species trees must be rooted, strictly binary, ultrametric in **generations**,
with an omitted/zero root stem. Gene trees must be rooted, strictly binary,
with unique tip IDs mapping to known species. Gene lengths are not scored.
H1 is one supplied non-root polyploid clade, not inferred here. H2 can be
restricted using the existing selectors. A configuration file supplies all
scientific assumptions and simulation budgets; no biological rates/bounds
are silently invented.

Example model JSON for the pilot chronogram `((A:1,X:1):1,B:2);`:

```json
{
  "schema": "nwkit-mul-locus-mc-model-v1",
  "copy_role": "distinct-loci",
  "root_locus_count": 1,
  "species_time_unit": "generations",
  "ancestral_stem": 0.5,
  "detection": {"A": 0.9, "X": 0.9, "B": 0.9},
  "max_observed_tips": 4,
  "parameter_grid": [
    {"duplication": 0.03, "loss": 0.03, "ne": 0.3, "hybridization_age": 0.3},
    {"duplication": 0.08, "loss": 0.03, "ne": 1.0, "hybridization_age": 0.7}
  ],
  "samples": 10000,
  "integration": "hybrid-rb",
  "seed": 20261026,
  "confidence": 0.99,
  "max_attempts": 1000000,
  "max_locus_nodes": 10000,
  "max_coalescent_states": 100000
}
```

All shown keys except `integration` are required; unknown/duplicate keys and
duplicate grid points fail. Numeric scientific fields require finite JSON
numbers, not booleans or numeric strings; the grid must be a nonempty array
of parameter objects. Combined rates and the diploid `2*Ne` scale must also
be finite. Omitted `integration` retains `selected-histogram`.
Detection must cover every input species with known probabilities in `(0,1]`.
Rates/stem are finite nonnegative, Ne is finite positive, and alternative
attachment ages must lie inside the fixed H1 stem. Incompatible H2/grid
pairs are explicit exclusions. At least one supported alternative must remain.
Redundant null age points are deduplicated because attachment age does not
exist under no polyploidy; null parameter records show that field as missing.

```bash
nwkit mul-reconcile -i genes.nwk --species-tree species.nwk \
  --species-regex '.*_([^_]+)$' --score-model locus-mc --h1 X --h2 'A B' \
  --locus-model model.json --locus-bootstrap 99 \
  --locus-null-calibration grid-supremum \
  -o scores.tsv --report genes.tsv --check-out checks.tsv \
  --model-out results.json --locus-calibration-out null-search.tsv
```

The example numbers are synthetic, not recommendations for biological data.
MSC-only options fail in this mode; locus-only options fail in D+L/MSC.
`--tree-out` is unsupported: a finite-grid MC winner is not an identified
dated biological point estimate. Only the primary table may use stdout.

## Integration Error and Search

With `selected-histogram`, every candidate/grid point receives `samples`
IID selected families. The alternative stratified raw-draw budget is
described below.
Continuous duplication/loss times, extinction, missing copies and gene
topology histories are integrated by sampling, not by discretizing event
times or enumerating only convenient reconciliations. Histogram probability
estimates are normalized over the selected species-colored observation
space. Log estimates can have sampling bias; their error bounds and raw
hit/sample counts must accompany interpretation.

Clopper-Pearson bounds apply to each histogram count. Bonferroni allocation
uses every bank and the complete finite observation-universe upper bound
`sum((2*n-3)!! * num_species**n, n=2..max_observed_tips)`. Thus the interval
allocation does not depend on selecting only patterns seen in the data or
bootstrap. Likelihood bounds propagate through multiplicities and the
finite-grid maxima. They bound **bank Monte Carlo error**, not gene-tree
error, biological uncertainty, model adequacy or global identifiability.

A zero histogram count has a zero point estimate and a positive upper bound,
not a pseudocount or proof of impossible biology. If no null/alternative grid
point has a finite score for the complete dataset, analysis fails explicitly;
increase simulation resources in a separately documented run. An insufficient
budget is not repaired by family dropping or smoothing. Node/state/selection
attempt caps abort the whole analysis; no overflowing history is discarded.
MC-overlapping grid points are recorded rather than called confident winners.
Serial and multiprocess execution use fixed candidate/grid seed streams.

Alternatively, explicit `integration: "ancestral-stratified"` spends `samples`
total **raw** draws, allocated evenly over the positive-prior ancestral
no-event and first-birth strata. At stem length T and combined rate r,
their original prior weights are `exp(-r*T)` and
`duplication/r*(1-exp(-r*T))`. A first loss on the one-locus stem cannot
produce a selected family. Conditional first-birth times follow the
truncated exponential; subsequent histories still use the complete DL
process. The observation estimate is the sum of prior-weighted raw pattern
frequencies divided by the sum of prior-weighted selection frequencies,
not an equal average of already-selected strata. This preserves the model
while sampling otherwise rare ancestral SSD histories more often.

For stratification, confidence allocation additionally covers both strata's
category and selection counts (a conservative `3*universe` bound per bank).
Numerator/denominator bounds propagate through their weighted ratio. Raw
stratum weights, samples, selection counts and histograms are saved. This
can reduce missing support but does not guarantee sufficient precision;
zero-hit and resource failures remain explicit. The same budget and method
apply to every hypothesis, and bootstrap draws use the original unstratified
generative model rather than treating the proposal as biological truth.
An underflowing positive stratum prior fails explicitly. Confidence precision
is checked as the observation-universe bound grows, before sampling or
constructing an unrepresentably large bound; this does not reduce the declared
family-size range. Work limits are checked before potentially dense
coalescent-transition allocations.

### Conditional Integration

Explicit `detection-rb` uses the same ancestral raw-draw stratification but
integrates all independent detection outcomes for each full genealogy.
`hybrid-rb` additionally integrates all gene genealogies for histories with
at most four **hidden extant loci**; larger histories use a genealogy draw
and exact detection integration. That rule is fixed before sampling, not a
retry after a resource failure. Neither method truncates hidden histories.
Observed counts above the declared limit are exactly aggregated into an
overflow state, preserving the selection denominator and total probability.

These modes require at least two draws per positive-prior stratum. Weighted
pattern/selection contributions save sparse Welford moments, including
implicit zero draws. Bounds invert `n*kl(mean,p) <= log(2/alpha)`, the
[two-sided Chernoff bound for IID values in [0,1]](https://arxiv.org/pdf/2205.07880),
not Clopper-Pearson on fractional hits. The same conservative stratum/universe
union budget and prior-weighted ratio propagation apply. Point likelihoods
are finite-MC ratios, not claimed unbiased. Raw variance reduction does not
guarantee precise normalized likelihoods or an identified parent.

The report retains the `hits` column but leaves it blank for weighted modes;
use `probability_estimate` and the saved `selected`/`patterns` moments, never
interpret fractional contributions as binomial successes. Each saved bank
identifies its `interval_method` (`chernoff-kl`), population and stratum
weights, enabling independent score replay. Older histogram outputs and
their seed streams are unchanged. Serial/parallel execution is byte-identical.

## Event Calibration

The statistic is `2*(max alternative log estimate - max null log estimate)`;
it may be negative because these models are not treated as nested. There is
no chi-square or AIC shortcut. Optional `--locus-bootstrap B` generates fresh
family sets under the numerically selected null grid point. Every replicate
reselects null parameters and searches all supplied alternative parents/grid
points, using the same frozen independent integration banks as the observed
analysis. The default `--locus-null-calibration plug-in` is a **plug-in
finite-grid parametric bootstrap**, not
uniform composite-null coverage or a posterior event probability.

The point P-value is `(1+exceedances)/(B+1)`, with ties included. A conservative
MC-overlap P-value interval counts definitely/possibly exceeding contrast
intervals. It conditions on the **numerical generating null grid point**:
it does not envelope uncertainty about which null point generates the
replicates, even if several null score intervals overlap. This interval is
about score-bank error and does not include bootstrap finite-sample
uncertainty; its resolution is limited by B. No automatic
production decision or scientific confidence interval is returned.
The interval does not by itself establish an event-test error guarantee:
bank coverage has its declared error probability, and plug-in composite-null
and model-misspecification errors remain outside it. Calibration P-values
and MC-overlap bounds are also printed to stderr when requested, even
without optional output files.

### Finite Null Grid Supremum

Explicit `--locus-null-calibration grid-supremum` instead generates B
replicates from **every distinct supplied null grid point**, including points
that are not numerical winners or whose MC score intervals do not overlap
the best null. Every replicate repeats the same full candidate/grid search.
No point is dropped using data-dependent screening. The P-value and each
MC-overlap endpoint are separately maximized over those generating points,
not averaged and not pooled into one mixture null. B is the number per
point; total work is B times the number of null points. This option requires
positive B. Omission preserves the plug-in behavior and its seed stream.

For fixed independent scoring banks, the numerical contrast is a fixed
statistic. If the true null is one of the supplied points and the data and
replicates have identical IID observation/selection laws, its rank P-value
with ties included is conservative. The maximum across all null points
cannot reject when that true-point P-value would not reject. This is a
finite-grid specialization of
[maximized Monte Carlo testing](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf),
not coverage of the continuous parameter space. It does not validate
off-grid rates/Ne, misspecified detection, unmodeled curation or empirical
gene-tree inference errors. Missing-support/resource failures abort rather
than supplying a partial P-value; analyses that completed cannot be used
as a selected-subset calibration guarantee.

The MC-overlap bounds still describe score-bank error, not finite-B
uncertainty, biological confidence or a posterior probability. Grid
maximization removes dependence on choosing just one numerical generating
null, but can substantially reduce power. The least-favorable point for
the numerical P-value need not be the point for its upper MC bound.
Grid streams use SeedSequence `[seed, 2, generating_grid, replicate]`;
bank streams remain unchanged and plug-in streams remain
`[seed, 1, replicate]` (indices are zero based).

The Python calibration helper snapshots one-shot bank/observation iterables
before reuse. A custom sampler must return exactly the requested family
count; missing or extra families abort before refitting, without filtering
or retrying. Numerical failures retain their exception type and identify
the replicate and, in grid mode, the generating null point. These checks
protect batch integrity; they cannot establish that a custom sampler's
observation/error law matches the data.

CLI calibration assumes the fixed rooted input gene topologies are observed
without an additional inference-error process. For estimated gene trees,
the null must repeat the actual sequence/tree/root pipeline; the research
pilot supplies such a sampler for its JC69/NJ procedure. Its results do not
calibrate IQ-TREE or other real-data pipelines. Gene-tree bootstrap stability,
block/family resampling and the event P-value are separate quantities.

## Outputs and Evidence

Primary schema/method: `nwkit-mul-locus-mc-v1` /
`linear-DL-daughter-bounded-MLC-finite-grid-MC-v1`. Rows retain every evaluated
candidate/grid and exclusion, including the no-polyploidy candidate 0.
Columns are `mul.tree, h2.node, grid, log_likelihood, mc_lower, mc_upper,
parameters, status`, with `reason` when exclusions exist. Missing/nonfinite
numeric values are blank in TSV and null in JSON: an unobserved-pattern
log estimate/lower endpoint means negative infinity, not zero likelihood
evidence; an unbounded contrast upper endpoint means positive infinity.
The complete JSON has model settings, the supplied H1/H2 scope, a dated
species tree, ordered species-colored observations, all pattern counts/bounds,
best numerical grid points, overlapping points, contrast intervals, bootstrap
searches, exclusions and limitations. Each bank also saves its dated population
tree, population-tip-to-species map and parameters. These are the candidate
models needed to reproduce the scores, not inferred biological point trees.
It is not a GRAMPA-compatible table.

Reports/checks contain per-family hit/sample counts, integration method,
probability estimates and bounds for every evaluated bank. For stratified
integration, aggregate `hits/samples` is a **raw proposal frequency**, not
the conditional likelihood; use `probability_estimate`, which includes the
prior weights and selection normalization saved in JSON. Calibration TSV
contains every null replicate's
selected grid/candidate and contrast bounds. Grid-supremum calibration adds
`generating_null_grid`, distinct from `null_grid` reselected during scoring.
JSON identifies method `finite-grid-supremum-Monte-Carlo-test`, replicates
per point, total replicates, every point's P/bounds/parameters, and the
least-favorable point for the P-value and upper MC bound separately.
All five possible files form
one staged bundle protected against every declared input, including model
JSON. Failures preserve prior files. Only the main table may use stdout.
The shared provenance audit hashes the model and every declared output and
rejects aliases between the audit file, model and companion outputs.

See [the independent study protocol](../../examples/mul-locus/README.md),
[validation record](../validation/MUL_LOCUS_MC_VALIDATION.md), and
[subsequent numerical/reproducibility audit](../validation/MUL_LOCUS_MC_AUDIT.md).
The [separate finite-null-grid protocol](../../examples/mul-locus/CALIBRATION.md)
expands independent evaluation across turnover, detection and ILS settings;
see its [implementation and evaluation record](../validation/MUL_LOCUS_GRID_CALIBRATION.md).
The [subsequent input and runner audit](../validation/MUL_LOCUS_GRID_AUDIT.md)
records batch-integrity, reproducibility and worker-failure repairs.
The [research-only conditional integration protocol](../../examples/mul-locus/INTEGRATION.md)
and [its evaluation](../validation/MUL_LOCUS_INTEGRATION_PROBE.md) record
the first prototype's wide intervals and independent rejection-sampler caps.
The [conditional-reference/KL follow-up](../../examples/mul-locus/CONDITIONAL.md)
and [integration record](../validation/MUL_LOCUS_CONDITIONAL_INTEGRATION.md)
support explicit CLI modes and optional GeneGalleon integration only.
D+L remains the default. This is not an established production WGD detector;
the follow-up still has wide MC event bounds and a parental error in high ILS.
