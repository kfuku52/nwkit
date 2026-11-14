# Focal-Lineage Ks Correction

`nwkit ksrate` calculates species-node distances on each focal species' Ks
scale from between-species ortholog comparisons. It is a native trio-based
correction, not the ksrates program, a WGD detector, or an absolute-age model.

```sh
nwkit ksrate -i species.nwk --ks-tsv ortholog-ks.tsv \
  -o node-ks.tsv --trios-out trios.tsv --model-out ks-model.json \
  --bootstrap 199 --seed 71
```

## Inputs And Formula

The rooted species tree must have unique named tips and no unary nodes.
Branch lengths are not used in the correction. The TSV has
`species_a,species_b,family_id,ks` and one representative ortholog distance per
unordered species pair and family. IDs are literal strings; duplicate
pair/family rows, tree-absent species, within-species pairs, negative/nonfinite
Ks and malformed tables are rejected. Select and audit representatives upstream
instead of giving more weight to high-copy families.

For focal F, sister S and an outgroup O outside their species-node clade:

```text
corrected node Ks = Ks(F,S) + Ks(F,O) - Ks(S,O)
```

For an additive distance, this is twice the focal-lineage synonymous distance
from the F/S ancestor to F. Each pair distance is the median of its family
observations; a node's correction is the median across complete sister/outgroup
trios. Medians need not be exactly additive. This estimator choice is explicit
and is not a claim of equivalence to ortholog-Ks mode estimators or a fitted
substitution tree.

By default every tip is focal, each ancestral node contributes sisters outside
the focal child clade, and outgroups come from the node's immediate parent
clade. `--outgroup-policy all` includes all tree tips outside the node.
`--focals` restricts focal species. A work bound rejects more than 100,000 trios.
The species root has no external outgroup and remains unresolved; neither its
age nor an ancient WGD anchor is guessed from the data.

## Uncertainty And Diagnostics

The bootstrap jointly resamples shared family IDs across all species pairs.
A family gets the same multiplicity in every pair, preserving across-pair
dependence in the collected data. Percentile intervals condition on the
species tree, inferred families, selected ortholog pairs and Ks estimates.
They do not capture all topology, alignment, substitution-model, saturation,
gene-conversion or representative-selection uncertainty. They are neither
posterior probabilities nor genome-event age confidence intervals.

The default `--ci-method family-bootstrap-percentile` is an approximate
finite-sample interval, not guaranteed nominal coverage. In particular,
percentile intervals for medians can under-cover with modest family counts.

For conservative event placement, use `--ci-method pair-median-bonferroni`.
This bounds every compared population pair median using noninterpolated
binomial order statistics, with the total error budget divided over all used
species pairs. Bounds propagate through the trio formula and the node median.
Under independent identically distributed family observations within each pair,
all reported population corrections are covered simultaneously with probability
at least `--ci-level`. Dependence between different species pairs is allowed.
This is a conditional sampling guarantee, not a guarantee that the biological
orthologs, substitution model or evolutionary interpretation are correct.

The simultaneous method can produce much wider intervals, especially for
long outgroup distances. Too few independent families make the entire node
interval unavailable; unbounded complete trios are not silently dropped.
Repeated copies from the same family do not increase the independent sample
size. It is available with `--bootstrap 0`, since it has no Monte Carlo step.
When bootstrap draws are requested, `bootstrap_ci_lower`, `bootstrap_ci_upper`,
and `bootstrap_interval_status` retain the original percentile diagnostics.
The order-statistic construction is described by
[NIST](https://itl.nist.gov/div898/software/dataplot/refman1/auxillar/mediancl.htm);
NWKIT uses the conservative integer ranks, not interpolated endpoints.

Missing species-pair comparisons remain missing. If any bootstrap draw cannot
estimate a node, its interval is unavailable and the number of estimable draws
is recorded; incomplete draws are not dropped to produce a narrower interval.
With the default method, `--bootstrap 0` returns points and diagnostics without
intervals. Missing bootstrap draws do not invalidate separately calculated
finite-sample simultaneous bounds, whose availability is reported independently.

Negative trio corrections, negative node medians, incomplete comparisons and
nonmonotone older-versus-younger node distances are explicit diagnostics.
`raw_corrected_ks` retains negative values; `corrected_ks` is unavailable for a
negative median. No clipping, isotonic projection or clock fit hides these
inconsistencies. Bootstrap endpoints refer to raw corrections and can be
negative, which also signals uncertainty in a nonnegative distance.

## Outputs

The primary TSV has focal species, stable `species_event_id`, input `branch_id`,
descendant taxa, raw/usable corrected Ks, interval endpoints and method, draw
counts, trio counts, monotonicity and status. The trio TSV retains all three
pairwise medians, per-comparison family counts and raw corrections. JSON records
the formula, estimators, bootstrap unit, seed and limitations.

Only the primary output may use stdout. All related file outputs are staged
together and protected against input/output path aliases. See
[CLI conventions](CLI_TSV_CONVENTIONS.md).

The correction is motivated by lineage-specific rate adjustment in
[Sensalari et al. (2022)](https://doi.org/10.1093/bioinformatics/btab602).
Hand-computable asymmetric-rate examples, shared-family bootstrap dependence,
unresolved roots, missing draws, negative/nonmonotone diagnostics and real CLI
execution are checked in `tests/test_ksrate.py`.
