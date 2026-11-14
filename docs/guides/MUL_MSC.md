# Conditional Allopolyploid MUL/MSC Prototype

This guide describes conditional MSC only. The separate experimental
[locus DL + ILS model](MUL_LOCUS_MC.md) includes ordinary duplication/loss
and a normalized null comparison with explicit Monte Carlo error.

`mul-reconcile --score-model msc` compares **second-parent candidates, assuming
one disomic allotetraploid event already occurred**. Incomplete lineage sorting
(ILS) is explained by a multispecies coalescent (MSC), not automatically counted
as duplications or losses. This is a first research prototype, not a WGD test,
not full AlloppNET, and not a duplication-loss-coalescent likelihood model.
The default [D+L mode](MUL_RECONCILE.md) and GeneGalleon defaults are unchanged.

## Example

Species input, with fixed node ages and a shared population scale:

```text
(((A:2,X:2):1,B:3):2,C:5);
```

One observed gene family:

```text
((a_A,x1_X),(b_B,x2_X));
```

Here `X` is assumed allotetraploid. Its stem exists between ages 0 and 2.
The supplied hybridization age is 1. Candidate H2 branches `A`, `B` and `C`
all exist at that age. Input/output paths in this example should be distinct:

```bash
nwkit mul-reconcile -i genes.nwk --species-tree species.nwk \
  --species-regex '.*_([^_]+)$' --score-model msc --h1 X \
  --species-time-unit coalescent --hybridization-age 1 \
  -o likelihoods.tsv --report gene_likelihoods.tsv \
  --model-out msc.json --tree-out best_dated_mul.nwk
```

`--h2` can restrict the candidate list using the usual selectors. H1 must be
one fixed non-root clade; unrestricted H1 scans and arbitrary supplied
`--multree yes` topologies are deliberately unsupported in this mode.

Optional [bounded parameter estimation](MUL_MSC_FIT.md) uses `--msc-fit
age|ne|joint` and a separate fitted-output schema. The following describes
the default fixed-parameter mode.

## Model and Units

All species lengths must be explicit, finite, nonnegative, and ultrametric.
Tips are contemporaneous at age 0; small floating discrepancies are checked
using an absolute tolerance of `1e-10 * maximum_input_branch_length`. The root
stem must be omitted/zero; an infinite ancestral population above the root
allows the remaining lineages to coalesce. Both gene and species trees are
rooted and strictly binary. Gene branch lengths and internal labels are ignored:
this model uses gene **topology probabilities**, integrating coalescence times.

Time units are never guessed:

- `--species-time-unit coalescent`: every input length and hybridization age
  is in one **shared** scale of `2*Ne` generations. Do not supply Ne again.
- `--species-time-unit generations`: additionally supply
  `--effective-population-size N`; all branches use the fixed positive diploid
  or subgenome Ne `N`, and durations are divided by `2*N`. Hybridization age
  stays in generations. Calendar years and substitutions/site are not accepted
  as generations without an explicit upstream conversion.

In fixed mode no population size, node age, hybridization age, or mutation rate is estimated.
The duplicated descendant subtrees share their node ages and the fixed Ne.
The hybridization age must lie strictly inside H1's stem. H2 is the **direct
second-parent branch** donating at that same age; its branch must exist then.
The existing H1 stem represents the first parental lineage. The candidate is
unfolded as H1's copied subtree attached alongside H2, at the supplied age.
This does not add a separately dated unsampled donor that split from H2 earlier.
The constant population scale means no founder bottleneck or population jump
at hybridization is modeled. These restrictions are narrower than AlloppNET.

Baseline no-polyploidy and autopolyploid candidates remain in the primary table
as `excluded`, with reasons, but receive **no likelihood**. Time-incompatible
H2 candidates are also explicit exclusions. If no candidate can be evaluated,
the analysis fails rather than emitting an empty successful ranking.

## Sampling and Missing Copies

Each gene family must contain at most one copy from each diploid species and
at most two distinct homoeologs from each H1 species. One sampled copy per
subgenome is assumed. Multiple individuals/alleles, ordinary small-scale
duplicates (SSD), isoforms, and collapsed/chimeric sequences are not supported;
the software cannot infer their identity from copy counts. Matching is through
the shared species parser, regex, or species mapping TSV. Excess or unmatched
copies fail before scoring; no candidate-specific family filtering occurs.

When homoeolog identity is unknown, all injective species-preserving mappings
are considered. Two copies from one polyploid species must occupy different
subgenomes; they cannot both be assigned to the easier subgenome. The prior is
uniform over the permissible mappings **within each family**, normalized by
their count. For family `g` with `M_g` mappings:

```text
P(G_g | H, parameters, observed copies)
  = (1 / M_g) * sum_a P_MSC(G_g | dated_MUL(H), assignment a, parameters)

log_likelihood(H) = sum_g log P(G_g | H, parameters, observed copies)
```

This is a mean, not a maximum assignment probability or unnormalized sum.
The model conditions on observed copies: absent species/homoeologs are not
scored as losses or evidence against any candidate. Biased loss/detection can
invalidate the uniform mapping assumption. Families are treated as independent;
linked loci and synteny-derived joint assignments are not modeled.

## Computation and Outputs

An ancestral configuration is an antichain of nodes of the supplied gene tree.
Only mergers consistent with that tree are retained; incompatible mergers
lose probability mass rather than becoming duplication events. For a branch
starting with `k` lineages, `k*(k-1)/2` is the total merger rate; each pair has
rate 1. The algorithm counts all compatible ordered paths to each configuration,
multiplies by the integrated pure-death transition probability, and sums
configurations/mappings in log space. Above the root, all possible compatible
merger orders are integrated analytically to infinity.

Small-time lineage transitions use a positive uniformization series with an
omitted-tail bound of `1e-15` relative to the accumulated probability. Other
transitions use a rate-shifted matrix exponential so very small probabilities
remain finite log probabilities (e.g. a discordant topology after 10,000
coalescent units). At very long times, the slowest hypoexponential term is used
only if the sum of the absolute remaining terms is bounded below `1e-15` of
that term; this avoids overflowing a matrix exponential at extreme durations.
Positive log probabilities up to `1e-10` from floating
roundoff are clamped at probability 1; larger violations fail. There is no beam
search or selected-history cutoff.

`--max-coalescent-states` (default 100000) bounds per-family/candidate work:
configuration visits, merger transitions, joins, propagation updates and the
squared dimension of each new lineage-transition calculation.
`--max-coalescent-assignments` (default 10000) bounds family assignments.
Exceeding either fails the whole run without truncation or replacing results.
Legacy `--max-state-pairs`/`--max-maps` are DL-only.

The primary TSV schema is `nwkit-mul-msc-likelihood-v1`:

- `mul.tree, h1.node, h2.node`: run-local hypothesis identifiers/selectors.
- `log_likelihood`: larger is better; blank for excluded candidates.
- `delta_log_likelihood`: best minus candidate; not a P-value or probability.
- `dated.tree`: generated MUL-tree with lengths in the declared input unit.
- `status, reason`: evaluated or explicitly excluded.

`--report` and `--check-out` both contain all evaluated gene/candidate pairs:
`mul.tree, gene.tree, log_likelihood, num_assignments, coalescent_states`.
The last field is the above work counter, **not** the number of unique states.
They export no duplication/loss counts or GRAMPA map serialization.
`--tree-out` contains the highest-likelihood candidate with dated lengths,
explicit rooted declaration and 17-significant-digit branch serialization.
Exact floating-score ties prefer the lowest run-local candidate ID.
`--model-out` records method `conditional-direct-parent-MUL-MSC-topology-v1`,
schema, fixed units/parameters, normalized mapping policy, all exclusions,
gene-level contributions, and limitations. Only the primary TSV can use stdout.
All files use NWKIT's staged bundle and input-alias protection.

These likelihood outputs must not be passed to GeneGalleon's legacy GRAMPA
summary reader. No GeneGalleon mode/default has been changed to consume them.

## Evidence and Next Stages

See [prototype validation](../validation/MUL_MSC_VALIDATION.md). Analytic,
independent exhaustive, and known-parameter pilot evidence is distinct from
biological accuracy on estimated gene trees. Scores are not posterior candidate
probabilities, confidence intervals, or calibrated tests of polyploidy.

Remaining stages are a duplication-loss-locus layer with normalized history
integration and observation/ascertainment rules; fair no-polyploidy comparisons;
general parameter estimation and calibrated uncertainty; broader estimated-tree,
sequence and misspecified-model benchmarks; and opt-in GeneGalleon consumers/runtime
validation. Bounded shared-Ne/attachment estimation and a small estimated-tree pilot
are documented separately above; they do not complete these production gates.
Autopolyploid inheritance, multiple events, subgenome exchange, and gene flow
need additional models. A lower D+L score is not interchangeable with a higher
conditional MSC likelihood.

Scientific foundations: [Jones et al. (2013), AlloppNET/AlloppMUL](https://pubmed.ncbi.nlm.nih.gov/23427289/),
[Degnan and Salter (2005), gene topology distributions](https://onlinelibrary.wiley.com/doi/10.1111/j.0014-3820.2005.tb00891.x),
and [Rasmussen and Kellis (2012), DLCoal](https://genome.cshlp.org/content/early/2012/01/23/gr123901111).
The likelihood implementation here is independent and more restricted than
the first paper; it does not claim to implement that paper's Bayesian analysis.
