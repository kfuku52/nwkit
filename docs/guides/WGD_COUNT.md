# Native gene-count DL/WGM analysis

`nwkit wgd-count` searches non-root species-tree branches for one retained
genome multiplication at a time. It is an experimental, native count model,
not a substitute for synteny, gene-tree attribution, or an allopolyploid model.
No known WGD branch is required. Multiplicity is supplied, not inferred.

```sh
nwkit wgd-count -i species.nwk --counts families.tsv \
  --bootstrap 199 --seed 71 -o events.tsv --model-out count-model.json
```

With the default `--bootstrap 0`, results are exploratory candidates: p-values
are `NA` and `count_support=not_calibrated`. Increasing the number of bootstrap
datasets improves Monte Carlo precision, not model adequacy. A p-value is never
an event posterior probability. The command can be computationally expensive.

See the [scientific audit](../validation/WGD_SCIENTIFIC_VALIDATION.md) for
independent numerical and simulation evidence. In its small deliberately
misspecified local-SSD pilot, the branch-burst AIC gate still admitted conditional
WGD support in 1/12 datasets. This is not calibrated protection against SSD, and
the pilot does not establish population error rates.

## Inputs

The rooted species tree needs unique named tips and all non-root lengths,
expressed in consistent units. Polytomies are allowed; unary nodes are not.
Rooting follows the shared [CLI conventions](CLI_TSV_CONVENTIONS.md).

The count TSV has `family_id` plus exactly one column per species-tree tip.
Rows have unique nonempty family IDs and nonnegative integer counts. Missing
markers follow `--missing-values`; absence of a gene is `0`, unavailable data is
missing. Every included family must have at least one observed copy.

```text
family_id	A	B	C
f001	2	2	1
f002	1	NA	0
```

The observation/ascertainment contract matters:

- Families must descend from an ancestral family present at the species root.
  Later de novo family origins, horizontal transfer, arbitrary gene-family
  splitting/merging, and filtering to selected count patterns are not modeled.
- The likelihood is conditional on at least one observed copy, with separate
  conditioning for each missingness pattern. `--ascertainment root-clades`
  instead requires observation in every immediate root-child clade and uses
  that selection probability in the likelihood. This conservative selection
  excludes lineage-restricted families without treating the excluded rows as
  losses. Every supplied row must satisfy the selected condition. Arbitrary
  additional count-pattern filters are not implemented.
  Exact conditioning for a multifurcating root has exponential cost in the
  number of root-child clades; a highly unresolved root can be computationally
  prohibitive. Do not resolve it arbitrarily to accelerate inference.
- `--detection-tsv` accepts `leaf_name,detection_probability` for known binomial
  detection probabilities in `(0,1]`. It must exactly cover the tips. These
  probabilities are supplied, not estimated, and are not substitutes for
  modeling correlated annotation failures or uncollapsed isoforms/haplotypes.

Background rates default to separately estimated duplication/loss rates on
terminal and internal branches. `--rate-model homogeneous` uses one pair.
`--rate-groups-tsv` overrides these with `branch_id,regime` rows covering every
non-root branch. IDs are the input IDs from `nwkit nwk2table`, not row numbers.
Too many unconstrained rate groups can make a genome event unidentifiable.

Family rates default to a four-category, mean-one gamma quantile approximation
with fixed shape 1. `--family-gamma-shape none` selects homogeneous families;
otherwise vary the fixed shape/category count to check sensitivity. It is not
continuous gamma integration or inference of the gamma shape.

## Probability Model

Each copy follows a linear birth-death process with group-specific duplication
rate lambda and loss rate mu. For a single copy, extinction probability is `a`
and the positive descendant count is geometric with parameter `b`. Independent
ancestral copies are convolved. These branch probabilities include paths that
temporarily exceed the represented state range and then return.

At multiplication factor `m` and retention `q`, `n` copies become
`n + Binomial(n*(m-1),q)`. This retains the original copy at the event, with
subsequent loss in the branch DL process. `q=0` exactly recovers the background.
This is a retained-copy pulse, not an instantaneous complete chromosome-doubling
and correlated-loss model. The positive geometric root mean is estimated.

Pruning sums latent copy counts and averages family rate categories. Survival
conditioning is applied **after** averaging categories. Missing tips have a
flat observation likelihood. Zero observed copies have a binomial observation
likelihood, not the same treatment as missing data.

State probabilities and root tails are never renormalized to hide omitted
mass. Every fitted model is checked by doubling the maximum count state and
requiring the maximum per-family log-likelihood change to be below
`--state-tolerance`. Failed optimization or truncation checks abort the scan.
This is a numerical convergence diagnostic, not a mathematical uniform error
bound over all parameter values. Numerical rate bounds and root-mean upper
bounds are recorded in JSON; inspect `nuisance_bound_reached`.
The TSV reports this event-fit flag separately from
`background_nuisance_bound_reached` and `branch_burst_nuisance_bound_reached`.
Inspect all three: an event fit away from numerical bounds does not certify
its fitted null or SSD competitor. These diagnostics do not change the native
support status or provide calibration against SSD heterogeneity.

## Search And Calibration

Default event fractions `.25,.5,.75` are a finite grid measured from the branch
parent. Each candidate is separately refitted with one genome event. This is
not a simultaneous multiple-event search and does not produce dated-event
confidence intervals. `--candidate-branches` restricts the search; bootstrap
calibration repeats that exact branch set and fraction grid.

The statistic is twice the event-versus-background log-likelihood improvement.
The bootstrap simulates from the fitted null, preserves tip detection and each
family's missingness mask, conditions on observation, refits the null, and
repeats the complete candidate search. A candidate's corrected p-value uses
the simulated **maximum** statistic over all searched branches. The plus-one
estimate is `(1 + exceedances)/(1 + draws)`. Failed simulated refits abort;
draws are not silently discarded. This handles the finite search and retention
boundary but remains a plug-in, model-conditional calibration, not a supremum
test over nuisance parameters or guaranteed composite-null error control.

Each candidate also has an independently fitted branch-burst alternative:
the branch gets its own continuous DL rates without a genome pulse. The
`burst_minus_event_aic` comparison is a conditional diagnostic, not a calibrated
test of arbitrary SSD heterogeneity. It counts the fitted continuous parameters,
not the searched event-fraction grid. The bootstrap, not AIC, handles grid search.

## Outputs And Interpretation

The primary TSV uses stable `species_event_id` clade hashes, input `branch_id`,
descendant taxa, event fraction, retention, likelihoods, branch-burst diagnostic,
bootstrap p-value/Monte Carlo SE, state error and explicit support status.
`count_supported_conditional` requires a bootstrap p-value at/below `--alpha`
and a positive branch-burst diagnostic. It is **not** a WGD conclusion: synteny,
family-definition quality, detection and rate sensitivity remain necessary.

JSON records rates in inverse input branch-length units, conditioning, root
prior, family categories, branch groups, state checks, bootstrap statistics,
seed and limitations. Each candidate also records its branch-burst group
assignments, so those fitted rates can be mapped back to branches. Related file
outputs are staged together and cannot
replace declared inputs. Only the primary TSV may use standard output.

Count likelihoods cannot classify individual gene-tree duplication nodes.
They also cannot distinguish allopolyploid donor divergence from the genome
multiplication time. Use unresolved labels when those distinctions are needed.

## Verification

`tests/test_wgd_count_model.py` compares transitions with an independently
constructed large-state CTMC matrix exponential, brute-force pruning,
probability normalization, limiting cases, units and unbounded simulations.
`tests/test_wgd_count_fit.py` checks candidate recovery and branch-burst
comparisons on small examples. These checks validate implementation behavior,
not empirical false-positive rates across real genome-analysis conditions.

The scientific foundation is the gene-count birth-death/retained-WGD class of
[Rabier, Ta and Ane (2014)](https://doi.org/10.1093/molbev/mst263), with the
importance of rate heterogeneity and event-search calibration illustrated by
[Zwaenepoel and Van de Peer (2020)](https://doi.org/10.1093/molbev/msaa111).
This command is an independent, limited implementation, not a claim of software
or inferential equivalence to those studies.
