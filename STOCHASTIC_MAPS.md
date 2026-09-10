# Conditional stochastic maps

`nwkit asr` can retain complete discrete-state histories, summarize time spent
in each state, aggregate changes through time, and draw state probabilities
along branches. All requested outputs use the same conditional map draws.
The existing `--stochastic-map-out` transition-count table keeps its schema and
seeded results when these outputs are added.

```bash
nwkit asr --infile examples/stochastic_maps/tree.nwk \
  --trait examples/stochastic_maps/traits.tsv --state-column habitat \
  --model ER --rate 0.65 --states forest,grassland,wetland \
  --n-sim 500 --seed 42 --threads 2 \
  --outfile nodes.tsv --model-out model.tsv \
  --stochastic-map-out counts.tsv --map-history-out histories.tsv \
  --map-summary-out durations.tsv --map-time-out time.tsv \
  --map-probabilities-out probabilities.tsv --map-figure-out probabilities.png
```

![Conditional state probabilities](examples/stochastic_maps/probabilities.png)

Each colored ribbon has constant total thickness; the thickness of each color
estimates its state's conditional probability at that position. This is a
summary of many histories, not one sampled history. Gray vertical connectors
show topology and do not represent evolutionary time. Missing tip states remain
uncertain. PNG, PDF, and SVG are supported; figure dimensions follow tip count.

## Models and sampling

Histories support one input tree and the discrete CTMC models ER, SYM, ARD, F81,
GTR, CUSTOM, MK-DESIGN, MK-REGIME, HRM, COVARION, PAGEL-INDEPENDENT, and
PAGEL-DEPENDENT. MK-MIXTURE, THRESHOLD, continuous models, and tree ensembles
are rejected. Existing transition-graph constraints apply. HRM and COVARION
histories report observed states: hidden-only changes are removed and adjacent
segments with the same observed state are merged. Pagel models retain their
JSON-encoded joint-state labels.

Node states are sampled jointly conditional on tip evidence. A uniformization
bridge is sampled on each branch conditional on its endpoints, then its event
times are sorted uniform draws conditional on the uniformization event count.
Virtual self transitions are removed. The uniformization cutoff bounds the
omitted Poisson mass relative to the smallest positive branch transition
probability (target `1e-12`), so a rare endpoint does not amplify an
unconditional tail cutoff. This also applies to count-only mapping.
Event times use a separate random stream,
so retaining them does not alter the existing transition-count random stream.
Shared node endpoints agree across all incident branches. The same seed and
simulation count produce identical results across worker counts.

These draws condition on the input tree, the fitted or supplied rate matrices,
and the selected root-prior model. They do **not** integrate parameter or tree
uncertainty. Duration quantiles describe variation among conditional histories;
they are not confidence intervals for fitted rates. Probability `mc_se` measures
Monte Carlo sampling error only, using `sqrt(p * (1 - p) / n)`; zero at an
observed frequency of zero or one is not proof of certainty.

## Outputs and coordinates

All TSV files have a header and use the original tree's `branch_id` and `parent`
identifiers. Branch records also contain `node_class` and `name`. A branch is
identified by its child node. The root has no incoming branch and no map rows.
Auxiliary outputs require distinct file paths; only the primary `--outfile` may
be `-`. When any new map output is requested, all ASR file outputs are staged
and committed together. A fitting, aggregation, or export failure preserves
existing outputs and does not emit a partial primary table to standard output.

Time is in input branch-length units, increasing from the input root (zero).
There is no implicit present date, ultrametric conversion, or calibration.
If lengths are substitutions, the time coordinates are substitutions as well.

| Option | Rows and fields beyond branch metadata |
| --- | --- |
| `--map-history-out` | One row per draw/branch/segment: `simulation`, `segment` (both one-based), `state`, local `start`, `end`, `duration`, `start_from_root`, `end_from_root`. |
| `--map-summary-out` | Every branch/state, including zero occupancy: `state`, `branch_length`, `total_duration`, `mean_duration`, `duration_sd`, `duration_q025`, `duration_q975`, `mean_fraction`, `num_simulations`. |
| `--map-probabilities-out` | Every branch/grid point/state: `position` (fraction from 0 to 1), local `distance`, `time_from_root`, `state`, `probability`, `mc_se`, `num_simulations`. |
| `--stochastic-map-out` | Existing branch/directed-state-pair counts: `from_state`, `to_state`, `total_count`, `mean_count`, `posterior_frequency`, `num_simulations`. |

History segments cover each entire branch, including branches with no changes.
A transition occurs at the start of every segment after the first; its source
state is the previous segment's state. Segments are left-closed/right-open,
with the final endpoint included. A zero-length branch has one zero-duration
segment. Its mean occupancy fraction is blank; probability rows repeat the
shared endpoint state distribution at every grid position. `duration_sd` is
the sample SD (denominator `n - 1`) and is blank for one draw. Quantiles use
linear interpolation of empirical draw durations. Duration statistics are
computed in branch-length-scaled units to avoid variance overflow or underflow.
Unrepresentable totals or lineage-time rates are rejected rather than exported
as infinite values. All totals sum over draws,
while means divide by `num_simulations`.

`--map-time-out` aggregates over equal-width bins spanning zero to the greatest
node depth. Its columns are `bin` (one-based), `time_start`, `time_end`,
`lineage_time`, `quantity`, `state`, `other_state`, `total`, `mean`,
`per_lineage_time`, and `num_simulations`.

- `quantity=duration` reports each state's summed occupancy across overlapping
  branches; `other_state` is blank.
- `quantity=transitions` reports each directed pair, excluding self changes;
  `state` is the source and `other_state` is the destination.
- `lineage_time` is the sum of branch overlap lengths in that bin, counted once
  per tree, not once per draw. `per_lineage_time` is `mean / lineage_time` (blank
  if exposure is zero). For duration rows it is a proportion of lineage time;
  for transition rows it is a realized count per lineage-time unit, not a
  fitted generator entry or a state-specific hazard estimate.

Bins are left-closed/right-open, except the final right endpoint is included.
Events exactly on a bin boundary enter the following bin. Bins with no observed
changes still have explicit zero rows. Summing duration means over states and
bins recovers total branch length; summing transition totals recovers the
transition-count table. A tree without a positive time span cannot produce
`--map-time-out`; its other tables may be empty or contain zero-length branches.

## Controls and limits

Any map output enables `--n-sim` (default 100), `--seed`, and `--threads`
(default 1). Counts output is optional. `--map-time-bins` requires time output
(default 20; range 1–10,000). `--map-grid-points` requires probabilities or a
map figure (default 51; range 2–10,001, including both branch endpoints).
More draws reduce Monte Carlo error; a finer grid only refines the sampled
probability display. Underscore aliases are available for the new options.

Histories are retained in memory. Preflight limits include 100,000 draw/node
records, a 256 MiB bound on the dense simulation arrays before histories,
a conservative 1,000,000 uniformization-segment bound, and the existing
CTMC mapping work limits. Probability, duration, and time outputs are limited
to 500,000 rows each; time aggregation additionally allows at most 20,000,000
segment/bin intersections. Probability evaluation is limited to 20,000,000
draw/grid-point pairs across branches. Figures support at most 20 states and
200 tips; state and tip labels are rendered literally, including dollar signs.
Reduce draws, grid points, bins, or input size when the applicable limit is
exceeded. Branch depths that lose more than `1e-8` relative branch-length
precision on the root-time axis are rejected. History and probability exports
also reject distinct local times that collapse to the same root-time value.
Duration and time-bin summaries use local branch coordinates to avoid
subtraction of large root depths. Reduce the disparity between root depths and
branch/event durations for coordinate exports; multiplying every branch by the
same factor does not generally resolve that disparity.

See the [runnable example](examples/stochastic_maps/README.md) and
[ASR guide](ASR.md) for input and model conventions.
