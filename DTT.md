# Disparity through time

`nwkit dtt` measures how continuous-trait diversity is distributed within
surviving clades from the crown root to the present. One selected column gives
scalar DTT; multiple columns give a joint squared-Euclidean DTT curve.

```sh
nwkit dtt -i examples/dtt/tree.nwk --trait examples/dtt/traits.tsv \
  --columns size,shape --n-sim 999 --seed 42 --threads 2 \
  -o dtt.tsv --summary-out summary.tsv --figure-out dtt.png \
  --model-out model.json
```

![Example observed DTT and Brownian comparison](examples/dtt/dtt.png)

## Inputs and definition

Use a rooted ultrametric tree with finite, nonnegative branch lengths, positive
crown age and at least three complete tips. Relative tip-depth spread must be
at most `1e-8`. Event times are normalized from accumulated input depths to
keep all split times within `[0,1]`, including zero-length present-day splits. The tip-keyed TSV uses `leaf_name`; shared missing-token and
unmatched-row policies are described in [CLI conventions](CLI_TSV_CONVENTIONS.md).
`--missing error` is the default. `--missing drop` omits incomplete tips jointly
and rebases time at their retained most recent common ancestor. The original
tree must be ultrametric even when tips are dropped. Original branch IDs are
preserved in the clade table.

Disparity is the mean squared Euclidean distance between pairs of tips:
`2 * sum_i ||x_i - mean(x)||² / (n - 1)`. Each clade's disparity is divided by
the retained crown clade's disparity. At each split, DTT is the unweighted mean
of these ratios over active clades containing at least two retained tips.
Singleton lineages are excluded. If no such clades remain, DTT is zero.
Ratios can exceed one. Polytomies and zero-length internal branches are allowed;
simultaneous split events are combined. Unary nodes do not introduce split events.

The curve begins at one immediately before the crown split, then records its
post-split value at the same time zero. Relative time increases toward one at
the present. A final zero point is included even when the last split is earlier.
These conventions reproduce `geiger::dtt(index="avg.sq")` after combining its
repeated event times and adding the zero tail. See the
[geiger reference](https://search.r-project.org/CRAN/refmans/geiger/html/dtt.html)
and the runnable [R comparison](examples/dtt/reference.R).

`--scale raw` preserves squared-Euclidean geometry, so input units affect the
relative importance of traits. `--scale standardize` divides each trait by its
observed sample SD before analysis. Numerical centering and scaling are recorded
in the model JSON. Every selected trait must be nonconstant.

## Brownian comparison and MDI

By default, 999 simulations evolve joint traits along the retained tree under
Brownian motion (BM). The fitted rate covariance uses GLS residual cross-products
with divisor `n - 1`, equivalent to independent-contrast rate estimation. The
simulation starts at a fixed root; root translation does not affect disparity.
Standardization, when selected, is held fixed across simulations.

The shaded band is a pointwise BM null simulation envelope, conditional on the
fitted rates, scaling and tree. It is **not a confidence band** for the observed
curve. It excludes rate-estimation and tree uncertainty. `--ci-level` sets its
coverage fraction (default `.95`). Both BM mean and median are exported; the
figure and MDI use the pointwise median.

The morphological disparity index (MDI) integrates observed DTT minus the BM
median using linear interpolation and trapezoids. Positive MDI indicates greater
within-clade disparity than the BM reference; negative MDI indicates less.
`--mdi-range start,end` restricts integration in relative time (default `0,1`),
with interpolated boundary values. Exported BM-null MDI quantiles and the right
histogram integrate each simulated curve against the same pointwise median.
They are conditional simulation summaries, not MDI confidence intervals or
calibrated hypothesis-test p-values.

BM simulations require a nonsingular tip covariance, more tips than traits and
linearly independent trait columns. `--n-sim 0` permits descriptive DTT without
these BM restrictions and leaves BM/MDI fields empty. A star tree or a tree with
no split strictly between crown and present returns `uninformative_topology`:
its DTT timing contains no internal diversification information.

## Outputs and reproducibility

The primary TSV (`-o`, stdout by default) has `time_index`, `phase`,
`relative_time`, `time_from_crown`, `relative_disparity`, `num_clades`,
`bm_mean`, `bm_median`, `bm_lower`, `bm_upper`, and `num_simulations`.
`phase` distinguishes the two time-zero rows. `num_clades` counts only the
non-singleton clades entering the average.

Optional outputs:

- `--summary-out`: one-row TSV with dimensions, scaling, crown age, integration
  interval, observed/BM median areas, MDI, BM-null MDI quantiles and status.
- `--clades-out`: original branch and parent IDs, name, node class, retained tip
  count and relative disparity for every retained node (tips have zero).
- `--simulations-out`: all simulated curves in long form, keyed by one-based
  `simulation` and `time_index`.
- `--model-out`: JSON with used/excluded taxa, transformations, fitted rate
  covariance, seed and statistical conventions. Rates use transformed trait
  units per unit relative crown time. To recover original trait units per input
  tree-time unit, multiply entry `(i,j)` by `data_scale[i]*data_scale[j]/crown_age`.
- `--figure-out`: PNG, PDF or SVG with retained-tip traits above DTT and, when
  simulated, the MDI histogram. The tree and DTT share an aligned time axis.
  Gray branches show chronology, without ancestral-state reconstruction.

### Compact trait display

`--figure-layout heatmap` (default) always uses a single tree plus heatmap,
including when only one trait is displayed. `auto` is an alias for this same
layout. Use `--figure-layout trees` to request individual trait trees.

`--figure-columns` selects an ordered subset of `--columns` for display only.
It does not change the traits entering DTT, missing-tip filtering, fitted BM
covariance, or MDI. The figure reports both analyzed and displayed trait counts.
For example, calculate three traits but display two:

```sh
nwkit dtt -i examples/dtt/tree.nwk --trait examples/dtt/traits.tsv \
  --columns size,shape,performance --figure-columns size,performance \
  --figure-layout heatmap --figure-out heatmap.png
```

In a heatmap, rows follow the tree's tip order and columns follow the requested
trait order. `--figure-scale standardize` (default) colors each column by its
observed z score: subtract its mean over retained tips and divide by its sample
SD. Blue indicates below-mean values and red above-mean values. This display
transformation is independent of the analysis option `--scale`.
`--figure-scale raw` uses a shared color scale in original units, suitable when
traits have comparable units. Individual-tree mode always shows original values
with separate color bars; `--figure-scale` applies only to heatmaps.

Tip names are labeled up to 40 tips (with numeric values in individual-tree
mode). Larger figures retain all tip rows and colors without tip text labels.
Unicode labels require an installed font covering their characters; missing
coverage produces an error instead of silently replacing text with empty glyphs.
Raw color scales use dimensionless internal coordinates with original-unit
labels, preserving contrasts even in very small or very large units.

Figures exceeding 32 million pixels are rejected; use a heatmap or fewer
`--figure-columns` to reduce display size without changing the analysis.

All file outputs are committed together; failed computation or drawing preserves
existing files and emits no partial primary stdout table. Output paths must be
distinct and cannot replace inputs. Multiword flags also accept underscore aliases.

`--seed` defaults to 1. Each simulation has its own indexed random stream, so
changing `--threads` (1–32 worker processes) preserves results in the same software
environment. Simulation counts may be 0 or 2–10,000. Input is limited to 2,000
tips, 4,000 nodes and 64 traits. Preflight limits cap retained simulation/time
values at 2 million, exported simulation rows at 500,000 and estimated simulation
work at 100 million operations; reduce simulation count or data dimensions if a
limit is exceeded. These are work guards, not statistical adequacy criteria.
