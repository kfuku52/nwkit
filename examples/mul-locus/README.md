# DL + ILS null-comparison pilot

For the subsequent, separately frozen finite-null-grid study, see
[calibration protocol](CALIBRATION.md). The pilot protocol below is unchanged.
For the later research-only integration comparison, see
[conditional integration protocol](INTEGRATION.md) and its linked evaluation.

This research study evaluates the experimental `locus-mc` model, not a
production WGD detector. Its complete generative observation distribution
includes continuous linear locus birth/death, one ancestral origin locus on
a finite DL stem, daughter-bounded coalescence, independent detection, and
the explicit family selection `2 <= observed tips <= 4`. No DL histories
are silently truncated. Resource-cap failures remain failures.

## Frozen Protocol

Species chronogram `((A:1,X:1):1,B:2);`, ancestral DL stem 0.5 generations,
known detection 0.9 per locus. Independent forward simulation uses a global
SSA over simultaneous loci (not the production depth-first branch sampler).
msprime simulates full-locus ancestry, accepting a genealogy only if each
SSD daughter has coalesced by its birth time. Undetected loci remain sampled
during that conditioning and are removed afterwards. Root and tip counts
are not altered to make topology comparisons agree.

Cases: no polyploidy with (duplication, Ne) `(0.03,0.3)` and `(0.08,1)`;
direct-parent allopolyploid H2 A or B with duplication 0.05, Ne 1,
attachment age 0.5. Loss rate is 0.03 in every case. The supplied search grid
crosses duplication `0.03/0.08`, Ne `0.3/1`, age `0.3/0.7`, and loss 0.03.
True allopolyploid duplication/age are not on that grid; fitting never
receives a per-dataset true parameter or parent. Null grids omit redundant
attachment ages, which have no meaning without polyploidy.

Defaults: seed `20261025`; 5 independent replicates x 30 families per case,
20 datasets total. Fit banks use 10,000 selected simulations per candidate/
grid point, seed `20261026`, with simultaneous 99% Monte Carlo bounds.
The finite-grid maximum is not continuous MLE, and MC intervals are not
biological confidence intervals. Estimated trees use JC69 600-site sequences,
rate 0.01/site/generation, NJ and midpoint rooting without true-root access.

Both views use the identical family sets and model/grid/candidate search.
True-tree calibration runs 39 fresh null replicates. Estimated-tree
calibration runs 19 and repeats sequence simulation/NJ/rooting under the
fitted null; it does not substitute a true-gene-tree null distribution.
Every bootstrap reselects null parameters and searches every alternative
parent/grid point. P-values use the finite-sample `(1+exceedances)/(B+1)`
correction; estimates and conservative MC-overlap bounds are separate.
Gene-tree/bootstrap stability is not converted into event-test support.

Report null false positives and their denominator, parent recovery, reported
versus MC-unresolved decisions, tree error, and all failures. At alpha 0.05
an event is reported only if the upper MC bound on the bootstrap P-value
is <=0.05. A tiny pilot cannot establish general false-positive control;
the plug-in bootstrap is conditional on a finite model/grid, not uniform
composite-null coverage. Do not change this protocol after inspecting
results to force successful calibration or power.

## Integration Resource Follow-Up

The first frozen run (`/tmp/nwkit-locus-pilot-20261002`) completed only 17/40
analyses; 23 had missing MC support, mostly inside null replicates rather
than in observed datasets. This is an integration limitation, not evidence
that an observed topology is biologically impossible.

The follow-up retains the identical scientific model, evaluation dataset
seed, family sets and thresholds, and uses an independent bank/bootstrap
seed `20261027`, a fixed 100,000 **raw** draws per bank, and explicit
`ancestral-stratified` integration. Before rerunning it, tests compare the
stratum-weighted ancestor/selection distribution to the independent linear
BD closed form. No-event and first-birth ancestral histories are sampled
separately, then reweighted by their original prior AND their selection
mass; initial-stem first-loss histories have zero selected probability.
Simultaneous MC bounds include both strata and selection denominators.
The first failed experiment is retained; it is not replaced by this run.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/mul-locus/pilot.py \
  --output /tmp/new-stratified-locus-pilot --samples 100000 \
  --bank-seed 20261027 --integration ancestral-stratified
```

## Separate Gene-Tree Stability Protocol

After the complete paired event-calibration study, `stability.py` resamples
600 alignment columns within each family, using the same column indices for
all sequences in a family, then repeats JC69 distance/NJ/midpoint rooting.
It keeps all original families and searches every original null/alternative
grid bank. Default seed `20261029`, 19 replicates for each of all 20 datasets;
there is no case selection or resource/threshold tuning from these results.
Every missing-support or tree-inference failure remains a failed replicate.
Input alignment/bank/source hashes and every full fit are retained.

This is **descriptive parent-ranking stability**, not an event P-value,
posterior support or biological confidence interval. Its numerical parent
frequencies, MC-overlapping parents and failures are reported separately
from the plug-in search-wide null calibration. No stability frequency is
converted to a WGD decision, and no independence between bootstrap fits is
claimed beyond their conditional resampling construction.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/mul-locus/stability.py \
  --study /tmp/nwkit-locus-stratified-20261002 \
  --output /tmp/new-locus-site-bootstrap
```

## Run and Outputs

Research dependencies `msprime` and `biopython` belong in an isolated
environment; no new NWKIT runtime dependency is required. From the repo root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/mul-locus/pilot.py --output /tmp/new-locus-pilot
```

The output must be a fresh directory. Protocol/source hashes, model JSON,
independent locus histories and acceptance audits, alignments, gene trees,
all fit/bootstrap results and a complete summary are retained. Each failed
analysis is retained and counted; the study exits nonzero if any fail.
Frozen banks are independent of evaluation datasets and reusable only with
their identical model/source snapshot. Docker/SIF, GeneGalleon adoption,
external empirical analyses, ML gene-tree pipelines, arbitrary copy roles,
population heterogeneity and other event modes are later work.
