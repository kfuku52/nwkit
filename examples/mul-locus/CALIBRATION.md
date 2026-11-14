# Finite Null Grid Calibration Study

This is a new fixed protocol, not a revision or replacement of the
[original pilot](README.md). It evaluates `grid-supremum` calibration using
independent SSA/msprime datasets and fresh integration banks. The original
pilot, stability experiment, scripts, outputs, and source-hash records remain
historical evidence. These new results do not calibrate empirical pipelines.
See the [implementation and evaluation record](../../docs/validation/MUL_LOCUS_GRID_CALIBRATION.md).

## Frozen Before Evaluation

Species chronogram `((A:1,X:1):1,B:2);`, H1 X, H2 A/B, one ancestral origin
locus, DL stem 0.5 generations, known species detection, disomic direct-parent
allopolyploidy, and selection `2 <= observed tips <= 4` are retained.
The search grid crosses the two listed duplication rates, two diploid Ne
values, and attachment ages 0.3/0.7; the loss rate and detection are known
within each scenario. No per-dataset true parameters or parent are passed
to the fitter.

| Scenario | Duplication grid | Loss | Ne grid | Detection |
| --- | --- | --- | --- | --- |
| baseline | 0.03 / 0.08 | 0.03 | 0.3 / 1 | 0.9 |
| turnover | 0.08 / 0.16 | 0.12 | 0.3 / 1 | 0.9 |
| missing-ils | 0.03 / 0.08 | 0.03 | 1 / 3 | 0.6 |

Each scenario has seven truth cases: all four distinct on-grid null
parameter points; an off-grid null with midpoint duplication and maximum
grid Ne; and allopolyploid H2 A/B with that same midpoint duplication,
maximum Ne, and age 0.5. The last three cases are off-grid. Off-grid null
results are misspecification diagnostics, not tests of the finite-grid
validity claim.

Five independent replicates per truth case, 30 families each: **105 datasets,
75 null and 30 allopolyploid**, 3,150 selected families. Independent generation
uses the existing global SSA, full msprime genealogies with joint daughter
rejection, and detection. Full histories, acceptance audits, 600-site JC69
alignments, true trees, and NJ/midpoint estimated trees are retained.
Only true rooted topologies are calibrated in this phase. Estimated trees
are preserved, not substituted into a true-tree event null distribution.

Each scenario builds 20 fresh integration banks with **100,000 raw draws
per bank**, using `ancestral-stratified`, simultaneous 99% score-bank bounds,
and the unchanged node/state/selection caps. No zero-hit smoothing or
family/history truncation. Failed banks prevent scoring that scenario;
every affected planned dataset is recorded as failed.

Seeds: integration 20261030 + scenario index; evaluation SeedSequence
`[20261033, scenario_index, case_index, replicate_index]`; calibration seed
is the uint32 SeedSequence output from the same index tuple headed by
20261034. The bank seed remains separately recorded. Every dataset uses
fresh calibration draws rather than a shared null pool. Worker count
changes scheduling only.

For every observed dataset, all four null grid points generate 99 fresh
replicates of 30 selected families. Each replicate repeats the full
null/alternative parameter and parent search with the independently frozen
integration banks. Point P-values use `(1 + exceedances)/(99 + 1)`, including
ties. Supremum P-values and each MC-overlap endpoint are the separate maxima
over the four generating null points. The numerical fitted-null result
from these **same draws** is recorded as a paired comparison; it is not a
replay of the old plug-in seed stream.

Decision threshold is fixed at 0.05: an event is reported only when the
supremum **upper score-MC P bound** is <=0.05. Point-P decisions, fitted-null
decisions, least-favorable generating null points, parent recovery, and
estimated-tree accuracy are separate descriptive outputs.

## Denominators and Limits

All generation, integration, support, and calibration failures remain in
the planned counts. Reports distinguish completed and failed counts.
The operational report rate counts failures as no report, not as a
scientific negative. Failure-aware rate ranges allow either outcome for
each failed analysis. Their pointwise 95% binomial envelope takes the lower
endpoint with zero failed positives and the upper endpoint with all failed
positives. The intervals are descriptive, conditional on shared frozen
banks, not simultaneous confidence bounds or a general biological guarantee.
Parent counts use completed alternative analyses, with failures reported
alongside rather than removed from the planned denominator.

Finite-grid supremum rank-test validity requires a true null point on the
supplied grid, identical generative observation/selection laws for data
and null replicates, independently frozen scoring banks, and complete IID
replicates. It does not cover continuous off-grid nuisance values,
misspecified detection, empirical gene-tree errors, or broader event modes.
Five datasets per truth point cannot establish a uniform empirical 5%
false-positive guarantee; report uncertainty and failures without tuning
the protocol from results.

The statistical construction is a finite-grid specialization of maximized
Monte Carlo testing, not a new proof for the unrestricted biological model.
See [Dufour (2006)](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf).

## Run

Use the existing isolated research environment with msprime and Biopython.
The output directory must be fresh; no automatic retries with larger
resources, changed thresholds, or substituted datasets are performed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/mul-locus/calibration.py \
  --output /tmp/new-locus-grid-calibration --cpus 4
```

`protocol.json` is written before simulation. It records the complete
configuration, cases, seeds, arguments, dependencies, and source hashes.
Scenario banks, every dataset's observations/fit/calibration or failure,
`summary.tsv`, and failure-aware `rates.tsv` are retained. The program
finishes remaining independent scenarios and exits nonzero if any analysis
fails. Resource/dependency or incomplete-study failures must not be reported
as a successful accuracy study.

## Post-Evaluation Runner Audit

The [subsequent audit](../../docs/validation/MUL_LOCUS_GRID_AUDIT.md)
preserves the frozen results above. Future runs use fixed named scenario
indices (baseline 0, turnover 1, missing-ils 2), including for subsets or
reordered requests, and reject calibration-seed collisions before writing
the protocol or simulating. The original default full-study streams remain
unchanged; scientific settings and thresholds were not revised.

Bank failures retain a full trial record and failure file for every planned
dataset. A parallel-worker failure retains partial banks, recovers valid
saved trial summaries/failures and records unresolved trials as
`worker-failed`, with their original metadata and seeds. Completed workers
save `summary.json` before returning their row. Unreadable/incomplete records
remain preserved; a companion `worker-failure.json` records failure without
overwriting existing evidence. No resampling or automatic retry is used.
Main-process termination or an unavailable output filesystem can still leave
an incomplete study; such artifacts are not a successful validation gate.
