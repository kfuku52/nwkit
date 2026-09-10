# OU shift simulation validation

**Historical pBIC results:** the pBIC numbers below were generated with the old
kfl1ou 3.0.9 scoring implementation. The [pBIC correction](SHIFT_PBIC.md) changes
its coefficient determinant and can change selections; these results do not
validate the corrected criterion. Original evidence is retained. BIC is
unaffected by the correction.

This is a small, reproducible statistical pilot, not a calibration guarantee.
It tests whether `shift` recovers independently generated changes and shared
optima, how often it selects shifts when none exist, and whether bootstrap
support and failures are reported honestly. See [SHIFT.md](SHIFT.md) for the
model and command contract.

## Pilot results (2026-09-09)

**Assessment: share with caveats.** The software pipeline completed the eight-tip
pilot, but selection was unreliable in several known-truth settings. Successful
execution and bootstrap completion do not establish statistical calibration.
Model averaging should not be used to conceal this uncertainty.

The following counts pool the four root/SE settings, with 20 datasets per
scenario. Recovery compares the complete partition of extant tips, not just
whether at least one true shift was found. All 320 point fits
completed (80 datasets × two criteria × two point-estimation modes).

| Eight-tip metric | BIC | pBIC (command default) |
| --- | ---: | ---: |
| Any shift selected under no-shift truth | 18/20 | 9/20 |
| Single-shift ancestry partition recovered | 2/20 | 12/20 |
| Distinct two-shift ancestry partition recovered | 18/20 | 7/20 |
| Convergent two-shift ancestry partition recovered | 6/20 | 3/20 |
| Convergent shared-optimum partition recovered after convergence | 6/20 | 3/20 |
| Any merge in convergent cases | 11/20 | 4/20 |

pBIC reduced false shift selections in this paired pilot, but also missed more
of the distinct and convergent two-shift patterns. These observations do not
justify declaring either criterion generally superior or assigning a nominal
false-positive rate. Five replicates per cell leave large Monte Carlo uncertainty.
The [full table](examples/shift/validation/table.md) also reports expected-tip
RMSE; the [cell summaries](examples/shift/validation/cell-summary.json) retain
root and SE differences and Wilson intervals.

For `null`, `single`, and `distinct` combined, BIC produced 12/60 merge events
and pBIC 0/60. Those are pipeline events, not all proven false equality claims:
most arise after incorrectly selecting the shift locations. Within the subset
with exactly correct shift edges, neither criterion merged these nonconvergent
truths (BIC 0/22, pBIC 0/30). These selected subsets are too small and conditional
to establish a low false-convergence rate.

All 160 outer bootstrap runs returned results, and all 480 inner refits
succeeded. With only three inner draws per fit, support is intentionally coarse.
The mean bootstrap frequency of the **true shared-optimum partition** was:

| Truth | BIC | pBIC |
| --- | ---: | ---: |
| No shift | 0.017 | 0.367 |
| Single shift | 0.050 | 0.550 |
| Two distinct shifts | 0.783 | 0.317 |
| Two convergent shifts | 0.183 | 0.183 |

These are means over 20 available outer fits per table entry, not
posterior probabilities, confidence coverage, or the probability that the
observed selected model is correct. No missing or failed fit was dropped.

### Sixteen-tip extension

All 192 additional point fits completed. pBIC selected a shift in 8/24 no-shift
truth datasets. Splitting these independent null draws by their inactive effect
label gives 6/12 and 2/12; that difference is Monte Carlo variation, not a signal
strength effect on a no-shift model.

| Sixteen-tip pBIC recovery metric | Effect 0.5 | Effect 2 |
| --- | ---: | ---: |
| Single-shift ancestry partition | 1/12 | 11/12 |
| Distinct two-shift ancestry partition | 0/12 | 7/12 |
| Convergent two-shift ancestry partition | 0/12 | 2/12 |
| Convergent shared-optimum partition after convergence | 0/12 | 4/12 |

The stronger single shifts were often recovered, while convergence remained
hard even at the stronger effect. Shared-optimum recovery can succeed using a
different ancestral shift partition, hence 4/12 shared recoveries versus 2/12
ancestry recoveries. The [full table](examples/shift/validation-16tip/table.md)
and [cell summaries](examples/shift/validation-16tip/cell-summary.json) retain
all outcomes. Only three replicates per root/SE/effect cell were used.

Across both studies, all 672 outer runs completed on 176 distinct simulated
datasets; the eight-tip datasets were reused across criteria and methods.
No outer run or available inner bootstrap refit failed. This measures returned
results under this grid, not guaranteed numerical reliability or inferential
accuracy in other datasets.

### Implication for the next implementation

Keep this interface experimental. Increase outer replication on an independent,
predeclared grid before interpreting selection or bootstrap support as reliable
error control. Include larger and imbalanced trees, matched total tree heights,
weaker/stronger OU pull, heterogeneous SEs and nested shifts. Bootstrap requires
substantially more inner draws before support calibration can be assessed.

Backward convergence only merges regimes found by the initial unconstrained
search; it cannot restore a shift that discovery excluded. The
[small-tree joint reference](SHIFT_JOINT.md) now compares both procedures using
held-out seeds: it found lower BIC scores in 2/12 datasets without increasing
shared-partition recovery. It also identified an incompatible pBIC penalty
between equivalent free and constrained representations.
Model averaging should follow that evaluation and preserve each candidate's
shared-optimum constraints.

## Design and independent truth

`tools/shift_simulation_cases.py` generates an eight-tip balanced ultrametric
tree with unit branches. Each edge is simulated directly as

```
X_child = exp(-alpha*t) * X_parent
          + (1-exp(-alpha*t)) * theta_edge + Normal(0, q)
q = sigma2 * (1-exp(-2*alpha*t)) / (2*alpha)
```

The generator does not call `kfl1ou`, a fitted model, or the NWKIT Gaussian
simulator. It uses alpha = 0.7, sigma2 = 0.25, baseline optimum/root mean = 0,
and effect magnitude = 2. A fixed root is exactly zero; a random root is drawn
from Normal(0, sigma2/(2*alpha)). Independent measurement noise with known
SE = 0 or 0.2 is added to every tip and supplied to the fit.

The four scenarios are `null` (no shifts), `single` (one quarter-clade shifts
to +2), `distinct` (disconnected quarter-clades shift to +2 and -2), and
`convergent` (those two clades independently shift to +2). Five independent
replicates of each scenario × root treatment × SE cell give 80 datasets.
BIC and the command's default pBIC use exactly the same datasets, providing a
paired comparison. Search is exhaustive over at most two shifts. The generator
records tip expectations after OU attenuation, rather than treating tip optima
as expected observations.

The master seed is 20260909. Each cell and replicate has its own deterministic
NumPy SeedSequence stream. The manifest records the ordered grid, versions and
SHA-256 hashes of the runner, generator and shift implementation. Reordering or
changing the grid can change cell indices and therefore streams; use archived
inputs for exact comparisons across different grids. The two tree-size grids
reuse some random seeds, so their observations must not be pooled as mutually
independent trials. Comparisons between those grids here are descriptive; use
a different master seed for a genuinely independent follow-up study.

Each dataset is run three times: unconstrained shift discovery, discovery plus
backward convergence, and the complete convergence pipeline with three
parametric bootstrap replicates. Bootstrap uses its fitted model's simulations,
whereas the outer datasets use independent known truth. Separate point fits
ensure an all-failed bootstrap does not erase a valid point estimate.

A secondary pBIC grid uses 16 tips, effect magnitudes 0.5 and 2, both roots,
both SE levels and three replicates per cell: 96 further datasets. It evaluates
point estimates with and without convergence; bootstrap is disabled in this
extension. Its inputs and results are archived separately under
[validation-16tip](examples/shift/validation-16tip/records.json). Effect magnitude
has no effect in a `null` generating model; the two null cells use independent
seeds and are kept separately in the summaries. With unit branches, total tree
height is 3 for eight tips and 4 for sixteen tips; the shifted clades also have
more time to approach their optima. This extension does not isolate tip count
while holding tree height or the induced tip-mean effect fixed.

## Metrics and denominators

- `any_shift` under `null` is a false shift selection. It is not a hypothesis
  test with a nominal 5% level; information criteria do not specify that level.
- `exact_edges` compares selected branch IDs with generating shift locations.
  `ancestry_recovered` compares label-independent partitions of extant tips by
  their most recent shift. Equivalent representations can disagree on edges.
- `shared_recovered` compares the entire inferred partition of tip optima with
  the generating equality groups. This is stricter than observing any merge.
  It requires positive fitted alpha; at the BM boundary no OU optima are
  identifiable. In unconstrained fits, different shift regimes remain separate
  labels. All fits in this pilot had positive alpha.
- `any_merge` records whether the backward procedure joined regimes. A merge
  can involve spurious discoveries, so this alone is **not** a false-convergence
  rate. `merges_given_correct_edges` additionally conditions on recovering the
  exact generating shift configuration; under `null`, `single`, or `distinct`,
  a merge in that subset is unwarranted.
- Tip-mean RMSE compares fitted expected tip values with the generating
  expectations in original trait units; it does not compare to noisy data.
- Conditional rates use completed fits. The `*_returned` recovery rates use
  **all attempted outer fits**, counting failures as no correct result returned.
  Every failed run retains its case, mode, seed and error text.
- Cell-level binary rates include Wilson 95% Monte Carlo intervals. With five
  replicates, even 0/5 has an upper bound of approximately 43%; 5/5 has a lower
  bound near 57%. The pooled table is descriptive across four root/SE settings.
  Use the individual cells for uncertainty; pooled binomial intervals in the
  machine-readable summary are only approximations for this heterogeneous mix.
- Bootstrap truth-partition frequencies divide by successful inner refits and
  are averaged over outer fits with an available bootstrap result. They are
  selection frequencies, not posterior probabilities or confidence coverage.
  Inner counts cover returned bootstrap results only. A failed outer bootstrap
  has unknown inner counts and is reported separately, never as zero failures.
  Three inner replicates provide only a coarse diagnostic.

The archived [inputs](examples/shift/validation/inputs.json),
[per-run records](examples/shift/validation/records.json), and
[cell summaries](examples/shift/validation/cell-summary.json) permit independent
recalculation. The CSV is an alternative export: when reading it with pandas,
use `keep_default_na=False` so the scenario name `null` stays a string.

## Reproduction

Install NWKIT's numerical dependencies and R with kfl1ou >= 3.0.9. From the
checkout, choose new output directories (existing directories are refused):

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/validate_shift_simulations.py \
  --output /tmp/shift-bic --criterion BIC --replicates 5 --bootstrap 3 \
  --rscript /path/to/Rscript
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/validate_shift_simulations.py \
  --output /tmp/shift-pbic --criterion pBIC --replicates 5 --bootstrap 3 \
  --rscript /path/to/Rscript
PYTHONPATH=. python tools/summarize_shift_validation.py \
  /tmp/shift-bic /tmp/shift-pbic --output /tmp/shift-evidence
```

Raw directories contain every Newick tree, trait table, truth record, fit JSON,
regime map, execution log and failure record. The compact evidence exporter
retains all input truth, per-run metrics and checksum-verified source snapshots
without duplicating large fit JSONs.
It requires completed summaries, rejects duplicate runs and checks paired input
identity. No fit is retried or omitted because of its result.

The secondary grid is reproduced with:

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/validate_shift_simulations.py \
  --output /tmp/shift-16tip --criterion pBIC --tips 16 --effects 0.5 2 \
  --replicates 3 --bootstrap 0 --rscript /path/to/Rscript
PYTHONPATH=. python tools/summarize_shift_validation.py \
  /tmp/shift-16tip --output /tmp/shift-16tip-evidence
```

For a larger study use `--tips 8 16 32`, `--effects 0.5 1 2`, larger
`--replicates` and `--bootstrap`, and preserve all cell definitions and failures.
Exhaustive search becomes expensive with tree size. This pilot does not cover
imbalanced trees, nested shifts, misspecified root/error models, missing tips,
or more than two shifts; those need separate studies.

## Verification

Tests compare the generator's exact means and covariance against the independent
Gaussian reference for all four scenarios and both roots, and check empirical
innovation moments using 2,000 draws per root. They also verify label-invariant
partitions, distinct ancestry/shared-optimum metrics, invalid inputs, and failure
denominators. Related integration tests execute the real R backend, including
convergence constraints, known SEs, bootstrap and same-parameter likelihoods.


## Evidence audit (2026-09-10)

The exporter now verifies the declared grid and seed against regenerated truth,
checks the actual trait table and tree, and links each trial record to its saved
result and fitted model. It checks model settings, tip coverage, branch topology,
observations and SEs, then recalculates selection metrics and cell summaries.
Floating-point comparisons allow rounding across numerical environments; changed
labels, missing keys, invalid coordinates and stale summaries are rejected.
This is a consistency audit, not authentication against coordinated alteration.

All 672 existing fitted models and 176 inputs passed this audit. Re-exporting
changed none of the trial metrics or aggregate results above. The original
simulation source snapshots are retained alongside the current exporter and
audit code. No statistical replicate was rerun or selectively replaced.

Evidence is assembled in a temporary sibling directory and published after all
validation and serialization succeeds. Existing output bundles are refused;
a failed write leaves no partial published bundle. Simulation grids reject
duplicate axis entries and out-of-range bootstrap counts before creating a run,
so configuration mistakes do not become misleading refit-failure statistics.

Regression tests inject changed truth, trait rows, condition labels, missing
error messages, stale summaries, contradictory models, cyclic branch parents
and disk-write failures. The focused suite, including real kfl1ou integration,
passed 161 tests; lint, formatting and type checks also passed.
