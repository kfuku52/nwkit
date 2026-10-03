# Conditional Integration Probe

This is a new exploratory protocol, not a retuning or replacement of the
[finite-grid study](CALIBRATION.md). Its historical banks, datasets, failures
and decisions remain unchanged. Do not interpret this low-budget probe as
empirical 5% error control or a power study.

## Fixed First Probe

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/integration.py \
  --output /tmp/nwkit-locus-integration-probe-20261002 \
  --source /tmp/nwkit-locus-grid-calibration-20261002
```

Before diagnosis, simulation or scoring, the script saves arguments, model
configurations, source hashes and runtime information to `protocol.json`.
The output directory must not exist. Fixed defaults are 1,000 raw DL draws
per candidate/grid, 10 families per dataset, one dataset at each of the seven
truth points in each of the three original scenarios, and 19 independent
null datasets per null grid point. Thus there are 21 paired datasets and
63 planned method analyses. Rates/Ne, detection, species ages, parent scope
and observed-family selection (2-4 tips) are unchanged. Data use the
independent global SSA/msprime sampler, including during null calibration.
No estimated-tree error process is included in this probe.

Data, bank and calibration base seeds are respectively 20261101, 20261102
and 20261103. Canonical scenario indices are baseline 0, turnover 1,
missing-ils 2, regardless of request order. DL/gene/detection integration
namespaces are 3/4/5; independent data/calibration namespaces are 6/7.
Calibration seeds retain four uint32 words as a Python integer, not a single
uint32 compression. Every method shares the same complete observed and null
data streams. No failure triggers a retry or another seed.

## Compared Estimators

- `histogram`: original ancestral-stratified estimator and Clopper-Pearson
  bounds, with the new explicitly paired RNG streams (not the historical
  bank RNG streams).
- `detection`: sample the complete daughter-bounded genealogy and integrate
  all detection outcomes using a tree DP.
- `hybrid`: for histories with at most **four extant hidden loci**, integrate
  both the full genealogy and detection; for larger histories, use the same
  genealogy draw/detection integration as `detection`. This rule is fixed
  before data and independent of observed patterns or work-limit failures.

The full hidden locus tree, including undetected copies, is used for every
daughter constraint. The observed four-tip limit is not a hidden-history
limit. No large history is discarded. Exceeding the 100,000-work cap aborts
the paired bank, retains the failure and marks all planned analyses as
bank-failed. A failure is never converted into a negative event or a
partial-bank result. Genealogy/detection DP mass failures also abort.

The first execution exhausted the detection-DP cap in all three scenarios
before scoring any dataset. All 63 planned analyses remain bank-failed in
`/tmp/nwkit-locus-integration-probe-20261002`, with its original source archive.
The follow-up changes only computation: detection outcomes above four
observed tips are exactly aggregated into an overflow probability at every
gene-tree node. Adding the other child cannot reduce a completed child's
detected-tip count, so these states can never become selected later. No
hidden gene tip or daughter constraint is removed. Independent tests compare
the aggregated distribution with complete detection-mask enumeration.
All original scientific settings, seeds, budgets and decisions stay fixed;
the separate output is `/tmp/nwkit-locus-integration-overflow-probe-20261003`.

That follow-up constructed all 20 baseline banks, but stopped during the
second dataset's calibration because the original reference helper always
performed sequence/NJ inference, even when returning true trees. Biopython
raised `UnboundLocalError` while midpoint-rooting a zero-distance estimated
tree. That incomplete execution and its original source archive remain
preserved; it is not a successful probe or a completed 63-analysis study.

The final operational follow-up uses the explicit `reference_true_selected`
path, which returns the detected true genealogy before mutation/NJ steps.
The original reference-family/default estimated-tree behavior and RNG stream
are unchanged. The true-only observation law is unchanged, but removing the
unused draws changes the new study's realized random stream; this is not
an exact replay of the incomplete follow-up. Unexpected per-bank/generation/
calibration exceptions are now retained too. Scientific parameters, draw
budgets, base seeds and decision thresholds remain fixed, with separate output
`/tmp/nwkit-locus-integration-true-probe-20261003`.

The exact forest DP combines count transition probabilities with uniform
random-pair merger paths, and uses the full-history daughter normalizers in
log space. Independent small-history tests instead exponentiate an all-pair
forest CTMC and enumerate every detection mask, including nested daughters,
losses, repeated species colors and undetected tips.

Each conditional contribution is a bounded variable in [0,1]. Raw pattern
and selection estimates are unbiased within each IID stratum; the final
ratio of estimated pattern to selection mass is not claimed unbiased.
Prior weighting and selection normalization are identical across methods.
Conditional expectation reduces raw-contribution variance, but does not
guarantee narrower finite-sample intervals or solve unsampled DL histories.

## MC Bounds And Decisions

Weighted contributions do **not** use Clopper-Pearson. With sample variance
`s^2`, n >=2 and per-variable error budget alpha, the prototype uses

```text
radius = sqrt(2*s^2*log(4/alpha)/n) + 7*log(4/alpha)/(3*(n-1))
interval = [max(0, mean-radius), min(1, mean+radius)]
```

This applies [Maurer and Pontil (2009), Theorem 4](https://www.cs.mcgill.ca/~colt2009/papers/012.pdf)
to both tails. The same existing simultaneous budget, `(1-confidence) /
(3 * number_of_banks * category_bound)`, covers at most two active strata:
two intervals per possible pattern plus two selection-denominator intervals
is at most three times the observation-universe bound. Coverage is per
method, not a simultaneous claim across the three compared estimators.
Intervals are prior weighted and ratio bounded before taking logs. A zero
estimate/lower bound remains zero; there are no pseudocounts or epsilon
log-likelihood repairs. A zero lower selection endpoint yields upper bound 1.

All methods use the same full-search finite-grid supremum calibration.
An event requires maximum upper score-MC P <=0.05, unchanged. Failure rates,
point and upper-bound decisions, parent recovery, probability interval
widths and complete planned denominators are retained. Contrast infinities
are JSON nulls, not finite evidence. Per-method elapsed calibration and
combined paired-bank time are operational costs, not a speed comparison;
methods can stop at different failure points.

## Historical Diagnosis

`diagnosis.json` reconstructs all nine historical failing null family sets
from their recorded seeds, generating grid and replicate. It records every
observed pattern against every bank, and separately totals hits across
all topologies with the same species copy-count vector. Zero copy-vector
hits identify **DL-history or detection support**, not necessarily a
genealogy problem. Without the historical raw hidden histories, those two
causes cannot be separated further. Positive copy-vector hits but zero
pattern hits identify a topology-support boundary in that bank. These are
sample-support descriptions, not proof of zero biological probability.
Different banks are never pooled to fabricate a finite candidate score.

Results and adoption judgment are recorded separately in
[the validation record](../../docs/validation/MUL_LOCUS_INTEGRATION_PROBE.md).
