# Native GY94/profile validation with external sequences

This study tests the model/interval settings used by GeneGalleon, without
adding species-age uncertainty propagation or changing the estimator. It is
not a demonstration of universal 95% coverage. The public sequence profile
method remains experimental.

## Design and evidence

- [Primary protocol](protocol.json): 200 root-duplication, 100 high-rate-SD,
  and 100 internal-duplication families.
- [Stress protocol](stress-protocol.json): 50 families each with longer
  alignments, more tips, loss plus a copy-wide rate shift, one wrong tip mapping,
  and low rate SD.
- [Pilot protocol](pilot-protocol.json): 12 separate development families,
  excluded from the 650-family validation study.
- `results/summary.tsv` and `results/coverage.png`: audited numerical results.
- Per-family JSONL and metadata retain failures, source/input hashes, seeds,
  actual estimator selection and unavailable-interval reasons.

The external [AliSim sequence generator](https://iqtree.github.io/doc/AliSim)
uses `GY{0.5,2}+FQ+G4{0.7}`: omega 0.5, kappa 2, uniform codon frequencies,
four gamma categories and shape 0.7. AliSim's `--length` is the nucleotide
length for codon input; the runner verifies the returned codon count and absence
of standard-code stops. A preliminary plumbing run caught the length convention
before inference; those development cases are not validation observations.

The chronology fixtures do not invoke NWKIT inference. Independent lognormal
edge rates have median 0.01 expected nucleotide changes per codon per time unit.
The species-root age is 10; duplication truth is 20 (root) or 7.5 (internal).
AliSim simulates the sequences using those branches. Native inference refits
branch lengths, frequencies, kappa, omega and gamma shape from each alignment.
The true substitution lengths are starting values, not fixed observations in
the sequence analysis. True rooted topology is supplied, with LCA reconciliation;
this is not an end-to-end test of gene-tree inference or GeneRax searching.

Each family contributes one prespecified duplication. Fixed species-node ages
do not count as successful age estimates. The wrong-map condition assesses the
original root duplication only; it does not validate new or missing inferred
events. Loss and copy-wide rate variation are combined in one stress condition
and cannot be separated causally from that condition alone.

Two runs record point-only success separately from profile success. Auto can
refit point ages during a profile calculation; both point estimates and actual
estimators are retained. Bias and RMSE use the final profile point when returned,
otherwise the completed point-only result, and are conditional on point success.
Timeouts and missing intervals remain in all-trial denominators. Calibration-
limited returned intervals remain visible and are included in returned-interval
coverage. A strict-clock boundary can yield a completed point but no interval.

The 95% Wilson intervals describe Monte Carlo uncertainty in the observed
fractions; they are not gene-age intervals. A confidence interval including
0.95 does not establish calibration. Neither a wide interval covering truth
nor an unavailable interval is silently replaced by a narrow interval.

## Completed pre-fix results

| Condition | Trials | Points | Intervals | Covered / intervals | Correct return / trials |
|---|---:|---:|---:|---:|---:|
| root | 200 | 200 | 189 | 176/189 (93.1%) | 88.0% |
| high_sd | 100 | 100 | 98 | 85/98 (86.7%) | 85.0% |
| internal | 100 | 99 | 99 | 95/99 (96.0%) | 95.0% |
| long_alignment | 50 | 50 | 49 | 41/49 (83.7%) | 82.0% |
| larger_family | 50 | 50 | 48 | 48/48 (100.0%) | 96.0% |
| loss_rate_shift | 50 | 50 | 50 | 49/50 (98.0%) | 98.0% |
| wrong_mapping | 50 | 50 | 43 | 41/43 (95.3%) | 82.0% |
| low_sd | 50 | 50 | 37 | 35/37 (94.6%) | 70.0% |

All 650 input/result records passed the independent audit. High rate variation
and longer alignments show undercoverage even among returned intervals. Low
rate variation returned intervals in only 37/50 trials (13 strict-clock limits).
The larger-family condition retained two 180-second profile timeouts. Wrong
mapping returned 43 calibration-limited intervals and seven strict-clock limits;
its 95.3% conditional coverage is only an 82% correct-return fraction over all
trials and does not validate the incorrect reconciliation. Bias, RMSE, interval
width and all status counts are in [the audited summary](results/summary.json).

The numerical initialization repair below is separate from these pre-fix
calibration results. Undercoverage remains unresolved; no empirical correction
factor or new default interval method was introduced.

## Reproduce

Run in a GeneGalleon container with IQ-TREE 3 and the current NWKIT checkout
on `PYTHONPATH`. Use new output directories; existing evidence is not overwritten.

```bash
export PYTHONPATH=/nwkit
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python tools/validate_radte_default_profile.py \
  --protocol examples/radte/default-profile-validation/protocol.json \
  --output /results/primary --workers 6 --timeout 180
python tools/validate_radte_default_profile.py \
  --protocol examples/radte/default-profile-validation/stress-protocol.json \
  --output /results/stress --workers 3 --timeout 180
python tools/report_radte_default_profile.py \
  --results /results/primary /results/stress --output /results/report
```

The audit checks complete, unique family sets, per-family JSON/JSONL agreement,
input hashes, truth ages against the generating chronology, coverage indicators,
point-error metrics, denominators and Wilson limits. Each raw
family directory retains the input files, actual command arrays, process logs,
point/profile manifests, complete output tables and the result JSON. Hashes
identify the tested sources rather than assuming a package version identifies
uncommitted source changes.

The validation runtime is `local/genegalleon:standard-iqtree-dev`, arm64 image
`sha256:deca052d1d71b2c9530c92a33b053be06a1d6e86088325de3ce0ed5f7eee7912`,
with local NWKIT sources mounted on `PYTHONPATH`. This is Docker validation with
a source overlay, not verification of the image's installed NWKIT snapshot or
SIF compatibility. Image/revision identities here are experiment records, not
new dependency defaults. The tests started from NWKIT `e2465b048f7529b7dad84cb6722278406bb49bbe`
and GeneGalleon `021c5b62c551d073016ce14676a56ec9768f75c8`.

## Interpretation limits and diagnostics

The data do not include ILS, uncertainty in species ages, empirical topology
errors, gene conversion, indels or protein models. The true topology and simple
stationary generating substitution model make this more favorable than many
real families. No best-looking interval method is selected after inspecting
coverage, and no scale factor is fitted to these outcomes.

`tools/diagnose_radte_default_profile.py` replays the first 50 primary families
with the generating rate SD supplied. This is an oracle diagnostic, not a
practical method or 50 additional independent families. Supplying SD can also
change auto's estimator selection, so this comparison alone cannot attribute
every change to variance-estimation uncertainty.

In the [paired diagnostic](oracle/summary.json), original fits covered 43/48
returned intervals (43/50 all trials); supplying the generating SD covered
47/50, with all 50 intervals returned. Median widths were 13.60 and 13.17 age
units respectively. This supports sensitivity to rate-SD treatment, not a claim
that uniformly widening intervals solves the problem. The small paired sample
and changed auto path do not identify a unique cause or establish calibration.

## Numerical initialization repair and separate checks

Primary internal family 99 failed before any profile interval could be computed.
A branch-only warm start collapsed one duration to the chronology boundary,
producing an exact-sequence age gradient about 2.96e9. The fix resets that
near-minimum-duration warm start to the existing chronology interior and records
`sequence_initial_ages_reset_from_duration_boundary`. The objective, rate-SD fit,
constraints and final minimum duration are unchanged. The retained fixture
`tests/data/radte-initial-boundary/` reproduces the original failure; independent
restarts of the unchanged objective agree after repair. The repaired age is
9.3410, interval [7.3875, 10], truth 7.5. The original failure is still counted
in the pre-fix study.

The separate [post-fix check](postfix/summary.json) uses 20 new internal families:
20 point estimates and intervals returned, 19 covering truth. None triggered
the reset, so these are general regression observations rather than 20 examples
of the boundary repair. They are excluded from the 650-family pre-fix study.
The numerical/profile regression suite passed 119 tests; all nine GeneGalleon
native/IQ-TREE dating integration scenarios passed against the repaired source.

The plot/report changes make the exploratory/conditional scope visible. Statistical calibration improvements belong in NWKIT;
GeneGalleon must not inflate, clip, or replace intervals to conceal a failure.

## Commit review checks

The final source snapshot passed 354 dating/plot/CLI/example tests, full-source
Ruff/format and mypy checks, and the maintainability complexity limits. The
GeneGalleon snapshot passed 256 integration/summary/static tests and five smoke
tests. Distribution contents and reproducibility checks passed. These scoped
local checks do not represent the complete hosted CI matrix or SIF execution.

Nonfinite returned ages are now recorded as failures without terminating the
validation study; an invalid profile result retains a valid point-only fit.
New regression tests cover both runner paths and the oracle path. The 650-family
results remain unchanged. The boundary initializer was subsequently extracted
into a small helper with the same decision rule to satisfy the complexity limit;
the archived repair patch/source hashes describe the original experimental fix.
