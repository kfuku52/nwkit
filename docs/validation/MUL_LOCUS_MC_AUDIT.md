# Locus MC Numerical and Reproducibility Audit

Executed 2026-10-02 on the uncommitted locus/MSC implementation based on
`0a4cea2881cf0559f124f68b03d4b0203759bb67` (checkout version `0.43.37`).
This is a subsequent audit, not a replacement for the historical
[independent study record](MUL_LOCUS_MC_VALIDATION.md). That record, its
study artifacts, and its recorded source hashes were not rewritten.

## Findings and Repairs

1. **Scientific JSON values could be misinterpreted.** Python treated
   booleans as numerical rates, Ne, ages, stem lengths, or detection
   probabilities. Other malformed JSON structures produced incidental
   `TypeError`/`AttributeError`, and very large integers caused conversion
   overflows. Validation now requires a model object, detection object,
   nonempty array of parameter objects, and finite real scientific numbers
   excluding booleans/strings. The combined duplication/loss rate and
   diploid `2*Ne` scale are checked too, including integer overflow cases.
2. **The MC precision guard came after an unbounded universe calculation.**
   `max_observed_tips=1000` raised `OverflowError` while converting the
   enormous category bound to float; larger settings could spend excessive
   work constructing it before rejecting precision. The monotonically
   increasing bound is now checked as it grows. Tests at 100, 1000, and
   one billion tips verify rejection before bank sampling and before
   creating output files. This changes neither accepted small-family
   bounds nor the declared scientific family-size selection.
3. **The result artifact did not retain its full scored model.** Candidate
   population chronograms, population-tip species labels, supplied H1/H2,
   and ordered observations were absent. These are now additive JSON
   fields. For both integration methods, tests reconstruct the population
   trees and banks from JSON, regenerate their histograms/strata with the
   same seeds, and reproduce all saved scores exactly.

The initial pre-repair regression run demonstrated **13 failures and two
passes**. It did not execute the billion-tip case on the defective code.
Further boundary inspection reproduced the integer combined-rate/`2*Ne`
overflow cases before their repair. Invalid configurations fail explicitly;
there is no replacement value, pseudocount, family dropping, or weaker cap.

Two interpretation risks were also made explicit:

- Stratified report `hits/samples` is an aggregate raw proposal frequency,
  not the selected-family probability. Reports now identify the integration
  method; the guide directs readers to the prior-weighted,
  selection-normalized `probability_estimate` and the saved strata.
- Bootstrap MC-overlap bounds condition on the numerical fitted null used
  to generate replicates. They do not envelope uncertainty about which
  null grid point should generate those replicates, or finite-bootstrap
  sampling error. This restriction is stated in both JSON and the guide.
  No uniform composite-null error-control claim was added.

## Independent Numerical Check

A new fixed-history reference includes two nested SSD daughter bounds,
two losses, five surviving loci (including repeated species labels),
species-dependent non-detection, and selection of two to four observed
tips. The reference tracks every pair merger using the independent full
forest CTMC and matrix exponential. It imposes daughter bounds jointly
without the production count DP or local normalizations.

The exact joint bound probability agrees with the product of the production
DP's conditional normalizers to absolute tolerance `2e-13`. Complete
gene-topology probabilities sum to one. Independent enumeration of every
detection subset, followed by tree pruning, gives the expected colored
observation distribution. In 12,000 full-genealogy draws, both selection
frequency and all conditional observation probabilities lie within the
predeclared per-comparison binomial intervals (`alpha=0.0001`). No
tolerances, rates, or seeds were changed to obtain agreement.

Existing zero-DL MSC/full-forest checks for both parents, daughter-bounded
checks, birth/death closed forms, weighted-stratum normalization,
full-search bootstrap reselection, serial/process equality, input/output
alias protection, and bundle rollback checks remain active.

## Executed Checks

Runtime: `/tmp/nwkit-wgd-dev-20261001/bin/python`, CPython 3.12.14 on
macOS 27 ARM64; NumPy 2.5.3, SciPy 1.18.1. The independent study replay
used msprime 1.4.4, tskit 1.0.3, and Biopython 1.88. Compiled-library
imports and `pip check` passed before reuse.

```sh
PY=/tmp/nwkit-wgd-dev-20261001/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  "$PY" tools/check.py quick -- \
  tests/test_mul_locus.py tests/test_mul_coalescent.py \
  tests/test_mul_msc.py tests/test_mul_msc_fit.py \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli.py tests/test_cli_contracts.py \
  tests/test_interface_conventions.py tests/test_output_transaction.py \
  tests/test_rooting_state.py tests/test_provenance.py \
  tests/test_util_tree_io.py tests/test_tree_outputs.py -x -rs
"$PY" tools/check.py test -- -m slow \
  tests/test_mul_locus.py tests/test_mul_msc_fit.py tests/test_mul_msc.py -rs
"$PY" tools/check.py quick -- tests/test_mul_locus.py -k saved_population -rs
"$PY" tools/check_maintainability.py
"$PY" -m nwkit mul-reconcile --help
git diff --check
```

- Broad quick lane: **524 passed**, 13 slow tests deselected, no skips.
  Ruff lint/format checked 562 files and mypy checked 249 source files.
  All four original-GRAMPA reference comparisons executed.
- Explicit slow lane: **13 passed**, 149 other tests deselected, no skips.
  Together these cover 537 distinct tests, not the entire NWKIT suite.
- After strengthening JSON round-trip cases with histogram/stratum
  regeneration, both affected cases passed again with static checks.
- Maintainability: 3,605 functions, average complexity 6.92, maximum 50;
  all hard limits passed. Existing baseline-comparison warnings remain.
- CLI help agrees with the guide; documentation links and README command
  alphabetic order were checked. No root README command ordering changed.

## Frozen Study Replay

A read-only replay used the original banks and model from
`/tmp/nwkit-locus-stratified-20261002`. It reparsed all 40 saved true/estimated
gene-tree collections and reran the complete null calibration with the
original 39/19 replicate counts, seeds, and independent sequence/NJ/root
sampler for estimated trees. All 40 fit records and **1,160 bootstrap
searches** matched exactly, including probabilities, likelihood/contrast
bounds, selected grids/parents, P-values, and attempt counts.

Only the deliberately clarified explanatory `calibration.limitations`
string was excluded from JSON equality. No numerical field or other
calibration field was excluded. This validates unchanged scoring/sampling
with frozen banks; it is not a new two-million-draw bank rebuild or a new
accuracy/calibration study. The original study scripts were unchanged.

Audited source SHA256:

```text
nwkit/mul_locus.py: 5f2a22170472d452e8f057ed14151c06ce66643ad1f541b18ee769b4dc41bb96
nwkit/mul_locus_mc.py: 58c932eee7082d0ee8339eb0e21743d58dbc0ffe5036496e003c5078ccc8823e
nwkit/mul_locus_cli.py: 97653fa109c90a134677ea88fcdd371dfa072fabc28caeff7b67c55e1a7d8799
tests/test_mul_locus.py: 388a60ef8fc25502048e962eb5bbb10db1790641ec691e065446f70964d2b119
historical validation record: ee29afe198e962418f05b0f83e7834364b7ca774d83eb01853182a0c2520dc00
```

## Remaining Boundaries

No defect was found in the checked fixed-history bounded-coalescent law,
stratum weighting, or frozen-study score calculations. That is not proof
of general estimation accuracy or calibrated false-positive control.
Uncertainty in the generating null choice, broader composite-null error
control, unknown detection/species ages, alternate ascertainment, larger
families, model misspecification, and empirical gene-tree pipelines remain
research boundaries. The pilot's JC69/NJ calibration does not calibrate
IQ-TREE or any other empirical inference/rooting pipeline.

Full release/security/coverage/distribution gates, Python 3.10, other
platforms, empirical datasets, and GeneGalleon Docker/SIF integration were
not validated here. No GeneGalleon files or curated inputs were changed;
no release, version bump, commit/push, or Wiki deployment was performed.
