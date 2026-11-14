# Finite Null Grid Calibration and Independent Evaluation

Executed 2026-10-02 on the uncommitted locus/MSC implementation based on
`0a4cea2881cf0559f124f68b03d4b0203759bb67`, checkout version `0.43.37`.
This is new evidence, not a rewrite of the
[original pilot](MUL_LOCUS_MC_VALIDATION.md) or
[subsequent numerical audit](MUL_LOCUS_MC_AUDIT.md).
Their protocols, artifacts and source-hash records remain historical.

## Implemented Method

Explicit `--locus-null-calibration grid-supremum` generates B complete null
datasets at every distinct supplied no-polyploidy parameter point. Each
dataset repeats the full null/alternative parameter and parent search,
using the same independently frozen scoring banks as the observed dataset.
The reported P-value is the maximum of the per-point rank P-values,
`(1 + number of ties or exceedances)/(B + 1)`. Both score-MC overlap
endpoints are independently maximized too. Results retain every generating
point, every refitted point/parent, B per point, total work and the two
potentially different least-favorable points.

No data-dependent null screening, mixture-null averaging, pseudocounts,
family dropping or partial-result calibration is used. A missing-support
or resource failure aborts the CLI output bundle and preserves previous
files. The default remains `plug-in`, including its original seed stream
and calibration schema. The new option is rejected outside `locus-mc`.

With fixed independent banks and identical data/null observation and
selection laws, the true-point rank P-value is conservative if the true
null belongs to the supplied finite grid. Taking the maximum cannot reject
when that true-point P-value would not reject. This is a finite-grid
specialization of [Dufour's maximized Monte Carlo testing
construction](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf),
not a guarantee over off-grid continuous parameters, model
misspecification, failed analyses or empirical gene-tree pipelines.
Score-MC bounds are not finite-B sampling confidence intervals.

## Executed Checks

Runtime: `/tmp/nwkit-wgd-dev-20261001/bin/python`, CPython 3.12.14,
macOS 27 ARM64, NumPy 2.5.3, SciPy 1.18.1, msprime 1.4.4,
tskit 1.0.3 and Biopython 1.88. Compiled-library imports and `pip check`
passed before reuse. msprime/Biopython are research dependencies, not new
required NWKIT runtime dependencies.

```sh
PY=/tmp/nwkit-wgd-dev-20261001/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  "$PY" tools/check.py quick -- \
  tests/test_mul_locus.py tests/test_mul_locus_calibration.py \
  tests/test_mul_locus_study.py tests/test_mul_coalescent.py \
  tests/test_mul_msc.py tests/test_mul_msc_fit.py \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli.py tests/test_cli_contracts.py \
  tests/test_interface_conventions.py tests/test_output_transaction.py \
  tests/test_rooting_state.py tests/test_provenance.py \
  tests/test_util_tree_io.py tests/test_tree_outputs.py -x -rs
"$PY" tools/check.py test -- -m slow \
  tests/test_mul_locus.py tests/test_mul_msc_fit.py tests/test_mul_msc.py -rs
"$PY" tools/check_maintainability.py
"$PY" -m nwkit mul-reconcile --help
"$PY" tools/check.py quick -- tests/test_mul_locus_calibration.py -k cli_grid -rs
git diff --check
```

- Broad quick lane: **541 passed**, 13 slow cases deselected, no skips.
  Ruff lint/format checked 564 files; mypy checked 249 source files.
  All four original-GRAMPA reference comparisons executed.
- Explicit slow lane: **13 passed**, 149 other cases deselected, no skips.
  Together these are **554 distinct selected tests**, not the full suite.
- New cases cover all-null-point evaluation, full search/refitting,
  maximum rather than average aggregation, ties, separate least-favorable
  points, exact seed namespaces, reordered banks, serial/process equality,
  per-replicate output fields, invalid options and bundle preservation.
  An exhaustive discrete Bernoulli rank-test example enumerates every
  three-replicate outcome under both supplied true null points; its
  rejection probability is below the predeclared level. This checks the
  finite-grid construction, not general biological accuracy.
- Study tests retain every planned denominator after bank/calibration
  failures and verify failure-aware ranges and all on-grid truth points.
- Maintainability: 3,608 functions, mean complexity 6.92, maximum 50;
  all hard limits passed, with existing baseline-comparison warnings.
- After the final help-only clarification of `--locus-bootstrap`, actual
  command help, the static checks and both affected CLI cases passed
  again (2 passed, 12 deselected). These are repeat checks, not additional
  distinct cases in the total above.

## Default-Mode Compatibility

A read-only replay used the original model and frozen banks in
`/tmp/nwkit-locus-stratified-20261002`. All 40 saved true/estimated
gene-tree datasets were reparsed and calibrated with their original
39/19 replicate counts, seeds and independent sequence/NJ/root sampler.
All 40 complete fit records and **1,160 bootstrap searches** matched
exactly after the new calibration refactor, including selected
parameters/parents, probability and contrast bounds, P-values and attempts.

Only the previously clarified explanatory `calibration.limitations` string
was excluded from equality. No numerical or other calibration field was
excluded. This is unchanged default behavior with historical banks, not
a new bank rebuild or independent accuracy study.

## Independent Study Protocol

The [separate frozen protocol](../../examples/mul-locus/CALIBRATION.md)
uses independent global SSA/msprime data and fresh integration banks in
three fixed DL/ILS/detection scenarios. Each scenario includes every one
of its four null grid points, an off-grid null diagnostic and both
allopolyploid parents, with five independent datasets per truth case.
Thirty families per dataset give 105 planned datasets and 3,150 selected
families. Each observed dataset requests 99 fresh null replicates per
generating point and full-search calibration. Decisions use the maximum
upper score-MC P bound at a fixed 0.05 threshold.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/calibration.py \
  --output /tmp/nwkit-locus-grid-calibration-20261002 --cpus 4
```

The protocol was written before evaluation, including configurations,
seeds, dependency versions and ten used-source hashes. The output keeps
all planned cases, successful fits and complete calibrations, failures,
histories, true trees, 600-site JC69 alignments and NJ/midpoint trees.
Only the true rooted topologies are event-calibrated here. Estimated
trees are retained for later pipeline-aware work.

A preceding small operational run used a separate output directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/calibration.py \
  --output /tmp/nwkit-locus-grid-smoke-20261002 --scenarios baseline \
  --samples 1000 --replicates 1 --families 2 --bootstrap 3 --cpus 2
```

It retained all seven planned cases: six completed and one failed during
calibration because no finite class score remained. The study continued
other cases and exited 1, as intended. This low-budget run checks failure
bookkeeping, not accuracy; no scientific settings in the main protocol
were changed in response to its outcomes.

## Main Study Results

All 105 planned datasets were processed: **96 completed and 9
calibration failures**. The program exited 1 because failures remained;
this is not a fully successful accuracy-validation gate. All 60 banks
completed, consuming 6,000,000 raw draws. All 3,150 evaluation families
and their paired sequences/trees were retained. Completed calibrations
retain 38,016 full null searches (96 times 4 times 99); this count does
not include partial work before the nine failed calibrations.

| Scenario | Planned | Completed | Calibration Failed |
| --- | ---: | ---: | ---: |
| baseline | 35 | 28 | 7 |
| turnover | 35 | 33 | 2 |
| missing-ils | 35 | 35 | 0 |

The following ranges allow every failed analysis either decision. They
bound the realized report proportion in this fixed study, not the unknown
population false-positive rate or power. Groups contain different truth
points; no pooled IID-binomial confidence claim is made.

| Scenario | Truth Class | Planned | Completed | Failed | Event Reports | Failure-Aware Report Range | Correct Alternative Parents |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| baseline | on-grid null | 20 | 18 | 2 | 0 | 0-10% | - |
| baseline | off-grid null | 5 | 4 | 1 | 0 | 0-20% | - |
| baseline | allopolyploid | 10 | 6 | 4 | 6 | 60-100% | 6/6 completed; 4 failed |
| turnover | on-grid null | 20 | 19 | 1 | 0 | 0-5% | - |
| turnover | off-grid null | 5 | 5 | 0 | 0 | 0-0% | - |
| turnover | allopolyploid | 10 | 9 | 1 | 8 | 80-90% | 9/9 completed; 1 failed |
| missing-ils | on-grid null | 20 | 20 | 0 | 0 | 0-0% | - |
| missing-ils | off-grid null | 5 | 5 | 0 | 0 | 0-0% | - |
| missing-ils | allopolyploid | 10 | 10 | 0 | 2 | 20-20% | 9/10 completed; 0 failed |

All 57 completed on-grid nulls and 14 completed off-grid diagnostics had
no event report, but four failed null analyses remain unresolved. Neither
zero observed reports nor the off-grid diagnostics establish general
5% calibration. The three failed on-grid null analyses alone leave the
realized pooled on-grid report proportion between 0 and 3/60; that is a
failure envelope, not a 95% error-control result.

Per-truth-point uncertainty is large. The following pointwise 95% binomial
envelopes use the planned denominator of five and allow each failed
analysis either outcome. They are conditional on shared independently
frozen banks, not simultaneous bounds across cases.

| Scenario | Truth Case | Completed/5 | Failed | Event Reports/5 | Failure-Aware 95% Envelope |
| --- | --- | ---: | ---: | ---: | --- |
| baseline | null-d0-ne0 | 5 | 0 | 0 | 0-52.18% |
| baseline | null-d0-ne1 | 4 | 1 | 0 | 0-71.64% |
| baseline | null-d1-ne0 | 5 | 0 | 0 | 0-52.18% |
| baseline | null-d1-ne1 | 4 | 1 | 0 | 0-71.64% |
| baseline | null-off-grid | 4 | 1 | 0 | 0-71.64% |
| baseline | allop-A | 3 | 2 | 3 | 14.66-100% |
| baseline | allop-B | 3 | 2 | 3 | 14.66-100% |
| turnover | null-d0-ne0 | 4 | 1 | 0 | 0-71.64% |
| turnover | null-d0-ne1 | 5 | 0 | 0 | 0-52.18% |
| turnover | null-d1-ne0 | 5 | 0 | 0 | 0-52.18% |
| turnover | null-d1-ne1 | 5 | 0 | 0 | 0-52.18% |
| turnover | null-off-grid | 5 | 0 | 0 | 0-52.18% |
| turnover | allop-A | 5 | 0 | 5 | 47.82-100% |
| turnover | allop-B | 4 | 1 | 3 | 14.66-99.49% |
| missing-ils | null-d0-ne0 | 5 | 0 | 0 | 0-52.18% |
| missing-ils | null-d0-ne1 | 5 | 0 | 0 | 0-52.18% |
| missing-ils | null-d1-ne0 | 5 | 0 | 0 | 0-52.18% |
| missing-ils | null-d1-ne1 | 5 | 0 | 0 | 0-52.18% |
| missing-ils | null-off-grid | 5 | 0 | 0 | 0-52.18% |
| missing-ils | allop-A | 5 | 0 | 1 | 0.51-71.64% |
| missing-ils | allop-B | 5 | 0 | 1 | 0.51-71.64% |

### Detection and Parent Recovery

All 25 completed allopolyploid datasets had point P <=0.05, but only
16 passed the predeclared **upper score-MC bound** rule. Under missing-ils,
the point rule reported 10/10 while the upper-bound rule reported 2/10.
Thus limited integration precision is a substantial detection-power
boundary in this study, not evidence that the numerical point statistic
failed to separate these alternatives. The decision rule was not relaxed.

The fitted-null upper bound from the same grid-mode draws would report
17 events versus the supremum's 16. The differing case was missing-ils,
allop-B replicate 2: fitted-point upper P=0.04, supremum upper P=0.07,
supremum point P=0.01. Across all 96 completions, the supremum point P
was strictly larger than the fitted-point P in 50 cases; the upper bound
was strictly larger in four. This paired comparison is not a fresh test
of the original plug-in seed stream.

Alternative-parent selection was correct in 24/25 completed alternative
datasets, with five further alternative analyses failed. The sole wrong
parent was missing-ils allop-A replicate 4, which selected B but did not
report an event (upper P=0.33). All 16 reported events selected the correct
parent in this small study. These conditional counts do not establish
parental identifiability or biological confidence.

### Failure Diagnosis

All nine failures occurred during null calibration, not integration-bank
construction, evaluation-family generation or observed-data fitting.
Replaying each failing selected family set with its recorded dataset
seed, generating grid and replicate reproduced missing support for at
least one entire hypothesis class. Finite-bank counts below refer to
banks that represented every pattern in that complete 30-family draw.

| Scenario | Dataset | Generating Grid | Replicate | Finite Null Banks | Finite Alternative Banks |
| --- | --- | ---: | ---: | ---: | ---: |
| baseline | null-d0-ne1-r4 | 6 | 17 | 4 | 0 |
| baseline | null-d1-ne1-r2 | 6 | 20 | 4 | 0 |
| baseline | null-off-grid-r2 | 4 | 3 | 4 | 0 |
| baseline | allop-A-r3 | 6 | 47 | 4 | 0 |
| baseline | allop-A-r5 | 6 | 14 | 0 | 16 |
| baseline | allop-B-r3 | 6 | 24 | 0 | 0 |
| baseline | allop-B-r4 | 6 | 67 | 1 | 0 |
| turnover | null-d0-ne0-r4 | 6 | 5 | 4 | 0 |
| turnover | allop-B-r3 | 4 | 80 | 4 | 0 |

For example, the baseline draw for allop-A-r5 contained the four-X
ladder topology once: its total null-bank hit count was zero while the
alternative-bank total was 55. The allop-B-r3 draw contained a four-A
ladder with zero hits across both classes. Some other failures had
patterns represented in different banks, but no single alternative
bank represented the complete dataset. Pooling those banks would change
the hypothesis and is not a repair.

These are finite-histogram support failures, not proof that the
biological model assigns zero probability. The chosen 100,000-raw-draw
budget per bank and current stratification do not suffice throughout
calibration. No failed draw was discarded, retried with a different
seed, smoothed or converted into a negative event. Improving rare-pattern
integration and narrowing score-MC bounds is the next research priority;
any revised integration/budget study needs its own fixed protocol.

## Artifact Audit and Source Snapshot

A read-only independent artifact check passed: all 105 planned case keys
and calibration seeds were unique and present; all family, sequence and
tree counts matched; all 60 banks retained the declared raw/stratum
counts; every completed calibration had four times 99 ordered draw
records. Every per-point P-value, supremum endpoint, least-favorable
tie-break, event/parent field and per-case failure-aware binomial envelope
was recomputed from stored results. All nine failure records remained
present. All ten protocol source hashes still matched after evaluation.

Artifacts are local in `/tmp/nwkit-locus-grid-calibration-20261002`;
the scripts and this record are portable, but this document does not
claim that temporary outputs are archived or publicly available.

```text
protocol.json: 8bdc51f7163437d0ae4b3dfbcf8a570ea96b29511c5a223dddb580fdd3a559cc
summary.tsv: 5ce9cdbf0dab085bdcb311b194db7b192439a8e059861c39b5560e9db6a9f1ab
rates.tsv: bc28d51ca2aa8da2a398084a500dc4d058396576a86f5c42b2f13421c40fa890
examples/mul-locus/calibration.py: 6866fcf7e3f581776a25e642a1f0e79b8274675cf58b9a084d8e94f66a114a9b
examples/mul-locus/pilot.py: 858048382b20e2a9a3c68f70320e766e02d38c8e8bd352e664484d35f5c27772
nwkit/mul_locus.py: 5f2a22170472d452e8f057ed14151c06ce66643ad1f541b18ee769b4dc41bb96
nwkit/mul_locus_mc.py: 09070220ffbdd6cc164c44ea1580ac8e27ec6b87ff2b1dbf2459d932f6964d89
nwkit/mul_locus_cli.py: f8919aa79298567ce45c84ad519c6147fa33f3577dd6bd3621f31165ea0dae89
nwkit/mul_msc_model.py: fd8d37c00ab3c981328f4dd85b785750047b84781ed02be811a0af5aa681bab9
nwkit/mul_coalescent.py: cabb89e6bed59c4e96c1c26615993310f8c84e51b58d967c19c0770d609a4173
nwkit/mul_msc_fit.py: 29ea640c73469096ed7c22b88166988007866d943a7b3b3467afe8cb396c7c5e
nwkit/species_parser.py: 99659e0e57a4ad6b6b7387cd5a130c7cbc05df163fd8d2ab3884cb4e74af664d
nwkit/util.py: 9495ff0f3070348fcc79c2cf6b23919fdaf320e5735f63415b227ce30e9c4c72
nwkit/mul_reconcile.py: 590ddc5c15cc5e687b7dc2327fa9dc2bd34e101d333d8875b4045a18f1cb2c41
nwkit/mul_reconcile_cli.py: 6b396dbc979a5dc30baf2055796b8f83cff9dd4ec554f9ffb0a1de8a87a6041a
tests/test_mul_locus_calibration.py: 7696bb909292dfb0124a97b65b60e6005381e83b2c937fe12ceb5b9aceb3285e
tests/test_mul_locus_study.py: 350daa887c81097700fcabc6812c3d423b437e17fb982a7b07660d79681ed65b
historical validation record: ee29afe198e962418f05b0f83e7834364b7ca774d83eb01853182a0c2520dc00
historical audit record: 23277580a4f6d0358ec1de67487af613f4de1f44cb9a29c360d2176ae9d466be
```

## Remaining Boundaries

Five datasets per truth point cannot establish a uniform empirical 5%
false-positive guarantee. Completed-only rates cannot remove support
failures from that assessment. Off-grid nulls are misspecification
diagnostics outside the finite-grid validity claim. Parent recovery is
selection within the alternative class, not evidence of an event when
the event test does not reject.

Unknown detection or species ages, other ascertainment, larger families,
multiple events, other polyploid modes and empirical inference/rooting
pipelines remain unvalidated. The original pilot's JC69/NJ calibration
does not calibrate IQ-TREE or other real-data workflows.

Full release/security/coverage/distribution gates, Python 3.10, other
platforms and GeneGalleon Docker/SIF integration were not run in this
phase. No GeneGalleon files, curated inputs or historical study records
were changed. No version bump, commit, push, release or Wiki deployment
was performed. This experimental mode is not adopted for production.
