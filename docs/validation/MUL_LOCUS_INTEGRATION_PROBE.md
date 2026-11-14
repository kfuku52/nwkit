# Locus Conditional Integration Probe

Executed 2026-10-02 through 2026-10-03 on checkout version 0.43.37, based
on `0a4cea2881cf0559f124f68b03d4b0203759bb67` with the existing uncommitted
locus/MSC work preserved. This is exploratory evidence, not a production
adoption or a rewrite of the original finite-grid study.

See [the fixed protocol](../../examples/mul-locus/INTEGRATION.md) and
[historical study](MUL_LOCUS_GRID_CALIBRATION.md).

## Diagnosis And Implementation

All nine historical failing 30-family null draws reproduced from their
saved dataset seeds, generating points and replicate indices. The diagnosis
records every pattern/bank hit and the summed hits for the same species
copy-count vector, without pooling different candidate/grid banks for scoring.

Two failures contain a copy vector absent throughout an entire failed
hypothesis class: baseline allop-A-r5 has four X copies absent from all null
banks; baseline allop-B-r3 has four A copies absent from all alternative
banks. Genealogy integration alone cannot manufacture an unsampled hidden
copy configuration. These absences can reflect DL-history **or detection**
sampling; historical raw hidden histories were not retained, so those two
causes cannot be distinguished further. Baseline null-d0-ne1-r4 has only
topology absences; other draws combine topology and copy-vector absences in
individual banks. Absence is not proof of zero biological probability.

The research module `nwkit/mul_locus_integral.py` compares a paired histogram,
detection integration, and predeclared hybrid genealogy/detection integration.
It retains full hidden histories, uses complete-history daughter bounds,
and keeps the observed-size selection denominator. Histories with more than
four extant hidden loci use a genealogy draw, not a truncated history or a
retry after a resource failure. Weighted contributions use two-sided
empirical Bernstein bounds, not Clopper-Pearson on fractional counts.

The shared search/calibration kernel gained an optional scorer, used for
every observed and null search. Omission preserves the original scorer,
seeds, search scope, calibration and CLI output schemas. No new CLI
integration mode is exposed. The independent reference helper gained a
true-only opt-in path; its default sequence/NJ behavior is unchanged.

## Operational Failures Preserved

The first execution in `/tmp/nwkit-locus-integration-probe-20261002` reached
the detection-DP work cap in all three scenarios, before dataset scoring.
All 63 planned analyses are retained as bank-failed. Its nine source files
were archived and every archived SHA256 matched its frozen protocol.

The next execution in `/tmp/nwkit-locus-integration-overflow-probe-20261003`
exactly aggregated detection outcomes above four observed tips into an
overflow mass. It built all 20 baseline banks but stopped in the second
dataset when unnecessary NJ/midpoint inference raised Biopython's
`UnboundLocalError` on a zero-distance tree. This incomplete execution and
its matching source archive remain preserved, not counted as a completed
study or repaired by inventing an estimated root.

The final execution uses true-only independent data/null sampling and
retains unexpected per-bank/generation/calibration exceptions. It leaves
all scientific settings, budgets and base seeds unchanged. Removing unused
mutation/NJ draws changes realized data streams relative to the incomplete
execution, but not the true-genealogy observation law. This distinction is
explicit in the follow-up protocol.

## Independent Paired Results

```sh
PY=/tmp/nwkit-wgd-dev-20261001/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
"$PY" examples/mul-locus/integration.py \
  --output /tmp/nwkit-locus-integration-true-probe-20261003 \
  --source /tmp/nwkit-locus-grid-calibration-20261002
```

The program processed all 63 planned analyses: 21 independently generated
datasets, 10 families each, compared across three methods with identical
observed/null streams. All 60 paired banks completed, with 60,000 raw locus
histories. The hybrid integrated 25,359 small histories and retained 34,641
larger histories using their full genealogy/detection path. No histories
or failed analyses were dropped. The program exited **1** because failures
remain; this is not a fully successful accuracy-validation gate.

| Method | Planned | Completed | Score Support Failures | Independent Sampler Cap Failures |
| --- | ---: | ---: | ---: | ---: |
| histogram | 21 | 0 | 21 | 0 |
| detection RB | 21 | 18 | 0 | 3 |
| hybrid RB | 21 | 18 | 0 | 3 |

Both conditional methods completed all seven baseline and all seven turnover
datasets, and four of seven missing-ils datasets. Their three failing cases
are missing-ils null-d0-ne1, null-d1-ne1 and allop-B. Each reaches the same
independent bounded-coalescent rejection cap in the shared null stream.
The histogram stops earlier on score support, so absence of recorded sampler
failures in that method is not evidence of sampler robustness.

The 36 completed calibrations retain **2,736** complete null searches
(four points times 19 draws per analysis). Every stored point P-value,
MC-overlap endpoint and grid supremum was independently recomputed. All 36
observed fit records exactly replayed from the serialized weighted moments
and observations. All 63 case/method identities and 21 dataset seeds were
unique, and the final archived source hashes matched the protocol.

For each conditional method, all five completed alternative datasets had
point P <=0.05 and selected the correct parent; one further alternative
dataset failed. **No alternative passed the upper-MC decision rule**:
upper P was 1 in every completion. All contrast lower bounds were unbounded;
35 of 36 upper bounds were also unbounded. The one finite upper endpoint
was baseline null-d1-ne0/hybrid, 22.082247051022236.

No completed null reported an event, but two of 12 planned on-grid nulls
failed for each conditional method. Their realized failure-aware on-grid
report proportion remains 0-2/12. Neither that envelope nor zero observed
reports establishes an empirical 5% guarantee. One dataset per truth point
and 19 null draws per point cannot establish power or parental identifiability.

| Method | Evaluation Pattern/Bank Cells | Zero Point Estimates | Median Probability Interval Width |
| --- | ---: | ---: | ---: |
| histogram | 2,700 | 265 | 0.127662 |
| detection RB | 2,700 | 0 | 0.224506 |
| hybrid RB | 2,700 | 0 | 0.217397 |

These are paired, dependent cells from observed evaluation patterns, not
independent trials or an audit of all possible null patterns. Conditional
integration improved represented support in this probe, but **did not narrow
the MC intervals**. The generic bounded-variable interval's additive term
and the simultaneous observation-universe budget are substantial at 500
draws per stratum. No confidence level or scientific threshold was relaxed.
This low-budget comparison does not imply that the nine failures in the
earlier 100,000-draw study have been eliminated.

## Computational Cost

Combined paired-bank construction took 94.891 seconds. Summed method
analysis times were histogram 39.348 s, detection 252.767 s and hybrid
252.197 s. These are operational totals, **not a fair speed comparison**:
methods stop at different failures, and independent rejection sampling
dominates the three shared cap failures.

A separate isolated detection-DP benchmark used balanced 8/12-tip genealogies,
species colors A/X/B, detection 0.6, the same 100,000-work cap, one warmup,
three repetitions of 20 calls, and no concurrent numerical workload:

```sh
"$PY" examples/mul-locus/benchmark_detection.py
```

| Hidden Tips | Full Detection Median | Overflow Median | Full Python Peak | Overflow Python Peak |
| --- | ---: | ---: | ---: | ---: |
| 8 | 0.143437 ms | 0.249608 ms | 19,208 B | 9,520 B |
| 12 | 0.824365 ms | 1.383579 ms | 78,752 B | 11,576 B |

Overflow reduces allocation/state growth but is slower on these small/
moderate examples; no overall speedup is claimed. Target pattern and
overflow masses differed by at most 3.33e-15. A separate replay against
the first execution's archived pre-overflow module confirmed the same
probability equivalence and runtime/memory tradeoff. Memory here is peak
traced Python allocation, **not process RSS**. A 20-X-copy binomial detection
test completes within 1,000 charged operations with exact overflow mass.

## Verification

Runtime: CPython 3.12.14, macOS 27 ARM64, NumPy 2.5.3, SciPy 1.18.1,
msprime 1.4.4 and Biopython 1.88 in the isolated Python environment above.
Compiled-library imports and `pip check` passed. No required NWKIT
dependency was added.

```sh
NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  "$PY" tools/check.py quick -- \
  tests/test_mul_locus_integral.py tests/test_mul_locus_integration_study.py \
  tests/test_mul_locus.py tests/test_mul_locus_calibration.py \
  tests/test_mul_locus_study.py tests/test_mul_coalescent.py \
  tests/test_mul_msc.py tests/test_mul_msc_fit.py \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli.py tests/test_cli_contracts.py \
  tests/test_interface_conventions.py tests/test_output_transaction.py \
  tests/test_rooting_state.py tests/test_provenance.py \
  tests/test_util_tree_io.py tests/test_tree_outputs.py -x -rs
"$PY" tools/check.py quick -- tests/test_mul_locus_integral.py \
  tests/test_mul_locus_integration_study.py tests/test_mul_locus_study.py -x -rs
"$PY" tools/check.py test -- -m slow tests/test_mul_locus.py \
  tests/test_mul_locus_study.py tests/test_mul_msc_fit.py tests/test_mul_msc.py -rs
"$PY" -m ruff check examples/mul-locus/pilot.py \
  examples/mul-locus/integration.py examples/mul-locus/benchmark_detection.py
"$PY" tools/check_maintainability.py
git diff --check
```

- Broad quick lane: 608 passed, 15 slow cases deselected, no skips; Ruff
  checked/formatted 567 source/test files and mypy checked 250 source files.
- After the opt-in true-only repair, affected quick lane: 61 passed,
  two slow cases deselected, no skips. This adds two distinct cases, not 61
  new cases, to the broad lane's count. Static checks passed again.
- Explicit slow lane: 15 passed, 162 other cases deselected, no skips,
  including actual bank/evaluation worker death. **625 distinct selected
  tests passed**, not the full suite.
- Independent CTMC/mask oracles cover all candidate classes at three Ne
  values, nested daughter bounds, loss, invisible copies, detection 0/1,
  zero-duration branches, tiny daughter normalizers and repeated colors.
  An exhaustive non-Bernoulli IID example checks interval failure mass
  against the declared error budget across all 64-draw count configurations.
- Four old/default reference-family evaluations, null/allopolyploid and
  true/estimated views, exactly matched the archived implementation's
  three trees, sequence/history audits and final RNG state. This is a
  focused default-compatibility check, not a rerun of the historical study.
- Maintainability: 3,630 functions, mean 6.91, maximum 50; all hard limits
  respected, with existing baseline-comparison warnings.

## Evidence And Adoption Decision

Local final artifacts: `/tmp/nwkit-locus-integration-true-probe-20261003`.
Temporary local output is not claimed archived publicly or deployed.

```text
protocol.json: 7e3175897e2a88b50168f299d07976bd85110f0ac7313ab23714c117d63a6614
diagnosis.json: 6b433e3f5a7a12770bb5d750625a665295b0a8d4d9700a04c4d850696b5ee1d5
summary.tsv: 168e93282f1272486b23fc1ac5af2323295bfc4481e95e710c045c634de14803
source-snapshot.tar: 517ed4dee31cff12e2a48890ae22695fbeadf87aa102754c3d116ae587a96c4b
```

**Do not adopt as the CLI default.** Support representation improved in the
new low-budget probe, but MC interval precision did not. The independent
reference simulator also needs a separately validated rare-bound sampling
method before a larger held-out evaluation can reliably finish. Next work
should sharpen bounded-contribution intervals and investigate rare DL/
detection strata, with a fresh fixed protocol and no tuning of this study.

Full release/security/coverage/distribution gates, Python 3.10, other
platforms, GeneGalleon Docker/SIF integration and empirical gene-tree
pipelines were not run. No GeneGalleon files, curated inputs, default CLI
integration, model thresholds, historical results, versions or dependency
constraints were changed. No commit, push, release or Wiki deployment.
