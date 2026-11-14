# Finite Null Grid Input and Runner Audit

Executed 2026-10-02 on the uncommitted implementation based on
`0a4cea2881cf0559f124f68b03d4b0203759bb67`, checkout version `0.43.37`.
This audit follows the [finite-grid evaluation](MUL_LOCUS_GRID_CALIBRATION.md).
That record, the original pilot/audit records and all frozen study artifacts
were preserved, including their original source hashes and failure counts.

## Verified Findings and Repairs

1. **One-shot observations could change scores between hypotheses.**
   `search_banks` repeatedly consumed the same iterator, so later banks
   scored an empty dataset. Calibration then failed when taking its length.
   A one-shot bank iterator was also exhausted before null selection or
   subsequent searches. Observations and banks are now snapshotted before
   reuse. Regression tests compare complete iterator/list results in both
   calibration modes. Normal list inputs and their ordering are unchanged.
2. **A custom null sampler could change family count silently.**
   Returning too many families was accepted; an empty result reached the
   scorer instead of a batch-integrity check. Returned iterables are now
   snapshotted and their count must equal the observed count before any
   refit. No family is added, dropped or retried. These checks do not certify
   an arbitrary sampler's IID law or inference-error model.
3. **Numerical failure diagnostics omitted generating context.**
   Only `ValueError` received replicate/grid context; `ArithmeticError`
   escaped without it. Numerical errors now preserve their exception type
   and receive the same context, with explicit chaining. Failures still
   abort calibration and preserve the CLI output bundle.
   Invalid-replicate diagnostics also distinguish the nonnegative plug-in
   count from grid-supremum's positive-count requirement.
4. **Scenario selection changed reproducibility.**
   Scenario indices and bank seeds depended on request order, so a standalone
   missing-ils run did not reproduce that scenario in the default full run.
   Named indices are now fixed and recorded in the protocol. Subsets and
   reversed requests retain the full-protocol streams. All 105 default
   trial metadata/seeds and all three numerical model configurations match
   the frozen full study exactly.
5. **Worker death could lose planned denominators and metadata.**
   A broken process pool escaped the runner, leaving later trials/scenarios
   unrecorded. Bank-failure rows also omitted family counts, parameters and
   calibration seeds. The runner now retains returned partial banks, gives
   every affected trial complete metadata and a failure record, and continues
   independent scenarios. Evaluation workers persist a trial summary before
   returning; a pool failure recovers valid saved results/failures and marks
   unresolved trials `worker-failed`. Malformed or incomplete saved records
   remain preserved rather than counted as completions. No failed trial is
   resimulated. An unavailable output filesystem or main-process death can
   still prevent recording; incomplete outputs are not a passing gate.
6. **32-bit calibration-seed collisions were not checked.**
   The supplied scheme really collides across baseline case/replicate
   indices `(3, 11575)` and `(5, 13904)`: both produce `3590509910` with
   base seed 20261034. Such a large run would reuse a null stream, contrary
   to the fresh-stream protocol. A preflight now checks all planned trial
   seeds and fails with both identities before protocol writing/simulation.
   It does not silently choose replacement seeds. The original 105-trial
   protocol passes this check and retains its original numbers.

The first pre-repair targeted run demonstrated **15 failures**. Two additional
one-shot-bank regressions failed before that repair. A new recovery boundary
test also rejected treating an identity-only saved record as completed;
the recovery guard was strengthened accordingly. No scientific tolerances,
model rates, family-selection rules, thresholds or integration budgets were
changed to obtain passing checks.

## Executed Verification

Runtime: `/tmp/nwkit-wgd-dev-20261001/bin/python`, CPython 3.12.14 on
macOS 27 ARM64; NumPy 2.5.3, SciPy 1.18.1, msprime 1.4.4,
tskit 1.0.3 and Biopython 1.88. Compiled-library imports and `pip check`
passed before reuse. Research dependencies remain optional, not new required
NWKIT dependencies.

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
"$PY" tools/check.py quick -- -m slow tests/test_mul_locus_study.py -rs
"$PY" -m ruff check examples/mul-locus/calibration.py
"$PY" tools/check_maintainability.py
"$PY" -m nwkit mul-reconcile --help
"$PY" tools/check.py quick -- tests/test_mul_locus_calibration.py -k invalid_calibration -rs
git diff --check
```

- Broad quick lane: **562 passed**, no skips. All four original-GRAMPA
  reference comparisons executed. Ruff lint/format checked 564 files and
  mypy checked 249 source files.
- Existing explicit scientific slow lane: **13 passed**, no skips.
- Two newly added slow tests terminate actual child processes with
  `os._exit(17)`, during bank construction and evaluation. Both passed,
  retaining all 14 planned trial/failure records in each test. Mock-pool
  tests additionally cover partial returns and recovering a saved completion
  not delivered by the broken pool. These are controlled failure tests,
  not full scientific datasets.
- Combined coverage is **577 distinct selected tests**, not the full NWKIT
  suite. All applicable static checks passed after the final test additions.
- Maintainability hard limits passed: 3,608 functions, mean complexity 6.92,
  maximum 50, with existing baseline warnings. The check covers its configured
  source scope; the study script also passed explicit Ruff lint/format.
- After clarifying the final invalid-replicate message, all four affected
  validation cases and the static checks passed again. These are repeat
  cases, not additions to the distinct-test total.

### Frozen Result Compatibility

Read-only replay with the original stratified-pilot banks reproduced all
40 true/estimated fits and **1,160 plug-in calibration searches** exactly.
Only the previously clarified explanatory `calibration.limitations` string
was excluded, not any numerical or other calibration field. Estimated-tree
replays retained the original independent sequence/NJ/root sampler.

Four complete finite-grid calibrations also matched every stored fit and
calibration field, including **1,584 null-search rows**:

- baseline/null-d0-ne0-r1;
- turnover/allop-B-r4;
- missing-ils/allop-B-r2 (the changed decision versus fitted-null upper P);
- missing-ils/allop-A-r4 (the nonreported wrong-parent case).

The complete failure path for baseline/allop-A-r5 reproduced its original
exception/context exactly. All 105 default trial metadata records/seeds,
the three model configurations and all 21 original aggregate rows matched.
No six-million-draw bank rebuild or full 105-dataset accuracy rerun was
performed; these are compatibility replays using independently frozen banks.

### Operational Study Replay

The same separate low-budget protocol was run in fresh directories with
`--cpus 1` and `--cpus 2`:

```sh
"$PY" examples/mul-locus/calibration.py \
  --output /tmp/nwkit-locus-grid-audit-serial-20261002 --scenarios baseline \
  --samples 1000 --replicates 1 --families 2 --bootstrap 3 --cpus 1
"$PY" examples/mul-locus/calibration.py \
  --output /tmp/nwkit-locus-grid-audit-parallel-20261002 --scenarios baseline \
  --samples 1000 --replicates 1 --families 2 --bootstrap 3 --cpus 2
```

Both retained all seven planned trials: **six completed and one calibration
support failure**, and correctly exited 1. All 39 non-protocol files were
byte-identical; protocols differed only in requested output/CPU arguments.
All 33 non-protocol files of the original low-budget smoke run were exactly
reproduced; six per-trial `summary.json` files are additive. This checks
execution and bookkeeping, not false-positive control or detection power.

## Source Snapshot

```text
nwkit/mul_locus_mc.py: b54cf478a5a4ceb63292305bc44b68e06a6e73068571859835b4b6db9dc1904f
examples/mul-locus/calibration.py: 8893256ef2d4d1ca690a170917f01825c9629ce848b23f999bf858b1f3c532dd
tests/test_mul_locus_calibration.py: 4d61d4a9eb31f48713b161a68c385301ccd2cc54ec1b6c7de7083fce2ac0818c
tests/test_mul_locus_study.py: 6c1188f7fd9c19001379f8c624d15d826a37c55cac4cba635d79e0e4fe56a31b
docs/validation/MUL_LOCUS_GRID_CALIBRATION.md: 45507972b2c252173bd49247975d33f2248357f686a1e646a0ad5b8dc3b6c2aa
docs/validation/MUL_LOCUS_MC_AUDIT.md: 23277580a4f6d0358ec1de67487af613f4de1f44cb9a29c360d2176ae9d466be
docs/validation/MUL_LOCUS_MC_VALIDATION.md: ee29afe198e962418f05b0f83e7834364b7ca774d83eb01853182a0c2520dc00
```

## Scientific Boundaries

The maximum-P construction, full candidate/grid refitting, tie-inclusive
rank formula and independent frozen-bank assumptions were rechecked against
[Dufour (2006)](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf).
The scope remains the supplied finite null grid with matching observation
laws and complete replicates. It is not continuous-nuisance coverage or
empirical inference-pipeline calibration.

The original study still has **nine MC-support failures**, and detection
under missing-ils remains **2/10** with the predeclared upper-MC rule.
These repairs do not solve rare-pattern integration or wide score-MC bounds,
and do not establish uniform empirical accuracy. There is still no
production-adoption claim or GeneGalleon locus-mode integration.

Full release/security/coverage/distribution gates, Python 3.10, other
platforms and GeneGalleon Docker/SIF validation were not run in this phase.
No GeneGalleon files or curated inputs were changed. No version bump,
commit, push, release or Wiki deployment was performed.
