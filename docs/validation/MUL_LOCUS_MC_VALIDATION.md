# Experimental Locus DL + ILS Validation

Executed 2026-10-02 against uncommitted NWKIT source based on
`0a4cea2881cf0559f124f68b03d4b0203759bb67`, version 0.43.37.
This milestone implements an explicit small-family generative locus model
and finite-grid no-polyploidy/one-allopolyploid comparison. It completes a
research implementation and pilot, **not general false-positive calibration
or production WGD detection**. Default D+L and conditional MSC retain their
separate contracts. Historical [fixed MSC](MUL_MSC_VALIDATION.md) and
[MSC fitting](MUL_MSC_FIT_VALIDATION.md) records were not regenerated.

See the [model/CLI guide](../guides/MUL_LOCUS_MC.md) and the
[frozen study protocol](../../examples/mul-locus/README.md).

## Executed Checks

Host: macOS 27.0 ARM64, CPython 3.12.14, NumPy 2.5.3, SciPy 1.18.1;
interpreter `/tmp/nwkit-wgd-dev-20261001/bin/python`. Compiled-library import
preflight and `pip check` passed. Optional isolated research dependencies
were msprime 1.4.4, tskit 1.0.3 and Biopython 1.88; NWKIT dependency metadata
was not changed. Commands ran from the NWKIT repository root:

```bash
env NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python tools/check.py quick -- \
  tests/test_mul_locus.py tests/test_mul_coalescent.py tests/test_mul_msc.py \
  tests/test_mul_msc_fit.py tests/test_mul_reconcile.py \
  tests/test_mul_reconcile_reference.py tests/test_cli.py \
  tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_output_transaction.py tests/test_rooting_state.py \
  tests/test_provenance.py -x -rs

env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python tools/check.py test -- \
  -m slow tests/test_mul_locus.py tests/test_mul_msc_fit.py \
  tests/test_mul_msc.py -s -rs

/tmp/nwkit-wgd-dev-20261001/bin/python tools/check_maintainability.py
/tmp/nwkit-wgd-dev-20261001/bin/python -m pip check
/tmp/nwkit-wgd-dev-20261001/bin/python -m ruff check examples/mul-locus/
/tmp/nwkit-wgd-dev-20261001/bin/python -m ruff format --check examples/mul-locus/
```

Final focused quick gate: **434 passed**, 12 slow cases deselected, no skips,
47.28 seconds. Whole-repository Ruff/format (562 files) and incremental mypy
(249 source files) passed. All four original-GRAMPA reference cases executed.
Final explicit slow gate: **12 passed**, 110 non-slow cases deselected, no
skips, 11.36 seconds. These are **446 distinct selected tests**, not a full
suite or release result. Maintainability checked 3,604 functions, mean
complexity 6.92, all hard limits respected; unrelated baseline warnings
remain. No tolerance, biological parameter or complexity baseline was
weakened to pass. Both research scripts passed separate lint/format checks.
CLI help was checked against the guide; README command order remains ABC.

Intermediate gates found missing Counter annotations and a stale MSC
exclusion-message assertion. These were fixed and the final gates rerun.
An initial stability preflight against the older first-run bank file failed
because that historical record lacked the new explicit sample-count field;
the completed follow-up's full bank records were used for the actual study.
Shared-column resampling and identical seeded full-search replay then passed.

## Numerical and Contract Evidence

- The continuous linear locus birth/death process includes the ancestral
  stem, arbitrary SSD histories and extinction. Pure-birth, critical and
  subcritical ancestor counts (5,000 draws each) matched independent BD
  closed forms inside the declared binomial intervals.
- Count DP implements normalized daughter-bounded coalescence, following
  the [DLCoal construction](https://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html),
  with a separately specified experimental direct-parent MUL extension.
  Two-tip and three-tip daughter normalizers matched a full all-pair forest
  oracle within `2e-13`; 8,000 three-tip daughter genealogies matched its
  normalized topology law. This is not an implementation-equivalence claim
  for the complete DLCoal reconstruction program.
- At zero DL, 8,000 sampled null genealogies and 8,000 genealogies for each
  of both allopolyploid parents matched the independent full-forest CTMC.
  The allopolyploid checks sum over repeated-species colored topology
  classes, rather than matching a single arbitrary homoeolog assignment.
- Genealogy generation precedes detection; undetected extant loci remain
  in daughter conditioning. Both hypotheses use exactly the same family
  selection, `2 <= observed tips <= K`, and never drop input families.
- A weighted two-stratum example verifies that raw pattern probabilities
  and selection masses are normalized together, not averaged after separate
  selection. Stratified ancestral counts also matched the conditional BD
  closed form, including total observation mass one.
- Simultaneous binomial MC bounds cover the complete declared observation
  universe, not only observed/bootstrap-selected patterns. Zero estimates
  remain zero; missing-support failures have no pseudocount fallback.
  Underflowing positive stratum priors fail instead of deleting a stratum.
- Work caps now precede potentially dense lineage-transition allocations.
  Unrepresentable confidence precision fails before bank simulation.
  Node/state/selection caps abort, never discard a problematic history.
- A mock null bootstrap reselected different null parameter grids and
  alternative parents across replicates. Every real replicate repeats the
  complete supplied search; redundant null attachment ages are deduplicated.
- Both integration methods produced byte-identical serial/process output
  bundles. Late writer failures and path aliases preserve all five old
  files. Model JSON and every output are included in provenance hashes;
  `--audit` cannot alias the model, checks or calibration output. This fixes
  missing common path registrations, not just a command-local check.

## Independent Event-Comparison Pilot

The protocol was written before execution. Twenty datasets contain 30
families each: **600 independent family genealogies and 40 paired
true/estimated-tree analyses**. The two views are paired, not 40 independent
datasets. Evaluation seed `20261025`; H1 X and parent candidates A/B are
supplied. Species ages, detection 0.9, loss 0.03 and the ancestral stem 0.5
are fixed. The allopolyploid truth's duplication 0.05 and age 0.5 lie outside
the fitted discrete grid. The null has four unique parameter points and the
alternative has 16 candidate/grid points.

The generator uses a global simultaneous-locus SSA, independent of the
production depth-first sampler. Full ancestry uses msprime's standard
coalescent with diploid ploidy 2 and one haploid sample per distinct extant
locus, as specified by [msprime's ancestry API](https://tskit.dev/msprime/docs/latest/ancestry.html#ploidy).
Conditional genealogy rejection enforces every daughter bound without
removing undetected extant loci. All raw histories and acceptance counts
are retained. Estimated trees use 600-site JC69 alignments, rate 0.01,
NJ and midpoint rooting without access to the true root.

### First Run: Missing Integration Support

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/pilot.py \
  --output /tmp/nwkit-locus-pilot-20261002
```

Ten thousand IID selected draws per bank, bank/bootstrap seed `20261026`:
the study exited nonzero with **17 completed and 23 failed analyses**.
Only three observed-tree analyses lacked finite support under either null or
alternative; the other 20 failures occurred inside null bootstraps. The
complete 40-row summary and all failed cases were preserved. Successful
null analyses were 5/10 true and 3/10 estimated; successful allopolyploid
analyses were 9/10 true and 0/10 estimated. These incomplete denominators
must not be reported as full-study false-positive or accuracy estimates.

### Fixed Resource Follow-Up

Before rerunning, the protocol separately declared an independent
bank/bootstrap seed `20261027` and 100,000 raw draws per bank using
ancestral first-event stratification. The scientific model, evaluation
seed, family sets and decision threshold were unchanged. The 20 complete
`families.jsonl` files, including genealogies, RNG seeds and raw alignments,
are byte-identical to the first run. The original failed study remains.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/pilot.py \
  --output /tmp/nwkit-locus-stratified-20261002 --samples 100000 \
  --bank-seed 20261027 --integration ancestral-stratified
```

This follow-up exited zero: **40/40 completed, no omitted family or failed
analysis**. The 20 banks contain two million total raw draws. True-tree null
calibration has 39 replicates per analysis; estimated-tree calibration has
19 and repeats the independent sequence/NJ/rooting pipeline under the fitted
null. All **1,160 null bootstrap searches** completed and reselected null
parameters and all alternative parents/grid points.

The protocol reports an event only when the MC-overlap upper P-value is
at most 0.05; no decision is made from parent stability frequencies.

| Input View | Completed | Null Reported Events | Allopolyploid Reported Events | Numerical Parent Correct |
| --- | ---: | ---: | ---: | ---: |
| True gene trees | 20/20 | 0/10 | 10/10 | 9/10 |
| Estimated gene trees | 20/20 | 0/10 | 10/10 | 10/10 |

All positive-case P-values reached the minimum resolution: 0.025 for true
trees, 0.05 for estimated trees. This means no bootstrap exceedance at these
budgets, **not exact P-values or strong tail precision**. All null upper MC
P-value bounds were 1; null intervals remain broad. Five replicates per null
condition cannot establish a 5% false-positive guarantee or uniform
composite-null calibration. Estimated-tree results only calibrate this
JC69/NJ/root procedure, not IQ-TREE or real-data curation.

The wrong numerical parent in `allop-A-r1/true` had log estimate -88.1763
versus -88.2038 for the best true-parent point. Their MC intervals overlap;
this is not a confidently wrong identified parent. Both views leave parents
A/B MC-overlapping in two of five parent-A datasets; the remaining three
parent-A and all five parent-B datasets resolve one parent. Matching
numerical winners must not be conflated with resolved biological estimates.
Colored rooted topology accuracy was 0.7333-0.9667 across datasets.

## Separate Site-Bootstrap Stability

The additional protocol was declared before its execution: 19 within-family
column bootstraps per original dataset, seed `20261029`, no family dropping.
All sequences in a family share the resampled column indices. Each bootstrap
repeats NJ/rooting and searches the same complete frozen integration banks.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python examples/mul-locus/stability.py \
  --study /tmp/nwkit-locus-stratified-20261002 \
  --output /tmp/nwkit-locus-site-bootstrap-20261002
```

All **380/380 bootstrap analyses completed**. Parent-A numerical recovery
counts were 18, 16, 19, 19 and 19 of 19 in its five datasets; parent-B was
19/19 in all five. Across all cases, MC-parent overlap was A/B in 240/380,
only A in 47/380 and only B in 93/380. Null datasets also receive an
alternative-parent ranking by construction; that ranking does not report
an event. These are descriptive frequencies conditional on alignments and
banks, not event P-values, posterior probabilities or confidence intervals.

## Performance and Provenance

A singleton conditional-choice shortcut preserves the uniform draw consumed
by NumPy choice. Before/after measurements used the same 2,000 selected
families, H2 B, duplication 0.06, loss 0.03, Ne/age 0.5, seed `20261020`,
100-family warm-up, three repetitions and tracemalloc on this same host.
Wall seconds were 10.668/9.934/9.778 before and 5.366/4.801/4.766 after;
steady Python-traced peaks remained about 594 KB. Every run had 2,822
attempts and identical histogram SHA256
`ac96d2e6e7d12f26c7afe8edcf9aa659c9a70f95c1dd562f7efe131df13c0bf8`.
These are instrumented local measurements, not a general throughput claim.

Each study's `protocol.json` records its actual source hashes. Final
hardening after the follow-up started added annotations, strict root-count
type/underflow checks, earlier resource/precision guards, stderr calibration
and clearer exclusions/path registration. Its scientific sampling law and
score calculations were not changed. Frozen study source hashes therefore
must not be presented as final-source hashes, nor reused banks as newly
resimulated final-source banks.

After final hardening, a read-only replay loaded the saved follow-up banks,
reparsed all 40 saved gene-tree collections, and reran `calibrate` with the
original 39/19 bootstrap counts and seeds. Estimated-tree replicates again
used the independent sequence/NJ/rooting sampler. Every full fit and
calibration JSON value matched exactly, including all 1,160 bootstrap
searches. This validates current rescoring/bootstrap behavior with frozen
banks; it is not a new two-million-draw integration-bank rebuild.

Study protocol / complete summary SHA256:

```text
first protocol: 38153ec042f5a0539a395ad70f60b6a32bb40ac74b33b7db00393a39bce31692
first summary:  5641f1d4ad55876ddce344dc0bfb337a5af932c05a69f48ac8578740cb27c1bf
follow protocol: e3ef963da1397db32d9aad201f187e80081f3542ba2a564fb74ecc2d55838a38
follow summary: eac932ef4e8a104f600f5b29592efe711cd22a11bf128ea568f3a306bdc832a0
stability protocol: 96a0e56f6db2cb06002db90c8b530e672bf95b66bfb23151dd51fca793ba264c
stability summary: 39fcdad239a3e672726db945b4bcb39a07c44ea3354cf8d15556b90679aa7bb8
```

Final checked core/research SHA256:

```text
nwkit/mul_locus.py: 9b008f428f730465ca5d9e9dbc9c01c662c4739c19343ce8ba516b68792ece51
nwkit/mul_locus_mc.py: 1dbcda1d5f241ec5b84cdf989eb479cc72844e6512eff4a6e743cb6f589717ea
nwkit/mul_locus_cli.py: aa255d3de8227d33a05e0ceba70c8b2f9d8f2032db3d5fa48c73e5baea82f024
nwkit/mul_msc_model.py: fd8d37c00ab3c981328f4dd85b785750047b84781ed02be811a0af5aa681bab9
nwkit/mul_coalescent.py: cabb89e6bed59c4e96c1c26615993310f8c84e51b58d967c19c0770d609a4173
nwkit/provenance.py: c460b89a9cf6afb5d9d4e2d0f3c862e5f64b3cc2c4e1b576463bed907a64e739
tests/test_mul_locus.py: 315320109370520f70788d1bb30749c1527874e19d72eb544e9c3d2d94d52603
examples/mul-locus/pilot.py: 858048382b20e2a9a3c68f70320e766e02d38c8e8bd352e664484d35f5c27772
examples/mul-locus/stability.py: acb9652735189673b373600b0b5727418be46eebc6ba78476d845f54b55537cc
```

## Remaining Boundaries

This milestone does not establish continuous MLE, biological parameter
identifiability/intervals, general calibrated composite-null error control,
arbitrary family ascertainment, uncertain detection/species ages, heterogeneous
population rates, autopolyploidy, multiple events, hemiplasy or exchange.
Larger families and rare multi-SSD histories still need substantially more
integration precision. Inference-error-aware calibration must reproduce the
actual empirical sequence/tree/root pipeline, not reuse this pilot's error law.

Full release/security/coverage/distribution gates, minimum Python 3.10,
other platforms, empirical datasets and Docker/SIF runtime validation were
not run for this milestone. GeneGalleon integration, release/version bump,
commit/push and Wiki deployment are step 6 and were not performed. Existing
GeneGalleon worktree changes and curated research data were left intact.
