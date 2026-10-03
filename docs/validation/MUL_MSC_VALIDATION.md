# Conditional MUL/MSC Prototype Validation

Executed 2026-10-02 against uncommitted source based on NWKIT
`0a4cea2881cf0559f124f68b03d4b0203759bb67` (version 0.43.37).
This is the first milestone: a fixed-parameter MSC kernel and conditional
second-parent comparison for one disomic allotetraploid event. It is not a
WGD test, a full AlloppNET implementation, or a DLCoal likelihood.
See the [model assumptions and CLI guide](../guides/MUL_MSC.md).

## Executed Checks

Host runtime: macOS 27.0 ARM64, CPython 3.12.14, NumPy 2.5.3, SciPy 1.18.1,
using `/tmp/nwkit-wgd-dev-20261001/bin/python`. The documented import preflight
and `pip check` passed. Commands ran from the NWKIT repository root:

```bash
env NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/nwkit-wgd-dev-20261001/bin/python tools/check.py quick -- \
  tests/test_mul_coalescent.py tests/test_mul_msc.py \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli.py tests/test_cli_contracts.py \
  tests/test_interface_conventions.py tests/test_output_transaction.py \
  tests/test_rooting_state.py -x -rs

/tmp/nwkit-wgd-dev-20261001/bin/python tools/check.py test -- \
  -m slow tests/test_mul_msc.py -s -rs

/tmp/nwkit-wgd-dev-20261001/bin/python tools/check_maintainability.py
```

The focused quick gate passed **325 tests**, with three slow pilot cases
deselected and no skips. Whole-repository Ruff, formatting (555 files) and mypy
(244 source files) passed. All four original-GRAMPA reference tests were
explicitly enabled and executed. The separate final pilot run passed those
**three tests**, with 48 non-slow cases deselected and no skips. These are
328 distinct selected tests, not a full-suite or release result.
Maintainability respected every hard limit; existing unrelated baseline
warnings remain. No thresholds, baselines or scientific models were loosened
to obtain passing results.

## Numerical Reference Evidence

- Pure-death lineage transitions matched an independent hypoexponential
  closed form using adaptive-precision `Decimal`: 216 combinations of
  `k=2..7`, `j=1..k`, and times `1e-100, 1e-12, .001, .1, .5, 2, 20, 1000`.
  Each complete transition row normalized to one.
- Three-taxon gene topology probabilities matched
  `1 - (2/3)*exp(-t)` for the concordant topology and `exp(-t)/3` for each
  discordant topology, including zero, tiny, and 10,000-unit branches.
- All 150 topologies across three four-tip population trees (15 each) and
  one five-tip tree (105) matched an independent full gene-forest CTMC.
  That reference retains every pair merger, including histories incompatible
  with the target topology; it does not reuse the production configuration
  recurrence. Complete topology distributions normalized to one.
- An internal H2 donor with two polyploid descendant species and four latent
  homoeolog assignments matched the independent forest reference over all
  15 four-tip topologies. The 105 five-tip unknown-homoeolog probabilities
  also normalized to one. Assignment probabilities are averaged, not summed
  without normalization or maximized.
- Missing population tips, child order, and generations/Ne versus coalescent
  units satisfied their invariance checks. A family with no observed
  polyploid copies acquired no artificial candidate-specific signal.
- Zero-length and extreme-time calculations passed; probabilities outside
  the representable log range failed clearly. Work limits failed before
  excessive matrix allocation and never silently truncated histories.

## Known-Parameter Parent-Only Pilot

Seed `20261005`; NumPy `default_rng` was reset for each time scale. Truth
probabilities for all 105 rooted five-tip gene topologies came from the
independent all-pair forest CTMC, averaging the two homoeolog assignments.
For each true H2 and scale, 50 replicates each sampled 500 independent
families by multinomial counts. Both scoring methods saw identical counts
and the same three parent candidates. Species/hybridization durations and
population scale were fixed to their true values for MSC.

| Shared time multiplier | True H2 | MSC correct / 50 | D+L parent-only correct / 50 |
| --- | --- | --- | --- |
| 0.1 | A | 50 | 0 |
| 0.1 | B | 50 | 9 |
| 0.1 | C | 50 | 50 |
| 0.5 | A | 50 | 49 |
| 0.5 | B | 50 | 50 |
| 0.5 | C | 50 | 50 |
| 2.0 | A | 50 | 50 |
| 2.0 | B | 50 | 50 |
| 2.0 | C | 50 | 50 |

This toy pilot recovered 450/450 parents with MSC and 358/450 with D+L.
It is **not an equal-information accuracy comparison**: MSC receives true
coalescent durations while D+L ignores branch lengths. The D+L comparison
excludes baseline and autopolyploid candidates; it is not the complete default
workflow. No unknown parameters, estimated gene trees, root uncertainty,
linkage, losses, missing copies or small-scale duplication were simulated.
These numbers cannot establish real-data accuracy or WGD detection rates.

Tests require independent probability agreement and the correct
population-expectation likelihood ordering; no finite-sample recovery
threshold was selected or tuned as a CI requirement.

## CLI, Compatibility and Output Evidence

The [guide example](../guides/MUL_MSC.md#example) ran successfully in
`/tmp/nwkit-mul-msc-20261002.R9pM98`, using the checkout through `PYTHONPATH`.
For its single gene family, the evaluated parent log likelihoods were:
H2 B `-1.0610520887581791`, A `-2.808494409935022`, and C
`-3.7940079775623525`. Baseline, autopolyploid and time-incompatible candidates
had explicit exclusions and blank likelihoods. The best dated tree was:

```text
[&R](((A:2,X+:2)<1>:1,(B:1,X*:1)<2>:2)<3>:2,C:5)<4>:0;
```

Serial and two-process CLI runs produced identical output bundles. Tests
covered invalid/ambiguous sampling, missing units/Ne, nonzero root stems,
empty inputs, unrooted or polytomous trees, input/output aliases and resource
limits. Failures before scoring preserved existing score/report files. The
injected late tree-output failure preserved all five previous bundle files
with no staging files left behind.

A serialization audit found that ETE's default six-digit lengths lost dated
precision. MSC output now uses the existing RADTE-style 17-significant-digit
writer without mutating global parser settings. Exact branch-float round trips
and quoted node names passed. Unrequested gene diagnostics are not collected.

Default and explicit `--score-model dl` retain the existing D+L solver and
output schemas; MSC-only arguments fail in that mode rather than being ignored.
MSC has its own `nwkit-mul-msc-likelihood-v1` schema and cannot feed the legacy
GeneGalleon GRAMPA summary consumer. This milestone changed no GeneGalleon
files/defaults and preserved its existing worktree changes.

## Provenance and Remaining Gates

Final tested source SHA-256:

```text
mul_coalescent.py       cabb89e6bed59c4e96c1c26615993310f8c84e51b58d967c19c0770d609a4173
mul_msc_model.py        52ffc8b47f36fd524a72920a2030cdef52357cf2df9eb6f84b5ed3296615e9e6
mul_msc.py              f56867ef7cbe75c552b76d34fe7cdd166aa5642c72f3a700197f2c508a50da0b
mul_reconcile.py       dd96de7b8307ebfd82c654ce7a1221e3ccd0311f7c71019e8fed9ce6115110ba
mul_reconcile_cli.py   49a13bc3af75b6725365ce69cf9ef8809cf52eee8efa3c1fc2950edecbf8c55b
```

Full/release/distribution gates, minimum-Python and other-platform runs,
Docker/SIF validation, and production biological analyses were not executed
for this host-only research milestone. There was no commit, version bump,
push, Wiki update or deployment. Historical D+L/GeneGalleon validation is
separate evidence, not validation of this new MSC mode.

The next gate is a normalized duplication-loss-locus-history model, including
ancestral copy priors, extinction and sampling/ascertainment. Only then can
no-polyploidy and allopolyploid hypotheses receive comparable likelihoods.
Adding an ILS penalty to D+L, or transferring DLCoal to an arbitrary MUL-tree
without those histories and observation rules, would not supply that model.
Parameter identifiability/estimation, estimated-tree and misspecification
benchmarks, calibrated support, and opt-in GeneGalleon consumers/runtime checks
remain before adopting MSC as a production default.
