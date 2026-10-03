# Conditional Reference, KL Intervals And GeneGalleon Integration

Executed 2026-10-03 on local NWKIT 0.43.37 with the prior uncommitted
locus/MSC work preserved. This is opt-in research integration, not a
production-model adoption or a rewrite of previous validation records.

## Repairs And Contracts

- Added two-sided unit-interval Chernoff-KL inversion for weighted moments.
  The independent non-Bernoulli finite-sample enumeration respects its
  error budget. Conservative outward brackets handle extreme means;
  non-integer/nonfinite sample counts fail rather than return spurious bounds.
  The last sample-count validation guard was added after the study freeze;
  all stored study moments have integer counts and every completed weighted
  observed fit still replays exactly with the final implementation.
- Added explicit CLI JSON integration `detection-rb` and `hybrid-rb`.
  Full hidden daughter conditioning, prior-weighted selection normalization,
  whole-universe/stratum error allocation and resource failures are retained.
  Weighted reports leave `hits` blank; saved stratum moments, populations
  and interval methods reproduce scores. Both modes give identical serial
  and multiprocess output bytes. No pseudocounts or failed-family dropping.
- Built a separate topology-only conditional reference with dense pure-death
  matrix exponentials and uniform pair histories. It uses the independent
  global SSA locus simulator but imports no production coalescent or detection
  helpers. Joint full-forest and rare-bound tests precede the study.
  It cannot simulate sequence/tree-inference error or branch lengths.
- The first execution exposed tiny negative roundoff in impossible
  lineage-count increases from `expm`. Only exact pure-death structural zeros
  are restored; a negative feasible transition still fails. The first complete
  63 failure records and source snapshot are preserved, not repaired in place.

Bounds follow Theorem 3 of [Foong, Bruinsma and Burt (2022)](https://arxiv.org/pdf/2205.07880).
The daughter-conditioned construction follows the
[DLCoal model](https://compbio.mit.edu/dlcoal/pub/dlcoal/doc/dlcoal-manual.html).
See [the predeclared protocol](../../examples/mul-locus/CONDITIONAL.md).

## Frozen Follow-Up

Directory: `/tmp/nwkit-locus-conditional-probe-v2-20261003`.
Three scenarios, seven truth cases each, one dataset per case, ten families,
2,000 raw draws per bank, 19 null replicates per each of four null grid points.
Hidden exact-integration limit four; work cap 100,000. Original biological
settings/decision rule were not tuned. Seeds and sources froze before sampling.

| Method | Completed / Planned | Score-Support Failures | Reference Failures |
| --- | ---: | ---: | ---: |
| Histogram | 1 / 21 | 20 | 0 |
| Detection RB | 21 / 21 | 0 | 0 |
| Hybrid RB | 21 / 21 | 0 | 0 |

The runner exits 1 because histogram failures remain. It retains all 63 planned
method analyses. All 3,268 complete null searches (43 analyses times 76)
had their point P, MC-overlap endpoints and grid supremum independently
recomputed. All 42 weighted observed fits exactly replay from saved moments.
Both first/second source archives match their frozen source hashes.

For each weighted method, all six alternative datasets have point P <=0.05;
five select the correct parent. Missing-ils allop-B selects the wrong parent.
No alternative passes the upper-MC decision rule: upper P is 1 for every
completed analysis. No on-grid null reports an event, but only one dataset
per truth point is available: no empirical 5% guarantee or power claim.
This new low-budget probe does not reproduce or eliminate the nine failures
of the historical 100,000-draw experiment.

The following comparisons recompute **both** intervals on the same stored
banks, observations and alpha. They isolate interval construction from the
changed seeds/budget/reference; the point probabilities are exactly equal.
Cells are dependent, not independent trials.

| Scenario / Method | Pattern-Bank Cells | Median KL Width | Median Prior EB Width |
| --- | ---: | ---: | ---: |
| Baseline / Detection | 820 | 0.081268 | 0.117592 |
| Baseline / Hybrid | 820 | 0.081221 | 0.109775 |
| Turnover / Detection | 980 | 0.083388 | 0.134307 |
| Turnover / Hybrid | 980 | 0.081660 | 0.122555 |
| Missing-ILS / Detection | 960 | 0.107970 | 0.128763 |
| Missing-ILS / Hybrid | 960 | 0.107229 | 0.119389 |

These interval improvements do not establish biological accuracy, parent
identifiability or a resolved event decision. D+L remains the default.

## GeneGalleon Integration

All three existing BUSCO DNA/protein/orthogroup stage contracts optionally
add the locus analysis via `grampa_locus_model`, with explicit generation-scaled
`grampa_locus_species_tree`, donor scope and null-calibration settings.
The generation tree must match the original species/rooted topology; child
order is aligned to preserve numeric postorder selectors. Model detection
names are normalized with the species/gene labels, without modifying inputs.
Both model/tree must be outside the result directory.

Legacy D+L files/columns are unchanged. The separate `locus/` directory and
existing nine outputs publish in one recoverable bundle. Invalid inputs,
unsupported observed support, missing generation units, configuration errors
and both early/late publication failures preserve previous files. A cleared
option retains older experiments with a warning, not as a new locus analysis.
Source/model/tree/config/output changes participate in artifact invalidation.

No scientific thresholds, biological rates, curated inputs, upstream branch
defaults or default model were changed. Workflow curation, size defaults
(5-50) and IQ-TREE/rooting error are not silently accommodated by the model.
See [the GeneGalleon guide](https://github.com/kfuku52/genegalleon/blob/main/docs/grampa-replacement.md).

## Executed Checks

NWKIT runtime: isolated CPython 3.12.14/macOS ARM64; compiled SciPy imports,
ETE and package preflight pass, with no broken requirements. Broad quick lane
includes all locus/MSC/coalescent/reconciliation/reference tests and CLI,
rooting, tree I/O, provenance and output-transaction consumers:

```sh
PY=/tmp/nwkit-wgd-dev-20261001/bin/python
NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  "$PY" tools/check.py quick -- tests/test_mul_locus*.py \
  tests/test_mul_coalescent.py tests/test_mul_msc.py tests/test_mul_msc_fit.py \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py \
  tests/test_output_transaction.py tests/test_rooting_state.py \
  tests/test_provenance.py tests/test_util_tree_io.py tests/test_tree_outputs.py -x -rs
"$PY" tools/check.py quick -- tests/test_mul_locus_integral.py -x -rs
"$PY" tools/check.py test -- -m slow tests/test_mul_locus.py \
  tests/test_mul_locus_study.py tests/test_mul_msc_fit.py tests/test_mul_msc.py -rs
```

- Broad quick: 628 passed, 15 slow deselected, no skips. Ruff/format and mypy
  (250 sources) pass.
- Final bounded-count guard: focused quick 59 passed, adding ten distinct
  boundary cases to the broad selection, with static checks passing again.
- Explicit slow lane: 15 passed, 164 other tests deselected, no skips.
  **653 distinct selected tests**, not the full suite.
- Research scripts pass Ruff; maintainability reports 3,638 functions,
  mean 6.91, maximum 50 with existing exception; all hard limits pass.

GeneGalleon runtime: Linux ARM64, CPython 3.12.14, Docker local-source overlay
`local/genegalleon:locus-integration-20261003`, immutable final image
`sha256:808188429bc65df5adbc802122425b81bc21d1768dae99141ef07ccbcfbfca53`.
It inherits the current verified upstream/container snapshot and records the
local NWKIT source hashes separately; default upstream SHAs were not pinned.
The normal wrapper's daily freshness check passes, never disabled.

```sh
GG_TEST_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:locus-integration-20261003 \
  bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_native_mul_reconcile_runtime.py -x -rs
GG_TEST_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:locus-integration-20261003 \
  bash ./dev check static workflow/tests/test_shell_static_safety.py \
  workflow/tests/test_shell_entrypoint_static.py \
  workflow/tests/test_artifact_provenance.py workflow/tests/test_parse_grampa.py \
  workflow/tests/test_species_labeling_helpers.py -x -rs
bash ./dev config-check
```

The static selection has 311 passes, no skips. Native wrapper validation has
27 passes, no skips before the additional late-publication case; the final
affected selection has 12 passes, 16 deselected and no skips, adding one
distinct late-publication case. **339 distinct GeneGalleon tests pass**.
Container `bash ./dev lint` also
passes with a per-command Git safe-directory setting. Host lint is blocked
by macOS Bash 3.2; the first container attempt lacked Ruff, which was added
only to the local validation overlay. No global Git trust was modified.

Full release/security/coverage/distribution gates, Python 3.10, other platforms,
SIF, full scientific workflow runs and empirical inference-error calibration
were not executed. Local README/Wiki feature descriptions were updated and
command lists remain alphabetical. No commit, push or Wiki deployment.

```text
protocol.json: c3ff8357eb4c155b443aed1785bb9e687929cefdbd1acff250221a4a1177aeb0
summary.tsv: 852e2bd53c277da77a26e4603acbdaeaa6bbad23b2628057eedc0925dd59b5ea
source-snapshot.tar: 63800cbe7c69f0656787c010dae400b332d57f5e0da006c38e824ca6f8e90e42
```
