# BMI1 targeted-search validation, 2026-10-06

The standalone [`gene-tree-search`](../guides/GENE_TREE_SEARCH.md) implementation
was validated using the Maruyama query2family data received on 2026-10-06.
This study tests detection and conditional topology comparison. The supplied
gene tree is not known biological ground truth.

The numerical search results below were recorded before the publication review.
They retain their original search ordering and final-round protocol; they have
not been silently regenerated after fixes to canonical candidate metadata and
best-round retention. The current command reports the best joint fit per
topology and can change finite-budget candidate ordering. The fixture, models
and all-tip comparison contract are unchanged.

## Inputs and protocol

The committed regression fixture is
[`tests/data/gene_tree_search/bmi1`](../../tests/data/gene_tree_search/bmi1).
`PROVENANCE.json` identifies the source members and SHA-256 hashes:

- Gene tree: `generax_nwk/PRC1BMI1_generax.nwk` (79 tips).
- Species tree: `parameters/undated_species_tree.pruned.nwk` (27 tips).
- Alignment: `clipkit/PRC1BMI1_cds.clipkit.fa.gz` (79 sequences, 3,273 nucleotide sites).
- Received-run model: GTR+G4, UndatedDL, MAD rooting.

Comparisons below use fresh GeneRax EVAL fits: GTR+G4, UndatedDL, independently
optimized per-family rates, reconciliation weight 1, optimized reconciliation
roots, seed 12345, and **two identical refit rounds per topology**. They compare
the input topology against candidates under the same EVAL conditions, rather
than against an archived MAD-weighted score. Unrooted-equivalent candidates
are deduplicated; fixed-root and DTL modes were checked separately.

Scientific execution used GeneRax 2.1.3 in
`local/genegalleon:busco-guide-20261006`, image ID
`sha256:59d846d984ffdf31111f586d9b990926125154974242cc8215a14da1d86c04b7`,
with the current NWKIT source mounted as the working directory. The GeneGalleon
freshness wrapper rejected this image manifest with `Unknown source in runtime
manifest: rapidnj`. The explicitly identified snapshot was therefore used for
standalone scientific checks; GeneGalleon workflow integration, daily freshness
and SIF compatibility are not established by these results.

## Detection and complete-tip searches

No gene identifier is specified to the detector. Its highest-ranked set is the
single tip `Nymphaea_colorata_GeneID116261861`. Input-root LCA counts are 34
duplications and 39 losses; diagnostic removal of this tip decreases their sum
by 6. Renaming every gene tip, with explicit unchanged species mappings, gives
the same leading structural anomaly. Removal is never used to evaluate a
likelihood: every evaluated topology retains all 79 sequences and all sites.

The multiset run retained 32 detected sets and evaluated 64 unique unrooted
topologies. Detection considered 8,278 cover states without exhausting the
20,000-state budget; 162 admissible sets were discovered and 130 fell outside
the selected-set budget. Larger sets were excluded by the eight-tip ceiling.
The regraft search generated 27,072 states, discarded 26,320 at beam boundaries,
found 238 endpoint topologies, and retained 64. These are substantial search
restrictions, explicitly reported by the tool.

The conditional exhaustive run selected the **automatically top-ranked set**,
without naming Nymphaea, and considered every backbone attachment position.
It generated 155 rooted regraft states and retained/evaluated 153 unique
unrooted topologies, including the baseline. No beam or endpoint truncation
occurred for this set. Other detected sets remain outside this conditional
search; this is not an exhaustive search over all gene trees or all move sets.

## Results

The independently refitted baseline has sequence log likelihood -83,999.0,
reconciliation log likelihood -185.853 and joint log likelihood -84,184.853.
LCA counts below refer to the fitted output root.

| Search/candidate | Sequence log likelihood | Reconciliation log likelihood | Joint gain over baseline | D / L |
|---|---:|---:|---:|---:|
| Multiset search, Anthoceros single-tip move | -84,003.6 | -178.674 | +2.579 | 33 / 36 |
| All Nymphaea placements, best local move | -84,000.2 | -184.099 | +0.554 | 34 / 38 |
| Best Nymphaea attachment to an angiosperm-only sister subtree | -84,036.9 | -176.480 | -28.527 | 34 / 36 |

The multiset best moves
`Anthoceros_agrestis_AagrBONN_evm.TU.Sc2ySwM_228.4798`.
Its gain is a conditional model comparison, not evidence of a known biological
error or of a globally optimal replacement. The fitted species-overlap count
actually changes from 31 to 32 despite the improved D+L/joint score; reducing
species-overlap events is not equivalent to maximizing the reconciliation
likelihood or establishing orthology.

The Nymphaea-only best leaves the tip associated with the gymnosperm subtree.
Its joint gain is small. The tool accounts for displayed score rounding:
the baseline bound is +/-0.5005 and this candidate's bound is +/-0.0505, so the
gain only narrowly exceeds that conservative bound. No statistical support
claim follows from it.

For the descriptive angiosperm check, sister-subtree taxa in each **proposed**
rooted tree were mapped to the species-tree clade containing Amborella and
Arabidopsis. Even the best placement with an entirely angiosperm sister subtree
has lower joint likelihood by 28.527. This input/model comparison therefore
does **not** support simply moving the highlighted Nymphaea gene into that
group. The correct biological placement remains unresolved; alignment/model
assessment and independent evidence are needed before adopting such a change.
Neither sequence exclusion nor a paralog-to-ortholog label change was used.

## Reproduce and inspect

From the NWKIT root, create a new output directory, then run the following
inside the identified container or an equivalent environment with GeneRax:

```sh
python -m nwkit gene-tree-search \
  -i tests/data/gene_tree_search/bmi1/gene.nwk \
  --species-tree tests/data/gene_tree_search/bmi1/species.nwk \
  --species-parser taxonomic --evaluation generax \
  --alignment tests/data/gene_tree_search/bmi1/alignment.fa.gz \
  --subst-model GTR+G4 --rec-model UndatedDL \
  --generax-command 'mpiexec -np 4 generax' --eval-rounds 2 \
  --max-proposals 32 --beam-width 8 --max-candidates 64 --max-evaluations 64 \
  --workdir NEW_OUTPUT/generax --tree-out NEW_OUTPUT/best.nwk \
  --sets-out NEW_OUTPUT/sets.tsv --candidates-out NEW_OUTPUT/candidates.tsv \
  --report-out NEW_OUTPUT/report.json -o NEW_OUTPUT/scores.tsv
```

For all attachments of the automatic top set, replace the four search budgets
with `--max-proposals 1 --beam-width 256 --max-candidates 256
--max-evaluations 256` and use a new work directory. The actual conditional run
used two MPI ranks; all candidates and the baseline used the same launcher.
For Docker root execution, set `OMPI_MCA_plm=isolated`,
`OMPI_ALLOW_RUN_AS_ROOT=1` and `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` in the container
environment. Mount the NWKIT source at its absolute host path and use that path
as the container working directory so `python -m nwkit` runs the current source.

The local evidence directory is
`output/gene-tree-search-20261006-144613/`:

- `final-scores.tsv`, `final-sets.tsv`, `final-candidates.tsv`, `final-best.nwk`,
  and `final-report.json`: final multiset run.
- `final-all-target-scores.tsv`, corresponding candidates/sets/best/report
  files: final 153-topology conditional exhaustive run.
- `all-target-placement-summary.tsv`: descriptive sister-subtree classification.
- `final-multiset-v2/` and `final-all-target/`: complete per-round GeneRax
  inputs, logs, trees and stats.
- `validation-environment.json`: container identity, source hashes and check
  results. Earlier exploratory directories retain preliminary and interrupted
  runs; the `final-*` files above are the reported evidence.

## Executed checks

Host checks used an isolated Python 3.12.14 environment, preserving the existing
checkout environment. A same-version ETE 4.4.0 source rebuild repaired an
incompatible wheel in the isolated environment. Import preflight and `pip check`
passed. No dependency defaults or scientific model settings were changed.

```sh
python tools/check.py quick -- tests/test_gene_tree_search.py tests/test_cli.py \
  tests/test_cli_contracts.py tests/test_interface_conventions.py tests/test_reconcile.py -x
python tools/check.py test -- tests/test_provenance.py tests/test_output_transaction.py -rs
python tools/check.py test -- tests/test_gene_tree_search.py -m slow -rs
```

- Host quick checks: Ruff lint/format and mypy passed; **163 tests passed**,
  four real-backend cases deselected by the documented quick lane.
- Host provenance/output-transaction checks: **58 tests passed**.
- Docker real GeneRax: **four cases passed**, testing DL/DTL crossed with
  kept/optimized roots, with two refit rounds each; no optional backend skips.
- The five-tip reference enumerates all **105 rooted / 15 unrooted** topologies,
  verifies search coverage, and matches shared reconciliation D+L exports.
- New-module Bandit checks, maintainability limits and `git diff --check` passed.

At the initial scientific-validation stage, the full repository/release suite,
performance benchmarking and SIF validation were not run. These focused checks
cover the changed command, scientific
counting/search invariants, CLI contracts, atomic output behavior and the actual
GeneRax adapter; they are not a release or calibrated misplacement classifier.

## Publication review, 0.43.40

The review corrected best-fit retention across rounds, selection among resolved
score improvements, canonical move metadata, protein-asparagine validation,
inline/stdin tree hashes, output/audit collision protection, and local POSIX
worker cleanup on timeout, cancellation or launcher failure. GeneRax family
files and native arguments now use round-relative
paths with a recorded `cwd`; directories containing spaces, quotes and `#` were
tested with the real backend. These fixes do not change scientific thresholds,
fixture sequences or substitution/reconciliation models.

Validation used isolated Python 3.12.14 in a Docker development snapshot
`sha256:0f6f281409c34c45d7183e6fca642cba828f9bc0e1274727a451a3732c92e696`,
derived from the GeneGalleon image identified above, with the checkout mounted
at `/repo`. It preserved the existing checkout environment and prior build
artifacts. The focused quick command initially passed Ruff, mypy and **239 tests**;
four real-backend cases were deselected by this lane. After deep-tree and
launcher-failure fixes, the selection was extended with MUL and count-WGM model,
CLI, fitting and scientific-reference tests: **400 passed, 11 slow cases
deselected** in 93.90 seconds. It includes a 1,201-tip serialization/copy check
and a child-ready handshake that reproduces orphaned workers after nonzero
launcher exit before the fix. The explicit slow command passed **all four**
DL/DTL × kept/optimized-root cases
with the unusual directory names and two refit rounds each.

The final isolated numerical runtime uses NumPy 2.5.3 and SciPy 1.17.1. Release
review reproduced a finite-input SciPy 1.18.0/1.18.1 exponential hang outside
NWKIT; those two versions are excluded, with a reproducer and removal condition
in [DEVELOPMENT.md](../../DEVELOPMENT.md). The scientific search bounds and
models remain intact. Review also fixed count-WGM selection probabilities
exceeding one through complement roundoff, using equivalent algebra rather
than clipping probabilities. Independent probability-one tests failed in eight
cases before this correction and now pass, alongside the existing rare-Yule
and dense-transition references.

The final `python tools/check.py release` gate passed on this Docker runtime:
Ruff checks/formatting, uncached mypy, dependency consistency, Bandit,
dependency vulnerability audit, **4,846 tests passed / 82 skipped** in
4,327.30 seconds, **86% combined statement/branch coverage**, maintainability
limits, and reproducible wheel/sdist contents. The four real GeneRax cases
also passed within this full run. Skips cover optional R/kfl1ou and independent
research/reference backends, plus one case-sensitive-filesystem check on the
mounted checkout. Eight warnings report Python's multiprocessing `fork()`
deprecation; no test failed. The final validation-document update was followed
by another `python tools/check.py dist` check of the publication bytes.
SIF validation and a complete external CI platform matrix were not run locally.

A separate Linux x86-64 / Python 3.11.17 container with NumPy 2.4.6 and SciPy
1.17.1 passed **146 affected tests, four deselected** in 19.96 seconds, including
the deep MUL-tree annotation that failed in the preceding public CI run.
Three preceding CI failures were also checked in isolation: deep annotation,
rare-pulse likelihood and frozen SHIFT replay all passed on this runtime.
The recorded SHIFT replay failure was not reproduced; no historical archive
or replay tolerance was changed.

A separate real CLI check used the committed BMI1 inputs, GTR+G4/UndatedDL,
optimized roots, seed 12345, two rounds and budgets of one proposal, beam width
two, three candidates and three evaluations. It retained all 79 tips/3,273
sites and selected the baseline: the two D+L-favored candidates had joint gains
-223.789 and -59.589. Their rejection verifies the likelihood adoption gate;
this deliberately narrow scan does not revise the earlier exhaustive results.
The fitted output topology, input hashes, two-round metadata and all five
companion-output audit records were independently checked.

A real BMI1 timeout check also returned a failure, published no score table,
retained its diagnostic work directory and left no local GeneRax/MPI workers
running after cleanup. The synthetic child-process regression separately
verifies that killing only a launcher cannot leave its worker active.

Per-round commands must now be replayed from their recorded `cwd`. The current
score is the best joint fit across all rounds, not necessarily the final round.
The earlier numerical tables remain historical results with their stated
protocol rather than being presented as a regenerated release experiment.
