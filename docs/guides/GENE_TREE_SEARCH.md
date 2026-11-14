# Targeted gene-tree search

`nwkit gene-tree-search` proposes complete-tip gene trees from automatically
detected sets of potentially anomalous tips. It optionally compares them with
an independently refitted input topology using GeneRax EVAL. A detection is a
structural hypothesis, not a diagnosis of a misplaced gene. Genuine ancient
duplications, differential loss, ILS, species-tree error and sequence/model
error can also produce the detected patterns.

## Detection and joint moves

Both inputs must be rooted, strictly binary trees with unique nonempty tip
labels. Species labels use the shared parser or `--species-map-tsv`; every gene
species must be represented in the species tree.

For a gene node whose two child subtrees share species, the detector finds the
complete set of copies of each overlapping species on either side. It combines
one side choice for **every** overlapping species. Thus, a two-species overlap
directly generates two-tip or larger covers; a failed singleton search cannot
exclude them. Neither side is presumed erroneous.

LCA duplications also produce cuts separating the two species branches below
the mapped ancestor. These detect non-overlap duplications that species overlap
alone misses. Whole duplication child clades and unions of related constraints
(shared event or shared species) provide additional proposals. The inferred tip
count is determined by each set. `--max-moved-tips` is a computational ceiling,
not a requested move count or a biological threshold.

For each set, maximal monophyletic selected components are detached together.
Every tip is then reattached: one two-tip clade requires one pruned component,
whereas two separated tips require two. Candidate endpoints can cross a valley
where an individual move would make the score worse. A beam ranks partial
regrafts by rooted LCA duplication-plus-loss count without an improvement gate.
The score counts losses below the gene root; unobserved ancestral stem losses
are excluded. Temporary removal is diagnostic only and never changes the data
used for likelihood comparison.

The topology induced by each detached component's original tips remains fixed.
A monophyletic two-tip set therefore moves as one clade; it does not also generate
independent moves of those two tips. With multiple components, later insertions
may attach within an earlier regraft; that component need not remain monophyletic
with respect to the other inserted tips. Coupling considers pairwise unions of the best
`--max-proposals` seed sets, not recursive unions of all subsets. The backbone
must retain at least two tips. These restrictions apply even with very large
beam/evaluation budgets; broad-budget sensitivity is not full tree enumeration.

Detected sets are ranked by diagnostic D+L reduction per removed tip. Candidate
retention and evaluation reserve an endpoint for each represented set before
filling by D+L rank. Limits can omit sets, attachments and combinations; the JSON
report records oversized sets, cover enumeration truncation, beam discards,
coupling seed exclusions, insufficient backbones, candidate discards and
unevaluated retained trees. Increase budgets to inspect
sensitivity. No budgeted result establishes a globally optimal topology.

## Propose and inspect

```sh
nwkit gene-tree-search -i gene.nwk --species-tree species.nwk \
  --species-parser taxonomic \
  -o scores.tsv --sets-out sets.tsv --candidates-out candidates.tsv \
  --report-out report.json --tree-out best.nwk
```

Without `--evaluation generax`, the best output retains the input topology.
D+L reduction alone cannot authorize a repair. The score table includes the
baseline, moved tips as JSON arrays, tip/component counts and structural scores;
likelihood fields remain empty and `evaluation_status=not_evaluated`.

`sets.tsv` records each set, the stable constraint clade IDs, proposal sources
and diagnostic pruning gain. `candidates.tsv` records `candidate_id` and Newick.
Only the primary score table accepts stdout. Companion outputs are staged
together and cannot overwrite any input. Output directories must already
exist. Output trees carry tip labels, branch lengths and rooting declarations;
internal names, old event annotations and support values are omitted. Recompute
support and reconciliation annotations for an adopted topology.

Gene/species trees also accept inline Newick or `-` for stdin; only one input
may own stdin. File hashes are captured before search, and inline/stdin tree
hashes describe the text actually parsed. Add `--audit audit.jsonl` to record
all outputs and a stdin species-map hash; the command JSON flags that map hash
as unavailable without the shared audit record.

## Fit the same full-tip alignment

```sh
nwkit gene-tree-search -i gene.nwk --species-tree species.nwk \
  --species-parser taxonomic --evaluation generax \
  --alignment alignment.fa.gz --subst-model GTR+G4 --rec-model UndatedDL \
  --generax-command 'mpiexec -np 4 generax' --workdir new-eval-directory \
  --eval-rounds 2 --root-policy optimize \
  -o scores.tsv --sets-out sets.tsv --candidates-out candidates.tsv \
  --report-out report.json --tree-out best.nwk
```

The substitution model is mandatory. FASTA IDs must exactly match every tree
tip; extra/missing IDs, duplicates, unequal lengths and recognized all-missing
sequences are rejected. One normalized FASTA is retained for all evaluations.
GeneRax taxon-label restrictions apply to this backend.

Every candidate uses the same alignment, species tree, substitution model,
reconciliation model and optimization protocol. Each is presented as a separate
family with `--per-family-rates`; rates are independently optimized, never
shared across artificial copies of the family. Sequence and reconciliation log
likelihoods are added with reconciliation weight 1. Models `UndatedDL` and
`UndatedDTL` are supported; the **proposal heuristic remains D+L** for both.
The DTL model therefore needs broader-budget sensitivity checks when transfer
may account for a discordance.

EVAL preserves unrooted topology. `--root-policy optimize` permits the GeneRax
reconciliation root fit and deduplicates unrooted-equivalent proposals;
`keep` adds `--enforce-gene-tree-root` and compares rooted proposals. Rooting
metadata and explicit overrides follow the shared
[CLI conventions](CLI_TSV_CONVENTIONS.md). This does not reproduce a historical
MAD-weighted optimization objective: both baseline and candidates are freshly
fit under the selected EVAL policy. See the
[GeneRax documentation](https://github.com/BenoitMorel/GeneRax/wiki/GeneRax).

The default two EVAL rounds give every topology the same refit opportunity.
The best joint fit across all its rounds is retained, including for the baseline;
subsequent rounds start with those best fitted branches. Inputs, commands, logs,
results and scores for every round are retained. These are GeneRax local fits,
not a proof of globally optimized nuisance parameters. Inspect the round scores
and increase `--eval-rounds` if the baseline/candidate ranking is unstable.

GeneRax stats can be rounded. The score table includes a conservative bound
from the displayed decimal places. A topology replaces the baseline only when
its score interval lies wholly above the baseline interval. A small unresolved
difference keeps the baseline. Among candidates with resolved improvements,
the largest joint score is selected, even if another candidate has a higher but
unresolved point estimate. The `fitted_*` event-count columns refer to the
optimized output root; the unprefixed structural columns refer to the proposed
root. A selected candidate is best within the evaluated set and the chosen
model; its gain is not a P-value or a probability of misplacement.

The work directory must be new and independent of input/report paths. The
generated family file uses relative paths from each recorded round's working
directory, so spaces, quotes or `#` in ancestor directories do not corrupt its
syntax. Commands in the JSON report must be replayed from their recorded `cwd`.
The launcher is parsed without a shell. Backend failure, timeout, missing results,
nonfinite likelihoods or changed topology/tips abort the operation and preserve
existing report outputs; diagnostics remain in the work directory. Container
MPI settings belong to the invocation environment. For the GeneGalleon Docker
runtime, use the documented isolated launcher environment (`OMPI_MCA_plm=isolated`)
and the Open MPI root permission environment when running as root.

Use an inline substitution-model string, such as `GTR+G4` or `LG+G4`.
FASTA IDs must match tips exactly, without descriptions; equal nonzero sequence
lengths, duplicate IDs and missing/extra tips are checked before evaluation.
An all-unknown DNA sequence is rejected for recognized nucleotide models;
protein `N` represents observed asparagine. GeneRax performs the remaining
model/alphabet checks. Evaluation requires at least three gene tips. Fitted
trees must be rooted, strictly binary, preserve the chosen topology and tips,
and have finite nonnegative branch lengths. Local POSIX worker processes are
stopped with the launcher on failure, timeout or cancellation; remote MPI jobs and
non-POSIX worker cleanup remain the launcher's responsibility.

## BMI1 regression input

[`tests/data/gene_tree_search/bmi1`](../../tests/data/gene_tree_search/bmi1)
contains the supplied GeneRax gene tree, species tree and trimmed CDS alignment:
79 tips, 27 species-tree tips and 3,273 sites. Its `PROVENANCE.json` records the
source members and SHA-256 hashes. It is a real-input fixture with an unknown
correct topology, not ground truth. The detector identifies the previously
discussed Nymphaea tip without an ID rule; renaming all gene tips gives the same
structural detection.

For a conditional exhaustive scan of the **automatically top-ranked set**, use
the fixture inputs above and `--max-proposals 1 --beam-width 256
--max-candidates 256 --max-evaluations 256`. This covers all single-component
attachment positions for BMI1, but intentionally leaves other detected sets
out. The default multiset search includes coupled moves from the outset.

Focused checks:

```sh
python tools/check.py quick -- tests/test_gene_tree_search.py tests/test_cli.py \
  tests/test_cli_contracts.py tests/test_interface_conventions.py tests/test_reconcile.py
python tools/check.py test -- tests/test_gene_tree_search.py -m slow -rs
```

The slow check runs real GeneRax when available and otherwise reports a skip.
The small reference enumerates all 105 rooted five-tip topologies, compares
proposal coverage, and checks every D+L count against the shared reconciliation
export. Additional tests cover joint sets, one-component versus multi-component
moves, true ancient duplication, deterministic detection, input protection,
failure rollback, likelihood gating and score rounding.
