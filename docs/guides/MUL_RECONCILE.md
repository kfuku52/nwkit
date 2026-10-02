# Native MUL-Tree Reconciliation

`nwkit mul-reconcile` replaces the GRAMPA-style search over one duplicated
species clade (H1) and a second-parent attachment branch (H2). It evaluates
the singly labeled baseline and auto-/allopolyploid MUL-tree hypotheses using
exact minimum duplication-plus-loss (D+L) reconciliation. This is a parsimony
ranking, not a posterior probability or calibrated polyploidy test.

```bash
nwkit mul-reconcile -i rooted_genes.nwk --species-tree species.nwk \
  --species-regex '.*_([^_]+)$' --h1 'x,y,z' \
  -o scores.tsv --report mappings.tsv --check-out checks.tsv \
  --tree-out best.nwk --model-out model.json --cpus 4
```

Gene input is a Newick collection. Both gene and species trees must be rooted
and strictly binary; gene tips and singly labeled species tips must be unique.
The shared species parser, regex and mapping TSV associate gene tips with
species. Branch lengths and input internal labels do not affect D+L scores.
Explicit unrooted declarations fail; an unmarked binary root retains the
usual NWKIT rooted interpretation. No tree is silently pruned or resolved.

## Hypotheses

H1's subtree is copied as a new sibling of H2. H2 cannot be strictly below H1;
H1=H2 represents an autopolyploid placement. H1/H2 accept exact monophyletic
comma-separated species sets, tip names, or one-based postorder internal
numbers, matching GRAMPA's closing-parenthesis numbering. Space separates
multiple selectors; omission searches all admissible nodes. Baseline ID 0 is
always included. Numeric species labels, `<n>` tip labels and species labels
ending in `+` or `*` are reserved in search mode.

The original H1 occurrences receive `+` and copied ones `*` in output trees;
gene tips match their original species identity. Candidate IDs are run-local,
not stable biological identifiers. Equal-score candidates sort by candidate
ID, with the baseline preferred when it ties. Every tied candidate remains in
the score table and `best_hypotheses` metadata. The detailed report exports
every optimal assignment for the **first** best candidate only.

`--multree yes` instead reconciles against a supplied rooted binary tree with
repeated, identical species tip names. It disables H1/H2 generation and allows
more general multiple-copy topologies; no copy labels should be added to those
input species names. This is not an inferred multi-event history.

## Exact Objective

For a fixed leaf-occurrence assignment, an internal gene node maps to the LCA
of its child maps. It is a duplication if that map equals either child map.
Each child edge contributes `depth(child_map) - depth(parent_map) - 1 + dup`
losses. The mapped gene-root depth adds unobserved ancestral lineages, as in
GRAMPA. Thus this objective differs from `nwkit reconcile`'s below-gene-root
loss count. All edges have unit topological length for scoring.

Dynamic-program states retain the minimum subtree score for each possible
species-node map, tie counts, and backpointers. The recurrence combines every
pair of child states; higher-cost alternatives at the same map cannot improve
any ancestral score. Traceback reports all minimizing leaf assignments without
grouping distinct gene tips or fixing assignments using sister heuristics.
Traceback is iterative, including deeply nested binary trees. Detailed mapping
rows are written incrementally into staged output rather than accumulated in RAM.
All hypotheses use all input gene trees. Neither missing species nor numerous
ambiguous copies triggers candidate-dependent filtering.

`--max-candidates`, `--max-state-pairs` (per gene/candidate), and `--max-maps`
(per reported gene) bound resources. Exceeding a bound fails the analysis before
replacing outputs; maps are never silently truncated. With no detailed report,
`--max-maps` does not limit scoring. These limits are not scientific filters.

## Outputs and Compatibility

The score and detail tables retain GRAMPA's modern column names. `maps` retains
the annotated-Newick `gene-node[species-map-duplication]` serialization.
Nonportable characters in the bracketed species map are percent-escaped to
avoid breaking Newick comments; JSON retains the exact original labels.
Details add `node.maps` and `node.maps.format=nwkit-node-map-json-v1`: JSON
node rows contain gene and MUL-tree indices/labels, duplications, child-edge
losses and gene-root losses. Node
indices are zero-based postorder for genes and preorder for MUL-trees. Summing the
duplication and both loss fields reconstructs the reported D+L totals.

Check-table `groups` means ambiguous tips, `fixed` is zero, and `combinations`
is the exact number of all leaf-occurrence assignments. These are diagnostics,
not GRAMPA's collapsed-group counts. `optimal.mappings` and `state.pairs` add
exact solver diagnostics. The model JSON records method, tie handling,
filtering policy, score table, and limitations. Primary TSV alone can use
stdout; related files are staged together and cannot replace input paths.

## Interpretation and Validation

The model conditions on supplied roots/topologies and omits sequence likelihood,
ILS, HGT and rooting uncertainty. Lower D+L alone does not establish WGD or
identify a real donor. It is separate from `wgd-count`/`wgd-tree`; those likelihood
models still do not fit an allopolyploid donor network.

The implementation follows the D+L objective and MUL placement described by
[Thomas et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28419377/) and the
[GRAMPA documentation](https://gwct.bio/grampa/readme.html), not its grouping
heuristics. See [executed reference and exhaustive checks](../validation/MUL_RECONCILE_VALIDATION.md):
an identified whole-root auto example scores 13 under complete enumeration
but 15 in the grouped original program. This difference is intentional and
verified against the original program's own LCA scoring function.
