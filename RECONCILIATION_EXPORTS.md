# Reusable reconciliation exports

`nwkit root --method reconciliation --candidates-out roots.nwk` writes all
exactly equally optimal physical root edges as one Newick tree per line.
`--outfile` still receives the selected canonical root. The collection uses the
same output format, name quoting, branch annotations and pairwise distances as
the selected tree. Rooting preserves the existing position when its split
already matches; rerooted edges use the usual midpoint placement. Candidate order follows the deterministic canonical
ordering of reconciliation rooting. There is no candidate-count cap.

```sh
nwkit root --method reconciliation --infile gene.nwk \
  --species-tree species.nwk --duplication-cost 1.5 --loss-cost 1 \
  --outfile selected.nwk --candidates-out roots.nwk
nwkit reconcile --infile selected.nwk --species-tree species.nwk \
  --event-source lca --unmatched error --outfile reconciliation.tsv
```

Costs retain their NWKIT defaults of 1 and 1 unless explicitly set.
`--candidates-out` requires a file path and reconciliation rooting. File outputs
are published as a recoverable bundle: validation or handled write failures
leave prior files intact. When `--outfile -` is used, the candidate file is
published before the selected tree is printed. Input/output collisions are
rejected, including collisions with the species tree or mapping TSV.

`reconcile` adds the nullable integer column `implied_losses`. At each mapped
internal LCA event, a duplication contributes the entire topological distance
to each child placement; a speciation contributes each distance minus one.
The sum counts missing species lineages below the gene root. No losses are
inferred on the path from the species root to the gene root. Leaves contribute
zero. Other event sources and unmapped events have undefined loss counts.
Existing clade IDs, event definitions and RADTE-compatible columns are retained.
