# Reconciliation and RADTE figure examples

The small example in the parent directory has two gene speciation nodes, S1
and S2, sharing the species event AB at age 10. Its symmetric-rate duplication
has age 20. It is a deterministic illustration, not an accuracy benchmark.
Commands for drawing both reports are in [RECONCILIATION_PLOTS.md](../../../RECONCILIATION_PLOTS.md).

`generax-gene.nhx` and `generax-species.nwk` are copied from
`kfuku52/RADTE/data/example_generax_01/gene_tree.nhx` and `species_tree.nwk`.
They contain Arabidopsis thaliana, Capsella rubella, and Utricularia gibba gene
identifiers and GeneRax annotations. No historical RADTE expected output is
used as a numerical reference; the following commands generate fresh NWKIT
results from those inputs.

Run from the NWKIT repository root:

```sh
nwkit reconcile --infile examples/radte/visualization/generax-gene.nhx \
  --species-tree examples/radte/visualization/generax-species.nwk \
  --event-source nhx --outfile output/pdf/reconcile-generax.tsv \
  --figure-out output/pdf/reconcile-generax.pdf

nwkit radte --generax-nhx examples/radte/visualization/generax-gene.nhx \
  --species-tree examples/radte/visualization/generax-species.nwk \
  --max-age 1000 --uncertainty profile --out-prefix output/pdf/radte-generax \
  --figure-out output/pdf/radte-generax.pdf

nwkit draw --radte-prefix output/pdf/radte-generax \
  --species-tree examples/radte/visualization/generax-species.nwk \
  --outfile output/pdf/radte-generax.svg

nwkit draw --infile examples/radte/visualization/generax-gene.nhx \
  --species-tree examples/radte/visualization/generax-species.nwk \
  --reconciliation output/pdf/reconcile-generax.tsv \
  --outfile output/pdf/reconcile-generax.svg
```

The figures visualize the current experimental estimator and its recorded
diagnostics. They do not establish accuracy or interval coverage on biological
data. The input time units are preserved without assigning a different unit.

The biological example also exercises profile continuation with feasible initial
ages. Its intervals are conditional on the supplied tree and calibrations.
