# NWKIT documentation

The repository documentation is split by purpose so the root stays focused on
project entry points and release metadata.

## Guides

User-facing command, data-format, model, and mathematical documentation lives
in [`guides/`](guides/):

- [All optimal reconciliation roots and loss exports](guides/RECONCILIATION_EXPORTS.md)
- [ASR and ancestral-trait models](guides/ASR.md)
- [Bounded conditional MUL/MSC estimation](guides/MUL_MSC_FIT.md)
- [CLI and TSV conventions](guides/CLI_TSV_CONVENTIONS.md)
- [Conditional allopolyploid MUL/MSC prototype](guides/MUL_MSC.md)
- [Discrete stochastic maps](guides/STOCHASTIC_MAPS.md)
- [Disparity through time](guides/DTT.md)
- [Experimental locus DL + ILS null comparison](guides/MUL_LOCUS_MC.md)
- [Fixed gene-topology DL/WGD likelihood](guides/WGD_TREE.md)
- [Focal-lineage Ks correction](guides/KSRATE.md)
- [Installation check and first analysis](guides/QUICK_START.md)
- [Native gene-count DL/WGM candidates](guides/WGD_COUNT.md)
- [Native MUL-tree reconciliation](guides/MUL_RECONCILE.md)
- [Phylogenetic PCA](guides/PCA.md)
- [Phylogenetic regression](guides/PHYLOGENETIC_REGRESSION.md)
- [Phylogenetic signal](guides/SIGNAL.md)
- [RADTE dating](guides/RADTE.md)
- [SHIFT inference](guides/SHIFT.md)
- [Tree format conversion](guides/CONVERT.md)

The complete guide list is available in the directory; filenames retain the
historical names used by the CLI and examples.

## Validation and research notes

Reproducibility studies, calibration experiments, performance measurements, and
adoption decisions live in [`validation/`](validation/). These documents record
the evidence and limits for experimental features; they are not a replacement
for the corresponding user guide.

- [Native WGD and Ks scientific audit](validation/WGD_SCIENTIFIC_VALIDATION.md)
- [MUL-tree exact and original-GRAMPA validation](validation/MUL_RECONCILE_VALIDATION.md)
- [Conditional allopolyploid MUL/MSC prototype validation](validation/MUL_MSC_VALIDATION.md)
- [Conditional MUL/MSC estimation and estimated-tree pilot](validation/MUL_MSC_FIT_VALIDATION.md)
- [Locus DL + ILS comparison and search-wide calibration](validation/MUL_LOCUS_MC_VALIDATION.md)
- [Locus MC numerical and reproducibility audit](validation/MUL_LOCUS_MC_AUDIT.md)
- [Finite null grid calibration and independent evaluation](validation/MUL_LOCUS_GRID_CALIBRATION.md)
- [Finite null grid input and runner audit](validation/MUL_LOCUS_GRID_AUDIT.md)
- [Conditional locus integration support and precision probe](validation/MUL_LOCUS_INTEGRATION_PROBE.md)
- [Conditional reference, KL intervals and GeneGalleon integration](validation/MUL_LOCUS_CONDITIONAL_INTEGRATION.md)
- [Node diagnostics and fixed-data MC budget follow-up](validation/MUL_NODE_DIAGNOSTICS_AND_MC_BUDGET.md)

## Project operations

- [Development checks](../DEVELOPMENT.md)
- [Release checklist](../RELEASING.md)
- [Change history](../CHANGELOG.md)
