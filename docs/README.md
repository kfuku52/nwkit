# NWKIT documentation

The repository documentation is split by purpose so the root stays focused on
project entry points and release metadata.

## Guides

User-facing command, data-format, model, and mathematical documentation lives
in [`guides/`](guides/):

- [ASR and ancestral-trait models](guides/ASR.md)
- [CLI and TSV conventions](guides/CLI_TSV_CONVENTIONS.md)
- [Phylogenetic regression](guides/PHYLOGENETIC_REGRESSION.md)
- [RADTE dating](guides/RADTE.md)
- [SHIFT inference](guides/SHIFT.md)
- [DTT, PCA, signal, and stochastic maps](guides/DTT.md)
- [All optimal reconciliation roots and loss exports](guides/RECONCILIATION_EXPORTS.md)

The complete guide list is available in the directory; filenames retain the
historical names used by the CLI and examples.

## Validation and research notes

Reproducibility studies, calibration experiments, performance measurements, and
adoption decisions live in [`validation/`](validation/). These documents record
the evidence and limits for experimental features; they are not a replacement
for the corresponding user guide.

## Project operations

- [Development checks](../DEVELOPMENT.md)
- [Release checklist](../RELEASING.md)
- [Change history](../CHANGELOG.md)
